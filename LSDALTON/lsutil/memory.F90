MODULE memory_handling
use,intrinsic :: iso_c_binding,only:c_ptr,c_f_pointer,c_associated,c_null_ptr,c_loc
use infpar_module
use background_buffer_module
use TYPEDEFTYPE, only : mpi_realk
use AO_typetype
use OD_typetype
use molecule_typetype
use LSmatrix_type
use Matrix_module
use precision
use LSTENSOR_TYPETYPE
use basis_typetype
use dec_typedef_module
use OverlapType
use tensor_type_def_module
#ifdef MOD_UNRELEASED
use lattice_type
#endif
#ifdef VAR_MPI
use lsmpi_module
#endif
private
public DetectHeapMemoryInit,DetectHeapMemory
public Set_PrintSCFmemory,MemInGB,stats_globalmem
public Set_MemModParamPrintMemory, Print_Memory_info
public MemModParamPrintMemorylupri,MemModParamPrintMemory
public Set_PrintMemoryLowerLimit
public get_available_memory
public init_globalmemvar
public init_threadmemvar
public collect_thread_memory
public mem_TurnOffThread_Memory,mem_TurnONThread_Memory
public stats_mem,stats_mpi_mem,stats_mem_tp,scf_stats_debug_mem
public set_mem_xccalc,debug_mem_stats,print_maxmem
!  print_mem_alloc
public mem_InsideOMPsection
public mem_alloc
public mem_dealloc
public mem_pseudo_alloc
public mem_pseudo_dealloc
public mem_init_background_alloc
public mem_change_background_alloc
public mem_free_background_alloc
public mem_is_background_buf_init,mem_get_bg_buf_n,mem_get_bg_buf_free
public mem_allocated_mem_real, mem_deallocated_mem_real
public mem_allocated_global,mem_allocated_type_matrix
!parameters
public mem_realsize,mem_logicalsize,mem_complexsize,mem_intsize,mem_shortintsize
! special routines to keep track of mem in special types
public mem_allocated_mem_intwork
public mem_deallocated_mem_intwork
public mem_allocated_mem_overlap
public mem_deallocated_mem_overlap
public mem_allocated_mem_overlapT
public mem_deallocated_mem_overlapT
public mem_allocated_mem_linkshell
public mem_deallocated_mem_linkshell
public mem_allocated_mem_integralitem
public mem_deallocated_mem_integralitem
public mem_allocated_mem_integrand
public mem_deallocated_mem_integrand
public mem_allocated_mem_ODitem
public mem_deallocated_mem_ODitem
public mem_allocated_mem_FMM
public mem_deallocated_mem_FMM
public mem_allocated_mem_lstensor
public mem_deallocated_mem_lstensor
public mem_allocated_lstensor
public mem_allocated_mem_type_matrix
public mem_deallocated_mem_type_matrix
public copy_from_mem_stats
public copy_to_mem_stats
public max_mem_used_global
#ifdef MOD_UNRELEASED
public mem_allocated_lvec_data
public mem_allocated_mem_lvec_data
public mem_deallocated_mem_lvec_data
public mem_allocated_lattice_cell
public mem_allocated_mem_lattice_cell
public mem_deallocated_mem_lattice_cell
#endif
public mem_add_external_memory
public longintbuffersize
logical,save :: PrintSCFmemory
logical,save :: MemModParamPrintMemory
integer,save :: MemModParamPrintMemorylupri
integer(kind=8),save :: PrintMemoryLowerLimit
integer,save :: longintbuffersize
!Monitor memory for integral code and possibly other parts!
!GLOBAL VARIABLES
integer(KIND=long),save :: mem_allocated_global, max_mem_used_global         !Count all memory
integer(KIND=long),save :: mem_allocated_Heap_global, max_mem_used_Heap_global !Count all memory to determine Heap Memory Usage
integer(KIND=long),save :: mem_allocated_type_matrix, max_mem_used_type_matrix !Count memory, density opt code
integer(KIND=long),save :: mem_allocated_type_matrix_MPIFULL, max_mem_used_type_matrix_MPIFULL !Count memory across MPI nodes
integer(KIND=long),save :: mem_allocated_real, max_mem_used_real             !Count 'real' memory, integral code
integer(KIND=long),save :: mem_allocated_mpi, max_mem_used_mpi             !Count memory,allocated with the MPI_MEM_ALLOC routine 
integer(KIND=long),save :: mem_allocated_complex, max_mem_used_complex             !Count 'real' memory, integral code
integer(KIND=long),save :: mem_allocated_integer, max_mem_used_integer       !Count 'integer' memory, integral code
integer(KIND=long),save :: mem_allocated_character, max_mem_used_character   !Count 'character' memory, integral code
integer(KIND=long),save :: mem_allocated_logical, max_mem_used_logical       !Count 'logical' memory, integral code

integer(KIND=long),save :: mem_allocated_AOBATCH, max_mem_used_AOBATCH       !Count 'AOBATCH' memory, integral code
integer(KIND=long),save :: mem_allocated_XCcalc, max_mem_used_XCcalc       !Count XC memory, integral code
integer(KIND=long),save :: mem_allocated_DECORBITAL, max_mem_used_DECORBITAL       !Count 'DECORBITAL' memory, deccc code
integer(KIND=long),save :: mem_allocated_DECFRAG, max_mem_used_DECFRAG       !Count 'DECFRAG' memory, deccc code
integer(KIND=long),save :: mem_allocated_BATCHTOORB, max_mem_used_BATCHTOORB       !Count 'BATCHTOORB' memory, deccc code
integer(KIND=long),save :: mem_allocated_DecAObatchinfo, max_mem_used_DecAObatchinfo !Count 'DecAObatchinfo' memory, deccc code
integer(KIND=long),save :: mem_allocated_MYPOINTER, max_mem_used_MYPOINTER       !Count 'MYPOINTER' memory, deccc code
integer(KIND=long),save :: mem_allocated_fragmentAOS, max_mem_used_fragmentAOS       !Count 'fragmentAOS' memory, deccc code
integer(KIND=long),save :: mem_allocated_ARRAY2, max_mem_used_ARRAY2       !Count 'ARRAY2' memory, deccc code
integer(KIND=long),save :: mem_allocated_ARRAY4, max_mem_used_ARRAY4       !Count 'ARRAY4' memory, deccc code
integer(KIND=long),save :: mem_allocated_ARRAY, max_mem_used_ARRAY         !Count 'ARRAY ' memory, deccc code
integer(KIND=long),save :: mem_allocated_PNOSPACEINFO, max_mem_used_PNOSPACEINFO         !Count 'PNOSpaceInfo ' memory, deccc code
integer(KIND=long),save :: mem_allocated_MP2DENS, max_mem_used_MP2DENS       !Count 'MP2DENS' memory, deccc code
integer(KIND=long),save :: mem_allocated_TRACEBACK, max_mem_used_TRACEBACK       !Count 'TRACEBACK' memory, deccc code
integer(KIND=long),save :: mem_allocated_MP2GRAD, max_mem_used_MP2GRAD       !Count 'MP2GRAD' memory, deccc code
integer(KIND=long),save :: mem_allocated_ODBATCH, max_mem_used_ODBATCH       !Count 'ODBATCH' memory, integral code
integer(KIND=long),save :: mem_allocated_LSAOTENSOR, max_mem_used_LSAOTENSOR       !Count 'LSAOTENSOR' memory, integral code
integer(KIND=long),save :: mem_allocated_SLSAOTENSOR, max_mem_used_SLSAOTENSOR       !Count 'SLSAOTENSOR' memory, integral code
integer(KIND=long),save :: mem_allocated_GLOBALLSAOTENSOR, max_mem_used_GLOBALLSAOTENSOR       !Count 'GLOBALLSAOTENSOR' memory, integral code
integer(KIND=long),save :: mem_allocated_ATOMTYPEITEM, max_mem_used_ATOMTYPEITEM       !Count 'ATOMTYPEITEM' memory, integral code
integer(KIND=long),save :: mem_allocated_ATOMITEM, max_mem_used_ATOMITEM       !Count 'ATOMITEM' memory, integral code
integer(KIND=long),save :: mem_allocated_LSMATRIX, max_mem_used_LSMATRIX       !Count 'LSMATRIX' memory, integral code
#ifdef MOD_UNRELEASED
integer(KIND=long),save :: mem_allocated_lvec_data, max_mem_used_lvec_data !Count memory, density opt code
integer(KIND=long),save :: mem_allocated_lattice_cell, max_mem_used_lattice_cell !Count memory, density opt code
#endif

!Memory distributed on types:
integer(KIND=long),save :: mem_allocated_linkshell, max_mem_used_linkshell         !Count memory, type linkshell
integer(KIND=long),save :: mem_allocated_integralitem, max_mem_used_integralitem   !Count memory, type integralitem
integer(KIND=long),save :: mem_allocated_integrand, max_mem_used_integrand         !Count memory, type integrand
integer(KIND=long),save :: mem_allocated_overlap, max_mem_used_overlap             !Count memory, type Overlap
integer(KIND=long),save :: mem_allocated_overlapT, max_mem_used_overlapT           !Count memory, type Overlap
integer(KIND=long),save :: mem_allocated_intwork, max_mem_used_intwork             !Count memory, integral work array
integer(KIND=long),save :: mem_allocated_ODitem, max_mem_used_ODitem               !Count memory, type ODitem
integer(KIND=long),save :: mem_allocated_FMM, max_mem_used_FMM
!Count memory, type FMM
integer(KIND=long),save :: mem_allocated_lstensor, max_mem_used_lstensor
logical,save :: mem_InsideOMPsection
!THREADPRIVATE LOCAL VARIABLES
!
! The subroutine init_threadmemvar should be called before an OMP section to init these variables and
! after the OMP section the subroutine collect_thread_memory should be called. Note that the
! max_mem_used will here be an estimate and not the absolute correct value. TK
!
integer(KIND=long),save :: mem_tp_allocated_global, max_mem_tp_used_global         !Count all memory
integer(KIND=long),save :: mem_tp_allocated_Heap_global, max_mem_tp_used_Heap_global !Count all memory to determine Heap Memory Usage
integer(KIND=long),save :: mem_tp_allocated_type_matrix, max_mem_tp_used_type_matrix !Count memory, density opt code
integer(KIND=long),save :: mem_tp_allocated_type_matrix_MPIFULL, max_mem_tp_used_type_matrix_MPIFULL
integer(KIND=long),save :: mem_tp_allocated_real, max_mem_tp_used_real             !Count 'real' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_mpi, max_mem_tp_used_mpi             !Count memory, MPI_ALLOC_MEM
integer(KIND=long),save :: mem_tp_allocated_complex, max_mem_tp_used_complex             !Count 'complex' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_integer, max_mem_tp_used_integer       !Count 'integer' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_character, max_mem_tp_used_character   !Count 'character' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_logical, max_mem_tp_used_logical       !Count 'logical' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_AOBATCH, max_mem_tp_used_AOBATCH       !Count 'AOBATCH' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_DECORBITAL, max_mem_tp_used_DECORBITAL       !Count 'DECORBITAL' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_DECFRAG, max_mem_tp_used_DECFRAG       !Count 'DECFRAG' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_BATCHTOORB, max_mem_tp_used_BATCHTOORB       !Count 'BATCHTOORB' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_DecAObatchinfo, max_mem_tp_used_DecAObatchinfo !Count 'DecAObatchinfo' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_MYPOINTER, max_mem_tp_used_MYPOINTER       !Count 'MYPOINTER' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_fragmentAOS, max_mem_tp_used_fragmentAOS       !Count 'fragmentAOS' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_ARRAY2, max_mem_tp_used_ARRAY2       !Count 'ARRAY2' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_ARRAY4, max_mem_tp_used_ARRAY4       !Count 'ARRAY4' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_ARRAY, max_mem_tp_used_ARRAY         !Count 'ARRAY' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_PNOSPACEINFO, max_mem_tp_used_PNOSPACEINFO         !Count 'PNOSPACEINFO' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_MP2DENS, max_mem_tp_used_MP2DENS       !Count 'MP2DENS' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_TRACEBACK, max_mem_tp_used_TRACEBACK       !Count 'TRACEBACK' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_MP2GRAD, max_mem_tp_used_MP2GRAD       !Count 'MP2GRAD' memory, deccc code
integer(KIND=long),save :: mem_tp_allocated_ODBATCH, max_mem_tp_used_ODBATCH       !Count 'ODBATCH' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_LSAOTENSOR, max_mem_tp_used_LSAOTENSOR       !Count 'LSAOTENSOR' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_SLSAOTENSOR, max_mem_tp_used_SLSAOTENSOR       !Count 'SLSAOTENSOR' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_GLOBALLSAOTENSOR, max_mem_tp_used_GLOBALLSAOTENSOR       !Count 'GLOBALLSAOTENSOR' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_ATOMTYPEITEM, max_mem_tp_used_ATOMTYPEITEM       !Count 'ATOMTYPEITEM' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_ATOMITEM, max_mem_tp_used_ATOMITEM       !Count 'ATOMITEM' memory, integral code
integer(KIND=long),save :: mem_tp_allocated_LSMATRIX, max_mem_tp_used_LSMATRIX       !Count 'LSMATRIX' memory, integral code
!Memory distributed on types:
integer(KIND=long),save :: mem_tp_allocated_linkshell, max_mem_tp_used_linkshell         !Count memory, type linkshell
integer(KIND=long),save :: mem_tp_allocated_integralitem, max_mem_tp_used_integralitem   !Count memory, type integralitem
integer(KIND=long),save :: mem_tp_allocated_integrand, max_mem_tp_used_integrand         !Count memory, type integrand
integer(KIND=long),save :: mem_tp_allocated_overlap, max_mem_tp_used_overlap             !Count memory, type Overlap
integer(KIND=long),save :: mem_tp_allocated_overlapT, max_mem_tp_used_overlapT           !Count memory, type OverlapT
integer(KIND=long),save :: mem_tp_allocated_intwork, max_mem_tp_used_intwork             !Count memory, type Overlap
integer(KIND=long),save :: mem_tp_allocated_ODitem, max_mem_tp_used_ODitem               !Count memory, type ODitem
!Count memory, type FMM
integer(KIND=long),save :: mem_tp_allocated_FMM, max_mem_tp_used_FMM
integer(KIND=long),save :: mem_tp_allocated_lstensor, max_mem_tp_used_lstensor
! integers used to collect OMP thread mem !NOT THREAD PRIVATE - BUT SHARED
integer(KIND=long),save :: max_mem_used_global_tmp,max_mem_used_type_matrix_tmp
integer(KIND=long),save :: max_mem_used_Heap_global_tmp
integer(KIND=long),save :: max_mem_used_type_matrix_MPIFULL_tmp,max_mem_used_real_tmp
integer(KIND=long),save :: max_mem_used_complex_tmp
integer(KIND=long),save :: max_mem_used_mpi_tmp
integer(KIND=long),save :: max_mem_used_integer_tmp,max_mem_used_character_tmp
integer(KIND=long),save :: max_mem_used_logical_tmp,max_mem_used_linkshell_tmp
integer(KIND=long),save :: max_mem_used_AOBATCH_tmp,max_mem_used_ODBATCH_tmp
integer(KIND=long),save :: max_mem_used_DECORBITAL_tmp
integer(KIND=long),save :: max_mem_used_DECFRAG_tmp
integer(KIND=long),save :: max_mem_used_BATCHTOORB_tmp
integer(KIND=long),save :: max_mem_used_DecAObatchinfo_tmp
integer(KIND=long),save :: max_mem_used_MYPOINTER_tmp
integer(KIND=long),save :: max_mem_used_fragmentAOS_tmp
integer(KIND=long),save :: max_mem_used_ARRAY2_tmp
integer(KIND=long),save :: max_mem_used_ARRAY4_tmp
integer(KIND=long),save :: max_mem_used_tensor_tmp
integer(KIND=long),save :: max_mem_used_PNOSPACEINFO_tmp
integer(KIND=long),save :: max_mem_used_MP2DENS_tmp
integer(KIND=long),save :: max_mem_used_TRACEBACK_tmp
integer(KIND=long),save :: max_mem_used_MP2GRAD_tmp
integer(KIND=long),save :: max_mem_used_LSAOTENSOR_tmp,max_mem_used_SLSAOTENSOR_tmp
integer(KIND=long),save :: max_mem_used_GLOBALLSAOTENSOR_tmp,max_mem_used_ATOMTYPEITEM_tmp
integer(KIND=long),save :: max_mem_used_ATOMITEM_tmp,max_mem_used_LSMATRIX_tmp

integer(KIND=long),save :: max_mem_used_integralitem_tmp,max_mem_used_integrand_tmp
integer(KIND=long),save :: max_mem_used_overlap_tmp,max_mem_used_ODitem_tmp
integer(KIND=long),save :: max_mem_used_FMM_tmp,max_mem_used_lstensor_tmp
integer(KIND=long),save :: max_mem_used_intwork_tmp,max_mem_used_overlapT_tmp
#ifdef MOD_UNRELEASED
integer(KIND=long),save :: mem_tp_allocated_lvec_data, max_mem_tp_used_lvec_data
integer(KIND=long),save :: max_mem_used_lvec_data_tmp
integer(KIND=long),save :: mem_tp_allocated_lattice_cell, max_mem_tp_used_lattice_cell
integer(KIND=long),save :: max_mem_used_lattice_cell_tmp
#endif
!Memory PARAMETERS
!sizes of types found by SIZEOF() !which only works for some compilers
!so these numbers are then hardcoded which requires that when
!the types change these numbers should be changes accordingly.
integer(KIND=long),save :: mem_OVERLAPsize
integer(KIND=long),save :: mem_AOBATCHsize
integer(KIND=long),save :: mem_DECORBITALsize
integer(KIND=long),save :: mem_DECFRAGsize
integer(KIND=long),save :: mem_BATCHTOORBsize
integer(KIND=long),save :: mem_DecAObatchinfosize
integer(KIND=long),save :: mem_MYPOINTERsize
integer(KIND=long),save :: mem_fragmentAOSsize
integer(KIND=long),save :: mem_ARRAY2size
integer(KIND=long),save :: mem_ARRAY4size
integer(KIND=long),save :: mem_ARRAYsize
integer(KIND=long),save :: mem_PNOSPACEINFOsize
integer(KIND=long),save :: mem_MP2DENSsize
integer(KIND=long),save :: mem_TRACEBACKsize
integer(KIND=long),save :: mem_MP2GRADsize
integer(KIND=long),save :: mem_ODBATCHsize
integer(KIND=long),save :: mem_LSAOTENSORsize
integer(KIND=long),save :: mem_SLSAOTENSORsize
integer(KIND=long),save :: mem_GLOBALLSAOTENSORsize
integer(KIND=long),save :: mem_ATOMTYPEITEMsize
integer(KIND=long),save :: mem_ATOMITEMsize
integer(KIND=long),save :: mem_LSMATRIXsize
integer(KIND=long),save :: mem_MATRIXsize
integer(KIND=long),parameter :: mem_realsize=8
integer(KIND=long),parameter :: mem_logicalsize=4
integer(KIND=long),parameter :: mem_complexsize=16
!integer(KIND=long),parameter :: mem_charsize=4 It is not obvious how we
!should count character memory!
!integer(KIND=long),parameter :: mem_pointersize=8
#if VAR_INT64
integer(KIND=long),parameter :: mem_intsize=8_long
#else
integer(KIND=long),parameter :: mem_intsize=4_long
#endif
integer(KIND=long),parameter :: mem_int4size=4_long
integer(KIND=long),parameter :: mem_int8size=8_long
integer(KIND=long),parameter :: mem_shortintsize=1_long
#ifdef MOD_UNRELEASED
integer(KIND=long),save :: mem_lvec_datasize
integer(KIND=long),save :: mem_lattice_cellsize
#endif
!$OMP THREADPRIVATE(mem_tp_allocated_global, max_mem_tp_used_global,&
!$OMP mem_tp_allocated_Heap_global, max_mem_tp_used_Heap_global,&
!$OMP mem_tp_allocated_type_matrix, max_mem_tp_used_type_matrix,&
!$OMP mem_tp_allocated_type_matrix_MPIFULL, max_mem_tp_used_type_matrix_MPIFULL,&
!$OMP mem_tp_allocated_real, max_mem_tp_used_real,&
!$OMP mem_tp_allocated_complex, max_mem_tp_used_complex,&
!$OMP mem_tp_allocated_integer, max_mem_tp_used_integer,&
!$OMP mem_tp_allocated_character, max_mem_tp_used_character,&
!$OMP mem_tp_allocated_logical, max_mem_tp_used_logical,&
!$OMP mem_tp_allocated_AOBATCH, max_mem_tp_used_AOBATCH,&
!$OMP mem_tp_allocated_DECORBITAL, max_mem_tp_used_DECORBITAL,&
!$OMP mem_tp_allocated_DECFRAG, max_mem_tp_used_DECFRAG,&
!$OMP mem_tp_allocated_BATCHTOORB, max_mem_tp_used_BATCHTOORB,&
!$OMP mem_tp_allocated_DecAObatchinfo, max_mem_tp_used_DecAObatchinfo,&
!$OMP mem_tp_allocated_MYPOINTER, max_mem_tp_used_MYPOINTER,&
!$OMP mem_tp_allocated_fragmentAOS, max_mem_tp_used_fragmentAOS,&
!$OMP mem_tp_allocated_ARRAY2, max_mem_tp_used_ARRAY2,&
!$OMP mem_tp_allocated_ARRAY4, max_mem_tp_used_ARRAY4,&
!$OMP mem_tp_allocated_ARRAY, max_mem_tp_used_ARRAY,&
!$OMP mem_tp_allocated_PNOSPACEINFO, max_mem_tp_used_PNOSPACEINFO,&
!$OMP mem_tp_allocated_mpi, max_mem_tp_used_mpi, &
!$OMP mem_tp_allocated_MP2DENS, max_mem_tp_used_MP2DENS,&
!$OMP mem_tp_allocated_TRACEBACK, max_mem_tp_used_TRACEBACK,&
!$OMP mem_tp_allocated_MP2GRAD, max_mem_tp_used_MP2GRAD,&
!$OMP mem_tp_allocated_ODBATCH, max_mem_tp_used_ODBATCH,&
!$OMP mem_tp_allocated_LSAOTENSOR, max_mem_tp_used_LSAOTENSOR,&
!$OMP mem_tp_allocated_SLSAOTENSOR, max_mem_tp_used_SLSAOTENSOR,&
!$OMP mem_tp_allocated_GLOBALLSAOTENSOR, max_mem_tp_used_GLOBALLSAOTENSOR,&
!$OMP mem_tp_allocated_ATOMTYPEITEM, max_mem_tp_used_ATOMTYPEITEM,&
!$OMP mem_tp_allocated_ATOMITEM, max_mem_tp_used_ATOMITEM,&
!$OMP mem_tp_allocated_LSMATRIX, max_mem_tp_used_LSMATRIX,&
!$OMP mem_tp_allocated_linkshell, max_mem_tp_used_linkshell,&
!$OMP mem_tp_allocated_integralitem, max_mem_tp_used_integralitem,&
!$OMP mem_tp_allocated_integrand, max_mem_tp_used_integrand,&
!$OMP mem_tp_allocated_intwork, max_mem_tp_used_intwork,&
!$OMP mem_tp_allocated_overlap, max_mem_tp_used_overlap,&
!$OMP mem_tp_allocated_overlapT, max_mem_tp_used_overlapT,&  
!$OMP mem_tp_allocated_ODitem, max_mem_tp_used_ODitem,&
!$OMP mem_tp_allocated_FMM, max_mem_tp_used_FMM,&
#ifdef MOD_UNRELEASED
!$OMP mem_tp_allocated_lvec_data, max_mem_tp_used_lvec_data,&
!$OMP mem_tp_allocated_lattice_cell, max_mem_tp_used_lattice_cell,&
#endif
!$OMP mem_tp_allocated_lstensor, max_mem_tp_used_lstensor)

!Interfaces for allocating/deallocating pointers
INTERFACE mem_alloc
   MODULE PROCEDURE real_allocate_1dim, real_allocate_1dim_int64, &
      &             real_allocate_1dim_sp, real_allocate_2dim, &
      &             real_allocate_2dim_sp, real_allocate_2dim_zero, real_allocate_3dim, &
      &             real_allocate_3dim_sp, real_allocate_3dim_zero, real_allocate_4dim, &
      &             real_allocate_5dim, real_allocate_5dim_zero, &
      &             real_allocate_7dim_zero, &
      &             complex_allocate_1dim, complex_allocate_2dim, &
      &             intS_allocate_1dim,intS_allocate_1dim_wrapper4, &
      &             int4_allocate_1dim,int4_allocate_1dim_wrapper4, &
      &             int8_allocate_1dim,int8_allocate_1dim_wrapper4, &
      &             int4_allocate_2dim,int4_allocate_2dim_wrapper4, &
      &             int8_allocate_2dim,int8_allocate_2dim_wrapper4, &
      &             int_allocate_3dim, &
      &             int4_allocate_4dim,int4_allocate_4dim_wrapper4,&
      &             int8_allocate_4dim,int8_allocate_4dim_wrapper4,&
      &             int_allocate_1dim_zero, int_allocate_2dim_zero,  &
      &             int_allocate_3dim_zero, int_allocate_4dim_zero, &
      &             char_allocate_1dim_wrapper4,char_allocate_1dim,shortint_allocate_2dim, &
      !     &             shortint_allocate_1dim, &
   &             logic4_allocate_1dim,logic4_allocate_1dim_wrapper4,logic4_allocate_1dim_zero, &
      &             logic4_allocate_2dim, logic4_allocate_3dim,&
      &             logic8_allocate_1dim,logic8_allocate_1dim_wrapper4,logic8_allocate_1dim_zero, &
      &             logic8_allocate_2dim, logic8_allocate_3dim,&
      &             AOBATCH_allocate_1dim, ODBATCH_allocate_1dim, DECORBITAL_allocate_1dim, &
      &             LSAOTENSOR_allocate_1dim, SLSAOTENSOR_allocate_1dim, &
      &             GLOBALLSAOTENSOR_allocate_1dim, ATOMTYPEITEM_allocate_1dim, &
      &             ATOMITEM_allocate_1dim, LSMATRIX_allocate_1dim,&! LSMATRIXP_allocate_1dim, &
      &             MATRIX_allocate_1dim, MATRIXP_allocate_1dim, DECFRAG_allocate_1dim, &
      &             BATCHTOORB_allocate_1dim,MYPOINTER_allocate_1dim, MYPOINTER_allocate_2dim, &
      &             fragmentAOS_allocate_1dim, &
      &             ARRAY2_allocate_1dim,ARRAY4_allocate_1dim,MP2DENS_allocate_1dim, &
      &             TRACEBACK_allocate_1dim,MP2GRAD_allocate_1dim,PNOSPACEINFO_allocate_1dim, &
      &             OVERLAPT_allocate_1dim,tensor_allocate_1dim, lsmpi_allocate_i8V, lsmpi_allocate_i4V,&
      &             lsmpi_allocate_dV4,lsmpi_allocate_dV8, lsmpi_local_allocate_dV8, &
      &             lsmpi_local_allocate_I8V8,lsmpi_local_allocate_I4V4,&
      &             lsmpi_allocate_d,DECAOBATCHINFO_allocate_1dim,&
#ifdef MOD_UNRELEASED
      &             lvec_data_allocate_1dim, lattice_cell_allocate_1dim, &
#endif
      &             lsmpi_local_allocate_dV4
   END INTERFACE
   !
   INTERFACE mem_dealloc
      MODULE PROCEDURE real_deallocate_1dim, real_deallocate_2dim,  &
         &             real_deallocate_1dim_sp, real_deallocate_2dim_sp, real_deallocate_3dim_sp,&
         &             real_deallocate_3dim, real_deallocate_4dim, &
         &             real_deallocate_5dim, real_deallocate_7dim, &
         &             complex_deallocate_1dim, complex_deallocate_2dim, &
         &             int4_deallocate_1dim, int8_deallocate_1dim,  &
         &             int4_deallocate_2dim, int8_deallocate_2dim,  &
         &             int4_deallocate_3dim, int8_deallocate_3dim, &
         &             int4_deallocate_4dim, int8_deallocate_4dim, &
         &             char_deallocate_1dim, shortint_deallocate_2dim, &
         &             shortint_deallocate_1dim, &
         &             logic4_deallocate_1dim, logic4_deallocate_2dim,logic4_deallocate_3dim, &
         &             logic8_deallocate_1dim, logic8_deallocate_2dim,logic8_deallocate_3dim, &
         &             AOBATCH_deallocate_1dim, ODBATCH_deallocate_1dim, DECORBITAL_deallocate_1dim, &
         &             LSAOTENSOR_deallocate_1dim, SLSAOTENSOR_deallocate_1dim, &
         &             GLOBALLSAOTENSOR_deallocate_1dim, ATOMTYPEITEM_deallocate_1dim, &
         &             ATOMITEM_deallocate_1dim, LSMATRIX_deallocate_1dim, &!LSMATRIXP_deallocate_1dim, &
         &             MATRIX_deallocate_1dim, MATRIXP_deallocate_1dim,DECFRAG_deallocate_1dim, &
         &             BATCHTOORB_deallocate_1dim,MYPOINTER_deallocate_1dim,MYPOINTER_deallocate_2dim, &
         &             fragmentAOS_deallocate_1dim, &
         &             ARRAY2_deallocate_1dim,ARRAY4_deallocate_1dim,MP2DENS_deallocate_1dim, &
         &             TRACEBACK_deallocate_1dim,MP2GRAD_deallocate_1dim, &
         &             OVERLAPT_deallocate_1dim,tensor_deallocate_1dim,PNOSPACEINFO_deallocate_1dim,&
         &             lsmpi_local_deallocate_I4V,lsmpi_local_deallocate_I8V,&
         &             lsmpi_deallocate_d,DECAOBATCHINFO_deallocate_1dim,&
#ifdef MOD_UNRELEASED
         &             lvec_data_deallocate_1dim,lattice_cell_deallocate_1dim, &
#endif
         &             lsmpi_deallocate_i8V,lsmpi_deallocate_i4V,lsmpi_deallocate_dV,lsmpi_local_deallocate_dV
   END INTERFACE

   INTERFACE MemInGB
      MODULE PROCEDURE Mem4dimInGB,Mem4dimInGBint8,&
           & Mem2dimInGB,Mem2dimInGBint8,&
           & Mem1dimInGB,Mem1dimInGBint8
   END INTERFACE MemInGB
   interface mem_pseudo_alloc
      module procedure mem_pseudo_alloc_realk,mem_pseudo_alloc_realk2,mem_pseudo_alloc_realk3,mem_pseudo_alloc_mpirealk
   end interface mem_pseudo_alloc
   interface mem_pseudo_dealloc
      module procedure mem_pseudo_dealloc_realk,mem_pseudo_dealloc_realk2,mem_pseudo_dealloc_realk3,mem_pseudo_dealloc_mpirealk
   end interface mem_pseudo_dealloc


   CONTAINS

     real(realk) function Mem4dimInGB(n1,n2,n3,n4)
       implicit none
       integer(kind=4),intent(in)  :: n1,n2,n3,n4
       integer(kind=8) :: Ln1,Ln2,Ln3,Ln4
       Ln1 = n1; Ln2 = n2; Ln3 = n3; Ln4 = n4
       Mem4dimInGB = Mem4dimInGBint8(Ln1,Ln2,Ln3,Ln4)
     end function Mem4dimInGB

     real(realk) function Mem4dimInGBint8(n1,n2,n3,n4)
       implicit none
       integer(kind=8),intent(in)  :: n1,n2,n3,n4
       Mem4dimInGBint8 = n1*n2*n3*n4*8.0E-9_realk
     end function Mem4dimInGBint8

     real(realk) function Mem2dimInGB(n1,n2)
       implicit none
       integer(kind=4),intent(in)  :: n1,n2
       integer(kind=8) :: Ln1,Ln2
       Ln1 = n1; Ln2 = n2
       Mem2dimInGB = Mem2dimInGBint8(Ln1,Ln2)
     end function Mem2dimInGB

     real(realk) function Mem2dimInGBint8(n1,n2)
       implicit none
       integer(kind=8),intent(in)  :: n1,n2
       Mem2dimInGBint8 = n1*n2*8.0E-9_realk
     end function Mem2dimInGBint8

     real(realk) function Mem1dimInGB(n1)
       implicit none
       integer(kind=4),intent(in)  :: n1
       integer(kind=8) :: Ln1
       Ln1 = n1
       Mem1dimInGB = Mem1dimInGBint8(Ln1)
     end function Mem1dimInGB

     real(realk) function Mem1dimInGBint8(n1)
       implicit none
       integer(kind=8),intent(in)  :: n1
       Mem1dimInGBint8 = n1*8.0E-9_realk
     end function Mem1dimInGBint8

   subroutine Set_PrintSCFmemory(inputPrintSCFmemory)
      implicit none
      logical,intent(in) :: inputPrintSCFmemory
      PrintSCFmemory = inputPrintSCFmemory
   end subroutine Set_PrintSCFmemory

   subroutine Set_MemModParamPrintMemory(inputMemModParamPrintMemory,lupri)
      implicit none
      logical,intent(in) :: inputMemModParamPrintMemory
      integer,intent(in) :: lupri
      MemModParamPrintMemory = inputMemModParamPrintMemory
      IF(MemModParamPrintMemory)THEN
         call Set_PrintSCFmemory(inputMemModParamPrintMemory)
      ENDIF
      MemModParamPrintMemorylupri = lupri
      PrintMemoryLowerLimit = 0
   end subroutine Set_MemModParamPrintMemory

   Subroutine Set_PrintMemoryLowerLimit(inputPrintMemoryLowerLimit)
      implicit none
      integer(kind=8),intent(in) :: inputPrintMemoryLowerLimit
      IF(MemModParamPrintMemory)THEN
         PrintMemoryLowerLimit = inputPrintMemoryLowerLimit
      ENDIF
   End Subroutine Set_PrintMemoryLowerLimit

   subroutine set_sizes_of_types()
      implicit none
      TYPE(OVERLAP) :: OVERLAPitem
      TYPE(AOBATCH) :: AOBATCHitem
      TYPE(DECORBITAL) :: DECORBITALitem
      TYPE(DECFRAG) :: DECFRAGitem
      TYPE(BATCHTOORB) :: BATCHTOORBitem
      TYPE(DECAOBATCHINFO) :: DECAOBATCHINFOitem
      TYPE(MYPOINTER) :: MYPOINTERitem
      TYPE(fragmentAOS) :: fragmentAOSitem
      TYPE(ARRAY2) :: ARRAY2item
      TYPE(ARRAY4) :: ARRAY4item
      type(tensor) :: ARRAYitem
      TYPE(PNOSpaceInfo) :: PNOSpaceitem
      TYPE(MP2DENS) :: MP2DENSitem
      TYPE(TRACEBACK) :: TRACEBACKitem
      TYPE(MP2GRAD) :: MP2GRADitem
      TYPE(ODBATCH) :: ODBATCHitem
      TYPE(LSAOTENSOR) :: LSAOTENSORitem
      TYPE(SLSAOTENSOR) :: SLSAOTENSORitem
      TYPE(GLOBALLSAOTENSOR) :: GLOBALLSAOTENSORitem
      TYPE(ATOMTYPEITEM) :: ATOMTYPEITEMitem
      TYPE(ATOMITEM) :: ATOMITEMitem
      TYPE(LSMATRIX) :: LSMATRIXitem
      TYPE(MATRIX) :: MATRIXitem
#ifdef MOD_UNRELEASED
      TYPE(lvec_data_t) :: lvec_dataitem
      TYPE(lattice_cell_info_t) :: lattice_cellitem
#endif
      ! Size of buffer handling for long integer buffer
      longintbuffersize = 81

#if defined (VAR_XLF) || defined (VAR_G95) || defined (VAR_CRAY)
#ifdef VAR_LSDEBUG
      print*,'DEBUG: Warning set sizes of Types Manual!'
      print*,'DEBUG: This is error prone - verify that the hardcoded sizes are up to date!'
#endif
      mem_AOBATCHsize=496
      mem_DECORBITALsize=64
      mem_DECFRAGsize=9016
      mem_BATCHTOORBsize=56
      mem_DECAOBATCHINFOsize=20
      mem_MYPOINTERsize=72
      mem_fragmentAOSsize=72
      mem_ARRAY2size=80
      mem_ARRAY4size=384
      mem_ARRAYsize=1360
      mem_PNOSpaceInfosize=1360
      mem_MP2DENSsize=504
      mem_TRACEBACKsize=12
      mem_MP2GRADsize=848
      mem_ODBATCHsize=88
      mem_LSAOTENSORsize=304
      mem_SLSAOTENSORsize=184
      mem_GLOBALLSAOTENSORsize=20
      mem_ATOMTYPEITEMsize=38264
      mem_ATOMITEMsize=216
      mem_LSMATRIXsize=56
      mem_MATRIXsize=368
      mem_OVERLAPsize=2272
#ifdef MOD_UNRELEASED
      mem_lvec_datasize=1944
      mem_lattice_cellsize=104
#endif
#else
      !implemented for VAR_PGF90 VAR_GFORTRAN VAR_IFORT we think!
      !we assume that all other compilers work with sizeof()
      mem_AOBATCHsize=sizeof(AOBATCHitem)
      mem_DECORBITALsize=sizeof(DECORBITALitem)
      mem_DECFRAGsize=sizeof(DECFRAGitem)
      mem_BATCHTOORBsize=sizeof(BATCHTOORBitem)
      mem_DECAOBATCHINFOsize=sizeof(DECAOBATCHINFOitem)
      mem_MYPOINTERsize=sizeof(MYPOINTERitem)
      mem_fragmentAOSsize=sizeof(fragmentAOSitem)
      mem_ARRAY2size=sizeof(ARRAY2item)
      mem_ARRAY4size=sizeof(ARRAY4item)
      mem_ARRAYsize=sizeof(ARRAYitem)
      mem_PNOSpaceInfosize=sizeof(PNOSpaceitem)
      mem_MP2DENSsize=sizeof(MP2DENSitem)
      mem_TRACEBACKsize=sizeof(TRACEBACKitem)
      mem_MP2GRADsize=sizeof(MP2GRADitem)
      mem_ODBATCHsize=sizeof(ODBATCHitem)
      mem_LSAOTENSORsize=sizeof(LSAOTENSORitem)
      mem_SLSAOTENSORsize=sizeof(SLSAOTENSORitem)
      mem_GLOBALLSAOTENSORsize=sizeof(GLOBALLSAOTENSORitem)
      mem_ATOMTYPEITEMsize=sizeof(ATOMTYPEITEMitem)
      mem_ATOMITEMsize=sizeof(ATOMITEMitem)
      mem_LSMATRIXsize=sizeof(LSMATRIXitem)
      mem_MATRIXsize=sizeof(MATRIXitem)
      mem_OVERLAPsize=sizeof(OVERLAPitem)
#ifdef MOD_UNRELEASED
      mem_lvec_datasize=sizeof(lvec_dataitem)
      mem_lattice_cellsize=sizeof(lattice_cellitem)

#endif

      !print *,sizeof(AOBATCHitem)
      !print *,sizeof(DECORBITALitem)
      !print *,sizeof(DECFRAGitem)
      !print *,sizeof(BATCHTOORBitem)
      !print *,sizeof(DECAOBATCHINFOitem)
      !print *,sizeof(MYPOINTERitem)
      !print *,sizeof(fragmentAOSitem)
      !print *,sizeof(ARRAY2item)
      !print *,sizeof(ARRAY4item)
      !print *,sizeof(ARRAYitem)
      !print *,sizeof(PNOSpaceitem)
      !print *,sizeof(MP2DENSitem)
      !print *,sizeof(TRACEBACKitem)
      !print *,sizeof(MP2GRADitem)
      !print *,sizeof(ODBATCHitem)
      !print *,sizeof(LSAOTENSORitem)
      !print *,sizeof(SLSAOTENSORitem)
      !print *,sizeof(GLOBALLSAOTENSORitem)
      !print *,sizeof(ATOMTYPEITEMitem)
      !print *,sizeof(ATOMITEMitem)
      !print *,sizeof(LSMATRIXitem)
      !print *,sizeof(MATRIXitem)
      !print *,sizeof(OVERLAPitem)
      !#ifdef MOD_UNRELEASED
      !print *,sizeof(lvec_dataitem)
      !print *,sizeof(lattice_cellitem)
      !#endif

#endif
   end subroutine set_sizes_of_types

   subroutine mem_TurnOffThread_Memory()
      implicit none
      mem_InsideOMPsection = .FALSE.

      IF(MemModParamPrintMemory)THEN
         IF(max_mem_used_global_tmp.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_TurnOffThread_Memory'
               call print_mem_alloc(MemModParamPrintMemorylupri,max_mem_used_global_tmp,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_TurnOffThread_Memory'
               call print_mem_alloc(6,max_mem_used_global_tmp,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_TurnOffThread_Memory')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,max_mem_used_global_tmp)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,max_mem_used_Heap_global_tmp)
      max_mem_used_type_matrix = MAX(max_mem_used_type_matrix,max_mem_used_type_matrix_tmp)
      max_mem_used_type_matrix_MPIFULL = MAX(max_mem_used_type_matrix_MPIFULL,max_mem_used_type_matrix_MPIFULL_tmp)
      max_mem_used_real = MAX(max_mem_used_real,max_mem_used_real_tmp)
      max_mem_used_mpi = MAX(max_mem_used_mpi,max_mem_used_mpi_tmp)
      max_mem_used_complex = MAX(max_mem_used_complex,max_mem_used_complex_tmp)
      max_mem_used_integer = MAX(max_mem_used_integer,max_mem_used_integer_tmp)
      max_mem_used_character = MAX(max_mem_used_character,max_mem_used_character_tmp)
      max_mem_used_logical = MAX(max_mem_used_logical,max_mem_used_logical_tmp)
      max_mem_used_AOBATCH = MAX(max_mem_used_AOBATCH,max_mem_used_AOBATCH_tmp)
      max_mem_used_DECORBITAL = MAX(max_mem_used_DECORBITAL,max_mem_used_DECORBITAL_tmp)
      max_mem_used_DECFRAG = MAX(max_mem_used_DECFRAG,max_mem_used_DECFRAG_tmp)
      max_mem_used_BATCHTOORB = MAX(max_mem_used_BATCHTOORB,max_mem_used_BATCHTOORB_tmp)
      max_mem_used_DECAOBATCHINFO = MAX(max_mem_used_DECAOBATCHINFO,max_mem_used_DECAOBATCHINFO_tmp)
      max_mem_used_MYPOINTER = MAX(max_mem_used_MYPOINTER,max_mem_used_MYPOINTER_tmp)
      max_mem_used_fragmentAOS = MAX(max_mem_used_fragmentAOS,max_mem_used_fragmentAOS_tmp)
      max_mem_used_ARRAY2 = MAX(max_mem_used_ARRAY2,max_mem_used_ARRAY2_tmp)
      max_mem_used_ARRAY4 = MAX(max_mem_used_ARRAY4,max_mem_used_ARRAY4_tmp)
      max_mem_used_ARRAY = MAX(max_mem_used_ARRAY,max_mem_used_tensor_tmp)
      max_mem_used_PNOSpaceInfo = MAX(max_mem_used_PNOSpaceInfo,max_mem_used_PNOSpaceInfo_tmp)
      max_mem_used_MP2DENS = MAX(max_mem_used_MP2DENS,max_mem_used_MP2DENS_tmp)
      max_mem_used_TRACEBACK = MAX(max_mem_used_TRACEBACK,max_mem_used_TRACEBACK_tmp)
      max_mem_used_MP2GRAD = MAX(max_mem_used_MP2GRAD,max_mem_used_MP2GRAD_tmp)
      max_mem_used_ODBATCH = MAX(max_mem_used_ODBATCH,max_mem_used_ODBATCH_tmp)
      max_mem_used_LSAOTENSOR = MAX(max_mem_used_LSAOTENSOR,max_mem_used_LSAOTENSOR_tmp)
      max_mem_used_SLSAOTENSOR = MAX(max_mem_used_SLSAOTENSOR,max_mem_used_SLSAOTENSOR_tmp)
      max_mem_used_GLOBALLSAOTENSOR = MAX(max_mem_used_GLOBALLSAOTENSOR,max_mem_used_GLOBALLSAOTENSOR_tmp)
      max_mem_used_ATOMTYPEITEM = MAX(max_mem_used_ATOMTYPEITEM,max_mem_used_ATOMTYPEITEM_tmp)
      max_mem_used_ATOMITEM = MAX(max_mem_used_ATOMITEM,max_mem_used_ATOMITEM_tmp)
      max_mem_used_LSMATRIX = MAX(max_mem_used_LSMATRIX,max_mem_used_LSMATRIX_tmp)

      max_mem_used_linkshell = MAX(max_mem_used_linkshell,max_mem_used_linkshell_tmp)
      max_mem_used_integralitem = MAX(max_mem_used_integralitem,max_mem_used_integralitem_tmp)
      max_mem_used_integrand = MAX(max_mem_used_integrand,max_mem_used_integrand_tmp)
      max_mem_used_intwork = MAX(max_mem_used_intwork,max_mem_used_intwork_tmp)
      max_mem_used_overlap = MAX(max_mem_used_overlap,max_mem_used_overlap_tmp)
      max_mem_used_overlapT = MAX(max_mem_used_overlapT,max_mem_used_overlapT_tmp)
      max_mem_used_ODitem = MAX(max_mem_used_ODitem,max_mem_used_ODitem_tmp)
      max_mem_used_FMM = MAX(max_mem_used_FMM,max_mem_used_FMM_tmp)
      max_mem_used_lstensor = MAX(max_mem_used_lstensor,max_mem_used_lstensor_tmp)
#ifdef MOD_UNRELEASED
      max_mem_used_lvec_data = MAX(max_mem_used_lvec_data,max_mem_used_lvec_data_tmp)
      max_mem_used_lattice_cell = MAX(max_mem_used_lattice_cell,max_mem_used_lattice_cell_tmp)
#endif

   end subroutine mem_TurnOffThread_Memory

   subroutine mem_TurnONThread_Memory()
      implicit none
      mem_InsideOMPsection = .TRUE.
      max_mem_used_global_tmp = 0
      max_mem_used_Heap_global_tmp = 0
      max_mem_used_type_matrix_tmp = 0
      max_mem_used_type_matrix_MPIFULL_tmp = 0
      max_mem_used_real_tmp = 0
      max_mem_used_mpi_tmp = 0
      max_mem_used_complex_tmp = 0
      max_mem_used_integer_tmp = 0
      max_mem_used_character_tmp = 0
      max_mem_used_logical_tmp = 0
      max_mem_used_AOBATCH_tmp = 0
      max_mem_used_DECORBITAL_tmp = 0
      max_mem_used_DECFRAG_tmp = 0
      max_mem_used_BATCHTOORB_tmp = 0
      max_mem_used_DECAOBATCHINFO_tmp = 0
      max_mem_used_MYPOINTER_tmp = 0
      max_mem_used_fragmentAOS_tmp = 0
      max_mem_used_ARRAY2_tmp = 0
      max_mem_used_ARRAY4_tmp = 0
      max_mem_used_tensor_tmp = 0
      max_mem_used_PNOSpaceInfo_tmp = 0
      max_mem_used_MP2DENS_tmp = 0
      max_mem_used_TRACEBACK_tmp = 0
      max_mem_used_MP2GRAD_tmp = 0
      max_mem_used_ODBATCH_tmp = 0
      max_mem_used_LSAOTENSOR_tmp = 0
      max_mem_used_SLSAOTENSOR_tmp = 0
      max_mem_used_GLOBALLSAOTENSOR_tmp = 0
      max_mem_used_ATOMTYPEITEM_tmp = 0
      max_mem_used_ATOMITEM_tmp = 0
      max_mem_used_LSMATRIX_tmp = 0
      max_mem_used_overlapT_tmp = 0

      max_mem_used_linkshell_tmp = 0
      max_mem_used_integralitem_tmp = 0
      max_mem_used_integrand_tmp = 0
      max_mem_used_intwork_tmp = 0
      max_mem_used_overlap_tmp = 0
      max_mem_used_ODitem_tmp = 0
      max_mem_used_FMM_tmp = 0
      max_mem_used_lstensor_tmp = 0
#ifdef MOD_UNRELEASED
      max_mem_used_lvec_data_tmp = 0
      max_mem_used_lattice_cell_tmp = 0
#endif
   end subroutine mem_TurnONThread_Memory

   subroutine init_globalmemvar()
      implicit none
      MemModParamPrintMemory = .FALSE.
      MemModParamPrintMemorylupri = 6
      PrintMemoryLowerLimit = 0
      call Set_PrintSCFmemory(.FALSE.)
      call set_sizes_of_types()
      mem_InsideOMPsection = .FALSE.
      mem_allocated_global = 0
      max_mem_used_global = 0
      mem_allocated_Heap_global = 0
      max_mem_used_Heap_global = 0
      mem_allocated_type_matrix = 0
      max_mem_used_type_matrix = 0
      mem_allocated_type_matrix_MPIFULL = 0
      max_mem_used_type_matrix_MPIFULL = 0
      mem_allocated_real = 0
      max_mem_used_real = 0
      mem_allocated_mpi = 0
      max_mem_used_mpi = 0
      mem_allocated_complex = 0
      max_mem_used_complex = 0
      mem_allocated_integer = 0
      max_mem_used_integer = 0
      mem_allocated_character = 0
      max_mem_used_character = 0
      mem_allocated_logical = 0
      max_mem_used_logical = 0
      mem_allocated_AOBATCH = 0
      max_mem_used_AOBATCH = 0
      mem_allocated_XCcalc = 0
      max_mem_used_XCcalc = 0 
      mem_allocated_DECORBITAL = 0
      max_mem_used_DECORBITAL = 0
      mem_allocated_DECFRAG = 0
      mem_allocated_BATCHTOORB = 0
      mem_allocated_DECAOBATCHINFO = 0
      mem_allocated_MYPOINTER = 0
      mem_allocated_fragmentAOS = 0
      mem_allocated_ARRAY2 = 0
      mem_allocated_ARRAY4 = 0
      mem_allocated_ARRAY = 0
      mem_allocated_PNOSpaceInfo = 0
      mem_allocated_MP2DENS = 0
      mem_allocated_TRACEBACK = 0
      mem_allocated_MP2GRAD = 0
      max_mem_used_DECFRAG = 0
      max_mem_used_BATCHTOORB = 0
      max_mem_used_DECAOBATCHINFO = 0
      max_mem_used_MYPOINTER = 0
      max_mem_used_fragmentAOS = 0
      max_mem_used_ARRAY2 = 0
      max_mem_used_ARRAY4 = 0
      max_mem_used_ARRAY = 0
      max_mem_used_PNOSpaceInfo = 0
      max_mem_used_MP2DENS = 0
      max_mem_used_TRACEBACK = 0
      max_mem_used_MP2GRAD = 0
      mem_allocated_ODBATCH = 0
      max_mem_used_ODBATCH = 0
      mem_allocated_LSAOTENSOR = 0
      max_mem_used_LSAOTENSOR = 0
      mem_allocated_SLSAOTENSOR = 0
      max_mem_used_SLSAOTENSOR = 0
      mem_allocated_GLOBALLSAOTENSOR = 0
      max_mem_used_GLOBALLSAOTENSOR = 0
      mem_allocated_ATOMTYPEITEM = 0
      max_mem_used_ATOMTYPEITEM = 0
      mem_allocated_ATOMITEM = 0
      max_mem_used_ATOMITEM = 0
      mem_allocated_LSMATRIX = 0
      max_mem_used_LSMATRIX = 0

      mem_allocated_linkshell = 0
      max_mem_used_linkshell = 0
      mem_allocated_integralitem = 0
      max_mem_used_integralitem = 0
      mem_allocated_integrand = 0
      max_mem_used_integrand = 0
      mem_allocated_intwork = 0
      max_mem_used_intwork = 0
      mem_allocated_overlap = 0
      max_mem_used_overlap = 0
      mem_allocated_overlapT = 0
      max_mem_used_overlapT = 0
      mem_allocated_ODitem = 0
      max_mem_used_ODitem = 0
      mem_allocated_FMM = 0
      max_mem_used_FMM = 0
      mem_allocated_lstensor = 0
      max_mem_used_lstensor = 0
#ifdef MOD_UNRELEASED
      mem_allocated_lvec_data = 0
      max_mem_used_lvec_data = 0
      mem_allocated_lattice_cell = 0
      max_mem_used_lattice_cell = 0
#endif
      call init_threadmemvar()
   end subroutine init_globalmemvar

   subroutine init_threadmemvar()
      implicit none

      mem_tp_allocated_global = 0
      max_mem_tp_used_global = 0
      mem_tp_allocated_Heap_global = 0
      max_mem_tp_used_Heap_global = 0
      mem_tp_allocated_type_matrix = 0
      max_mem_tp_used_type_matrix = 0
      mem_tp_allocated_type_matrix_MPIFULL = 0
      max_mem_tp_used_type_matrix_MPIFULL = 0
      mem_tp_allocated_real = 0
      max_mem_tp_used_real = 0
      mem_tp_allocated_mpi = 0
      max_mem_tp_used_mpi = 0
      mem_tp_allocated_complex = 0
      max_mem_tp_used_complex = 0
      mem_tp_allocated_integer = 0
      max_mem_tp_used_integer = 0
      mem_tp_allocated_character = 0
      max_mem_tp_used_character = 0
      mem_tp_allocated_logical = 0
      max_mem_tp_used_logical = 0
      mem_tp_allocated_AOBATCH = 0
      max_mem_tp_used_AOBATCH = 0
      mem_tp_allocated_DECORBITAL = 0
      max_mem_tp_used_DECORBITAL = 0
      mem_tp_allocated_DECFRAG = 0
      mem_tp_allocated_BATCHTOORB = 0
      mem_tp_allocated_DECAOBATCHINFO = 0
      mem_tp_allocated_MYPOINTER = 0
      mem_tp_allocated_fragmentAOS = 0
      mem_tp_allocated_ARRAY2 = 0
      mem_tp_allocated_ARRAY4 = 0
      mem_tp_allocated_ARRAY = 0
      mem_tp_allocated_PNOSpaceInfo = 0
      mem_tp_allocated_MP2DENS = 0
      mem_tp_allocated_TRACEBACK = 0
      mem_tp_allocated_MP2GRAD = 0
      max_mem_tp_used_DECFRAG = 0
      max_mem_tp_used_BATCHTOORB = 0
      max_mem_tp_used_DECAOBATCHINFO = 0
      max_mem_tp_used_MYPOINTER = 0
      max_mem_tp_used_fragmentAOS = 0
      max_mem_tp_used_ARRAY2 = 0
      max_mem_tp_used_ARRAY4 = 0
      max_mem_tp_used_ARRAY = 0
      max_mem_tp_used_PNOSpaceInfo = 0
      max_mem_tp_used_MP2DENS = 0
      max_mem_tp_used_TRACEBACK = 0
      max_mem_tp_used_MP2GRAD = 0
      mem_tp_allocated_ODBATCH = 0
      max_mem_tp_used_ODBATCH = 0
      mem_tp_allocated_LSAOTENSOR = 0
      max_mem_tp_used_LSAOTENSOR = 0
      mem_tp_allocated_SLSAOTENSOR = 0
      max_mem_tp_used_SLSAOTENSOR = 0
      mem_tp_allocated_GLOBALLSAOTENSOR = 0
      max_mem_tp_used_GLOBALLSAOTENSOR = 0
      mem_tp_allocated_ATOMTYPEITEM = 0
      max_mem_tp_used_ATOMTYPEITEM = 0
      mem_tp_allocated_ATOMITEM = 0
      max_mem_tp_used_ATOMITEM = 0
      mem_tp_allocated_LSMATRIX = 0
      max_mem_tp_used_LSMATRIX = 0

      mem_tp_allocated_linkshell = 0
      max_mem_tp_used_linkshell = 0
      mem_tp_allocated_integralitem = 0
      max_mem_tp_used_integralitem = 0
      mem_tp_allocated_integrand = 0
      max_mem_tp_used_integrand = 0
      mem_tp_allocated_intwork = 0
      max_mem_tp_used_intwork = 0
      mem_tp_allocated_overlap = 0
      max_mem_tp_used_overlap = 0
      mem_tp_allocated_overlapT = 0
      max_mem_tp_used_overlapT = 0
      mem_tp_allocated_ODitem = 0
      max_mem_tp_used_ODitem = 0
      mem_tp_allocated_FMM = 0
      max_mem_tp_used_FMM = 0
      mem_tp_allocated_lstensor = 0
      max_mem_tp_used_lstensor = 0
#ifdef MOD_UNRELEASED
      mem_tp_allocated_lvec_data = 0
      max_mem_tp_used_lvec_data = 0
      mem_tp_allocated_lattice_cell = 0
      max_mem_tp_used_lattice_cell = 0
#endif
   end subroutine init_threadmemvar

   subroutine collect_thread_memory()
      implicit none
      !$OMP CRITICAL
      mem_allocated_global = mem_allocated_global+mem_tp_allocated_global
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in collect_thread_memory'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in collect_thread_memory'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in collect_thread_memory')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global_tmp = max_mem_used_global_tmp+max_mem_tp_used_global
      max_mem_used_Heap_global_tmp = max_mem_used_Heap_global_tmp+max_mem_tp_used_Heap_global
      mem_allocated_type_matrix = mem_allocated_type_matrix+mem_tp_allocated_type_matrix
      max_mem_used_type_matrix_tmp = max_mem_used_type_matrix_tmp+max_mem_tp_used_type_matrix
      mem_allocated_type_matrix_MPIFULL = mem_allocated_type_matrix_MPIFULL+mem_tp_allocated_type_matrix_MPIFULL
      max_mem_used_type_matrix_MPIFULL_tmp = max_mem_used_type_matrix_MPIFULL_tmp+max_mem_tp_used_type_matrix_MPIFULL
      mem_allocated_real = mem_allocated_real+mem_tp_allocated_real
      max_mem_used_real_tmp = max_mem_used_real_tmp+max_mem_tp_used_real
      mem_allocated_mpi = mem_allocated_mpi+mem_tp_allocated_mpi
      max_mem_used_mpi_tmp = max_mem_used_mpi_tmp+max_mem_tp_used_mpi
      mem_allocated_complex = mem_allocated_complex+mem_tp_allocated_complex
      max_mem_used_complex_tmp = max_mem_used_complex_tmp+max_mem_tp_used_complex
      mem_allocated_integer = mem_allocated_integer+mem_tp_allocated_integer
      max_mem_used_integer_tmp = max_mem_used_integer_tmp+max_mem_tp_used_integer
      mem_allocated_character = mem_allocated_character+mem_tp_allocated_character
      max_mem_used_character_tmp = max_mem_used_character_tmp+max_mem_tp_used_character
      mem_allocated_logical = mem_allocated_logical+mem_tp_allocated_logical
      max_mem_used_logical_tmp = max_mem_used_logical_tmp+max_mem_tp_used_logical
      mem_allocated_AOBATCH = mem_allocated_AOBATCH+mem_tp_allocated_AOBATCH
      max_mem_used_AOBATCH_tmp = max_mem_used_AOBATCH_tmp+max_mem_tp_used_AOBATCH
      mem_allocated_DECORBITAL = mem_allocated_DECORBITAL+mem_tp_allocated_DECORBITAL
      max_mem_used_DECORBITAL_tmp = max_mem_used_DECORBITAL_tmp+max_mem_tp_used_DECORBITAL
      mem_allocated_DECFRAG = mem_allocated_DECFRAG+mem_tp_allocated_DECFRAG
      mem_allocated_BATCHTOORB = mem_allocated_BATCHTOORB+mem_tp_allocated_BATCHTOORB
      mem_allocated_DECAOBATCHINFO = mem_allocated_DECAOBATCHINFO+mem_tp_allocated_DECAOBATCHINFO
      mem_allocated_MYPOINTER = mem_allocated_MYPOINTER+mem_tp_allocated_MYPOINTER
      mem_allocated_fragmentAOS = mem_allocated_fragmentAOS+mem_tp_allocated_fragmentAOS
      mem_allocated_ARRAY2 = mem_allocated_ARRAY2+mem_tp_allocated_ARRAY2
      mem_allocated_ARRAY4 = mem_allocated_ARRAY4+mem_tp_allocated_ARRAY4
      mem_allocated_ARRAY = mem_allocated_ARRAY+mem_tp_allocated_ARRAY
      mem_allocated_PNOSpaceInfo = mem_allocated_PNOSpaceInfo+mem_tp_allocated_PNOSpaceInfo
      mem_allocated_MP2DENS = mem_allocated_MP2DENS+mem_tp_allocated_MP2DENS
      mem_allocated_TRACEBACK = mem_allocated_TRACEBACK+mem_tp_allocated_TRACEBACK
      mem_allocated_MP2GRAD = mem_allocated_MP2GRAD+mem_tp_allocated_MP2GRAD
      max_mem_used_DECFRAG_tmp = max_mem_used_DECFRAG_tmp+max_mem_tp_used_DECFRAG
      max_mem_used_BATCHTOORB_tmp = max_mem_used_BATCHTOORB_tmp+max_mem_tp_used_BATCHTOORB
      max_mem_used_DECAOBATCHINFO_tmp = max_mem_used_DECAOBATCHINFO_tmp+max_mem_tp_used_DECAOBATCHINFO
      max_mem_used_MYPOINTER_tmp = max_mem_used_MYPOINTER_tmp+max_mem_tp_used_MYPOINTER
      max_mem_used_fragmentAOS_tmp = max_mem_used_fragmentAOS_tmp+max_mem_tp_used_fragmentAOS
      max_mem_used_ARRAY2_tmp = max_mem_used_ARRAY2_tmp+max_mem_tp_used_ARRAY2
      max_mem_used_ARRAY4_tmp = max_mem_used_ARRAY4_tmp+max_mem_tp_used_ARRAY4
      max_mem_used_tensor_tmp = max_mem_used_tensor_tmp+max_mem_tp_used_ARRAY
      max_mem_used_PNOSpaceInfo_tmp = max_mem_used_PNOSpaceInfo_tmp+max_mem_tp_used_PNOSpaceInfo
      max_mem_used_MP2DENS_tmp = max_mem_used_MP2DENS_tmp+max_mem_tp_used_MP2DENS
      max_mem_used_TRACEBACK_tmp = max_mem_used_TRACEBACK_tmp+max_mem_tp_used_TRACEBACK
      max_mem_used_MP2GRAD_tmp = max_mem_used_MP2GRAD_tmp+max_mem_tp_used_MP2GRAD
      mem_allocated_ODBATCH = mem_allocated_ODBATCH+mem_tp_allocated_ODBATCH
      max_mem_used_ODBATCH_tmp = max_mem_used_ODBATCH_tmp+max_mem_tp_used_ODBATCH
      mem_allocated_LSAOTENSOR = mem_allocated_LSAOTENSOR+mem_tp_allocated_LSAOTENSOR
      max_mem_used_LSAOTENSOR_tmp = max_mem_used_LSAOTENSOR_tmp+max_mem_tp_used_LSAOTENSOR
      mem_allocated_SLSAOTENSOR = mem_allocated_SLSAOTENSOR+mem_tp_allocated_SLSAOTENSOR
      max_mem_used_SLSAOTENSOR_tmp = max_mem_used_SLSAOTENSOR_tmp+max_mem_tp_used_SLSAOTENSOR
      mem_allocated_GLOBALLSAOTENSOR = mem_allocated_GLOBALLSAOTENSOR+mem_tp_allocated_GLOBALLSAOTENSOR
      max_mem_used_GLOBALLSAOTENSOR_tmp = max_mem_used_GLOBALLSAOTENSOR_tmp+max_mem_tp_used_GLOBALLSAOTENSOR
      mem_allocated_ATOMTYPEITEM = mem_allocated_ATOMTYPEITEM+mem_tp_allocated_ATOMTYPEITEM
      max_mem_used_ATOMTYPEITEM_tmp = max_mem_used_ATOMTYPEITEM_tmp+max_mem_tp_used_ATOMTYPEITEM
      mem_allocated_ATOMITEM = mem_allocated_ATOMITEM+mem_tp_allocated_ATOMITEM
      max_mem_used_ATOMITEM_tmp = max_mem_used_ATOMITEM_tmp+max_mem_tp_used_ATOMITEM
      mem_allocated_LSMATRIX = mem_allocated_LSMATRIX+mem_tp_allocated_LSMATRIX
      max_mem_used_LSMATRIX_tmp = max_mem_used_LSMATRIX_tmp+max_mem_tp_used_LSMATRIX

      mem_allocated_linkshell = mem_allocated_linkshell+mem_tp_allocated_linkshell
      max_mem_used_linkshell_tmp = max_mem_used_linkshell_tmp+max_mem_tp_used_linkshell
      mem_allocated_integralitem = mem_allocated_integralitem+mem_tp_allocated_integralitem
      max_mem_used_integralitem_tmp = max_mem_used_integralitem_tmp+max_mem_tp_used_integralitem
      mem_allocated_integrand = mem_allocated_integrand+mem_tp_allocated_integrand
      max_mem_used_integrand_tmp = max_mem_used_integrand_tmp+max_mem_tp_used_integrand
      mem_allocated_intwork = mem_allocated_intwork+mem_tp_allocated_intwork
      max_mem_used_intwork_tmp = max_mem_used_intwork_tmp+max_mem_tp_used_intwork
      mem_allocated_overlap = mem_allocated_overlap+mem_tp_allocated_overlap
      max_mem_used_overlap_tmp = max_mem_used_overlap_tmp+max_mem_tp_used_overlap
      mem_allocated_overlapT = mem_allocated_overlapT+mem_tp_allocated_overlapT
      max_mem_used_overlapT_tmp = max_mem_used_overlapT_tmp+max_mem_tp_used_overlapT
      mem_allocated_ODitem = mem_allocated_ODitem+mem_tp_allocated_ODitem
      max_mem_used_ODitem_tmp = max_mem_used_ODitem_tmp+max_mem_tp_used_ODitem
      mem_allocated_FMM = mem_allocated_FMM+mem_tp_allocated_FMM
      max_mem_used_FMM_tmp = max_mem_used_FMM_tmp+max_mem_tp_used_FMM
      mem_allocated_lstensor = mem_allocated_lstensor+mem_tp_allocated_lstensor
      max_mem_used_lstensor_tmp = max_mem_used_lstensor_tmp+max_mem_tp_used_lstensor
#ifdef MOD_UNRELEASED
      mem_allocated_lvec_data = mem_allocated_lvec_data+mem_tp_allocated_lvec_data
      max_mem_used_lvec_data_tmp = max_mem_used_lvec_data_tmp+max_mem_tp_used_lvec_data
      mem_allocated_lattice_cell = mem_allocated_lattice_cell+mem_tp_allocated_lattice_cell
      max_mem_used_lattice_cell_tmp = max_mem_used_lattice_cell_tmp+max_mem_tp_used_lattice_cell
#endif
      !$OMP END CRITICAL
   end subroutine collect_thread_memory

   !> \brief Print global mem
   !> \author T. Kjaergaard
   !> \date 2014
   subroutine stats_globalmem(lupri)
      implicit none
      !> Logical unit number for output file.
      integer,intent(in) :: lupri
      call print_mem_alloc(lupri,mem_allocated_global,'TOTAL')
      call print_maxmem(lupri,max_mem_used_global,'TOTAL')
    end subroutine stats_globalmem

   !> \brief Tool to detect Heap Memory usage in part of code
   !> \author T. Kjaergaard
   !> \date 2014
   subroutine DetectHeapMemoryInit()
      implicit none
      mem_allocated_Heap_global = 0
      max_mem_used_Heap_global = 0
    end subroutine DetectHeapMemoryInit

   !> \brief Tool to detect Heap Memory usage in part of code
   !> \author T. Kjaergaard
   !> \date 2014
   subroutine DetectHeapMemory(mem_allocated,Max_mem_used,String,lupri)
      implicit none
      integer(kind=long),intent(inout) :: mem_allocated,Max_mem_used
      Character*(*),intent(in) ::  STRING
      integer,intent(in) :: lupri
      !
      character(len=100) :: trimstring
      trimstring = adjustl(string)
      mem_allocated = mem_allocated_Heap_global
      Max_mem_used = max_mem_used_Heap_global
!      if (max_mem_used_Heap_global < 100) then !Divide by 1 to typecast real
!         WRITE(LUPRI,'(A,A,A,f9.3,A)')'Max Memory Usage in ',TRIM(trimstring),' = ', &
!              & max_mem_used_Heap_global/1.0E0_realk,' Byte'
!      else if (max_mem_used_Heap_global < 1000000) then
!         WRITE(LUPRI,'(A,A,A,f9.3,A)')'Max Memory Usage in ',TRIM(trimstring),' = ', &
!              & max_mem_used_Heap_global/1000.0E0_realk,' kB'
!      else if (max_mem_used_Heap_global < 1000000000) then
!         WRITE(LUPRI,'(A,A,A,f9.3,A)')'Max Memory Usage in ',TRIM(trimstring),' = ', &
!              & max_mem_used_Heap_global/(1000.0E0_realk*1000.0E0_realk),' MB'
!      else
         WRITE(LUPRI,'(A,A,A,f9.3,A)')'Max Memory Usage in ',TRIM(trimstring),' = ', &
              & max_mem_used_Heap_global/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk),' GB'
!      endif
    end subroutine DetectHeapMemory

   !> \brief Print current and max. amount of memory allocated for different data types.
   !> \author S. Host
   !> \date 2009
   subroutine stats_mem(lupri)
      implicit none
      !> Logical unit number for output file.
      integer,intent(in) :: lupri
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("                  Memory statistics          ")')
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("  Allocated memory (TOTAL):             ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_global
      WRITE(LUPRI,'("  Allocated memory (type(matrix)):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_type_matrix
      WRITE(LUPRI,'("  Allocated memory (real):              ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_real
      WRITE(LUPRI,'("  Allocated memory (MPI):               ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_mpi
      WRITE(LUPRI,'("  Allocated memory (complex):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_complex
      WRITE(LUPRI,'("  Allocated memory (integer):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integer
      WRITE(LUPRI,'("  Allocated memory (logical):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_logical
      WRITE(LUPRI,'("  Allocated memory (character):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_character
      WRITE(LUPRI,'("  Allocated memory (AOBATCH):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_AOBATCH
      WRITE(LUPRI,'("  Allocated memory (ODBATCH):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ODBATCH
      WRITE(LUPRI,'("  Allocated memory (LSAOTENSOR):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_LSAOTENSOR
      WRITE(LUPRI,'("  Allocated memory (SLSAOTENSOR):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_SLSAOTENSOR
      WRITE(LUPRI,'("  Allocated memory (GLOBALLSAOTENSOR):  ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_GLOBALLSAOTENSOR
      WRITE(LUPRI,'("  Allocated memory (ATOMTYPEITEM):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ATOMTYPEITEM
      WRITE(LUPRI,'("  Allocated memory (ATOMITEM):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ATOMITEM
      WRITE(LUPRI,'("  Allocated memory (LSMATRIX):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_LSMATRIX
      WRITE(LUPRI,'("  Allocated memory (DECORBITAL):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_DECORBITAL
      WRITE(LUPRI,'("  Allocated memory (DECFRAG):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_DECFRAG
      WRITE(LUPRI,'("  Allocated memory (overlapType):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_overlapT
      WRITE(LUPRI,'("  Allocated memory (BATCHTOORB):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_BATCHTOORB
      WRITE(LUPRI,'("  Allocated memory (DECAOBATCHINFO):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_DECAOBATCHINFO
      WRITE(LUPRI,'("  Allocated memory (MYPOINTER):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_MYPOINTER
      WRITE(LUPRI,'("  Allocated memory (fragmentAOS):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_fragmentAOS
      WRITE(LUPRI,'("  Allocated memory (ARRAY2):            ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ARRAY2
      WRITE(LUPRI,'("  Allocated memory (ARRAY4):            ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ARRAY4
      WRITE(LUPRI,'("  Allocated memory (ARRAY):             ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ARRAY
      WRITE(LUPRI,'("  Allocated memory (PNOSpaceInfo):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_PNOSpaceInfo
      WRITE(LUPRI,'("  Allocated memory (MP2DENS):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_MP2DENS
      WRITE(LUPRI,'("  Allocated memory (TRACEBACK):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_TRACEBACK
      WRITE(LUPRI,'("  Allocated memory (MP2GRAD):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_MP2GRAD
#ifdef MOD_UNRELEASED
      WRITE(LUPRI,'("  Allocated memory (type(lvec_data)):   ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_lvec_data
      WRITE(LUPRI,'("  Allocated memory (type(lattice_cell)):",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_lattice_cell
#endif
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("                  Additional Memory information          ")')
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("  Allocated memory (linkshell):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_linkshell
      WRITE(LUPRI,'("  Allocated memory (integrand):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integrand
      WRITE(LUPRI,'("  Allocated memory (integralitem):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integralitem
      WRITE(LUPRI,'("  Allocated memory (intwork):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_intwork
      WRITE(LUPRI,'("  Allocated memory (overlap):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_overlap
      WRITE(LUPRI,'("  Allocated memory (ODitem):            ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ODitem
      WRITE(LUPRI,'("  Allocated memory (lstensor):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_lstensor
      WRITE(LUPRI,'("  Allocated memory (FMM   ):            ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_FMM
      WRITE(LUPRI,'("  Allocated memory (XC    ):            ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_XCcalc

      call print_maxmem(lupri,max_mem_used_global,'TOTAL')
#ifdef VAR_SCALAPACK
      call print_maxmem(lupri,max_mem_used_type_matrix,'Type(matrix) on Master node')
      call print_maxmem(lupri,max_mem_used_type_matrix_MPIFULL,'Type(matrix) on all nodes')
#else
      call print_maxmem(lupri,max_mem_used_type_matrix,'type(matrix)')
#endif
      CALL print_maxmem(lupri,max_mem_used_real,'real(realk)')
      CALL print_maxmem(lupri,max_mem_used_mpi,'MPI')
      CALL print_maxmem(lupri,max_mem_used_complex,'complex(complexk)')
      CALL print_maxmem(lupri,max_mem_used_integer,'integer')
      CALL print_maxmem(lupri,max_mem_used_logical,'logical')
      CALL print_maxmem(lupri,max_mem_used_character,'character')
      CALL print_maxmem(lupri,max_mem_used_AOBATCH,'AOBATCH')
      CALL print_maxmem(lupri,max_mem_used_DECORBITAL,'DECORBITAL')
      CALL print_maxmem(lupri,max_mem_used_DECFRAG,'DECFRAG')
      CALL print_maxmem(lupri,max_mem_used_BATCHTOORB,'BATCHTOORB')
      CALL print_maxmem(lupri,max_mem_used_DECAOBATCHINFO,'DECAOBATCHINFO')
      CALL print_maxmem(lupri,max_mem_used_MYPOINTER,'MYPOINTER')
      CALL print_maxmem(lupri,max_mem_used_fragmentAOS,'fragmentAOS')
      CALL print_maxmem(lupri,max_mem_used_ARRAY2,'ARRAY2')
      CALL print_maxmem(lupri,max_mem_used_ARRAY4,'ARRAY4')
      CALL print_maxmem(lupri,max_mem_used_ARRAY,'ARRAY')
      CALL print_maxmem(lupri,max_mem_used_PNOSpaceInfo,'PNOSpaceInfo')
      CALL print_maxmem(lupri,max_mem_used_MP2DENS,'MP2DENS') 
      CALL print_maxmem(lupri,max_mem_used_TRACEBACK,'TRACEBACK')
      CALL print_maxmem(lupri,max_mem_used_MP2GRAD,'MP2GRAD')
      CALL print_maxmem(lupri,max_mem_used_ODBATCH,'ODBATCH')
      CALL print_maxmem(lupri,max_mem_used_LSAOTENSOR,'LSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_used_SLSAOTENSOR,'SLSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_used_GLOBALLSAOTENSOR,'GLOBALLSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_used_ATOMTYPEITEM,'ATOMTYPEITEM')
      CALL print_maxmem(lupri,max_mem_used_ATOMITEM,'ATOMITEM')
      CALL print_maxmem(lupri,max_mem_used_LSMATRIX,'LSMATRIX')
      CALL print_maxmem(lupri,max_mem_used_overlapT,'OverlapT')


      CALL print_maxmem(lupri,max_mem_used_linkshell,'linkshell')
      CALL print_maxmem(lupri,max_mem_used_integrand,'integrand')
      CALL print_maxmem(lupri,max_mem_used_integralitem,'integralitem')
      CALL print_maxmem(lupri,max_mem_used_intwork,'IntWork')
      CALL print_maxmem(lupri,max_mem_used_overlap,'Overlap')
      CALL print_maxmem(lupri,max_mem_used_ODitem,'ODitem')
      CALL print_maxmem(lupri,max_mem_used_lstensor,'LStensor')
      CALL print_maxmem(lupri,max_mem_used_FMM,'FMM')
      CALL print_maxmem(lupri,max_mem_used_XCcalc,'XC')
#ifdef MOD_UNRELEASED
      CALL print_maxmem(lupri,max_mem_used_lvec_data,'Lvec_data')
      CALL print_maxmem(lupri,max_mem_used_lattice_cell,'Lattice_cell')
#endif

      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,*)
   end subroutine stats_mem

   !> \brief Print current and max. amount of memory allocated for different data types.
   !> \brief In order to test MPI and NON MPI seperately we have seperate printout
   !> \author S. Host
   !> \date 2009
   subroutine stats_mpi_mem(lupri)
      implicit none
      !> Logical unit number for output file.
      integer,intent(in) :: lupri
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("                  Memory statistics          ")')
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("  Allocated MPI memory (TOTAL):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_global
      WRITE(LUPRI,'("  Allocated MPI memory (type(matrix)):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_type_matrix
      WRITE(LUPRI,'("  Allocated MPI memory (real):            ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_real
      WRITE(LUPRI,'("  Allocated MPI memory (MPI):             ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_mpi
      WRITE(LUPRI,'("  Allocated MPI memory (complex):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_complex
      WRITE(LUPRI,'("  Allocated MPI memory (integer):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integer
      WRITE(LUPRI,'("  Allocated MPI memory (logical):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_logical
      WRITE(LUPRI,'("  Allocated MPI memory (character):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_character
      WRITE(LUPRI,'("  Allocated MPI memory (AOBATCH):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_AOBATCH
      WRITE(LUPRI,'("  Allocated MPI memory (ODBATCH):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ODBATCH
      WRITE(LUPRI,'("  Allocated MPI memory (LSAOTENSOR):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_LSAOTENSOR
      WRITE(LUPRI,'("  Allocated MPI memory (SLSAOTENSOR):     ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_SLSAOTENSOR
      WRITE(LUPRI,'("  Allocated MPI memory (GLOBALLSAOTENSOR):",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_GLOBALLSAOTENSOR
      WRITE(LUPRI,'("  Allocated MPI memory (ATOMTYPEITEM):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ATOMTYPEITEM
      WRITE(LUPRI,'("  Allocated MPI memory (ATOMITEM):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ATOMITEM
      WRITE(LUPRI,'("  Allocated MPI memory (LSMATRIX):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_LSMATRIX
      WRITE(LUPRI,'("  Allocated MPI memory (DECORBITAL):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_DECORBITAL
      WRITE(LUPRI,'("  Allocated MPI memory (DECFRAG):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_DECFRAG
      WRITE(LUPRI,'("  Allocated MPI memory (overlapType):     ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_overlapT
      WRITE(LUPRI,'("  Allocated MPI memory (BATCHTOORB):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_BATCHTOORB
      WRITE(LUPRI,'("  Allocated MPI memory (DECAOBATCHINFO):  ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_DECAOBATCHINFO
      WRITE(LUPRI,'("  Allocated MPI memory (MYPOINTER):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_MYPOINTER
      WRITE(LUPRI,'("  Allocated MPI memory (fragmentAOS):     ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_fragmentAOS
      WRITE(LUPRI,'("  Allocated MPI memory (ARRAY2):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ARRAY2
      WRITE(LUPRI,'("  Allocated MPI memory (ARRAY4):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ARRAY4
      WRITE(LUPRI,'("  Allocated MPI memory (ARRAY):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ARRAY
      WRITE(LUPRI,'("  Allocated MPI memory (PNOSpaceInfo):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_PNOSpaceInfo
      WRITE(LUPRI,'("  Allocated MPI memory (MP2DENS):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_MP2DENS
      WRITE(LUPRI,'("  Allocated MPI memory (TRACEBACK):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_TRACEBACK
      WRITE(LUPRI,'("  Allocated MPI memory (MP2GRAD):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_MP2GRAD
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("                  Additional Memory information          ")')
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("  Allocated MPI memory (linkshell):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_linkshell
      WRITE(LUPRI,'("  Allocated MPI memory (integrand):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integrand
      WRITE(LUPRI,'("  Allocated MPI memory (integralitem):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_integralitem
      WRITE(LUPRI,'("  Allocated MPI memory (IntWork):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_IntWork
      WRITE(LUPRI,'("  Allocated MPI memory (overlap):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_overlap
      WRITE(LUPRI,'("  Allocated MPI memory (ODitem):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_ODitem
      WRITE(LUPRI,'("  Allocated MPI memory (lstensor):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_lstensor
      WRITE(LUPRI,'("  Allocated MPI memory (FMM   ):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_FMM
      WRITE(LUPRI,'("  Allocated MPI memory (XC    ):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_XCcalc
#ifdef MOD_UNRELEASED
      WRITE(LUPRI,'("  Allocated MPI memory (lvec_data):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_lvec_data
      WRITE(LUPRI,'("  Allocated MPI memory (lattice_cell):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_allocated_lattice_cell
#endif

      call print_maxmem(lupri,max_mem_used_global,'TOTAL')
#ifdef VAR_SCALAPACK
      call print_maxmem(lupri,max_mem_used_type_matrix,'Type(matrix) on Master node')
      call print_maxmem(lupri,max_mem_used_type_matrix_MPIFULL,'Type(matrix) on all nodes')
#else
      call print_maxmem(lupri,max_mem_used_type_matrix,'type(matrix)')
#endif
      CALL print_maxmem(lupri,max_mem_used_real,'real(realk)')
      CALL print_maxmem(lupri,max_mem_used_mpi,'MPI')
      CALL print_maxmem(lupri,max_mem_used_complex,'complex(complexk)')
      CALL print_maxmem(lupri,max_mem_used_integer,'integer')
      CALL print_maxmem(lupri,max_mem_used_logical,'logical')
      CALL print_maxmem(lupri,max_mem_used_character,'character')
      CALL print_maxmem(lupri,max_mem_used_AOBATCH,'AOBATCH')
      CALL print_maxmem(lupri,max_mem_used_DECORBITAL,'DECORBITAL')
      CALL print_maxmem(lupri,max_mem_used_DECFRAG,'DECFRAG')
      CALL print_maxmem(lupri,max_mem_used_BATCHTOORB,'BATCHTOORB')
      CALL print_maxmem(lupri,max_mem_used_DECAOBATCHINFO,'DECAOBATCHINFO')
      CALL print_maxmem(lupri,max_mem_used_MYPOINTER,'MYPOINTER')
      CALL print_maxmem(lupri,max_mem_used_fragmentAOS,'fragmentAOS')
      CALL print_maxmem(lupri,max_mem_used_ARRAY2,'ARRAY2')
      CALL print_maxmem(lupri,max_mem_used_ARRAY4,'ARRAY4')
      CALL print_maxmem(lupri,max_mem_used_ARRAY,'ARRAY')
      CALL print_maxmem(lupri,max_mem_used_PNOSpaceInfo,'PNOSpaceInfo')
      CALL print_maxmem(lupri,max_mem_used_MP2DENS,'MP2DENS')
      CALL print_maxmem(lupri,max_mem_used_TRACEBACK,'TRACEBACK')
      CALL print_maxmem(lupri,max_mem_used_MP2GRAD,'MP2GRAD')
      CALL print_maxmem(lupri,max_mem_used_ODBATCH,'ODBATCH')
      CALL print_maxmem(lupri,max_mem_used_LSAOTENSOR,'LSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_used_SLSAOTENSOR,'SLSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_used_GLOBALLSAOTENSOR,'GLOBALLSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_used_ATOMTYPEITEM,'ATOMTYPEITEM')
      CALL print_maxmem(lupri,max_mem_used_ATOMITEM,'ATOMITEM')
      CALL print_maxmem(lupri,max_mem_used_LSMATRIX,'LSMATRIX')
      CALL print_maxmem(lupri,max_mem_used_overlapT,'OverlapType')


      CALL print_maxmem(lupri,max_mem_used_linkshell,'linkshell')
      CALL print_maxmem(lupri,max_mem_used_integrand,'integrand')
      CALL print_maxmem(lupri,max_mem_used_integralitem,'integralitem')
      CALL print_maxmem(lupri,max_mem_used_IntWork,'IntWork')
      CALL print_maxmem(lupri,max_mem_used_overlap,'Overlap')
      CALL print_maxmem(lupri,max_mem_used_ODitem,'ODitem')
      CALL print_maxmem(lupri,max_mem_used_lstensor,'LStensor')
      CALL print_maxmem(lupri,max_mem_used_FMM,'FMM')
      CALL print_maxmem(lupri,max_mem_used_XCcalc,'XC ')
#ifdef MOD_UNRELEASED
      CALL print_maxmem(lupri,max_mem_used_lvec_data,'Lvec_data')
      CALL print_maxmem(lupri,max_mem_used_lattice_cell,'Lattice_cell')
#endif

      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,*)
   end subroutine stats_mpi_mem

   !> \brief Print current and max. amount of memory allocated for different data types.
   !> \author S. Host
   !> \date 2009
   subroutine stats_mem_tp(lupri)
      implicit none
      !> Logical unit number for output file.
      integer,intent(in) :: lupri
      !$OMP CRITICAL
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("  Tread Private             Memory statistics          ")')
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("  Allocated memory (TOTAL):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_global
      WRITE(LUPRI,'("  Allocated memory (type(matrix)):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_type_matrix
      WRITE(LUPRI,'("  Allocated memory (real):            ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_real
      WRITE(LUPRI,'("  Allocated memory (MPI):             ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_mpi
      WRITE(LUPRI,'("  Allocated memory (complex):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_complex
      WRITE(LUPRI,'("  Allocated memory (integer):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_integer
      WRITE(LUPRI,'("  Allocated memory (logical):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_logical
      WRITE(LUPRI,'("  Allocated memory (character):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_character
      WRITE(LUPRI,'("  Allocated memory (AOBATCH):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_AOBATCH
      WRITE(LUPRI,'("  Allocated memory (ODBATCH):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_ODBATCH
      WRITE(LUPRI,'("  Allocated memory (LSAOTENSOR):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_LSAOTENSOR
      WRITE(LUPRI,'("  Allocated memory (SLSAOTENSOR):     ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_SLSAOTENSOR
      WRITE(LUPRI,'("  Allocated memory (GLOBALLSAOTENSOR):",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_GLOBALLSAOTENSOR
      WRITE(LUPRI,'("  Allocated memory (ATOMTYPEITEM):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_ATOMTYPEITEM
      WRITE(LUPRI,'("  Allocated memory (ATOMITEM):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_ATOMITEM
      WRITE(LUPRI,'("  Allocated memory (LSMATRIX):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_LSMATRIX
      WRITE(LUPRI,'("  Allocated memory (DECORBITAL):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_DECORBITAL
      WRITE(LUPRI,'("  Allocated memory (DECFRAG):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_DECFRAG
      WRITE(LUPRI,'("  Allocated memory (overlapType):     ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_overlapT
      WRITE(LUPRI,'("  Allocated memory (BATCHTOORB):      ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_BATCHTOORB
      WRITE(LUPRI,'("  Allocated memory (DECAOBATCHINFO):  ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_DECAOBATCHINFO
      WRITE(LUPRI,'("  Allocated memory (MYPOINTER):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_MYPOINTER
      WRITE(LUPRI,'("  Allocated memory (fragmentAOS):     ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_fragmentAOS
      WRITE(LUPRI,'("  Allocated memory (ARRAY2):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_ARRAY2
      WRITE(LUPRI,'("  Allocated memory (ARRAY4):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_ARRAY4
      WRITE(LUPRI,'("  Allocated memory (ARRAY):           ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_ARRAY
      WRITE(LUPRI,'("  Allocated memory (PNOSpaceInfo):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_PNOSpaceInfo
      WRITE(LUPRI,'("  Allocated memory (MP2DENS):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_MP2DENS
      WRITE(LUPRI,'("  Allocated memory (TRACEBACK):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_TRACEBACK
      WRITE(LUPRI,'("  Allocated memory (MP2GRAD):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_MP2GRAD
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("                  Additional Memory information          ")')
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("  Allocated memory (linkshell):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_linkshell
      WRITE(LUPRI,'("  Allocated memory (integrand):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_integrand
      WRITE(LUPRI,'("  Allocated memory (integralitem):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_integralitem
      WRITE(LUPRI,'("  Allocated memory (IntWork):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_IntWork
      WRITE(LUPRI,'("  Allocated memory (overlap):         ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_overlap
      WRITE(LUPRI,'("  Allocated memory (ODitem):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_ODitem
      WRITE(LUPRI,'("  Allocated memory (lstensor):        ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_lstensor
      WRITE(LUPRI,'("  Allocated memory (FMM   ):          ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_FMM
#ifdef MOD_UNRELEASED
      WRITE(LUPRI,'("  Allocated memory (lvec_data):       ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_lvec_data
      WRITE(LUPRI,'("  Allocated memory (lattice_cell):    ",i9," byte  &
         &- Should be zero - otherwise a leakage is present")') mem_tp_allocated_lattice_cell
#endif

      call print_maxmem(lupri,max_mem_tp_used_global,'TOTAL')
#ifdef VAR_SCALAPACK
      call print_maxmem(lupri,max_mem_tp_used_type_matrix,'type(matrix) on Master node')
      call print_maxmem(lupri,max_mem_tp_used_type_matrix_MPIFULL,'Type(matrix) across all MPI nodes')
#else
      call print_maxmem(lupri,max_mem_tp_used_type_matrix,'type(matrix)')
#endif
      CALL print_maxmem(lupri,max_mem_tp_used_real,'real(realk)')
      CALL print_maxmem(lupri,max_mem_tp_used_mpi,'real(mpi)')
      CALL print_maxmem(lupri,max_mem_tp_used_complex,'complex(complexk)')
      CALL print_maxmem(lupri,max_mem_tp_used_integer,'integer')
      CALL print_maxmem(lupri,max_mem_tp_used_logical,'logical')
      CALL print_maxmem(lupri,max_mem_tp_used_character,'character')
      CALL print_maxmem(lupri,max_mem_tp_used_AOBATCH,'AOBATCH')
      CALL print_maxmem(lupri,max_mem_tp_used_DECORBITAL,'DECORBITAL')
      CALL print_maxmem(lupri,max_mem_tp_used_DECFRAG,'DECFRAG')
      CALL print_maxmem(lupri,max_mem_tp_used_BATCHTOORB,'BATCHTOORB')
      CALL print_maxmem(lupri,max_mem_tp_used_DECAOBATCHINFO,'DECAOBATCHINFO')
      CALL print_maxmem(lupri,max_mem_tp_used_MYPOINTER,'MYPOINTER')
      CALL print_maxmem(lupri,max_mem_tp_used_fragmentAOS,'fragmentAOS')
      CALL print_maxmem(lupri,max_mem_tp_used_ARRAY2,'ARRAY2')
      CALL print_maxmem(lupri,max_mem_tp_used_ARRAY4,'ARRAY4')
      CALL print_maxmem(lupri,max_mem_tp_used_ARRAY,'ARRAY')
      CALL print_maxmem(lupri,max_mem_tp_used_PNOSpaceInfo,'PNOSpaceInfo')
      CALL print_maxmem(lupri,max_mem_tp_used_MP2DENS,'MP2DENS')
      CALL print_maxmem(lupri,max_mem_tp_used_TRACEBACK,'TRACEBACK')
      CALL print_maxmem(lupri,max_mem_tp_used_MP2GRAD,'MP2GRAD')
      CALL print_maxmem(lupri,max_mem_tp_used_ODBATCH,'ODBATCH')
      CALL print_maxmem(lupri,max_mem_tp_used_LSAOTENSOR,'LSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_tp_used_SLSAOTENSOR,'SLSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_tp_used_GLOBALLSAOTENSOR,'GLOBALLSAOTENSOR')
      CALL print_maxmem(lupri,max_mem_tp_used_ATOMTYPEITEM,'ATOMTYPEITEM')
      CALL print_maxmem(lupri,max_mem_tp_used_ATOMITEM,'ATOMITEM')
      CALL print_maxmem(lupri,max_mem_tp_used_LSMATRIX,'LSMATRIX')
      CALL print_maxmem(lupri,max_mem_tp_used_overlapT,'OverlapType')


      CALL print_maxmem(lupri,max_mem_tp_used_linkshell,'linkshell')
      CALL print_maxmem(lupri,max_mem_tp_used_integrand,'integrand')
      CALL print_maxmem(lupri,max_mem_tp_used_integralitem,'integralitem')
      CALL print_maxmem(lupri,max_mem_tp_used_intwork,'IntWork')
      CALL print_maxmem(lupri,max_mem_tp_used_overlap,'Overlap')
      CALL print_maxmem(lupri,max_mem_tp_used_ODitem,'ODitem')
      CALL print_maxmem(lupri,max_mem_tp_used_lstensor,'LStensor')
      CALL print_maxmem(lupri,max_mem_tp_used_FMM,'FMM')
#ifdef MOD_UNRELEASED
      CALL print_maxmem(lupri,max_mem_tp_used_lvec_data,'Lvec_data')
      CALL print_maxmem(lupri,max_mem_tp_used_lattice_cell,'Lattice_cell')
#endif

      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,*)
      !$OMP END CRITICAL
   end subroutine stats_mem_tp

   subroutine Print_Memory_info(lupri,String)
      implicit none
      integer       :: lupri    
      Character*(*) ::  STRING
      character(len=4) :: Unit
      real(realk) :: Mem_with_correct_unit
      IF(MemModParamPrintMemory)THEN
         IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
            call GetMemwithUnit(Unit,max_mem_used_global,Mem_with_correct_unit)
            WRITE(lupri,'(A)')' '
            WRITE(lupri,'(3A,F6.2,1X,A)')'The maximum memory used ',TRIM(STRING),':     ',Mem_with_correct_unit,Unit
            WRITE(*,'(3A,F6.2,1X,A)')'The maximum memory used ',TRIM(STRING),':     ',Mem_with_correct_unit,Unit
            call GetMemwithUnit(Unit,mem_allocated_global,Mem_with_correct_unit)
            WRITE(lupri,'(3A,F6.2,1X,A)')'The current allocated memory ',TRIM(STRING),': ',Mem_with_correct_unit,Unit
            WRITE(*,'(3A,F6.2,1X,A)')'The current allocated memory ',TRIM(STRING),': ',Mem_with_correct_unit,Unit
            WRITE(*,'(A)')'-----------------------------------------------------------------------------------'
            if (mem_allocated_type_matrix > 0_long) call print_mem_alloc(lupri,mem_allocated_type_matrix,'type(matrix)')
         if (mem_allocated_real > 0_long) call print_mem_alloc(lupri,mem_allocated_real,'real')
         if (mem_allocated_mpi > 0_long) call print_mem_alloc(lupri,mem_allocated_mpi,'MPI')
         if (mem_allocated_complex > 0_long) call print_mem_alloc(lupri,mem_allocated_complex,'complex')
         if (mem_allocated_integer > 0_long) call print_mem_alloc(lupri,mem_allocated_integer,'integer')
         if (mem_allocated_logical > 0_long) call print_mem_alloc(lupri,mem_allocated_logical,'logical')
         if (mem_allocated_character > 0_long) call print_mem_alloc(lupri,mem_allocated_character,'character')
         if (mem_allocated_AOBATCH > 0_long) call print_mem_alloc(lupri,mem_allocated_AOBATCH,'AOBATCH')
         if (mem_allocated_DECORBITAL > 0_long) call print_mem_alloc(lupri,mem_allocated_DECORBITAL,'DECORBITAL')
         if (mem_allocated_DECFRAG > 0_long) call print_mem_alloc(lupri,mem_allocated_DECFRAG,'DECFRAG')
         if (mem_allocated_BATCHTOORB > 0_long) call print_mem_alloc(lupri,mem_allocated_BATCHTOORB,'BATCHTOORB')
         if (mem_allocated_DECAOBATCHINFO > 0_long) call print_mem_alloc(lupri,mem_allocated_DECAOBATCHINFO,'DECAOBATCHINFO')
         if (mem_allocated_MYPOINTER > 0_long) call print_mem_alloc(lupri,mem_allocated_MYPOINTER,'MYPOINTER')
         if (mem_allocated_fragmentAOS > 0_long) call print_mem_alloc(lupri,mem_allocated_fragmentAOS,'fragmentAOS')
         if (mem_allocated_ARRAY2 > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY2,'ARRAY2')
         if (mem_allocated_ARRAY4 > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY4,'ARRAY4')
         if (mem_allocated_ARRAY > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY,'ARRAY')
         if (mem_allocated_MP2DENS > 0_long) call print_mem_alloc(lupri,mem_allocated_MP2DENS,'MP2DENS')
         if (mem_allocated_TRACEBACK > 0_long) call print_mem_alloc(lupri,mem_allocated_TRACEBACK,'TRACEBACK')
         if (mem_allocated_MP2GRAD > 0_long) call print_mem_alloc(lupri,mem_allocated_MP2GRAD,'MP2GRAD')
         if (mem_allocated_ODBATCH > 0_long) call print_mem_alloc(lupri,mem_allocated_ODBATCH,'ODBATCH')
         if (mem_allocated_LSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_LSAOTENSOR,'LSAOTENSOR')
         if (mem_allocated_SLSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_SLSAOTENSOR,'SLSAOTENSOR')
         if (mem_allocated_GLOBALLSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_GLOBALLSAOTENSOR,'GLOBALLSAOTENSOR')
         if (mem_allocated_ATOMTYPEITEM > 0_long) call print_mem_alloc(lupri,mem_allocated_ATOMTYPEITEM,'ATOMTYPEITEM')
         if (mem_allocated_ATOMITEM > 0_long) call print_mem_alloc(lupri,mem_allocated_ATOMITEM,'ATOMITEM')
         if (mem_allocated_LSMATRIX > 0_long) call print_mem_alloc(lupri,mem_allocated_LSMATRIX,'LSMATRIX')
         if (mem_allocated_overlapT > 0_long) call print_mem_alloc(lupri,mem_allocated_overlapT,'overlapType')
         if (mem_allocated_FMM > 0_long) call print_mem_alloc(lupri,mem_allocated_FMM,'FMM   ')
         WRITE(lupri,'(A)')' '
         call LsTraceBack(STRING)
      ENDIF
   ENDIF
end subroutine Print_Memory_info

subroutine GetMemWithUnit(Unit,Mem,Mem_with_correct_unit)
   implicit none
   character(len=4) :: Unit    
   integer(KIND=long) :: Mem
   real(realk) :: Mem_with_correct_unit

   if (mem < 100) then !Divide by 1 to typecast real
      Unit = 'Byte' 
      Mem_with_correct_unit = mem/1.0E0_realk
   else if (mem < 1000000) then
      Unit = 'kB  ' 
      Mem_with_correct_unit = mem/1000.0E0_realk
   else if (mem < 1000000000) then
      Unit = 'MB  ' 
      Mem_with_correct_unit = mem/(1000.0E0_realk*1000.0E0_realk)
   else
      Unit = 'GB  ' 
      Mem_with_correct_unit = mem/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
   endif
end subroutine GetMemWithUnit

!> \brief For given SCF it, print amount of memory allocated for different data types.
!> \author S. Host
!> \date 2009
subroutine scf_stats_debug_mem(lupri,it)
   implicit none
   !> Current SCF iteration
   integer, intent(in) :: it
   !> Logical unit number for output file
   integer,intent(in) :: lupri
   IF(PrintSCFmemory)THEN
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,'("                  Memory statistics, iteration", i4)') it
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      call print_maxmem(lupri,max_mem_used_global,'TOTAL')
      if (max_mem_used_type_matrix > 0_long) call print_maxmem(lupri,max_mem_used_type_matrix,'type(matrix)')
      if (max_mem_used_real > 0_long) call print_maxmem(lupri,max_mem_used_real,'real')
      if (max_mem_used_mpi > 0_long) call print_maxmem(lupri,max_mem_used_mpi,'MPI')
      if (max_mem_used_complex > 0_long) call print_maxmem(lupri,max_mem_used_complex,'complex')
      if (max_mem_used_integer > 0_long) call print_maxmem(lupri,max_mem_used_integer,'integer')
      if (max_mem_used_logical > 0_long) call print_maxmem(lupri,max_mem_used_logical,'logical')

      if (max_mem_used_character > 0_long) call print_maxmem(lupri,max_mem_used_character,'character')

      if (max_mem_used_AOBATCH > 0_long) call print_maxmem(lupri,max_mem_used_AOBATCH,'AOBATCH')
      if (max_mem_used_DECORBITAL > 0_long) call print_maxmem(lupri,max_mem_used_DECORBITAL,'DECORBITAL')
      if (max_mem_used_DECFRAG > 0_long) call print_maxmem(lupri,max_mem_used_DECFRAG,'DECFRAG')
      if (max_mem_used_BATCHTOORB > 0_long) call print_maxmem(lupri,max_mem_used_BATCHTOORB,'BATCHTOORB')
      if (max_mem_used_DECAOBATCHINFO > 0_long) call print_maxmem(lupri,max_mem_used_DECAOBATCHINFO,'DECAOBATCHINFO')
      if (max_mem_used_MYPOINTER > 0_long) call print_maxmem(lupri,max_mem_used_MYPOINTER,'MYPOINTER')
      if (max_mem_used_fragmentAOS > 0_long) call print_maxmem(lupri,max_mem_used_fragmentAOS,'fragmentAOS')
      if (max_mem_used_ARRAY2 > 0_long) call print_maxmem(lupri,max_mem_used_ARRAY2,'ARRAY2')
      if (max_mem_used_ARRAY4 > 0_long) call print_maxmem(lupri,max_mem_used_ARRAY4,'ARRAY4')
      if (max_mem_used_ARRAY > 0_long) call print_maxmem(lupri,max_mem_used_ARRAY,'ARRAY')
      if (max_mem_used_PNOSpaceInfo > 0_long) call print_maxmem(lupri,max_mem_used_PNOSpaceInfo,'PNOSpaceInfo')
      if (max_mem_used_MP2DENS > 0_long) call print_maxmem(lupri,max_mem_used_MP2DENS,'MP2DENS')
      if (max_mem_used_TRACEBACK > 0_long) call print_maxmem(lupri,max_mem_used_TRACEBACK,'TRACEBACK')
      if (max_mem_used_MP2GRAD > 0_long) call print_maxmem(lupri,max_mem_used_MP2GRAD,'MP2GRAD')
      if (max_mem_used_ODBATCH > 0_long) call print_maxmem(lupri,max_mem_used_ODBATCH,'ODBATCH')
      if (max_mem_used_LSAOTENSOR > 0_long) call print_maxmem(lupri,max_mem_used_LSAOTENSOR,'LSAOTENSOR')
      if (max_mem_used_SLSAOTENSOR > 0_long) call print_maxmem(lupri,max_mem_used_SLSAOTENSOR,'SLSAOTENSOR')
      if (max_mem_used_GLOBALLSAOTENSOR > 0_long) call print_maxmem(lupri,max_mem_used_GLOBALLSAOTENSOR,'GLOBALLSAOTENSOR')
      if (max_mem_used_ATOMTYPEITEM > 0_long) call print_maxmem(lupri,max_mem_used_ATOMTYPEITEM,'ATOMTYPEITEM')
      if (max_mem_used_ATOMITEM > 0_long) call print_maxmem(lupri,max_mem_used_ATOMITEM,'ATOMITEM')
      if (max_mem_used_LSMATRIX > 0_long) call print_maxmem(lupri,max_mem_used_LSMATRIX,'LSMATRIX')
      if (max_mem_used_overlapT > 0_long) call print_maxmem(lupri,max_mem_used_overlapT,'overlapType')


      if (max_mem_used_linkshell > 0_long) call print_maxmem(lupri,max_mem_used_linkshell,'linkshell')
      if (max_mem_used_integrand > 0_long) call print_maxmem(lupri,max_mem_used_integrand,'integrand')
      if (max_mem_used_integralitem > 0_long) call print_maxmem(lupri,max_mem_used_integralitem,'integralitem')
      if (max_mem_used_IntWork > 0_long) call print_maxmem(lupri,max_mem_used_IntWork,'IntWork')
      if (max_mem_used_overlap > 0_long) call print_maxmem(lupri,max_mem_used_overlap,'overlap')
      if (max_mem_used_ODitem > 0_long) call print_maxmem(lupri,max_mem_used_ODitem,'ODitem')
      if (max_mem_used_lstensor > 0_long) call print_maxmem(lupri,max_mem_used_lstensor,'lstensor')
      if (max_mem_used_FMM > 0_long) call print_maxmem(lupri,max_mem_used_FMM,'FMM    ')
#ifdef MOD_UNRELEASED
      if (max_mem_used_lvec_data > 0_long) call print_maxmem(lupri,max_mem_used_lvec_data,'lvec_data')
      if (max_mem_used_lattice_cell > 0_long) call print_maxmem(lupri,max_mem_used_lattice_cell,'lattice_cell')
#endif
      WRITE(LUPRI,*)
      call print_mem_alloc(lupri,mem_allocated_global,'TOTAL')
      if (mem_allocated_type_matrix > 0_long) call print_mem_alloc(lupri,mem_allocated_type_matrix,'type(matrix)')
      if (mem_allocated_real > 0_long) call print_mem_alloc(lupri,mem_allocated_real,'real')
      if (mem_allocated_mpi > 0_long) call print_mem_alloc(lupri,mem_allocated_mpi,'MPI')
      if (mem_allocated_complex > 0_long) call print_mem_alloc(lupri,mem_allocated_complex,'complex')
      if (mem_allocated_integer > 0_long) call print_mem_alloc(lupri,mem_allocated_integer,'integer')
      if (mem_allocated_logical > 0_long) call print_mem_alloc(lupri,mem_allocated_logical,'logical')
      if (mem_allocated_character > 0_long) call print_mem_alloc(lupri,mem_allocated_character,'character')
      if (mem_allocated_AOBATCH > 0_long) call print_mem_alloc(lupri,mem_allocated_AOBATCH,'AOBATCH')
      if (mem_allocated_DECORBITAL > 0_long) call print_mem_alloc(lupri,mem_allocated_DECORBITAL,'DECORBITAL')
      if (mem_allocated_DECFRAG > 0_long) call print_mem_alloc(lupri,mem_allocated_DECFRAG,'DECFRAG')
      if (mem_allocated_BATCHTOORB > 0_long) call print_mem_alloc(lupri,mem_allocated_BATCHTOORB,'BATCHTOORB')
      if (mem_allocated_DECAOBATCHINFO > 0_long) call print_mem_alloc(lupri,mem_allocated_DECAOBATCHINFO,'DECAOBATCHINFO')
      if (mem_allocated_MYPOINTER > 0_long) call print_mem_alloc(lupri,mem_allocated_MYPOINTER,'MYPOINTER')
      if (mem_allocated_fragmentAOS > 0_long) call print_mem_alloc(lupri,mem_allocated_fragmentAOS,'fragmentAOS')
      if (mem_allocated_ARRAY2 > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY2,'ARRAY2')
      if (mem_allocated_ARRAY4 > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY4,'ARRAY4')
      if (mem_allocated_ARRAY > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY,'ARRAY')
      if (mem_allocated_PNOSpaceInfo > 0_long) call print_mem_alloc(lupri,mem_allocated_PNOSpaceInfo,'PNOSpaceInfo')
      if (mem_allocated_MP2DENS > 0_long) call print_mem_alloc(lupri,mem_allocated_MP2DENS,'MP2DENS')
      if (mem_allocated_TRACEBACK > 0_long) call print_mem_alloc(lupri,mem_allocated_TRACEBACK,'TRACEBACK')
      if (mem_allocated_MP2GRAD > 0_long) call print_mem_alloc(lupri,mem_allocated_MP2GRAD,'MP2GRAD')
      if (mem_allocated_ODBATCH > 0_long) call print_mem_alloc(lupri,mem_allocated_ODBATCH,'ODBATCH')
      if (mem_allocated_LSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_LSAOTENSOR,'LSAOTENSOR')
      if (mem_allocated_SLSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_SLSAOTENSOR,'SLSAOTENSOR')
      if (mem_allocated_GLOBALLSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_GLOBALLSAOTENSOR,'GLOBALLSAOTENSOR')
      if (mem_allocated_ATOMTYPEITEM > 0_long) call print_mem_alloc(lupri,mem_allocated_ATOMTYPEITEM,'ATOMTYPEITEM')
      if (mem_allocated_ATOMITEM > 0_long) call print_mem_alloc(lupri,mem_allocated_ATOMITEM,'ATOMITEM')
      if (mem_allocated_LSMATRIX > 0_long) call print_mem_alloc(lupri,mem_allocated_LSMATRIX,'LSMATRIX')
      if (mem_allocated_overlapT > 0_long) call print_mem_alloc(lupri,mem_allocated_overlapT,'overlapType')
      if (mem_allocated_FMM > 0_long) call print_mem_alloc(lupri,mem_allocated_FMM,'FMM   ')
      
      
      if (mem_allocated_linkshell > 0_long) call print_mem_alloc(lupri,mem_allocated_linkshell,'linkshell')
      if (mem_allocated_integrand > 0_long) call print_mem_alloc(lupri,mem_allocated_integrand,'integrand')
      if (mem_allocated_integralitem > 0_long) call print_mem_alloc(lupri,mem_allocated_integralitem,'integralitem')

      if (mem_allocated_overlap > 0_long) call print_mem_alloc(lupri,mem_allocated_overlap,'overlap')
      if (mem_allocated_ODitem > 0_long) call print_mem_alloc(lupri,mem_allocated_ODitem,'ODitem')
      if (mem_allocated_lstensor > 0_long) call print_mem_alloc(lupri,mem_allocated_lstensor,'lstensor')
#ifdef MOD_UNRELEASED
      if (mem_allocated_lvec_data > 0_long) call print_mem_alloc(lupri,mem_allocated_lvec_data,'lvec_data')
      if (mem_allocated_lattice_cell > 0_long) call print_mem_alloc(lupri,mem_allocated_lattice_cell,'lattice_cell')
#endif
      WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
      WRITE(LUPRI,*)
   ENDIF
end subroutine scf_stats_debug_mem

!> \brief status information printout. Print amount of memory allocated for different data types.
!> \author T. Kjaergaard
!> \date 2009
subroutine debug_mem_stats(lupri)
   implicit none
   !> Logical unit number for output file
   integer,intent(in) :: lupri

   WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
   WRITE(LUPRI,'("                  Debug Memory Statistics              ")')
   WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
   call print_maxmem(lupri,max_mem_used_global,'TOTAL')
   if (max_mem_used_type_matrix > 0_long) call print_maxmem(lupri,max_mem_used_type_matrix,'type(matrix)')
   if (max_mem_used_real > 0_long) call print_maxmem(lupri,max_mem_used_real,'real')
   if (max_mem_used_mpi > 0_long) call print_maxmem(lupri,max_mem_used_mpi,'MPI')
   if (max_mem_used_complex > 0_long) call print_maxmem(lupri,max_mem_used_complex,'complex')
   if (max_mem_used_integer > 0_long) call print_maxmem(lupri,max_mem_used_integer,'integer')
   if (max_mem_used_logical > 0_long) call print_maxmem(lupri,max_mem_used_logical,'logical')
   if (max_mem_used_character > 0_long) call print_maxmem(lupri,max_mem_used_character,'character')
   if (max_mem_used_AOBATCH > 0_long) call print_maxmem(lupri,max_mem_used_AOBATCH,'AOBATCH')
   if (max_mem_used_DECORBITAL > 0_long) call print_maxmem(lupri,max_mem_used_DECORBITAL,'DECORBITAL')
   if (max_mem_used_DECFRAG > 0_long) call print_maxmem(lupri,max_mem_used_DECFRAG,'DECFRAG')
   if (max_mem_used_BATCHTOORB > 0_long) call print_maxmem(lupri,max_mem_used_BATCHTOORB,'BATCHTOORB')
   if (max_mem_used_DECAOBATCHINFO > 0_long) call print_maxmem(lupri,max_mem_used_DECAOBATCHINFO,'DECAOBATCHINFO')
   if (max_mem_used_MYPOINTER > 0_long) call print_maxmem(lupri,max_mem_used_MYPOINTER,'MYPOINTER')
   if (max_mem_used_fragmentAOS > 0_long) call print_maxmem(lupri,max_mem_used_fragmentAOS,'fragmentAOS')
   if (max_mem_used_ARRAY2 > 0_long) call print_maxmem(lupri,max_mem_used_ARRAY2,'ARRAY2')
   if (max_mem_used_ARRAY4 > 0_long) call print_maxmem(lupri,max_mem_used_ARRAY4,'ARRAY4')
   if (max_mem_used_ARRAY > 0_long) call print_maxmem(lupri,max_mem_used_ARRAY,'ARRAY')
   if (max_mem_used_PNOSpaceInfo > 0_long) call print_maxmem(lupri,max_mem_used_PNOSpaceInfo,'PNOSpaceInfo')
   if (max_mem_used_MP2DENS > 0_long) call print_maxmem(lupri,max_mem_used_MP2DENS,'MP2DENS')
   if (max_mem_used_TRACEBACK > 0_long) call print_maxmem(lupri,max_mem_used_TRACEBACK,'TRACEBACK')
   if (max_mem_used_MP2GRAD > 0_long) call print_maxmem(lupri,max_mem_used_MP2GRAD,'MP2GRAD')
   if (max_mem_used_ODBATCH > 0_long) call print_maxmem(lupri,max_mem_used_ODBATCH,'ODBATCH')
   if (max_mem_used_LSAOTENSOR > 0_long) call print_maxmem(lupri,max_mem_used_LSAOTENSOR,'LSAOTENSOR')
   if (max_mem_used_SLSAOTENSOR > 0_long) call print_maxmem(lupri,max_mem_used_SLSAOTENSOR,'SLSAOTENSOR')
   if (max_mem_used_GLOBALLSAOTENSOR > 0_long) call print_maxmem(lupri,max_mem_used_GLOBALLSAOTENSOR,'GLOBALLSAOTENSOR')
   if (max_mem_used_ATOMTYPEITEM > 0_long) call print_maxmem(lupri,max_mem_used_ATOMTYPEITEM,'ATOMTYPEITEM')
   if (max_mem_used_ATOMITEM > 0_long) call print_maxmem(lupri,max_mem_used_ATOMITEM,'ATOMITEM')
   if (max_mem_used_LSMATRIX > 0_long) call print_maxmem(lupri,max_mem_used_LSMATRIX,'LSMATRIX')
   if (max_mem_used_overlapT > 0_long) call print_maxmem(lupri,max_mem_used_overlapT,'overlapType')
   
   if (max_mem_used_linkshell > 0_long) call print_maxmem(lupri,max_mem_used_linkshell,'linkshell')
   if (max_mem_used_integrand > 0_long) call print_maxmem(lupri,max_mem_used_integrand,'integrand')
   if (max_mem_used_integralitem > 0_long) call print_maxmem(lupri,max_mem_used_integralitem,'integralitem')
   if (max_mem_used_IntWork > 0_long) call print_maxmem(lupri,max_mem_used_IntWork,'IntWork')
   if (max_mem_used_overlap > 0_long) call print_maxmem(lupri,max_mem_used_overlap,'overlap')
   if (max_mem_used_ODitem > 0_long) call print_maxmem(lupri,max_mem_used_ODitem,'ODitem')
   if (max_mem_used_lstensor > 0_long) call print_maxmem(lupri,max_mem_used_lstensor,'lstensor')
   if (max_mem_used_FMM > 0_long) call print_maxmem(lupri,max_mem_used_FMM,'FMM    ')
#ifdef MOD_UNRELEASED
   if (max_mem_used_lvec_data > 0_long) call print_maxmem(lupri,max_mem_used_lvec_data,'lvec_data')
   if (max_mem_used_lattice_cell > 0_long) call print_maxmem(lupri,max_mem_used_lattice_cell,'lattice_cell')
#endif
   WRITE(LUPRI,*)
   call print_mem_alloc(lupri,mem_allocated_global,'TOTAL')
   if (mem_allocated_type_matrix > 0_long) call print_mem_alloc(lupri,mem_allocated_type_matrix,'type(matrix)')
   if (mem_allocated_real > 0_long) call print_mem_alloc(lupri,mem_allocated_real,'real')
   if (mem_allocated_mpi > 0_long) call print_mem_alloc(lupri,mem_allocated_mpi,'MPI')
   if (mem_allocated_complex > 0_long) call print_mem_alloc(lupri,mem_allocated_complex,'complex')
   if (mem_allocated_integer > 0_long) call print_mem_alloc(lupri,mem_allocated_integer,'integer')
   if (mem_allocated_logical > 0_long) call print_mem_alloc(lupri,mem_allocated_logical,'logical')
   if (mem_allocated_character > 0_long) call print_mem_alloc(lupri,mem_allocated_character,'character')
   if (mem_allocated_AOBATCH > 0_long) call print_mem_alloc(lupri,mem_allocated_AOBATCH,'AOBATCH')
   if (mem_allocated_DECORBITAL > 0_long) call print_mem_alloc(lupri,mem_allocated_DECORBITAL,'DECORBITAL')
   if (mem_allocated_DECFRAG > 0_long) call print_mem_alloc(lupri,mem_allocated_DECFRAG,'DECFRAG')
   if (mem_allocated_BATCHTOORB > 0_long) call print_mem_alloc(lupri,mem_allocated_BATCHTOORB,'BATCHTOORB')
   if (mem_allocated_DECAOBATCHINFO > 0_long) call print_mem_alloc(lupri,mem_allocated_DECAOBATCHINFO,'DECAOBATCHINFO')
   if (mem_allocated_MYPOINTER > 0_long) call print_mem_alloc(lupri,mem_allocated_MYPOINTER,'MYPOINTER')
   if (mem_allocated_fragmentAOS > 0_long) call print_mem_alloc(lupri,mem_allocated_fragmentAOS,'fragmentAOS')
   if (mem_allocated_ARRAY2 > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY2,'ARRAY2')
   if (mem_allocated_ARRAY4 > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY4,'ARRAY4')
   if (mem_allocated_ARRAY > 0_long) call print_mem_alloc(lupri,mem_allocated_ARRAY,'ARRAY')
   if (mem_allocated_PNOSpaceInfo > 0_long) call print_mem_alloc(lupri,mem_allocated_PNOSpaceInfo,'PNOSpaceInfo')
   if (mem_allocated_MP2DENS > 0_long) call print_mem_alloc(lupri,mem_allocated_MP2DENS,'MP2DENS')
   if (mem_allocated_TRACEBACK > 0_long) call print_mem_alloc(lupri,mem_allocated_TRACEBACK,'TRACEBACK')
   if (mem_allocated_MP2GRAD > 0_long) call print_mem_alloc(lupri,mem_allocated_MP2GRAD,'MP2GRAD')
   if (mem_allocated_ODBATCH > 0_long) call print_mem_alloc(lupri,mem_allocated_ODBATCH,'ODBATCH')
   if (mem_allocated_LSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_LSAOTENSOR,'LSAOTENSOR')
   if (mem_allocated_SLSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_SLSAOTENSOR,'SLSAOTENSOR')
   if (mem_allocated_GLOBALLSAOTENSOR > 0_long) call print_mem_alloc(lupri,mem_allocated_GLOBALLSAOTENSOR,'GLOBALLSAOTENSOR')
   if (mem_allocated_ATOMTYPEITEM > 0_long) call print_mem_alloc(lupri,mem_allocated_ATOMTYPEITEM,'ATOMTYPEITEM')
   if (mem_allocated_ATOMITEM > 0_long) call print_mem_alloc(lupri,mem_allocated_ATOMITEM,'ATOMITEM')
   if (mem_allocated_LSMATRIX > 0_long) call print_mem_alloc(lupri,mem_allocated_LSMATRIX,'LSMATRIX')
   if (mem_allocated_overlapT > 0_long) call print_mem_alloc(lupri,mem_allocated_overlapT,'overlapType')
   
   
   if (mem_allocated_linkshell > 0_long) call print_mem_alloc(lupri,mem_allocated_linkshell,'linkshell')
   if (mem_allocated_integrand > 0_long) call print_mem_alloc(lupri,mem_allocated_integrand,'integrand')
   if (mem_allocated_integralitem > 0_long) call print_mem_alloc(lupri,mem_allocated_integralitem,'integralitem')
   if (mem_allocated_IntWork > 0_long) call print_mem_alloc(lupri,mem_allocated_IntWork,'IntWork')
   if (mem_allocated_overlap > 0_long) call print_mem_alloc(lupri,mem_allocated_overlap,'overlap')
   if (mem_allocated_ODitem > 0_long) call print_mem_alloc(lupri,mem_allocated_ODitem,'ODitem')
   if (mem_allocated_lstensor > 0_long) call print_mem_alloc(lupri,mem_allocated_lstensor,'lstensor')
   if (mem_allocated_FMM > 0_long) call print_mem_alloc(lupri,mem_allocated_FMM,'FMM   ')
#ifdef MOD_UNRELEASED
   if (mem_allocated_lvec_data > 0_long) call print_mem_alloc(lupri,mem_allocated_lvec_data,'lvec_data')
   if (mem_allocated_lattice_cell > 0_long) call print_mem_alloc(lupri,mem_allocated_lattice_cell,'lattice_cell')
#endif
   WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
   WRITE(LUPRI,*)

  end subroutine debug_mem_stats

  !> \brief Print routine for memory statistics (max. allocated memory).
  !> \author S. Host
  !> \date 2009
  subroutine print_maxmem(lupri,max_mem_used,STRING)
     implicit none
     !> Logical unit number for output file
     integer :: lupri
     !> Amount of memory used
     integer(KIND=long) :: max_mem_used
     !> Name of data type for which memory usage is printed
     Character*(*)    ::  STRING
     character(len=27) :: trimstring
     !Character(len=3) ::  FORMATSTRING1
     !Character(len=23) ::  FORMATSTRING2
     trimstring = adjustl(string)
     if (max_mem_used < 100) then !Divide by 1 to typecast real
        WRITE(LUPRI,'("  Max allocated memory, ", a27, f9.3, " Byte")') trimstring, max_mem_used/1.0E0_realk
     else if (max_mem_used < 1000000) then
        WRITE(LUPRI,'("  Max allocated memory, ", a27, f9.3, " kB")') trimstring, max_mem_used/1000.0E0_realk
     else if (max_mem_used < 1000000000) then
        WRITE(LUPRI,'("  Max allocated memory, ", a27, f9.3, " MB")') trimstring, max_mem_used/(1000.0E0_realk*1000.0E0_realk)
     else
        WRITE(LUPRI,'("  Max allocated memory, ", a27, f9.3, " GB")') trimstring, &
           & max_mem_used/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
     endif
  end subroutine print_maxmem

  !> \brief Print routine for memory statistics (current allocated memory).
  !> \author S. Host
  !> \date 2009
  subroutine print_mem_alloc(lupri,mem_used,STRING)
     implicit none
     !> Logical unit number for output file
     integer :: lupri
     !> Amount of memory used
     integer(KIND=long) :: mem_used
     !> Name of data type for which memory usage is printed
     Character*(*)    ::  STRING
     character(len=20) :: trimstring
     trimstring = adjustl(string)
     if (mem_used < 100) then !Divide by 1 to typecast real
        WRITE(LUPRI,'("  Allocated memory, ", a24, f9.3, " Byte")')&
           &trimstring, mem_used/1.0E0_realk
     else if (mem_used < 1000000) then
        WRITE(LUPRI,'("  Allocated memory, ", a24, f9.3, " kB")')&
           &trimstring, mem_used/1000.0E0_realk
     else if (mem_used < 1000000000) then
        WRITE(LUPRI,'("  Allocated memory, ", a24, f9.3, " MB")')&
           &trimstring, mem_used/(1000.0E0_realk*1000.0E0_realk)
     else
        WRITE(LUPRI,'("  Allocated memory, ", a24, f9.3, " GB")')&
           &trimstring, mem_used/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
     endif
  end subroutine print_mem_alloc

  !SUBROUTINE mem_print(iunit)
  !implicit none
  !integer,intent(in) :: iunit
  !  call print_mem('Global',mem_allocated_global,iunit)
  !END SUBROUTINE mem_print
  !
  !SUBROUTINE mem_max_print(iunit)
  !implicit none
  !integer,intent(in) :: iunit
  !  call print_mem('Mamimum',max_mem_used_global,iunit)
  !END SUBROUTINE mem_max_print
  !
  !SUBROUTINE print_mem(id,mem,iunit)
  !implicit none
  !character*(*),intent(in) :: id
  !integer,intent(in) :: iunit
  !integer(kind=long), intent(in) :: mem
  !!
  !character(80) :: memory_string
  !
  !IF (mem.GT.(1000**3)) THEN
  !  WRITE(memory_string,'(1X,I10,A1,I3,1X,A6)') int(mem/1000E0_realk**3),'.', &
  !     &   int((mem-int(mem/1000E0_realk**3)*1000E0_realk**3)/1000E0_realk**3*1000),'Gbytes'
  !ELSE IF (mem.GT.(1000**2)) THEN
  !  WRITE(memory_string,'(1X,I10,A1,I3,1X,A6)') int(mem/1000E0_realk**2),'.', &
  !     &   int((mem-int(mem/1000E0_realk**2)*1000E0_realk**2)/1000E0_realk**2*1000),'Mbytes'
  !ELSE IF (mem.GT.(1000)) THEN
  !  WRITE(memory_string,'(1X,I10,A1,I3,1X,A6)') int(mem/1000E0_realk),'.', &
  !     &   int((mem-int(mem/1000E0_realk)*1000E0_realk)/1000E0_realk*1000),'kbytes'
  !ELSE
  !  WRITE(memory_string,'(1X,I14,1X,A6)') mem,' bytes'
  !ENDIF
  !WRITE(IUNIT,'(3X,A,X,A)') id,memory_string
  !WRITE(*,'(3X,A10,A14,A40)') id,' memory usage:',memory_string
  !END SUBROUTINE print_mem

  !> \brief When memory allocation fails, this subroutine should be called.
  !> It prints a memory statistics summary and quits LSDALTON.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine memory_error_quit(mylabel,error_size)
     implicit none
     !> Label for routine where memory allocation failed
     integer(kind=8),intent(in) :: error_size
     character(*),intent(in) :: mylabel
     !> Unit number for lsquit output 
     integer :: myoutput
     character(8) :: ERR,GLOB

     myoutput=6

     ! Memory statistics

     IF (error_size.LT.0) THEN
        write(ERR,'(A8)') ' < zero '
     ELSEIF (error_size.LT.1000) THEN
        write(ERR,'(F5.1,A3)') error_size*1E0_realk," B "
     ELSEIF (error_size.LT.1000000) THEN
        write(ERR,'(F5.1,A3)') error_size*1E-3_realk," kB"
     ELSEIF (error_size.LT.1000000000) THEN
        write(ERR,'(F5.1,A3)') error_size*1E-6_realk," MB"
#ifdef VAR_INT64
     ELSEIF (error_size.LT.1000000000000) THEN
        write(ERR,'(F5.1,A3)') error_size*1E-9_realk," GB"
     ELSEIF (error_size.LT.1000000000000000) THEN
        write(ERR,'(F5.1,A3)') error_size*1E-12_realk," TB"
     ELSEIF (error_size.LT.1000000000000000000) THEN
        write(ERR,'(F5.1,A3)') error_size*1E-15_realk," PB"
     ELSE
        write(ERR,'(F5.1,A3)') error_size*1E-18_realk," EB"
     ENDIF
#else
     ELSE
        write(ERR,'(F5.1,A3)') error_size*1E-9_realk," GB"
     ENDIF
#endif

     IF (mem_allocated_global.LT.0) THEN
        write(GLOB,'(A8)') ' < zero '
     ELSEIF (mem_allocated_global.LT.1000) THEN
        write(GLOB,'(F5.1,A3)') mem_allocated_global*1E0_realk," B "
     ELSEIF (mem_allocated_global.LT.1000000) THEN
        write(GLOB,'(F5.1,A3)') mem_allocated_global*1E-3_realk," kB"
     ELSEIF (mem_allocated_global.LT.1000000000) THEN
        write(GLOB,'(F5.1,A3)') mem_allocated_global*1E-6_realk," MB"
#ifdef VAR_INT64
     ELSEIF (mem_allocated_global.LT.1000000000000) THEN
        write(GLOB,'(F5.1,A3)') mem_allocated_global*1E-9_realk," GB"
     ELSEIF (mem_allocated_global.LT.1000000000000000) THEN
        write(GLOB,'(F5.1,A3)') mem_allocated_global*1E-12_realk," TB"
     ELSEIF (mem_allocated_global.LT.1000000000000000000) THEN
        write(GLOB,'(F5.1,A3)') mem_allocated_global*1E-15_realk," PB"
     ELSE
        write(GLOB,'(F5.1,A3)') mem_allocated_global*1E-18_realk," EB"
     ENDIF
#else
     ELSE
        write(GLOB,'(F5.1,A3)') mem_allocated_global*1E-9_realk," GB"
     ENDIF
#endif

     write(myoutput,*) 
     write(myoutput,*) 'LSDALTON is quitting because there is too little memory available!'
     write(myoutput,*) 'The program was trying to allocate ',ERR,&
        &' in addition to the ',GLOB,' already allocated.'
     write(myoutput,*) 'Increase available memory if possible (eg. through your submit script or you),'
     write(myoutput,*) 'may try to distribute memory over more nodes (eg. by using ScaLapack/PBLAS). See the'
     write(myoutput,*) 'LSDALTON manual or consult the Dalton Forum for details.'
     write(myoutput,*) 
     write(myoutput,*) 'LSDALTON provide no a priori estimate of the memory it needs, but the memory statistics'
     write(myoutput,*) 'from similar calculations on smaller systems may hopefully give you some clues as to what'
     write(myoutput,*) 'memory requirements are needed.'
     write(myoutput,*) 
     write(myoutput,*) 'Printing memory statistics for the current calculation before quitting LSDALTON...'
     
     call stats_mem(myoutput)
     flush(myoutput)

     ! Quit
     call lsquit(mylabel,-1)

   end subroutine memory_error_quit


!----- Background buffer handling -----!
subroutine mem_init_background_alloc(bytes)
   implicit none
   integer(kind=8), intent(in) :: bytes
   integer(kind=8) :: nelms

   if(buf_realk%init)then
      call lsquit("ERROR(mem_free_background_alloc): background buffer already initialized",-1)
   endif

   nelms = bytes/8_long
#ifdef VAR_MPI
   call mem_alloc(buf_realk%p,buf_realk%c,nelms)
#else
   call mem_alloc(buf_realk%p,nelms)
#endif

   buf_realk%init   = .true.
   buf_realk%offset = 0
   buf_realk%nmax   = nelms

   buf_realk%n      = 1

   buf_realk%c_addr = c_null_ptr
   buf_realk%c_mdel = c_null_ptr
   buf_realk%f_mdel = c_null_ptr
   buf_realk%e_mdel = 0
   buf_realk%n_mdel = 0
   buf_realk%l_mdel = .false.

end subroutine mem_init_background_alloc
subroutine mem_change_background_alloc(bytes,not_lazy)
   implicit none
   integer(kind=8), intent(in) :: bytes
   logical, intent(in), optional :: not_lazy
   integer :: i
   integer(kind=8) :: nelms
   logical :: change

   nelms = bytes

   !per default use lazy memory handling
   change = (nelms>buf_realk%nmax)
   if(present(not_lazy)) change = not_lazy

   if(buf_realk%n/=1)then
      do i = 1, buf_realk%n-1
         print *,"address not freed",i
      enddo
      call lsquit("ERROR(mem_change_background_alloc): pointers is still associated, not possible",-1)
   endif

   if(.not.buf_realk%init)then
      call lsquit("ERROR(mem_free_background_alloc): background buffer not initialized",-1)
   endif


   if( change )then
      print *,"mem_change_bg_alloc: Allocating ",bytes/(1024.0_realk**3),"GB"

#ifdef VAR_MPI
      call mem_dealloc(buf_realk%p,buf_realk%c)

      call mem_alloc(buf_realk%p,buf_realk%c,nelms)
#else
      call mem_dealloc(buf_realk%p)

      call mem_alloc(buf_realk%p,nelms)
#endif

      buf_realk%init   = .true.
      buf_realk%offset = 0
      buf_realk%nmax   = nelms

      buf_realk%n      = 1
      buf_realk%c_addr = c_null_ptr
      buf_realk%c_mdel = c_null_ptr
      buf_realk%f_mdel = c_null_ptr
      buf_realk%e_mdel = 0
      buf_realk%n_mdel = 0
      buf_realk%l_mdel = .false.
   endif

end subroutine mem_change_background_alloc
subroutine mem_free_background_alloc()
   implicit none
   integer :: i

   if(buf_realk%n/=1)then
      do i = 1, buf_realk%n-1
         print *,"address not freed",i
      enddo
      call lsquit("ERROR(mem_free_background_alloc): pointers is still associated",-1)
   endif

   if(.not.buf_realk%init)then
      call lsquit("ERROR(mem_free_background_alloc): background buffer not initialized",-1)
   endif

#ifdef VAR_MPI
   call mem_dealloc(buf_realk%p,buf_realk%c)
#else
   call mem_dealloc(buf_realk%p)
#endif

   buf_realk%init   = .false.
   buf_realk%offset = 0
   buf_realk%nmax   = 0

   buf_realk%n      = 1
   buf_realk%c_addr = c_null_ptr
   buf_realk%c_mdel = c_null_ptr
   buf_realk%f_mdel = c_null_ptr
   buf_realk%e_mdel = 0
   buf_realk%n_mdel = 0
   buf_realk%l_mdel = .false.

end subroutine mem_free_background_alloc


subroutine mem_clear_mdel_bg_buf()
   implicit none
   integer :: i
   if (buf_realk%l_mdel) then
      do i=1,buf_realk%n_mdel

         if(.not.c_associated(buf_realk%c_mdel(i),c_loc(buf_realk%p(buf_realk%f_addr(buf_realk%n-1)))))then
            call lsquit("ERROR(mem_clear_mdel_bg_buf): wrong sequence of&
               & deallocating, make sure you dealloc in the opposite seqence as&
               & allocating, also when using mark_deleted",-1)
         endif

         buf_realk%n = buf_realk%n - 1
         buf_realk%offset = buf_realk%offset - buf_realk%e_mdel(i)

      enddo
      
      buf_realk%c_mdel = c_null_ptr
      buf_realk%f_mdel = c_null_ptr
      buf_realk%e_mdel = 0
      buf_realk%n_mdel = 0
      buf_realk%l_mdel = .false.
   endif
end subroutine mem_clear_mdel_bg_buf

!ATTENTION, WHEN USING THIS STRUCTURE YOU NEED TO MAKE SURE THAT YOU DEALLOCATE
!IN THE OPPOSITE SEQUENCE OF ALLOCATING, OTHERWISE YOUR DATA WILL BE CORRUPTED
subroutine mem_pseudo_alloc_realk(p,n)
   implicit none
   real(realk), pointer :: p(:)
   integer(kind=8), intent(in) :: n

   if (buf_realk%offset+n > buf_realk%nmax)then
      print *,"Buffer Space (#elements):",buf_realk%nmax," Used:",buf_realk%offset," Requested:",n
      call lsquit("ERROR(mem_pseudo_alloc_realk): more requested than available",-1)
   endif

   p => buf_realk%p(buf_realk%offset+1:buf_realk%offset+n)

   buf_realk%f_addr(buf_realk%n) = buf_realk%offset+1
   buf_realk%c_addr(buf_realk%n) = c_loc(p(1))

   buf_realk%offset = buf_realk%offset+n


   buf_realk%n = buf_realk%n + 1

   if(buf_realk%n > max_n_pointers)then
      call lsquit("ERROR(mem_pseudo_alloc_realk): more pointers associated then currently supported, &
         &please change max_n_pointers in background_buffer.F90",-1)
   endif

end subroutine mem_pseudo_alloc_realk
subroutine mem_pseudo_alloc_realk2(p,n1,n2)
   implicit none
   real(realk), pointer :: p(:,:)
   integer(kind=8), intent(in) :: n1,n2
   integer(kind=8) :: nelms
   nelms = n1*n2

   if (buf_realk%offset+nelms > buf_realk%nmax)then
      print *,"Buffer Space (#elements):",buf_realk%nmax," Used:",buf_realk%offset," Requested:",nelms
      call lsquit("ERROR(mem_pseudo_alloc_realk2): more requested than available",-1)
   endif

#ifdef VAR_PTR_RESHAPE
   p(1:n1,1:n2) => buf_realk%p(buf_realk%offset+1:buf_realk%offset+nelms)
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
   call c_f_pointer(c_loc(buf_realk%p(buf_realk%offset+1)),p,[n1,n2])
#else
   call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

   buf_realk%f_addr(buf_realk%n) = buf_realk%offset+1
   buf_realk%c_addr(buf_realk%n) = c_loc(p(1,1))

   buf_realk%offset = buf_realk%offset+nelms


   buf_realk%n = buf_realk%n + 1

   if(buf_realk%n > max_n_pointers)then
      call lsquit("ERROR(mem_pseudo_alloc_realk2): more pointers associated then currently supported, &
         &please change max_n_pointers in background_buffer.F90",-1)
   endif

end subroutine mem_pseudo_alloc_realk2
subroutine mem_pseudo_alloc_realk3(p,n1,n2,n3)
   implicit none
   real(realk), pointer :: p(:,:,:)
   integer(kind=8), intent(in) :: n1,n2,n3
   integer(kind=8) :: nelms
   nelms = n1*n2*n3

   if (buf_realk%offset+nelms > buf_realk%nmax)then
      print *,"Buffer Space (#elements):",buf_realk%nmax," Used:",buf_realk%offset," Requested:",nelms
      call lsquit("ERROR(mem_pseudo_alloc_realk3): more requested than available",-1)
   endif

#ifdef VAR_PTR_RESHAPE
   p(1:n1,1:n2,1:n3) => buf_realk%p(buf_realk%offset+1:buf_realk%offset+nelms)
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
   call c_f_pointer(c_loc(buf_realk%p(buf_realk%offset+1)),p,[n1,n2,n3])
#else
   call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif


   buf_realk%f_addr(buf_realk%n) = buf_realk%offset+1
   buf_realk%c_addr(buf_realk%n) = c_loc(p(1,1,1))

   buf_realk%offset = buf_realk%offset+nelms
   buf_realk%n = buf_realk%n + 1

   if(buf_realk%n > max_n_pointers)then
      call lsquit("ERROR(mem_pseudo_alloc_realk2): more pointers associated then currently supported, &
         &please change max_n_pointers in background_buffer.F90",-1)
   endif

end subroutine mem_pseudo_alloc_realk3
subroutine mem_pseudo_dealloc_realk(p, mark_deleted)
   implicit none
   real(realk), pointer :: p(:)
   ! if a pointer is to be freed it may be marked as deleted, as long as there
   ! are pointers marked as deleted, no new pointers can be associated to the bg
   ! buffer
   logical, optional, intent(in) :: mark_deleted
   integer(kind=8) :: n
   logical :: md, last_assoc, del_assoc
   integer :: i
   type(c_ptr) :: loc1,loc2

   md=.false.
   if(present(mark_deleted))md=mark_deleted

   !no need to mark_deleted if it is the correct last element in the buffer
   last_assoc = c_associated(buf_realk%c_addr(buf_realk%n-1),c_loc(p(1)))
   del_assoc  = c_associated(buf_realk%f_mdel,               c_loc(p(1)))
   if(last_assoc) md = .false.

   if(md)then
      if(.not.buf_realk%l_mdel)then
         FindPos:do i=1,buf_realk%n
            loc1 = c_loc(p(1))
            loc2 = c_loc(buf_realk%p(buf_realk%f_addr(i)))
            if( c_associated(loc1,loc2))then
               buf_realk%f_mdel = c_loc(buf_realk%p(buf_realk%f_addr(i+1)))
               exit FindPos
            endif
         enddo FindPos
      endif
      buf_realk%l_mdel = .true.
      buf_realk%n_mdel = buf_realk%n_mdel + 1
   else
      buf_realk%n = buf_realk%n - 1
   endif

   if(buf_realk%n < 0)then
      call lsquit("ERROR(mem_pseudo_dealloc_realk): programming error, more&
      & pointers freed than associated",-1)
   endif

   if(.not.md.and..not.last_assoc)then
      call lsquit("ERROR(mem_pseudo_dealloc_realk): wrong sequence of &
      &deallocating, make sure you dealloc in the opposite seqence as allocating, &
      &otherwise you will corrupt your data",-1)
   endif

   n = size(p,kind=8)

   if (buf_realk%offset-n < 0)then
      call lsquit("ERROR(mem_pseudo_dealloc_realk): more freed than allocated",-1)
   endif


   if(md)then
      buf_realk%e_mdel(buf_realk%n_mdel) = n
      buf_realk%c_mdel(buf_realk%n_mdel) = c_loc(p(1))
   else
      buf_realk%offset = buf_realk%offset-n
   end if

   p => null()
   buf_realk%c_addr(buf_realk%n) = c_null_ptr

   !Right now, always clear the mark_deleted buffer at the next real
   !deallocation, otherwise bookeeping becomes more cumbersome
   !NOTE: the deassociation from the buffer happens after the real deallocation
   !and it has to be performed in the correct sequence to ensure the correctness
   !of the data in the buffer (without too much bookkeeping)
   if(del_assoc)then
      call mem_clear_mdel_bg_buf()
   endif

end subroutine mem_pseudo_dealloc_realk
subroutine mem_pseudo_dealloc_realk2(p)
   implicit none
   real(realk), pointer :: p(:,:)
   integer(kind=8) :: n

   buf_realk%n = buf_realk%n - 1

   if(buf_realk%n < 0)then
      call lsquit("ERROR(mem_pseudo_dealloc_realk2): programming error, more&
      & pointers freed than associated",-1)
   endif

   if( .not. c_associated(buf_realk%c_addr(buf_realk%n),c_loc(p(1,1))))then
      call lsquit("ERROR(mem_pseudo_dealloc_realk2): wrong sequence of &
      &deallocating, make sure you dealloc in the opposite seqence as allocating, &
      &otherwise you will corrupt your data",-1)
   endif

   n = size(p,kind=8)

   if (buf_realk%offset-n < 0)then
      call lsquit("ERROR(mem_pseudo_dealloc_realk2): more freed than allocated",-1)
   endif

   p => null()

   buf_realk%offset = buf_realk%offset-n

   buf_realk%c_addr(buf_realk%n) = c_null_ptr

end subroutine mem_pseudo_dealloc_realk2
subroutine mem_pseudo_dealloc_realk3(p)
   implicit none
   real(realk), pointer :: p(:,:,:)
   integer(kind=8) :: n

   buf_realk%n = buf_realk%n - 1

   if(buf_realk%n < 0)then
      call lsquit("ERROR(mem_pseudo_dealloc_realk3): programming error, more&
      & pointers freed than associated",-1)
   endif

   if( .not. c_associated(buf_realk%c_addr(buf_realk%n),c_loc(p(1,1,1))))then
      call lsquit("ERROR(mem_pseudo_dealloc_realk3): wrong sequence of &
      &deallocating, make sure you dealloc in the opposite seqence as allocating, &
      &otherwise you will corrupt your data",-1)
   endif

   n = size(p,kind=8)

   if (buf_realk%offset-n < 0)then
      call lsquit("ERROR(mem_pseudo_dealloc_realk3): more freed than allocated",-1)
   endif

   p => null()

   buf_realk%offset = buf_realk%offset-n

   buf_realk%c_addr(buf_realk%n) = c_null_ptr

end subroutine mem_pseudo_dealloc_realk3

function mem_is_background_buf_init() result(init)
   implicit none
   logical :: init
   init = buf_realk%init
end function mem_is_background_buf_init

function mem_get_bg_buf_n() result(n)
   implicit none
   integer(kind=8) :: n
   n = buf_realk%nmax
end function mem_get_bg_buf_n

function mem_get_bg_buf_free() result(n)
   implicit none
   integer(kind=8) :: n
   n = buf_realk%nmax-buf_realk%offset
end function mem_get_bg_buf_free

subroutine mem_pseudo_alloc_mpirealk(A,n,comm,local,simple) 
   implicit none
   integer(kind=8),intent(in)          :: n
   logical,intent(in),optional         :: local,simple
   type(mpi_realk)                     :: A
   integer(kind=ls_mpik),optional,intent(in)    ::comm
   integer (kind=ls_mpik)              :: IERR,info
   integer (kind=long)                 :: nsize
   character(120)                      :: errmsg
   logical                             :: loc,simp
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
#endif

   simp = .false.
   if(present(simple))simp = simple

   if(present(comm).and..not.simp)then
      call lsquit("ERROR(mem_pseudo_alloc_mpirealk): this feature is not available",-1)
      !in principle a simple call to MPI_WIN_CREATE should do the trick here
   else
      !if no communicator is given and simple option is chosen, do a normal alloc
      if(simp)then
         call mem_pseudo_alloc(A%d,n)
         A%c = c_null_ptr
         A%n = n
         A%w = 0
         A%t = 0
      else
         !if no communicator is given the default is to use MPI_alloc_mem with mpi
#ifdef VAR_MPI
         call mem_pseudo_alloc(A%d,n)
         A%c = c_loc(A%d(1))
         A%n = n
         A%w = 0
         A%t = 1
#else
         call mem_pseudo_alloc(A%d,n)
         A%c = c_null_ptr
         A%n = n
         A%w = 0
         A%t = 0
#endif
      endif
   endif
end subroutine mem_pseudo_alloc_mpirealk
subroutine mem_pseudo_dealloc_mpirealk(A)
   implicit none
   type(mpi_realk)                     :: A
   integer(kind=ls_mpik)               :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
#endif
   !check the allocation type. If 0 a normal allocation was used, normal
   !deallocation will be used, if 1 MPI_ALLOC_MEM was used and it will be deallocd 
   !accordingly,               if 2 MPI_ALLOC_WIN was used and freeing the window
   !will deallocate the pointer
   if(A%t==0)then
      call mem_pseudo_dealloc(A%d)
      A%n = 0
   else if(A%t==1)then
      call mem_pseudo_dealloc(A%d)
      A%c = c_null_ptr
      A%n = 0
   else if(A%t==2)then
      call lsquit("ERROR(mem_pseudo_dealloc_mpirealk): not yet implemented",-1)

   else
      call lsquit("ERROR(mem_pseudo_dealloc_mpirealk):wrong type",-1)
   endif

   !set back the default after deallocation
   A%t = -1
end subroutine mem_pseudo_dealloc_mpirealk

!----- ALLOCATE REAL POINTERS -----!

SUBROUTINE real_allocate_1dim(A,n)
   implicit none
   integer(kind=4),intent(in)  :: n
   REAL(REALK),pointer :: A(:)
   integer :: IERR
   integer (kind=long) :: nsize
   !$OMP CRITICAL
   nullify(A)
   ALLOCATE(A(n),STAT = IERR)
   !$OMP END CRITICAL
   nsize = size(A,KIND=long)*mem_realsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_1dim',IERR,n
      CALL memory_error_quit('Error in real_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_1dim

SUBROUTINE real_allocate_1dim_sp(A,n)  ! single precision
   implicit none
   integer,intent(in)  :: n
   REAL(4),pointer :: A(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n),STAT = IERR)
   nsize = size(A,KIND=long)*4
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_1dim_sp',IERR,n
      CALL memory_error_quit('Error in real_allocate_1dim_sp',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_1dim_sp


SUBROUTINE real_allocate_1dim_int64(A,n)
   implicit none
   integer(kind=long),intent(in)  :: n
   REAL(REALK),pointer :: A(:)
   integer :: IERR
   integer (kind=long) :: nsize
   real(realk) :: MemAvail
   logical :: memfound

   nullify(A)
   ALLOCATE(A(n),STAT = IERR)

   nsize = size(A,KIND=long)*mem_realsize

   IF (IERR.NE. 0) THEN
      MemAvail = 0.0E0_realk
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
      call get_available_memory(6,MemAvail,memfound)
#endif
      write(*,'("ERROR(real_allocate_1dim_int64),&
         & status=",I6," n=",I15," GBallocd=",f19.10," Sys free:",f19.10)')&
         & IERR,n,1.0E-9_realk*mem_allocated_global,MemAvail
      CALL memory_error_quit('Error in real_allocate_1dim_int64',nsize)
   ENDIF

   call mem_allocated_mem_real(nsize)

END SUBROUTINE real_allocate_1dim_int64

SUBROUTINE real_allocate_2dim(A,n1,n2)
   implicit none
   integer,intent(in)  :: n1, n2
   REAL(REALK),pointer :: A(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2),STAT = IERR)
   nsize = size(A,KIND=long)*mem_realsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_2dim',IERR,n1,n2
      call memory_error_quit('Error in real_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_2dim

SUBROUTINE real_allocate_2dim_sp(A,n1,n2)  ! single precision
   implicit none
   integer,intent(in)  :: n1, n2
   REAL(4),pointer :: A(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2),STAT = IERR)
   nsize = size(A,KIND=long)*4
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_2dim_sp',IERR,n1,n2
      call memory_error_quit('Error in real_allocate_2dim_sp',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_2dim_sp


SUBROUTINE real_allocate_2dim_zero(A,n1,n2,First,Second)
   ! Allocates 2d arrays starting from zero index
   ! for first,second or both dimensions
   implicit none
   integer,intent(in)  :: n1, n2
   REAL(REALK),pointer :: A(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   Logical :: First, Second
   !
   nullify(A)
   ! Both need zeroth element
   If (First .AND. Second) then
      ALLOCATE(A(0:n1,0:n2),STAT = IERR)
      IF (IERR.NE. 0) THEN
         write(*,*) 'Error1 in real_allocate_2dim_zero',IERR,n1,n2
         CALL MEMORY_ERROR_QUIT('Error1 in real_allocate_2dim_zero',nsize)
      ENDIF
   Else
      ! Only one
      If (First) then
         ALLOCATE(A(0:n1,n2),STAT = IERR)
         IF (IERR.NE. 0) THEN
            write(*,*) 'Error2 in real_allocate_2dim_zero',IERR,n1,n2
            CALL MEMORY_ERROR_QUIT('Error2 in real_allocate_2dim_zero',nsize)
         ENDIF
      Else
         If (Second) then
            ALLOCATE(A(n1,0:n2),STAT = IERR)
            IF (IERR.NE. 0) THEN
               write(*,*) 'Error2 in real_allocate_2dim_zero',IERR,n1,n2
               CALL MEMORY_ERROR_QUIT('Error2 in real_allocate_2dim_zero',nsize)
            ENDIF
         Else
            ! None :: an error, should be at least one.
            write(*,*) 'Error2 in real_allocate_2dim_zero'
            CALL MEMORY_ERROR_QUIT('Error2 in real_allocate_2dim_zero',nsize)
         Endif
      Endif
   Endif
   nsize = size(A,KIND=long)*mem_realsize
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_2dim_zero

SUBROUTINE real_allocate_3dim(A,n1,n2,n3)
   implicit none
   integer,intent(in)  :: n1, n2, n3
   REAL(REALK),pointer :: A(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2,n3),STAT = IERR)
   nsize = size(A,KIND=long)*mem_realsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_3dim',IERR,n1,n2,n3
      CALL MEMORY_ERROR_QUIT('Error in real_allocate_3dim',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_3dim

SUBROUTINE real_allocate_3dim_sp(A,n1,n2,n3)
   implicit none
   integer,intent(in)  :: n1, n2, n3
   REAL(4),pointer :: A(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2,n3),STAT = IERR)
   nsize = size(A,KIND=long)*4
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_3dim_sp',IERR,n1,n2,n3
      CALL MEMORY_ERROR_QUIT('Error in real_allocate_3dim_sp',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_3dim_sp


SUBROUTINE real_allocate_3dim_zero(A,n1,n2,n3,z1,z2,z3)
   implicit none
   integer,intent(in)  :: n1, n2, n3
   logical,intent(in)  :: z1,z2,z3
   REAL(REALK),pointer :: A(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   integer             :: i1,i2,i3
   nullify(A)
   i1=1
   IF(z1)i1=0
   i2=1
   IF(z2)i2=0
   i3=1
   IF(z3)i3=0
   ALLOCATE(A(i1:n1,i2:n2,i3:n3),STAT = IERR)
   nsize = size(A,KIND=long)*mem_realsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_3dim',IERR,n1,n2,n3
      CALL MEMORY_ERROR_QUIT('Error in real_allocate_3dim',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_3dim_zero

SUBROUTINE real_allocate_4dim(A,n1,n2,n3,n4)
   implicit none
   integer,intent(in)  :: n1,n2,n3,n4
   REAL(REALK),pointer :: A(:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2,n3,n4),STAT = IERR)
   nsize = size(A,KIND=long)*mem_realsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_4dim',IERR,n1,n2,n3,n4
      CALL MEMORY_ERROR_QUIT('Error in real_allocate_4dim',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_4dim

SUBROUTINE real_allocate_5dim(A,n1,n2,n3,n4,n5)
   implicit none
   integer,intent(in)  :: n1,n2,n3,n4,n5
   REAL(REALK),pointer :: A(:,:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2,n3,n4,n5),STAT = IERR)
   nsize = size(A,KIND=long)*mem_realsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_5dim',IERR,n1,n2,n3,n4,n5
      CALL MEMORY_ERROR_QUIT('Error in real_allocate_5dim',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_5dim

SUBROUTINE real_allocate_5dim_zero(A,n1,n2,n3,n4,n5,z1,z2,z3,z4,z5)
   ! Allocates 5d arrays starting from zero index
   ! for some or all dimensions
   implicit none
   integer,intent(in)  :: n1,n2,n3,n4,n5
   logical,intent(in)  :: z1,z2,z3,z4,z5
   REAL(REALK),pointer :: A(:,:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   integer             :: i1,i2,i3,i4,i5
   nullify(A)
   i1=1
   IF(z1)i1=0
   i2=1
   IF(z2)i2=0
   i3=1
   IF(z3)i3=0
   i4=1
   IF(z4)i4=0
   i5=1
   IF(z5)i5=0
   ALLOCATE(A(i1:n1,i2:n2,i3:n3,i4:n4,i5:n5),STAT = IERR)
   nsize = size(A,KIND=long)*mem_realsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_5dim_zero',IERR,n1,n2,n3,n4,n5
      CALL MEMORY_ERROR_QUIT('Error in real_allocate_5dim_zero',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_5dim_zero

SUBROUTINE real_allocate_7dim_zero(A,n1,n2,n3,n4,n5,n6,n7,z1,z2,z3,z4,z5,z6,z7)
   ! Allocates 5d arrays starting from zero index
   ! for some or all dimensions
   implicit none
   integer,intent(in)  :: n1,n2,n3,n4,n5,n6,n7
   logical,intent(in)  :: z1,z2,z3,z4,z5,z6,z7
   REAL(REALK),pointer :: A(:,:,:,:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   integer             :: i1,i2,i3,i4,i5,i6,i7
   nullify(A)
   i1=1
   IF(z1)i1=0
   i2=1
   IF(z2)i2=0
   i3=1
   IF(z3)i3=0
   i4=1
   IF(z4)i4=0
   i5=1
   IF(z5)i5=0
   i6=1
   IF(z6)i6=0
   i7=1
   IF(z7)i7=0
   ALLOCATE(A(i1:n1,i2:n2,i3:n3,i4:n4,i5:n5,i6:n6,i7:n7),STAT = IERR)
   nsize = size(A,KIND=long)*mem_realsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_allocate_7dim_zero',IERR,n1,n2,n3,n4,n5,n6,n7
      CALL MEMORY_ERROR_QUIT('Error in real_allocate_7dim_zero',nsize)
   ENDIF
   call mem_allocated_mem_real(nsize)
END SUBROUTINE real_allocate_7dim_zero

!----- DEALLOCATE REAL POINTERS -----!

SUBROUTINE real_deallocate_1dim(A)
   implicit none
   REAL(REALK),pointer :: A(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_1dim',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_1dim

SUBROUTINE real_deallocate_1dim_sp(A)
   implicit none
   REAL(4),pointer :: A(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*4
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_1dim_sp - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_1dim_sp',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_1dim_sp',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_1dim_sp


SUBROUTINE real_deallocate_2dim(A)
   implicit none
   REAL(REALK),pointer :: A(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_2dim - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_2dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_2dim',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_2dim

SUBROUTINE real_deallocate_2dim_sp(A)
   implicit none
   REAL(4),pointer :: A(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*4
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_2dim_sp - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_2dim_sp',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_2dim_sp',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_2dim_sp


SUBROUTINE real_deallocate_3dim(A)
   implicit none
   REAL(REALK),pointer :: A(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_3dim - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_3dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_3dim',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_3dim

SUBROUTINE real_deallocate_3dim_sp(A)
   implicit none
   REAL(4),pointer :: A(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*4
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_3dim_sp - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_3dim_sp',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_3dim_sp',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_3dim_sp

SUBROUTINE real_deallocate_4dim(A)
   implicit none
   REAL(REALK),pointer :: A(:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_4dim - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_4dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_4dim',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_4dim

SUBROUTINE real_deallocate_5dim(A)
   implicit none
   REAL(REALK),pointer :: A(:,:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_5dim - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_5dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_5dim',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_5dim

SUBROUTINE real_deallocate_7dim(A)
   implicit none
   REAL(REALK),pointer :: A(:,:,:,:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in real_deallocate_5dim - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in real_deallocate_7dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in real_deallocate_7dim',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_7dim
SUBROUTINE lsmpi_allocate_dV8(A,cip,n) 
   implicit none
   integer(kind=8),intent(in)  :: n
   real(realk),pointer :: A(:)
   type(c_ptr), intent(inout) :: cip
   integer (kind=ls_mpik) :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
   nullify(A)
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,lsmpi_len_realk,IERR)

   bytes =n*lsmpi_len_realk
   nsize = bytes

   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_allocate_dV8):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

   if(bytes<0)then
      print *,"calling MPI_ALLOC_MEM with",bytes,n,lsmpi_len_realk
      call lsquit('ERROR(mpi_allocate_dV8):something wrong in the subroutine, bytes<0',-1)
   endif
   call MPI_ALLOC_MEM(bytes,info,cip,IERR)

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_allocate_dV8):error in alloc",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   call c_f_pointer(cip,A,[n])

   call mem_allocated_mem_mpi(nsize)

#else
   call lsquit("ERROR(mpi_allocate_dV8):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_allocate_dV8
!allocate MPI memory
SUBROUTINE lsmpi_allocate_dV4(A,cip,n) 
   implicit none
   integer(kind=4),intent(in)  :: n
   real(realk),pointer :: A(:)
   type(c_ptr), intent(inout) :: cip
   integer (kind=ls_mpik) :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
   nullify(A)
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,lsmpi_len_realk,IERR)

   bytes =n*lsmpi_len_realk
   nsize = bytes

   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_allocate_dV4):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif


   if(bytes<0)then
      print *,"calling MPI_ALLOC_MEM with",bytes,n,lsmpi_len_realk
      call lsquit('ERROR(mpi_allocate_dV4):something wrong in the subroutine, bytes<0',-1)
   endif

   call MPI_ALLOC_MEM(bytes,info,cip,IERR)

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_allocate_dV4):error in alloc",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   call c_f_pointer(cip,A,[n])

   call mem_allocated_mem_mpi(nsize)

#else
   call lsquit("ERROR(mpi_allocate_dV4):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_allocate_dV4
SUBROUTINE lsmpi_local_allocate_I4V4(A,cip,n1,win,comm,n2) 
   implicit none
   integer(kind=4),intent(in)          :: n1
   integer(kind=4),intent(in),optional :: n2
   integer(kind=4),pointer             :: A(:)
   integer(kind=ls_mpik),intent(in)    ::comm
   integer(kind=ls_mpik),intent(inout) :: win
   type(c_ptr), intent(inout)          :: cip
   integer (kind=ls_mpik)              :: IERR,info
   integer (kind=long)                 :: nsize
   character(120)                      :: errmsg
   integer(kind=8)                     :: assoc
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_int4,lb,bytes
   nullify(A)
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_INTEGER4,lb,lsmpi_int4,IERR)
   assoc = n1
   if(present(n2))assoc=n2
   nsize = assoc*lsmpi_int4
   bytes =n1*lsmpi_int4

   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_local_allocate_dV8):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

#ifdef VAR_HAVE_MPI3
   call MPI_WIN_ALLOCATE_SHARED(bytes,lsmpi_int4,info,comm,cip,win,IERR)
#else
   call lsquit("ERROR(mpi_local_allocate_I4V4): not possible withot mpi3",-1)
#endif

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_local_allocate_I4V4):error in alloc",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF
   call c_f_pointer(cip,A,[assoc])

   call mem_allocated_mem_mpi(nsize)

#else
   call lsquit("ERROR(mpi_allocate_I4V4):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_local_allocate_I4V4
SUBROUTINE lsmpi_local_allocate_I8V8(A,cip,n1,win,comm,n2) 
   implicit none
   integer(kind=8),intent(in)          :: n1
   integer(kind=8),intent(in),optional :: n2
   integer(kind=8),pointer             :: A(:)
   integer(kind=ls_mpik),intent(in)    ::comm
   integer(kind=ls_mpik),intent(inout) :: win
   type(c_ptr), intent(inout)          :: cip
   integer (kind=ls_mpik)              :: IERR,info
   integer (kind=long)                 :: nsize
   character(120)                      :: errmsg
   integer(kind=8)                     :: assoc
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_int8,lb,bytes
   nullify(A)
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_INTEGER8,lb,lsmpi_int8,IERR)

   assoc = n1
   if(present(n2))assoc=n2
   nsize = assoc*lsmpi_int8
   bytes =n1*lsmpi_int8

   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_local_allocate_dV8):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

#ifdef VAR_HAVE_MPI3
   call MPI_WIN_ALLOCATE_SHARED(bytes,lsmpi_int8,info,comm,cip,win,IERR)
#else
   call lsquit("ERROR(mpi_local_allocate_I8V8): not possible withot mpi3",-1)
#endif

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_local_allocate_I8V8):error in alloc",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF
   call c_f_pointer(cip,A,[assoc])

   call mem_allocated_mem_mpi(nsize)

#else
   call lsquit("ERROR(mpi_allocate_I8V8):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_local_allocate_I8V8
SUBROUTINE lsmpi_allocate_d(A,n1,comm,local,simple) 
   implicit none
   integer(kind=8),intent(in)          :: n1
   logical,intent(in),optional         :: local,simple
   type(mpi_realk)                     :: A
   integer(kind=ls_mpik),optional,intent(in)    ::comm
   integer (kind=ls_mpik)              :: IERR,info
   integer (kind=long)                 :: nsize
   character(120)                      :: errmsg
   logical                             :: loc,simp
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
#endif

   simp = .false.
   if(present(simple))simp = simple

   if(present(comm).and..not.simp)then
#ifdef VAR_MPI
      loc = .false.
      if(present(local))loc=local

      nullify(A%d)
      info = MPI_INFO_NULL

      call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,lsmpi_len_realk,IERR)

      nsize = n1 * lsmpi_len_realk

      if(IERR/=0)then
         write (errmsg,'("ERROR(mpi_local_allocate_dV8):error in&
            & lsmpi_type_get_extent",I5)') IERR
         call memory_error_quit(errmsg,nsize)
      endif

#ifdef VAR_HAVE_MPI3
      if(loc) then
         bytes = int(0,kind=MPI_ADDRESS_KIND)
         if( infpar%pc_mynum == infpar%pc_nodtot - 1 ) bytes = n1 * lsmpi_len_realk

         if(bytes<0)then
            print *,"calling MPI_WIN_ALLOCATE with",bytes,n1,lsmpi_len_realk
            call lsquit('ERROR(lsmpi_allocate_d):something wrong in the subroutine, bytes<0',-1)
         endif

         call MPI_WIN_ALLOCATE_SHARED(bytes,lsmpi_len_realk,info,comm,A%c,A%w,IERR)

      else
         bytes = n1 * lsmpi_len_realk

         if(bytes<0)then
            print *,"calling MPI_WIN_ALLOCATE with",bytes,n1,lsmpi_len_realk
            call lsquit('ERROR(lsmpi_allocate_d):something wrong in the subroutine, bytes<0',-1)
         endif

         call MPI_WIN_ALLOCATE(bytes,lsmpi_len_realk,info,comm,A%c,A%w,IERR)
      endif
#else
      call lsquit("ERROR(lsmpi_allocate_d): not possible withot mpi3",-1)
#endif

      IF (IERR.NE. 0) THEN
         write (errmsg,'("ERROR(lsmpi_allocate_d):error in alloc",I5)') IERR
         CALL memory_error_quit(errmsg,nsize)
      ENDIF

      call c_f_pointer(A%c,A%d,[n1])
      call mem_allocated_mem_mpi(nsize)
      A%n = n1

      A%t = 2
#else
      call lsquit("ERROR(mpi_allocate_d):compiled without MPI, this is not&
         &available",-1)
#endif
   else
      !if no communicator is given and simple option is chosen, do a normal alloc
      if(simp)then
         call mem_alloc(A%d,n1)
         A%c = c_null_ptr
         A%n = n1
         A%w = 0
         A%t = 0
      else
         !if no communicator is given the default is to use MPI_alloc_mem with mpi
#ifdef VAR_MPI
         call mem_alloc(A%d,A%c,n1)
         A%n = n1
         A%w = 0
         A%t = 1
#else
         call mem_alloc(A%d,n1)
         A%c = c_null_ptr
         A%n = n1
         A%w = 0
         A%t = 0
#endif
      endif
   endif
END SUBROUTINE lsmpi_allocate_d
SUBROUTINE lsmpi_local_allocate_dV8(A,cip,n1,win,comm,n2) 
   implicit none
   integer(kind=8),intent(in)          :: n1
   integer(kind=8),intent(in),optional :: n2
   real(realk),pointer                 :: A(:)
   integer(kind=ls_mpik),intent(in)    ::comm
   integer(kind=ls_mpik),intent(inout) :: win
   type(c_ptr), intent(inout)          :: cip
   integer (kind=ls_mpik)              :: IERR,info
   integer (kind=long)                 :: nsize
   character(120)                      :: errmsg
   integer(kind=8)                     :: assoc
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
   nullify(A)
   ierr = 0
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,lsmpi_len_realk,IERR)

   bytes =n1*lsmpi_len_realk
   assoc = n1
   if(present(n2))assoc=n2
   nsize = assoc*lsmpi_len_realk
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_local_allocate_dV8):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

#ifdef VAR_HAVE_MPI3
   call MPI_WIN_ALLOCATE_SHARED(bytes,lsmpi_len_realk,info,comm,cip,win,IERR)
#else
   call lsquit("ERROR(mpi_local_allocate_dV8): not possible withot mpi3",-1)
#endif

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_local_allocate_dV8):error in alloc",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF
   call c_f_pointer(cip,A,[assoc])

   call mem_allocated_mem_mpi(nsize)

#else
   call lsquit("ERROR(mpi_allocate_dV8):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_local_allocate_dV8
SUBROUTINE lsmpi_local_allocate_dV4(A,cip,n1,win,comm,n2) 
   implicit none
   integer(kind=4),intent(in)          :: n1
   integer(kind=4),intent(in),optional :: n2
   real(realk),pointer                 :: A(:)
   integer(kind=ls_mpik),intent(in)    ::comm
   integer(kind=ls_mpik),intent(inout) :: win
   type(c_ptr), intent(inout)          :: cip
   integer (kind=ls_mpik)              :: IERR,info
   integer (kind=long)                 :: nsize
   integer(kind=4)                     :: assoc
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
   nullify(A)
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,lsmpi_len_realk,IERR)

   bytes =n1*lsmpi_len_realk
   assoc = n1
   if(present(n2))assoc=n2
   nsize = assoc*lsmpi_len_realk
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_local_allocate_dV4):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

#ifdef VAR_HAVE_MPI3
   call MPI_WIN_ALLOCATE_SHARED(bytes,lsmpi_len_realk,info,comm,cip,win,IERR)
#else
   call lsquit("ERROR(mpi_local_allocate_dV4): not possible withot mpi3",-1)
#endif

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_local_allocate_dV4):error in alloc",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   call c_f_pointer(cip,A,[assoc])

   call mem_allocated_mem_mpi(nsize)

#else
   call lsquit("ERROR(mpi_local_allocate_dV4):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_local_allocate_dV4

SUBROUTINE lsmpi_allocate_i8V(A,cip,n)  ! single precision
   implicit none
   integer,intent(in)  :: n
   integer(kind=8),pointer :: A(:)
   type(c_ptr), intent(inout) :: cip
   integer(kind=ls_mpik) :: IERR
#ifdef VAR_HAVE_MPI3
   !type(MPI_Info) :: info
   integer(kind=ls_mpik) :: info
#else
   integer(kind=ls_mpik) :: info
#endif
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_intlen,lb,bytes
   nullify(A)

   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_INTEGER8,lb,lsmpi_intlen,IERR)
   bytes =n*lsmpi_intlen
   nsize = bytes
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_allocate_i8V):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif


   if(bytes<0)then
      print *,"calling MPI_ALLOC_MEM with",bytes,n,lsmpi_intlen
      call lsquit('ERROR(mpi_allocate_i8V):something wrong in the subroutine, bytes<0',-1)
   endif

   call MPI_ALLOC_MEM(bytes,info,cip,IERR)

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_allocate_i8V):error in alloc",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF
   call c_f_pointer(cip,A,[n])

   call mem_allocated_mem_mpi(nsize)

#else
   call lsquit("ERROR(mpi_allocate_i8V):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_allocate_i8V
SUBROUTINE lsmpi_allocate_i4V(A,cip,n)  ! single precision
   implicit none
   integer,intent(in)  :: n
   integer(kind=4),pointer :: A(:)
   type(c_ptr), intent(inout) :: cip
   integer(kind=ls_mpik) :: IERR
#ifdef VAR_HAVE_MPI3
   !type(MPI_Info) :: info
   integer(kind=ls_mpik) :: info
#else
   integer(kind=ls_mpik) :: info
#endif
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_intlen,lb,bytes
   nullify(A)

   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_INTEGER4,lb,lsmpi_intlen,IERR)
   bytes =n*lsmpi_intlen
   nsize = bytes
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_allocate_i4V):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif


   if(bytes<0)then
      print *,"calling MPI_ALLOC_MEM with",bytes,n,lsmpi_intlen
      call lsquit('ERROR(mpi_allocate_i4V):something wrong in the subroutine, bytes<0',-1)
   endif

   call MPI_ALLOC_MEM(bytes,info,cip,IERR)

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_allocate_i4V):error in alloc",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF
   call c_f_pointer(cip,A,[n])

   call mem_allocated_mem_mpi(nsize)

#else
   call lsquit("ERROR(mpi_allocate_i4V):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_allocate_i4V

!deallcate MPI memory
SUBROUTINE lsmpi_deallocate_i8V(A,cip)
   implicit none
   integer(kind=8),pointer :: A(:)
   type(c_ptr), intent(inout) :: cip
   integer(kind=ls_mpik) :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_intlen,lb,bytes
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_INTEGER8,lb,lsmpi_intlen,IERR)
   nsize = size(A)*lsmpi_intlen
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_deallocate_i8V):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

   call mem_deallocated_mem_mpi(nsize)

   if (.not.ASSOCIATED(A).or..not.c_associated(cip)) then
      print *,'Memory previously released!!'
      call memory_error_quit('ERROR(mpi_deallocate_i8V): memory previously released',nsize)
   endif

   call MPI_FREE_MEM(A,IERR)
   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_deallocate_i8V):error in MPI_FREE_MEM",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   nullify(A)
   cip = c_null_ptr
#else
   call lsquit("ERROR(mpi_deallocate_i8V):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_deallocate_i8V
!deallcate MPI memory
SUBROUTINE lsmpi_deallocate_i4V(A,cip)
   implicit none
   integer(kind=4),pointer :: A(:)
   type(c_ptr), intent(inout) :: cip
   integer(kind=ls_mpik) :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_intlen,lb,bytes
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_INTEGER4,lb,lsmpi_intlen,IERR)
   nsize = size(A)*lsmpi_intlen
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_deallocate_i4V):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

   call mem_deallocated_mem_mpi(nsize)

   if (.not.ASSOCIATED(A).or..not.c_associated(cip)) then
      print *,'Memory previously released!!'
      call memory_error_quit('ERROR(mpi_deallocate_i4V): memory previously released',nsize)
   endif

   call MPI_FREE_MEM(A,IERR)
   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_deallocate_i4V):error in MPI_FREE_MEM",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   nullify(A)
   cip = c_null_ptr
#else
   call lsquit("ERROR(mpi_deallocate_i4V):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_deallocate_i4V

SUBROUTINE lsmpi_deallocate_dV(A,cip)
   implicit none
   real(realk),pointer :: A(:)
   type(c_ptr), intent(inout) :: cip
   integer(kind=ls_mpik) :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,lsmpi_len_realk,IERR)
   nsize = size(A)*lsmpi_len_realk
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_deallocate_dV):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

   call mem_deallocated_mem_mpi(nsize)

   if (.not.ASSOCIATED(A).or..not.c_associated(cip)) then
      print *,'Memory previously released!!'
      call memory_error_quit('ERROR(mpi_deallocate_dV): memory previously released',nsize)
   endif

   call MPI_FREE_MEM(A,IERR)

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_allocate_dV):error in MPI_FREE_MEM",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   nullify(A)
   cip = c_null_ptr
#else
   call lsquit("ERROR(mpi_deallocate_dV):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_deallocate_dV
SUBROUTINE lsmpi_deallocate_d(A)
   implicit none
   type(mpi_realk)                     :: A
   integer(kind=ls_mpik)               :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
#endif
   !check the allocation type. If 0 a normal allocation was used, normal
   !deallocation will be used, if 1 MPI_ALLOC_MEM was used and it will be deallocd 
   !accordingly,               if 2 MPI_ALLOC_WIN was used and freeing the window
   !will deallocate the pointer
   if(A%t==0)then
      call mem_dealloc(A%d)
      A%n = 0
      elseif(A%t==1)then
      call mem_dealloc(A%d,A%c)
      A%n = 0
      elseif(A%t==2)then
      if(c_associated(A%c).and.associated(A%d))then
#ifdef VAR_MPI
         info = MPI_INFO_NULL

         call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,lsmpi_len_realk,IERR)
         nsize = size(A%d)*lsmpi_len_realk

         if(IERR/=0)then
            write (errmsg,'("ERROR(mpi_deallocate_d):error in lsmpi_type_get_extent",I5)') IERR
            call memory_error_quit(errmsg,nsize)
         endif

         call mem_deallocated_mem_mpi(nsize)

         if (.not.associated(A%d)) then
            print *,'Memory previously released!!'
            call memory_error_quit('ERROR(mpi_deallocate_d): memory previously released',nsize)
         endif

         call MPI_WIN_FREE(A%w,IERR)

         IF (IERR.NE. 0) THEN
            write (errmsg,'("ERROR(mpi_deallocate_d):error in MPI_FREE_MEM",I5)') IERR
            CALL memory_error_quit(errmsg,nsize)
         ENDIF

         nullify(A%d)
         A%c = c_null_ptr

#else
         call lsquit("ERROR(mpi_deallocate_d):compiled without MPI, this is not&
            &available",-1)
#endif

      else
         call lsquit("ERROR(mpi_deallocate_d):pointer not allocated",-1)
      endif

   else
      call lsquit("ERROR(mpi_deallocate_d):wrong type",-1)
   endif

   !set back the default after deallocation
   A%t = -1
END SUBROUTINE lsmpi_deallocate_d

SUBROUTINE lsmpi_local_deallocate_dV(A,cip,win)
   implicit none
   real(realk),pointer                 :: A(:)
   type(c_ptr), intent(inout)          :: cip
   integer(kind=ls_mpik),intent(inout) :: win
   integer(kind=ls_mpik)               :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_len_realk,lb,bytes
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,lsmpi_len_realk,IERR)
   nsize = size(A)*lsmpi_len_realk
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_deallocate_dV):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

   call mem_deallocated_mem_mpi(nsize)

   if (.not.ASSOCIATED(A).or..not.c_associated(cip)) then
      print *,'Memory previously released!!'
      call memory_error_quit('ERROR(mpi_deallocate_dV): memory previously released',nsize)
   endif

   call MPI_WIN_FREE(win,IERR)

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_local_allocate_dV):error in MPI_FREE_MEM",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   nullify(A)
   cip = c_null_ptr
#else
   call lsquit("ERROR(mpi_deallocate_dV):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_local_deallocate_dV
SUBROUTINE lsmpi_local_deallocate_I8V(A,cip,win)
   implicit none
   integer(kind=8),pointer             :: A(:)
   type(c_ptr), intent(inout)          :: cip
   integer(kind=ls_mpik),intent(inout) :: win
   integer(kind=ls_mpik)               :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_int8,lb,bytes
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_INTEGER8,lb,lsmpi_int8,IERR)
   nsize = size(A)*lsmpi_int8
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_deallocate_I8V):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

   call mem_deallocated_mem_mpi(nsize)

   if (.not.ASSOCIATED(A).or..not.c_associated(cip)) then
      print *,'Memory previously released!!'
      call memory_error_quit('ERROR(mpi_deallocate_I8V): memory previously released',nsize)
   endif

   call MPI_WIN_FREE(win,IERR)

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_local_allocate_I8V):error in MPI_FREE_MEM",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   nullify(A)
   cip = c_null_ptr
#else
   call lsquit("ERROR(mpi_deallocate_I8V):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_local_deallocate_I8V
SUBROUTINE lsmpi_local_deallocate_I4V(A,cip,win)
   implicit none
   integer(kind=4),pointer             :: A(:)
   type(c_ptr), intent(inout)          :: cip
   integer(kind=ls_mpik),intent(inout) :: win
   integer(kind=ls_mpik)               :: IERR,info
   integer (kind=long) :: nsize
   character(120) :: errmsg
#ifdef VAR_MPI
   integer(kind=MPI_ADDRESS_KIND) :: lsmpi_int4,lb,bytes
   info = MPI_INFO_NULL

   call MPI_TYPE_GET_EXTENT(MPI_INTEGER4,lb,lsmpi_int4,IERR)
   nsize = size(A)*lsmpi_int4
   if(IERR/=0)then
      write (errmsg,'("ERROR(mpi_deallocate_I4V):error in&
         & lsmpi_type_get_extent",I5)') IERR
      call memory_error_quit(errmsg,nsize)
   endif

   call mem_deallocated_mem_mpi(nsize)

   if (.not.ASSOCIATED(A).or..not.c_associated(cip)) then
      print *,'Memory previously released!!'
      call memory_error_quit('ERROR(mpi_deallocate_I4V): memory previously released',nsize)
   endif

   call MPI_WIN_FREE(win,IERR)

   IF (IERR.NE. 0) THEN
      write (errmsg,'("ERROR(mpi_local_allocate_I4V):error in MPI_FREE_MEM",I5)') IERR
      CALL memory_error_quit(errmsg,nsize)
   ENDIF

   nullify(A)
   cip = c_null_ptr
#else
   call lsquit("ERROR(mpi_deallocate_I4V):compiled without MPI, this is not&
      &available",-1)
#endif
END SUBROUTINE lsmpi_local_deallocate_I4V


!ALlocate complex
SUBROUTINE complex_allocate_1dim(A,n)
   implicit none
   integer,intent(in)  :: n
   complex(COMPLEXK),pointer :: A(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n),STAT = IERR)
   nsize = size(A,KIND=long)*mem_complexsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in complex_allocate_1dim',IERR,n
      CALL memory_error_quit('Error in complex_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_complex(nsize)
END SUBROUTINE complex_allocate_1dim

SUBROUTINE complex_allocate_2dim(A,n1,n2)
   implicit none
   integer,intent(in)  :: n1, n2
   COMPLEX(COMPLEXK),pointer :: A(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(A)
   ALLOCATE(A(n1,n2),STAT = IERR)
   nsize = size(A,KIND=long)*mem_complexsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in complex_allocate_2dim',IERR,n1,n2
      CALL MEMORY_ERROR_QUIT('Error in complex_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_complex(nsize)
END SUBROUTINE complex_allocate_2dim

!----- DEALLOCATE COMPLEX POINTERS -----!

SUBROUTINE complex_deallocate_1dim(A)
   implicit none
   complex(COMPLEXK),pointer :: A(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_complexsize
   call mem_deallocated_mem_complex(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in complex_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in complex_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in complex_deallocate_1dim',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE complex_deallocate_1dim

SUBROUTINE complex_deallocate_2dim(A)
   implicit none
   complex(COMPLEXK),pointer :: A(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_complexsize
   call mem_deallocated_mem_complex(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in complex_deallocate_2dim - memory previously released',nsize)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in complex_deallocate_2dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in complex_deallocate_2dim',nsize)
   ENDIF
   nullify(A)
END SUBROUTINE complex_deallocate_2dim

!----- ALLOCATE INTEGER POINTERS -----!
SUBROUTINE int8_allocate_1dim_wrapper4(I,n)
   implicit none
   integer(kind=4),intent(in)  :: n
   INTEGER(kind=8),pointer     :: I(:)
   integer(kind=8)             :: n8
   n8=n
   call int8_allocate_1dim(I,n8)
END SUBROUTINE int8_allocate_1dim_wrapper4
SUBROUTINE int8_allocate_1dim(I,n)
   implicit none
   integer(kind=8),intent(in)  :: n
   INTEGER(kind=8),pointer     :: I(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n),STAT = IERR)
   nsize = size(I,KIND=long)*mem_int8size
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int8_allocate_1dim',IERR,n
      call memory_error_quit('Error in int8_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int8_allocate_1dim

SUBROUTINE int4_allocate_1dim_wrapper4(I,n)
   implicit none
   integer(kind=4),intent(in)  :: n
   INTEGER(kind=4),pointer     :: I(:)
   integer(kind=8)             :: n8
   n8=n
   call int4_allocate_1dim(I,n8)
END SUBROUTINE int4_allocate_1dim_wrapper4
SUBROUTINE int4_allocate_1dim(I,n)
   implicit none
   integer(kind=8),intent(in)  :: n
   INTEGER(kind=4),pointer     :: I(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n),STAT = IERR)
   nsize = size(I,KIND=long)*mem_int4size
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int4_allocate_1dim',IERR,n
      call memory_error_quit('Error in int4_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int4_allocate_1dim



SUBROUTINE intS_allocate_1dim_wrapper4(I,n)
   implicit none
   integer(kind=4),intent(in)  :: n
   INTEGER(kind=short),pointer :: I(:)
   integer(kind=8)             :: n8
   n8=n
   call intS_allocate_1dim(I,N8)
END SUBROUTINE intS_allocate_1dim_wrapper4
SUBROUTINE intS_allocate_1dim(I,n)
   implicit none
   integer(kind=8),intent(in)  :: n
   INTEGER(kind=short),pointer     :: I(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n),STAT = IERR)
   nsize = size(I,KIND=long)*mem_shortintsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int32_allocateS_1dim',IERR,n
      call memory_error_quit('Error in int32_allocateS_1dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE intS_allocate_1dim



SUBROUTINE int4_allocate_2dim_wrapper4(I,n1,n2)
   implicit none
   integer(kind=4),intent(in) :: n1,n2
   INTEGER(kind=4),pointer    :: I(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   integer(kind=8) :: n18, n28
   n18=n1;n28=n2
   call int4_allocate_2dim(I,n18,n28)
END SUBROUTINE int4_allocate_2dim_wrapper4
SUBROUTINE int4_allocate_2dim(I,n1,n2)
   implicit none
   integer(kind=8),intent(in) :: n1,n2
   INTEGER(kind=4),pointer    :: I(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n1,n2),STAT = IERR)
   nsize = size(I,KIND=long)*mem_int4size
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int4_allocate_2dim',IERR,n1,n2
      call memory_error_quit('Error in int4_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int4_allocate_2dim

SUBROUTINE int8_allocate_2dim_wrapper4(I,n1,n2)
   implicit none
   integer(kind=4),intent(in) :: n1,n2
   INTEGER(kind=8),pointer    :: I(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   integer(kind=8) :: n18, n28
   n18=n1;n28=n2
   call int8_allocate_2dim(I,n18,n28)
END SUBROUTINE int8_allocate_2dim_wrapper4
SUBROUTINE int8_allocate_2dim(I,n1,n2)
   implicit none
   integer(kind=8),intent(in) :: n1,n2
   INTEGER(kind=8),pointer    :: I(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n1,n2),STAT = IERR)
   nsize = size(I,KIND=long)*mem_int8size
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int8_allocate_2dim',IERR,n1,n2
      call memory_error_quit('Error in int8_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int8_allocate_2dim

SUBROUTINE shortint_allocate_2dim(I,n1,n2)
   implicit none
   integer,intent(in) :: n1,n2
   INTEGER(kind=short),pointer    :: I(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n1,n2),STAT = IERR)
   nsize = size(I,KIND=long)*mem_shortintsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in shortint_allocate_2dim',IERR,n1,n2
      call memory_error_quit('Error in shortint_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE shortint_allocate_2dim

SUBROUTINE int_allocate_3dim(I,n1,n2,n3)
   implicit none
   integer,intent(in) :: n1,n2,n3
   INTEGER,pointer    :: I(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n1,n2,n3),STAT = IERR)
   nsize = size(I,KIND=long)*mem_intsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int_allocate_3dim',IERR,n1,n2,n3
      call memory_error_quit('Error in int_allocate_3dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_3dim



SUBROUTINE int8_allocate_4dim_wrapper4(I,n1,n2,n3,n4)
   implicit none
   integer(kind=4),intent(in) :: n1,n2,n3,n4
   INTEGER(kind=8),pointer    :: I(:,:,:,:)
   integer(kind=8)            :: n18,n28,n38,n48
   n18=int(n1,kind=8)
   n28=int(n2,kind=8)
   n38=int(n3,kind=8)
   n48=int(n4,kind=8)
   call int8_allocate_4dim(I,n18,n28,n38,n48)
END SUBROUTINE int8_allocate_4dim_wrapper4
SUBROUTINE int8_allocate_4dim(I,n1,n2,n3,n4)
   implicit none
   integer(kind=8),intent(in) :: n1,n2,n3,n4
   INTEGER(kind=8),pointer    :: I(:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n1,n2,n3,n4),STAT = IERR)
   nsize = size(I,KIND=long)*mem_int8size
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int8_allocate_4dim',IERR,n1,n2,n3,n4
      call memory_error_quit('Error in int8_allocate_4dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int8_allocate_4dim

SUBROUTINE int4_allocate_4dim_wrapper4(I,n1,n2,n3,n4)
   implicit none
   integer(kind=4),intent(in) :: n1,n2,n3,n4
   INTEGER(kind=4),pointer    :: I(:,:,:,:)
   integer(kind=8)            :: n18,n28,n38,n48
   n18=n1;n28=n2;n38=n3;n48=n4
   call int4_allocate_4dim(I,n18,n28,n38,n48)
END SUBROUTINE int4_allocate_4dim_wrapper4
SUBROUTINE int4_allocate_4dim(I,n1,n2,n3,n4)
   implicit none
   integer(kind=8),intent(in) :: n1,n2,n3,n4
   INTEGER(kind=4),pointer    :: I(:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(I)
   ALLOCATE(I(n1,n2,n3,n4),STAT = IERR)
   nsize = size(I,KIND=long)*mem_int4size
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int4_allocate_4dim',IERR,n1,n2,n3,n4
      call memory_error_quit('Error in int4_allocate_4dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int4_allocate_4dim



SUBROUTINE int_allocate_1dim_zero(I,n,z1)
   implicit none
   integer,intent(in)  :: n
   INTEGER,pointer     :: I(:)
   logical,intent(in) :: z1
   integer :: IERR,i1
   integer (kind=long) :: nsize
   i1=1
   IF(z1)i1=0
   nullify(I)
   ALLOCATE(I(i1:n),STAT = IERR)
   nsize = size(I,KIND=long)*mem_intsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int_allocate_1dim',IERR,n
      call memory_error_quit('Error in int_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_1dim_zero

SUBROUTINE int_allocate_2dim_zero(I,n1,n2,z1,z2)
   implicit none
   integer,intent(in) :: n1,n2
   INTEGER,pointer    :: I(:,:)
   logical,intent(in) :: z1,z2
   integer :: IERR,i1,i2
   integer (kind=long) :: nsize
   i1=1
   IF(z1)i1=0
   i2=1
   IF(z2)i2=0
   nullify(I)
   ALLOCATE(I(i1:n1,i2:n2),STAT = IERR)
   nsize = size(I,KIND=long)*mem_intsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int_allocate_2dim',IERR,n1,n2
      call memory_error_quit('Error in int_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_2dim_zero

SUBROUTINE int_allocate_3dim_zero(I,n1,n2,n3,z1,z2,z3)
   implicit none
   integer,intent(in) :: n1,n2,n3
   logical,intent(in) :: z1,z2,z3
   INTEGER,pointer    :: I(:,:,:)
   integer :: IERR,i1,i2,i3
   integer (kind=long) :: nsize
   i1=1
   IF(z1)i1=0
   i2=1
   IF(z2)i2=0
   i3=1
   IF(z3)i3=0
   nullify(I)
   ALLOCATE(I(i1:n1,i2:n2,i3:n3),STAT = IERR)
   nsize = size(I,KIND=long)*mem_intsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int_allocate_3dim',IERR,n1,n2,n3
      call memory_error_quit('Error in int_allocate_3dim',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_3dim_zero

SUBROUTINE int_allocate_4dim_zero(I,n1,n2,n3,n4,z1,z2,z3,z4)
   implicit none
   integer,intent(in) :: n1,n2,n3,n4
   logical,intent(in) :: z1,z2,z3,z4
   INTEGER,pointer    :: I(:,:,:,:)
   integer :: IERR,i1,i2,i3,i4
   integer (kind=long) :: nsize
   i1=1
   IF(z1)i1=0
   i2=1
   IF(z2)i2=0
   i3=1
   IF(z2)i3=0
   i4=1
   IF(z2)i4=0
   nullify(I)
   ALLOCATE(I(i1:n1,i2:n2,i3:n3,i4:n4),STAT = IERR)
   nsize = size(I,KIND=long)*mem_intsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int_allocate_4dim_zero',IERR,n1,n2,n3,n4
      call memory_error_quit('Error in int_allocate_4dim_zero',nsize)
   ENDIF
   call mem_allocated_mem_integer(nsize)
END SUBROUTINE int_allocate_4dim_zero

!----- DEALLOCATE INTEGER POINTERS -----!
SUBROUTINE int8_deallocate_1dim(I)
   implicit none
   INTEGER(kind=8),pointer :: I(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_int8size
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in int8_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int8_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in int8_deallocate_1dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE int8_deallocate_1dim

SUBROUTINE int4_deallocate_1dim(I)
   implicit none
   INTEGER(kind=4),pointer :: I(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_int4size
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in int4_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int4_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in int4_deallocate_1dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE int4_deallocate_1dim

SUBROUTINE int4_deallocate_2dim(I)
   implicit none
   INTEGER(kind=4),pointer :: I(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_int4size
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in int4_deallocate_2dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int4_deallocate_2dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in int4_deallocate_2dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE int4_deallocate_2dim
SUBROUTINE int8_deallocate_2dim(I)
   implicit none
   INTEGER(kind=8),pointer :: I(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_int8size
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in int8_deallocate_2dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int8_deallocate_2dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in int8_deallocate_2dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE int8_deallocate_2dim

SUBROUTINE shortint_deallocate_1dim(I)
   implicit none
   INTEGER(kind=short),pointer :: I(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_shortintsize
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in shortint_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in shortint_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in shortint_deallocate_1dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE shortint_deallocate_1dim

SUBROUTINE shortint_deallocate_2dim(I)
   implicit none
   INTEGER(kind=short),pointer :: I(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_shortintsize
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in shortint_deallocate_2dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in shortint_deallocate_2dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in shortint_deallocate_2dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE shortint_deallocate_2dim

SUBROUTINE int4_deallocate_3dim(I)
   implicit none
   INTEGER(kind=4),pointer :: I(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_int4size
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in int4_deallocate_3dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int4_deallocate_3dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in int4_deallocate_3dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE int4_deallocate_3dim
SUBROUTINE int8_deallocate_3dim(I)
   implicit none
   INTEGER(kind=8),pointer :: I(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_int8size
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in int8_deallocate_3dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int8_deallocate_3dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in int8_deallocate_3dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE int8_deallocate_3dim

SUBROUTINE int4_deallocate_4dim(I)
   implicit none
   INTEGER(kind=4),pointer :: I(:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_int4size
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in int4_deallocate_4dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int4_deallocate_4dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in int4_deallocate_4dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE int4_deallocate_4dim
SUBROUTINE int8_deallocate_4dim(I)
   implicit none
   INTEGER(kind=8),pointer :: I(:,:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_int8size
   call mem_deallocated_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in int8_deallocate_4dim - memory previously released',nsize)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int8_deallocate_4dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in int8_deallocate_4dim',nsize)
   ENDIF
   nullify(I)
END SUBROUTINE int8_deallocate_4dim


!----- ALLOCATE CHARACTER POINTERS -----!

SUBROUTINE char_allocate_1dim_wrapper4(C,n)
   implicit none
   integer(kind=4),intent(in):: n
   CHARACTER(LEN=*),pointer  :: C(:)
   integer(kind=8)::n8
   n8=n
   call char_allocate_1dim(C,n8)
END SUBROUTINE char_allocate_1dim_wrapper4
SUBROUTINE char_allocate_1dim(C,n)
   implicit none
   integer(kind=8),intent(in)         :: n
   CHARACTER(LEN=*),pointer :: C(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(C)
   ALLOCATE(C(n),STAT = IERR)
   nsize = mem_complexsize*size(C,KIND=long)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in char_allocate_1dim_n8',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in char_allocate_1dim_n8',nsize)
   ENDIF
   call mem_allocated_mem_character(nsize)
END SUBROUTINE char_allocate_1dim

!----- DEALLOCATE CHARACTER POINTERS -----!

SUBROUTINE char_deallocate_1dim(C)
   implicit none
   CHARACTER(LEN=*),pointer :: C(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = mem_complexsize*size(C,KIND=long)
   call mem_deallocated_mem_character(nsize)
   if (.not.ASSOCIATED(C)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in char_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(C,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in char_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in char_deallocate_1dim',nsize)
   ENDIF
END SUBROUTINE char_deallocate_1dim
!----- ALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic4_allocate_1dim_wrapper4(L,n)
   implicit none
   integer(kind=4),intent(in) :: n
   LOGICAL(kind=4),pointer            :: L(:)
   integer(kind=8)            :: n8
   n8=n
   call logic4_allocate_1dim(L,n8)
END SUBROUTINE logic4_allocate_1dim_wrapper4

SUBROUTINE logic8_allocate_1dim_wrapper4(L,n)
   implicit none
   integer(kind=4),intent(in) :: n
   LOGICAL(kind=8),pointer            :: L(:)
   integer(kind=8)            :: n8
   n8=n
   call logic8_allocate_1dim(L,n8)
END SUBROUTINE logic8_allocate_1dim_wrapper4

SUBROUTINE logic4_allocate_1dim(L,n)
   implicit none
   integer(kind=8),intent(in) :: n
   LOGICAL(kind=4),pointer    :: L(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(L)
   ALLOCATE(L(n),STAT = IERR)
   nsize = size(L,KIND=long)*mem_logicalsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in logic_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic4_allocate_1dim

SUBROUTINE logic8_allocate_1dim(L,n)
   implicit none
   integer(kind=8),intent(in) :: n
   LOGICAL(kind=8),pointer    :: L(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(L)
   ALLOCATE(L(n),STAT = IERR)
   nsize = size(L,KIND=long)*mem_logicalsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in logic_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic8_allocate_1dim

SUBROUTINE logic4_allocate_1dim_zero(L,n,z1)
   implicit none
   integer,intent(in) :: n
   logical(kind=4),intent(in) :: z1
   LOGICAL(kind=4),pointer    :: L(:)
   integer :: IERR,i1
   integer (kind=long) :: nsize
   nullify(L)
   i1=1
   if(z1)i1=0
   ALLOCATE(L(i1:n),STAT = IERR)
   nsize = size(L,KIND=long)*mem_logicalsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in logic_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic4_allocate_1dim_zero

SUBROUTINE logic8_allocate_1dim_zero(L,n,z1)
   implicit none
   integer,intent(in) :: n
   logical(kind=8),intent(in) :: z1
   LOGICAL(kind=8),pointer    :: L(:)
   integer :: IERR,i1
   integer (kind=long) :: nsize
   nullify(L)
   i1=1
   if(z1)i1=0
   ALLOCATE(L(i1:n),STAT = IERR)
   nsize = size(L,KIND=long)*mem_logicalsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in logic_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic8_allocate_1dim_zero

SUBROUTINE logic4_allocate_2dim(L,n1,n2)
   implicit none
   integer,intent(in) :: n1,n2
   LOGICAL(kind=4),pointer    :: L(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(L)
   ALLOCATE(L(n1,n2),STAT = IERR)
   nsize = size(L,KIND=long)*mem_logicalsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_allocate_2dim',IERR,n1,n2
      CALL MEMORY_ERROR_QUIT('Error in logic_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic4_allocate_2dim
SUBROUTINE logic8_allocate_2dim(L,n1,n2)
   implicit none
   integer,intent(in) :: n1,n2
   LOGICAL(kind=8),pointer    :: L(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(L)
   ALLOCATE(L(n1,n2),STAT = IERR)
   nsize = size(L,KIND=long)*mem_logicalsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_allocate_2dim',IERR,n1,n2
      CALL MEMORY_ERROR_QUIT('Error in logic_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic8_allocate_2dim

SUBROUTINE logic4_allocate_3dim(L,n1,n2,n3)
   implicit none
   integer,intent(in) :: n1,n2,n3
   LOGICAL(kind=4),pointer    :: L(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(L)
   ALLOCATE(L(n1,n2,n3),STAT = IERR)
   nsize = size(L,KIND=long)*mem_logicalsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_allocate_3dim',IERR,n1,n2,n3
      CALL MEMORY_ERROR_QUIT('Error in logic_allocate_3dim',nsize)
   ENDIF
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic4_allocate_3dim
SUBROUTINE logic8_allocate_3dim(L,n1,n2,n3)
   implicit none
   integer,intent(in) :: n1,n2,n3
   LOGICAL(kind=8),pointer    :: L(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(L)
   ALLOCATE(L(n1,n2,n3),STAT = IERR)
   nsize = size(L,KIND=long)*mem_logicalsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_allocate_3dim',IERR,n1,n2,n3
      CALL MEMORY_ERROR_QUIT('Error in logic_allocate_3dim',nsize)
   ENDIF
   call mem_allocated_mem_logical(nsize)
END SUBROUTINE logic8_allocate_3dim


!----- DEALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic4_deallocate_1dim(L)
   implicit none
   LOGICAL(kind=4),pointer :: L(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(L,KIND=long)*mem_logicalsize
   call mem_deallocated_mem_logical(nsize)
   if (.not.ASSOCIATED(L)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in logic_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(L,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in logic_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(L)
END SUBROUTINE logic4_deallocate_1dim

SUBROUTINE logic8_deallocate_1dim(L)
   implicit none
   LOGICAL(kind=8),pointer :: L(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(L,KIND=long)*mem_logicalsize
   call mem_deallocated_mem_logical(nsize)
   if (.not.ASSOCIATED(L)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in logic_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(L,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in logic_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(L)
END SUBROUTINE logic8_deallocate_1dim

SUBROUTINE logic4_deallocate_2dim(L)
   implicit none
   LOGICAL(kind=4),pointer :: L(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(L,KIND=long)*mem_logicalsize
   call mem_deallocated_mem_logical(nsize)
   if (.not.ASSOCIATED(L)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in logic_deallocate_2dim - memory previously released',nsize)
   endif
   DEALLOCATE(L,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_deallocate_2dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in logic_deallocate_2dim',nsize)
   ENDIF
   NULLIFY(L)
END SUBROUTINE logic4_deallocate_2dim

SUBROUTINE logic8_deallocate_2dim(L)
   implicit none
   LOGICAL(kind=8),pointer :: L(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(L,KIND=long)*mem_logicalsize
   call mem_deallocated_mem_logical(nsize)
   if (.not.ASSOCIATED(L)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in logic_deallocate_2dim - memory previously released',nsize)
   endif
   DEALLOCATE(L,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_deallocate_2dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in logic_deallocate_2dim',nsize)
   ENDIF
   NULLIFY(L)
END SUBROUTINE logic8_deallocate_2dim

SUBROUTINE logic4_deallocate_3dim(L)
   implicit none
   LOGICAL(kind=4),pointer :: L(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(L,KIND=long)*mem_logicalsize
   call mem_deallocated_mem_logical(nsize)
   if (.not.ASSOCIATED(L)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in logic_deallocate_3dim - memory previously released',nsize)
   endif
   DEALLOCATE(L,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_deallocate_3dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in logic_deallocate_3dim',nsize)
   ENDIF
   NULLIFY(L)
END SUBROUTINE logic4_deallocate_3dim

SUBROUTINE logic8_deallocate_3dim(L)
   implicit none
   LOGICAL(kind=8),pointer :: L(:,:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(L,KIND=long)*mem_logicalsize
   call mem_deallocated_mem_logical(nsize)
   if (.not.ASSOCIATED(L)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in logic_deallocate_3dim - memory previously released',nsize)
   endif
   DEALLOCATE(L,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in logic_deallocate_3dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in logic_deallocate_3dim',nsize)
   ENDIF
   NULLIFY(L)
END SUBROUTINE logic8_deallocate_3dim


!----- ALLOCATE AOBATCH POINTERS -----!

SUBROUTINE AOBATCH_allocate_1dim(AOBATCHITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(AOBATCH),pointer    :: AOBATCHITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(AOBATCHITEM)
   ALLOCATE(AOBATCHITEM(n),STAT = IERR)
   nsize = size(AOBATCHITEM,KIND=long)*mem_AOBATCHsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in AOBATCH_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in AOBATCH_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_AOBATCH(nsize)
END SUBROUTINE AOBATCH_allocate_1dim

SUBROUTINE AOBATCH_deallocate_1dim(AOBATCHITEM)
   implicit none
   TYPE(AOBATCH),pointer :: AOBATCHITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(AOBATCHITEM,KIND=long)*mem_AOBATCHsize
   call mem_deallocated_mem_AOBATCH(nsize)
   if (.not.ASSOCIATED(AOBATCHITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in AOBATCH_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(AOBATCHITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in AOBATCH_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in AOBATCH_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(AOBATCHITEM)
END SUBROUTINE AOBATCH_deallocate_1dim

! Overlap Type

SUBROUTINE OVERLAPT_allocate_1dim(OVERLAPITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(OVERLAP),pointer    :: OVERLAPITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(OVERLAPITEM)
   ALLOCATE(OVERLAPITEM(n),STAT = IERR)
   nsize = size(OVERLAPITEM,KIND=long)*mem_OVERLAPsize 
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in OVERLAP_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in OVERLAP_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_OVERLAPT(nsize)
END SUBROUTINE OVERLAPT_allocate_1dim

SUBROUTINE OVERLAPT_deallocate_1dim(OVERLAPITEM)
   implicit none
   TYPE(OVERLAP),pointer :: OVERLAPITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(OVERLAPITEM,KIND=long)*mem_OVERLAPsize
   call mem_deallocated_mem_OVERLAPT(nsize)
   if (.not.ASSOCIATED(OVERLAPITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in OVERLAP_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(OVERLAPITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in OVERLAP_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in OVERLAP_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(OVERLAPITEM)
END SUBROUTINE OVERLAPT_deallocate_1dim

!----- ALLOCATE DECORBITAL POINTERS -----!

SUBROUTINE DECORBITAL_allocate_1dim(DECORBITALITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(DECORBITAL),pointer    :: DECORBITALITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(DECORBITALITEM)
   ALLOCATE(DECORBITALITEM(n),STAT = IERR)
   nsize = size(DECORBITALITEM,KIND=long)*mem_DECORBITALsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in DECORBITAL_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in DECORBITAL_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_DECORBITAL(nsize)
END SUBROUTINE DECORBITAL_allocate_1dim

SUBROUTINE DECORBITAL_deallocate_1dim(DECORBITALITEM)
   implicit none
   TYPE(DECORBITAL),pointer :: DECORBITALITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(DECORBITALITEM,KIND=long)*mem_DECORBITALsize
   call mem_deallocated_mem_DECORBITAL(nsize)
   if (.not.ASSOCIATED(DECORBITALITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in DECORBITAL_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(DECORBITALITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in DECORBITAL_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in DECORBITAL_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(DECORBITALITEM)
END SUBROUTINE DECORBITAL_deallocate_1dim


!----- ALLOCATE DECFRAG POINTERS -----!

SUBROUTINE DECFRAG_allocate_1dim(DECFRAGITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(DECFRAG),pointer    :: DECFRAGITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(DECFRAGITEM)
   ALLOCATE(DECFRAGITEM(n),STAT = IERR)
   nsize = size(DECFRAGITEM,KIND=long)*mem_DECFRAGsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in DECFRAG_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in DECFRAG_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_DECFRAG(nsize)
END SUBROUTINE DECFRAG_allocate_1dim

SUBROUTINE DECFRAG_deallocate_1dim(DECFRAGITEM)
   implicit none
   TYPE(DECFRAG),pointer :: DECFRAGITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(DECFRAGITEM,KIND=long)*mem_DECFRAGsize
   call mem_deallocated_mem_DECFRAG(nsize)
   if (.not.ASSOCIATED(DECFRAGITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in DECFRAG_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(DECFRAGITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in DECFRAG_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in DECFRAG_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(DECFRAGITEM)
END SUBROUTINE DECFRAG_deallocate_1dim


!----- ALLOCATE BATCHTOORB POINTERS -----!

SUBROUTINE BATCHTOORB_allocate_1dim(BATCHTOORBITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(BATCHTOORB),pointer    :: BATCHTOORBITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(BATCHTOORBITEM)
   ALLOCATE(BATCHTOORBITEM(n),STAT = IERR)
   nsize = size(BATCHTOORBITEM,KIND=long)*mem_BATCHTOORBsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in BATCHTOORB_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in BATCHTOORB_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_BATCHTOORB(nsize)
END SUBROUTINE BATCHTOORB_allocate_1dim

SUBROUTINE BATCHTOORB_deallocate_1dim(BATCHTOORBITEM)
   implicit none
   TYPE(BATCHTOORB),pointer :: BATCHTOORBITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(BATCHTOORBITEM,KIND=long)*mem_BATCHTOORBsize
   call mem_deallocated_mem_BATCHTOORB(nsize)
   if (.not.ASSOCIATED(BATCHTOORBITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in BATCHTOORB_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(BATCHTOORBITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in BATCHTOORB_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in BATCHTOORB_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(BATCHTOORBITEM)
END SUBROUTINE BATCHTOORB_deallocate_1dim

!----- ALLOCATE DECAOBATCHINFO POINTERS -----!

SUBROUTINE DECAOBATCHINFO_allocate_1dim(DECAOBATCHINFOITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(DECAOBATCHINFO),pointer    :: DECAOBATCHINFOITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(DECAOBATCHINFOITEM)
   ALLOCATE(DECAOBATCHINFOITEM(n),STAT = IERR)
   nsize = size(DECAOBATCHINFOITEM,KIND=long)*mem_DECAOBATCHINFOsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in DECAOBATCHINFO_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in DECAOBATCHINFO_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_DECAOBATCHINFO(nsize)
END SUBROUTINE DECAOBATCHINFO_allocate_1dim

SUBROUTINE DECAOBATCHINFO_deallocate_1dim(DECAOBATCHINFOITEM)
   implicit none
   TYPE(DECAOBATCHINFO),pointer :: DECAOBATCHINFOITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(DECAOBATCHINFOITEM,KIND=long)*mem_DECAOBATCHINFOsize
   call mem_deallocated_mem_DECAOBATCHINFO(nsize)
   if (.not.ASSOCIATED(DECAOBATCHINFOITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in DECAOBATCHINFO_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(DECAOBATCHINFOITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in DECAOBATCHINFO_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in DECAOBATCHINFO_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(DECAOBATCHINFOITEM)
END SUBROUTINE DECAOBATCHINFO_deallocate_1dim

!----- ALLOCATE MYPOINTER POINTERS -----!

SUBROUTINE MYPOINTER_allocate_1dim(MYPOINTERITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(MYPOINTER),pointer    :: MYPOINTERITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(MYPOINTERITEM)
   ALLOCATE(MYPOINTERITEM(n),STAT = IERR)
   nsize = size(MYPOINTERITEM,KIND=long)*mem_MYPOINTERsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MYPOINTER_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in MYPOINTER_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_MYPOINTER(nsize)
END SUBROUTINE MYPOINTER_allocate_1dim

SUBROUTINE MYPOINTER_deallocate_1dim(MYPOINTERITEM)
   implicit none
   TYPE(MYPOINTER),pointer :: MYPOINTERITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(MYPOINTERITEM,KIND=long)*mem_MYPOINTERsize
   call mem_deallocated_mem_MYPOINTER(nsize)
   if (.not.ASSOCIATED(MYPOINTERITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in MYPOINTER_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(MYPOINTERITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MYPOINTER_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in MYPOINTER_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(MYPOINTERITEM)
END SUBROUTINE MYPOINTER_deallocate_1dim


SUBROUTINE MYPOINTER_allocate_2dim(MYPOINTERITEM,n1,n2)
   implicit none
   integer,intent(in) :: n1,n2
   TYPE(MYPOINTER),pointer    :: MYPOINTERITEM(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(MYPOINTERITEM)
   ALLOCATE(MYPOINTERITEM(n1,n2),STAT = IERR)
   nsize = size(MYPOINTERITEM,KIND=long)*mem_MYPOINTERsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MYPOINTER_allocate_2dim',IERR,n1,n2
      CALL MEMORY_ERROR_QUIT('Error in MYPOINTER_allocate_2dim',nsize)
   ENDIF
   call mem_allocated_mem_MYPOINTER(nsize)
END SUBROUTINE MYPOINTER_allocate_2dim

SUBROUTINE MYPOINTER_deallocate_2dim(MYPOINTERITEM)
   implicit none
   TYPE(MYPOINTER),pointer :: MYPOINTERITEM(:,:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(MYPOINTERITEM,KIND=long)*mem_MYPOINTERsize
   call mem_deallocated_mem_MYPOINTER(nsize)
   if (.not.ASSOCIATED(MYPOINTERITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in MYPOINTER_deallocate_2dim - memory previously released',nsize)
   endif
   DEALLOCATE(MYPOINTERITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MYPOINTER_deallocate_2dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in MYPOINTER_deallocate_2dim',nsize)
   ENDIF
   NULLIFY(MYPOINTERITEM)
END SUBROUTINE MYPOINTER_deallocate_2dim


!----- ALLOCATE fragmentAOS POINTERS -----!

SUBROUTINE fragmentAOS_allocate_1dim(fragmentAOSITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(fragmentAOS),pointer    :: fragmentAOSITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize

   nullify(fragmentAOSITEM)
   if(n==0) then
      ! No allocation
      return
   else
      ! Allocate
      ALLOCATE(fragmentAOSITEM(n),STAT = IERR)
      nsize = size(fragmentAOSITEM,KIND=long)*mem_fragmentAOSsize
      IF (IERR.NE. 0) THEN
         write(*,*) 'Error in fragmentAOS_allocate_1dim',IERR,n
         CALL MEMORY_ERROR_QUIT('Error in fragmentAOS_allocate_1dim',nsize)
      ENDIF
      call mem_allocated_mem_fragmentAOS(nsize)
   end if

END SUBROUTINE fragmentAOS_allocate_1dim

SUBROUTINE fragmentAOS_deallocate_1dim(fragmentAOSITEM)
   implicit none
   TYPE(fragmentAOS),pointer :: fragmentAOSITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   
   if(.not. associated(fragmentAOSITEM)) then
      ! This means that the dimension is zero, and fragmentAOSITEM was therefore never
      ! allocated (see fragmentAOS_allocate_1dim).
      ! We therefore just return.
      return
   end if

   nsize = size(fragmentAOSITEM,KIND=long)*mem_fragmentAOSsize
   call mem_deallocated_mem_fragmentAOS(nsize)
   if (.not.ASSOCIATED(fragmentAOSITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in fragmentAOS_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(fragmentAOSITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in fragmentAOS_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in fragmentAOS_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(fragmentAOSITEM)

END SUBROUTINE fragmentAOS_deallocate_1dim



!----- ALLOCATE ARRAY2 POINTERS -----!

SUBROUTINE ARRAY2_allocate_1dim(ARRAY2ITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(ARRAY2),pointer    :: ARRAY2ITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(ARRAY2ITEM)
   ALLOCATE(ARRAY2ITEM(n),STAT = IERR)
   nsize = size(ARRAY2ITEM,KIND=long)*mem_ARRAY2size
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ARRAY2_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in ARRAY2_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_ARRAY2(nsize)
END SUBROUTINE ARRAY2_allocate_1dim

SUBROUTINE ARRAY2_deallocate_1dim(ARRAY2ITEM)
   implicit none
   TYPE(ARRAY2),pointer :: ARRAY2ITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(ARRAY2ITEM,KIND=long)*mem_ARRAY2size
   call mem_deallocated_mem_ARRAY2(nsize)
   if (.not.ASSOCIATED(ARRAY2ITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in ARRAY2_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(ARRAY2ITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ARRAY2_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in ARRAY2_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(ARRAY2ITEM)
END SUBROUTINE ARRAY2_deallocate_1dim

!----- ALLOCATE ARRAY4 POINTERS -----!

SUBROUTINE ARRAY4_allocate_1dim(ARRAY4ITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(ARRAY4),pointer    :: ARRAY4ITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(ARRAY4ITEM)
   ALLOCATE(ARRAY4ITEM(n),STAT = IERR)
   nsize = size(ARRAY4ITEM,KIND=long)*mem_ARRAY4size
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ARRAY4_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in ARRAY4_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_ARRAY4(nsize)
END SUBROUTINE ARRAY4_allocate_1dim

SUBROUTINE ARRAY4_deallocate_1dim(ARRAY4ITEM)
   implicit none
   TYPE(ARRAY4),pointer :: ARRAY4ITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(ARRAY4ITEM,KIND=long)*mem_ARRAY4size
   call mem_deallocated_mem_ARRAY4(nsize)
   if (.not.ASSOCIATED(ARRAY4ITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in ARRAY4_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(ARRAY4ITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ARRAY4_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in ARRAY4_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(ARRAY4ITEM)
END SUBROUTINE ARRAY4_deallocate_1dim

!----- ALLOCATE ARRAY POINTERS -----!

SUBROUTINE tensor_allocate_1dim(ARRAYITEM,n)
   implicit none
   integer,intent(in) :: n
   type(tensor),pointer    :: ARRAYITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(ARRAYITEM)
   ALLOCATE(ARRAYITEM(n),STAT = IERR)
   nsize = size(ARRAYITEM,KIND=long)*mem_ARRAYsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in tensor_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in tensor_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_ARRAY(nsize)
END SUBROUTINE tensor_allocate_1dim

SUBROUTINE tensor_deallocate_1dim(ARRAYITEM)
   implicit none
   type(tensor),pointer :: ARRAYITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(ARRAYITEM,KIND=long)*mem_ARRAYsize
   call mem_deallocated_mem_ARRAY(nsize)
   if (.not.ASSOCIATED(ARRAYITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in tensor_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(ARRAYITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in tensor_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in tensor_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(ARRAYITEM)
END SUBROUTINE tensor_deallocate_1dim

!----- ALLOCATE MP2DENS POINTERS -----!

SUBROUTINE MP2DENS_allocate_1dim(MP2DENSITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(MP2DENS),pointer    :: MP2DENSITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(MP2DENSITEM)
   ALLOCATE(MP2DENSITEM(n),STAT = IERR)
   nsize = size(MP2DENSITEM,KIND=long)*mem_MP2DENSsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MP2DENS_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in MP2DENS_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_MP2DENS(nsize)
END SUBROUTINE MP2DENS_allocate_1dim

SUBROUTINE MP2DENS_deallocate_1dim(MP2DENSITEM)
   implicit none
   TYPE(MP2DENS),pointer :: MP2DENSITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(MP2DENSITEM,KIND=long)*mem_MP2DENSsize
   call mem_deallocated_mem_MP2DENS(nsize)
   if (.not.ASSOCIATED(MP2DENSITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in MP2DENS_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(MP2DENSITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MP2DENS_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in MP2DENS_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(MP2DENSITEM)
END SUBROUTINE MP2DENS_deallocate_1dim

!----- ALLOCATE TRACEBACK POINTERS -----!

SUBROUTINE TRACEBACK_allocate_1dim(TRACEBACKITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(TRACEBACK),pointer    :: TRACEBACKITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(TRACEBACKITEM)
   ALLOCATE(TRACEBACKITEM(n),STAT = IERR)
   nsize = size(TRACEBACKITEM,KIND=long)*mem_TRACEBACKsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in TRACEBACK_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in TRACEBACK_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_TRACEBACK(nsize)
END SUBROUTINE TRACEBACK_allocate_1dim

SUBROUTINE TRACEBACK_deallocate_1dim(TRACEBACKITEM)
   implicit none
   TYPE(TRACEBACK),pointer :: TRACEBACKITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(TRACEBACKITEM,KIND=long)*mem_TRACEBACKsize
   call mem_deallocated_mem_TRACEBACK(nsize)
   if (.not.ASSOCIATED(TRACEBACKITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in TRACEBACK_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(TRACEBACKITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in TRACEBACK_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in TRACEBACK_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(TRACEBACKITEM)
END SUBROUTINE TRACEBACK_deallocate_1dim

!----- ALLOCATE MP2GRAD POINTERS -----!

SUBROUTINE MP2GRAD_allocate_1dim(MP2GRADITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(MP2GRAD),pointer    :: MP2GRADITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(MP2GRADITEM)
   ALLOCATE(MP2GRADITEM(n),STAT = IERR)
   nsize = size(MP2GRADITEM,KIND=long)*mem_MP2GRADsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MP2GRAD_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in MP2GRAD_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_MP2GRAD(nsize)
END SUBROUTINE MP2GRAD_allocate_1dim

SUBROUTINE MP2GRAD_deallocate_1dim(MP2GRADITEM)
   implicit none
   TYPE(MP2GRAD),pointer :: MP2GRADITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(MP2GRADITEM,KIND=long)*mem_MP2GRADsize
   call mem_deallocated_mem_MP2GRAD(nsize)
   if (.not.ASSOCIATED(MP2GRADITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in MP2GRAD_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(MP2GRADITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MP2GRAD_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in MP2GRAD_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(MP2GRADITEM)
END SUBROUTINE MP2GRAD_deallocate_1dim



!----- ALLOCATE PNOSpaceInfo POINTERS -----!
SUBROUTINE PNOSPACEINFO_allocate_1dim(SpaceItem,n)
  implicit none
  integer,intent(in) :: n
  type(PNOSpaceInfo), intent(inout), pointer :: SpaceItem(:)
  integer :: ierr
  integer(kind=long) :: nsize

  SpaceItem => null()

  allocate(SpaceItem(n),STAT = ierr )

  if (ierr.ne. 0) then
     WRITE(*,*) 'ERROR(PNOSPACEINFO_allocate_1dim) Status/nel',ierr,"/",n
     call memory_error_quit('ERROR(PNOSPACEINFO_allocate_1dim) status /= 0',nsize)
  endif

  nsize = size(SpaceItem,kind=long)*mem_PNOSpaceInfosize
  call mem_allocated_mem_PNOSpaceInfo(nsize)
 
END SUBROUTINE PNOSPACEINFO_allocate_1dim


SUBROUTINE PNOSpaceInfo_deallocate_1dim(PNOSpaceInfoITEM)
  implicit none
  TYPE(PNOSpaceInfo),pointer,intent(inout) :: PNOSpaceInfoITEM(:)
  integer :: IERR
  integer (kind=long) :: nsize
  nsize = size(PNOSpaceInfoITEM,KIND=long)*mem_PNOSpaceInfosize
  call mem_deallocated_mem_PNOSpaceInfo(nsize)

  if (.not.ASSOCIATED(PNOSpaceInfoITEM)) then
     print *,'Memory previously released!!'
     call memory_error_quit('Error in PNOSpaceInfo_deallocate_1dim - memory previously released',nsize)
  endif
  DEALLOCATE(PNOSpaceInfoITEM,STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in PNOSpaceInfo_deallocate_1dim',IERR
     CALL MEMORY_ERROR_QUIT('Error in PNOSpaceInfo_deallocate_1dim',nsize)
  ENDIF
  NULLIFY(PNOSpaceInfoITEM)
END SUBROUTINE PNOSpaceInfo_deallocate_1dim



!----- ALLOCATE ODBATCH POINTERS -----!

SUBROUTINE ODBATCH_allocate_1dim(ODBATCHITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(ODBATCH),pointer    :: ODBATCHITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(ODBATCHITEM)
   ALLOCATE(ODBATCHITEM(n),STAT = IERR)
   nsize = size(ODBATCHITEM,KIND=long)*mem_ODBATCHsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ODBATCH_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in ODBATCH_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_ODBATCH(nsize)
END SUBROUTINE ODBATCH_allocate_1dim

SUBROUTINE ODBATCH_deallocate_1dim(ODBATCHITEM)
   implicit none
   TYPE(ODBATCH),pointer :: ODBATCHITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(ODBATCHITEM,KIND=long)*mem_ODBATCHsize
   call mem_deallocated_mem_ODBATCH(nsize)
   if (.not.ASSOCIATED(ODBATCHITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in ODBATCH_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(ODBATCHITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ODBATCH_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in ODBATCH_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(ODBATCHITEM)
END SUBROUTINE ODBATCH_deallocate_1dim

!----- ALLOCATE LSAOTENSOR POINTERS -----!

SUBROUTINE LSAOTENSOR_allocate_1dim(LSAOTENSORITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(LSAOTENSOR),pointer    :: LSAOTENSORITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(LSAOTENSORITEM)
   ALLOCATE(LSAOTENSORITEM(n),STAT = IERR)
   nsize = size(LSAOTENSORITEM,KIND=long)*mem_LSAOTENSORsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in LSAOTENSOR_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in LSAOTENSOR_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_LSAOTENSOR(nsize)
END SUBROUTINE LSAOTENSOR_allocate_1dim

SUBROUTINE LSAOTENSOR_deallocate_1dim(LSAOTENSORITEM)
   implicit none
   TYPE(LSAOTENSOR),pointer :: LSAOTENSORITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(LSAOTENSORITEM,KIND=long)*mem_LSAOTENSORsize
   call mem_deallocated_mem_LSAOTENSOR(nsize)
   if (.not.ASSOCIATED(LSAOTENSORITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in LSAOTENSOR_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(LSAOTENSORITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in LSAOTENSOR_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in LSAOTENSOR_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(LSAOTENSORITEM)
END SUBROUTINE LSAOTENSOR_deallocate_1dim

!----- ALLOCATE SLSAOTENSOR POINTERS -----!

SUBROUTINE SLSAOTENSOR_allocate_1dim(SLSAOTENSORITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(SLSAOTENSOR),pointer    :: SLSAOTENSORITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(SLSAOTENSORITEM)
   ALLOCATE(SLSAOTENSORITEM(n),STAT = IERR)
   nsize = size(SLSAOTENSORITEM,KIND=long)*mem_SLSAOTENSORsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in SLSAOTENSOR_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in SLSAOTENSOR_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_SLSAOTENSOR(nsize)
END SUBROUTINE SLSAOTENSOR_allocate_1dim

SUBROUTINE SLSAOTENSOR_deallocate_1dim(SLSAOTENSORITEM)
   implicit none
   TYPE(SLSAOTENSOR),pointer :: SLSAOTENSORITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(SLSAOTENSORITEM,KIND=long)*mem_SLSAOTENSORsize
   call mem_deallocated_mem_SLSAOTENSOR(nsize)
   if (.not.ASSOCIATED(SLSAOTENSORITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in SLSAOTENSOR_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(SLSAOTENSORITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in SLSAOTENSOR_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in SLSAOTENSOR_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(SLSAOTENSORITEM)
END SUBROUTINE SLSAOTENSOR_deallocate_1dim

!----- ALLOCATE GLOBALLSAOTENSOR POINTERS -----!

SUBROUTINE GLOBALLSAOTENSOR_allocate_1dim(GLOBALLSAOTENSORITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(GLOBALLSAOTENSOR),pointer    :: GLOBALLSAOTENSORITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(GLOBALLSAOTENSORITEM)
   ALLOCATE(GLOBALLSAOTENSORITEM(n),STAT = IERR)
   nsize = size(GLOBALLSAOTENSORITEM,KIND=long)*mem_GLOBALLSAOTENSORsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in GLOBALLSAOTENSOR_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in GLOBALLSAOTENSOR_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_GLOBALLSAOTENSOR(nsize)
END SUBROUTINE GLOBALLSAOTENSOR_allocate_1dim

SUBROUTINE GLOBALLSAOTENSOR_deallocate_1dim(GLOBALLSAOTENSORITEM)
   implicit none
   TYPE(GLOBALLSAOTENSOR),pointer :: GLOBALLSAOTENSORITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(GLOBALLSAOTENSORITEM,KIND=long)*mem_GLOBALLSAOTENSORsize
   call mem_deallocated_mem_GLOBALLSAOTENSOR(nsize)
   if (.not.ASSOCIATED(GLOBALLSAOTENSORITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in GLOBALLSAOTENSOR_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(GLOBALLSAOTENSORITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in GLOBALLSAOTENSOR_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in GLOBALLSAOTENSOR_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(GLOBALLSAOTENSORITEM)
END SUBROUTINE GLOBALLSAOTENSOR_deallocate_1dim

!----- ALLOCATE ATOMTYPEITEM POINTERS -----!

SUBROUTINE ATOMTYPEITEM_allocate_1dim(ATOMTYPEITEMITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(ATOMTYPEITEM),pointer    :: ATOMTYPEITEMITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(ATOMTYPEITEMITEM)
   ALLOCATE(ATOMTYPEITEMITEM(n),STAT = IERR)
   nsize = size(ATOMTYPEITEMITEM,KIND=long)*mem_ATOMTYPEITEMsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ATOMTYPEITEM_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in ATOMTYPEITEM_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_ATOMTYPEITEM(nsize)
END SUBROUTINE ATOMTYPEITEM_allocate_1dim

SUBROUTINE ATOMTYPEITEM_deallocate_1dim(ATOMTYPEITEMITEM)
   implicit none
   TYPE(ATOMTYPEITEM),pointer :: ATOMTYPEITEMITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(ATOMTYPEITEMITEM,KIND=long)*mem_ATOMTYPEITEMsize
   call mem_deallocated_mem_ATOMTYPEITEM(nsize)
   if (.not.ASSOCIATED(ATOMTYPEITEMITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in ATOMTYPEITEM_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(ATOMTYPEITEMITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ATOMTYPEITEM_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in ATOMTYPEITEM_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(ATOMTYPEITEMITEM)
END SUBROUTINE ATOMTYPEITEM_deallocate_1dim

!----- ALLOCATE ATOMTYPEITEM POINTERS -----!

SUBROUTINE ATOMITEM_allocate_1dim(ATOMITEMITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(ATOMITEM),pointer    :: ATOMITEMITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(ATOMITEMITEM)
   ALLOCATE(ATOMITEMITEM(n),STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ATOMITEM_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in ATOMITEM_allocate_1dim',nsize)
   ENDIF
   nsize = size(ATOMITEMITEM,KIND=long)*mem_ATOMITEMsize
   call mem_allocated_mem_ATOMITEM(nsize)
END SUBROUTINE ATOMITEM_allocate_1dim

SUBROUTINE ATOMITEM_deallocate_1dim(ATOMITEMITEM)
   implicit none
   TYPE(ATOMITEM),pointer :: ATOMITEMITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(ATOMITEMITEM,KIND=long)*mem_ATOMITEMsize
   call mem_deallocated_mem_ATOMITEM(nsize)
   if (.not.ASSOCIATED(ATOMITEMITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in ATOMITEM_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(ATOMITEMITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in ATOMITEM_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in ATOMITEM_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(ATOMITEMITEM)
END SUBROUTINE ATOMITEM_deallocate_1dim

!----- ALLOCATE LSMATRIX POINTERS -----!

SUBROUTINE LSMATRIX_allocate_1dim(LSMATRIXITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(LSMATRIX),pointer    :: LSMATRIXITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(LSMATRIXITEM)
   ALLOCATE(LSMATRIXITEM(n),STAT = IERR)
   nsize = size(LSMATRIXITEM,KIND=long)*mem_LSMATRIXsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in LSMATRIX_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in LSMATRIX_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_LSMATRIX(nsize)
END SUBROUTINE LSMATRIX_allocate_1dim

SUBROUTINE LSMATRIX_deallocate_1dim(LSMATRIXITEM)
   implicit none
   TYPE(LSMATRIX),pointer :: LSMATRIXITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(LSMATRIXITEM,KIND=long)*mem_LSMATRIXsize
   call mem_deallocated_mem_LSMATRIX(nsize)
   if (.not.ASSOCIATED(LSMATRIXITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in LSMATRIX_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(LSMATRIXITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in LSMATRIX_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in LSMATRIX_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(LSMATRIXITEM)
END SUBROUTINE LSMATRIX_deallocate_1dim

!SUBROUTINE LSMATRIXP_allocate_1dim(LSMATRIXITEM,n)
!implicit none
!integer,intent(in) :: n
!TYPE(LSMATRIXpointer),pointer    :: LSMATRIXITEM(:)
!integer :: IERR
!integer (kind=long) :: nsize
!nullify(LSMATRIXITEM)
!ALLOCATE(LSMATRIXITEM(n),STAT = IERR)
!IF (IERR.NE. 0) THEN
!   write(*,*) 'Error in LSMATRIX_allocate_1dim',IERR,n
!   CALL MEMORY_ERROR_QUIT('Error in LSMATRIX_allocate_1dim',nsize)
!ENDIF
!nsize = size(LSMATRIXITEM,KIND=long)*8
!call mem_allocated_mem_LSMATRIX(nsize)
!END SUBROUTINE LSMATRIXP_allocate_1dim
!
!SUBROUTINE LSMATRIXP_deallocate_1dim(LSMATRIXITEM)
!implicit none
!TYPE(LSMATRIXpointer),pointer :: LSMATRIXITEM(:)
!integer :: IERR
!integer (kind=long) :: nsize
!   nsize = size(LSMATRIXITEM,KIND=long)*8
!   call mem_deallocated_mem_LSMATRIX(nsize)
!   if (.not.ASSOCIATED(LSMATRIXITEM)) then
!      print *,'Memory previously released!!'
!      call memory_error_quit('Error in LSMATRIX_deallocate_1dim - memory previously released',nsize)
!   endif
!   DEALLOCATE(LSMATRIXITEM,STAT = IERR)
!   IF (IERR.NE. 0) THEN
!      write(*,*) 'Error in LSMATRIX_deallocate_1dim',IERR
!      CALL MEMORY_ERROR_QUIT('Error in LSMATRIX_deallocate_1dim',nsize)
!   ENDIF
!   NULLIFY(LSMATRIXITEM)
!END SUBROUTINE LSMATRIXP_deallocate_1dim

!----- ALLOCATE MATRIX POINTERS -----!

SUBROUTINE MATRIX_allocate_1dim(MATRIXITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(MATRIX),pointer    :: MATRIXITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(MATRIXITEM)
   ALLOCATE(MATRIXITEM(n),STAT = IERR)
   nsize = size(MATRIXITEM,KIND=long)*mem_MATRIXsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MATRIX_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in MATRIX_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_type_matrix(nsize)
END SUBROUTINE MATRIX_allocate_1dim

SUBROUTINE MATRIX_deallocate_1dim(MATRIXITEM)
   implicit none
   TYPE(MATRIX),pointer :: MATRIXITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(MATRIXITEM,KIND=long)*mem_MATRIXsize
   call mem_deallocated_mem_type_matrix(nsize)
   if (.not.ASSOCIATED(MATRIXITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in MATRIX_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(MATRIXITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MATRIX_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in MATRIX_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(MATRIXITEM)
END SUBROUTINE MATRIX_deallocate_1dim

SUBROUTINE MATRIXP_allocate_1dim(MATRIXITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(MATRIXp),pointer    :: MATRIXITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(MATRIXITEM)
   ALLOCATE(MATRIXITEM(n),STAT = IERR)
   nsize = size(MATRIXITEM,KIND=long)*8
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MATRIX_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in MATRIX_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_type_matrix(nsize)
END SUBROUTINE MATRIXP_allocate_1dim

SUBROUTINE MATRIXP_deallocate_1dim(MATRIXITEM)
   implicit none
   TYPE(MATRIXp),pointer :: MATRIXITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(MATRIXITEM,KIND=long)*8
   call mem_deallocated_mem_type_matrix(nsize)
   if (.not.ASSOCIATED(MATRIXITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in MATRIX_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(MATRIXITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in MATRIX_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in MATRIX_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(MATRIXITEM)
END SUBROUTINE MATRIXP_deallocate_1dim

#ifdef MOD_UNRELEASED
SUBROUTINE Lvec_data_allocate_1dim(Lvec_dataITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(Lvec_data_t),pointer    :: Lvec_dataITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(Lvec_dataITEM)
   ALLOCATE(Lvec_dataITEM(n),STAT = IERR)
   nsize = size(Lvec_dataITEM,KIND=long)*mem_Lvec_datasize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in Lvec_data_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in Lvec_data_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_Lvec_data(nsize)
END SUBROUTINE Lvec_data_allocate_1dim

SUBROUTINE Lvec_data_deallocate_1dim(Lvec_dataITEM)
   implicit none
   TYPE(Lvec_data_t),pointer :: Lvec_dataITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(Lvec_dataITEM,KIND=long)*mem_Lvec_datasize
   call mem_deallocated_mem_Lvec_data(nsize)
   if (.not.ASSOCIATED(Lvec_dataITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in Lvec_data_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(Lvec_dataITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in Lvec_data_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in Lvec_data_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(Lvec_dataITEM)
END SUBROUTINE Lvec_data_deallocate_1dim

SUBROUTINE Lattice_cell_allocate_1dim(Lattice_cellITEM,n)
   implicit none
   integer,intent(in) :: n
   TYPE(Lattice_cell_info_t),pointer    :: Lattice_cellITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nullify(Lattice_cellITEM)
   ALLOCATE(Lattice_cellITEM(n),STAT = IERR)
   nsize = size(Lattice_cellITEM,KIND=long)*mem_Lattice_cellsize
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in Lattice_cell_allocate_1dim',IERR,n
      CALL MEMORY_ERROR_QUIT('Error in Lattice_cell_allocate_1dim',nsize)
   ENDIF
   call mem_allocated_mem_Lattice_cell(nsize)
END SUBROUTINE Lattice_cell_allocate_1dim

SUBROUTINE Lattice_cell_deallocate_1dim(Lattice_cellITEM)
   implicit none
   TYPE(Lattice_cell_info_t),pointer :: Lattice_cellITEM(:)
   integer :: IERR
   integer (kind=long) :: nsize
   nsize = size(Lattice_cellITEM,KIND=long)*mem_Lattice_cellsize
   call mem_deallocated_mem_Lattice_cell(nsize)
   if (.not.ASSOCIATED(Lattice_cellITEM)) then
      print *,'Memory previously released!!'
      call memory_error_quit('Error in Lattice_cell_deallocate_1dim - memory previously released',nsize)
   endif
   DEALLOCATE(Lattice_cellITEM,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in Lattice_cell_deallocate_1dim',IERR
      CALL MEMORY_ERROR_QUIT('Error in Lattice_cell_deallocate_1dim',nsize)
   ENDIF
   NULLIFY(Lattice_cellITEM)
END SUBROUTINE Lattice_cell_deallocate_1dim


#endif



!----- MEMORY HANDLING -----!

!subroutine that adds the memory used in external libraries. 
subroutine mem_add_external_memory(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global+nsize)
   max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global+nsize)
   if (mem_allocated_global  + nsize < 0) then
      write(*,*) 'Allocated Memory                               =', mem_allocated_global
      write(*,*) 'Added Memory                                   =', nsize
      write(*,*) 'Total memory negative! mem_add_external_memory =', mem_allocated_global + nsize
      IF(mem_allocated_global.GT.0.AND.nsize.GT.0)THEN
         call memory_error_quit('Error in mem_add_external_memory - probably integer overflow!',nsize)
      else
         call memory_error_quit('Error in mem_add_external_memory!',nsize)
      endif
   endif
end subroutine mem_add_external_memory

subroutine mem_allocated_mem_real(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_real = mem_tp_allocated_real + nsize
      if (mem_tp_allocated_real < 0) then
         write(*,*) 'Real memory negative! mem_tp_allocated_real =', mem_tp_allocated_real
         write(*,*) 'Real memory negative! nsize =', nsize
         call memory_error_quit('Error in mem_allocated_mem_real - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_real = MAX(max_mem_tp_used_real,mem_tp_allocated_real)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_allocated_mem_tp_real - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   else
      mem_allocated_real = mem_allocated_real + nsize
      if (mem_allocated_real < 0) then
         write(*,*) 'Real memory negative! mem_allocated_real =', mem_allocated_real
         write(*,*) 'Real memory negative! nsize =', nsize
         call memory_error_quit('Error in mem_allocated_mem_real - probably integer overflow!',nsize)
      endif
      max_mem_used_real = MAX(max_mem_used_real,mem_allocated_real)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_allocated_mem_real - probably integer overflow!',nsize)
      endif
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_real'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_real'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_real')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   endif
end subroutine mem_allocated_mem_real

subroutine mem_deallocated_mem_real(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize

   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_real = mem_tp_allocated_real - nsize
      if (mem_tp_allocated_real < 0) then
         write(*,*) 'Real memory negative! mem_tp_allocated_real =', mem_tp_allocated_real
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_real - something wrong with deallocation!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_real - something wrong with deallocation!',nsize)
      endif
   ELSE
      mem_allocated_real = mem_allocated_real - nsize
      if (mem_allocated_real < 0) then
         write(*,*) 'Real memory negative! mem_allocated_real =', mem_allocated_real
         call memory_error_quit('Error in mem_deallocated_mem_real - something wrong with deallocation!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_deallocated_mem_real - something wrong with deallocation!',nsize)
      endif
   ENDIF

end subroutine mem_deallocated_mem_real

subroutine mem_allocated_mem_mpi(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_mpi = mem_tp_allocated_mpi + nsize
      if (mem_tp_allocated_mpi < 0) then
         write(*,*) 'Real memory negative! mem_tp_allocated_mpi =', mem_tp_allocated_mpi
         write(*,*) 'Real memory negative! nsize =', nsize
         call memory_error_quit('Error in mem_allocated_mem_mpi - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_mpi = MAX(max_mem_tp_used_mpi,mem_tp_allocated_mpi)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_allocated_mem_tp_mpi - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   else
      mem_allocated_mpi = mem_allocated_mpi + nsize
      if (mem_allocated_mpi < 0) then
         write(*,*) 'Real memory negative! mem_allocated_mpi =', mem_allocated_mpi
         write(*,*) 'Real memory negative! nsize =', nsize
         call memory_error_quit('Error in mem_allocated_mem_mpi - probably integer overflow!',nsize)
      endif
      max_mem_used_mpi = MAX(max_mem_used_mpi,mem_allocated_mpi)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_allocated_mem_mpi - probably integer overflow!',nsize)
      endif
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_mpi'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_mpi'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_mpi')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   endif
end subroutine mem_allocated_mem_mpi

subroutine mem_deallocated_mem_mpi(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize

   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_mpi = mem_tp_allocated_mpi - nsize
      if (mem_tp_allocated_mpi < 0) then
         write(*,*) 'Real memory negative! mem_tp_allocated_mpi =', mem_tp_allocated_mpi
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_mpi - something wrong with deallocation!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_mpi - something wrong with deallocation!',nsize)
      endif
   ELSE
      mem_allocated_mpi = mem_allocated_mpi - nsize
      if (mem_allocated_mpi < 0) then
         write(*,*) 'Real memory negative! mem_allocated_mpi =', mem_allocated_mpi
         call memory_error_quit('Error in mem_deallocated_mem_mpi - something wrong with deallocation!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_deallocated_mem_mpi - something wrong with deallocation!',nsize)
      endif
   ENDIF

end subroutine mem_deallocated_mem_mpi


subroutine mem_allocated_mem_complex(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_complex = mem_tp_allocated_complex + nsize
      if (mem_tp_allocated_complex < 0) then
         write(*,*) 'Complex memory negative! mem_tp_allocated_complex =', mem_tp_allocated_complex
         write(*,*) 'Complex memory negative! nsize =', nsize
         call memory_error_quit('Error in mem_allocated_mem_complex - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_complex = MAX(max_mem_tp_used_complex,mem_tp_allocated_complex)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_allocated_mem_tp_complex - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   else
      mem_allocated_complex = mem_allocated_complex + nsize
      if (mem_allocated_complex < 0) then
         write(*,*) 'Complex memory negative! mem_allocated_Complex =', mem_allocated_complex
         write(*,*) 'Complex memory negative! nsize =', nsize
         call memory_error_quit('Error in mem_allocated_mem_Complex - probably integer overflow!',nsize)
      endif
      max_mem_used_Complex = MAX(max_mem_used_Complex,mem_allocated_Complex)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_allocated_mem_Complex - probably integer overflow!',nsize)
      endif
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_complex'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_complex'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_complex')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   endif
end subroutine mem_allocated_mem_complex

subroutine mem_deallocated_mem_complex(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize

   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_complex = mem_tp_allocated_complex - nsize
      if (mem_tp_allocated_complex < 0) then
         write(*,*) 'Complex memory negative! mem_tp_allocated_complex =', mem_tp_allocated_complex
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_complex - something wrong with deallocation!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_complex - something wrong with deallocation!',nsize)
      endif
   ELSE
      mem_allocated_complex = mem_allocated_complex - nsize
      if (mem_allocated_complex < 0) then
         write(*,*) 'Real memory negative! mem_allocated_complex =', mem_allocated_complex
         call memory_error_quit('Error in mem_deallocated_mem_complex - something wrong with deallocation!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_deallocated_mem_complex - something wrong with deallocation!',nsize)
      endif
   ENDIF

end subroutine mem_deallocated_mem_complex


subroutine mem_allocated_mem_integer(nsize)
   implicit none
   integer(kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_integer = mem_tp_allocated_integer + nsize
      max_mem_tp_used_integer = MAX(max_mem_tp_used_integer,mem_tp_allocated_integer)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_integer = mem_allocated_integer + nsize
      max_mem_used_integer = MAX(max_mem_used_integer,mem_allocated_integer)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_integer'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_integer'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_integer')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_integer

subroutine mem_deallocated_mem_integer(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
#ifdef VAR_OMP
   integer, external :: OMP_GET_THREAD_NUM
#endif
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_integer = mem_tp_allocated_integer - nsize
      if (mem_tp_allocated_integer < 0) then
#ifdef VAR_OMP
         print*,'Local: OMP THREAD ID:',OMP_GET_THREAD_NUM()
#endif
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_integer - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
#ifdef VAR_OMP
         print*,'Global: OMP THREAD ID:',OMP_GET_THREAD_NUM()
#endif
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_integer - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_integer = mem_allocated_integer - nsize
      if (mem_allocated_integer < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_integer1 - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_integer2 - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_integer

subroutine mem_allocated_mem_character(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_character = mem_tp_allocated_character + nsize
      max_mem_tp_used_character = MAX(max_mem_tp_used_character,mem_tp_allocated_character)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_character = mem_allocated_character + nsize
      max_mem_used_character = MAX(max_mem_used_character,mem_allocated_character)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_character'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_character'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_character')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_character

subroutine mem_deallocated_mem_character(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN !we add to thread private variables
      mem_tp_allocated_character = mem_tp_allocated_character - nsize
      if (mem_tp_allocated_character < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_character - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_character - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_character = mem_allocated_character - nsize
      if (mem_allocated_character < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_character - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_character - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_character

subroutine mem_allocated_mem_logical(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_logical = mem_tp_allocated_logical + nsize
      max_mem_tp_used_logical = MAX(max_mem_tp_used_logical,mem_tp_allocated_logical)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_logical = mem_allocated_logical + nsize
      max_mem_used_logical = MAX(max_mem_used_logical,mem_allocated_logical)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_logical'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_logical'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_logical')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_logical

subroutine mem_deallocated_mem_logical(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_logical = mem_tp_allocated_logical - nsize
      if (mem_tp_allocated_logical < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_logical - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_logical - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_logical = mem_allocated_logical - nsize
      if (mem_allocated_logical < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_logical - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_logical - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_logical

subroutine mem_allocated_mem_AOBATCH(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_AOBATCH = mem_tp_allocated_AOBATCH + nsize
      max_mem_tp_used_AOBATCH = MAX(max_mem_tp_used_AOBATCH,mem_tp_allocated_AOBATCH)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_AOBATCH = mem_allocated_AOBATCH + nsize
      max_mem_used_AOBATCH = MAX(max_mem_used_AOBATCH,mem_allocated_AOBATCH)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_AOBATCH'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_AOBATCH'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_AOBATCH')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_AOBATCH

subroutine mem_deallocated_mem_AOBATCH(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_AOBATCH = mem_tp_allocated_AOBATCH - nsize
      if (mem_tp_allocated_AOBATCH < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_AOBATCH - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_AOBATCH - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_AOBATCH = mem_allocated_AOBATCH - nsize
      if (mem_allocated_AOBATCH < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_AOBATCH - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_AOBATCH - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_AOBATCH

subroutine mem_allocated_mem_OVERLAPT(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_OVERLAPT = mem_tp_allocated_OVERLAPT + nsize
      max_mem_tp_used_OVERLAPT = MAX(max_mem_tp_used_OVERLAPT,mem_tp_allocated_OVERLAPT)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_OVERLAPT = mem_allocated_OVERLAPT + nsize
      max_mem_used_OVERLAPT = MAX(max_mem_used_OVERLAPT,mem_allocated_OVERLAPT)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_OVERLAPT'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_OVERLAPT'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_OVERLAPT')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_OVERLAPT

subroutine mem_deallocated_mem_OVERLAPT(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_OVERLAPT = mem_tp_allocated_OVERLAPT - nsize
      if (mem_tp_allocated_OVERLAPT < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_OVERLAPT - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_OVERLAPT - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_OVERLAPT = mem_allocated_OVERLAPT - nsize
      if (mem_allocated_OVERLAPT < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_OVERLAPT - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_OVERLAPT - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_OVERLAPT

subroutine mem_allocated_mem_DECORBITAL(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_DECORBITAL = mem_tp_allocated_DECORBITAL + nsize
      max_mem_tp_used_DECORBITAL = MAX(max_mem_tp_used_DECORBITAL,mem_tp_allocated_DECORBITAL)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_DECORBITAL = mem_allocated_DECORBITAL + nsize
      max_mem_used_DECORBITAL = MAX(max_mem_used_DECORBITAL,mem_allocated_DECORBITAL)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_DECORBITAL'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_DECORBITAL'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_DECORBITAL')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_DECORBITAL

subroutine mem_deallocated_mem_DECORBITAL(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_DECORBITAL = mem_tp_allocated_DECORBITAL - nsize
      if (mem_tp_allocated_DECORBITAL < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_DECORBITAL - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_DECORBITAL - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_DECORBITAL = mem_allocated_DECORBITAL - nsize
      if (mem_allocated_DECORBITAL < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_DECORBITAL - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      mem_allocated_global = mem_allocated_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_DECORBITAL - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_DECORBITAL





subroutine mem_allocated_mem_DECFRAG(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_DECFRAG = mem_tp_allocated_DECFRAG + nsize
      max_mem_tp_used_DECFRAG = MAX(max_mem_tp_used_DECFRAG,mem_tp_allocated_DECFRAG)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_DECFRAG = mem_allocated_DECFRAG + nsize
      max_mem_used_DECFRAG = MAX(max_mem_used_DECFRAG,mem_allocated_DECFRAG)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_DECFRAG'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_DECFRAG'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_DECFRAG')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_DECFRAG

subroutine mem_deallocated_mem_DECFRAG(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_DECFRAG = mem_tp_allocated_DECFRAG - nsize
      if (mem_tp_allocated_DECFRAG < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_DECFRAG - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_DECFRAG - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_DECFRAG = mem_allocated_DECFRAG - nsize
      if (mem_allocated_DECFRAG < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_DECFRAG - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_DECFRAG - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_DECFRAG


subroutine mem_allocated_mem_BATCHTOORB(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_BATCHTOORB = mem_tp_allocated_BATCHTOORB + nsize
      max_mem_tp_used_BATCHTOORB = MAX(max_mem_tp_used_BATCHTOORB,mem_tp_allocated_BATCHTOORB)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_BATCHTOORB = mem_allocated_BATCHTOORB + nsize
      max_mem_used_BATCHTOORB = MAX(max_mem_used_BATCHTOORB,mem_allocated_BATCHTOORB)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_BATCHTOORB'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_BATCHTOORB'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_BATCHTOORB')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_BATCHTOORB

subroutine mem_deallocated_mem_BATCHTOORB(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_BATCHTOORB = mem_tp_allocated_BATCHTOORB - nsize
      if (mem_tp_allocated_BATCHTOORB < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_BATCHTOORB - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_BATCHTOORB - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_BATCHTOORB = mem_allocated_BATCHTOORB - nsize
      if (mem_allocated_BATCHTOORB < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_BATCHTOORB - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_BATCHTOORB - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_BATCHTOORB

subroutine mem_allocated_mem_DECAOBATCHINFO(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_DECAOBATCHINFO = mem_tp_allocated_DECAOBATCHINFO + nsize
      max_mem_tp_used_DECAOBATCHINFO = MAX(max_mem_tp_used_DECAOBATCHINFO,mem_tp_allocated_DECAOBATCHINFO)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_DECAOBATCHINFO = mem_allocated_DECAOBATCHINFO + nsize
      max_mem_used_DECAOBATCHINFO = MAX(max_mem_used_DECAOBATCHINFO,mem_allocated_DECAOBATCHINFO)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_DECAOBATCHINFO'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_DECAOBATCHINFO'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_DECAOBATCHINFO')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_DECAOBATCHINFO

subroutine mem_deallocated_mem_DECAOBATCHINFO(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_DECAOBATCHINFO = mem_tp_allocated_DECAOBATCHINFO - nsize
      if (mem_tp_allocated_DECAOBATCHINFO < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_DECAOBATCHINFO - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_DECAOBATCHINFO - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_DECAOBATCHINFO = mem_allocated_DECAOBATCHINFO - nsize
      if (mem_allocated_DECAOBATCHINFO < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_DECAOBATCHINFO - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_DECAOBATCHINFO - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_DECAOBATCHINFO


subroutine mem_allocated_mem_MYPOINTER(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_MYPOINTER = mem_tp_allocated_MYPOINTER + nsize
      max_mem_tp_used_MYPOINTER = MAX(max_mem_tp_used_MYPOINTER,mem_tp_allocated_MYPOINTER)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_MYPOINTER = mem_allocated_MYPOINTER + nsize
      max_mem_used_MYPOINTER = MAX(max_mem_used_MYPOINTER,mem_allocated_MYPOINTER)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_MYPOINTER'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_MYPOINTER'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_MYPOINTER')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_MYPOINTER

subroutine mem_deallocated_mem_MYPOINTER(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_MYPOINTER = mem_tp_allocated_MYPOINTER - nsize
      if (mem_tp_allocated_MYPOINTER < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_MYPOINTER - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_MYPOINTER - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_MYPOINTER = mem_allocated_MYPOINTER - nsize
      if (mem_allocated_MYPOINTER < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_MYPOINTER - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_MYPOINTER - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_MYPOINTER



subroutine mem_allocated_mem_fragmentAOS(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_fragmentAOS = mem_tp_allocated_fragmentAOS + nsize
      max_mem_tp_used_fragmentAOS = MAX(max_mem_tp_used_fragmentAOS,mem_tp_allocated_fragmentAOS)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_fragmentAOS = mem_allocated_fragmentAOS + nsize
      max_mem_used_fragmentAOS = MAX(max_mem_used_fragmentAOS,mem_allocated_fragmentAOS)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_fragmentAOS'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_fragmentAOS'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_fragmentAOS')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_fragmentAOS

subroutine mem_deallocated_mem_fragmentAOS(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_fragmentAOS = mem_tp_allocated_fragmentAOS - nsize
      if (mem_tp_allocated_fragmentAOS < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_fragmentAOS - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_fragmentAOS - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_fragmentAOS = mem_allocated_fragmentAOS - nsize
      if (mem_allocated_fragmentAOS < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_fragmentAOS - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_fragmentAOS - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_fragmentAOS



subroutine mem_allocated_mem_ARRAY2(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ARRAY2 = mem_tp_allocated_ARRAY2 + nsize
      max_mem_tp_used_ARRAY2 = MAX(max_mem_tp_used_ARRAY2,mem_tp_allocated_ARRAY2)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_ARRAY2 = mem_allocated_ARRAY2 + nsize
      max_mem_used_ARRAY2 = MAX(max_mem_used_ARRAY2,mem_allocated_ARRAY2)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_ARRAY2'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_ARRAY2'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_ARRAY2')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
   ENDIF
end subroutine mem_allocated_mem_ARRAY2

subroutine mem_deallocated_mem_ARRAY2(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ARRAY2 = mem_tp_allocated_ARRAY2 - nsize
      if (mem_tp_allocated_ARRAY2 < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ARRAY2 - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ARRAY2 - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_ARRAY2 = mem_allocated_ARRAY2 - nsize
      if (mem_allocated_ARRAY2 < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ARRAY2 - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ARRAY2 - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_ARRAY2


subroutine mem_allocated_mem_ARRAY4(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ARRAY4 = mem_tp_allocated_ARRAY4 + nsize
      max_mem_tp_used_ARRAY4 = MAX(max_mem_tp_used_ARRAY4,mem_tp_allocated_ARRAY4)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_ARRAY4 = mem_allocated_ARRAY4 + nsize
      max_mem_used_ARRAY4 = MAX(max_mem_used_ARRAY4,mem_allocated_ARRAY4)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_ARRAY4'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_ARRAY4'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_ARRAY4')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_ARRAY4

subroutine mem_deallocated_mem_ARRAY4(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ARRAY4 = mem_tp_allocated_ARRAY4 - nsize
      if (mem_tp_allocated_ARRAY4 < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ARRAY4 - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ARRAY4 - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_ARRAY4 = mem_allocated_ARRAY4 - nsize
      if (mem_allocated_ARRAY4 < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ARRAY4 - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ARRAY4 - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_ARRAY4

subroutine mem_allocated_mem_ARRAY(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ARRAY = mem_tp_allocated_ARRAY + nsize
      max_mem_tp_used_ARRAY = MAX(max_mem_tp_used_ARRAY,mem_tp_allocated_ARRAY)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_ARRAY = mem_allocated_ARRAY + nsize
      max_mem_used_ARRAY = MAX(max_mem_used_ARRAY,mem_allocated_ARRAY)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_ARRAY'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_ARRAY'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_ARRAY')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_ARRAY

subroutine mem_deallocated_mem_ARRAY(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ARRAY = mem_tp_allocated_ARRAY - nsize
      if (mem_tp_allocated_ARRAY < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ARRAY - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ARRAY - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_ARRAY = mem_allocated_ARRAY - nsize
      if (mem_allocated_ARRAY < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ARRAY - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ARRAY - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_ARRAY

  subroutine mem_allocated_mem_PNOSpaceInfo(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
     IF(mem_InsideOMPsection)THEN!we add to thread private variables
        mem_tp_allocated_PNOSpaceInfo = mem_tp_allocated_PNOSpaceInfo + nsize
        max_mem_tp_used_PNOSpaceInfo = MAX(max_mem_tp_used_PNOSpaceInfo,mem_tp_allocated_PNOSpaceInfo)
        !Count also the total memory:
        mem_tp_allocated_global = mem_tp_allocated_global  + nsize
        mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
        max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
        max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
     ELSE
        mem_allocated_PNOSpaceInfo = mem_allocated_PNOSpaceInfo + nsize
        max_mem_used_PNOSpaceInfo = MAX(max_mem_used_PNOSpaceInfo,mem_allocated_PNOSpaceInfo)
        !Count also the total memory:
        mem_allocated_global = mem_allocated_global  + nsize
        mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
        max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
        max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
     ENDIF
   end subroutine mem_allocated_mem_PNOSpaceInfo

   subroutine mem_deallocated_mem_PNOSpaceInfo(nsize)
     implicit none
     integer (kind=long), intent(in) :: nsize
     IF(mem_InsideOMPsection)THEN!we add to thread private variables
        mem_tp_allocated_PNOSpaceInfo = mem_tp_allocated_PNOSpaceInfo - nsize
        if (mem_tp_allocated_PNOSpaceInfo < 0) then
           call memory_error_quit('Error in mem_tp_deallocated_mem_tp_PNOSpaceInfo - probably integer overflow!',nsize)
        endif
        !Count also the total memory:
        mem_tp_allocated_global = mem_tp_allocated_global - nsize
        mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
        if (mem_tp_allocated_global < 0) then
           call memory_error_quit('Error in mem_tp_deallocated_mem_tp_PNOSpaceInfo - probably integer overflow!',nsize)
        endif
     ELSE
        mem_allocated_PNOSpaceInfo = mem_allocated_PNOSpaceInfo - nsize
        if (mem_allocated_PNOSpaceInfo < 0) then
           call memory_error_quit('Error in mem_deallocated_mem_PNOSpaceInfo - probably integer overflow!',nsize)
        endif
        !Count also the total memory:
        mem_allocated_global = mem_allocated_global - nsize
        mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
        if (mem_allocated_global < 0) then
           call memory_error_quit('Error in mem_deallocated_mem_PNOSpaceInfo - probably integer overflow!',nsize)
        endif
     ENDIF
   end subroutine mem_deallocated_mem_PNOSpaceInfo
subroutine mem_allocated_mem_MP2DENS(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_MP2DENS = mem_tp_allocated_MP2DENS + nsize
      max_mem_tp_used_MP2DENS = MAX(max_mem_tp_used_MP2DENS,mem_tp_allocated_MP2DENS)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_MP2DENS = mem_allocated_MP2DENS + nsize
      max_mem_used_MP2DENS = MAX(max_mem_used_MP2DENS,mem_allocated_MP2DENS)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_MP2DENS'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_MP2DENS'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_MP2DENS')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_MP2DENS

subroutine mem_deallocated_mem_MP2DENS(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_MP2DENS = mem_tp_allocated_MP2DENS - nsize
      if (mem_tp_allocated_MP2DENS < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_MP2DENS - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_MP2DENS - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_MP2DENS = mem_allocated_MP2DENS - nsize
      if (mem_allocated_MP2DENS < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_MP2DENS - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_MP2DENS - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_MP2DENS



subroutine mem_allocated_mem_TRACEBACK(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_TRACEBACK = mem_tp_allocated_TRACEBACK + nsize
      max_mem_tp_used_TRACEBACK = MAX(max_mem_tp_used_TRACEBACK,mem_tp_allocated_TRACEBACK)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_TRACEBACK = mem_allocated_TRACEBACK + nsize
      max_mem_used_TRACEBACK = MAX(max_mem_used_TRACEBACK,mem_allocated_TRACEBACK)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_TRACEBACK'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_TRACEBACK'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_TRACEBACK')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_TRACEBACK

subroutine mem_deallocated_mem_TRACEBACK(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_TRACEBACK = mem_tp_allocated_TRACEBACK - nsize
      if (mem_tp_allocated_TRACEBACK < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_TRACEBACK - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_TRACEBACK - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_TRACEBACK = mem_allocated_TRACEBACK - nsize
      if (mem_allocated_TRACEBACK < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_TRACEBACK - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_TRACEBACK - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_TRACEBACK


subroutine mem_allocated_mem_MP2GRAD(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_MP2GRAD = mem_tp_allocated_MP2GRAD + nsize
      max_mem_tp_used_MP2GRAD = MAX(max_mem_tp_used_MP2GRAD,mem_tp_allocated_MP2GRAD)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_MP2GRAD = mem_allocated_MP2GRAD + nsize
      max_mem_used_MP2GRAD = MAX(max_mem_used_MP2GRAD,mem_allocated_MP2GRAD)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_MP2GRAD'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_MP2GRAD'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_MP2GRAD')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_MP2GRAD

subroutine mem_deallocated_mem_MP2GRAD(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_MP2GRAD = mem_tp_allocated_MP2GRAD - nsize
      if (mem_tp_allocated_MP2GRAD < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_MP2GRAD - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_MP2GRAD - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_MP2GRAD = mem_allocated_MP2GRAD - nsize
      if (mem_allocated_MP2GRAD < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_MP2GRAD - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_MP2GRAD - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_MP2GRAD



subroutine mem_allocated_mem_ODBATCH(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ODBATCH = mem_tp_allocated_ODBATCH + nsize
      max_mem_tp_used_ODBATCH = MAX(max_mem_tp_used_ODBATCH,mem_tp_allocated_ODBATCH)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_ODBATCH = mem_allocated_ODBATCH + nsize
      max_mem_used_ODBATCH = MAX(max_mem_used_ODBATCH,mem_allocated_ODBATCH)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_ODBATCH'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_ODBATCH'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_ODBATCH')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_ODBATCH

subroutine mem_deallocated_mem_ODBATCH(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ODBATCH = mem_tp_allocated_ODBATCH - nsize
      if (mem_tp_allocated_ODBATCH < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ODBATCH - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ODBATCH - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_ODBATCH = mem_allocated_ODBATCH - nsize
      if (mem_allocated_ODBATCH < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ODBATCH - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ODBATCH - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_ODBATCH

subroutine mem_allocated_mem_LSAOTENSOR(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_LSAOTENSOR = mem_tp_allocated_LSAOTENSOR + nsize
      max_mem_tp_used_LSAOTENSOR = MAX(max_mem_tp_used_LSAOTENSOR,mem_tp_allocated_LSAOTENSOR)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_LSAOTENSOR = mem_allocated_LSAOTENSOR + nsize
      max_mem_used_LSAOTENSOR = MAX(max_mem_used_LSAOTENSOR,mem_allocated_LSAOTENSOR)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_LSAOTENSOR'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_LSAOTENSOR'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_LSAOTENSOR')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_LSAOTENSOR

subroutine mem_deallocated_mem_LSAOTENSOR(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_LSAOTENSOR = mem_tp_allocated_LSAOTENSOR - nsize
      if (mem_tp_allocated_LSAOTENSOR < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_LSAOTENSOR - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_LSAOTENSOR - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_LSAOTENSOR = mem_allocated_LSAOTENSOR - nsize
      if (mem_allocated_LSAOTENSOR < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_LSAOTENSOR - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_LSAOTENSOR - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_LSAOTENSOR

subroutine mem_allocated_mem_SLSAOTENSOR(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_SLSAOTENSOR = mem_tp_allocated_SLSAOTENSOR + nsize
      max_mem_tp_used_SLSAOTENSOR = MAX(max_mem_tp_used_SLSAOTENSOR,mem_tp_allocated_SLSAOTENSOR)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_SLSAOTENSOR = mem_allocated_SLSAOTENSOR + nsize
      max_mem_used_SLSAOTENSOR = MAX(max_mem_used_SLSAOTENSOR,mem_allocated_SLSAOTENSOR)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_SLSAOTENSOR'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_SLSAOTENSOR'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_SLSAOTENSOR')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_SLSAOTENSOR

subroutine mem_deallocated_mem_SLSAOTENSOR(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_SLSAOTENSOR = mem_tp_allocated_SLSAOTENSOR - nsize
      if (mem_tp_allocated_SLSAOTENSOR < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_SLSAOTENSOR - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_SLSAOTENSOR - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_SLSAOTENSOR = mem_allocated_SLSAOTENSOR - nsize
      if (mem_allocated_SLSAOTENSOR < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_SLSAOTENSOR - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_SLSAOTENSOR - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_SLSAOTENSOR

subroutine mem_allocated_mem_GLOBALLSAOTENSOR(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_GLOBALLSAOTENSOR = mem_tp_allocated_GLOBALLSAOTENSOR + nsize
      max_mem_tp_used_GLOBALLSAOTENSOR = MAX(max_mem_tp_used_GLOBALLSAOTENSOR,mem_tp_allocated_GLOBALLSAOTENSOR)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_GLOBALLSAOTENSOR = mem_allocated_GLOBALLSAOTENSOR + nsize
      max_mem_used_GLOBALLSAOTENSOR = MAX(max_mem_used_GLOBALLSAOTENSOR,mem_allocated_GLOBALLSAOTENSOR)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_GLOBALLSAOTENSOR'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_GLOBALLSAOTENSOR'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_GLOBALLSAOTENSOR')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_GLOBALLSAOTENSOR

subroutine mem_deallocated_mem_GLOBALLSAOTENSOR(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_GLOBALLSAOTENSOR = mem_tp_allocated_GLOBALLSAOTENSOR - nsize
      if (mem_tp_allocated_GLOBALLSAOTENSOR < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_GLOBALLSAOTENSOR - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_GLOBALLSAOTENSOR - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_GLOBALLSAOTENSOR = mem_allocated_GLOBALLSAOTENSOR - nsize
      if (mem_allocated_GLOBALLSAOTENSOR < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_GLOBALLSAOTENSOR - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_GLOBALLSAOTENSOR - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_GLOBALLSAOTENSOR

subroutine mem_allocated_mem_ATOMTYPEITEM(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ATOMTYPEITEM = mem_tp_allocated_ATOMTYPEITEM + nsize
      max_mem_tp_used_ATOMTYPEITEM = MAX(max_mem_tp_used_ATOMTYPEITEM,mem_tp_allocated_ATOMTYPEITEM)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_ATOMTYPEITEM = mem_allocated_ATOMTYPEITEM + nsize
      max_mem_used_ATOMTYPEITEM = MAX(max_mem_used_ATOMTYPEITEM,mem_allocated_ATOMTYPEITEM)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_ATOMTYPEITEM'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_ATOMTYPEITEM'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_ATOMTYPEITEM')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_ATOMTYPEITEM

subroutine mem_deallocated_mem_ATOMTYPEITEM(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ATOMTYPEITEM = mem_tp_allocated_ATOMTYPEITEM - nsize
      if (mem_tp_allocated_ATOMTYPEITEM < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ATOMTYPEITEM - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ATOMTYPEITEM - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_ATOMTYPEITEM = mem_allocated_ATOMTYPEITEM - nsize
      if (mem_allocated_ATOMTYPEITEM < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ATOMTYPEITEM - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ATOMTYPEITEM - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_ATOMTYPEITEM

subroutine mem_allocated_mem_ATOMITEM(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ATOMITEM = mem_tp_allocated_ATOMITEM + nsize
      max_mem_tp_used_ATOMITEM = MAX(max_mem_tp_used_ATOMITEM,mem_tp_allocated_ATOMITEM)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_ATOMITEM = mem_allocated_ATOMITEM + nsize
      max_mem_used_ATOMITEM = MAX(max_mem_used_ATOMITEM,mem_allocated_ATOMITEM)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_ATOMITEM'
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_ATOMITEM'
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_ATOMITEM')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_ATOMITEM

subroutine mem_deallocated_mem_ATOMITEM(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ATOMITEM = mem_tp_allocated_ATOMITEM - nsize
      if (mem_tp_allocated_ATOMITEM < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ATOMITEM - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ATOMITEM - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_ATOMITEM = mem_allocated_ATOMITEM - nsize
      if (mem_allocated_ATOMITEM < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ATOMITEM - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ATOMITEM - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_ATOMITEM

subroutine mem_allocated_mem_LSMATRIX(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_LSMATRIX = mem_tp_allocated_LSMATRIX + nsize
      max_mem_tp_used_LSMATRIX = MAX(max_mem_tp_used_LSMATRIX,mem_tp_allocated_LSMATRIX)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_LSMATRIX = mem_allocated_LSMATRIX + nsize
      max_mem_used_LSMATRIX = MAX(max_mem_used_LSMATRIX,mem_allocated_LSMATRIX)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_LSMATRIX '
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_LSMATRIX '
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_LSMATRIX')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_LSMATRIX

subroutine mem_deallocated_mem_LSMATRIX(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_LSMATRIX = mem_tp_allocated_LSMATRIX - nsize
      if (mem_tp_allocated_LSMATRIX < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_LSMATRIX - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_LSMATRIX - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_LSMATRIX = mem_allocated_LSMATRIX - nsize
      if (mem_allocated_LSMATRIX < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_LSMATRIX - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_LSMATRIX - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_LSMATRIX

!==================================================================================
!
!    The next subroutines are special subroutines to keep track of memory used
!    in special structures
!
!==================================================================================

! to keep track of memory used in the overlap structure
subroutine mem_allocated_mem_overlap(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_overlap = mem_tp_allocated_overlap + nsize
      max_mem_tp_used_overlap = MAX(max_mem_tp_used_overlap,mem_tp_allocated_overlap)
   ELSE
      mem_allocated_overlap = mem_allocated_overlap + nsize
      max_mem_used_overlap = MAX(max_mem_used_overlap,mem_allocated_overlap)
   ENDIF
end subroutine mem_allocated_mem_overlap

subroutine mem_deallocated_mem_overlap(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_overlap = mem_tp_allocated_overlap - nsize
      if (mem_tp_allocated_overlap < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_overlap - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_overlap = mem_allocated_overlap - nsize
      if (mem_allocated_overlap < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_overlap - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_overlap

! to keep track of memory used in the IntWork structure
subroutine mem_allocated_mem_IntWork(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_IntWork = mem_tp_allocated_IntWork + nsize
      max_mem_tp_used_IntWork = MAX(max_mem_tp_used_IntWork,mem_tp_allocated_IntWork)
   ELSE
      mem_allocated_IntWork = mem_allocated_IntWork + nsize
      max_mem_used_IntWork = MAX(max_mem_used_IntWork,mem_allocated_IntWork)
   ENDIF
   call mem_allocated_mem_real(nsize)
end subroutine mem_allocated_mem_IntWork

subroutine mem_deallocated_mem_IntWork(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_IntWork = mem_tp_allocated_IntWork - nsize
      if (mem_tp_allocated_IntWork < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_IntWork - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_IntWork = mem_allocated_IntWork - nsize
      if (mem_allocated_IntWork < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_IntWork - probably integer overflow!',nsize)
      endif
   ENDIF
   call mem_deallocated_mem_real(nsize)
end subroutine mem_deallocated_mem_IntWork

!2 to keep track of memory used in the linkshell structure
subroutine mem_allocated_mem_linkshell(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_linkshell = mem_tp_allocated_linkshell + nsize
      max_mem_tp_used_linkshell = MAX(max_mem_tp_used_linkshell,mem_tp_allocated_linkshell)
   ELSE
      mem_allocated_linkshell = mem_allocated_linkshell + nsize
      max_mem_used_linkshell = MAX(max_mem_used_linkshell,mem_allocated_linkshell)
   ENDIF
end subroutine mem_allocated_mem_linkshell

subroutine mem_deallocated_mem_linkshell(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_linkshell = mem_tp_allocated_linkshell - nsize
      if (mem_tp_allocated_linkshell < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_linkshell - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_linkshell = mem_allocated_linkshell - nsize
      if (mem_allocated_linkshell < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_linkshell - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_linkshell

subroutine set_mem_XCcalc(mem_allocated_dft,max_mem_used_dft)
   implicit none
   integer(KIND=long) :: mem_allocated_dft,max_mem_used_dft
   mem_allocated_XCcalc = mem_allocated_dft
   max_mem_used_XCcalc = max_mem_used_dft
end subroutine set_mem_XCcalc

!3 to keep track of memory used in the integralitem structure
subroutine mem_allocated_mem_integralitem(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_integralitem = mem_tp_allocated_integralitem + nsize
      max_mem_tp_used_integralitem = MAX(max_mem_tp_used_integralitem,mem_tp_allocated_integralitem)
   ELSE
      mem_allocated_integralitem = mem_allocated_integralitem + nsize
      max_mem_used_integralitem = MAX(max_mem_used_integralitem,mem_allocated_integralitem)
   ENDIF
end subroutine mem_allocated_mem_integralitem

subroutine mem_deallocated_mem_integralitem(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_integralitem = mem_tp_allocated_integralitem - nsize
      if (mem_tp_allocated_integralitem < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_integralitem - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_integralitem = mem_allocated_integralitem - nsize
      if (mem_allocated_integralitem < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_integralitem - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_integralitem
!4 to keep track of memory used in the integrand structure
subroutine mem_allocated_mem_integrand(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_integrand = mem_tp_allocated_integrand + nsize
      max_mem_tp_used_integrand = MAX(max_mem_tp_used_integrand,mem_tp_allocated_integrand)
   ELSE
      mem_allocated_integrand = mem_allocated_integrand + nsize
      max_mem_used_integrand = MAX(max_mem_used_integrand,mem_allocated_integrand)
   ENDIF
end subroutine mem_allocated_mem_integrand

subroutine mem_deallocated_mem_integrand(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_integrand = mem_tp_allocated_integrand - nsize
      if (mem_tp_allocated_integrand < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_integrand - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_integrand = mem_allocated_integrand - nsize
      if (mem_allocated_integrand < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_integrand - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_integrand
!7 to keep track of memory used in the ODitem structure
subroutine mem_allocated_mem_ODitem(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ODitem = mem_tp_allocated_ODitem + nsize
      max_mem_tp_used_ODitem = MAX(max_mem_tp_used_ODitem,mem_tp_allocated_ODitem)
   ELSE
      mem_allocated_ODitem = mem_allocated_ODitem + nsize
      max_mem_used_ODitem = MAX(max_mem_used_ODitem,mem_allocated_ODitem)
   ENDIF
end subroutine mem_allocated_mem_ODitem

subroutine mem_deallocated_mem_ODitem(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_ODitem = mem_tp_allocated_ODitem - nsize
      if (mem_tp_allocated_ODitem < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_ODitem - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_ODitem = mem_allocated_ODitem - nsize
      if (mem_allocated_ODitem < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_ODitem - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_ODitem
!9 to keep track of memory used in the FMM structure
subroutine mem_allocated_mem_FMM(nsize,global)
   implicit none
   logical :: global
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_FMM = mem_tp_allocated_FMM + nsize
      max_mem_tp_used_FMM = MAX(max_mem_tp_used_FMM,mem_tp_allocated_FMM)
      IF(global)THEN
         !Count also the total memory:
         mem_tp_allocated_global = mem_tp_allocated_global  + nsize
         mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
         if (mem_tp_allocated_global < 0) then
            write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
            call memory_error_quit('Error in mem_tp_allocated_mem_tp_real - probably integer overflow!',nsize)
         endif
         max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
         max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
      ENDIF
   ELSE
      mem_allocated_FMM = mem_allocated_FMM + nsize
      max_mem_used_FMM = MAX(max_mem_used_FMM,mem_allocated_FMM)
      IF(global)THEN
         !Count also the total memory:
         mem_allocated_global = mem_allocated_global  + nsize
         mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
         if (mem_allocated_global < 0) then
            write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
            call memory_error_quit('Error in mem_allocated_mem_real - probably integer overflow!',nsize)
         endif
         IF(MemModParamPrintMemory)THEN
            IF(mem_allocated_global.GT.max_mem_used_global)THEN              
               IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
                  WRITE(MemModParamPrintMemorylupri,'(A)')&
                     & 'Increased Maximum Memory Usage in mem_allocated_mem_FMM '
                  call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
                  WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_FMM '
                  call print_mem_alloc(6,mem_allocated_global,'TOTAL')
                  call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_FMM')
               ENDIF
            ENDIF
         ENDIF
         max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
         max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
      ENDIF
   ENDIF
end subroutine mem_allocated_mem_FMM

subroutine mem_deallocated_mem_FMM(nsize,global)
   implicit none
   logical :: global
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_FMM = mem_tp_allocated_FMM - nsize
      if (mem_tp_allocated_FMM < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_FMM - probably integer overflow!',nsize)
      endif
      IF(global)THEN
         !Count also the total memory:
         mem_tp_allocated_global = mem_tp_allocated_global - nsize
         mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
         if (mem_tp_allocated_global < 0) then
            call memory_error_quit('Error in mem_tp_deallocated_mem_tp_logical - probably integer overflow!',nsize)
         endif
      ENDIF
   ELSE
      mem_allocated_FMM = mem_allocated_FMM - nsize
      if (mem_allocated_FMM < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_FMM - probably integer overflow!',nsize)
      endif
      IF(global)THEN
         !Count also the total memory:
         mem_allocated_global = mem_allocated_global - nsize
         mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
         if (mem_allocated_global < 0) then
            call memory_error_quit('Error in mem_deallocated_mem_logical - probably integer overflow!',nsize)
         endif
      ENDIF
   ENDIF

end subroutine mem_deallocated_mem_FMM
!10 to keep track of memory used in the lstensor structure
subroutine mem_allocated_mem_lstensor(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_lstensor = mem_tp_allocated_lstensor + nsize
      max_mem_tp_used_lstensor = MAX(max_mem_tp_used_lstensor,mem_tp_allocated_lstensor)
   ELSE
      mem_allocated_lstensor = mem_allocated_lstensor + nsize
      max_mem_used_lstensor = MAX(max_mem_used_lstensor,mem_allocated_lstensor)
   ENDIF
end subroutine mem_allocated_mem_lstensor

subroutine mem_deallocated_mem_lstensor(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_lstensor = mem_tp_allocated_lstensor - nsize
      if (mem_tp_allocated_lstensor < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_lstensor - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_lstensor = mem_allocated_lstensor - nsize
      if (mem_allocated_lstensor < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_lstensor - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_lstensor

!11 to keep track of memory used in the type_matrix structure
subroutine mem_allocated_mem_type_matrix(nsize,nsizeFULL)
   implicit none
   integer (kind=long), intent(in) :: nsize
   integer (kind=long), optional,intent(in) :: nsizeFULL
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_type_matrix = mem_tp_allocated_type_matrix + nsize
      max_mem_tp_used_type_matrix = MAX(max_mem_tp_used_type_matrix,mem_tp_allocated_type_matrix)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_allocated_mem_type_matrix - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      mem_allocated_type_matrix = mem_allocated_type_matrix + nsize
      max_mem_used_type_matrix = MAX(max_mem_used_type_matrix,mem_allocated_type_matrix)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_allocated_mem_type_matrix - probably integer overflow!',nsize)
      endif
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_MATRIX '
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_MATRIX '
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_MATRIX')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
   IF(present(nsizeFULL))THEN
#ifdef VAR_SCALAPACK
      IF(mem_InsideOMPsection)THEN!we add to thread private variables
         mem_tp_allocated_type_matrix_MPIFULL = mem_tp_allocated_type_matrix_MPIFULL + nsizeFULL
         max_mem_tp_used_type_matrix_MPIFULL = MAX(max_mem_tp_used_type_matrix_MPIFULL,mem_tp_allocated_type_matrix_MPIFULL)
      ELSE
         mem_allocated_type_matrix_MPIFULL = mem_allocated_type_matrix_MPIFULL + nsizeFULL
         max_mem_used_type_matrix_MPIFULL = MAX(max_mem_used_type_matrix_MPIFULL,mem_allocated_type_matrix_MPIFULL)
      ENDIF
#endif
   ENDIF
end subroutine mem_allocated_mem_type_matrix

subroutine mem_deallocated_mem_type_matrix(nsize,nsizeFULL)
   implicit none
   integer (kind=long), intent(in) :: nsize
   integer (kind=long), optional,intent(in) :: nsizeFULL
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_type_matrix = mem_tp_allocated_type_matrix - nsize
      if (mem_tp_allocated_type_matrix < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_type_matrix - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_type_matrix - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_type_matrix = mem_allocated_type_matrix - nsize
      if (mem_allocated_type_matrix < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_type_matrix - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_logical - probably integer overflow!',nsize)
      endif
   ENDIF
   IF(present(nsizeFULL))THEN
#ifdef VAR_SCALAPACK
      IF(mem_InsideOMPsection)THEN!we add to thread private variables
         mem_tp_allocated_type_matrix_MPIFULL = mem_tp_allocated_type_matrix_MPIFULL - nsizeFULL
      ELSE
         mem_allocated_type_matrix_MPIFULL = mem_allocated_type_matrix_MPIFULL - nsizeFULL
      ENDIF
#endif
   ENDIF

end subroutine mem_deallocated_mem_type_matrix

#ifdef MOD_UNRELEASED
!13 to keep track of memory used in the lvec_data structure
subroutine mem_allocated_mem_lvec_data(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_lvec_data = mem_tp_allocated_lvec_data + nsize
      max_mem_tp_used_lvec_data = MAX(max_mem_tp_used_lvec_data,mem_tp_allocated_lvec_data)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_allocated_mem_tp_real - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      !        CALL print_maxmem(6,nsize,'add Lvec_data')
      mem_allocated_lvec_data = mem_allocated_lvec_data + nsize
      max_mem_used_lvec_data = MAX(max_mem_used_lvec_data,mem_allocated_lvec_data)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_allocated_mem_real - probably integer overflow!',nsize)
      endif
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_lvec_data '
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_lvec_data '
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_lvec_data')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_lvec_data

subroutine mem_deallocated_mem_lvec_data(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_lvec_data = mem_tp_allocated_lvec_data - nsize
      if (mem_tp_allocated_lvec_data < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_lvec_data - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_logical - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_lvec_data = mem_allocated_lvec_data - nsize
      if (mem_allocated_lvec_data < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_lvec_data - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_logical - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_lvec_data
!14 to keep track of memory used in the lattice_cell_info structure
subroutine mem_allocated_mem_lattice_cell(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_lattice_cell = mem_tp_allocated_lattice_cell + nsize
      max_mem_tp_used_lattice_cell = MAX(max_mem_tp_used_lattice_cell,mem_tp_allocated_lattice_cell)
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global  + nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global  + nsize
      if (mem_tp_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_tp_allocated_global =', mem_tp_allocated_global
         call memory_error_quit('Error in mem_tp_allocated_mem_tp_real - probably integer overflow!',nsize)
      endif
      max_mem_tp_used_global = MAX(max_mem_tp_used_global,mem_tp_allocated_global)
      max_mem_tp_used_Heap_global = MAX(max_mem_tp_used_Heap_global,mem_tp_allocated_Heap_global)
   ELSE
      !        CALL print_maxmem(6,nsize,'add Lattice_cell')
      mem_allocated_lattice_cell = mem_allocated_lattice_cell + nsize
      max_mem_used_lattice_cell = MAX(max_mem_used_lattice_cell,mem_allocated_lattice_cell)
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global  + nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global  + nsize
      if (mem_allocated_global < 0) then
         write(*,*) 'Total memory negative! mem_allocated_global =', mem_allocated_global
         call memory_error_quit('Error in mem_allocated_mem_real - probably integer overflow!',nsize)
      endif
      IF(MemModParamPrintMemory)THEN
         IF(mem_allocated_global.GT.max_mem_used_global)THEN              
            IF(max_mem_used_global.GT.PrintMemoryLowerLimit)THEN
               WRITE(MemModParamPrintMemorylupri,'(A)')&
                  & 'Increased Maximum Memory Usage in mem_allocated_mem_lattice_cell '
               call print_mem_alloc(MemModParamPrintMemorylupri,mem_allocated_global,'TOTAL')
               WRITE(*,'(A)')'Increased Maximum Memory Usage in mem_allocated_mem_lattice_cell '
               call print_mem_alloc(6,mem_allocated_global,'TOTAL')
               call LsTraceBack('Increased Maximum Memory Usage in mem_allocated_mem_lattice_cell')
            ENDIF
         ENDIF
      ENDIF
      max_mem_used_global = MAX(max_mem_used_global,mem_allocated_global)
      max_mem_used_Heap_global = MAX(max_mem_used_Heap_global,mem_allocated_Heap_global)
   ENDIF
end subroutine mem_allocated_mem_lattice_cell

subroutine mem_deallocated_mem_lattice_cell(nsize)
   implicit none
   integer (kind=long), intent(in) :: nsize
   IF(mem_InsideOMPsection)THEN!we add to thread private variables
      mem_tp_allocated_lattice_cell = mem_tp_allocated_lattice_cell - nsize
      if (mem_tp_allocated_lattice_cell < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_lattice_cell - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_tp_allocated_global = mem_tp_allocated_global - nsize
      mem_tp_allocated_Heap_global = mem_tp_allocated_Heap_global - nsize
      if (mem_tp_allocated_global < 0) then
         call memory_error_quit('Error in mem_tp_deallocated_mem_tp_logical - probably integer overflow!',nsize)
      endif
   ELSE
      mem_allocated_lattice_cell = mem_allocated_lattice_cell - nsize
      if (mem_allocated_lattice_cell < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_lattice_cell - probably integer overflow!',nsize)
      endif
      !Count also the total memory:
      mem_allocated_global = mem_allocated_global - nsize
      mem_allocated_Heap_global = mem_allocated_Heap_global - nsize
      if (mem_allocated_global < 0) then
         call memory_error_quit('Error in mem_deallocated_mem_logical - probably integer overflow!',nsize)
      endif
   ENDIF
end subroutine mem_deallocated_mem_lattice_cell



#endif


!> \brief Get how much memory is currently available.
!> Note: This uses a system call and checks whether Linux or Max
!> is used. If neither of these two systems are used, this routine
!> will return MemFound=.FALSE. and the MemoryAvailable will be set to
!> 1 GB be default.
!> \author Kasper Kristensen, modified by Patrick Ettenhuber
!> \date January 2012
subroutine get_available_memory(lupri,MemoryAvailable,memfound,suppress_print)
   !> Logical unit number for output file
   integer,intent(in) :: lupri
   !> Available memory measured in GB (using 1GB = 1000000000 bytes)
   real(realk),intent(inout) :: MemoryAvailable
   !> Was the memory information found?
   logical,intent(inout) :: MemFound
   logical,intent(in),optional :: suppress_print
   !> check for /proc/meminfo
   logical :: meminfo_found,sp

   ! If System will not be identified, use default values
   memfound=.false.
   MemoryAvailable=1.0E0_realk

   sp = .false.
   if(present(suppress_print))sp = suppress_print

   INQUIRE(FILE="/proc/meminfo", EXIST=meminfo_found)

   ! Is this a LINUX system?
   if(meminfo_found) then
      call get_available_memory_specific('MENFO',MemoryAvailable,memfound)
      if(memfound.and..not.sp) then
         write(lupri,*) 'get_available_memory: System identified to be LINUX!'
         write(lupri,'(1X,a,g16.5)') 'Available Memory (GB) = ', MemoryAvailable
      end if
   endif

   ! meminfo not found or mem in meminfo not found  --> Is this a MAC system?
   if((.not. meminfo_found).or.(.not.memfound)) then
      call get_available_memory_specific('MAC  ',MemoryAvailable,memfound)
      if(memfound.and..not.sp) then
         write(lupri,*) 'get_available_memory: System identified to be MAC!'
         write(lupri,'(1X,a,g16.5)') 'Available Memory (GB) = ', MemoryAvailable
      end if
   endif

   ! Still no result --> Is this a LINUX system without a  readable /proc/memninfo?
   if(.not.memfound) then
      call get_available_memory_specific('LINUX',MemoryAvailable,memfound)
      if(memfound.and..not.sp) then
         write(lupri,*) 'get_available_memory: System identified to be LINUX!(using the top command)'
         write(lupri,'(1X,a,g16.5)') 'Available Memory (GB) = ', MemoryAvailable
      endif
   endif


   if(.not.memfound)then
      write(lupri,*) '******************** WARNING WARNING WARNING ***********************'
      write(lupri,*) 'get_available_memory: System type NOT identified!'
      write(lupri,'(1X,a,g16.5)') 'Default setting: Available Memory (GB) = ', MemoryAvailable
   endif


end subroutine get_available_memory


!> \brief Get how much memory is currently available using the systemtype
!> defined in the input.
!> \author Kasper Kristensen, modified by Patrick Ettenhuber
!> \date January 2012
subroutine get_available_memory_specific(systemtype,MemoryAvailable,memfound)
   !> Which system - currently only "MAC  " and "LINUX" are allowed, and in principle all systems that have the /proc/meminfo file with "MENFO"
   character(len=5),intent(in) :: systemtype
   !> Available memory measured in GB (using 1GB = 1000000000 bytes)
   real(realk),intent(inout) :: MemoryAvailable
   !> Was the memory information found?
   logical,intent(inout) :: MemFound
   character(len=90) :: command,string
   integer :: funit,i,length, endpos, startpos, ios
   integer(kind=8) :: memInt,intZero
   logical :: file_exist,doublecheck
   real(realk) :: ConvertToGB
   string=''
   ios = -1
   funit = 0

   ! Default values
   doublecheck=.false.
   memfound=.false.
   MemoryAvailable=1.0E0_realk
   ConvertToGB=0.0_realk   ! this is just set to catch possible errors
   funit=200

   if(systemtype=="MENFO")then

      open(unit=funit,file='/proc/meminfo',status='old',form='formatted',iostat=ios)
      if(ios/=0) then
         ! Something went wrong with opening file and we use the default values
         MemoryAvailable=1.0E0_realk
         memfound=.false.
         close(funit)
         return
      end if

      ! Read from file expected to contain memory info
      read(funit,'(a90)',iostat=ios) string
      read(funit,'(a90)',iostat=ios) string

      if(ios/=0) then
         ! Something went wrong with reading from file and we use the default values
         MemoryAvailable=1.0E0_realk
         memfound=.false.
         close(funit)
         return
      end if

      if(string(1:8)=="MemFree:")then
         MemFound=.true.
         length = LEN(string)
         do i=length,1,-1
            if(string(i-1:i)=='MB') then
               ConvertToGB = 1.0E-3_realk
               doublecheck=.true.
               endpos=i-2
               exit
               elseif(string(i-1:i)=='kB')then
               ConvertToGB = 1.0E-6_realk
               doublecheck=.true.
               endpos=i-2
               exit
            endif
         enddo

         if(.not. doublecheck)then
            print *,"not the correct format of /proc/meminfo --> check and update get_available_memory_specific"
            print *,"in lsutil/memory.f90"
            intZero = 0
            call memory_error_quit('not the correct format of /proc/meminfo',intZero)
         endif

         do i=endpos-1,1,-1
            if( (string(i:i)==' ') .or. (string(i:i)==',') ) then
               startpos=i+1
               exit
            end if
         end do

         ! Set available memory measured in GB
         read(string(startpos:endpos),*) memInt  ! memory in kB or MB as an integer
         MemoryAvailable = ConvertToGB*real(memInt) ! memory in GB as a real
      endif

   else
      ! Create file 'memoryfile' containing the relevant memory info
      call get_memoryfile(systemtype)

      ! Check if the file was constructed at all
      inquire(file='memoryfile',exist=file_exist)

      TheFileExist: if(file_exist) then

         ! Open file expected to contain memory info
         open(unit=funit,file='memoryfile',status='old',form='formatted',iostat=ios)
         if(ios/=0) then
            ! Something went wrong with opening file and we use the default values
            MemoryAvailable=1.0E0_realk
            memfound=.false.
            close(funit)
            return
         end if

         ! Read from file expected to contain memory info
         read(funit,'(a90)',iostat=ios) string
         if(ios/=0) then
            ! Something went wrong with reading from file and we use the default values
            MemoryAvailable=1.0E0_realk
            memfound=.false.
            close(funit)
            return
         end if

         ! Check if memory information is found in string
         call is_memory_info_found(string,memfound)

      end if TheFileExist


      MemInfoPresentInFile: if(memfound) then

         !> Using MAC the file looks like:
         !> PhysMem:  445M wired,  533M active,  226M inactive, 1471M used, 2615M free.
         !> Using Linux it looks like:
         !> Mem:  16505788k total, 16139716k used,   366072k free,     3348k buffers

         ! String element defining last number in free memory
         ! For the MAC example above this would be the position of "5" in 2615
         ! For the Linux example this would be the position of "2" in 366072
         length = LEN(string)
         do i=length,1,-1
            if( (string(i-5:i)=='M free') .or. (string(i-5:i)=='k free') .or.  (string(i-5:i)=='M unus') ) then

               ! This is a double-check that the memory file has the correct format
               doublecheck=.true.

               ! position of last digit in free memory
               endpos=i-6

               ! Conversion factor (kB to GB or MB to GB?)
               if(string(i-5:i-5)=='M') ConvertToGB = 1.0E-3_realk
               if(string(i-5:i-5)=='k') ConvertToGB = 1.0E-6_realk
               exit

            end if
         end do

         ! String element defining first number in free memory
         ! For the MAC example above this would be the position of "2" in 2615
         ! For the Linux example this would be the position of "3" in 366072
         do i=endpos,1,-1
            if( (string(i:i)==' ') .or. (string(i:i)==',') ) then
               startpos=i+1
               exit
            end if
         end do

         ! Set available memory measured in GB
         read(string(startpos:endpos),*) memInt  ! memory in kB or MB as an integer
         MemoryAvailable = ConvertToGB*real(memInt) ! memory in GB as a real

      end if MemInfoPresentInFile


      ! Safety double check
      if(.not. doublecheck ) then
         ! Something went wrong and we use the default values instead of trusting what we read
         MemoryAvailable=1.0E0_realk
         memfound=.false.
      end if

   endif
   close(funit)

end subroutine get_available_memory_specific


!> \brief Create file 'memoryfile' containing memory info. The file will be a single line.
!> Using MAC it looks like:
!> PhysMem:  445M wired,  533M active,  226M inactive, 1471M used, 2615M free.
!> Using Linux it looks like:
!> Mem:  16505788k total, 16139716k used,   366072k free,     3348k buffers
!> \author Kasper Kristensen
!> \date January 2012
subroutine get_memoryfile(systemtype)
   !> Which system - currently only "MAC  " and "LINUX" are allowed
   character(len=5),intent(in) :: systemtype
   character(len=90) :: command

   ! Call the system command to get a file containing the free memory info
   if(systemtype=='MAC  ') then
      write(command,'(a)') 'top -l 1 | grep "Mem:" > memoryfile'
      elseif(systemtype=='LINUX') then
      write(command,'(a)') 'top -b -n1 | grep "Mem:" > memoryfile'
   else
      stop 'get_memory_command: &
         & Currently only implemented for MAC and LINUX'
   end if
   call system(command)

end subroutine get_memoryfile


!> \brief Find out whether memory information string determined in get_available_memory
!> contains the seeked memory information, i.e. whether it contains
!> a string bit called "Mem:".
!> find out which system (MAC,Linux, or default) is being used.
!> \author Kasper Kristensen
!> \date January 2012
subroutine is_memory_info_found(string,memfound)
   !> String which may contain memory information
   character(len=90),intent(in) :: string
   !> Was the memory information found correctly?
   logical,intent(inout) :: memfound
   integer :: i,length

   length = len(string)

   memfound=.false.
   do i=1,length-3
      if(string(i:(i+3))=="Mem:") then
         memfound=.true.
         exit
      end if
   end do

end subroutine is_memory_info_found

! NOTE: If you add stuff here, remember to change
! longintbuffersize accordingly!
subroutine copy_from_mem_stats(longintbufferInt)
   implicit none
   integer(kind=long) :: longintbufferInt(longintbuffersize)
   longintbufferInt(1) = mem_allocated_global
   longintbufferInt(2) = mem_allocated_type_matrix
   longintbufferInt(3) = mem_allocated_real
   longintbufferInt(4) = mem_allocated_integer
   longintbufferInt(5) = mem_allocated_logical
   longintbufferInt(6) = mem_allocated_linkshell
   longintbufferInt(7) = mem_allocated_integrand
   longintbufferInt(8) = mem_allocated_integralitem
   longintbufferInt(9) = mem_allocated_overlap
   longintbufferInt(12) = mem_allocated_ODitem
   longintbufferInt(13) = mem_allocated_lstensor
   longintbufferInt(14) = mem_allocated_FMM
   longintbufferInt(15) = mem_allocated_AOBATCH
   longintbufferInt(16) = mem_allocated_ODBATCH
   longintbufferInt(17) = mem_allocated_LSAOTENSOR
   longintbufferInt(18) = mem_allocated_SLSAOTENSOR
   longintbufferInt(19) = mem_allocated_GLOBALLSAOTENSOR
   longintbufferInt(20) = mem_allocated_ATOMTYPEITEM
   longintbufferInt(21) = mem_allocated_ATOMITEM
   longintbufferInt(22) = mem_allocated_LSMATRIX
   longintbufferInt(23) = max_mem_used_global
   longintbufferInt(24) = max_mem_used_type_matrix
   longintbufferInt(25) = max_mem_used_real
   longintbufferInt(26) = max_mem_used_integer
   longintbufferInt(27) = max_mem_used_logical
   longintbufferInt(28) = max_mem_used_linkshell
   longintbufferInt(29) = max_mem_used_integrand
   longintbufferInt(30) = max_mem_used_integralitem
   longintbufferInt(31) = max_mem_used_overlap
   longintbufferInt(34) = max_mem_used_ODitem
   longintbufferInt(35) = max_mem_used_lstensor
   longintbufferInt(36) = max_mem_used_FMM
   longintbufferInt(37) = max_mem_used_AOBATCH
   longintbufferInt(38) = max_mem_used_ODBATCH
   longintbufferInt(39) = max_mem_used_LSAOTENSOR
   longintbufferInt(40) = max_mem_used_SLSAOTENSOR
   longintbufferInt(41) = max_mem_used_GLOBALLSAOTENSOR
   longintbufferInt(42) = max_mem_used_ATOMTYPEITEM
   longintbufferInt(43) = max_mem_used_ATOMITEM
   longintbufferInt(44) = max_mem_used_LSMATRIX
   longintbufferInt(45) = mem_allocated_DECORBITAL
   longintbufferInt(46) = max_mem_used_DECORBITAL
   longintbufferInt(47) = mem_allocated_DECFRAG
   longintbufferInt(48) = mem_allocated_BATCHTOORB
   longintbufferInt(49) = mem_allocated_MYPOINTER
   longintbufferInt(50) = mem_allocated_ARRAY2
   longintbufferInt(51) = mem_allocated_ARRAY4
   longintbufferInt(52) = mem_allocated_MP2DENS
   longintbufferInt(53) = mem_allocated_TRACEBACK
   longintbufferInt(54) = mem_allocated_MP2GRAD
   longintbufferInt(55) = max_mem_used_DECFRAG
   longintbufferInt(56) = max_mem_used_BATCHTOORB
   longintbufferInt(57) = max_mem_used_MYPOINTER
   longintbufferInt(58) = max_mem_used_ARRAY2
   longintbufferInt(59) = max_mem_used_ARRAY4
   longintbufferInt(60) = max_mem_used_MP2DENS
   longintbufferInt(61) = max_mem_used_TRACEBACK
   longintbufferInt(62) = max_mem_used_MP2GRAD
   longintbufferInt(63) = max_mem_used_overlapT
   longintbufferInt(64) = mem_allocated_overlapT
   longintbufferInt(65) = mem_allocated_complex
   longintbufferInt(66) = max_mem_used_complex
   longintbufferInt(67) = max_mem_used_character
   longintbufferInt(68) = mem_allocated_character
   longintbufferInt(69) = mem_allocated_ARRAY
   longintbufferInt(70) = max_mem_used_ARRAY
   longintbufferInt(71) = mem_allocated_mpi
   longintbufferInt(72) = max_mem_used_mpi
   longintbufferInt(73) = mem_allocated_DECAOBATCHINFO
   longintbufferInt(74) = max_mem_used_DECAOBATCHINFO
   longintbufferInt(75) = mem_allocated_PNOSpaceInfo
   longintbufferInt(76) = max_mem_used_PNOSpaceInfo
   longintbufferInt(77) = max_mem_used_type_matrix_MPIFULL
#ifdef MOD_UNRELEASED
   longintbufferInt(78) = mem_allocated_lvec_data
   longintbufferInt(79) = mem_allocated_lattice_cell
#endif
   longintbufferInt(80) = mem_allocated_fragmentAOS
   longintbufferInt(81) = max_mem_used_fragmentAOS
   ! NOTE: If you add stuff here, remember to change
   ! longintbuffersize accordingly!
end subroutine copy_from_mem_stats

! NOTE: If you add stuff here, remember to change
! longintbuffersize accordingly!
subroutine copy_to_mem_stats(longintbufferInt)
   implicit none
   integer(kind=long) :: longintbufferInt(longintbuffersize)
   mem_allocated_global = longintbufferInt(1)
   mem_allocated_type_matrix = longintbufferInt(2)
   mem_allocated_real = longintbufferInt(3)
   mem_allocated_integer = longintbufferInt(4)
   mem_allocated_logical = longintbufferInt(5)
   mem_allocated_linkshell = longintbufferInt(6)
   mem_allocated_integrand = longintbufferInt(7)
   mem_allocated_integralitem = longintbufferInt(8)
   mem_allocated_overlap = longintbufferInt(9)
   mem_allocated_ODitem = longintbufferInt(12)
   mem_allocated_lstensor = longintbufferInt(13)
   mem_allocated_FMM = longintbufferInt(14)
   mem_allocated_AOBATCH = longintbufferInt(15)
   mem_allocated_ODBATCH = longintbufferInt(16)
   mem_allocated_LSAOTENSOR = longintbufferInt(17)
   mem_allocated_SLSAOTENSOR = longintbufferInt(18)
   mem_allocated_GLOBALLSAOTENSOR = longintbufferInt(19)
   mem_allocated_ATOMTYPEITEM = longintbufferInt(20)
   mem_allocated_ATOMITEM = longintbufferInt(21)
   mem_allocated_LSMATRIX = longintbufferInt(22)
   max_mem_used_global = longintbufferInt(23)
   max_mem_used_type_matrix = longintbufferInt(24)
   max_mem_used_real = longintbufferInt(25)
   max_mem_used_integer = longintbufferInt(26)
   max_mem_used_logical = longintbufferInt(27)
   max_mem_used_linkshell = longintbufferInt(28)
   max_mem_used_integrand = longintbufferInt(29)
   max_mem_used_integralitem = longintbufferInt(30)
   max_mem_used_overlap = longintbufferInt(31)
   max_mem_used_ODitem = longintbufferInt(34)
   max_mem_used_lstensor = longintbufferInt(35)
   max_mem_used_FMM = longintbufferInt(36)
   max_mem_used_AOBATCH = longintbufferInt(37)
   max_mem_used_ODBATCH = longintbufferInt(38)
   max_mem_used_LSAOTENSOR = longintbufferInt(39)
   max_mem_used_SLSAOTENSOR = longintbufferInt(40)
   max_mem_used_GLOBALLSAOTENSOR = longintbufferInt(41)
   max_mem_used_ATOMTYPEITEM = longintbufferInt(42)
   max_mem_used_ATOMITEM = longintbufferInt(43)
   max_mem_used_LSMATRIX = longintbufferInt(44)
   mem_allocated_DECORBITAL = longintbufferInt(45)
   max_mem_used_DECORBITAL = longintbufferInt(46)
   mem_allocated_DECFRAG = longintbufferInt(47)
   mem_allocated_BATCHTOORB = longintbufferInt(48)
   mem_allocated_MYPOINTER = longintbufferInt(49)
   mem_allocated_ARRAY2 = longintbufferInt(50)
   mem_allocated_ARRAY4 = longintbufferInt(51)
   mem_allocated_MP2DENS = longintbufferInt(52)
   mem_allocated_TRACEBACK = longintbufferInt(53)
   mem_allocated_MP2GRAD = longintbufferInt(54)
   max_mem_used_DECFRAG = longintbufferInt(55)
   max_mem_used_BATCHTOORB = longintbufferInt(56)
   max_mem_used_MYPOINTER = longintbufferInt(57)
   max_mem_used_ARRAY2 = longintbufferInt(58)
   max_mem_used_ARRAY4 = longintbufferInt(59)
   max_mem_used_MP2DENS = longintbufferInt(60)
   max_mem_used_TRACEBACK = longintbufferInt(61)
   max_mem_used_MP2GRAD = longintbufferInt(62)
   max_mem_used_overlapT = longintbufferInt(63)
   mem_allocated_overlapT = longintbufferInt(64)
   mem_allocated_complex = longintbufferInt(65)
   max_mem_used_complex = longintbufferInt(66)
   max_mem_used_character = longintbufferInt(67)
   mem_allocated_character = longintbufferInt(68)
   mem_allocated_ARRAY = longintbufferInt(69)
   max_mem_used_ARRAY = longintbufferInt(70)
   mem_allocated_mpi = longintbufferInt(71)
   max_mem_used_mpi = longintbufferInt(72)
   mem_allocated_DECAOBATCHINFO = longintbufferInt(73)
   max_mem_used_DECAOBATCHINFO = longintbufferInt(74)
   mem_allocated_PNOSpaceInfo = longintbufferInt(75)
   max_mem_used_PNOSpaceInfo = longintbufferInt(76)
   max_mem_used_type_matrix_MPIFULL = longintbufferInt(77)
#ifdef MOD_UNRELEASED
   mem_allocated_lvec_data = longintbufferInt(78)
   mem_allocated_lattice_cell = longintbufferInt(79)
#endif
   mem_allocated_fragmentAOS = longintbufferInt(80)
   max_mem_used_fragmentAOS = longintbufferInt(81)
   ! NOTE: If you add stuff here, remember to change
   ! longintbuffersize accordingly!
end subroutine copy_to_mem_stats

END MODULE memory_handling

