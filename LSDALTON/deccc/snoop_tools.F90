!> @file
!> Utils for DEC subroutines
!> \author Marcin Ziolkowski (modified by Kasper Kristensen)
module snoop_tools_module

  use fundamental
  use precision
  use lstiming
  use ls_util!,only: dgemm_ts
  use typedeftype!,only: lsitem
  use molecule_module!, only: get_geometry
  use files!,only:lsopen,lsclose
  use DALTONINFO!, only: ls_free
  use dec_typedef_module
  use memory_handling!, only: mem_alloc, mem_dealloc, mem_allocated_global,&
  !       & stats_mem, get_avaiLable_memory
  use,intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
  use matrix_module!, only:matrix
  use matrix_operations
  use matrix_util
  use matrix_operations_aux
  use tensor_interface_module
  use array4_simple_operations
  use IntegralInterfaceMOD!, only: ii_get_h1, ii_get_nucpot
  use BUILDAOBATCH
  use DALTONINFO

#ifdef VAR_MPI
  use infpar_module
  use lsmpi_op
#endif

#ifdef VAR_PAPI
  use papi_module
#endif


contains

  subroutine testme()
    print *, 'Testing'

  end subroutine testme


end module snoop_tools_module
