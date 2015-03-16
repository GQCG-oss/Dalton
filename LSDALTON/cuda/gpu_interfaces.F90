!> @file
!> GPU interfaces
!> \brief: gpu interface module
!> \author: Janus Juul Eriksen
!> \date: 2015, Aarhus
module gpu_interfaces

#ifdef VAR_OPENACC
#ifdef VAR_CUBLAS

  interface

     ! cublasCreate
     integer (C_INT) function cublasCreate_v2 ( handle ) bind (C, name="cublasCreate_v2")
       use iso_c_binding
       implicit none
       type (C_PTR) :: handle
     end function cublasCreate_v2

     ! cublasDestroy
     integer (C_INT) function cublasDestroy_v2 ( handle ) bind (C, name="cublasDestroy_v2")
       use iso_c_binding
       implicit none
       type (C_PTR), value :: handle
     end function cublasDestroy_v2

    ! cublasDgemm_v2
    integer (C_INT) function cublasDgemm_v2 ( handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc ) bind (C, name="cublasDgemm_v2")
      use iso_c_binding
      implicit none
      type (C_PTR), value :: handle
      type (C_PTR), value :: A, B, C
      integer (C_INT), value :: m, n, k, lda, ldb, ldc
      integer (C_INT), value :: transa, transb
      real (C_DOUBLE) :: alpha, beta
    end function cublasDgemm_v2

    ! cublasSetStream_v2
    integer (C_INT) function cublasSetStream_v2 ( handle, stream ) bind (C, name="cublasSetStream_v2")
      use iso_c_binding
      implicit none
      type (C_PTR), value :: handle
      type (C_PTR), value :: stream
    end function cublasSetStream_v2

  end interface

#endif
#endif

#ifdef VAR_CUDA

  interface

    subroutine get_dev_mem( total , free ) bind(C, name="get_dev_mem")

       use iso_c_binding

       integer (C_SIZE_T) :: total,free

    end subroutine get_dev_mem

  end interface

#endif

contains

  subroutine gpu_interfaces_dummy_routine()

  end subroutine gpu_interfaces_dummy_routine

end module gpu_interfaces
