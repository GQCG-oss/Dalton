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

#ifdef VAR_OPENACC
  subroutine ls_dgemm_acc(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc,acc_handle,cublas_handle)

       use openacc
       use iso_c_binding

       character, intent(inout) :: transa,transb
       integer, intent(inout) :: m,n,k,lda,ldb,ldc
       real(realk), dimension(lda,*), intent(inout) :: a
       real(realk), dimension(ldb,*), intent(inout) :: b
       real(realk), dimension(ldc,*), intent(inout) :: c
       real(realk), intent(inout) :: alpha,beta
       integer(kind=acc_handle_kind), intent(inout), optional :: acc_handle
       type(c_ptr), intent(inout), optional :: cublas_handle
       logical :: async
       integer(kind=acc_handle_kind) :: handle_2
       integer :: transa_2,transb_2
#ifdef VAR_PGF90
       integer*4, external :: acc_get_cuda_stream
#endif

       async = .false.

       if (present(acc_handle)) async = .true.

       if (present(cublas_handle)) then

          handle_2 = acc_get_cuda_stream(cublas_handle)

          if (handle_2 .ne. acc_async_sync) async = .true.

       endif

       if ((transa .ne. 'n') .or. (transa .ne. 'N') .or. (transa .ne. 't') .or. (transa .ne. 'T') then

          call lsquit('wrong argument to transa in ls_dgemm_acc',DECinfo%output)

       elseif ((transb .ne. 'n') .or. (transb .ne. 'N') .or. (transb .ne. 't') .or. (transb .ne. 'T') then

          call lsquit('wrong argument to transb in ls_dgemm_acc',DECinfo%output)

       endif

#ifdef defined(VAR_CUBLAS)

       if ((transa .eq. 'n') .or. (transa .eq. 'N')) then

          transa_2 = 0 

       else

          transa_2 = 1

       endif
       if ((transb .eq. 'n') .or. (transb .eq. 'N')) then

          transb_2 = 0

       else

          transb_2 = 1

       endif

#endif

#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)

       if (async) then

!$acc host_data use_device(a,b,c)
          call dgemm_acc_openacc_async(handle,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$acc end host_data

       else

!$acc host_data use_device(a,b,c)
          call dgemm_acc(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$acc end host_data

       endif

#elif defined(VAR_CUBLAS)

       if (.not. present(cublas_handle)) call lsquit('defined(VAR_CUBLAS), but no cublas handle.',DECinfo%output) 

!$acc host_data use_device(a,b,c)
       stat = cublasDgemm_v2(cublas_handle,int(transa_2,kind=4),int(transb_2,kind=4),int(m,kind=4),int(n,kind=4),int(k,kind=4),&
                             & alpha,c_loc(a),int(lda,kind=4),c_loc(b),int(ldb,kind=4),&
                             & beta,c_loc(c),int(ldc,kind=4))
!$acc end host_data

#endif

  end subroutine ls_dgemm_acc
#endif

  subroutine gpu_interfaces_dummy_routine()

  end subroutine gpu_interfaces_dummy_routine

end module gpu_interfaces
