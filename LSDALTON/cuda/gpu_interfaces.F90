!> @file
!> GPU interfaces
!> \brief: gpu interface module
!> \author: Janus Juul Eriksen
!> \date: 2015, Aarhus
module gpu_interfaces

  use iso_c_binding
  use precision

  !> module variable to count the FLOPs done on the GPU
  real(realk), save :: FLOPonGPU

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
       implicit none
       integer (C_SIZE_T) :: total,free
    end subroutine get_dev_mem

  end interface

#endif

contains

  subroutine ls_dgemm_acc(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc,na,nb,nc,acc_handle,cublas_handle)

       use precision
       use iso_c_binding
#ifdef VAR_OPENACC
       use openacc
#endif

       character(len=1) :: transa,transb
       integer :: m,n,k,lda,ldb,ldc
       integer(kind=8) :: na,nb,nc
       real(realk), dimension(na)  :: a
       real(realk), dimension(nb)  :: b
       real(realk), dimension(nc)  :: c
       real(realk)  :: alpha,beta
#ifdef VAR_OPENACC
       integer(kind=acc_handle_kind) :: acc_handle
#else
       integer :: acc_handle
#endif
       type(c_ptr) :: cublas_handle ! NOTE: To use cublas, make sure you've called cublasCreate_v2 beforehand!!!
       logical :: async,false_arg1,false_arg2
       integer :: transa_2,transb_2

#ifdef VAR_OPENACC

       async = .false.

       if (acc_handle .ne. acc_async_sync) async = .true.

       false_arg1 = .true.; false_arg2 = .true.

       if ((transa .eq. 'n') .or. (transa .eq. 'N') .or. (transa .eq. 't') .or. (transa .eq. 'T')) false_arg1 = .false.
       if ((transb .eq. 'n') .or. (transb .eq. 'N') .or. (transb .eq. 't') .or. (transb .eq. 'T')) false_arg2 = .false.

       if (false_arg1) call lsquit('wrong argument to transa in ls_dgemm_acc',-1)
       if (false_arg2) call lsquit('wrong argument to transb in ls_dgemm_acc',-1)

#ifdef VAR_CUBLAS

       transa_2 = 0; transb_2 = 0

       if ((transa .eq. 't') .or. (transa .eq. 'T')) transa_2 = 1
       if ((transb .eq. 't') .or. (transb .eq. 'T')) transb_2 = 1

#endif

#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)

       if (async) then

!$acc host_data use_device(a,b,c)
          call dgemm_acc_openacc_async(acc_handle,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$acc end host_data

       else

!$acc host_data use_device(a,b,c)
          call dgemm_acc(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
!$acc end host_data

       endif

#elif defined(VAR_CUBLAS)

!$acc host_data use_device(a,b,c)
       stat = cublasDgemm_v2(cublas_handle,int(transa_2,kind=4),int(transb_2,kind=4),int(m,kind=4),int(n,kind=4),int(k,kind=4),&
                             & alpha,c_loc(a),int(lda,kind=4),c_loc(b),int(ldb,kind=4),&
                             & beta,c_loc(c),int(ldc,kind=4))
!$acc end host_data

#endif

       ! calculate the gpu flop count
       call addDGEMM_FLOPonGPUaccouting(m,n,k,beta)

#else

       ! call ordinary cpu dgemm
       call dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

#endif

  end subroutine ls_dgemm_acc

  subroutine AddFLOP_FLOPonGPUaccouting(inputFLOPonGPU)

    implicit none

    real(realk), intent(in) :: inputFLOPonGPU

    FLOPonGPU = FLOPonGPU + inputFLOPonGPU

  end subroutine AddFLOP_FLOPonGPUaccouting

  subroutine addDGEMM_FLOPonGPUaccouting(M,N,K,beta)

    implicit none

    integer, intent(in) :: M,N,K
    real(realk), intent(in) :: beta

    IF(ABS(beta).GT.1.0E-10_realk)THEN
       FLOPonGPU = FLOPonGPU + 2*M*N*K
    ELSE
       FLOPonGPU = FLOPonGPU + M*N*(2*K-1)
    ENDIF

  end subroutine AddDGEMM_FLOPonGPUaccouting

  subroutine init_FLOPonGPUaccouting()

    implicit none

    FLOPonGPU = 0.0E0_realk

  end subroutine Init_FLOPonGPUaccouting

  subroutine extract_FLOPonGPUaccouting(outFLOPonGPU)

    implicit none
    real(realk), intent(inout) :: outFLOPonGPU

    outFLOPonGPU = FLOPonGPU

  end subroutine Extract_FLOPonGPUaccouting
 
end module gpu_interfaces
