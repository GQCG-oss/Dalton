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

    ! cublasSgemm_v2
    integer (C_INT) function cublasSgemm_v2 ( handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc ) bind (C, name="cublasSgemm_v2")
      use iso_c_binding
      implicit none
      type (C_PTR), value :: handle
      type (C_PTR), value :: A, B, C
      integer (C_INT), value :: m, n, k, lda, ldb, ldc
      integer (C_INT), value :: transa, transb
      real (C_FLOAT) :: alpha, beta
    end function cublasSgemm_v2

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

       character(len=1), intent(in) :: transa,transb
       integer, intent(in) :: m,n,k,lda,ldb,ldc
       integer(kind=8), intent(in) :: na,nb,nc
       real(realk), dimension(na), intent(in), target :: a
       real(realk), dimension(nb), intent(in), target :: b
       real(realk), dimension(nc), intent(inout), target :: c
       real(realk), intent(in)  :: alpha,beta
#ifdef VAR_OPENACC
       integer(kind=acc_handle_kind), intent(in) :: acc_handle
#else
       integer, intent(in) :: acc_handle
#endif
       type(c_ptr), intent(in) :: cublas_handle ! NOTE: To use cublas, make sure you've called cublasCreate_v2 beforehand!!!
       logical :: async,false_arg1,false_arg2
       integer :: transa_2,transb_2

!       print *,'transa,transb     = ',transa,transb
!       print *,'m,n,k,lda,ldb,ldc = ',m,n,k,lda,ldb,ldc
!       print *,'alpha,beta        = ',alpha,beta
!       print *,'acc_handle        = ',acc_handle
!       print *,'cublas_handle     = ',cublas_handle

       false_arg1 = .true.; false_arg2 = .true.
       if ((transa .eq. 'n') .or. (transa .eq. 'N') .or. (transa .eq. 't') .or. (transa .eq. 'T')) false_arg1 = .false.
       if ((transb .eq. 'n') .or. (transb .eq. 'N') .or. (transb .eq. 't') .or. (transb .eq. 'T')) false_arg2 = .false.

       if (false_arg1) call lsquit('wrong argument to transa in ls_dgemm_acc',-1)
       if (false_arg2) call lsquit('wrong argument to transb in ls_dgemm_acc',-1)

#ifdef VAR_OPENACC

       async = .false.
       if (acc_handle .ne. acc_async_sync) async = .true.

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
       call addGEMM_FLOPonGPUaccouting(int(m,kind=8),int(n,kind=8),int(k,kind=8),beta)

#else

       ! call ordinary cpu dgemm
       call dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

#endif

  end subroutine ls_dgemm_acc


  subroutine ls_sgemm_acc(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc,na,nb,nc,acc_handle,cublas_handle)

       use precision
       use iso_c_binding
#ifdef VAR_OPENACC
       use openacc
#endif

       character(len=1), intent(in) :: transa,transb
       integer, intent(in) :: m,n,k,lda,ldb,ldc
       integer(kind=8), intent(in) :: na,nb,nc
       real(realk), dimension(na), intent(in), target :: a
       real(realk), dimension(nb), intent(in), target :: b
       real(kind=4), dimension(nc), intent(inout), target :: c
       real(realk), intent(in)  :: alpha,beta
#ifdef VAR_OPENACC
       integer(kind=acc_handle_kind), intent(in) :: acc_handle
#else
       integer, intent(in) :: acc_handle
#endif
       type(c_ptr), intent(in) :: cublas_handle ! NOTE: To use cublas, make sure you've called cublasCreate_v2 beforehand!!!
       logical :: async,false_arg1,false_arg2
       integer :: transa_2,transb_2

!       print *,'transa,transb     = ',transa,transb
!       print *,'m,n,k,lda,ldb,ldc = ',m,n,k,lda,ldb,ldc
!       print *,'alpha,beta        = ',alpha,beta
!       print *,'acc_handle        = ',acc_handle
!       print *,'cublas_handle     = ',cublas_handle

       false_arg1 = .true.; false_arg2 = .true.
       if ((transa .eq. 'n') .or. (transa .eq. 'N') .or. (transa .eq. 't') .or. (transa .eq. 'T')) false_arg1 = .false.
       if ((transb .eq. 'n') .or. (transb .eq. 'N') .or. (transb .eq. 't') .or. (transb .eq. 'T')) false_arg2 = .false.

       if (false_arg1) call lsquit('wrong argument to transa in ls_sgemm_acc',-1)
       if (false_arg2) call lsquit('wrong argument to transb in ls_sgemm_acc',-1)

#ifdef VAR_OPENACC

       async = .false.
       if (acc_handle .ne. acc_async_sync) async = .true.

#ifdef VAR_CUBLAS

       transa_2 = 0; transb_2 = 0
       if ((transa .eq. 't') .or. (transa .eq. 'T')) transa_2 = 1
       if ((transb .eq. 't') .or. (transb .eq. 'T')) transb_2 = 1

#endif

#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)

       if (async) then

!$acc host_data use_device(a,b,c)
          call sgemm_acc_openacc_async(acc_handle,transa,transb,m,n,k,&
                   & real(alpha,kind=4),real(a,kind=4),lda,real(b,kind=4),ldb,real(beta,kind=4),real(c,kind=4),ldc)
!$acc end host_data

       else

!$acc host_data use_device(a,b,c)
          call sgemm_acc(transa,transb,m,n,k,&
                   & real(alpha,kind=4),real(a,kind=4),lda,real(b,kind=4),ldb,real(beta,kind=4),real(c,kind=4),ldc)
!$acc end host_data

       endif

#elif defined(VAR_CUBLAS)

!!$acc host_data use_device(a,b,c)
!       stat = cublasSgemm_v2(cublas_handle,int(transa_2,kind=4),int(transb_2,kind=4),int(m,kind=4),int(n,kind=4),int(k,kind=4),&
!                             & real(alpha,kind=4),c_loc(real(a,kind=4)),int(lda,kind=4),c_loc(real(b,kind=4)),int(ldb,kind=4),&
!                             & real(beta,kind=4),c_loc(real(c,kind=4)),int(ldc,kind=4))
!!$acc end host_data

#endif

       ! calculate the gpu flop count
       call addGEMM_FLOPonGPUaccouting(int(m,kind=8),int(n,kind=8),int(k,kind=8),beta)

#else

       ! call ordinary cpu sgemm
       call sgemm(transa,transb,m,n,k,&
                & real(alpha,kind=4),real(a,kind=4),lda,real(b,kind=4),ldb,real(beta,kind=4),real(c,kind=4),ldc)

#endif

  end subroutine ls_sgemm_acc


  subroutine ls_ddot_acc(n,a,inca,b,incb,alpha,c,acc_handle,cublas_handle)

       use precision
       use iso_c_binding
#ifdef VAR_OPENACC
       use openacc
#endif

       integer(kind=8), intent(in) :: n
       integer, intent(in) :: inca,incb
       real(realk), dimension(n), intent(in), target :: a,b
       real(realk), intent(inout), target :: c
       real(realk), intent(in)  :: alpha
#ifdef VAR_OPENACC
       integer(kind=acc_handle_kind), intent(in) :: acc_handle
#else
       integer, intent(in) :: acc_handle
#endif
       type(c_ptr), intent(in) :: cublas_handle ! NOTE: To use cublas, make sure you've called cublasCreate_v2 beforehand!!!
       logical :: async
       !> ddot
       real(realk), external :: ddot

!       print *,'n,inca,incb       = ',n,inca,incb
!       print *,'alpha             = ',alpha,beta
!       print *,'acc_handle        = ',acc_handle
!       print *,'cublas_handle     = ',cublas_handle

#ifdef VAR_OPENACC

       async = .false.
       if (acc_handle .ne. acc_async_sync) async = .true.

#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)

       if (async) then

!$acc host_data use_device(a,b,c)
          call dgemm_acc_openacc_async(acc_handle,'n','n',1,1,n,alpha,a,1,b,n,1.0E0_realk,c,1)
!$acc end host_data

       else

!$acc host_data use_device(a,b,c)
!          call dgemm_acc('n','n',1,1,n,alpha,a,1,b,n,1.0E0_realk,c,1)
          c = c + alpha * ddot_acc(n,a,1,b,1)
!$acc end host_data

       endif

#elif defined(VAR_CUBLAS)

!$acc host_data use_device(a,b,c)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(n,kind=4),&
                             & alpha,c_loc(a),int(1,kind=4),c_loc(b),int(n,kind=4),&
                             & 1.0E0_realk,c_loc(c),int(1,kind=4))
!$acc end host_data

#endif

       ! calculate the gpu flop count
       call addDOT_FLOPonGPUaccouting(n)

#else

       ! call ordinary cpu ddot
       c = c + alpha * ddot(n,a,1,b,1)

#endif

  end subroutine ls_ddot_acc


  subroutine AddFLOP_FLOPonGPUaccouting(inputFLOPonGPU)

    implicit none

    real(realk), intent(in) :: inputFLOPonGPU

    FLOPonGPU = FLOPonGPU + inputFLOPonGPU

  end subroutine AddFLOP_FLOPonGPUaccouting

  subroutine addGEMM_FLOPonGPUaccouting(M,N,K,beta)

    implicit none

    integer(kind=8), intent(in) :: M,N,K
    real(realk), intent(in) :: beta

    IF(ABS(beta).GT.1.0E-10_realk)THEN
       FLOPonGPU = FLOPonGPU + 2*M*N*K
    ELSE
       FLOPonGPU = FLOPonGPU + M*N*(2*K-1)
    ENDIF

  end subroutine AddGEMM_FLOPonGPUaccouting

  subroutine addDOT_FLOPonGPUaccouting(K)

    implicit none

    integer(kind=8), intent(in) :: K

    ! (2*K - 1) operations [K multiplications + (K - 1) additions] + 1 final multiplication + 1 final addition
    FLOPonGPU = FLOPonGPU + (2*K + 1)

  end subroutine AddDOT_FLOPonGPUaccouting

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
