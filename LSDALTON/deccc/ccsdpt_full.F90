!> @file
!> DEC-CCSD(T) kernels
!> \brief: ccsd(t) module
!> \author: Janus Juul Eriksen
!> \date: 2012-2014, Aarhus
module ccsdpt_full_module
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer, c_size_t

  use precision
  use Fundamental, only: bohr_to_angstrom
#ifdef VAR_OPENACC
  use openacc
#endif
  use gpu_interfaces

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use reorder_frontend_module
  use dec_workarounds_module
  
#ifdef MOD_UNRELEASED
  public :: trip_generator_ijk_case1
  public :: trip_generator_ijk_case2
  public :: trip_generator_ijk_case3
  public :: trip_generator_abc_case1
  public :: trip_generator_abc_case2
  public :: trip_generator_abc_case3
  public :: trip_denom_ijk
  public :: trip_denom_abc
  public :: ccsdpt_energy_full_ijk_case1
  public :: ccsdpt_energy_full_ijk_case2
  public :: ccsdpt_energy_full_ijk_case3
  public :: ccsdpt_energy_full_abc_case1
  public :: ccsdpt_energy_full_abc_case2
  public :: ccsdpt_energy_full_abc_case3
#endif

  private

contains

#ifdef MOD_UNRELEASED


  !> \brief: generator for triples amplitudes, case(1)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_ijk_case1(oindex1,oindex3,no,nv,ccsd_doubles_1,ccsd_doubles_3,&
                                & vvvo_tile_1,vvvo_tile_3,ovoo_tile_11,&
                                & ovoo_tile_13,ovoo_tile_31,trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex3,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(real_pt), dimension(nv,nv,no) :: ccsd_doubles_1, ccsd_doubles_3
    !> tiles of ovoo 2-el integrals
    real(real_pt), dimension(no,nv) :: ovoo_tile_11, ovoo_tile_13, ovoo_tile_31
    !> tiles of vvvo 2-el integrals
    real(real_pt), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_1, vvvo_tile_3
    !> triples amplitude and work array
    real(real_pt), dimension(nv,nv,nv) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs 
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! iik,iki
    call trip_amplitudes_ijk_virt(oindex1,oindex1,oindex3,no,nv,ccsd_doubles_1(:,:,oindex1),&
                            & vvvo_tile_3,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex3,oindex1,no,nv,ccsd_doubles_1,&
                            & ovoo_tile_13,trip_tmp,handle,cublas_handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_ampl,'=',trip_tmp,(i8*nv)*(i8*nv**2))
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_ampl = trip_tmp
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
#endif

    ! kii,iik
    call trip_amplitudes_ijk_virt(oindex3,oindex1,oindex1,no,nv,ccsd_doubles_3(:,:,oindex1),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex1,oindex3,no,nv,ccsd_doubles_1,&
                            & ovoo_tile_31,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

    ! iki.kii
    call trip_amplitudes_ijk_virt(oindex1,oindex3,oindex1,no,nv,ccsd_doubles_1(:,:,oindex3),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex3,oindex1,oindex1,no,nv,ccsd_doubles_3,&
                            & ovoo_tile_11,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_ijk_case1


  !> \brief: generator for triples amplitudes, case(1)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_abc_case1(vindex1,vindex3,no,nv,ccsd_doubles_1,ccsd_doubles_3,&
                                & ooov_tile_1,ooov_tile_3,vovv_11,vovv_13,vovv_31,&
                                & trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> a, b, a, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex3,no,nv
    !> ccsd doubles
    real(real_pt), dimension(no,no,nv), intent(inout)  :: ccsd_doubles_1,ccsd_doubles_3
    !> vovv integrals
    real(real_pt), dimension(nv,no), intent(inout)  :: vovv_11,vovv_13,vovv_31
    !> tiles of ooov 2-el integrals
    real(real_pt), dimension(no,no,no), intent(inout) :: ooov_tile_1, ooov_tile_3
    !> triples amplitude and work array
    real(real_pt), dimension(no,no,no) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! aac,aca
    call trip_amplitudes_abc_occ(vindex1,vindex1,vindex3,no,nv,ccsd_doubles_1(:,:,vindex1),&
                            & ooov_tile_3,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex1,vindex3,vindex1,no,nv,ccsd_doubles_1,&
                            & vovv_13,trip_tmp,handle,cublas_handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_ampl,'=',trip_tmp,(i8*no)*no**2)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_ampl = trip_tmp
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
#endif

    ! caa,aac
    call trip_amplitudes_abc_occ(vindex3,vindex1,vindex1,no,nv,ccsd_doubles_3(:,:,vindex1),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex1,vindex1,vindex3,no,nv,ccsd_doubles_1,&
                            & vovv_31,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

    ! aca.caa
    call trip_amplitudes_abc_occ(vindex1,vindex3,vindex1,no,nv,ccsd_doubles_1(:,:,vindex3),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex3,vindex1,vindex1,no,nv,ccsd_doubles_3,&
                            & vovv_11,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_abc_case1


  !> \brief: generator for triples amplitudes, case(2)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_ijk_case2(oindex1,oindex2,no,nv,ccsd_doubles_1,ccsd_doubles_2,&
                                & vvvo_tile_1,vvvo_tile_2,ovoo_tile_12,&
                                & ovoo_tile_21,ovoo_tile_22,trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(real_pt), dimension(nv,nv,no) :: ccsd_doubles_1, ccsd_doubles_2
    !> tiles of ovoo 2-el integrals
    real(real_pt), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_21, ovoo_tile_22
    !> tiles of vvvo 2-el integrals
    real(real_pt), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_1, vvvo_tile_2
    !> triples amplitude and work array
    real(real_pt), dimension(nv,nv,nv) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! ijj.jji
    call trip_amplitudes_ijk_virt(oindex1,oindex2,oindex2,no,nv,ccsd_doubles_1(:,:,oindex2),&
                            & vvvo_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex2,oindex2,oindex1,no,nv,ccsd_doubles_2,&
                            & ovoo_tile_12,trip_tmp,handle,cublas_handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_ampl,'=',trip_tmp,(i8*nv)*(i8*nv**2))
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_ampl = trip_tmp
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
#endif

    ! jij,ijj
    call trip_amplitudes_ijk_virt(oindex2,oindex1,oindex2,no,nv,ccsd_doubles_2(:,:,oindex1),&
                            & vvvo_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex2,oindex2,no,nv,ccsd_doubles_1,&
                            & ovoo_tile_22,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
#endif 

    ! jji,jij
    call trip_amplitudes_ijk_virt(oindex2,oindex2,oindex1,no,nv,ccsd_doubles_2(:,:,oindex2),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex2,oindex1,oindex2,no,nv,ccsd_doubles_2,&
                            & ovoo_tile_21,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_ijk_case2


  !> \brief: generator for triples amplitudes, case(2)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_abc_case2(vindex1,vindex2,no,nv,ccsd_doubles_1,ccsd_doubles_2,&
                                & ooov_tile_1,ooov_tile_2,vovv_12,vovv_21,vovv_22,&
                                & trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> a, b, a, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,no,nv
    real(real_pt), dimension(no,no,nv), intent(inout)  :: ccsd_doubles_1,ccsd_doubles_2
    !> vovv integrals
    real(real_pt), dimension(nv,no), intent(inout)  :: vovv_12,vovv_21,vovv_22
    !> tiles of ooov 2-el integrals
    real(real_pt), dimension(no,no,no), intent(inout) :: ooov_tile_1, ooov_tile_2
    !> triples amplitude and work array
    real(real_pt), dimension(no,no,no) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! abb.bba
    call trip_amplitudes_abc_occ(vindex1,vindex2,vindex2,no,nv,ccsd_doubles_1(:,:,vindex2),&
                            & ooov_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex2,vindex2,vindex1,no,nv,ccsd_doubles_2,&
                            & vovv_12,trip_tmp,handle,cublas_handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_ampl,'=',trip_tmp,(i8*no)*no**2)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_ampl = trip_tmp
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
#endif

    ! bab,abb
    call trip_amplitudes_abc_occ(vindex2,vindex1,vindex2,no,nv,ccsd_doubles_2(:,:,vindex1),&
                            & ooov_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex1,vindex2,vindex2,no,nv,ccsd_doubles_1,&
                            & vovv_22,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
#endif

    ! bba,bab
    call trip_amplitudes_abc_occ(vindex2,vindex2,vindex1,no,nv,ccsd_doubles_2(:,:,vindex2),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex2,vindex1,vindex2,no,nv,ccsd_doubles_2,&
                            & vovv_21,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_abc_case2


  !> \brief: generator for triples amplitudes, case(3)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_ijk_case3(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_1,ccsd_doubles_2,ccsd_doubles_3,&
                                & vvvo_tile_1,vvvo_tile_2,vvvo_tile_3,&
                                & ovoo_tile_12,ovoo_tile_13,ovoo_tile_21,ovoo_tile_23,ovoo_tile_31,ovoo_tile_32,&
                                & trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,oindex3,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(real_pt), dimension(nv,nv,no) :: ccsd_doubles_1, ccsd_doubles_2, ccsd_doubles_3
    !> tiles of ovoo 2-el integrals
    real(real_pt), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_13, ovoo_tile_21
    real(real_pt), dimension(no,nv) :: ovoo_tile_23, ovoo_tile_31, ovoo_tile_32
    !> tiles of vvvo 2-el integrals 
    real(real_pt), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_1, vvvo_tile_2, vvvo_tile_3
    !> triples amplitude and work array
    real(real_pt), dimension(nv,nv,nv) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! ijk.jki
    call trip_amplitudes_ijk_virt(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_1(:,:,oindex2),&
                            & vvvo_tile_3,trip_ampl,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex2,oindex3,oindex1,no,nv,ccsd_doubles_2,&
                            & ovoo_tile_13,trip_ampl,handle,cublas_handle)

    ! jik,ikj
    call trip_amplitudes_ijk_virt(oindex2,oindex1,oindex3,no,nv,ccsd_doubles_2(:,:,oindex1),&
                            & vvvo_tile_3,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex3,oindex2,no,nv,ccsd_doubles_1,&
                            & ovoo_tile_23,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
#endif

    ! kij,ijk
    call trip_amplitudes_ijk_virt(oindex3,oindex1,oindex2,no,nv,ccsd_doubles_3(:,:,oindex1),&
                            & vvvo_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_1,&
                            & ovoo_tile_32,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
#endif

    ! jki,kij
    call trip_amplitudes_ijk_virt(oindex2,oindex3,oindex1,no,nv,ccsd_doubles_2(:,:,oindex3),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex3,oindex1,oindex2,no,nv,ccsd_doubles_3,&
                            & ovoo_tile_21,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
#endif

    ! ikj,kji
    call trip_amplitudes_ijk_virt(oindex1,oindex3,oindex2,no,nv,ccsd_doubles_1(:,:,oindex3),&
                            & vvvo_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex3,oindex2,oindex1,no,nv,ccsd_doubles_3,&
                            & ovoo_tile_12,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
#endif 

    ! kji,jik
    call trip_amplitudes_ijk_virt(oindex3,oindex2,oindex1,no,nv,ccsd_doubles_3(:,:,oindex2),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex2,oindex1,oindex3,no,nv,ccsd_doubles_2,&
                            & ovoo_tile_31,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_ijk_case3


  !> \brief: generator for triples amplitudes, case(3)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_abc_case3(vindex1,vindex2,vindex3,no,nv,ccsd_doubles_1,ccsd_doubles_2,ccsd_doubles_3,&
                                & ooov_tile_1,ooov_tile_2,ooov_tile_3,&
                                & vovv_12,vovv_13,vovv_21,vovv_23,vovv_31,vovv_32,&
                                & trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> a, b, a, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,vindex3,no,nv
    real(real_pt), dimension(no,no,nv), intent(inout)  :: ccsd_doubles_1,ccsd_doubles_2,ccsd_doubles_3
    !> vovv integrals
    real(real_pt), dimension(nv,no), intent(inout)  :: vovv_12,vovv_13,vovv_21,vovv_23,vovv_31,vovv_32
    !> tiles of ooov 2-el integrals 
    real(real_pt), dimension(no,no,no), intent(inout) :: ooov_tile_1, ooov_tile_2, ooov_tile_3
    !> triples amplitude and work array
    real(real_pt), dimension(no,no,no) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! abc.bca
    call trip_amplitudes_abc_occ(vindex1,vindex2,vindex3,no,nv,ccsd_doubles_1(:,:,vindex2),&
                            & ooov_tile_3,trip_ampl,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex2,vindex3,vindex1,no,nv,ccsd_doubles_2,&
                            & vovv_13,trip_ampl,handle,cublas_handle)

    ! cab,abc
    call trip_amplitudes_abc_occ(vindex3,vindex1,vindex2,no,nv,ccsd_doubles_3(:,:,vindex1),&
                            & ooov_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex1,vindex2,vindex3,no,nv,ccsd_doubles_1,&
                            & vovv_32,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
#endif

    ! bca,cab
    call trip_amplitudes_abc_occ(vindex2,vindex3,vindex1,no,nv,ccsd_doubles_2(:,:,vindex3),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex3,vindex1,vindex2,no,nv,ccsd_doubles_3,&
                            & vovv_21,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
#endif

    ! acb,cba
    call trip_amplitudes_abc_occ(vindex1,vindex3,vindex2,no,nv,ccsd_doubles_1(:,:,vindex3),&
                            & ooov_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex3,vindex2,vindex1,no,nv,ccsd_doubles_3,&
                            & vovv_12,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
#endif

    ! bac,acb
    call trip_amplitudes_abc_occ(vindex2,vindex1,vindex3,no,nv,ccsd_doubles_2(:,:,vindex1),&
                            & ooov_tile_3,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex1,vindex3,vindex2,no,nv,ccsd_doubles_1,&
                            & vovv_23,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
#endif

    ! cba,bac
    call trip_amplitudes_abc_occ(vindex3,vindex2,vindex1,no,nv,ccsd_doubles_3(:,:,vindex2),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_abc_virt(vindex2,vindex1,vindex3,no,nv,ccsd_doubles_2,&
                            & vovv_31,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_abc_case3


  !> \brief: create VIRTUAL part of a triples amplitude ([a,b,c] tuple) for a fixed [i,j,k] tuple, that is, t^{***}_{ijk}
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  subroutine trip_amplitudes_ijk_virt(oindex1,oindex2,oindex3,no,nv,doub_ampl_v2,int_virt_tile,trip,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(real_pt), dimension(nv,nv,nv), target, intent(in) :: int_virt_tile
    real(real_pt), dimension(nv,nv), target, intent(in) :: doub_ampl_v2
    real(real_pt), dimension(nv,nv,nv), target, intent(inout) :: trip
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    call ls_dgemm_acc('n','n',nv,nv**2,nv,1.0E0_realk,doub_ampl_v2,nv,int_virt_tile,nv,0.0E0_realk,trip,nv,&
                    & int(i8*nv**2,kind=8),int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),&
                    & async_idx,cublas_handle)

  end subroutine trip_amplitudes_ijk_virt


  !> \brief: create OCCUPIED part of a triples amplitude ([i,j,k] tuple) for a fixed [a,b,c] tuple, that is, t^{***}_{abc}
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine trip_amplitudes_abc_occ(vindex1,vindex2,vindex3,no,nv,doub_ampl_o2,int_occ_tile,trip,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(real_pt), dimension(no,no,no), target, intent(in) :: int_occ_tile
    real(real_pt), dimension(no,no), target, intent(in) :: doub_ampl_o2
    real(real_pt), dimension(no,no,no), target, intent(inout) :: trip
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    call ls_dgemm_acc('n','n',no,no**2,no,-1.0E0_realk,doub_ampl_o2,no,int_occ_tile,no,0.0E0_realk,trip,no,&
                    & int(i8*no**2,kind=8),int((i8*no**2)*no,kind=8),int((i8*no**2)*no,kind=8),&
                    & async_idx,cublas_handle)

  end subroutine trip_amplitudes_abc_occ


  !> \brief: create OCCUPIED part of a triples amplitude ([a,b,c] tuple) for a fixed [i,j,k] tuple, that is, t^{***}_{ijk}
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  subroutine trip_amplitudes_ijk_occ(oindex1,oindex2,oindex3,no,nv,doub_ampl_ov2,int_occ_portion,trip,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(real_pt), dimension(no,nv), target, intent(in) :: int_occ_portion
    real(real_pt), dimension(nv,nv,no), target, intent(in) :: doub_ampl_ov2
    real(real_pt), dimension(nv,nv,nv), target, intent(inout) :: trip
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    call ls_dgemm_acc('t','t',nv,nv**2,no,-1.0E0_realk,int_occ_portion,no,doub_ampl_ov2,nv**2,1.0E0_realk,trip,nv,&
                    & int(i8*nv*no,kind=8),int((i8*nv**2)*no,kind=8),int((i8*nv**2)*nv,kind=8),&
                    & async_idx,cublas_handle)

  end subroutine trip_amplitudes_ijk_occ


  !> \brief: create VIRTUAL part of a triples amplitude ([i,j,k] tuple) for a fixed [a,b,c] tuple, that is, t^{***}_{abc}
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine trip_amplitudes_abc_virt(vindex1,vindex2,vindex3,no,nv,doub_ampl_vo2,int_virt_portion,trip,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(real_pt), dimension(nv,no), target, intent(in) :: int_virt_portion
    real(real_pt), dimension(no,no,nv), target, intent(in) :: doub_ampl_vo2
    real(real_pt), dimension(no,no,no), target, intent(inout) :: trip
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    call ls_dgemm_acc('t','t',no,no**2,nv,1.0E0_realk,int_virt_portion,nv,doub_ampl_vo2,no**2,1.0E0_realk,trip,no,&
                    & int(i8*nv*no,kind=8),int((i8*no**2)*nv,kind=8),int((i8*no**2)*no,kind=8),&
                    & async_idx,cublas_handle)

  end subroutine trip_amplitudes_abc_virt

  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  !> \param: oindex1, oindex2, and oindex3 are the three occupied indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: eigenocc and eigenvirt are vectors containing occupied and virtual orbital energies, respectively
  !> \param: amplitudes are the final triples amplitude tuple [a,b,c], that is, of the size (virt)³ kept in memory
  subroutine trip_denom_ijk(oindex1,oindex2,oindex3,no,nv,eigenocc,eigenvirt,trip,async_idx)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    real(real_pt), dimension(nv,nv,nv), intent(inout) :: trip
    !> temporary quantities
    integer :: a, b, c
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    real(realk) :: e_orb_occ,e_denom

    ! at first, calculate the sum of the three participating occupied orbital energies, as this
    ! is a constant for the three incomming occupied indices

    e_orb_occ = eigenocc(oindex1) + eigenocc(oindex2) + eigenocc(oindex3)

#ifdef VAR_OPENACC
!$acc parallel present(trip,eigenvirt) firstprivate(nv,e_orb_occ) &
!$acc& private(a,b,c,e_denom) async(async_idx)
#else
!$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,b,c,e_denom),SHARED(nv,trip,eigenvirt,e_orb_occ)
#endif
#ifdef VAR_OPENACC
!$acc loop gang
#endif
    do a=1,nv
#ifdef VAR_OPENACC
!$acc loop worker
#endif
       do b=1,nv
#ifdef VAR_OPENACC
!$acc loop vector
#endif
          do c=1,nv

                  e_denom = e_orb_occ - eigenvirt(a) - eigenvirt(b) - eigenvirt(c)

#ifdef VAR_REAL_SP
                  trip(c,b,a) = trip(c,b,a) / real(e_denom,kind=4)
#else
                  trip(c,b,a) = trip(c,b,a) / e_denom
#endif

          end do
#ifdef VAR_OPENACC
!$acc end loop
#endif
       end do
#ifdef VAR_OPENACC
!$acc end loop
#endif
    end do
#ifdef VAR_OPENACC
!$acc end loop
!$acc end parallel
#else
!$OMP END PARALLEL DO
#endif

  end subroutine trip_denom_ijk


  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  !> \param: vindex1, vindex2, and vindex3 are the three virtual indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: eigenocc and eigenvirt are vectors containing occupied and virtual orbital energies, respectively
  !> \param: amplitudes are the final triples amplitude tuple [i,j,k], that is, of the size (occ)³ kept in memory
  subroutine trip_denom_abc(vindex1,vindex2,vindex3,no,nv,eigenocc,eigenvirt,trip,async_idx)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    real(real_pt), dimension(no,no,no), intent(inout) :: trip
    !> temporary quantities
    integer :: i, j, k
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    real(realk) :: e_orb_virt,e_denom

    ! at first, calculate the sum of the three participating virtual orbital energies, as this
    ! is a constant for the three incomming occupied indices

    e_orb_virt = eigenvirt(vindex1) + eigenvirt(vindex2) + eigenvirt(vindex3)

#ifdef VAR_OPENACC
!$acc parallel present(trip,eigenocc) firstprivate(no,e_orb_virt)&
!$acc& private(i,j,k,e_denom) async(async_idx)
#else
!$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,j,k,e_denom),SHARED(no,trip,eigenocc,e_orb_virt)
#endif
#ifdef VAR_OPENACC
!$acc loop gang
#endif
    do i=1,no
#ifdef VAR_OPENACC
!$acc loop worker
#endif
       do j=1,no
#ifdef VAR_OPENACC
!$acc loop vector
#endif
          do k=1,no

                  e_denom = eigenocc(i) + eigenocc(j) + eigenocc(k) - e_orb_virt

#ifdef VAR_REAL_SP
                  trip(k,j,i) = trip(k,j,i) / real(e_denom,kind=4)
#else
                  trip(k,j,i) = trip(k,j,i) / e_denom
#endif

          end do
#ifdef VAR_OPENACC
!$acc end loop
#endif
       end do
#ifdef VAR_OPENACC
!$acc end loop
#endif
    end do
#ifdef VAR_OPENACC
!$acc end loop
!$acc end parallel
#else
!$OMP END PARALLEL DO
#endif

  end subroutine trip_denom_abc


  subroutine ccsdpt_energy_full_ijk_case1(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & vvoo_tile_12,vvoo_tile_13,vvoo_tile_31,&
                                         & e4,e5,tmp_res,t1_1,t1_3,&
                                         & async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: o1,o2,o3,no,nv
    !> trip arrays
    real(real_pt), dimension(nv,nv,nv), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4,e5
    !> ccsd(t) singles amplitudes
    real(real_pt), dimension(nv) :: tmp_res,t1_1,t1_3
    !> tiles of vvoo integrals
    real(real_pt), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_13, vvoo_tile_31
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp energy
    real(real_pt) :: e4_tmp(1),e5_tmp(1)

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

!$acc enter data create(e4_tmp,e5_tmp) async(handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,(i8*nv)*(i8*nv**2))
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 2.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o1,o1,o3,nv,no,vvoo_tile_13,vvoo_tile_31,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o1,o1,o3,nv,no,vvoo_tile_12,vvoo_tile_12,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - dble(e4_tmp(1))
!$acc end kernels

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o3,o1,o1,nv,no,vvoo_tile_12,vvoo_tile_12,&
                 & tmp_res,t1_3,trip_ampl,.true.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o3,o1,o1,nv,no,vvoo_tile_13,vvoo_tile_31,&
                 & tmp_res,t1_1,trip_ampl,.true.,e5,e5_tmp,handle,cublas_handle)

!$acc exit data delete(e4_tmp,e5_tmp) async(handle)

  end subroutine ccsdpt_energy_full_ijk_case1


  subroutine ccsdpt_energy_full_abc_case1(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & oovv_tile_12,oovv_tile_13,oovv_tile_31,&
                                         & e4,e5,tmp_res,t1_1,t1_3,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: v1,v2,v3,no,nv
    !> trip arrays
    real(real_pt), dimension(no,no,no), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4,e5
    !> ccsd(t) singles amplitudes
    real(real_pt), dimension(no) :: tmp_res,t1_1,t1_3
    !> tiles of vvoo integrals
    real(real_pt), dimension(no,no) :: oovv_tile_12, oovv_tile_13, oovv_tile_31
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp energy
    real(real_pt) :: e4_tmp(1),e5_tmp(1)

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

!$acc enter data create(e4_tmp,e5_tmp) async(handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,(i8*no)*no**2)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 2.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v1,v1,v3,nv,no,oovv_tile_13,oovv_tile_31,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v1,v1,v3,nv,no,oovv_tile_12,oovv_tile_12,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - dble(e4_tmp(1))
!$acc end kernels

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v3,v1,v1,nv,no,oovv_tile_12,oovv_tile_12,&
                 & tmp_res,t1_3,trip_ampl,.true.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v3,v1,v1,nv,no,oovv_tile_13,oovv_tile_31,&
                 & tmp_res,t1_1,trip_ampl,.true.,e5,e5_tmp,handle,cublas_handle)

!$acc exit data delete(e4_tmp,e5_tmp) async(handle)

  end subroutine ccsdpt_energy_full_abc_case1


  subroutine ccsdpt_energy_full_ijk_case2(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & vvoo_tile_12,vvoo_tile_21,vvoo_tile_23,&
                                         & e4,e5,tmp_res,t1_1,t1_2,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: o1,o2,o3,no,nv
    !> trip arrays
    real(real_pt), dimension(nv,nv,nv), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4,e5
    !> ccsd(t) singles amplitudes
    real(real_pt), dimension(nv) :: tmp_res,t1_1,t1_2
    !> tiles of vvoo integrals
    real(real_pt), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_21, vvoo_tile_23
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp energy
    real(real_pt) :: e4_tmp(1),e5_tmp(1)

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

!$acc enter data create(e4_tmp,e5_tmp) async(handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,(i8*nv)*(i8*nv**2))
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 2.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o1,o2,o2,nv,no,vvoo_tile_23,vvoo_tile_23,&
                 & tmp_res,t1_1,trip_ampl,.true.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o1,o2,o2,nv,no,vvoo_tile_21,vvoo_tile_12,&
                 & tmp_res,t1_2,trip_ampl,.true.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o2,o2,o1,nv,no,vvoo_tile_21,vvoo_tile_12,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o2,o2,o1,nv,no,vvoo_tile_23,vvoo_tile_23,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - dble(e4_tmp(1))
!$acc end kernels

!$acc exit data delete(e4_tmp,e5_tmp) async(handle)

  end subroutine ccsdpt_energy_full_ijk_case2


  subroutine ccsdpt_energy_full_abc_case2(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & oovv_tile_12,oovv_tile_21,oovv_tile_23,&
                                         & e4,e5,tmp_res,t1_1,t1_2,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: v1,v2,v3,no,nv
    !> trip arrays
    real(real_pt), dimension(no,no,no), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4,e5
    !> ccsd(t) singles amplitudes
    real(real_pt), dimension(no) :: tmp_res,t1_1,t1_2
    !> tiles of vvoo integrals
    real(real_pt), dimension(no,no) :: oovv_tile_12, oovv_tile_21, oovv_tile_23
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp energy
    real(real_pt) :: e4_tmp(1),e5_tmp(1)

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

!$acc enter data create(e4_tmp,e5_tmp) async(handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,(i8*no)*no**2)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 2.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v1,v2,v2,nv,no,oovv_tile_23,oovv_tile_23,&
                 & tmp_res,t1_1,trip_ampl,.true.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v1,v2,v2,nv,no,oovv_tile_21,oovv_tile_12,&
                 & tmp_res,t1_2,trip_ampl,.true.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v2,v2,v1,nv,no,oovv_tile_21,oovv_tile_12,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v2,v2,v1,nv,no,oovv_tile_23,oovv_tile_23,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - dble(e4_tmp(1))
!$acc end kernels

!$acc exit data delete(e4_tmp,e5_tmp) async(handle)

  end subroutine ccsdpt_energy_full_abc_case2


  subroutine ccsdpt_energy_full_ijk_case3(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & vvoo_tile_12,vvoo_tile_13,vvoo_tile_21,vvoo_tile_23,vvoo_tile_31,vvoo_tile_32,&
                                         & e4,e5,tmp_res,t1_1,t1_2,t1_3,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: o1,o2,o3,no,nv
    !> trip arrays
    real(real_pt), dimension(nv,nv,nv), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4,e5
    !> ccsd(t) singles amplitudes
    real(real_pt), dimension(nv) :: tmp_res,t1_1,t1_2,t1_3
    !> tiles of vvoo integrals
    real(real_pt), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_13, vvoo_tile_21
    real(real_pt), dimension(nv,nv) :: vvoo_tile_23, vvoo_tile_31, vvoo_tile_32
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp energy
    real(real_pt) :: e4_tmp(1),e5_tmp(1)

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

!$acc enter data create(e4_tmp,e5_tmp) async(handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,(i8*nv)*(i8*nv**2))
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 8.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o1,o2,o3,nv,no,vvoo_tile_23,vvoo_tile_32,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o1,o2,o3,nv,no,vvoo_tile_21,vvoo_tile_12,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 2.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o2,o3,o1,nv,no,vvoo_tile_31,vvoo_tile_13,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o2,o3,o1,nv,no,vvoo_tile_32,vvoo_tile_23,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 2.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o3,o1,o2,nv,no,vvoo_tile_12,vvoo_tile_21,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o3,o1,o2,nv,no,vvoo_tile_13,vvoo_tile_31,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[3,2,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[3,2,1],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - 4.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o3,o2,o1,nv,no,vvoo_tile_21,vvoo_tile_12,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o3,o2,o1,nv,no,vvoo_tile_23,vvoo_tile_32,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[1,3,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[1,3,2],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - 4.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o1,o3,o2,nv,no,vvoo_tile_32,vvoo_tile_23,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o1,o3,o2,nv,no,vvoo_tile_31,vvoo_tile_13,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[2,1,3],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[2,1,3],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_ijk(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*nv**2)*nv,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - 4.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_ijk_11_full(o2,o1,o3,nv,no,vvoo_tile_13,vvoo_tile_31,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_ijk_12_full(o2,o1,o3,nv,no,vvoo_tile_12,vvoo_tile_21,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

!$acc exit data delete(e4_tmp,e5_tmp) async(handle)

  end subroutine ccsdpt_energy_full_ijk_case3


  subroutine ccsdpt_energy_full_abc_case3(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & oovv_tile_12,oovv_tile_13,oovv_tile_21,oovv_tile_23,oovv_tile_31,oovv_tile_32,&
                                         & e4,e5,tmp_res,t1_1,t1_2,t1_3,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: v1,v2,v3,no,nv
    !> trip arrays
    real(real_pt), dimension(no,no,no), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4,e5
    !> ccsd(t) singles amplitudes
    real(real_pt), dimension(no) :: tmp_res,t1_1,t1_2,t1_3
    !> tiles of vvoo integrals
    real(real_pt), dimension(no,no) :: oovv_tile_12,oovv_tile_13,oovv_tile_21,oovv_tile_23,oovv_tile_31,oovv_tile_32
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp energy
    real(real_pt) :: e4_tmp(1),e5_tmp(1)

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

!$acc enter data create(e4_tmp,e5_tmp) async(handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC) && !defined(VAR_REAL_SP)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,(i8*no)*no**2)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 8.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v1,v2,v3,nv,no,oovv_tile_23,oovv_tile_32,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v1,v2,v3,nv,no,oovv_tile_21,oovv_tile_12,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 2.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v2,v3,v1,nv,no,oovv_tile_31,oovv_tile_13,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v2,v3,v1,nv,no,oovv_tile_32,oovv_tile_23,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 + 2.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v3,v1,v2,nv,no,oovv_tile_12,oovv_tile_21,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v3,v1,v2,nv,no,oovv_tile_13,oovv_tile_31,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[3,2,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[3,2,1],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - 4.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v3,v2,v1,nv,no,oovv_tile_21,oovv_tile_12,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v3,v2,v1,nv,no,oovv_tile_23,oovv_tile_32,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[1,3,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[1,3,2],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - 4.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v1,v3,v2,nv,no,oovv_tile_32,oovv_tile_23,&
                 & tmp_res,t1_1,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v1,v3,v2,nv,no,oovv_tile_31,oovv_tile_13,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[2,1,3],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[2,1,3],0.0E0_realk,trip_ampl)
#endif

    call trip_denom_abc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)

    call ls_ddot_acc(int((i8*no**2)*no,kind=8),trip_tmp,1,trip_ampl,1,e4_tmp,handle,cublas_handle)
!$acc kernels present(e4,e4_tmp) async(handle)
    e4 = e4 - 4.0E0_realk * dble(e4_tmp(1))
!$acc end kernels

    call ccsdpt_contract_abc_11_full(v2,v1,v3,nv,no,oovv_tile_13,oovv_tile_31,&
                 & tmp_res,t1_2,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)
    call ccsdpt_contract_abc_12_full(v2,v1,v3,nv,no,oovv_tile_12,oovv_tile_21,&
                 & tmp_res,t1_3,trip_ampl,.false.,e5,e5_tmp,handle,cublas_handle)

!$acc exit data delete(e4_tmp,e5_tmp) async(handle)

  end subroutine ccsdpt_energy_full_abc_case3

  subroutine ccsdpt_contract_ijk_11_full(oindex1,oindex2,oindex3,nv,no,int_normal_23,int_normal_32,tmp_res,t1,&
                              & trip_ampl,special,e5,e5_tmp,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(real_pt), dimension(nv), target :: tmp_res,t1
    real(real_pt), dimension(nv,nv), target :: int_normal_23, int_normal_32
    real(real_pt), dimension(nv,nv,nv), target :: trip_ampl
    real(real_pt), target :: e5_tmp(1)
    real(realk) :: e5
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_ijk_11_full: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal and we thus only need to contract with the coulumb part.
       call ls_dgemm_acc('n','n',nv,1,nv**2,1.0E0_realk,trip_ampl,nv,int_normal_23,nv**2,0.0E0_realk,tmp_res,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*nv,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
       call ls_dgemm_acc('n','n',nv,1,nv**2,2.0E0_realk,trip_ampl,nv,int_normal_23,nv**2,0.0E0_realk,tmp_res,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*nv,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels
       call ls_dgemm_acc('n','n',nv,1,nv**2,-1.0E0_realk,trip_ampl,nv,int_normal_32,nv**2,0.0E0_realk,tmp_res,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*nv,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels

    end select TypeofContraction_ijk_11_full

  end subroutine ccsdpt_contract_ijk_11_full

  subroutine ccsdpt_contract_abc_11_full(vindex1,vindex2,vindex3,nv,no,int_normal_23,int_normal_32,tmp_res,t1,&
                              & trip_ampl,special,e5,e5_tmp,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(real_pt), dimension(no), target :: tmp_res,t1
    real(real_pt), dimension(no,no), target :: int_normal_23, int_normal_32
    real(real_pt), dimension(no,no,no), target :: trip_ampl
    real(real_pt), target :: e5_tmp(1)
    real(realk) :: e5
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_abc_11_full: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal and we thus only need to contract with the coulumb part. 
       call ls_dgemm_acc('n','n',no,1,no**2,1.0E0_realk,trip_ampl,no,int_normal_23,no**2,0.0E0_realk,tmp_res,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*no,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
       call ls_dgemm_acc('n','n',no,1,no**2,2.0E0_realk,trip_ampl,no,int_normal_23,no**2,0.0E0_realk,tmp_res,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*no,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels
       call ls_dgemm_acc('n','n',no,1,no**2,-1.0E0_realk,trip_ampl,no,int_normal_32,no**2,0.0E0_realk,tmp_res,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*no,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels

    end select TypeofContraction_abc_11_full

  end subroutine ccsdpt_contract_abc_11_full

  subroutine ccsdpt_contract_ijk_12_full(oindex1,oindex2,oindex3,nv,no,int_normal_21,int_normal_12,tmp_res,t1,&
                              & trip_ampl,special,e5,e5_tmp,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(real_pt), dimension(nv), target :: tmp_res,t1
    real(real_pt), dimension(nv,nv), target :: int_normal_21, int_normal_12
    real(real_pt), dimension(nv,nv,nv), target :: trip_ampl
    real(real_pt), target :: e5_tmp(1)
    real(realk) :: e5
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_ijk_12_full: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal
       ! and we thus only need to contract with (-1)*coulumb part. 
       ! now contract coulumb term over both indices
       call ls_dgemm_acc('n','n',nv,1,nv**2,-1.0E0_realk,trip_ampl,nv,int_normal_21,nv**2,0.0E0_realk,tmp_res,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*nv,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
       call ls_dgemm_acc('n','n',nv,1,nv**2,-2.0E0_realk,trip_ampl,nv,int_normal_21,nv**2,0.0E0_realk,tmp_res,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*nv,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels
       call ls_dgemm_acc('n','n',nv,1,nv**2,1.0E0_realk,trip_ampl,nv,int_normal_12,nv**2,0.0E0_realk,tmp_res,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*nv,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels

    end select TypeofContraction_ijk_12_full

  end subroutine ccsdpt_contract_ijk_12_full

  subroutine ccsdpt_contract_abc_12_full(vindex1,vindex2,vindex3,nv,no,int_normal_21,int_normal_12,tmp_res,t1,&
                              & trip_ampl,special,e5,e5_tmp,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(real_pt), dimension(no), target :: tmp_res,t1
    real(real_pt), dimension(no,no), target :: int_normal_21, int_normal_12
    real(real_pt), dimension(no,no,no), target :: trip_ampl
    real(real_pt), target :: e5_tmp(1)
    real(realk) :: e5
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_abc_12_full: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal
       ! and we thus only need to contract with (-1)*coulumb part. 

       ! now contract coulumb term over both indices
       call ls_dgemm_acc('n','n',no,1,no**2,-1.0E0_realk,trip_ampl,no,int_normal_21,no**2,0.0E0_realk,tmp_res,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*no,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
       call ls_dgemm_acc('n','n',no,1,no**2,-2.0E0_realk,trip_ampl,no,int_normal_21,no**2,0.0E0_realk,tmp_res,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*no,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels
       call ls_dgemm_acc('n','n',no,1,no**2,1.0E0_realk,trip_ampl,no,int_normal_12,no**2,0.0E0_realk,tmp_res,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)
       call ls_ddot_acc(int(i8*no,kind=8),tmp_res,1,t1,1,e5_tmp,async_idx,cublas_handle)
!$acc kernels present(e5,e5_tmp) async(async_idx)
       e5 = e5 + 2.0E0_realk * dble(e5_tmp(1))
!$acc end kernels

    end select TypeofContraction_abc_12_full

  end subroutine ccsdpt_contract_abc_12_full

!endif mod_unreleased
#endif

  subroutine dummy_ccsdpt_full_routine()

  end subroutine dummy_ccsdpt_full_routine

end module ccsdpt_full_module
