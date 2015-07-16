!> @file
!> DEC-CCSD(T) kernels
!> \brief: ccsd(t) module
!> \author: Janus Juul Eriksen
!> \date: 2012-2014, Aarhus
module ccsdpt_dec_module
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
  public :: ccsdpt_driver_ijk_case1_par
  public :: ccsdpt_driver_ijk_case1_ser
  public :: ccsdpt_driver_ijk_case2_par
  public :: ccsdpt_driver_ijk_case2_ser
  public :: ccsdpt_driver_ijk_case3_par
  public :: ccsdpt_driver_ijk_case3_ser
  public :: ccsdpt_driver_abc_case1_par
  public :: ccsdpt_driver_abc_case1_ser
  public :: ccsdpt_driver_abc_case2_par
  public :: ccsdpt_driver_abc_case2_ser
  public :: ccsdpt_driver_abc_case3_par
  public :: ccsdpt_driver_abc_case3_ser
#endif

contains

#ifdef MOD_UNRELEASED


#ifdef VAR_MPI
  !> \brief: driver routine for contractions in case(1) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case1_par(oindex1,oindex3,no,nv,vvoo_tile_12,vvoo_tile_13,vvoo_tile_31,&
                            & ovoo_tile_12,ovoo_tile_13,ovoo_tile_31,&
                            & vvvo_tile_o1,vvvo_tile_o3,&
                            & ccsdpt_singles_1,ccsdpt_singles_3,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_3,wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & i,k,tile_size_i,tile_size_k)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex3,no,nv
    integer, intent(in) :: i,k,tile_size_i,tile_size_k
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no) :: ccsdpt_doubles_1,ccsdpt_doubles_3
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_3
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_13, ovoo_tile_31
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv,tile_size_i,tile_size_i) :: vvoo_tile_12
    real(realk), dimension(nv,nv,tile_size_i,tile_size_k) :: vvoo_tile_13
    real(realk), dimension(nv,nv,tile_size_k,tile_size_i) :: vvoo_tile_31
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex3
    real(realk), dimension(nv,nv,nv,tile_size_i), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv,tile_size_k), intent(inout) :: vvvo_tile_o3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : iik --132--> iki --231--> kii
    ! in 211/212 : kii ........ kii ........ iik
    ! in 221/222 : iki ........ iik ........ iki

    do idx = 1,3

       if (idx .eq. 1) then ! iik

          ! calculate contribution to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex1,oindex3,nv,no,vvoo_tile_13(:,:,i,k),vvoo_tile_31(:,:,k,i),&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex1,oindex3,nv,no,vvoo_tile_12(:,:,i,i),vvoo_tile_12(:,:,i,i),&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex1,oindex1,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex3),&
                           & wrk_3d,trip,vvvo_tile_o1(:,:,:,i),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex1,oindex1,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o3(:,:,:,k),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex3,oindex1,no,nv,ovoo_tile_31,ovoo_tile_13,&
                           & ccsdpt_doubles_1,trip,wrk_3d,.true.,handle,cublas_handle)

       else if (idx .eq. 2) then ! kii

          ! iki: this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_iki trip amplitudes and therefore end up
          ! canceling each other when added to ccsdpt_singles

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex1),&
                           & trip,wrk_3d,vvvo_tile_o3(:,:,:,k),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex3),&
                           & trip,wrk_3d,vvvo_tile_o1(:,:,:,i),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex1,oindex1,no,nv,ovoo_tile_12,ovoo_tile_12,&
                           & ccsdpt_doubles_3,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex1,oindex1,no,nv,ovoo_tile_31,&
                           & ccsdpt_doubles_1,wrk_3d,trip,handle,cublas_handle)

       else if (idx .eq. 3) then ! iki

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,3,1],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,3,1],0.0E0_realk,trip)
#endif

          call ccsdpt_contract_ijk_11(oindex3,oindex1,oindex1,nv,no,vvoo_tile_12(:,:,i,i),vvoo_tile_12(:,:,i,i),&
                       & ccsdpt_singles_3,&
                       & trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex1,oindex1,nv,no,vvoo_tile_13(:,:,i,k),vvoo_tile_31(:,:,k,i),&
                       & ccsdpt_singles_1,&
                       & trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o1(:,:,:,i),.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex1,oindex3,no,nv,ovoo_tile_13,ovoo_tile_31,&
                           & ccsdpt_doubles_1,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex1,oindex3,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_3,trip,wrk_3d,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case1_par
#endif


  !> \brief: driver routine for contractions in case(1) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case1_ser(oindex1,oindex3,no,nv,vvoo_tile_12,vvoo_tile_13,vvoo_tile_31,&
                            & ovoo_tile_12,ovoo_tile_13,ovoo_tile_31,&
                            & vvvo_tile_o1,vvvo_tile_o3,&
                            & ccsdpt_singles_1,ccsdpt_singles_3,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_3,wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex3,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no) :: ccsdpt_doubles_1,ccsdpt_doubles_3
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_3
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_13, ovoo_tile_31
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_13, vvoo_tile_31
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex3
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : iik --132--> iki --231--> kii
    ! in 211/212 : kii ........ kii ........ iik
    ! in 221/222 : iki ........ iik ........ iki

    do idx = 1,3

       if (idx .eq. 1) then ! iik

          ! calculate contribution to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex1,oindex3,nv,no,vvoo_tile_13,vvoo_tile_31,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex1,oindex3,nv,no,vvoo_tile_12,vvoo_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex1,oindex1,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex3),&
                           & wrk_3d,trip,vvvo_tile_o1,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex1,oindex1,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o3,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex3,oindex1,no,nv,ovoo_tile_31,ovoo_tile_13,&
                           & ccsdpt_doubles_1,trip,wrk_3d,.true.,handle,cublas_handle)

       else if (idx .eq. 2) then ! kii

          ! iki: this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_iki trip amplitudes and therefore end up
          ! canceling each other when added to ccsdpt_singles

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex1),&
                           & trip,wrk_3d,vvvo_tile_o3,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex3),&
                           & trip,wrk_3d,vvvo_tile_o1,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex1,oindex1,no,nv,ovoo_tile_12,ovoo_tile_12,&
                           & ccsdpt_doubles_3,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex1,oindex1,no,nv,ovoo_tile_31,&
                           & ccsdpt_doubles_1,wrk_3d,trip,handle,cublas_handle)

       else if (idx .eq. 3) then ! iki

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,3,1],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,3,1],0.0E0_realk,trip)
#endif

          call ccsdpt_contract_ijk_11(oindex3,oindex1,oindex1,nv,no,vvoo_tile_12,vvoo_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex1,oindex1,nv,no,vvoo_tile_13,vvoo_tile_31,&
                       & ccsdpt_singles_1,&
                       & trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o1,.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex1,oindex3,no,nv,ovoo_tile_13,ovoo_tile_31,&
                           & ccsdpt_doubles_1,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex1,oindex3,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_3,trip,wrk_3d,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case1_ser


#ifdef VAR_MPI
  !> \brief: driver routine for contractions in case(1) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case1_par(vindex1,vindex3,no,nv,vovv_tile_1,vovv_tile_3,&
                            & oovv_tile_12,oovv_tile_13,oovv_tile_31,&
                            & ooov_tile_v1,ooov_tile_v3,&
                            & ccsdpt_singles_1,ccsdpt_singles_3,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_3,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & a,c,tile_size_a,tile_size_c)

    implicit none

    !> a, b, c, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex3,no,nv
    integer, intent(in) :: a,c,tile_size_a,tile_size_c
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no,nv) :: ccsdpt_doubles_1,ccsdpt_doubles_3
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_3
    !> tiles of vovv integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout) :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_c), intent(inout) :: vovv_tile_3 
    !> tiles of oovv integrals
    real(realk), dimension(no,no,tile_size_a,tile_size_a), intent(inout) :: oovv_tile_12
    real(realk), dimension(no,no,tile_size_a,tile_size_c), intent(inout) :: oovv_tile_13
    real(realk), dimension(no,no,tile_size_c,tile_size_a), intent(inout) :: oovv_tile_31
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex3
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v3
    !> triples amplitude and work array
    real(realk), dimension(nv,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : aac --132--> aca --231--> caa
    ! in 211/212 : caa ........ caa ........ aac
    ! in 221/222 : aca ........ aac ........ aca

    do idx = 1,3

       if (idx .eq. 1) then ! aac

          ! calculate contribution to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex1,vindex3,nv,no,&
                       & oovv_tile_13(:,:,a,c),oovv_tile_31(:,:,c,a),&
                       & ccsdpt_singles_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex1,vindex3,nv,no,&
                       & oovv_tile_12(:,:,a,a),oovv_tile_12(:,:,a,a),&
                       & ccsdpt_singles_3,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex3,vindex1,no,nv,&
                           & vovv_tile_1(:,:,vindex3,a),vovv_tile_3(:,:,vindex1,c),&
                           & ccsdpt_doubles_1,trip,wrk_3d,.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex1,vindex1,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex3),&
                           & wrk_3d,trip,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex1,vindex1,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 2) then ! caa

          ! aca: this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_aca trip amplitudes and therefore end up
          ! canceling each other when added to ccsdpt_singles

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex3,vindex1,vindex1,no,nv,&
                           & vovv_tile_1(:,:,vindex1,a),vovv_tile_1(:,:,vindex1,a),&
                           & ccsdpt_doubles_3,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex1,vindex1,no,nv,vovv_tile_1(:,:,vindex3,a),&
                           & ccsdpt_doubles_1,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex1),&
                           & trip,wrk_3d,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex3),&
                           & trip,wrk_3d,ooov_tile_v1,handle,cublas_handle)

       else if (idx .eq. 3) then ! aca

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,3,1],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,3,1],0.0E0_realk,trip)
#endif

          call ccsdpt_contract_abc_11(vindex3,vindex1,vindex1,nv,no,&
                       & oovv_tile_12(:,:,a,a),oovv_tile_12(:,:,a,a),&
                       & ccsdpt_singles_3,trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex1,vindex1,nv,no,&
                       & oovv_tile_13(:,:,a,c),oovv_tile_31(:,:,c,a),&
                       & ccsdpt_singles_1,trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex1,vindex3,no,nv,&
                           & vovv_tile_3(:,:,vindex1,c),vovv_tile_1(:,:,vindex3,a),&
                           & ccsdpt_doubles_1,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex1,vindex3,no,nv,vovv_tile_1(:,:,vindex1,a),&
                           & ccsdpt_doubles_3,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v1,.true.,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case1_par
#endif


  !> \brief: driver routine for contractions in case(1) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case1_ser(vindex1,vindex3,no,nv,oovv,vovv,&
                            & ooov_tile_v1,ooov_tile_v3,&
                            & ccsdpt_singles_1,ccsdpt_singles_3,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_3,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> a, b, c, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex3,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no,nv) :: ccsdpt_doubles_1,ccsdpt_doubles_3
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_3
    !> vovv integrals
    real(realk), dimension(nv,no,nv,nv), intent(inout)  :: vovv
    !> oovv integrals
    real(realk), dimension(no,no,nv,nv), intent(inout)  :: oovv
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex3
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v3
    !> triples amplitude and work array
    real(realk), dimension(nv,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : aac --132--> aca --231--> caa
    ! in 211/212 : caa ........ caa ........ aac
    ! in 221/222 : aca ........ aac ........ aca

    do idx = 1,3

       if (idx .eq. 1) then ! aac

          ! calculate contribution to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex1,vindex3,nv,no,&
                       & oovv(:,:,vindex1,vindex3),oovv(:,:,vindex3,vindex1),&
                       & ccsdpt_singles_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex1,vindex3,nv,no,&
                       & oovv(:,:,vindex1,vindex1),oovv(:,:,vindex1,vindex1),&
                       & ccsdpt_singles_3,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex3,vindex1,no,nv,&
                           & vovv(:,:,vindex3,vindex1),vovv(:,:,vindex1,vindex3),&
                           & ccsdpt_doubles_1,trip,wrk_3d,.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex1,vindex1,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex3),&
                           & wrk_3d,trip,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex1,vindex1,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 2) then ! caa

          ! aca: this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_aca trip amplitudes and therefore end up
          ! canceling each other when added to ccsdpt_singles

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex3,vindex1,vindex1,no,nv,&
                           & vovv(:,:,vindex1,vindex1),vovv(:,:,vindex1,vindex1),&
                           & ccsdpt_doubles_3,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex1,vindex1,no,nv,vovv(:,:,vindex3,vindex1),&
                           & ccsdpt_doubles_1,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex1),&
                           & trip,wrk_3d,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex3),&
                           & trip,wrk_3d,ooov_tile_v1,handle,cublas_handle)

       else if (idx .eq. 3) then ! aca

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,3,1],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,3,1],0.0E0_realk,trip)
#endif

          call ccsdpt_contract_abc_11(vindex3,vindex1,vindex1,nv,no,&
                       & oovv(:,:,vindex1,vindex1),oovv(:,:,vindex1,vindex1),&
                       & ccsdpt_singles_3,trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex1,vindex1,nv,no,&
                       & oovv(:,:,vindex1,vindex3),oovv(:,:,vindex3,vindex1),&
                       & ccsdpt_singles_1,trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex1,vindex3,no,nv,&
                           & vovv(:,:,vindex1,vindex3),vovv(:,:,vindex3,vindex1),&
                           & ccsdpt_doubles_1,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex1,vindex3,no,nv,vovv(:,:,vindex1,vindex1),&
                           & ccsdpt_doubles_3,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v1,.true.,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case1_ser


#ifdef VAR_MPI
  !> \brief: driver routine for contractions in case(2) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case2_par(oindex1,oindex2,no,nv,vvoo_tile_12,vvoo_tile_21,vvoo_tile_23,&
                            & ovoo_tile_12,ovoo_tile_21,ovoo_tile_23,&
                            & vvvo_tile_o1,vvvo_tile_o2,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_2,wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & i,j,tile_size_i,tile_size_j)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,no,nv
    integer, intent(in) :: i,j,tile_size_i,tile_size_j
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no) :: ccsdpt_doubles_1,ccsdpt_doubles_2
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_2
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_21, ovoo_tile_23
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv,tile_size_i,tile_size_j) :: vvoo_tile_12
    real(realk), dimension(nv,nv,tile_size_j,tile_size_i) :: vvoo_tile_21
    real(realk), dimension(nv,nv,tile_size_j,tile_size_j) :: vvoo_tile_23
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex2
    real(realk), dimension(nv,nv,nv,tile_size_i), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv,tile_size_j), intent(inout) :: vvvo_tile_o2
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : ijj --312--> jij --312--> jji
    ! in 211/212 : jij ........ jji ........ ijj
    ! in 221/222 : jji ........ ijj ........ jij

    do idx = 1,3

       if (idx .eq. 1) then
  
          ! calculate contributions to ccsdpt_singles:
 
          call ccsdpt_contract_ijk_11(oindex1,oindex2,oindex2,nv,no,vvoo_tile_23(:,:,j,j),vvoo_tile_23(:,:,j,j),&
                       & ccsdpt_singles_1,&
                       & trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex2,oindex2,nv,no,vvoo_tile_21(:,:,j,i),vvoo_tile_12(:,:,i,j),&
                       & ccsdpt_singles_2,&
                       & trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex2),&
                           & wrk_3d,trip,vvvo_tile_o2(:,:,:,j),.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex2,oindex1,no,nv,ovoo_tile_21,ovoo_tile_12,&
                           & ccsdpt_doubles_2,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex2,oindex1,no,nv,ovoo_tile_23,&
                           & ccsdpt_doubles_1,trip,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 2) then
   
          ! this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_jij trip amplitudes and therefore end up
          ! canceling each other when added to T_star

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d)
#endif 

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex2),&
                           & trip,wrk_3d,vvvo_tile_o1(:,:,:,i),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex1),&
                           & trip,wrk_3d,vvvo_tile_o2(:,:,:,j),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex2,oindex2,no,nv,ovoo_tile_23,ovoo_tile_23,&
                           & ccsdpt_doubles_1,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex2,oindex2,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_2,wrk_3d,trip,handle,cublas_handle)

       else if (idx .eq. 3) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:
   
          call ccsdpt_contract_ijk_11(oindex2,oindex2,oindex1,nv,no,vvoo_tile_21(:,:,j,i),vvoo_tile_12(:,:,i,j),&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex2,oindex1,nv,no,vvoo_tile_23(:,:,j,j),vvoo_tile_23(:,:,j,j),&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
   
          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex2,oindex2,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o2(:,:,:,j),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex2,oindex2,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex2),&
                           & wrk_3d,trip,vvvo_tile_o1(:,:,:,i),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex1,oindex2,no,nv,ovoo_tile_12,ovoo_tile_21,&
                           & ccsdpt_doubles_2,trip,wrk_3d,.true.,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case2_par
#endif


  !> \brief: driver routine for contractions in case(2) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case2_ser(oindex1,oindex2,no,nv,vvoo_tile_12,vvoo_tile_21,vvoo_tile_23,&
                            & ovoo_tile_12,ovoo_tile_21,ovoo_tile_23,&
                            & vvvo_tile_o1,vvvo_tile_o2,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_2,wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no) :: ccsdpt_doubles_1,ccsdpt_doubles_2
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_2
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_21, ovoo_tile_23
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_21, vvoo_tile_23
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex2
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o2
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : ijj --312--> jij --312--> jji
    ! in 211/212 : jij ........ jji ........ ijj
    ! in 221/222 : jji ........ ijj ........ jij

    do idx = 1,3

       if (idx .eq. 1) then
  
          ! calculate contributions to ccsdpt_singles:
 
          call ccsdpt_contract_ijk_11(oindex1,oindex2,oindex2,nv,no,vvoo_tile_23,vvoo_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex2,oindex2,nv,no,vvoo_tile_21,vvoo_tile_12,&
                       & ccsdpt_singles_2,&
                       & trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex2),&
                           & wrk_3d,trip,vvvo_tile_o2,.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex2,oindex1,no,nv,ovoo_tile_21,ovoo_tile_12,&
                           & ccsdpt_doubles_2,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex2,oindex1,no,nv,ovoo_tile_23,&
                           & ccsdpt_doubles_1,trip,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 2) then
   
          ! this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_jij trip amplitudes and therefore end up
          ! canceling each other when added to T_star

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d)
#endif 

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex2),&
                           & trip,wrk_3d,vvvo_tile_o1,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex1),&
                           & trip,wrk_3d,vvvo_tile_o2,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex2,oindex2,no,nv,ovoo_tile_23,ovoo_tile_23,&
                           & ccsdpt_doubles_1,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex2,oindex2,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_2,wrk_3d,trip,handle,cublas_handle)

       else if (idx .eq. 3) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:
   
          call ccsdpt_contract_ijk_11(oindex2,oindex2,oindex1,nv,no,vvoo_tile_21,vvoo_tile_12,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex2,oindex1,nv,no,vvoo_tile_23,vvoo_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
   
          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex2,oindex2,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o2,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex2,oindex2,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex2),&
                           & wrk_3d,trip,vvvo_tile_o1,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex1,oindex2,no,nv,ovoo_tile_12,ovoo_tile_21,&
                           & ccsdpt_doubles_2,trip,wrk_3d,.true.,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case2_ser


#ifdef VAR_MPI
  !> \brief: driver routine for contractions in case(2) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case2_par(vindex1,vindex2,no,nv,vovv_tile_1,vovv_tile_2,&
                            & oovv_tile_12,oovv_tile_21,oovv_tile_23,&
                            & ooov_tile_v1,ooov_tile_v2,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_2,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & a,b,tile_size_a,tile_size_b)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,no,nv
    integer, intent(in) :: a,b,tile_size_a,tile_size_b
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no,nv) :: ccsdpt_doubles_1,ccsdpt_doubles_2
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_2
    !> tiles of vovv integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout) :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_b), intent(inout) :: vovv_tile_2
    !> tiles of oovv integrals
    real(realk), dimension(no,no,tile_size_a,tile_size_b), intent(inout) :: oovv_tile_12
    real(realk), dimension(no,no,tile_size_b,tile_size_a), intent(inout) :: oovv_tile_21
    real(realk), dimension(no,no,tile_size_b,tile_size_b), intent(inout) :: oovv_tile_23
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex2
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v2
    !> triples amplitude and work array
    real(realk), dimension(nv,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : abb --312--> bab --312--> bba
    ! in 211/212 : bab ........ bba ........ abb
    ! in 221/222 : bba ........ abb ........ bab

    do idx = 1,3

       if (idx .eq. 1) then
  
          ! calculate contributions to ccsdpt_singles:
 
          call ccsdpt_contract_abc_11(vindex1,vindex2,vindex2,nv,no,&
                       & oovv_tile_23(:,:,b,b),oovv_tile_23(:,:,b,b),&
                       & ccsdpt_singles_1,trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex2,vindex2,nv,no,&
                       & oovv_tile_21(:,:,b,a),oovv_tile_12(:,:,a,b),&
                       & ccsdpt_singles_2,trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex2,vindex2,vindex1,no,nv,&
                           & vovv_tile_1(:,:,vindex2,a),vovv_tile_2(:,:,vindex1,b),&
                           & ccsdpt_doubles_2,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex2,vindex1,no,nv,vovv_tile_2(:,:,vindex2,b),&
                           & ccsdpt_doubles_1,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex2),&
                           & wrk_3d,trip,ooov_tile_v2,.true.,handle,cublas_handle)

       else if (idx .eq. 2) then
   
          ! this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_bab trip amplitudes and therefore end up
          ! canceling each other when added to T_star

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[3,1,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[3,1,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex2,vindex2,no,nv,&
                           & vovv_tile_2(:,:,vindex2,b),vovv_tile_2(:,:,vindex2,b),&
                           & ccsdpt_doubles_1,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex2,vindex2,no,nv,vovv_tile_2(:,:,vindex1,b),&
                           & ccsdpt_doubles_2,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex2),&
                           & trip,wrk_3d,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex1),&
                           & trip,wrk_3d,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 3) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[3,1,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[3,1,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:
   
          call ccsdpt_contract_abc_11(vindex2,vindex2,vindex1,nv,no,&
                       & oovv_tile_21(:,:,b,a),oovv_tile_12(:,:,a,b),&
                       & ccsdpt_singles_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex2,vindex1,nv,no,&
                       & oovv_tile_23(:,:,b,b),oovv_tile_23(:,:,b,b),&
                       & ccsdpt_singles_1,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex2,vindex1,vindex2,no,nv,&
                           & vovv_tile_2(:,:,vindex1,b),vovv_tile_1(:,:,vindex2,a),&
                           & ccsdpt_doubles_2,trip,wrk_3d,.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex2,vindex2,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex2,vindex2,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex2),&
                           & wrk_3d,trip,ooov_tile_v1,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case2_par
#endif


  !> \brief: driver routine for contractions in case(2) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case2_ser(vindex1,vindex2,no,nv,oovv,vovv,&
                            & ooov_tile_v1,ooov_tile_v2,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_2,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no,nv) :: ccsdpt_doubles_1,ccsdpt_doubles_2
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_2
    !> vovv integrals
    real(realk), dimension(nv,no,nv,nv), intent(inout)  :: vovv
    !> oovv integrals
    real(realk), dimension(no,no,nv,nv), intent(inout)  :: oovv
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex2
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v2
    !> triples amplitude and work array
    real(realk), dimension(nv,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : abb --312--> bab --312--> bba
    ! in 211/212 : bab ........ bba ........ abb
    ! in 221/222 : bba ........ abb ........ bab

    do idx = 1,3

       if (idx .eq. 1) then
  
          ! calculate contributions to ccsdpt_singles:
 
          call ccsdpt_contract_abc_11(vindex1,vindex2,vindex2,nv,no,&
                       & oovv(:,:,vindex2,vindex2),oovv(:,:,vindex2,vindex2),&
                       & ccsdpt_singles_1,trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex2,vindex2,nv,no,&
                       & oovv(:,:,vindex2,vindex1),oovv(:,:,vindex1,vindex2),&
                       & ccsdpt_singles_2,trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex2,vindex2,vindex1,no,nv,&
                           & vovv(:,:,vindex2,vindex1),vovv(:,:,vindex1,vindex2),&
                           & ccsdpt_doubles_2,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex2,vindex1,no,nv,vovv(:,:,vindex2,vindex2),&
                           & ccsdpt_doubles_1,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex2),&
                           & wrk_3d,trip,ooov_tile_v2,.true.,handle,cublas_handle)

       else if (idx .eq. 2) then
   
          ! this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_bab trip amplitudes and therefore end up
          ! canceling each other when added to T_star

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[3,1,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[3,1,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex2,vindex2,no,nv,&
                           & vovv(:,:,vindex2,vindex2),vovv(:,:,vindex2,vindex2),&
                           & ccsdpt_doubles_1,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex2,vindex2,no,nv,vovv(:,:,vindex1,vindex2),&
                           & ccsdpt_doubles_2,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex2),&
                           & trip,wrk_3d,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex1),&
                           & trip,wrk_3d,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 3) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[3,1,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[3,1,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:
   
          call ccsdpt_contract_abc_11(vindex2,vindex2,vindex1,nv,no,&
                       & oovv(:,:,vindex2,vindex1),oovv(:,:,vindex1,vindex2),&
                       & ccsdpt_singles_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex2,vindex1,nv,no,&
                       & oovv(:,:,vindex2,vindex2),oovv(:,:,vindex2,vindex2),&
                       & ccsdpt_singles_1,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex2,vindex1,vindex2,no,nv,&
                           & vovv(:,:,vindex1,vindex2),vovv(:,:,vindex2,vindex1),&
                           & ccsdpt_doubles_2,trip,wrk_3d,.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex2,vindex2,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex2,vindex2,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex2),&
                           & wrk_3d,trip,ooov_tile_v1,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case2_ser


#ifdef VAR_MPI
  !> \brief: driver routine for contractions in case(3) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case3_par(oindex1,oindex2,oindex3,no,nv,&
                            & vvoo_tile_12, vvoo_tile_13, vvoo_tile_21,&
                            & vvoo_tile_23, vvoo_tile_31, vvoo_tile_32,&
                            & ovoo_tile_12, ovoo_tile_13, ovoo_tile_21,&
                            & ovoo_tile_23, ovoo_tile_31, ovoo_tile_32,&
                            & vvvo_tile_o1,vvvo_tile_o2,vvvo_tile_o3,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_2,ccsdpt_doubles_3,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & i,j,k,tile_size_i,tile_size_j,tile_size_k)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,oindex3,no,nv
    integer, intent(in) :: i,j,k,tile_size_i,tile_size_j,tile_size_k
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no) :: ccsdpt_doubles_1,ccsdpt_doubles_2,ccsdpt_doubles_3
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_13, ovoo_tile_21
    real(realk), dimension(no,nv) :: ovoo_tile_23, ovoo_tile_31, ovoo_tile_32
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv,tile_size_i,tile_size_j) :: vvoo_tile_12
    real(realk), dimension(nv,nv,tile_size_i,tile_size_k) :: vvoo_tile_13
    real(realk), dimension(nv,nv,tile_size_j,tile_size_i) :: vvoo_tile_21
    real(realk), dimension(nv,nv,tile_size_j,tile_size_k) :: vvoo_tile_23
    real(realk), dimension(nv,nv,tile_size_k,tile_size_i) :: vvoo_tile_31
    real(realk), dimension(nv,nv,tile_size_k,tile_size_j) :: vvoo_tile_32
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex2,oindex3
    real(realk), dimension(nv,nv,nv,tile_size_i), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv,tile_size_j), intent(inout) :: vvvo_tile_o2
    real(realk), dimension(nv,nv,nv,tile_size_k), intent(inout) :: vvvo_tile_o3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. 
    ! in order to minimize the number of reorderings needed to be performed, 
    ! in order to take optimal advantage of the symmetry of the amplitudes, 
    ! AND to finish all work that involves int_virt_tile_o3 first such that this may be updated (gpu),
    ! we carry out the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : ijk --213--> jik --132--> jki --321--> ikj --213--> kij --132--> kji
    ! in 211/212 : kij ........ ikj ........ ijk ........ kji ........ jki ........ jik
    ! in 221/222 : jki ........ kji ........ kij ........ jik ........ ijk ........ ikj

    do idx = 1,6

       if (idx .eq. 1) then ! ijk

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex2,oindex3,nv,no,vvoo_tile_23(:,:,j,k),vvoo_tile_32(:,:,k,j),&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex2,oindex3,nv,no,vvoo_tile_21(:,:,j,i),vvoo_tile_12(:,:,i,j),&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex3),&
                           & wrk_3d,trip,vvvo_tile_o2(:,:,:,j),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex2),&
                           & wrk_3d,trip,vvvo_tile_o3(:,:,:,k),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex3,oindex1,no,nv,ovoo_tile_31,ovoo_tile_13,&
                           & ccsdpt_doubles_2,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex3,oindex1,no,nv,ovoo_tile_23,&
                           & ccsdpt_doubles_1,trip,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 2) then ! kij

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex2,oindex1,oindex3,nv,no,vvoo_tile_13(:,:,i,k),vvoo_tile_31(:,:,k,i),&
                       & ccsdpt_singles_2,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex1,oindex3,nv,no,vvoo_tile_12(:,:,i,j),vvoo_tile_21(:,:,j,i),&
                       & ccsdpt_singles_3,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex3),&
                           & trip,wrk_3d,vvvo_tile_o1(:,:,:,i),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex1),&
                           & trip,wrk_3d,vvvo_tile_o3(:,:,:,k),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex3,oindex2,no,nv,ovoo_tile_32,ovoo_tile_23,&
                           & ccsdpt_doubles_1,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex3,oindex2,no,nv,ovoo_tile_13,&
                           & ccsdpt_doubles_2,wrk_3d,trip,handle,cublas_handle)

       else if (idx .eq. 3) then ! jki

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex2,oindex3,oindex1,nv,no,vvoo_tile_31(:,:,k,i),vvoo_tile_13(:,:,i,k),&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex3,oindex1,nv,no,vvoo_tile_32(:,:,k,j),vvoo_tile_23(:,:,j,k),&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex2,oindex3,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o3(:,:,:,k),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex2,oindex3,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex3),&
                           & wrk_3d,trip,vvvo_tile_o1(:,:,:,i),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex1,oindex2,no,nv,ovoo_tile_12,ovoo_tile_21,&
                           & ccsdpt_doubles_3,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex1,oindex2,no,nv,ovoo_tile_31,&
                           & ccsdpt_doubles_2,trip,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 4) then ! kji

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,2,1],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,2,1],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex3,oindex2,nv,no,vvoo_tile_32(:,:,k,j),vvoo_tile_23(:,:,j,k),&
                       & ccsdpt_singles_1,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex3,oindex2,nv,no,vvoo_tile_31(:,:,k,i),vvoo_tile_13(:,:,i,k),&
                       & ccsdpt_singles_2,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex2),&
                           & trip,wrk_3d,vvvo_tile_o3(:,:,:,k),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex3),&
                           & trip,wrk_3d,vvvo_tile_o2(:,:,:,j),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex2,oindex1,no,nv,ovoo_tile_21,ovoo_tile_12,&
                           & ccsdpt_doubles_3,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex2,oindex1,no,nv,ovoo_tile_32,&
                           & ccsdpt_doubles_1,wrk_3d,trip,handle,cublas_handle)

       else if (idx .eq. 5) then ! ikj

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex3,oindex1,oindex2,nv,no,vvoo_tile_12(:,:,i,j),vvoo_tile_21(:,:,j,i),&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex1,oindex2,nv,no,vvoo_tile_13(:,:,i,k),vvoo_tile_31(:,:,k,i),&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex2),&
                           & wrk_3d,trip,vvvo_tile_o1(:,:,:,i),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o2(:,:,:,j),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex2,oindex3,no,nv,ovoo_tile_23,ovoo_tile_32,&
                           & ccsdpt_doubles_1,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex2,oindex3,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_3,trip,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 6) then ! jik

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex3,oindex2,oindex1,nv,no,vvoo_tile_21(:,:,j,i),vvoo_tile_12(:,:,i,j),&
                       & ccsdpt_singles_3,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex2,oindex1,nv,no,vvoo_tile_23(:,:,j,k),vvoo_tile_32(:,:,k,j),&
                       & ccsdpt_singles_1,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex3,oindex2,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex1),&
                           & trip,wrk_3d,vvvo_tile_o2(:,:,:,j),.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex3,oindex2,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex2),&
                           & trip,wrk_3d,vvvo_tile_o1(:,:,:,i),handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex1,oindex3,no,nv,ovoo_tile_13,ovoo_tile_31,&
                           & ccsdpt_doubles_2,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex1,oindex3,no,nv,ovoo_tile_21,&
                           & ccsdpt_doubles_3,wrk_3d,trip,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case3_par
#endif


  !> \brief: driver routine for contractions in case(3) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case3_ser(oindex1,oindex2,oindex3,no,nv,&
                            & vvoo_tile_12, vvoo_tile_13, vvoo_tile_21,&
                            & vvoo_tile_23, vvoo_tile_31, vvoo_tile_32,&
                            & ovoo_tile_12, ovoo_tile_13, ovoo_tile_21,&
                            & ovoo_tile_23, ovoo_tile_31, ovoo_tile_32,&
                            & vvvo_tile_o1,vvvo_tile_o2,vvvo_tile_o3,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_2,ccsdpt_doubles_3,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,oindex3,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no) :: ccsdpt_doubles_1,ccsdpt_doubles_2,ccsdpt_doubles_3
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_13, ovoo_tile_21
    real(realk), dimension(no,nv) :: ovoo_tile_23, ovoo_tile_31, ovoo_tile_32
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_13, vvoo_tile_21
    real(realk), dimension(nv,nv) :: vvoo_tile_23, vvoo_tile_31, vvoo_tile_32
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex2,oindex3
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o2
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. 
    ! in order to minimize the number of reorderings needed to be performed, 
    ! in order to take optimal advantage of the symmetry of the amplitudes, 
    ! AND to finish all work that involves int_virt_tile_o3 first such that this may be updated (gpu),
    ! we carry out the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : ijk --213--> jik --132--> jki --321--> ikj --213--> kij --132--> kji
    ! in 211/212 : kij ........ ikj ........ ijk ........ kji ........ jki ........ jik
    ! in 221/222 : jki ........ kji ........ kij ........ jik ........ ijk ........ ikj

    do idx = 1,6

       if (idx .eq. 1) then ! ijk

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex2,oindex3,nv,no,vvoo_tile_23,vvoo_tile_32,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex2,oindex3,nv,no,vvoo_tile_21,vvoo_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex3),&
                           & wrk_3d,trip,vvvo_tile_o2,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex2),&
                           & wrk_3d,trip,vvvo_tile_o3,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex3,oindex1,no,nv,ovoo_tile_31,ovoo_tile_13,&
                           & ccsdpt_doubles_2,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex3,oindex1,no,nv,ovoo_tile_23,&
                           & ccsdpt_doubles_1,trip,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 2) then ! kij

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex2,oindex1,oindex3,nv,no,vvoo_tile_13,vvoo_tile_31,&
                       & ccsdpt_singles_2,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex1,oindex3,nv,no,vvoo_tile_12,vvoo_tile_21,&
                       & ccsdpt_singles_3,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex3),&
                           & trip,wrk_3d,vvvo_tile_o1,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex1),&
                           & trip,wrk_3d,vvvo_tile_o3,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex3,oindex2,no,nv,ovoo_tile_32,ovoo_tile_23,&
                           & ccsdpt_doubles_1,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex3,oindex2,no,nv,ovoo_tile_13,&
                           & ccsdpt_doubles_2,wrk_3d,trip,handle,cublas_handle)

       else if (idx .eq. 3) then ! jki

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex2,oindex3,oindex1,nv,no,vvoo_tile_31,vvoo_tile_13,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex3,oindex1,nv,no,vvoo_tile_32,vvoo_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex2,oindex3,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o3,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex2,oindex3,nv,no,&
                           & ccsdpt_doubles_2(:,:,oindex3),&
                           & wrk_3d,trip,vvvo_tile_o1,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex1,oindex2,no,nv,ovoo_tile_12,ovoo_tile_21,&
                           & ccsdpt_doubles_3,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex1,oindex2,no,nv,ovoo_tile_31,&
                           & ccsdpt_doubles_2,trip,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 4) then ! kji

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,2,1],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,2,1],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex3,oindex2,nv,no,vvoo_tile_32,vvoo_tile_23,&
                       & ccsdpt_singles_1,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex3,oindex2,nv,no,vvoo_tile_31,vvoo_tile_13,&
                       & ccsdpt_singles_2,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex2),&
                           & trip,wrk_3d,vvvo_tile_o3,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,oindex3),&
                           & trip,wrk_3d,vvvo_tile_o2,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex2,oindex1,no,nv,ovoo_tile_21,ovoo_tile_12,&
                           & ccsdpt_doubles_3,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex2,oindex1,no,nv,ovoo_tile_32,&
                           & ccsdpt_doubles_1,wrk_3d,trip,handle,cublas_handle)

       else if (idx .eq. 5) then ! ikj

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex3,oindex1,oindex2,nv,no,vvoo_tile_12,vvoo_tile_21,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex1,oindex2,nv,no,vvoo_tile_13,vvoo_tile_31,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex2),&
                           & wrk_3d,trip,vvvo_tile_o1,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex1),&
                           & wrk_3d,trip,vvvo_tile_o2,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex2,oindex3,no,nv,ovoo_tile_23,ovoo_tile_32,&
                           & ccsdpt_doubles_1,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex2,oindex3,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_3,trip,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 6) then ! jik

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex3,oindex2,oindex1,nv,no,vvoo_tile_21,vvoo_tile_12,&
                       & ccsdpt_singles_3,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex2,oindex1,nv,no,vvoo_tile_23,vvoo_tile_32,&
                       & ccsdpt_singles_1,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex3,oindex2,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex1),&
                           & trip,wrk_3d,vvvo_tile_o2,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex3,oindex2,nv,no,&
                           & ccsdpt_doubles_3(:,:,oindex2),&
                           & trip,wrk_3d,vvvo_tile_o1,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex1,oindex3,no,nv,ovoo_tile_13,ovoo_tile_31,&
                           & ccsdpt_doubles_2,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex1,oindex3,no,nv,ovoo_tile_21,&
                           & ccsdpt_doubles_3,wrk_3d,trip,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case3_ser


#ifdef VAR_MPI
  !> \brief: driver routine for contractions in case(3) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case3_par(vindex1,vindex2,vindex3,no,nv,&
                            & vovv_tile_1,vovv_tile_2,vovv_tile_3,&
                            & oovv_tile_12, oovv_tile_13, oovv_tile_21,&
                            & oovv_tile_23, oovv_tile_31, oovv_tile_32,&
                            & ooov_tile_v1,ooov_tile_v2,ooov_tile_v3,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_2,ccsdpt_doubles_3,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & a,b,c,tile_size_a,tile_size_b,tile_size_c)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,vindex3,no,nv
    integer, intent(in) :: a,b,c,tile_size_a,tile_size_b,tile_size_c
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no,nv) :: ccsdpt_doubles_1,ccsdpt_doubles_2,ccsdpt_doubles_3
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3
    !> tiles of vovv integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout) :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_b), intent(inout) :: vovv_tile_2
    real(realk), dimension(nv,no,nv,tile_size_c), intent(inout) :: vovv_tile_3
    !> tiles of oovv integrals
    real(realk), dimension(no,no,tile_size_a,tile_size_b), intent(inout) :: oovv_tile_12
    real(realk), dimension(no,no,tile_size_a,tile_size_c), intent(inout) :: oovv_tile_13
    real(realk), dimension(no,no,tile_size_b,tile_size_a), intent(inout) :: oovv_tile_21
    real(realk), dimension(no,no,tile_size_b,tile_size_c), intent(inout) :: oovv_tile_23
    real(realk), dimension(no,no,tile_size_c,tile_size_a), intent(inout) :: oovv_tile_31
    real(realk), dimension(no,no,tile_size_c,tile_size_b), intent(inout) :: oovv_tile_32
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex2,vindex3
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v2
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v3
    !> triples amplitude and work array
    real(realk), dimension(nv,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. 
    ! in order to minimize the number of reorderings needed to be performed, 
    ! in order to take optimal advantage of the symmetry of the amplitudes, 
    ! AND to finish all work that involves int_occ_tile_v3 first such that this may be updated (gpu),
    ! we carry out the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : abc --213--> bac --132--> bca --321--> acb --213--> cab --132--> cba
    ! in 211/212 : cab ........ acb ........ abc ........ cba ........ bca ........ bac
    ! in 221/222 : bca ........ cba ........ cab ........ bac ........ abc ........ acb

    do idx = 1,6

       if (idx .eq. 1) then ! abc

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex2,vindex3,nv,no,&
                       & oovv_tile_23(:,:,b,c),oovv_tile_32(:,:,c,b),&
                       & ccsdpt_singles_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex2,vindex3,nv,no,&
                       & oovv_tile_21(:,:,b,a),oovv_tile_12(:,:,a,b),&
                       & ccsdpt_singles_3,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex2,vindex3,vindex1,no,nv,&
                           & vovv_tile_1(:,:,vindex3,a),vovv_tile_3(:,:,vindex1,c),&
                           & ccsdpt_doubles_2,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex3,vindex1,no,nv,vovv_tile_3(:,:,vindex2,c),&
                           & ccsdpt_doubles_1,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex3),&
                           & wrk_3d,trip,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex2),&
                           & wrk_3d,trip,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 2) then ! cab

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[2,1,3],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[2,1,3],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex2,vindex1,vindex3,nv,no,&
                       & oovv_tile_13(:,:,a,c),oovv_tile_31(:,:,c,a),&
                       & ccsdpt_singles_2,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex1,vindex3,nv,no,&
                       & oovv_tile_12(:,:,a,b),oovv_tile_21(:,:,b,a),&
                       & ccsdpt_singles_3,wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex3,vindex2,no,nv,&
                           & vovv_tile_2(:,:,vindex3,b),vovv_tile_3(:,:,vindex2,c),&
                           & ccsdpt_doubles_1,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex3,vindex2,no,nv,vovv_tile_3(:,:,vindex1,c),&
                           & ccsdpt_doubles_2,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex3),&
                           & trip,wrk_3d,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex1),&
                           & trip,wrk_3d,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 3) then ! bca

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[1,3,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[1,3,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex2,vindex3,vindex1,nv,no,&
                       & oovv_tile_31(:,:,c,a),oovv_tile_13(:,:,a,c),&
                       & ccsdpt_singles_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex3,vindex1,nv,no,&
                       & oovv_tile_32(:,:,c,b),oovv_tile_23(:,:,b,c),&
                       & ccsdpt_singles_1,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex3,vindex1,vindex2,no,nv,&
                           & vovv_tile_2(:,:,vindex1,b),vovv_tile_1(:,:,vindex2,a),&
                           & ccsdpt_doubles_3,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex1,vindex2,no,nv,vovv_tile_1(:,:,vindex3,a),&
                           & ccsdpt_doubles_2,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex2,vindex3,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex2,vindex3,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex3),&
                           & wrk_3d,trip,ooov_tile_v1,handle,cublas_handle)

       else if (idx .eq. 4) then ! cba

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[3,2,1],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[3,2,1],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex3,vindex2,nv,no,&
                       & oovv_tile_32(:,:,c,b),oovv_tile_23(:,:,b,c),&
                       & ccsdpt_singles_1,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex3,vindex2,nv,no,&
                       & oovv_tile_31(:,:,c,a),oovv_tile_13(:,:,a,c),&
                       & ccsdpt_singles_2,wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex3,vindex2,vindex1,no,nv,&
                           & vovv_tile_1(:,:,vindex2,a),vovv_tile_2(:,:,vindex1,b),&
                           & ccsdpt_doubles_3,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex2,vindex1,no,nv,vovv_tile_2(:,:,vindex3,b),&
                           & ccsdpt_doubles_1,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex2),&
                           & trip,wrk_3d,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex3),&
                           & trip,wrk_3d,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 5) then ! acb

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,1,3],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,1,3],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex3,vindex1,vindex2,nv,no,&
                       & oovv_tile_12(:,:,a,b),oovv_tile_21(:,:,b,a),&
                       & ccsdpt_singles_3,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex1,vindex2,nv,no,&
                       & oovv_tile_13(:,:,a,c),oovv_tile_31(:,:,c,a),&
                       & ccsdpt_singles_2,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex2,vindex3,no,nv,&
                           & vovv_tile_3(:,:,vindex2,c),vovv_tile_2(:,:,vindex3,b),&
                           & ccsdpt_doubles_1,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex2,vindex3,no,nv,vovv_tile_2(:,:,vindex1,b),&
                           & ccsdpt_doubles_3,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex2),&
                           & wrk_3d,trip,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 6) then ! bac

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex3,vindex2,vindex1,nv,no,&
                       & oovv_tile_21(:,:,b,a),oovv_tile_12(:,:,a,b),&
                       & ccsdpt_singles_3,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex2,vindex1,nv,no,&
                       & oovv_tile_23(:,:,b,c),oovv_tile_32(:,:,c,b),&
                       & ccsdpt_singles_1,wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex2,vindex1,vindex3,no,nv,&
                           & vovv_tile_3(:,:,vindex1,c),vovv_tile_1(:,:,vindex3,a),&
                           & ccsdpt_doubles_2,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex1,vindex3,no,nv,vovv_tile_1(:,:,vindex2,a),&
                           & ccsdpt_doubles_3,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex3,vindex2,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex1),&
                           & trip,wrk_3d,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex3,vindex2,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex2),&
                           & trip,wrk_3d,ooov_tile_v1,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case3_par
#endif


  !> \brief: driver routine for contractions in case(3) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case3_ser(vindex1,vindex2,vindex3,no,nv,oovv,vovv,&
                            & ooov_tile_v1,ooov_tile_v2,ooov_tile_v3,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3,&
                            & ccsdpt_doubles_1,ccsdpt_doubles_2,ccsdpt_doubles_3,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,vindex3,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no,nv) :: ccsdpt_doubles_1,ccsdpt_doubles_2,ccsdpt_doubles_3
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3
    !> vovv integrals
    real(realk), dimension(nv,no,nv,nv), intent(inout)  :: vovv
    !> oovv integrals
    real(realk), dimension(no,no,nv,nv), intent(inout)  :: oovv
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex2,vindex3
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v2
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v3
    !> triples amplitude and work array
    real(realk), dimension(nv,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
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

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. 
    ! in order to minimize the number of reorderings needed to be performed, 
    ! in order to take optimal advantage of the symmetry of the amplitudes, 
    ! AND to finish all work that involves int_occ_tile_v3 first such that this may be updated (gpu),
    ! we carry out the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : abc --213--> bac --132--> bca --321--> acb --213--> cab --132--> cba
    ! in 211/212 : cab ........ acb ........ abc ........ cba ........ bca ........ bac
    ! in 221/222 : bca ........ cba ........ cab ........ bac ........ abc ........ acb

    do idx = 1,6

       if (idx .eq. 1) then ! abc

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex2,vindex3,nv,no,&
                       & oovv(:,:,vindex2,vindex3),oovv(:,:,vindex3,vindex2),&
                       & ccsdpt_singles_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex2,vindex3,nv,no,&
                       & oovv(:,:,vindex2,vindex1),oovv(:,:,vindex1,vindex2),&
                       & ccsdpt_singles_3,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex2,vindex3,vindex1,no,nv,&
                           & vovv(:,:,vindex3,vindex1),vovv(:,:,vindex1,vindex3),&
                           & ccsdpt_doubles_2,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex3,vindex1,no,nv,vovv(:,:,vindex2,vindex3),&
                           & ccsdpt_doubles_1,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex3),&
                           & wrk_3d,trip,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex2),&
                           & wrk_3d,trip,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 2) then ! cab

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[2,1,3],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[2,1,3],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex2,vindex1,vindex3,nv,no,&
                       & oovv(:,:,vindex1,vindex3),oovv(:,:,vindex3,vindex1),&
                       & ccsdpt_singles_2,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex1,vindex3,nv,no,&
                       & oovv(:,:,vindex1,vindex2),oovv(:,:,vindex2,vindex1),&
                       & ccsdpt_singles_3,wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex3,vindex2,no,nv,&
                           & vovv(:,:,vindex3,vindex2),vovv(:,:,vindex2,vindex3),&
                           & ccsdpt_doubles_1,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex3,vindex2,no,nv,vovv(:,:,vindex1,vindex3),&
                           & ccsdpt_doubles_2,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex3),&
                           & trip,wrk_3d,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex1),&
                           & trip,wrk_3d,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 3) then ! bca

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[1,3,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[1,3,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex2,vindex3,vindex1,nv,no,&
                       & oovv(:,:,vindex3,vindex1),oovv(:,:,vindex1,vindex3),&
                       & ccsdpt_singles_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex3,vindex1,nv,no,&
                       & oovv(:,:,vindex3,vindex2),oovv(:,:,vindex2,vindex3),&
                       & ccsdpt_singles_1,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex3,vindex1,vindex2,no,nv,&
                           & vovv(:,:,vindex1,vindex2),vovv(:,:,vindex2,vindex1),&
                           & ccsdpt_doubles_3,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex1,vindex2,no,nv,vovv(:,:,vindex3,vindex1),&
                           & ccsdpt_doubles_2,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex2,vindex3,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex2,vindex3,nv,no,&
                           & ccsdpt_doubles_2(:,:,vindex3),&
                           & wrk_3d,trip,ooov_tile_v1,handle,cublas_handle)

       else if (idx .eq. 4) then ! cba

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[3,2,1],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[3,2,1],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex3,vindex2,nv,no,&
                       & oovv(:,:,vindex3,vindex2),oovv(:,:,vindex2,vindex3),&
                       & ccsdpt_singles_1,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex3,vindex2,nv,no,&
                       & oovv(:,:,vindex3,vindex1),oovv(:,:,vindex1,vindex3),&
                       & ccsdpt_singles_2,wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex3,vindex2,vindex1,no,nv,&
                           & vovv(:,:,vindex2,vindex1),vovv(:,:,vindex1,vindex2),&
                           & ccsdpt_doubles_3,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex2,vindex1,no,nv,vovv(:,:,vindex3,vindex2),&
                           & ccsdpt_doubles_1,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex2),&
                           & trip,wrk_3d,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_1(:,:,vindex3),&
                           & trip,wrk_3d,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 5) then ! acb

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,1,3],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,1,3],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex3,vindex1,vindex2,nv,no,&
                       & oovv(:,:,vindex1,vindex2),oovv(:,:,vindex2,vindex1),&
                       & ccsdpt_singles_3,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex1,vindex2,nv,no,&
                       & oovv(:,:,vindex1,vindex3),oovv(:,:,vindex3,vindex1),&
                       & ccsdpt_singles_2,trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex1,vindex2,vindex3,no,nv,&
                           & vovv(:,:,vindex2,vindex3),vovv(:,:,vindex3,vindex2),&
                           & ccsdpt_doubles_1,trip,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex2,vindex3,no,nv,vovv(:,:,vindex1,vindex2),&
                           & ccsdpt_doubles_3,trip,wrk_3d,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex2),&
                           & wrk_3d,trip,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex1),&
                           & wrk_3d,trip,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 6) then ! bac

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex3,vindex2,vindex1,nv,no,&
                       & oovv(:,:,vindex2,vindex1),oovv(:,:,vindex1,vindex2),&
                       & ccsdpt_singles_3,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex2,vindex1,nv,no,&
                       & oovv(:,:,vindex2,vindex3),oovv(:,:,vindex3,vindex2),&
                       & ccsdpt_singles_1,wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_abc_211(vindex2,vindex1,vindex3,no,nv,&
                           & vovv(:,:,vindex1,vindex3),vovv(:,:,vindex3,vindex1),&
                           & ccsdpt_doubles_2,wrk_3d,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex1,vindex3,no,nv,vovv(:,:,vindex2,vindex1),&
                           & ccsdpt_doubles_3,wrk_3d,trip,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex3,vindex2,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex1),&
                           & trip,wrk_3d,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex3,vindex2,nv,no,&
                           & ccsdpt_doubles_3(:,:,vindex2),&
                           & trip,wrk_3d,ooov_tile_v1,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case3_ser


  !> brief: do the first of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine. int_normal is abij of driver.
  !> nv is nvirt and T_star is ccsdpt_singles of driver. trip_ampl is the triples amplitude array.
  subroutine ccsdpt_contract_ijk_11(oindex1,oindex2,oindex3,nv,no,int_normal_23,int_normal_32,T_star_o1,&
                              & trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv), target :: T_star_o1 ! T_star(:,oinedx1)
    real(realk), dimension(nv,nv), target :: int_normal_23, int_normal_32
    real(realk), dimension(nv,nv,nv), target :: trip_ampl
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

    TypeofContraction_ijk_11: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal and we thus only need to contract with the coulumb part. 
       call ls_dgemm_acc('n','n',nv,1,nv**2,1.0E0_realk,trip_ampl,nv,int_normal_23,nv**2,1.0E0_realk,T_star_o1,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
       call ls_dgemm_acc('n','n',nv,1,nv**2,2.0E0_realk,trip_ampl,nv,int_normal_23,nv**2,1.0E0_realk,T_star_o1,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_dgemm_acc('n','n',nv,1,nv**2,-1.0E0_realk,trip_ampl,nv,int_normal_32,nv**2,1.0E0_realk,T_star_o1,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)

    end select TypeofContraction_ijk_11

  end subroutine ccsdpt_contract_ijk_11


  !> brief: do the first of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_11(vindex1,vindex2,vindex3,nv,no,int_normal_23,int_normal_32,T_star_v1,&
                              & trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(realk), dimension(no), target :: T_star_v1 ! T_star(:,oinedx1)
    real(realk), dimension(no,no), target :: int_normal_23, int_normal_32
    real(realk), dimension(no,no,no), target :: trip_ampl
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

    TypeofContraction_abc_11: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal and we thus only need to contract with the coulumb part. 
       call ls_dgemm_acc('n','n',no,1,no**2,1.0E0_realk,trip_ampl,no,int_normal_23,no**2,1.0E0_realk,T_star_v1,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
       call ls_dgemm_acc('n','n',no,1,no**2,2.0E0_realk,trip_ampl,no,int_normal_23,no**2,1.0E0_realk,T_star_v1,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)
       call ls_dgemm_acc('n','n',no,1,no**2,-1.0E0_realk,trip_ampl,no,int_normal_32,no**2,1.0E0_realk,T_star_v1,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)

    end select TypeofContraction_abc_11

  end subroutine ccsdpt_contract_abc_11


  !> brief: do the second of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine. int_normal is abij of driver.
  !> nv is nvirt and T_star is ccsdpt_singles of driver. trip_ampl is the triples amplitude array.
  subroutine ccsdpt_contract_ijk_12(oindex1,oindex2,oindex3,nv,no,int_normal_21,int_normal_12,T_star_o3,&
                              & trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv), target :: T_star_o3 ! T_star(:,oinedx3)
    real(realk), dimension(nv,nv), target :: int_normal_21, int_normal_12
    real(realk), dimension(nv,nv,nv), target :: trip_ampl
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

    TypeofContraction_ijk_12: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal
       ! and we thus only need to contract with (-1)*coulumb part. 

       ! now contract coulumb term over both indices
       call ls_dgemm_acc('n','n',nv,1,nv**2,-1.0E0_realk,trip_ampl,nv,int_normal_21,nv**2,1.0E0_realk,T_star_o3,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
       call ls_dgemm_acc('n','n',nv,1,nv**2,-2.0E0_realk,trip_ampl,nv,int_normal_21,nv**2,1.0E0_realk,T_star_o3,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_dgemm_acc('n','n',nv,1,nv**2,1.0E0_realk,trip_ampl,nv,int_normal_12,nv**2,1.0E0_realk,T_star_o3,nv,&
                       & int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),int(i8*nv,kind=8),&
                       & async_idx,cublas_handle)

    end select TypeofContraction_ijk_12

  end subroutine ccsdpt_contract_ijk_12


  !> brief: do the second of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_12(vindex1,vindex2,vindex3,nv,no,int_normal_21,int_normal_12,T_star_v3,&
                              & trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(realk), dimension(no), target :: T_star_v3 ! T_star(:,oinedx3)
    real(realk), dimension(no,no), target :: int_normal_21, int_normal_12
    real(realk), dimension(no,no,no), target :: trip_ampl
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

    TypeofContraction_abc_12: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal
       ! and we thus only need to contract with (-1)*coulumb part. 

       ! now contract coulumb term over both indices
       call ls_dgemm_acc('n','n',no,1,no**2,-1.0E0_realk,trip_ampl,no,int_normal_21,no**2,1.0E0_realk,T_star_v3,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
       call ls_dgemm_acc('n','n',no,1,no**2,-2.0E0_realk,trip_ampl,no,int_normal_21,no**2,1.0E0_realk,T_star_v3,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)
       call ls_dgemm_acc('n','n',no,1,no**2,1.0E0_realk,trip_ampl,no,int_normal_12,no**2,1.0E0_realk,T_star_v3,no,&
                       & int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),int(i8*no,kind=8),&
                       & async_idx,cublas_handle)

    end select TypeofContraction_abc_12

  end subroutine ccsdpt_contract_abc_12


  !> brief: do the first of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_1 of driver. trip_ampl is the triples amplitude array.
  !> int_virt_tile is a v^3 tile determined by driver occ index
  !> tmp_g is a 3d work array
  subroutine ccsdpt_contract_ijk_211(oindex1,oindex2,oindex3,nv,no,&
       & T_star_o1o2,tmp_g,trip_ampl,int_virt_tile,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv,nv), target :: T_star_o1o2 ! T_star(:,:,oindex1,oindex2)
    real(realk), dimension(nv,nv,nv), target :: tmp_g,trip_ampl
    real(realk), dimension(nv,nv,nv), target, intent(in) :: int_virt_tile
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type, i,j,k
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

    TypeofContraction_211: select case(contraction_type)

    case(0)

       ! note: here we collect contract over L_{dkbc} and g_{dkbc} in one go.

       ! reorder to obtain coulumb term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)
#endif

       ! now contract coulumb term over 2 first indices
       call ls_dgemm_acc('t','n',nv,nv,nv**2,1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv,&
                       & int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),&
                       & async_idx,cublas_handle)

       ! reorder to obtain exchange term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g)
#endif

       ! now contract exchange term over 2 first indices2
       call ls_dgemm_acc('t','n',nv,nv,nv**2,-1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv,&
                       & int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),&
                       & async_idx,cublas_handle)

    case(1)

       ! note: here we contract over L_{dkbc}.

       ! reorder to obtain coulumb term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)
#endif

       ! now contract coulumb term over 2 first indices
       call ls_dgemm_acc('t','n',nv,nv,nv**2,2.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv,&
                       & int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),&
                       & async_idx,cublas_handle)

       ! reorder to obtain exchange term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g)
#endif

       ! now contract exchange term over 2 first indices
       call ls_dgemm_acc('t','n',nv,nv,nv**2,-1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv,&
                       & int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),&
                       & async_idx,cublas_handle)

    end select TypeofContraction_211

  end subroutine ccsdpt_contract_ijk_211


  !> brief: do the first of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_211(vindex1,vindex2,vindex3,no,nv,&
                               & int_virt_23,int_virt_32,T_star_v1,trip_ampl,wrk_3d,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), dimension(no,no,nv), target :: T_star_v1 ! T_star(:,:,:,vindex1)
    real(realk), dimension(nv,no), target :: int_virt_23, int_virt_32
    real(realk), dimension(nv,no,no), target :: trip_ampl
    real(realk), dimension(nv,no,no), target :: wrk_3d
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
    ! is this a special contraction, i.e., can we handle 221 and 222 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 221 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_abc_211: select case(contraction_type)

    case(0)

       ! now contract coulumb term over first index, then contract exchange term over first index 
       ! for this special case, we only have to subtract one coulumb term
       call ls_dgemm_acc('n','n',nv,no**2,no,1.0E0_realk,int_virt_32,nv,trip_ampl,no,0.0E0_realk,wrk_3d,nv,&
                       & int(i8*nv*no,kind=8),int((i8*no**2)*nv,kind=8),int((i8*no**2)*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_dgemm_acc('n','n',nv,no**2,no,-1.0E0_realk,int_virt_23,nv,trip_ampl,no,1.0E0_realk,wrk_3d,nv,&
                       & int(i8*nv*no,kind=8),int((i8*no**2)*nv,kind=8),int((i8*no**2)*nv,kind=8),&
                       & async_idx,cublas_handle)

    case(1)
 
       ! now contract coulumb term over first index, next contract exchange term over first index
       call ls_dgemm_acc('n','n',nv,no**2,no,2.0E0_realk,int_virt_32,nv,trip_ampl,no,0.0E0_realk,wrk_3d,nv,&
                       & int(i8*nv*no,kind=8),int((i8*no**2)*nv,kind=8),int((i8*no**2)*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_dgemm_acc('n','n',nv,no**2,no,-1.0E0_realk,int_virt_23,nv,trip_ampl,no,1.0E0_realk,wrk_3d,nv,&
                       & int(i8*nv*no,kind=8),int((i8*no**2)*nv,kind=8),int((i8*no**2)*nv,kind=8),&
                       & async_idx,cublas_handle)

    end select TypeofContraction_abc_211

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,no,no,[3,2,1],1.0E0_realk,T_star_v1,async_idx)
#else
    call array_reorder_3d(1.0E0_realk,wrk_3d,nv,no,no,[3,2,1],1.0E0_realk,T_star_v1)
#endif

  end subroutine ccsdpt_contract_abc_211


  !> brief: do the second of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_1 of driver. trip_ampl is the triples amplitude array.
  !> int_virt_tile is a v^3 tile determined by driver occ index
  !> tmp_g is a 3d work array
  subroutine ccsdpt_contract_ijk_212(oindex1,oindex2,oindex3,nv,no,&
       & T_star_o3o2,tmp_g,trip_ampl,int_virt_tile,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv,nv), target :: T_star_o3o2 ! T_star(:,:,oindex3,oindex2)
    real(realk), dimension(nv,nv,nv), target, intent(in) :: int_virt_tile
    real(realk), dimension(nv,nv,nv), target :: tmp_g,trip_ampl
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! reorder to obtain coulumb term 
#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
    call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)
#endif

    ! now contract coulumb term over 2 first indices
    call ls_dgemm_acc('t','n',nv,nv,nv**2,-1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o3o2,nv,&
                    & int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),int(i8*nv**2,kind=8),&
                    & async_idx,cublas_handle)

  end subroutine ccsdpt_contract_ijk_212


  !> brief: do the second of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_212(vindex1,vindex2,vindex3,no,nv,&
                              & int_virt_12,T_star_v3,trip_ampl,wrk_3d,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), dimension(no,no,nv), target :: T_star_v3 ! T_star(:,:,:,vindex3)
    real(realk), dimension(nv,no), target :: int_virt_12
    real(realk), dimension(nv,no,no), target :: trip_ampl
    real(realk), dimension(nv,no,no), target :: wrk_3d
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    call ls_dgemm_acc('n','n',nv,no**2,no,-1.0E0_realk,int_virt_12,nv,trip_ampl,no,0.0E0_realk,wrk_3d,nv,&
                    & int(i8*no*nv,kind=8),int((i8*no**2)*nv,kind=8),int((i8*no**2)*nv,kind=8),&
                    & async_idx,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,no,no,[3,2,1],1.0E0_realk,T_star_v3,async_idx)
#else
    call array_reorder_3d(1.0E0_realk,wrk_3d,nv,no,no,[3,2,1],1.0E0_realk,T_star_v3)
#endif

  end subroutine ccsdpt_contract_abc_212


  !> brief: do the first of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_2 of driver. trip_ampl is the triples amplitud array.
  subroutine ccsdpt_contract_ijk_221(oindex1,oindex2,oindex3,no,nv,&
                               & int_occ_23,int_occ_32,T_star_o1,trip_ampl,wrk_3d,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(nv,nv,no), target :: T_star_o1 ! T_star(:,:,:,oindex1)
    real(realk), dimension(no,nv), target :: int_occ_23, int_occ_32
    real(realk), dimension(nv,nv,nv), target :: trip_ampl,wrk_3d
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
    ! is this a special contraction, i.e., can we handle 221 and 222 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 221 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_221: select case(contraction_type)

    case(0)

       ! now contract coulumb term over first index, then contract exchange term over first index 
       ! for this special case, we only have to subtract one coulumb term
       call ls_dgemm_acc('n','n',no,nv**2,nv,-1.0E0_realk,int_occ_32,no,trip_ampl,nv,0.0E0_realk,wrk_3d,no,&
                       & int(i8*no*nv,kind=8),int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_dgemm_acc('n','n',no,nv**2,nv,1.0E0_realk,int_occ_23,no,trip_ampl,nv,1.0E0_realk,wrk_3d,no,&
                       & int(i8*no*nv,kind=8),int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),&
                       & async_idx,cublas_handle)

    case(1)
 
       ! now contract coulumb term over first index, next contract exchange term over first index
       call ls_dgemm_acc('n','n',no,nv**2,nv,-2.0E0_realk,int_occ_32,no,trip_ampl,nv,0.0E0_realk,wrk_3d,no,&
                       & int(i8*no*nv,kind=8),int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),&
                       & async_idx,cublas_handle)
       call ls_dgemm_acc('n','n',no,nv**2,nv,1.0E0_realk,int_occ_23,no,trip_ampl,nv,1.0E0_realk,wrk_3d,no,&
                       & int(i8*no*nv,kind=8),int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),&
                       & async_idx,cublas_handle)

    end select TypeofContraction_221

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,nv,nv,[3,2,1],1.0E0_realk,T_star_o1,async_idx)
#else
    call array_reorder_3d(1.0E0_realk,wrk_3d,no,nv,nv,[3,2,1],1.0E0_realk,T_star_o1) 
#endif

  end subroutine ccsdpt_contract_ijk_221


  !> brief: do the first of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_221(vindex1,vindex2,vindex3,nv,no,&
       & T_star_v1v2,tmp_g,trip_ampl,int_occ_tile,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(realk), dimension(no,no), target :: T_star_v1v2 ! T_star(:,:,vindex1,vindex2)
    real(realk), dimension(no,no,no), target :: tmp_g,trip_ampl
    real(realk), dimension(no,no,no), target, intent(in) :: int_occ_tile
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

    TypeofContraction_abc_221: select case(contraction_type)

    case(0)

       ! reorder to obtain coulumb term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g)
#endif

       call ls_dgemm_acc('t','n',no,no,no**2,-1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no,&
                       & int((i8*no**2)*no,kind=8),int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),&
                       & async_idx,cublas_handle)

       ! reorder to obtain exchange term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[3,1,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[3,1,2],0.0E0_realk,tmp_g)
#endif

       call ls_dgemm_acc('t','n',no,no,no**2,1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no,&
                       & int((i8*no**2)*no,kind=8),int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),&
                       & async_idx,cublas_handle)

    case(1)

       ! reorder to obtain coulumb term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g)
#endif

       call ls_dgemm_acc('t','n',no,no,no**2,-2.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no,&
                       & int((i8*no**2)*no,kind=8),int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),&
                       & async_idx,cublas_handle)

       ! reorder to obtain exchange term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[3,1,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[3,1,2],0.0E0_realk,tmp_g)
#endif

       call ls_dgemm_acc('t','n',no,no,no**2,1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no,&
                       & int((i8*no**2)*no,kind=8),int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),&
                       & async_idx,cublas_handle)

    end select TypeofContraction_abc_221

  end subroutine ccsdpt_contract_abc_221


  !> brief: do the second of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_2 of driver. trip_ampl is the triples amplitud array.
  subroutine ccsdpt_contract_ijk_222(oindex1,oindex2,oindex3,no,nv,int_occ_12,T_star_o3,trip_ampl,wrk_3d,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(nv,nv,no), target :: T_star_o3 ! T_star(:,:,:,oindex3)
    real(realk), dimension(no,nv), target :: int_occ_12
    real(realk), dimension(nv,nv,nv), target :: trip_ampl,wrk_3d
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! contract coulumb term over first index
    call ls_dgemm_acc('n','n',no,nv**2,nv,1.0E0_realk,int_occ_12,no,trip_ampl,nv,0.0E0_realk,wrk_3d,no,&
                    & int(i8*nv*no,kind=8),int((i8*nv**2)*nv,kind=8),int((i8*nv**2)*nv,kind=8),&
                    & async_idx,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,nv,nv,[3,2,1],1.0E0_realk,T_star_o3,async_idx)
#else
    call array_reorder_3d(1.0E0_realk,wrk_3d,no,nv,nv,[3,2,1],1.0E0_realk,T_star_o3)
#endif

  end subroutine ccsdpt_contract_ijk_222


  !> brief: do the second of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_222(vindex1,vindex2,vindex3,nv,no,&
       & T_star_v3v2,tmp_g,trip_ampl,int_occ_tile,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(realk), dimension(no,no), target :: T_star_v3v2 ! T_star(:,:,vindex3,vindex2)
    real(realk), dimension(no,no,no), target, intent(in) :: int_occ_tile
    real(realk), dimension(no,no,no), target :: tmp_g,trip_ampl
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! reorder
#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
    call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g)
#endif

    call ls_dgemm_acc('t','n',no,no,no**2,1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v3v2,no,&
                    & int((i8*no**2)*no,kind=8),int((i8*no**2)*no,kind=8),int(i8*no**2,kind=8),&
                    & async_idx,cublas_handle)

  end subroutine ccsdpt_contract_abc_222


!endif mod_unreleased
#endif

  subroutine dummy_ccsdpt_dec_routine()

  end subroutine dummy_ccsdpt_dec_routine

end module ccsdpt_dec_module
