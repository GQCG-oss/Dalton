!> @file
!> DEC-CCSD(T) kernels
!> \brief: ccsd(t) module
!> \author: Janus Juul Eriksen
!> \date: 2012-2014, Aarhus
module ccsdpt_tools_module

  use precision
#ifdef VAR_OPENACC
  use openacc
#endif
  use memory_handling
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  use tensor_interface_module
  use cc_tools_module

#ifdef MOD_UNRELEASED
  public :: ptr_init_ijk_pt,ptr_init_abc_pt,ptr_final_ijk_pt,ptr_final_abc_pt,&
          & ptr_init_ijk_par,ptr_init_ijk_ser,ptr_init_abc_par,ptr_init_abc_ser,&
          & ptr_final_ijk_par,ptr_final_ijk_ser,ptr_final_abc_par,ptr_final_abc_ser,&
          & ptr_aliasing_ijk_par,ptr_aliasing_ijk_ser,ptr_aliasing_abc_par,ptr_aliasing_abc_ser
#ifdef VAR_REAL_SP
  public :: sp_ptr_init_ijk_pt,sp_ptr_init_abc_pt,sp_ptr_final_ijk_pt,sp_ptr_final_abc_pt,&
          & sp_ptr_init_ijk_par,sp_ptr_init_ijk_ser,sp_ptr_init_abc_par,sp_ptr_init_abc_ser,&
          & sp_ptr_final_ijk_par,sp_ptr_final_ijk_ser,sp_ptr_final_abc_par,sp_ptr_final_abc_ser,&
          & sp_ptr_aliasing_ijk_par,sp_ptr_aliasing_ijk_ser,sp_ptr_aliasing_abc_par,sp_ptr_aliasing_abc_ser
#endif
  public :: preload_tiles_in_bg_buf,create_comp_array_ccsdpt,job_distrib_ccsdpt  
#else
  public :: dummy_ccsdpt_tools_routine
#endif

  private

contains

#ifdef MOD_UNRELEASED

  subroutine ptr_init_ijk_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), dimension(nvirt,nocc), target :: ccsdpt_singles
    real(realk), dimension(nvirt,nvirt,nocc,nocc), target :: ccsdpt_doubles
    real(realk), pointer, dimension(:,:) :: pt_1
    real(realk), pointer, dimension(:,:,:,:) :: pt_2

    pt_1 => ccsdpt_singles
    pt_2 => ccsdpt_doubles

  end subroutine ptr_init_ijk_pt

  subroutine ptr_init_abc_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), dimension(nocc,nvirt), target :: ccsdpt_singles
    real(realk), dimension(nocc,nocc,nvirt,nvirt), target :: ccsdpt_doubles
    real(realk), pointer, dimension(:,:) :: pt_1
    real(realk), pointer, dimension(:,:,:,:) :: pt_2

    pt_1 => ccsdpt_singles
    pt_2 => ccsdpt_doubles

  end subroutine ptr_init_abc_pt

#ifdef VAR_REAL_SP
  subroutine sp_ptr_init_ijk_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), dimension(nvirt,nocc) :: ccsdpt_singles
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: ccsdpt_doubles
    real(real_sp), pointer, dimension(:,:) :: pt_1
    real(real_sp), pointer, dimension(:,:,:,:) :: pt_2

    call mem_alloc(pt_1,nvirt,nocc)
    call mem_alloc(pt_2,nvirt,nvirt,nocc,nocc)
    pt_1 = real(ccsdpt_singles,kind=4)
    pt_2 = real(ccsdpt_doubles,kind=4)

  end subroutine sp_ptr_init_ijk_pt
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_init_abc_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), dimension(nocc,nvirt) :: ccsdpt_singles
    real(realk), dimension(nocc,nocc,nvirt,nvirt) :: ccsdpt_doubles
    real(real_sp), pointer, dimension(:,:) :: pt_1
    real(real_sp), pointer, dimension(:,:,:,:) :: pt_2

    call mem_alloc(pt_1,nocc,nvirt)
    call mem_alloc(pt_2,nocc,nocc,nvirt,nvirt)
    pt_1 = real(ccsdpt_singles,kind=4)
    pt_2 = real(ccsdpt_doubles,kind=4)

  end subroutine sp_ptr_init_abc_pt
#endif

  subroutine ptr_final_ijk_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), dimension(nvirt,nocc) :: ccsdpt_singles
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: ccsdpt_doubles
    real(realk), pointer, dimension(:,:) :: pt_1
    real(realk), pointer, dimension(:,:,:,:) :: pt_2

    ! empty dummy routine

  end subroutine ptr_final_ijk_pt

  subroutine ptr_final_abc_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), dimension(nocc,nvirt) :: ccsdpt_singles
    real(realk), dimension(nocc,nocc,nvirt,nvirt) :: ccsdpt_doubles
    real(realk), pointer, dimension(:,:) :: pt_1
    real(realk), pointer, dimension(:,:,:,:) :: pt_2

    ! empty dummy routine

  end subroutine ptr_final_abc_pt

#ifdef VAR_REAL_SP
  subroutine sp_ptr_final_ijk_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), dimension(nvirt,nocc) :: ccsdpt_singles
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: ccsdpt_doubles
    real(real_sp), pointer, dimension(:,:) :: pt_1
    real(real_sp), pointer, dimension(:,:,:,:) :: pt_2

    call mem_dealloc(pt_1)
    call mem_dealloc(pt_2)

  end subroutine sp_ptr_final_ijk_pt
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_final_abc_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), dimension(nocc,nvirt) :: ccsdpt_singles
    real(realk), dimension(nocc,nocc,nvirt,nvirt) :: ccsdpt_doubles
    real(real_sp), pointer, dimension(:,:) :: pt_1
    real(real_sp), pointer, dimension(:,:,:,:) :: pt_2

    call mem_dealloc(pt_1)
    call mem_dealloc(pt_2)

  end subroutine sp_ptr_final_abc_pt
#endif

! #########################

  subroutine ptr_init_ijk_par(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), pointer, dimension(:) :: ccsd_i,ccsd_j,ccsd_k
    real(realk), pointer, dimension(:) :: vvvo_i,vvvo_j,vvvo_k
    real(realk), pointer, dimension(:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(realk), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), dimension(nvirt,nocc), target, optional :: t1
    real(realk), pointer, dimension(:,:), optional :: t1_ptr

    if ((present(t1) .and. (.not. present(t1_ptr))) .or. &
          & ((.not. present(t1)) .and. present(t1_ptr))) call lsquit("missing args in ptr_init_ijk_par",-1)
    if (present(t1) .and. present(t1_ptr)) t1_ptr => t1

  end subroutine ptr_init_ijk_par

  subroutine ptr_init_ijk_ser(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), pointer, dimension(:,:,:) :: ccsd_i,ccsd_j,ccsd_k
    real(realk), pointer, dimension(:,:,:) :: vvvo_i,vvvo_j,vvvo_k
    real(realk), pointer, dimension(:,:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(realk), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), dimension(nvirt,nocc), target, optional :: t1
    real(realk), pointer, dimension(:,:), optional :: t1_ptr

    if ((present(t1) .and. (.not. present(t1_ptr))) .or. &
          & ((.not. present(t1)) .and. present(t1_ptr))) call lsquit("missing args in ptr_init_ijk_ser",-1)
    if (present(t1) .and. present(t1_ptr)) t1_ptr => t1

  end subroutine ptr_init_ijk_ser

  subroutine ptr_init_abc_par(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                   & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                   & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), pointer, dimension(:) :: ccsd_a,ccsd_b,ccsd_c
    real(realk), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(realk), pointer, dimension(:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(realk), pointer, dimension(:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), dimension(nocc,nvirt), target, optional :: t1
    real(realk), pointer, dimension(:,:), optional :: t1_ptr

    if ((present(t1) .and. (.not. present(t1_ptr))) .or. &
          & ((.not. present(t1)) .and. present(t1_ptr))) call lsquit("missing args in ptr_init_abc_par",-1)
    if (present(t1) .and. present(t1_ptr)) t1_ptr => t1

  end subroutine ptr_init_abc_par

  subroutine ptr_init_abc_ser(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                   & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                   & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), pointer, dimension(:,:,:) :: ccsd_a,ccsd_b,ccsd_c
    real(realk), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(realk), pointer, dimension(:,:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(realk), pointer, dimension(:,:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), dimension(nocc,nvirt), target, optional :: t1
    real(realk), pointer, dimension(:,:), optional :: t1_ptr

    if ((present(t1) .and. (.not. present(t1_ptr))) .or. &
          & ((.not. present(t1)) .and. present(t1_ptr))) call lsquit("missing args in ptr_init_abc_ser",-1)
    if (present(t1) .and. present(t1_ptr)) t1_ptr => t1

  end subroutine ptr_init_abc_ser

  subroutine ptr_final_ijk_par(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), pointer, dimension(:) :: ccsd_i,ccsd_j,ccsd_k
    real(realk), pointer, dimension(:) :: vvvo_i,vvvo_j,vvvo_k
    real(realk), pointer, dimension(:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(realk), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), pointer, dimension(:,:), optional :: t1_ptr

    ! empty dummy routine

  end subroutine ptr_final_ijk_par

  subroutine ptr_final_ijk_ser(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), pointer, dimension(:,:,:) :: ccsd_i,ccsd_j,ccsd_k
    real(realk), pointer, dimension(:,:,:) :: vvvo_i,vvvo_j,vvvo_k
    real(realk), pointer, dimension(:,:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(realk), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), pointer, dimension(:,:), optional :: t1_ptr

    ! empty dummy routine

  end subroutine ptr_final_ijk_ser

  subroutine ptr_final_abc_par(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                   & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                   & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), pointer, dimension(:) :: ccsd_a,ccsd_b,ccsd_c
    real(realk), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(realk), pointer, dimension(:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(realk), pointer, dimension(:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), pointer, dimension(:,:), optional :: t1_ptr

    ! empty dummy routine

  end subroutine ptr_final_abc_par

  subroutine ptr_final_abc_ser(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                   & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                   & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(realk), pointer, dimension(:,:,:) :: ccsd_a,ccsd_b,ccsd_c
    real(realk), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(realk), pointer, dimension(:,:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(realk), pointer, dimension(:,:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), pointer, dimension(:,:), optional :: t1_ptr

    ! empty dummy routine

  end subroutine ptr_final_abc_ser

#ifdef VAR_REAL_SP
  subroutine sp_ptr_init_ijk_par(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(real_sp), pointer, dimension(:) :: ccsd_i,ccsd_j,ccsd_k
    real(real_sp), pointer, dimension(:) :: vvvo_i,vvvo_j,vvvo_k
    real(real_sp), pointer, dimension(:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(real_sp), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), dimension(nvirt,nocc), target, optional :: t1
    real(real_sp), pointer, dimension(:,:), optional :: t1_ptr

    if ((present(t1) .and. (.not. present(t1_ptr))) .or. &
          & ((.not. present(t1)) .and. present(t1_ptr))) call lsquit("missing args in sp_ptr_init_ijk_par",-1)

    call mem_alloc(t1_ptr,nvirt,nocc)
    t1_ptr = real(t1,kind=4)
    call mem_alloc(ccsd_i,(i8*nocc)*(i8*nvirt**2))
    call mem_alloc(ccsd_j,(i8*nocc)*(i8*nvirt**2))
    call mem_alloc(ccsd_k,(i8*nocc)*(i8*nvirt**2))
    call mem_alloc(vvvo_i,(i8*nvirt)*(i8*nvirt**2))
    call mem_alloc(vvvo_j,(i8*nvirt)*(i8*nvirt**2))
    call mem_alloc(vvvo_k,(i8*nvirt)*(i8*nvirt**2))
    call mem_alloc(vvoo_ij,i8*nvirt**2);call mem_alloc(vvoo_ik,i8*nvirt**2);call mem_alloc(vvoo_ji,i8*nvirt**2)
    call mem_alloc(vvoo_jk,i8*nvirt**2);call mem_alloc(vvoo_ki,i8*nvirt**2);call mem_alloc(vvoo_kj,i8*nvirt**2)
    call mem_alloc(ovoo_ij,nocc,nvirt);call mem_alloc(ovoo_ik,nocc,nvirt);call mem_alloc(ovoo_ji,nocc,nvirt)
    call mem_alloc(ovoo_jk,nocc,nvirt);call mem_alloc(ovoo_ki,nocc,nvirt);call mem_alloc(ovoo_kj,nocc,nvirt)

  end subroutine sp_ptr_init_ijk_par
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_init_ijk_ser(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(real_sp), pointer, dimension(:,:,:) :: ccsd_i,ccsd_j,ccsd_k
    real(real_sp), pointer, dimension(:,:,:) :: vvvo_i,vvvo_j,vvvo_k
    real(real_sp), pointer, dimension(:,:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(real_sp), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), dimension(nvirt,nocc), target, optional :: t1
    real(real_sp), pointer, dimension(:,:), optional :: t1_ptr

    if ((present(t1) .and. (.not. present(t1_ptr))) .or. &
          & ((.not. present(t1)) .and. present(t1_ptr))) call lsquit("missing args in sp_ptr_init_ijk_ser",-1)

    call mem_alloc(t1_ptr,nvirt,nocc)
    t1_ptr = real(t1,kind=4)
    call mem_alloc(ccsd_i,nvirt,nvirt,nocc);call mem_alloc(ccsd_j,nvirt,nvirt,nocc);call mem_alloc(ccsd_k,nvirt,nvirt,nocc)
    call mem_alloc(vvvo_i,nvirt,nvirt,nvirt);call mem_alloc(vvvo_j,nvirt,nvirt,nvirt);call mem_alloc(vvvo_k,nvirt,nvirt,nvirt)
    call mem_alloc(vvoo_ij,nvirt,nvirt);call mem_alloc(vvoo_ik,nvirt,nvirt);call mem_alloc(vvoo_ji,nvirt,nvirt)
    call mem_alloc(vvoo_jk,nvirt,nvirt);call mem_alloc(vvoo_ki,nvirt,nvirt);call mem_alloc(vvoo_kj,nvirt,nvirt)
    call mem_alloc(ovoo_ij,nocc,nvirt);call mem_alloc(ovoo_ik,nocc,nvirt);call mem_alloc(ovoo_ji,nocc,nvirt)
    call mem_alloc(ovoo_jk,nocc,nvirt);call mem_alloc(ovoo_ki,nocc,nvirt);call mem_alloc(ovoo_kj,nocc,nvirt)

  end subroutine sp_ptr_init_ijk_ser
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_init_abc_par(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                   & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                   & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(real_sp), pointer, dimension(:) :: ccsd_a,ccsd_b,ccsd_c
    real(real_sp), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(real_sp), pointer, dimension(:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(real_sp), pointer, dimension(:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), dimension(nocc,nvirt), optional :: t1
    real(real_sp), pointer, dimension(:,:), optional :: t1_ptr

    if ((present(t1) .and. (.not. present(t1_ptr))) .or. &
          & ((.not. present(t1)) .and. present(t1_ptr))) call lsquit("missing args in sp_ptr_init_abc_par",-1)

    call mem_alloc(t1_ptr,nocc,nvirt)
    t1_ptr = real(t1,kind=4)
    call mem_alloc(ccsd_a,(i8*nocc**2)*(i8*nvirt))
    call mem_alloc(ccsd_b,(i8*nocc**2)*(i8*nvirt))
    call mem_alloc(ccsd_c,(i8*nocc**2)*(i8*nvirt))
    call mem_alloc(vovv_ab,i8*nvirt*nocc);call mem_alloc(vovv_ac,i8*nvirt*nocc);call mem_alloc(vovv_ba,i8*nvirt*nocc)
    call mem_alloc(vovv_bc,i8*nvirt*nocc);call mem_alloc(vovv_ca,i8*nvirt*nocc);call mem_alloc(vovv_cb,i8*nvirt*nocc)
    call mem_alloc(oovv_ab,i8*nocc**2);call mem_alloc(oovv_ac,i8*nocc**2);call mem_alloc(oovv_ba,i8*nocc**2)
    call mem_alloc(oovv_bc,i8*nocc**2);call mem_alloc(oovv_ca,i8*nocc**2);call mem_alloc(oovv_cb,i8*nocc**2)
    call mem_alloc(ooov_a,nocc,nocc,nocc);call mem_alloc(ooov_b,nocc,nocc,nocc);call mem_alloc(ooov_c,nocc,nocc,nocc)

  end subroutine sp_ptr_init_abc_par
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_init_abc_ser(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                   & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                   & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(real_sp), pointer, dimension(:,:,:) :: ccsd_a,ccsd_b,ccsd_c
    real(real_sp), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(real_sp), pointer, dimension(:,:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(real_sp), pointer, dimension(:,:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), dimension(nocc,nvirt), optional :: t1
    real(real_sp), pointer, dimension(:,:), optional :: t1_ptr

    if ((present(t1) .and. (.not. present(t1_ptr))) .or. &
          & ((.not. present(t1)) .and. present(t1_ptr))) call lsquit("missing args in sp_ptr_init_abc_ser",-1)

    call mem_alloc(t1_ptr,nocc,nvirt)
    t1_ptr = real(t1,kind=4)
    call mem_alloc(ccsd_a,nocc,nocc,nvirt);call mem_alloc(ccsd_b,nocc,nocc,nvirt);call mem_alloc(ccsd_c,nocc,nocc,nvirt)
    call mem_alloc(ooov_a,nocc,nocc,nocc);call mem_alloc(ooov_b,nocc,nocc,nocc);call mem_alloc(ooov_c,nocc,nocc,nocc)
    call mem_alloc(oovv_ab,nocc,nocc);call mem_alloc(oovv_ac,nocc,nocc);call mem_alloc(oovv_ba,nocc,nocc)
    call mem_alloc(oovv_bc,nocc,nocc);call mem_alloc(oovv_ca,nocc,nocc);call mem_alloc(oovv_cb,nocc,nocc)
    call mem_alloc(vovv_ab,nvirt,nocc);call mem_alloc(vovv_ac,nvirt,nocc);call mem_alloc(vovv_ba,nvirt,nocc)
    call mem_alloc(vovv_bc,nvirt,nocc);call mem_alloc(vovv_ca,nvirt,nocc);call mem_alloc(vovv_cb,nvirt,nocc)

  end subroutine sp_ptr_init_abc_ser
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_final_ijk_par(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(real_sp), pointer, dimension(:) :: ccsd_i,ccsd_j,ccsd_k
    real(real_sp), pointer, dimension(:) :: vvvo_i,vvvo_j,vvvo_k
    real(real_sp), pointer, dimension(:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(real_sp), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(real_sp), pointer, dimension(:,:), optional :: t1_ptr

    if (.not. present(t1_ptr)) call lsquit("missing arg in sp_ptr_final_ijk_par",-1)

    call mem_dealloc(t1_ptr)
    call mem_dealloc(ccsd_i); call mem_dealloc(ccsd_j); call mem_dealloc(ccsd_k)
    call mem_dealloc(vvvo_i); call mem_dealloc(vvvo_j); call mem_dealloc(vvvo_k)
    call mem_dealloc(vvoo_ij); call mem_dealloc(vvoo_ik); call mem_dealloc(vvoo_ji)
    call mem_dealloc(vvoo_jk); call mem_dealloc(vvoo_ki); call mem_dealloc(vvoo_kj)
    call mem_dealloc(ovoo_ij); call mem_dealloc(ovoo_ik); call mem_dealloc(ovoo_ji)
    call mem_dealloc(ovoo_jk); call mem_dealloc(ovoo_ki); call mem_dealloc(ovoo_kj)

  end subroutine sp_ptr_final_ijk_par
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_final_ijk_ser(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(real_sp), pointer, dimension(:,:,:) :: ccsd_i,ccsd_j,ccsd_k
    real(real_sp), pointer, dimension(:,:,:) :: vvvo_i,vvvo_j,vvvo_k
    real(real_sp), pointer, dimension(:,:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(real_sp), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(real_sp), pointer, dimension(:,:), optional :: t1_ptr

    if (.not. present(t1_ptr)) call lsquit("missing arg in sp_ptr_final_ijk_ser",-1)

    call mem_dealloc(t1_ptr)
    call mem_dealloc(ccsd_i); call mem_dealloc(ccsd_j); call mem_dealloc(ccsd_k)
    call mem_dealloc(vvvo_i); call mem_dealloc(vvvo_j); call mem_dealloc(vvvo_k)
    call mem_dealloc(vvoo_ij); call mem_dealloc(vvoo_ik); call mem_dealloc(vvoo_ji)
    call mem_dealloc(vvoo_jk); call mem_dealloc(vvoo_ki); call mem_dealloc(vvoo_kj)
    call mem_dealloc(ovoo_ij); call mem_dealloc(ovoo_ik); call mem_dealloc(ovoo_ji)
    call mem_dealloc(ovoo_jk); call mem_dealloc(ovoo_ki); call mem_dealloc(ovoo_kj)

  end subroutine sp_ptr_final_ijk_ser
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_final_abc_par(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                   & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                   & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(real_sp), pointer, dimension(:) :: ccsd_a,ccsd_b,ccsd_c
    real(real_sp), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(real_sp), pointer, dimension(:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(real_sp), pointer, dimension(:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(real_sp), pointer, dimension(:,:), optional :: t1_ptr

    if (.not. present(t1_ptr)) call lsquit("missing arg in sp_ptr_final_abc_par",-1)

    call mem_dealloc(t1_ptr)
    call mem_dealloc(ccsd_a); call mem_dealloc(ccsd_b); call mem_dealloc(ccsd_c)
    call mem_dealloc(ooov_a); call mem_dealloc(ooov_b); call mem_dealloc(ooov_c)
    call mem_dealloc(vovv_ab); call mem_dealloc(vovv_ac); call mem_dealloc(vovv_ba)
    call mem_dealloc(vovv_bc); call mem_dealloc(vovv_ca); call mem_dealloc(vovv_cb)
    call mem_dealloc(oovv_ab); call mem_dealloc(oovv_ac); call mem_dealloc(oovv_ba)
    call mem_dealloc(oovv_bc); call mem_dealloc(oovv_ca); call mem_dealloc(oovv_cb)

  end subroutine sp_ptr_final_abc_par
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_final_abc_ser(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                   & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                   & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1_ptr)

    implicit none

    !> pointers
    integer :: nvirt,nocc
    real(real_sp), pointer, dimension(:,:,:) :: ccsd_a,ccsd_b,ccsd_c
    real(real_sp), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(real_sp), pointer, dimension(:,:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(real_sp), pointer, dimension(:,:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(real_sp), pointer, dimension(:,:), optional :: t1_ptr

    if (.not. present(t1_ptr)) call lsquit("missing arg in sp_ptr_final_abc_ser",-1)

    call mem_dealloc(t1_ptr)
    call mem_dealloc(ccsd_a); call mem_dealloc(ccsd_b); call mem_dealloc(ccsd_c)
    call mem_dealloc(ooov_a); call mem_dealloc(ooov_b); call mem_dealloc(ooov_c)
    call mem_dealloc(vovv_ab); call mem_dealloc(vovv_ac); call mem_dealloc(vovv_ba)
    call mem_dealloc(vovv_bc); call mem_dealloc(vovv_ca); call mem_dealloc(vovv_cb)
    call mem_dealloc(oovv_ab); call mem_dealloc(oovv_ac); call mem_dealloc(oovv_ba)
    call mem_dealloc(oovv_bc); call mem_dealloc(oovv_ca); call mem_dealloc(oovv_cb)

  end subroutine sp_ptr_final_abc_ser
#endif

  subroutine ptr_aliasing_ijk_par(nvirt,nocc,i,j,k,i_count,j_count,k_count,ts_i,ts_j,ts_k,&
                   & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                   & ccsd_pdm_i,ccsd_pdm_j,ccsd_pdm_k,&
                   & vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k,&
                   & vvoo_pdm_ij,vvoo_pdm_ik,vvoo_pdm_ji,vvoo_pdm_jk,vvoo_pdm_ki,vvoo_pdm_kj,&
                   & ovoo,async_id,num_ids,case_num)

    implicit none

    !> pointers
    integer :: nvirt,nocc,i,j,k,i_count,j_count,k_count,ts_i,ts_j,ts_k
    real(realk), pointer, dimension(:) :: ccsd_i,ccsd_j,ccsd_k
    real(realk), pointer, dimension(:) :: vvvo_i,vvvo_j,vvvo_k
    real(realk), pointer, dimension(:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(realk), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), pointer, dimension(:) :: ccsd_pdm_i,ccsd_pdm_j,ccsd_pdm_k
    real(realk), pointer, dimension(:) :: vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k
    real(realk), pointer, dimension(:) :: vvoo_pdm_ij,vvoo_pdm_ik,vvoo_pdm_ji,vvoo_pdm_jk,vvoo_pdm_ki,vvoo_pdm_kj
    real(realk), dimension(nocc,nvirt,nocc,nocc), target :: ovoo
    integer :: num_ids
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_id(num_ids)
#else
    integer :: async_id(num_ids)
#endif
    integer :: case_num
    integer(kind=long) :: start_ccsd_i,start_ccsd_j,start_ccsd_k,end_ccsd_i,end_ccsd_j,end_ccsd_k
    integer(kind=long) :: start_vvvo_i,start_vvvo_j,start_vvvo_k,end_vvvo_i,end_vvvo_j,end_vvvo_k
    integer(kind=long) :: start_ij,start_ji,start_ik,start_ki,start_jk,start_kj
    integer(kind=long) :: end_ij,end_ji,end_ik,end_ki,end_jk,end_kj

    select case(case_num)

    case(1)

       start_ccsd_i = int(i8*i_count*(i8*nocc)*(i8*nvirt**2),kind=8)
       start_vvvo_i = int(i8*i_count*(i8*nvirt)*(i8*nvirt**2),kind=8)
       end_ccsd_i = start_ccsd_i + int(i8*nocc*(i8*nvirt**2),kind=8)
       end_vvvo_i = start_vvvo_i + int(i8*nvirt*(i8*nvirt**2),kind=8)

       ccsd_i => ccsd_pdm_i(start_ccsd_i+1:end_ccsd_i)
       vvvo_i => vvvo_pdm_i(start_vvvo_i+1:end_vvvo_i)

!$acc enter data copyin(ccsd_i,vvvo_i) async(async_id(1))

    case(2)

       start_ccsd_j = int(i8*j_count*(i8*nocc)*(i8*nvirt**2),kind=8)
       start_vvvo_j = int(i8*j_count*(i8*nvirt)*(i8*nvirt**2),kind=8)
       end_ccsd_j = start_ccsd_j + int(i8*nocc*(i8*nvirt**2),kind=8)
       end_vvvo_j = start_vvvo_j + int(i8*nvirt*(i8*nvirt**2),kind=8)
       start_ij = int((i8*((i_count-1)+j_count*ts_i))*(i8*nvirt**2),kind=8)
       start_ji = int((i8*(j_count+(i_count-1)*ts_j))*(i8*nvirt**2),kind=8)
       end_ij = start_ij + int(i8*nvirt**2,kind=8)
       end_ji = start_ji + int(i8*nvirt**2,kind=8)

       if (j .eq. i) then

          ovoo_ij => ovoo(:,:,i,j)
          vvoo_ij => vvoo_pdm_ij(start_ij+1:end_ij)

!$acc enter data copyin(ovoo_ij) async(async_id(1))
!$acc enter data copyin(vvoo_ij) async(async_id(3))

       else if (j .lt. i) then

          ccsd_j => ccsd_pdm_j(start_ccsd_j+1:end_ccsd_j)
          vvvo_j => vvvo_pdm_j(start_vvvo_j+1:end_vvvo_j)
          ovoo_ij => ovoo(:,:,i,j); ovoo_ji => ovoo(:,:,j,i)
          vvoo_ij => vvoo_pdm_ij(start_ij+1:end_ij); vvoo_ji => vvoo_pdm_ji(start_ji+1:end_ji)

!$acc enter data copyin(ccsd_j,vvvo_j,ovoo_ij,ovoo_ji) async(async_id(1))
!$acc enter data copyin(vvoo_ij,vvoo_ji) async(async_id(3))

       endif

    case(3)

       start_ccsd_k = int(i8*k_count*(i8*nocc)*(i8*nvirt**2),kind=8)
       start_vvvo_k = int(i8*k_count*(i8*nvirt)*(i8*nvirt**2),kind=8)
       end_ccsd_k = start_ccsd_k + int(i8*nocc*(i8*nvirt**2),kind=8)
       end_vvvo_k = start_vvvo_k + int(i8*nvirt*(i8*nvirt**2),kind=8)
       start_ik = int((i8*((i_count-1)+k_count*ts_i))*(i8*nvirt**2),kind=8)
       start_ki = int((i8*(k_count+(i_count-1)*ts_k))*(i8*nvirt**2),kind=8)
       start_jk = int((i8*((j_count-1)+k_count*ts_j))*(i8*nvirt**2),kind=8)
       start_kj = int((i8*(k_count+(j_count-1)*ts_k))*(i8*nvirt**2),kind=8)
       end_ik = start_ik + int(i8*nvirt**2,kind=8)
       end_ki = start_ki + int(i8*nvirt**2,kind=8)
       end_jk = start_jk + int(i8*nvirt**2,kind=8)
       end_kj = start_kj + int(i8*nvirt**2,kind=8)
 
       if ((i .eq. j) .and. (j .gt. k)) then

          ccsd_k => ccsd_pdm_k(start_ccsd_k+1:end_ccsd_k)
          vvvo_k => vvvo_pdm_k(start_vvvo_k+1:end_vvvo_k)
          ovoo_ik => ovoo(:,:,i,k); ovoo_ki => ovoo(:,:,k,i)
          vvoo_ik => vvoo_pdm_ik(start_ik+1:end_ik); vvoo_ki => vvoo_pdm_ki(start_ki+1:end_ki)

!$acc enter data copyin(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki) async(async_id(1))
!$acc enter data copyin(vvoo_ik,vvoo_ki) async(async_id(3))

       else if ((i .gt. j) .and. (j .eq. k)) then

          ovoo_jk => ovoo(:,:,j,k)
          vvoo_jk => vvoo_pdm_jk(start_jk+1:end_jk)

!$acc enter data copyin(ovoo_jk) async(async_id(1))
!$acc enter data copyin(vvoo_jk) async(async_id(3))

       else if ((i .gt. j) .and. (j .gt. k)) then

          ccsd_k => ccsd_pdm_k(start_ccsd_k+1:end_ccsd_k)
          vvvo_k => vvvo_pdm_k(start_vvvo_k+1:end_vvvo_k)
          ovoo_ik => ovoo(:,:,i,k); ovoo_ki => ovoo(:,:,k,i)
          ovoo_jk => ovoo(:,:,j,k); ovoo_kj => ovoo(:,:,k,j)
          vvoo_ik => vvoo_pdm_ik(start_ik+1:end_ik); vvoo_ki => vvoo_pdm_ki(start_ki+1:end_ki)
          vvoo_jk => vvoo_pdm_jk(start_jk+1:end_jk); vvoo_kj => vvoo_pdm_kj(start_kj+1:end_kj)

!$acc enter data copyin(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki,ovoo_jk,ovoo_kj) async(async_id(1))
!$acc enter data copyin(vvoo_ik,vvoo_ki,vvoo_jk,vvoo_kj) async(async_id(3))

       end if

    end select

  end subroutine ptr_aliasing_ijk_par

  subroutine ptr_aliasing_ijk_ser(nvirt,nocc,i,j,k,&
                   & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                   & ccsd_doubles,vvvo,vvoo,ovoo,async_id,num_ids,case_num)

    implicit none

    !> pointers
    integer :: nvirt,nocc,i,j,k
    real(realk), pointer, dimension(:,:,:) :: ccsd_i,ccsd_j,ccsd_k
    real(realk), pointer, dimension(:,:,:) :: vvvo_i,vvvo_j,vvvo_k
    real(realk), pointer, dimension(:,:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(realk), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), dimension(nvirt,nvirt,nocc,nocc), target :: ccsd_doubles
    real(realk), dimension(nvirt,nvirt,nvirt,nocc), target :: vvvo
    real(realk), dimension(nvirt,nvirt,nocc,nocc), target :: vvoo
    real(realk), dimension(nocc,nvirt,nocc,nocc), target :: ovoo
    integer :: num_ids
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_id(num_ids)
#else
    integer :: async_id(num_ids)
#endif
    integer :: case_num

    select case(case_num)

    case(1)

       ccsd_i => ccsd_doubles(:,:,:,i)
       vvvo_i => vvvo(:,:,:,i)

!$acc enter data copyin(ccsd_i,vvvo_i) async(async_id(1))

    case(2)

       if (j .eq. i) then

          ovoo_ij => ovoo(:,:,i,j)
          vvoo_ij => vvoo(:,:,i,j)

!$acc enter data copyin(ovoo_ij) async(async_id(1))
!$acc enter data copyin(vvoo_ij) async(async_id(3))

       else if (j .lt. i) then

          ccsd_j => ccsd_doubles(:,:,:,j)
          vvvo_j => vvvo(:,:,:,j)
          ovoo_ij => ovoo(:,:,i,j); ovoo_ji => ovoo(:,:,j,i)
          vvoo_ij => vvoo(:,:,i,j); vvoo_ji => vvoo(:,:,j,i)

!$acc enter data copyin(ccsd_j,vvvo_j,ovoo_ij,ovoo_ji) async(async_id(1))
!$acc enter data copyin(vvoo_ij,vvoo_ji) async(async_id(3))

       end if

    case(3)

       if ((i .eq. j) .and. (j .gt. k)) then

          ccsd_k => ccsd_doubles(:,:,:,k)
          vvvo_k => vvvo(:,:,:,k)
          ovoo_ik => ovoo(:,:,i,k); ovoo_ki => ovoo(:,:,k,i)
          vvoo_ik => vvoo(:,:,i,k); vvoo_ki => vvoo(:,:,k,i)

!$acc enter data copyin(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki) async(async_id(1))
!$acc enter data copyin(vvoo_ik,vvoo_ki) async(async_id(3))

       else if ((i .gt. j) .and. (j .eq. k)) then

          ovoo_jk => ovoo(:,:,j,k)
          vvoo_jk => vvoo(:,:,j,k)

!$acc enter data copyin(ovoo_jk) async(async_id(1))
!$acc enter data copyin(vvoo_jk) async(async_id(3))

       else if ((i .gt. j) .and. (j .gt. k)) then

          ccsd_k => ccsd_doubles(:,:,:,k)
          vvvo_k => vvvo(:,:,:,k)
          ovoo_ik => ovoo(:,:,i,k); ovoo_ki => ovoo(:,:,k,i)
          ovoo_jk => ovoo(:,:,j,k); ovoo_kj => ovoo(:,:,k,j)
          vvoo_ik => vvoo(:,:,i,k); vvoo_ki => vvoo(:,:,k,i)
          vvoo_jk => vvoo(:,:,j,k); vvoo_kj => vvoo(:,:,k,j)

!$acc enter data copyin(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki,ovoo_jk,ovoo_kj) async(async_id(1))
!$acc enter data copyin(vvoo_ik,vvoo_ki,vvoo_jk,vvoo_kj) async(async_id(3))

       end if

    end select

  end subroutine ptr_aliasing_ijk_ser

  subroutine ptr_aliasing_abc_par(nvirt,nocc,a,b,c,a_count,b_count,c_count,ts_a,ts_b,ts_c,&
                     & ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                     & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                     & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                     & ccsd_pdm_a,ccsd_pdm_b,ccsd_pdm_c,ooov,&
                     & vovv_pdm_a,vovv_pdm_b,vovv_pdm_c,&
                     & oovv_pdm_ab,oovv_pdm_ba,oovv_pdm_ac,oovv_pdm_ca,oovv_pdm_bc,oovv_pdm_cb,async_id,num_ids,case_num)

    implicit none

    integer :: nvirt,nocc,a,b,c,a_count,b_count,c_count,ts_a,ts_b,ts_c,case_num
    real(realk), pointer, dimension(:) :: ccsd_a,ccsd_b,ccsd_c
    real(realk), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(realk), pointer, dimension(:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(realk), pointer, dimension(:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), pointer, dimension(:) :: ccsd_pdm_a,ccsd_pdm_b,ccsd_pdm_c
    real(realk), dimension(nocc,nocc,nocc,nvirt), target :: ooov
    real(realk), pointer, dimension(:) :: vovv_pdm_a,vovv_pdm_b,vovv_pdm_c
    real(realk), pointer, dimension(:) :: oovv_pdm_ab,oovv_pdm_ba,oovv_pdm_ac,oovv_pdm_ca,oovv_pdm_bc,oovv_pdm_cb
    integer :: num_ids
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_id(num_ids)
#else
    integer :: async_id(num_ids)
#endif
    integer(kind=long) :: start_ccsd_a,start_ccsd_b,start_ccsd_c,end_ccsd_a,end_ccsd_b,end_ccsd_c
    integer(kind=long) :: start_vovv_ab,start_vovv_ba,start_vovv_ac,start_vovv_ca,start_vovv_bc,start_vovv_cb
    integer(kind=long) :: start_oovv_ab,start_oovv_ba,start_oovv_ac,start_oovv_ca,start_oovv_bc,start_oovv_cb
    integer(kind=long) :: end_vovv_ab,end_vovv_ba,end_vovv_ac,end_vovv_ca,end_vovv_bc,end_vovv_cb
    integer(kind=long) :: end_oovv_ab,end_oovv_ba,end_oovv_ac,end_oovv_ca,end_oovv_bc,end_oovv_cb

    select case(case_num)

    case(1)

       start_ccsd_a = int(i8*a_count*(i8*nvirt)*(i8*nocc**2),kind=8)
       end_ccsd_a = start_ccsd_a + int(i8*nvirt*(i8*nocc**2),kind=8)
   
       ccsd_a => ccsd_pdm_a(start_ccsd_a+1:end_ccsd_a)
       ooov_a => ooov(:,:,:,a)

!$acc enter data copyin(ccsd_a,ooov_a) async(async_id(1))

    case(2)

       start_ccsd_b = int(i8*b_count*(i8*nvirt)*(i8*nocc**2),kind=8)
       end_ccsd_b = start_ccsd_b + int(i8*nvirt*(i8*nocc**2),kind=8)
       start_vovv_ab = int((i8*((a-1)+b_count*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_vovv_ba = int((i8*((b-1)+(a_count-1)*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_oovv_ab = int((i8*((a_count-1)+b_count*ts_a))*(i8*nocc**2),kind=8)
       start_oovv_ba = int((i8*(b_count+(a_count-1)*ts_b))*(i8*nocc**2),kind=8)
       end_vovv_ab = start_vovv_ab + int((i8*nvirt)*(i8*nocc),kind=8)
       end_vovv_ba = start_vovv_ba + int((i8*nvirt)*(i8*nocc),kind=8)
       end_oovv_ab = start_oovv_ab + int(i8*nocc**2,kind=8)
       end_oovv_ba = start_oovv_ba + int(i8*nocc**2,kind=8)

       if (b .eq. a) then

          vovv_ab => vovv_pdm_b(start_vovv_ab+1:end_vovv_ab)
          oovv_ab => oovv_pdm_ab(start_oovv_ab+1:end_oovv_ab)

!$acc enter data copyin(vovv_ab) async(async_id(1))
!$acc enter data copyin(oovv_ab) async(async_id(3))

       else if (b .lt. a) then

          ccsd_b => ccsd_pdm_b(start_ccsd_b+1:end_ccsd_b)
          ooov_b => ooov(:,:,:,b)
          vovv_ab => vovv_pdm_b(start_vovv_ab+1:end_vovv_ab)
          vovv_ba => vovv_pdm_a(start_vovv_ba+1:end_vovv_ba)
          oovv_ab => oovv_pdm_ab(start_oovv_ab+1:end_oovv_ab)
          oovv_ba => oovv_pdm_ba(start_oovv_ba+1:end_oovv_ba)

!$acc enter data copyin(ccsd_b,ooov_b,vovv_ab,vovv_ba) async(async_id(1))
!$acc enter data copyin(oovv_ab,oovv_ba) async(async_id(3))

       endif

    case(3)

       start_ccsd_c = int(i8*c_count*(i8*nvirt)*(i8*nocc**2),kind=8)
       end_ccsd_c = start_ccsd_c + int(i8*nvirt*(i8*nocc**2),kind=8)
       start_vovv_ac = int((i8*((a-1)+c_count*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_vovv_ca = int((i8*((c-1)+(a_count-1)*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_vovv_bc = int((i8*((b-1)+c_count*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_vovv_cb = int((i8*((c-1)+(b_count-1)*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_oovv_ac = int((i8*((a_count-1)+c_count*ts_a))*(i8*nocc**2),kind=8)
       start_oovv_ca = int((i8*(c_count+(a_count-1)*ts_c))*(i8*nocc**2),kind=8)
       start_oovv_bc = int((i8*((b_count-1)+c_count*ts_b))*(i8*nocc**2),kind=8)
       start_oovv_cb = int((i8*(c_count+(b_count-1)*ts_c))*(i8*nocc**2),kind=8)
       end_vovv_ac = start_vovv_ac + int((i8*nvirt)*(i8*nocc),kind=8)
       end_vovv_ca = start_vovv_ca + int((i8*nvirt)*(i8*nocc),kind=8)
       end_vovv_bc = start_vovv_bc + int((i8*nvirt)*(i8*nocc),kind=8)
       end_vovv_cb = start_vovv_cb + int((i8*nvirt)*(i8*nocc),kind=8)
       end_oovv_ac = start_oovv_ac + int(i8*nocc**2,kind=8)
       end_oovv_ca = start_oovv_ca + int(i8*nocc**2,kind=8)
       end_oovv_bc = start_oovv_bc + int(i8*nocc**2,kind=8)
       end_oovv_cb = start_oovv_cb + int(i8*nocc**2,kind=8)

       if ((a .eq. b) .and. (b .gt. c)) then
 
          ccsd_c => ccsd_pdm_c(start_ccsd_c+1:end_ccsd_c)
          ooov_c => ooov(:,:,:,c)
          vovv_ac => vovv_pdm_c(start_vovv_ac+1:end_vovv_ac)
          vovv_ca => vovv_pdm_a(start_vovv_ca+1:end_vovv_ca)
          oovv_ac => oovv_pdm_ac(start_oovv_ac+1:end_oovv_ac)
          oovv_ca => oovv_pdm_ca(start_oovv_ca+1:end_oovv_ca)

!$acc enter data copyin(ccsd_c,ooov_c,vovv_ac,vovv_ca) async(async_id(1))
!$acc enter data copyin(oovv_ac,oovv_ca) async(async_id(3))
 
       else if ((a .gt. b) .and. (b .eq. c)) then
 
          vovv_bc => vovv_pdm_c(start_vovv_bc+1:end_vovv_bc)
          oovv_bc => oovv_pdm_bc(start_oovv_bc+1:end_oovv_bc)
 
!$acc enter data copyin(vovv_bc) async(async_id(1))
!$acc enter data copyin(oovv_bc) async(async_id(3))

       else if ((a .gt. b) .and. (b .gt. c)) then
   
          ccsd_c => ccsd_pdm_c(start_ccsd_c+1:end_ccsd_c)
          ooov_c => ooov(:,:,:,c)
          vovv_ac => vovv_pdm_c(start_vovv_ac+1:end_vovv_ac)
          vovv_ca => vovv_pdm_a(start_vovv_ca+1:end_vovv_ca)
          vovv_bc => vovv_pdm_c(start_vovv_bc+1:end_vovv_bc)
          vovv_cb => vovv_pdm_b(start_vovv_cb+1:end_vovv_cb)
          oovv_ac => oovv_pdm_ac(start_oovv_ac+1:end_oovv_ac)
          oovv_ca => oovv_pdm_ca(start_oovv_ca+1:end_oovv_ca)
          oovv_bc => oovv_pdm_bc(start_oovv_bc+1:end_oovv_bc)
          oovv_cb => oovv_pdm_cb(start_oovv_cb+1:end_oovv_cb)

!$acc enter data copyin(ccsd_c,ooov_c,vovv_ac,vovv_ca,vovv_bc,vovv_cb) async(async_id(1))
!$acc enter data copyin(oovv_ac,oovv_ca,oovv_bc,oovv_cb) async(async_id(3))
   
       end if

    end select

  end subroutine ptr_aliasing_abc_par

  subroutine ptr_aliasing_abc_ser(nvirt,nocc,a,b,c,&
                     & ccsd_a,ccsd_b,ccsd_c,&
                     & ooov_a,ooov_b,ooov_c,&
                     & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                     & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                     & ccsd_doubles,ooov,vovv,oovv,async_id,num_ids,case_num)

    implicit none

    integer :: nvirt,nocc,a,b,c
    real(realk), pointer, dimension(:,:,:) :: ccsd_a,ccsd_b,ccsd_c
    real(realk), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(realk), pointer, dimension(:,:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(realk), pointer, dimension(:,:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), dimension(nocc,nocc,nvirt,nvirt), target, intent(inout) :: ccsd_doubles
    real(realk), dimension(nocc,nocc,nocc,nvirt), target, intent(inout) :: ooov
    real(realk), dimension(nvirt,nocc,nvirt,nvirt), target, intent(inout) :: vovv
    real(realk), dimension(nocc,nocc,nvirt,nvirt), target, intent(inout) :: oovv
    integer :: num_ids
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_id(num_ids)
#else
    integer :: async_id(num_ids)
#endif
    integer :: case_num

    select case(case_num)

    case(1)

       ccsd_a => ccsd_doubles(:,:,:,a)
       ooov_a => ooov(:,:,:,a)

!$acc enter data copyin(ccsd_a,ooov_a) async(async_id(1))

    case(2)

       if (b .eq. a) then
  
          vovv_ab => vovv(:,:,a,b)
          oovv_ab => oovv(:,:,a,b)
  
!$acc enter data copyin(vovv_ab) async(async_id(1))
!$acc enter data copyin(oovv_ab) async(async_id(3))

       else if (b .lt. a) then
  
          ccsd_b => ccsd_doubles(:,:,:,b)
          ooov_b => ooov(:,:,:,b)
          vovv_ab => vovv(:,:,a,b); vovv_ba => vovv(:,:,b,a)
          oovv_ab => oovv(:,:,a,b); oovv_ba => oovv(:,:,b,a)
  
!$acc enter data copyin(ccsd_b,ooov_b,vovv_ab,vovv_ba) async(async_id(1))
!$acc enter data copyin(oovv_ab,oovv_ba) async(async_id(3))

       endif

    case(3)

       if ((a .eq. b) .and. (b .gt. c)) then
 
          ccsd_c => ccsd_doubles(:,:,:,c)
          ooov_c => ooov(:,:,:,c)
          vovv_ac => vovv(:,:,a,c); vovv_ca => vovv(:,:,c,a)
          oovv_ac => oovv(:,:,a,c); oovv_ca => oovv(:,:,c,a)
 
!$acc enter data copyin(ccsd_c,ooov_c,vovv_ac,vovv_ca) async(async_id(1))
!$acc enter data copyin(oovv_ac,oovv_ca) async(async_id(3))

       else if ((a .gt. b) .and. (b .eq. c)) then
 
          vovv_bc => vovv(:,:,b,c)
          oovv_bc => oovv(:,:,b,c)
 
!$acc enter data copyin(vovv_bc) async(async_id(1))
!$acc enter data copyin(oovv_bc) async(async_id(3))

       else if ((a .gt. b) .and. (b .gt. c)) then
 
          ccsd_c => ccsd_doubles(:,:,:,c)
          ooov_c => ooov(:,:,:,c)
          vovv_ac => vovv(:,:,a,c); vovv_ca => vovv(:,:,c,a)
          vovv_bc => vovv(:,:,b,c); vovv_cb => vovv(:,:,c,b)
          oovv_ac => oovv(:,:,a,c); oovv_ca => oovv(:,:,c,a)
          oovv_bc => oovv(:,:,b,c); oovv_cb => oovv(:,:,c,b)
 
!$acc enter data copyin(ccsd_c,ooov_c,vovv_ac,vovv_ca,vovv_bc,vovv_cb) async(async_id(1))
!$acc enter data copyin(oovv_ac,oovv_ca,oovv_bc,oovv_cb) async(async_id(3))

       end if

    end select

  end subroutine ptr_aliasing_abc_ser

#ifdef VAR_REAL_SP
  subroutine sp_ptr_aliasing_ijk_par(nvirt,nocc,i,j,k,i_count,j_count,k_count,ts_i,ts_j,ts_k,&
                   & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                   & ccsd_pdm_i,ccsd_pdm_j,ccsd_pdm_k,&
                   & vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k,&
                   & vvoo_pdm_ij,vvoo_pdm_ik,vvoo_pdm_ji,vvoo_pdm_jk,vvoo_pdm_ki,vvoo_pdm_kj,&
                   & ovoo,async_id,num_ids,case_num)

    implicit none

    !> pointers
    integer :: nvirt,nocc,i,j,k,i_count,j_count,k_count,ts_i,ts_j,ts_k
    real(real_sp), pointer, dimension(:) :: ccsd_i,ccsd_j,ccsd_k
    real(real_sp), pointer, dimension(:) :: vvvo_i,vvvo_j,vvvo_k
    real(real_sp), pointer, dimension(:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(real_sp), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), pointer, dimension(:) :: ccsd_pdm_i,ccsd_pdm_j,ccsd_pdm_k
    real(realk), pointer, dimension(:) :: vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k
    real(realk), pointer, dimension(:) :: vvoo_pdm_ij,vvoo_pdm_ik,vvoo_pdm_ji,vvoo_pdm_jk,vvoo_pdm_ki,vvoo_pdm_kj
    real(realk), dimension(nocc,nvirt,nocc,nocc), target :: ovoo
    integer :: num_ids
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_id(num_ids)
#else
    integer :: async_id(num_ids)
#endif
    integer :: case_num
    integer(kind=long) :: start_ccsd_i,start_ccsd_j,start_ccsd_k,end_ccsd_i,end_ccsd_j,end_ccsd_k
    integer(kind=long) :: start_vvvo_i,start_vvvo_j,start_vvvo_k,end_vvvo_i,end_vvvo_j,end_vvvo_k
    integer(kind=long) :: start_ij,start_ji,start_ik,start_ki,start_jk,start_kj
    integer(kind=long) :: end_ij,end_ji,end_ik,end_ki,end_jk,end_kj

    select case(case_num)

    case(1)

       start_ccsd_i = int(i8*i_count*(i8*nocc)*(i8*nvirt**2),kind=8)
       start_vvvo_i = int(i8*i_count*(i8*nvirt)*(i8*nvirt**2),kind=8)
       end_ccsd_i = start_ccsd_i + int(i8*nocc*(i8*nvirt**2),kind=8)
       end_vvvo_i = start_vvvo_i + int(i8*nvirt*(i8*nvirt**2),kind=8)

       ccsd_i = real(ccsd_pdm_i(start_ccsd_i+1:end_ccsd_i),kind=4)
       vvvo_i = real(vvvo_pdm_i(start_vvvo_i+1:end_vvvo_i),kind=4)

!$acc enter data copyin(ccsd_i,vvvo_i) async(async_id(1))

    case(2)

       start_ccsd_j = int(i8*j_count*(i8*nocc)*(i8*nvirt**2),kind=8)
       start_vvvo_j = int(i8*j_count*(i8*nvirt)*(i8*nvirt**2),kind=8)
       end_ccsd_j = start_ccsd_j + int(i8*nocc*(i8*nvirt**2),kind=8)
       end_vvvo_j = start_vvvo_j + int(i8*nvirt*(i8*nvirt**2),kind=8)
       start_ij = int((i8*((i_count-1)+j_count*ts_i))*(i8*nvirt**2),kind=8)
       start_ji = int((i8*(j_count+(i_count-1)*ts_j))*(i8*nvirt**2),kind=8)
       end_ij = start_ij + int(i8*nvirt**2,kind=8)
       end_ji = start_ji + int(i8*nvirt**2,kind=8)

       if (j .eq. i) then

          ovoo_ij = real(ovoo(:,:,i,j),kind=4)
          vvoo_ij = real(vvoo_pdm_ij(start_ij+1:end_ij),kind=4)

!$acc enter data copyin(ovoo_ij) async(async_id(1))
!$acc enter data copyin(vvoo_ij) async(async_id(3))

       else if (j .lt. i) then

          ccsd_j = real(ccsd_pdm_j(start_ccsd_j+1:end_ccsd_j),kind=4)
          vvvo_j = real(vvvo_pdm_j(start_vvvo_j+1:end_vvvo_j),kind=4)
          ovoo_ij = real(ovoo(:,:,i,j),kind=4); ovoo_ji = real(ovoo(:,:,j,i),kind=4)
          vvoo_ij = real(vvoo_pdm_ij(start_ij+1:end_ij),kind=4); vvoo_ji = real(vvoo_pdm_ji(start_ji+1:end_ji),kind=4)

!$acc enter data copyin(ccsd_j,vvvo_j,ovoo_ij,ovoo_ji) async(async_id(1))
!$acc enter data copyin(vvoo_ij,vvoo_ji) async(async_id(3))

       endif

    case(3)

       start_ccsd_k = int(i8*k_count*(i8*nocc)*(i8*nvirt**2),kind=8)
       start_vvvo_k = int(i8*k_count*(i8*nvirt)*(i8*nvirt**2),kind=8)
       end_ccsd_k = start_ccsd_k + int(i8*nocc*(i8*nvirt**2),kind=8)
       end_vvvo_k = start_vvvo_k + int(i8*nvirt*(i8*nvirt**2),kind=8)
       start_ik = int((i8*((i_count-1)+k_count*ts_i))*(i8*nvirt**2),kind=8)
       start_ki = int((i8*(k_count+(i_count-1)*ts_k))*(i8*nvirt**2),kind=8)
       start_jk = int((i8*((j_count-1)+k_count*ts_j))*(i8*nvirt**2),kind=8)
       start_kj = int((i8*(k_count+(j_count-1)*ts_k))*(i8*nvirt**2),kind=8)
       end_ik = start_ik + int(i8*nvirt**2,kind=8)
       end_ki = start_ki + int(i8*nvirt**2,kind=8)
       end_jk = start_jk + int(i8*nvirt**2,kind=8)
       end_kj = start_kj + int(i8*nvirt**2,kind=8)
 
       if ((i .eq. j) .and. (j .gt. k)) then

          ccsd_k = real(ccsd_pdm_k(start_ccsd_k+1:end_ccsd_k),kind=4)
          vvvo_k = real(vvvo_pdm_k(start_vvvo_k+1:end_vvvo_k),kind=4)
          ovoo_ik = real(ovoo(:,:,i,k),kind=4); ovoo_ki = real(ovoo(:,:,k,i),kind=4)
          vvoo_ik=real(vvoo_pdm_ik(start_ik+1:end_ik),kind=4); vvoo_ki=real(vvoo_pdm_ki(start_ki+1:end_ki),kind=4)

!$acc enter data copyin(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki) async(async_id(1))
!$acc enter data copyin(vvoo_ik,vvoo_ki) async(async_id(3))

       else if ((i .gt. j) .and. (j .eq. k)) then

          ovoo_jk = real(ovoo(:,:,j,k),kind=4)
          vvoo_jk = real(vvoo_pdm_jk(start_jk+1:end_jk),kind=4)

!$acc enter data copyin(ovoo_jk) async(async_id(1))
!$acc enter data copyin(vvoo_jk) async(async_id(3))

       else if ((i .gt. j) .and. (j .gt. k)) then

          ccsd_k = real(ccsd_pdm_k(start_ccsd_k+1:end_ccsd_k),kind=4)
          vvvo_k = real(vvvo_pdm_k(start_vvvo_k+1:end_vvvo_k),kind=4)
          ovoo_ik = real(ovoo(:,:,i,k),kind=4); ovoo_ki = real(ovoo(:,:,k,i),kind=4)
          ovoo_jk = real(ovoo(:,:,j,k),kind=4); ovoo_kj = real(ovoo(:,:,k,j),kind=4)
          vvoo_ik=real(vvoo_pdm_ik(start_ik+1:end_ik),kind=4); vvoo_ki=real(vvoo_pdm_ki(start_ki+1:end_ki),kind=4)
          vvoo_jk=real(vvoo_pdm_jk(start_jk+1:end_jk),kind=4); vvoo_kj=real(vvoo_pdm_kj(start_kj+1:end_kj),kind=4)

!$acc enter data copyin(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki,ovoo_jk,ovoo_kj) async(async_id(1))
!$acc enter data copyin(vvoo_ik,vvoo_ki,vvoo_jk,vvoo_kj) async(async_id(3))

       end if

    end select

  end subroutine sp_ptr_aliasing_ijk_par
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_aliasing_ijk_ser(nvirt,nocc,i,j,k,&
                   & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                   & ccsd_doubles,vvvo,vvoo,ovoo,async_id,num_ids,case_num)

    implicit none

    !> pointers
    integer :: nvirt,nocc,i,j,k
    real(real_sp), pointer, dimension(:,:,:) :: ccsd_i,ccsd_j,ccsd_k
    real(real_sp), pointer, dimension(:,:,:) :: vvvo_i,vvvo_j,vvvo_k
    real(real_sp), pointer, dimension(:,:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(real_sp), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(realk), dimension(nvirt,nvirt,nocc,nocc), target :: ccsd_doubles
    real(realk), dimension(nvirt,nvirt,nvirt,nocc), target :: vvvo
    real(realk), dimension(nvirt,nvirt,nocc,nocc), target :: vvoo
    real(realk), dimension(nocc,nvirt,nocc,nocc), target :: ovoo
    integer :: num_ids
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_id(num_ids)
#else
    integer :: async_id(num_ids)
#endif
    integer :: case_num

    select case(case_num)

    case(1)

       ccsd_i = real(ccsd_doubles(:,:,:,i),kind=4)
       vvvo_i = real(vvvo(:,:,:,i),kind=4)

!$acc enter data copyin(ccsd_i,vvvo_i) async(async_id(1))

    case(2)

       if (j .eq. i) then

          ovoo_ij = real(ovoo(:,:,i,j),kind=4)
          vvoo_ij = real(vvoo(:,:,i,j),kind=4)

!$acc enter data copyin(ovoo_ij) async(async_id(1))
!$acc enter data copyin(vvoo_ij) async(async_id(3))

       else if (j .lt. i) then

          ccsd_j = real(ccsd_doubles(:,:,:,j),kind=4)
          vvvo_j = real(vvvo(:,:,:,j),kind=4)
          ovoo_ij = real(ovoo(:,:,i,j),kind=4); ovoo_ji = real(ovoo(:,:,j,i),kind=4)
          vvoo_ij = real(vvoo(:,:,i,j),kind=4); vvoo_ji = real(vvoo(:,:,j,i),kind=4)

!$acc enter data copyin(ccsd_j,vvvo_j,ovoo_ij,ovoo_ji) async(async_id(1))
!$acc enter data copyin(vvoo_ij,vvoo_ji) async(async_id(3))

       end if

    case(3)

       if ((i .eq. j) .and. (j .gt. k)) then

          ccsd_k = real(ccsd_doubles(:,:,:,k),kind=4)
          vvvo_k = real(vvvo(:,:,:,k),kind=4)
          ovoo_ik = real(ovoo(:,:,i,k),kind=4); ovoo_ki = real(ovoo(:,:,k,i),kind=4)
          vvoo_ik = real(vvoo(:,:,i,k),kind=4); vvoo_ki = real(vvoo(:,:,k,i),kind=4)

!$acc enter data copyin(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki) async(async_id(1))
!$acc enter data copyin(vvoo_ik,vvoo_ki) async(async_id(3))

       else if ((i .gt. j) .and. (j .eq. k)) then

          ovoo_jk = real(ovoo(:,:,j,k),kind=4)
          vvoo_jk = real(vvoo(:,:,j,k),kind=4)

!$acc enter data copyin(ovoo_jk) async(async_id(1))
!$acc enter data copyin(vvoo_jk) async(async_id(3))

       else if ((i .gt. j) .and. (j .gt. k)) then

          ccsd_k = real(ccsd_doubles(:,:,:,k),kind=4)
          vvvo_k = real(vvvo(:,:,:,k),kind=4)
          ovoo_ik = real(ovoo(:,:,i,k),kind=4); ovoo_ki = real(ovoo(:,:,k,i),kind=4)
          ovoo_jk = real(ovoo(:,:,j,k),kind=4); ovoo_kj = real(ovoo(:,:,k,j),kind=4)
          vvoo_ik = real(vvoo(:,:,i,k),kind=4); vvoo_ki = real(vvoo(:,:,k,i),kind=4)
          vvoo_jk = real(vvoo(:,:,j,k),kind=4); vvoo_kj = real(vvoo(:,:,k,j),kind=4)

!$acc enter data copyin(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki,ovoo_jk,ovoo_kj) async(async_id(1))
!$acc enter data copyin(vvoo_ik,vvoo_ki,vvoo_jk,vvoo_kj) async(async_id(3))

       end if

    end select

  end subroutine sp_ptr_aliasing_ijk_ser
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_aliasing_abc_par(nvirt,nocc,a,b,c,a_count,b_count,c_count,ts_a,ts_b,ts_c,&
                     & ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                     & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                     & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                     & ccsd_pdm_a,ccsd_pdm_b,ccsd_pdm_c,ooov,&
                     & vovv_pdm_a,vovv_pdm_b,vovv_pdm_c,&
                     & oovv_pdm_ab,oovv_pdm_ba,oovv_pdm_ac,oovv_pdm_ca,oovv_pdm_bc,oovv_pdm_cb,async_id,num_ids,case_num)

    implicit none

    integer :: nvirt,nocc,a,b,c,a_count,b_count,c_count,ts_a,ts_b,ts_c,case_num
    real(real_sp), pointer, dimension(:) :: ccsd_a,ccsd_b,ccsd_c
    real(real_sp), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(real_sp), pointer, dimension(:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(real_sp), pointer, dimension(:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), pointer, dimension(:) :: ccsd_pdm_a,ccsd_pdm_b,ccsd_pdm_c
    real(realk), dimension(nocc,nocc,nocc,nvirt), target :: ooov
    real(realk), pointer, dimension(:) :: vovv_pdm_a,vovv_pdm_b,vovv_pdm_c
    real(realk), pointer, dimension(:) :: oovv_pdm_ab,oovv_pdm_ba,oovv_pdm_ac,oovv_pdm_ca,oovv_pdm_bc,oovv_pdm_cb
    integer :: num_ids
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_id(num_ids)
#else
    integer :: async_id(num_ids)
#endif
    integer(kind=long) :: start_ccsd_a,start_ccsd_b,start_ccsd_c,end_ccsd_a,end_ccsd_b,end_ccsd_c
    integer(kind=long) :: start_vovv_ab,start_vovv_ba,start_vovv_ac,start_vovv_ca,start_vovv_bc,start_vovv_cb
    integer(kind=long) :: start_oovv_ab,start_oovv_ba,start_oovv_ac,start_oovv_ca,start_oovv_bc,start_oovv_cb
    integer(kind=long) :: end_vovv_ab,end_vovv_ba,end_vovv_ac,end_vovv_ca,end_vovv_bc,end_vovv_cb
    integer(kind=long) :: end_oovv_ab,end_oovv_ba,end_oovv_ac,end_oovv_ca,end_oovv_bc,end_oovv_cb

    select case(case_num)

    case(1)

       start_ccsd_a = int(i8*a_count*(i8*nvirt)*(i8*nocc**2),kind=8)
       end_ccsd_a = start_ccsd_a + int(i8*nvirt*(i8*nocc**2),kind=8)
   
       ccsd_a = real(ccsd_pdm_a(start_ccsd_a+1:end_ccsd_a),kind=4)
       ooov_a = real(ooov(:,:,:,a),kind=4)

!$acc enter data copyin(ccsd_a,ooov_a) async(async_id(1))

    case(2)

       start_ccsd_b = int(i8*b_count*(i8*nvirt)*(i8*nocc**2),kind=8)
       end_ccsd_b = start_ccsd_b + int(i8*nvirt*(i8*nocc**2),kind=8)
       start_vovv_ab = int((i8*((a-1)+b_count*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_vovv_ba = int((i8*((b-1)+(a_count-1)*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_oovv_ab = int((i8*((a_count-1)+b_count*ts_a))*(i8*nocc**2),kind=8)
       start_oovv_ba = int((i8*(b_count+(a_count-1)*ts_b))*(i8*nocc**2),kind=8)
       end_vovv_ab = start_vovv_ab + int((i8*nvirt)*(i8*nocc),kind=8)
       end_vovv_ba = start_vovv_ba + int((i8*nvirt)*(i8*nocc),kind=8)
       end_oovv_ab = start_oovv_ab + int(i8*nocc**2,kind=8)
       end_oovv_ba = start_oovv_ba + int(i8*nocc**2,kind=8)

       if (b .eq. a) then

          vovv_ab = real(vovv_pdm_b(start_vovv_ab+1:end_vovv_ab),kind=4)
          oovv_ab = real(oovv_pdm_ab(start_oovv_ab+1:end_oovv_ab),kind=4)

!$acc enter data copyin(vovv_ab) async(async_id(1))
!$acc enter data copyin(oovv_ab) async(async_id(3))

       else if (b .lt. a) then

          ccsd_b = real(ccsd_pdm_b(start_ccsd_b+1:end_ccsd_b),kind=4)
          ooov_b = real(ooov(:,:,:,b),kind=4)
          vovv_ab = real(vovv_pdm_b(start_vovv_ab+1:end_vovv_ab),kind=4)
          vovv_ba = real(vovv_pdm_a(start_vovv_ba+1:end_vovv_ba),kind=4)
          oovv_ab = real(oovv_pdm_ab(start_oovv_ab+1:end_oovv_ab),kind=4)
          oovv_ba = real(oovv_pdm_ba(start_oovv_ba+1:end_oovv_ba),kind=4)

!$acc enter data copyin(ccsd_b,ooov_b,vovv_ab,vovv_ba) async(async_id(1))
!$acc enter data copyin(oovv_ab,oovv_ba) async(async_id(3))

       endif

    case(3)

       start_ccsd_c = int(i8*c_count*(i8*nvirt)*(i8*nocc**2),kind=8)
       end_ccsd_c = start_ccsd_c + int(i8*nvirt*(i8*nocc**2),kind=8)
       start_vovv_ac = int((i8*((a-1)+c_count*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_vovv_ca = int((i8*((c-1)+(a_count-1)*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_vovv_bc = int((i8*((b-1)+c_count*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_vovv_cb = int((i8*((c-1)+(b_count-1)*nvirt))*(i8*nvirt)*(i8*nocc),kind=8)
       start_oovv_ac = int((i8*((a_count-1)+c_count*ts_a))*(i8*nocc**2),kind=8)
       start_oovv_ca = int((i8*(c_count+(a_count-1)*ts_c))*(i8*nocc**2),kind=8)
       start_oovv_bc = int((i8*((b_count-1)+c_count*ts_b))*(i8*nocc**2),kind=8)
       start_oovv_cb = int((i8*(c_count+(b_count-1)*ts_c))*(i8*nocc**2),kind=8)
       end_vovv_ac = start_vovv_ac + int((i8*nvirt)*(i8*nocc),kind=8)
       end_vovv_ca = start_vovv_ca + int((i8*nvirt)*(i8*nocc),kind=8)
       end_vovv_bc = start_vovv_bc + int((i8*nvirt)*(i8*nocc),kind=8)
       end_vovv_cb = start_vovv_cb + int((i8*nvirt)*(i8*nocc),kind=8)
       end_oovv_ac = start_oovv_ac + int(i8*nocc**2,kind=8)
       end_oovv_ca = start_oovv_ca + int(i8*nocc**2,kind=8)
       end_oovv_bc = start_oovv_bc + int(i8*nocc**2,kind=8)
       end_oovv_cb = start_oovv_cb + int(i8*nocc**2,kind=8)

       if ((a .eq. b) .and. (b .gt. c)) then
 
          ccsd_c = real(ccsd_pdm_c(start_ccsd_c+1:end_ccsd_c),kind=4)
          ooov_c = real(ooov(:,:,:,c),kind=4)
          vovv_ac = real(vovv_pdm_c(start_vovv_ac+1:end_vovv_ac),kind=4)
          vovv_ca = real(vovv_pdm_a(start_vovv_ca+1:end_vovv_ca),kind=4)
          oovv_ac = real(oovv_pdm_ac(start_oovv_ac+1:end_oovv_ac),kind=4)
          oovv_ca = real(oovv_pdm_ca(start_oovv_ca+1:end_oovv_ca),kind=4)
 
!$acc enter data copyin(ccsd_c,ooov_c,vovv_ac,vovv_ca) async(async_id(1))
!$acc enter data copyin(oovv_ac,oovv_ca) async(async_id(3))

       else if ((a .gt. b) .and. (b .eq. c)) then
 
          vovv_bc = real(vovv_pdm_c(start_vovv_bc+1:end_vovv_bc),kind=4)
          oovv_bc = real(oovv_pdm_bc(start_oovv_bc+1:end_oovv_bc),kind=4)
 
!$acc enter data copyin(vovv_bc) async(async_id(1))
!$acc enter data copyin(oovv_bc) async(async_id(3))

       else if ((a .gt. b) .and. (b .gt. c)) then  
 
          ccsd_c = real(ccsd_pdm_c(start_ccsd_c+1:end_ccsd_c),kind=4)
          ooov_c = real(ooov(:,:,:,c),kind=4)
          vovv_ac = real(vovv_pdm_c(start_vovv_ac+1:end_vovv_ac),kind=4)
          vovv_ca = real(vovv_pdm_a(start_vovv_ca+1:end_vovv_ca),kind=4)
          vovv_bc = real(vovv_pdm_c(start_vovv_bc+1:end_vovv_bc),kind=4)
          vovv_cb = real(vovv_pdm_b(start_vovv_cb+1:end_vovv_cb),kind=4)
          oovv_ac = real(oovv_pdm_ac(start_oovv_ac+1:end_oovv_ac),kind=4)
          oovv_ca = real(oovv_pdm_ca(start_oovv_ca+1:end_oovv_ca),kind=4)
          oovv_bc = real(oovv_pdm_bc(start_oovv_bc+1:end_oovv_bc),kind=4)
          oovv_cb = real(oovv_pdm_cb(start_oovv_cb+1:end_oovv_cb),kind=4)
   
!$acc enter data copyin(ccsd_c,ooov_c,vovv_ac,vovv_ca,vovv_bc,vovv_cb) async(async_id(1))
!$acc enter data copyin(oovv_ac,oovv_ca,oovv_bc,oovv_cb) async(async_id(3))

       end if

    end select

  end subroutine sp_ptr_aliasing_abc_par
#endif

#ifdef VAR_REAL_SP
  subroutine sp_ptr_aliasing_abc_ser(nvirt,nocc,a,b,c,&
                     & ccsd_a,ccsd_b,ccsd_c,&
                     & ooov_a,ooov_b,ooov_c,&
                     & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                     & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                     & ccsd_doubles,ooov,vovv,oovv,async_id,num_ids,case_num)

    implicit none

    integer :: nvirt,nocc,a,b,c
    real(real_sp), pointer, dimension(:,:,:) :: ccsd_a,ccsd_b,ccsd_c
    real(real_sp), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(real_sp), pointer, dimension(:,:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(real_sp), pointer, dimension(:,:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(realk), dimension(nocc,nocc,nvirt,nvirt), target, intent(inout) :: ccsd_doubles
    real(realk), dimension(nocc,nocc,nocc,nvirt), target, intent(inout) :: ooov
    real(realk), dimension(nvirt,nocc,nvirt,nvirt), target, intent(inout) :: vovv
    real(realk), dimension(nocc,nocc,nvirt,nvirt), target, intent(inout) :: oovv
    integer :: num_ids
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_id(num_ids)
#else
    integer :: async_id(num_ids)
#endif
    integer :: case_num

    select case(case_num)

    case(1)

        ccsd_a = real(ccsd_doubles(:,:,:,a),kind=4)
        ooov_a = real(ooov(:,:,:,a),kind=4)

!$acc enter data copyin(ccsd_a,ooov_a) async(async_id(1))

    case(2)

       if (b .eq. a) then

          vovv_ab = real(vovv(:,:,a,b),kind=4)
          oovv_ab = real(oovv(:,:,a,b),kind=4)

!$acc enter data copyin(vovv_ab) async(async_id(1))
!$acc enter data copyin(oovv_ab) async(async_id(3))

       else if (b .lt. a) then

          ccsd_b = real(ccsd_doubles(:,:,:,b),kind=4)
          ooov_b = real(ooov(:,:,:,b),kind=4)
          vovv_ab = real(vovv(:,:,a,b),kind=4); vovv_ba = real(vovv(:,:,b,a),kind=4)
          oovv_ab = real(oovv(:,:,a,b),kind=4); oovv_ba = real(oovv(:,:,b,a),kind=4)

!$acc enter data copyin(ccsd_b,ooov_b,vovv_ab,vovv_ba) async(async_id(1))
!$acc enter data copyin(oovv_ab,oovv_ba) async(async_id(3))

       endif

    case(3)

       if ((a .eq. b) .and. (b .gt. c)) then

          ccsd_c = real(ccsd_doubles(:,:,:,c),kind=4)
          ooov_c = real(ooov(:,:,:,c),kind=4)
          vovv_ac = real(vovv(:,:,a,c),kind=4); vovv_ca = real(vovv(:,:,c,a),kind=4)
          oovv_ac = real(oovv(:,:,a,c),kind=4); oovv_ca = real(oovv(:,:,c,a),kind=4)

!$acc enter data copyin(ccsd_c,ooov_c,vovv_ac,vovv_ca) async(async_id(1))
!$acc enter data copyin(oovv_ac,oovv_ca) async(async_id(3))

       else if ((a .gt. b) .and. (b .eq. c)) then

          vovv_bc = real(vovv(:,:,b,c),kind=4)
          oovv_bc = real(oovv(:,:,b,c),kind=4)

!$acc enter data copyin(vovv_bc) async(async_id(1))
!$acc enter data copyin(oovv_bc) async(async_id(3))

       else if ((a .gt. b) .and. (b .gt. c)) then

          ccsd_c = real(ccsd_doubles(:,:,:,c),kind=4)
          ooov_c = real(ooov(:,:,:,c),kind=4)
          vovv_ac = real(vovv(:,:,a,c),kind=4); vovv_ca = real(vovv(:,:,c,a),kind=4)
          vovv_bc = real(vovv(:,:,b,c),kind=4); vovv_cb = real(vovv(:,:,c,b),kind=4)
          oovv_ac = real(oovv(:,:,a,c),kind=4); oovv_ca = real(oovv(:,:,c,a),kind=4)
          oovv_bc = real(oovv(:,:,b,c),kind=4); oovv_cb = real(oovv(:,:,c,b),kind=4)

!$acc enter data copyin(ccsd_c,ooov_c,vovv_ac,vovv_ca,vovv_bc,vovv_cb) async(async_id(1))
!$acc enter data copyin(oovv_ac,oovv_ca,oovv_bc,oovv_cb) async(async_id(3))

       end if

    end select

  end subroutine sp_ptr_aliasing_abc_ser
#endif

  subroutine preload_tiles_in_bg_buf(array,jobs,b_size,nvirt,nocc,current_i,current_j,current_k,current_ij_count,nbuffs,needed,&
        &tiles_in_buf,array_pdm_buff,req,scheme,tdim,t1dim,dynamic_load,vovo_array)
     implicit none
     type(tensor), intent(inout) :: array
     integer, intent(in) :: b_size, current_i, current_j, current_k, current_ij_count, nbuffs, nvirt, nocc
     integer, intent(in) :: jobs(b_size+1)
     logical, intent(inout) :: needed(nbuffs)
     integer, intent(inout) :: tiles_in_buf(nbuffs)
     real(realk), pointer, intent(inout) :: array_pdm_buff(:,:)
     integer(kind=ls_mpik),intent(inout) :: req(nbuffs)
     integer, intent(in) :: tdim,t1dim
     logical,intent(in) :: scheme,dynamic_load ! if scheme == .true., then ijk, if .false., then abc
     logical, optional, intent(in) :: vovo_array
     integer :: i_test,j_test,k_test,i_search_buf
     integer :: ij_test,ji_test,ik_test,ki_test,jk_test,kj_test
     integer :: ibuf_test, jbuf_test, kbuf_test, ij_count_test, ij_new
     integer :: ijbuf_test,jibuf_test,ikbuf_test,kibuf_test,jkbuf_test,kjbuf_test
     logical :: keep_looping,ij_done, new_i_needed, new_j_needed, new_k_needed
     logical :: found,found_ij,found_ji,found_ik,found_ki,found_jk,found_kj
     logical :: new_ij_needed,new_ji_needed,new_ik_needed,new_ki_needed,new_jk_needed,new_kj_needed
     integer :: ts1,ts2
     integer(kind=ls_mpik) :: mode
     logical :: vovo
     integer :: tuple_type
#ifdef VAR_MPI
     mode = MPI_MODE_NOCHECK

     ! this routine is designed for the ijk scehem.
     ! however, by replacing i, j, k with a, b, and c, respectively, it also works for the abc scheme.

     vovo = .false.
     if (present(vovo_array)) vovo = vovo_array

     !set testing integers
     i_test        = current_i
     j_test        = current_j
     k_test        = current_k + 1 ! we manually increment k by 1
     ij_count_test = current_ij_count
     ij_done       = .false.
     keep_looping  = (count(needed) .lt. nbuffs)

     new_i_needed = .false.; new_i_needed = .false.; new_k_needed = .false.
     new_ij_needed = .false.; new_ji_needed = .false.
     new_ik_needed = .false.; new_ki_needed = .false.
     new_jk_needed = .false.; new_kj_needed = .false.
     found = .false.
     found_ij = .false.; found_ji = .false.
     found_ik = .false.; found_ki = .false.
     found_jk = .false.; found_kj = .false.

     if (scheme) then
        ts2 = nvirt**2*tdim**2 ! vovo ts
     else
        ts2 = nocc**2*tdim**2 ! vovo ts
     endif

     !load next bunch of tiles needed
     fill_buffer: do while(keep_looping)

        !break condition
        if (vovo) then

           keep_looping = (.not. (ij_done .and. (k_test .ge. j_test)))

        else

           keep_looping = (.not. (ij_done .and. (k_test .ge. j_test)) .and. count(needed) .lt. t1dim)

        endif

        ! return condition
        if (count(needed) .eq. nbuffs) return

        if (k_test < j_test) then

           tuple_type = -1

           if (vovo) then

              if ((i_test .eq. j_test) .and. (j_test .gt. k_test)) then

                 ik_test = (k_test-1)*t1dim+i_test; ki_test = (i_test-1)*t1dim+k_test
                 !Load the next ik tile
                 call check_if_new_instance_needed(ik_test,tiles_in_buf,nbuffs,new_ik_needed,set_needed=needed)
                 !Load the next ki tile
                 if (count(needed) .lt. nbuffs) then
                    call check_if_new_instance_needed(ki_test,tiles_in_buf,nbuffs,new_ki_needed,set_needed=needed)
                 else
                    new_ki_needed = .false.
                 endif

                 tuple_type = 1

              else if ((i_test .gt. j_test) .and. (j_test .eq. k_test)) then

                 jk_test = (k_test-1)*t1dim+j_test
                 !Load the next jk tile
                 call check_if_new_instance_needed(jk_test,tiles_in_buf,nbuffs,new_jk_needed,set_needed=needed)

                 tuple_type = 2

              else if ((i_test .gt. j_test) .and. (j_test .gt. k_test)) then

                 ik_test = (k_test-1)*t1dim+i_test; ki_test = (i_test-1)*t1dim+k_test
                 jk_test = (k_test-1)*t1dim+j_test; kj_test = (j_test-1)*t1dim+k_test
                 !Load the next ik tile
                 call check_if_new_instance_needed(ik_test,tiles_in_buf,nbuffs,new_ik_needed,set_needed=needed)
                 !Load the next ki tile
                 if (count(needed) .lt. nbuffs) then
                    call check_if_new_instance_needed(ki_test,tiles_in_buf,nbuffs,new_ki_needed,set_needed=needed)
                 else
                    new_ki_needed = .false.
                 endif
                 !Load the next jk tile
                 if (count(needed) .lt. nbuffs) then
                    call check_if_new_instance_needed(jk_test,tiles_in_buf,nbuffs,new_jk_needed,set_needed=needed)
                 else
                    new_jk_needed = .false.
                 endif
                 !Load the next kj tile
                 if (count(needed) .lt. nbuffs) then
                    call check_if_new_instance_needed(kj_test,tiles_in_buf,nbuffs,new_kj_needed,set_needed=needed)
                 else
                    new_kj_needed = .false.
                 endif

                 tuple_type = 3

              endif

           else

              !Load the next k tile
              call check_if_new_instance_needed(k_test,tiles_in_buf,nbuffs,new_k_needed,set_needed=needed)

           endif

           !load k
           if (new_k_needed .or. (tuple_type .ge. 1)) then

              if (vovo) then

                 if (tuple_type .eq. 1) then

                    !find pos in buff
                    if (new_ik_needed) then
                       call find_free_pos_in_buf(needed,nbuffs,ikbuf_test,found_ik)
                       if (found_ik) then
                          needed(ikbuf_test) = .true.
                          tiles_in_buf(ikbuf_test) = ik_test
                       endif
                    else
                       found_ik = .false.
                    endif
                    if (new_ki_needed .and. count(needed) .lt. nbuffs) then
                       call find_free_pos_in_buf(needed,nbuffs,kibuf_test,found_ki)
                       if (found_ki) then
                          needed(kibuf_test) = .true.
                          tiles_in_buf(kibuf_test) = ki_test
                       endif
                    else
                       found_ki = .false.
                    endif

                 else if (tuple_type .eq. 2) then

                    !find pos in buff
                    if (new_jk_needed) then
                       call find_free_pos_in_buf(needed,nbuffs,jkbuf_test,found_jk)
                       if (found_jk) then
                          needed(jkbuf_test) = .true.
                          tiles_in_buf(jkbuf_test) = jk_test
                       endif
                    else
                       found_jk = .false.
                    endif

                 else if (tuple_type .eq. 3) then

                    !find pos in buff
                    if (new_ik_needed) then
                       call find_free_pos_in_buf(needed,nbuffs,ikbuf_test,found_ik)
                       if (found_ik) then
                          needed(ikbuf_test) = .true.
                          tiles_in_buf(ikbuf_test) = ik_test
                       endif
                    else
                       found_ik = .false.
                    endif
                    if (new_ki_needed .and. count(needed) .lt. nbuffs) then
                       call find_free_pos_in_buf(needed,nbuffs,kibuf_test,found_ki)
                       if (found_ki) then
                          needed(kibuf_test) = .true.
                          tiles_in_buf(kibuf_test) = ki_test
                       endif
                    else
                       found_ki = .false.
                    endif
                    if (new_jk_needed .and. count(needed) .lt. nbuffs) then
                       call find_free_pos_in_buf(needed,nbuffs,jkbuf_test,found_jk)
                       if (found_jk) then
                          needed(jkbuf_test) = .true.
                          tiles_in_buf(jkbuf_test) = jk_test
                       endif
                    else
                       found_jk = .false.
                    endif
                    if (new_kj_needed .and. count(needed) .lt. nbuffs) then
                       call find_free_pos_in_buf(needed,nbuffs,kjbuf_test,found_kj)
                       if (found_kj) then
                          needed(kjbuf_test) = .true.
                          tiles_in_buf(kjbuf_test) = kj_test
                       endif
                    else
                       found_kj = .false.
                    endif

                 endif

                 if (found_ik .or. found_ki .or. found_jk .or. found_kj) found = .true.

              else

                 !find pos in buff
                 call find_free_pos_in_buf(needed,nbuffs,kbuf_test,found)

              endif

              if(found)then

                 if (vovo) then

                    if (tuple_type .eq. 1) then

                       if (.not. alloc_in_dummy) then
                          if (found_ik) call tensor_lock_win(array,ik_test,'s',assert=mode)
                          if (found_ki) call tensor_lock_win(array,ki_test,'s',assert=mode)
                       endif

                       if( alloc_in_dummy )then
                          if (found_ik) call tensor_get_tile(array,ik_test,array_pdm_buff(:,ikbuf_test),ts2,&
                             &lock_set=.true.,req=req(ikbuf_test))
                          if (found_ki) call tensor_get_tile(array,ki_test,array_pdm_buff(:,kibuf_test),ts2,&
                             &lock_set=.true.,req=req(kibuf_test))
                       else
                          if (found_ik) call tensor_get_tile(array,ik_test,array_pdm_buff(:,ikbuf_test),ts2,&
                             &lock_set=.true.,flush_it=.true.)
                          if (found_ki) call tensor_get_tile(array,ki_test,array_pdm_buff(:,kibuf_test),ts2,&
                             &lock_set=.true.,flush_it=.true.)
                       endif

                    else if (tuple_type .eq. 2) then

                       if (.not. alloc_in_dummy) call tensor_lock_win(array,jk_test,'s',assert=mode)

                       if( alloc_in_dummy )then
                          call tensor_get_tile(array,jk_test,array_pdm_buff(:,jkbuf_test),ts2,&
                             &lock_set=.true.,req=req(jkbuf_test))
                       else
                          call tensor_get_tile(array,jk_test,array_pdm_buff(:,jkbuf_test),ts2,&
                             &lock_set=.true.,flush_it=.true.)
                       endif

                    else if (tuple_type .eq. 3) then

                       if (.not. alloc_in_dummy) then
                          if (found_ik) call tensor_lock_win(array,ik_test,'s',assert=mode)
                          if (found_ki) call tensor_lock_win(array,ki_test,'s',assert=mode)
                          if (found_jk) call tensor_lock_win(array,jk_test,'s',assert=mode)
                          if (found_kj) call tensor_lock_win(array,kj_test,'s',assert=mode)
                       endif

                       if( alloc_in_dummy )then
                          if (found_ik) call tensor_get_tile(array,ik_test,array_pdm_buff(:,ikbuf_test),ts2,&
                             &lock_set=.true.,req=req(ikbuf_test))
                          if (found_ki) call tensor_get_tile(array,ki_test,array_pdm_buff(:,kibuf_test),ts2,&
                             &lock_set=.true.,req=req(kibuf_test))
                          if (found_jk) call tensor_get_tile(array,jk_test,array_pdm_buff(:,jkbuf_test),ts2,&
                             &lock_set=.true.,req=req(jkbuf_test))
                          if (found_kj) call tensor_get_tile(array,kj_test,array_pdm_buff(:,kjbuf_test),ts2,&
                             &lock_set=.true.,req=req(kjbuf_test))
                       else
                          if (found_ik) call tensor_get_tile(array,ik_test,array_pdm_buff(:,ikbuf_test),ts2,&
                             &lock_set=.true.,flush_it=.true.)
                          if (found_ki) call tensor_get_tile(array,ki_test,array_pdm_buff(:,kibuf_test),ts2,&
                             &lock_set=.true.,flush_it=.true.)
                          if (found_jk) call tensor_get_tile(array,jk_test,array_pdm_buff(:,jkbuf_test),ts2,&
                             &lock_set=.true.,flush_it=.true.)
                          if (found_kj) call tensor_get_tile(array,kj_test,array_pdm_buff(:,kjbuf_test),ts2,&
                             &lock_set=.true.,flush_it=.true.)
                       endif

                    endif

                 else

                    if( .not.alloc_in_dummy ) call tensor_lock_win(array,k_test,'s',assert=mode)
                    call get_tile_dim(ts1,array,k_test)
                    if( alloc_in_dummy )then
                       call tensor_get_tile(array,k_test,array_pdm_buff(:,kbuf_test),ts1,&
                          &lock_set=.true.,req=req(kbuf_test))
                    else
                       call tensor_get_tile(array,k_test,array_pdm_buff(:,kbuf_test),ts1,&
                          &lock_set=.true.,flush_it=.true.)
                    endif

                    needed(kbuf_test)       = .true.
                    tiles_in_buf(kbuf_test) = k_test

                 endif

              endif

           endif

        else ! k_test > j_test

           !Load the next i and j tiles
           if(.not.dynamic_load)then

              ij_count_test = ij_count_test + 1

              if(ij_count_test<=b_size)then

                 !is incremented by one at the end of the loop
                 !therefore we set it to 0 here
                 k_test = 0

                 ij_new = jobs(ij_count_test)

                 call calc_i_leq_j(ij_new,t1dim,i_test,j_test)

                 call check_if_new_instance_needed(j_test,tiles_in_buf,nbuffs,new_j_needed,set_needed=needed)

                 !load new j
                 if (new_j_needed .or. vovo) then

                    if (vovo) then

                       if (i_test .eq. j_test) then

                          ij_test = (j_test-1)*t1dim+i_test
                          !Load the next ij tile
                          call check_if_new_instance_needed(ij_test,tiles_in_buf,nbuffs,new_ij_needed,set_needed=needed)

                          !find pos in buff
                          if (new_ij_needed) then
                             call find_free_pos_in_buf(needed,nbuffs,ijbuf_test,found_ij)
                             if (found_ij) then
                                needed(ijbuf_test) = .true.
                                tiles_in_buf(ijbuf_test) = ij_test
                                found_ji = .false. 
                             endif
                          else
                             found_ij = .false.; found_ji = .false.
                          endif

                       else

                          ij_test = (j_test-1)*t1dim+i_test; ji_test = (i_test-1)*t1dim+j_test
                          !Load the next ij tile
                          call check_if_new_instance_needed(ij_test,tiles_in_buf,nbuffs,new_ij_needed,set_needed=needed)
                          !Load the next ji tile
                          if (count(needed) .lt. nbuffs) then
                             call check_if_new_instance_needed(ji_test,tiles_in_buf,nbuffs,new_ji_needed,set_needed=needed)
                          else
                             new_ji_needed = .false.
                          endif

                          !find pos in buff
                          if (new_ij_needed) then
                             call find_free_pos_in_buf(needed,nbuffs,ijbuf_test,found_ij)
                             if (found_ij) then
                                needed(ijbuf_test) = .true.
                                tiles_in_buf(ijbuf_test) = ij_test 
                             endif
                          else
                             found_ij = .false.
                          endif
                          if (new_ji_needed .and. count(needed) .lt. nbuffs) then
                             call find_free_pos_in_buf(needed,nbuffs,jibuf_test,found_ji)
                             if (found_ji) then
                                needed(jibuf_test) = .true.
                                tiles_in_buf(jibuf_test) = ji_test
                             endif
                          else
                             found_ji = .false.
                          endif

                       endif

                       if (found_ij .or. found_ji) found = .true.

                    else

                       !find pos in buff
                       call find_free_pos_in_buf(needed,nbuffs,jbuf_test,found)

                    endif

                    if(found)then

                       if (vovo) then

                          if (.not. alloc_in_dummy) then
                             if (found_ij) call tensor_lock_win(array,ij_test,'s',assert=mode)
                             if (found_ji) call tensor_lock_win(array,ji_test,'s',assert=mode)
                          endif
                          if( alloc_in_dummy )then
                             if (found_ij) call tensor_get_tile(array,ij_test,array_pdm_buff(:,ijbuf_test),ts2,&
                                &lock_set=.true.,req=req(ijbuf_test))
                             if (found_ji) call tensor_get_tile(array,ji_test,&
                                &array_pdm_buff(:,jibuf_test),ts2,lock_set=.true.,req=req(jibuf_test))
                          else
                             if (found_ij) call tensor_get_tile(array,ij_test,array_pdm_buff(:,ijbuf_test),ts2,&
                                &lock_set=.true.,flush_it=.true.)
                             if (found_ji) call tensor_get_tile(array,ji_test,&
                                &array_pdm_buff(:,jibuf_test),ts2,lock_set=.true.,flush_it=.true.)
                          endif

                       else

                          if( .not. alloc_in_dummy ) call tensor_lock_win(array,j_test,'s',assert=mode)
                          call get_tile_dim(ts1,array,j_test)
                          if(alloc_in_dummy)then
                             call tensor_get_tile(array,j_test,array_pdm_buff(:,jbuf_test),ts1,&
                                &lock_set=.true.,req=req(jbuf_test))
                          else
                             call tensor_get_tile(array,j_test,array_pdm_buff(:,jbuf_test),ts1,&
                                &lock_set=.true.,flush_it=.true.)
                          endif

                          needed(jbuf_test)       = .true.
                          tiles_in_buf(jbuf_test) = j_test

                       endif

                    endif

                 endif

                 if (.not. vovo) then

                    call check_if_new_instance_needed(i_test,tiles_in_buf,nbuffs,new_i_needed,set_needed=needed)

                    !load new i
                    if( new_i_needed )then

                       !find pos in buff
                       call find_free_pos_in_buf(needed,nbuffs,ibuf_test,found)

                       if(found)then

                          if( .not. alloc_in_dummy ) call tensor_lock_win(array,i_test,'s',assert=mode)
                          call get_tile_dim(ts1,array,i_test)
                          if(alloc_in_dummy )then
                             call tensor_get_tile(array,i_test,array_pdm_buff(:,ibuf_test),ts1,&
                                &lock_set=.true.,req=req(ibuf_test))
                          else
                             call tensor_get_tile(array,i_test,array_pdm_buff(:,ibuf_test),ts1,&
                                &lock_set=.true.,flush_it=.true.)
                          endif

                          needed(ibuf_test)       = .true.
                          tiles_in_buf(ibuf_test) = i_test

                       endif

                    endif

                 endif

              else

                 ij_done = .true.

              endif

           else

              ij_done = .true. 

           endif


        endif

        k_test = k_test + 1

     enddo fill_buffer
#endif
  end subroutine preload_tiles_in_bg_buf

  !> \brief: create ij_array for ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2013
  subroutine create_comp_array_ccsdpt(njobs,ts,array)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: njobs,ts
    !> ij_array
    integer, dimension(njobs), intent(inout) :: array
    !> integers
    integer :: counter,offset,fill_1,fill_2

    ! for ijk scheme, ts == nocc

    ! since i .ge. j, the composite ij indices will make up a lower triangular matrix.
    ! for each ij, k (where j .ge. k) jobs have to be carried out.
    ! thus, the largest jobs for a given i-value will be those that have the largest j-value,
    ! i.e. the largest jobs will be those for which the ij index appears near the diagonal.
    ! as the value of j specifies how large a given job is, we fill up the ij_array with jobs
    ! for j-values in descending order.

    ! the below is the lower triangular part of the ij (5*5) matrix written in row-major order

    ! ||   1              ||
    ! ||   2  3           ||
    ! ||   4  5  6        ||
    ! ||   7  8  9 10     ||
    ! ||  11 12 13 14 15  ||

    ! examples of ij --> i,j conversion
    ! - ij index 15 corresponds to (i,j)=(5,5) and thus to k=1,2,3,4,5
    ! - ij index 9  corresponds to (i,j)=(4,3) and thus to k=1,2,3
    ! - ij index 11  corresponds to (i,j)=(5,1) and thus to k=1

    ! we want ij_array to look like this
    ! (15 , 14 , 10 , 13 , 9 , 6 , 12 , 8 , 5 , 3 , 11 , 7 , 4 , 2 , 1)

    ! counter specifies the index of ij_array
    counter = 1

    do fill_1 = 0,ts-1

       ! zero the offset
       offset = 0

       if (fill_1 .eq. 0) then

          ! this is largest possible job, i.e., the (no,no)-th entry in the ij matrix
          array(counter) = njobs
          ! increment counter
          counter = counter + 1

       else

          do fill_2 = 0,fill_1

             if (fill_2 .eq. 0) then

                ! this is the largest i-value, for which we have to do k number of jobs, 
                ! that is, we are at the no'th row essentially moving from right towards left.
                array(counter) = njobs - fill_1
                ! increment counter
                counter = counter + 1

             else

                ! we loop through the i-values keeping the j-value (and k-range) fixed
                ! we thus loop from i == no up towards the diagonal of the lower triangular matrix
                offset = offset + (ts - fill_2)
                ! 'njobs - fill_1' gives the current column, while 'offset' moves us up through the rows
                ! while staying below or on the diagonal.(still row-major numbering)
                array(counter) = njobs - fill_1 - offset
                ! increment counter
                counter = counter + 1

             end if

          end do

       end if

    end do

  end subroutine create_comp_array_ccsdpt


  !> \brief: make job distribution list for ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2013
  subroutine job_distrib_ccsdpt(b_size,njobs,index_array,jobs)

    implicit none

    !> batch size (without remainder contribution) and njobs 
    integer, intent(in) :: b_size,njobs
    !> index_array
    integer, dimension(njobs), intent(inout) :: index_array
    !> jobs array
    integer, dimension(b_size+1), intent(inout) :: jobs
    !> integers
    integer :: nodtotal,fill,fill_sum

#ifdef VAR_MPI

    nodtotal = infpar%lg_nodtot

    ! fill the jobs array with composite index values stored in index_array.
    ! there are njobs jobs in total.

    ! the below algorithm distributes the jobs evenly among the nodes.

    do fill = 0,b_size

       fill_sum = infpar%lg_mynum + 1 + fill*nodtotal

       if (fill_sum .le. njobs) then

          jobs(fill + 1) = index_array(fill_sum) 

       else

          ! fill jobs array with negative number such that this number won't appear for any value of the composite index
          jobs(fill + 1) = -1

       end if

    end do

#endif

  end subroutine job_distrib_ccsdpt

!endif mod_unreleased
#endif

  subroutine dummy_ccsdpt_tools_routine()

  end subroutine dummy_ccsdpt_tools_routine

end module ccsdpt_tools_module
