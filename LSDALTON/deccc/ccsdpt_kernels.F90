!> @file
!> DEC-CCSD(T) kernels
!> \brief: ccsd(t) module
!> \author: Janus Juul Eriksen
!> \date: 2012-2014, Aarhus
module ccsdpt_kernels_module
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer, c_size_t

#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  use precision
  use dec_typedef_module
  use memory_handling
  use lstiming!, only: lstimer
  use Fundamental, only: bohr_to_angstrom
  use tensor_interface_module
  use reorder_frontend_module
  use background_buffer_module
#ifdef VAR_OPENACC
  use openacc
#endif
  use gpu_interfaces

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use ccsdpt_tools_module
  use ccsdpt_full_module
  use ccsdpt_dec_module
#ifdef VAR_REAL_SP
  use sp_ccsdpt_full_module
  use sp_ccsdpt_dec_module
#endif
  use decmpi_module
  use dec_workarounds_module
  use crop_tools_module
  use cc_tools_module
  use dec_fragment_utils

#ifdef MOD_UNRELEASED
  public :: ijk_loop_par
  public :: ijk_loop_ser
  public :: abc_loop_par
  public :: abc_loop_ser
#endif

  private

contains

#ifdef MOD_UNRELEASED

#ifdef VAR_MPI
  !> \brief: main ijk-loop (mpi version)
  !> \author: Janus Juul Eriksen
  !> \date: january 2014
  subroutine ijk_loop_par(nocc,nvirt,ovoo,vvoo,vvvo,ccsd_doubles,&
                        & eivalocc,eivalvirt,nodtotal,nbuffs,tile_size,ccsdpt_singles,&
                        & ccsdpt_doubles,e4,e5,t1)

    implicit none

    !> nocc,nvirt
    integer, intent(in) :: nocc,nvirt
    !> 2-el integrals
    real(realk), dimension(nocc,nvirt,nocc,nocc), target :: ovoo ! integrals (AI|JK) in the order (J,A,I,K)
    type(tensor), intent(inout) :: vvoo ! integrals (AI|BJ) in the order (A,B,I,J)
    real(realk), pointer, dimension(:) :: vvoo_pdm_ij,vvoo_pdm_ji,vvoo_pdm_ik,vvoo_pdm_ki ! v^2*tile_size tiles from vvoo
    real(realk), pointer, dimension(:) :: vvoo_pdm_jk,vvoo_pdm_kj ! v^2*tile_size tiles from vvoo
    real(realk), pointer, dimension(:,:) :: vvoo_pdm_buff      ! buffers to prefetch vvoo tiles
    type(tensor), intent(inout) :: vvvo ! integrals (AI|BC) in the order (C,B,A,I)
    real(realk), pointer, dimension(:) :: vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k ! v^3*tile_size tiles from vvvo
    real(realk), pointer, dimension(:,:) :: vvvo_pdm_buff      ! buffers to prefetch vvvo tiles
    integer, intent(inout) :: nodtotal, tile_size
    !> ccsd doubles amplitudes
    type(tensor), intent(inout)  :: ccsd_doubles
    real(realk), pointer, dimension(:) :: ccsd_pdm_i,ccsd_pdm_j,ccsd_pdm_k ! ov^2*tile_size tiles from ccsd_doubles
    real(realk), pointer, dimension(:,:) :: ccsd_pdm_buff ! buffers to prefetch ccsd_doubles tiles
    !> triples amplitudes and 3d work array
    real(real_pt), pointer, dimension(:,:,:) :: trip_tmp, trip_ampl
    !> ccsd(t) intermediates
    real(realk), dimension(nvirt,nocc), optional :: ccsdpt_singles
    real(realk), dimension(nvirt,nvirt,nocc,nocc), optional :: ccsdpt_doubles
    real(realk),optional :: e4,e5
    real(realk), dimension(nvirt,nocc), target, optional :: t1 
    logical :: full_no_frags
    real(real_pt), pointer, dimension(:) :: tmp_res_e5
    !> pointers
    real(real_pt), pointer, dimension(:) :: ccsd_i,ccsd_j,ccsd_k
    real(real_pt), pointer, dimension(:) :: vvvo_i,vvvo_j,vvvo_k
    real(real_pt), pointer, dimension(:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(real_pt), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(real_pt), pointer, dimension(:,:) :: t1_ptr
    !> tmp pointers
    real(real_pt), pointer :: pt_1(:,:),pt_2(:,:,:,:)
    !> orbital energiesi
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt) 
    !> loop integers
    integer :: b_size,njobs,ij_comp,ij_count,ij_max_count
    integer, pointer :: ij_array(:),jobs(:)
    integer :: i,j,k,tuple_type
    integer :: i_tile,j_tile,k_tile,i_pos,j_pos,k_pos,i_count,j_count,k_count
    integer :: i_buf_vvvo,i_buf_ccsd,j_buf_vvvo,j_buf_ccsd,k_buf_vvvo,k_buf_ccsd
    integer :: ijbuf,jibuf,ikbuf,kibuf,jkbuf,kjbuf
    integer :: ij,ji,ik,ki,jk,kj
    integer :: total_num_tiles_1,total_num_tiles_2,total_num_tiles,dim_ts
    integer :: nelms,tile_size_tmp_i,tile_size_tmp_j,tile_size_tmp_k
    !> preloading
    integer, intent(in) :: nbuffs
    integer,pointer, dimension(:) :: tiles_in_buf_vvvo,tiles_in_buf_ccsd,tiles_in_buf_vvoo
    integer(kind=ls_mpik), pointer, dimension(:) :: req_vvvo,req_ccsd,req_vvoo
    logical,pointer,dimension(:) :: needed_vvvo,needed_ccsd,needed_vvoo
    !> async handles
    integer :: num_ids,m
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
    integer(kind=acc_device_kind) :: acc_device_type
#else
    integer, pointer, dimension(:) :: async_id
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat
    ! timings
    real(realk) :: tcpu,twall,time_pt_ijk,time_pt_ijk_min,time_pt_ijk_max
    real(realk) :: time_trip,time_efull,time_driv,time_preload
    real(realk) :: time_trip_min,time_efull_min,time_driv_min,time_preload_min
    real(realk) :: time_trip_max,time_efull_max,time_driv_max,time_preload_max
    real(realk) :: time_trip_tot,time_efull_tot,time_driv_tot,time_preload_tot
    real(realk) :: phase_cntrs(nphases)
    real(realk) :: flushing_time, flushing_time_min, flushing_time_max
    real(realk) :: unlock_time,   unlock_time_min,   unlock_time_max
    real(realk) :: waiting_time,  waiting_time_min,  waiting_time_max
    real(realk) :: time_w_min, time_w_max
    real(realk) :: time_c_min, time_c_max
    real(realk) :: time_i_min, time_i_max
    logical :: use_bg_buf, dynamic_load
    ! remote counter and additional variables for dynamic load
    integer, pointer      :: dyn_i(:)
    type(c_ptr)           :: dyn_c
    integer(kind=ls_mpik) :: dyn_w,mode
    integer               :: plus_one

    mode = MPI_MODE_NOCHECK 

    ! init timings
    unlock_time   = time_lsmpi_win_unlock
    waiting_time  = time_lsmpi_wait
    flushing_time = time_lsmpi_win_flush
    time_trip_tot = 0.0E0_realk; time_preload_tot = 0.0E0_realk; time_efull_tot = 0.0E0_realk; time_driv_tot = 0.0E0_realk
    call time_start_phase( PHASE_WORK, twall = time_pt_ijk )
    call time_phases_get_current(current_wt=phase_cntrs)
    if (infpar%lg_mynum .eq. infpar%master) call LSTIMER('START',tcpu,twall,DECinfo%output)

    full_no_frags = .false.
    use_bg_buf    = mem_is_background_buf_init()
    dynamic_load  = DECinfo%dyn_load
    dynamic_load  = .false.
    plus_one      = 1

    if (present(e4) .and. present(e5)) full_no_frags = .true.

    if (full_no_frags) then

       call mem_alloc(tmp_res_e5,i8*nvirt)

       call ptr_init_ijk_par(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                      & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                      & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1,t1_ptr)

    else

       call ptr_init_ijk_par(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                      & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                      & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj)

       call ptr_init_ijk_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    endif

    call time_start_phase(PHASE_WORK)

    ! alloc and init stuff for preloading
    if(use_bg_buf)then
       call mem_pseudo_alloc(vvvo_pdm_buff,int((i8*nvirt)*(i8*nvirt**2)*tile_size,kind=8),int(i8*3*nbuffs,kind=8))
       call mem_pseudo_alloc(ccsd_pdm_buff,int(i8*nocc*(i8*nvirt**2)*tile_size,kind=8),int(i8*3*nbuffs,kind=8))
       call mem_pseudo_alloc(vvoo_pdm_buff,int((i8*nvirt**2)*tile_size**2,kind=8),int(i8*6*nbuffs,kind=8))
    else
       call mem_alloc(vvvo_pdm_buff,nvirt**3*tile_size,3*nbuffs)
       call mem_alloc(ccsd_pdm_buff,nocc*nvirt**2*tile_size,3*nbuffs)
       call mem_alloc(vvoo_pdm_buff,nvirt**2*tile_size**2,6*nbuffs)
    endif
    call mem_alloc(needed_vvvo,3*nbuffs)
    call mem_alloc(needed_ccsd,3*nbuffs)
    call mem_alloc(needed_vvoo,6*nbuffs)
    call mem_alloc(tiles_in_buf_vvvo,3*nbuffs)
    call mem_alloc(tiles_in_buf_ccsd,3*nbuffs)
    call mem_alloc(tiles_in_buf_vvoo,6*nbuffs)
    call mem_alloc(req_vvvo,3*nbuffs)
    call mem_alloc(req_ccsd,3*nbuffs)
    call mem_alloc(req_vvoo,6*nbuffs)
    if (alloc_in_dummy) then
       call tensor_lock_wins(vvvo,'s',all_nodes=.true.)
       call tensor_lock_wins(ccsd_doubles,'s',all_nodes=.true.)
       call tensor_lock_wins(vvoo,'s',all_nodes=.true.)
    endif
    needed_vvvo       = .false.
    needed_ccsd       = .false.
    needed_vvoo       = .false.
    tiles_in_buf_vvvo = -1
    tiles_in_buf_ccsd = -1
    tiles_in_buf_vvoo = -1

#ifndef VAR_REAL_SP
    if(use_bg_buf)then
       ! init triples tuples structure
       call mem_pseudo_alloc(trip_ampl,i8*nvirt,i8*nvirt,i8*nvirt)
       ! init 3d wrk array
       call mem_pseudo_alloc(trip_tmp,i8*nvirt,i8*nvirt,i8*nvirt)
    else
       ! init triples tuples structure
       call mem_alloc(trip_ampl,nvirt,nvirt,nvirt)
       ! init 3d wrk array
       call mem_alloc(trip_tmp,nvirt,nvirt,nvirt)
    endif
#else
    ! init triples tuples structure
    call mem_alloc(trip_ampl,nvirt,nvirt,nvirt)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nvirt,nvirt,nvirt)
#endif

    ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
    ! handle 1: ccsd_doubles and vvvo / ovoo integrals
    ! handle 3: vvoo integrals and ccsdpt_doubles intermediate
    ! handle 4: triples amplitudes
    ! handle 5: energy evaluation
    num_ids = 5
    call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

    if (DECinfo%acc_sync) then
       async_id = acc_async_sync
    else
       do m = 1,num_ids
          async_id(m) = int(m,kind=acc_handle_kind)
       enddo
    endif

#else

    if (DECinfo%acc_sync) then
       async_id = 0
    else
       do m = 1,num_ids
          async_id(m) = -m
       enddo
    endif

#endif

#ifdef VAR_CUBLAS

    ! initialize the CUBLAS context
    stat = cublasCreate_v2(cublas_handle)

#endif

    total_num_tiles_1 = vvvo%ntiles
    total_num_tiles_2 = ccsd_doubles%ntiles
    if (total_num_tiles_1 .ne. total_num_tiles_2) call lsquit('total_num_tiles_1 .ne. total_num_tiles_2 (ijk)',DECinfo%output) 
    total_num_tiles = total_num_tiles_1
    dim_ts = int(nocc / tile_size)
    if (mod(nocc,tile_size) .gt. 0) dim_ts = dim_ts + 1 

    if ((DECinfo%PL .gt. 2) .and. (infpar%lg_mynum .eq. infpar%master)) then
       print *,'nocc = ',nocc
       print *,'nvirt = ',nvirt
       print *,'total_num_tiles_1 = ',total_num_tiles_1
       print *,'total_num_tiles_2 = ',total_num_tiles_2
       print *,'tile_size = ',tile_size
       print *,'mod(nocc,tile_size) = ',mod(nocc,tile_size)
       print *,'dim_ts = ',dim_ts
    endif

    i_count = 0
    j_count = 0
    k_count = 0

    tile_size_tmp_i = 0
    tile_size_tmp_j = 0
    tile_size_tmp_k = 0

    ! create job distribution list
    ! first, determine common batch size from number of tasks and nodes
    ! in the ij matrix, njobs is the number of elements in the lower triangular matrix
    ! always an even number [ n(n+1) is always an even number ]
    njobs = int((dim_ts**2 + dim_ts)/2)
    b_size = int(njobs/nodtotal)

    if ((DECinfo%PL .gt. 2) .and. (infpar%lg_mynum .eq. infpar%master)) then
       print *,'nodtotal = ',nodtotal
       print *,'njobs = ',njobs
       print *,'b_size = ',b_size
    endif

    ! ij_array stores all jobs for composite ab indices in descending order
    call mem_alloc(ij_array,njobs)
    ! init list (one more than b_size since mod(njobs,nodtotal) is not necessearily zero
    call mem_alloc(jobs,b_size + 1)

    ! create ij_array
    call create_comp_array_ccsdpt(njobs,dim_ts,ij_array)
    ! fill the list
    call job_distrib_ccsdpt(b_size,njobs,ij_array,jobs)

    if(dynamic_load)then
       ij_count=0
       call mem_alloc( dyn_i, dyn_c, 1)

       dyn_i = 0
       if(infpar%lg_mynum == 0) dyn_i(1) = infpar%lg_nodtot

       call lsmpi_win_create(dyn_i,dyn_w,1,infpar%lg_comm)
#ifdef VAR_HAVE_MPI3
       call lsmpi_win_lock_all(dyn_w,ass=mode)
#else
       call lsquit("ERROR(ijk_loop_par): dynamic load only implemented for MPI3",-1)
#endif
       ij_max_count = njobs
    else
       ij_count=0
       ij_max_count = b_size + 1
    endif

    ! now follows the main loop, which is collapsed.

!$acc enter data create(trip_tmp,trip_ampl,tmp_res_e5)&
!$acc& copyin(eivalvirt,t1_ptr,e4,e5) if(full_no_frags)
!
!$acc enter data create(trip_tmp,trip_ampl)&
!$acc& copyin(eivalvirt,ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

    ijrun_par: do while (ij_count <= ij_max_count)

          !Get Job index
          if(dynamic_load)then

             if(ij_count == 0) then
                ij_count = infpar%lg_mynum + 1
             else
                call lsmpi_get_acc(plus_one,ij_count,infpar%master,1,dyn_w)
                call lsmpi_win_flush(dyn_w,local=.true.)
                ij_count = ij_count + 1
             endif

             if(ij_count <= njobs)then
                ij_comp  = ij_array(ij_count)
             else
                exit ijrun_par
             endif

          else
             ij_count = ij_count + 1
             if (ij_count > b_size + 1) exit ijrun_par
             ij_comp  = jobs(ij_count)
             ! no more jobs to be done? otherwise leave the loop
             if (ij_comp .lt. 0) exit ijrun_par
          endif

          ! calculate i and j from composite ij value
          call calc_i_leq_j(ij_comp,dim_ts,i_tile,j_tile)

          i_pos = (i_tile-1)*tile_size+1
          j_pos = (j_tile-1)*tile_size+1

          call get_tile_dim(nelms,vvvo,i_tile)
          tile_size_tmp_i = int(nelms/nvirt**3)
          call get_tile_dim(nelms,vvvo,j_tile)
          tile_size_tmp_j = int(nelms/nvirt**3)

          if((DECinfo%PL .gt. 2) .and. (ij_comp .gt. 0)) &
             & write(*,'("Rank ",I3," does PT ij job in ijk loop: (",I5"/",I5"/",I5"/",I5"/",I5"/",I5"/",I5"/",I5,")")') &
             &infpar%lg_mynum,ij_comp,njobs,i_tile,j_tile,tile_size_tmp_i,tile_size_tmp_j,i_pos,j_pos

          !FIND i and j in buffer
          call assoc_ptr_to_buf(i_tile,vvvo,3*nbuffs,tiles_in_buf_vvvo,needed_vvvo,&
                               & vvvo_pdm_i,vvvo_pdm_buff,i_buf_vvvo,req_vvvo)
          call assoc_ptr_to_buf(j_tile,vvvo,3*nbuffs,tiles_in_buf_vvvo,needed_vvvo,&
                               & vvvo_pdm_j,vvvo_pdm_buff,j_buf_vvvo,req_vvvo)
          call assoc_ptr_to_buf(i_tile,ccsd_doubles,3*nbuffs,tiles_in_buf_ccsd,needed_ccsd,&
                               & ccsd_pdm_i,ccsd_pdm_buff,i_buf_ccsd,req_ccsd)
          call assoc_ptr_to_buf(j_tile,ccsd_doubles,3*nbuffs,tiles_in_buf_ccsd,needed_ccsd,&
                               & ccsd_pdm_j,ccsd_pdm_buff,j_buf_ccsd,req_ccsd)
   
          !FIND ij and ji in buffer
          ij = (j_tile-1)*dim_ts+i_tile; ji = (i_tile-1)*dim_ts+j_tile
          call assoc_ptr_to_buf(ij,vvoo,6*nbuffs,tiles_in_buf_vvoo,needed_vvoo,&
                            & vvoo_pdm_ij,vvoo_pdm_buff,ijbuf,req_vvoo)
          call assoc_ptr_to_buf(ji,vvoo,6*nbuffs,tiles_in_buf_vvoo,needed_vvoo,&
                            & vvoo_pdm_ji,vvoo_pdm_buff,jibuf,req_vvoo)

          call time_start_phase(PHASE_COMM)
   
          if( alloc_in_dummy )then
   
             call lsmpi_wait(req_vvvo(i_buf_vvvo))
             call lsmpi_wait(req_vvvo(j_buf_vvvo))
             call lsmpi_wait(req_ccsd(i_buf_ccsd))
             call lsmpi_wait(req_ccsd(j_buf_ccsd))
             call lsmpi_wait(req_vvoo(ijbuf))
             call lsmpi_wait(req_vvoo(jibuf))
   
          else
   
             if(vvvo%lock_set(i_tile)) call tensor_unlock_win(vvvo,i_tile)
             if(vvvo%lock_set(j_tile)) call tensor_unlock_win(vvvo,j_tile)
             if(ccsd_doubles%lock_set(i_tile)) call tensor_unlock_win(ccsd_doubles,i_tile)
             if(ccsd_doubles%lock_set(j_tile)) call tensor_unlock_win(ccsd_doubles,j_tile)
             if(vvoo%lock_set(ij)) call tensor_unlock_win(vvoo,ij)
             if(vvoo%lock_set(ji)) call tensor_unlock_win(vvoo,ji)
   
          endif
   
          needed_vvvo(i_buf_vvvo) = .true.; needed_vvvo(j_buf_vvvo) = .true.
          needed_ccsd(i_buf_ccsd) = .true.; needed_ccsd(j_buf_ccsd) = .true.
          needed_vvoo(ijbuf) = .true.; needed_vvoo(jibuf) = .true.
        
          call time_start_phase(PHASE_WORK)

          do k_tile = 1,j_tile

             k_pos = (k_tile-1)*tile_size+1

             call get_tile_dim(nelms,vvvo,k_tile)
             tile_size_tmp_k = int(nelms/((i8*nvirt)*(i8*nvirt**2)))

             !FIND k in buffer
             call assoc_ptr_to_buf(k_tile,vvvo,3*nbuffs,tiles_in_buf_vvvo,needed_vvvo,&
                                  & vvvo_pdm_k,vvvo_pdm_buff,k_buf_vvvo,req_vvvo)
             call assoc_ptr_to_buf(k_tile,ccsd_doubles,3*nbuffs,tiles_in_buf_ccsd,needed_ccsd,&
                                  & ccsd_pdm_k,ccsd_pdm_buff,k_buf_ccsd,req_ccsd)

             !FIND ik, ki, jk, and kj in buffer
             ik = (k_tile-1)*dim_ts+i_tile; ki = (i_tile-1)*dim_ts+k_tile
             jk = (k_tile-1)*dim_ts+j_tile; kj = (j_tile-1)*dim_ts+k_tile
             call assoc_ptr_to_buf(ik,vvoo,6*nbuffs,tiles_in_buf_vvoo,needed_vvoo,&
                               & vvoo_pdm_ik,vvoo_pdm_buff,ikbuf,req_vvoo)
             call assoc_ptr_to_buf(ki,vvoo,6*nbuffs,tiles_in_buf_vvoo,needed_vvoo,&
                               & vvoo_pdm_ki,vvoo_pdm_buff,kibuf,req_vvoo)
             call assoc_ptr_to_buf(jk,vvoo,6*nbuffs,tiles_in_buf_vvoo,needed_vvoo,&
                               & vvoo_pdm_jk,vvoo_pdm_buff,jkbuf,req_vvoo)
             call assoc_ptr_to_buf(kj,vvoo,6*nbuffs,tiles_in_buf_vvoo,needed_vvoo,&
                               & vvoo_pdm_kj,vvoo_pdm_buff,kjbuf,req_vvoo)

             call time_start_phase(PHASE_COMM)
   
             if( alloc_in_dummy )then
   
                call lsmpi_wait(req_vvvo(k_buf_vvvo))
                call lsmpi_wait(req_ccsd(k_buf_ccsd))
                call lsmpi_wait(req_vvoo(ikbuf))
                call lsmpi_wait(req_vvoo(kibuf))
                call lsmpi_wait(req_vvoo(jkbuf))
                call lsmpi_wait(req_vvoo(kjbuf))
 
             else
  
                if(vvvo%lock_set(k_tile)) call tensor_unlock_win(vvvo,k_tile)
                if(ccsd_doubles%lock_set(k_tile)) call tensor_unlock_win(ccsd_doubles,k_tile)   
                if(vvoo%lock_set(ik)) call tensor_unlock_win(vvoo,ik)
                if(vvoo%lock_set(ki)) call tensor_unlock_win(vvoo,ki)
                if(vvoo%lock_set(jk)) call tensor_unlock_win(vvoo,jk)
                if(vvoo%lock_set(kj)) call tensor_unlock_win(vvoo,kj)

             endif

             needed_vvvo(k_buf_vvvo) = .true.
             needed_ccsd(k_buf_ccsd) = .true.
             needed_vvoo(ikbuf) = .true.; needed_vvoo(kibuf) = .true.
             needed_vvoo(jkbuf) = .true.; needed_vvoo(kjbuf) = .true.
 
             call time_start_phase(PHASE_WORK)

             call time_start_phase(PHASE_WORK, twall = time_preload )
             call preload_tiles_in_bg_buf(vvvo,jobs,b_size,nvirt,nocc,i_tile,j_tile,k_tile,ij_count,3*nbuffs,&
                                         & needed_vvvo,tiles_in_buf_vvvo,vvvo_pdm_buff,req_vvvo,&
                                         & .true.,tile_size,dim_ts,dynamic_load)
             call preload_tiles_in_bg_buf(ccsd_doubles,jobs,b_size,nvirt,nocc,i_tile,j_tile,k_tile,ij_count,3*nbuffs,&
                                         & needed_ccsd,tiles_in_buf_ccsd,ccsd_pdm_buff,req_ccsd,&
                                         & .true.,tile_size,dim_ts,dynamic_load)
             call preload_tiles_in_bg_buf(vvoo,jobs,b_size,nvirt,nocc,i_tile,j_tile,k_tile,ij_count,6*nbuffs,&
                                         & needed_vvoo,tiles_in_buf_vvoo,vvoo_pdm_buff,req_vvoo,&
                                         & .true.,tile_size,dim_ts,dynamic_load,vovo_array=.true.)
             call time_start_phase(PHASE_WORK, ttot = time_preload )
             time_preload_tot = time_preload_tot + time_preload

! ##########################

             do i = i_pos,i_pos+tile_size_tmp_i-1

                call ptr_aliasing_ijk_par(nvirt,nocc,i,j,k,i_count,j_count,k_count,&
                                 & tile_size_tmp_i,tile_size_tmp_j,tile_size_tmp_k,&
                                 & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                                 & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                 & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                 & ccsd_pdm_i,ccsd_pdm_j,ccsd_pdm_k,&
                                 & vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k,&
                                 & vvoo_pdm_ij,vvoo_pdm_ik,vvoo_pdm_ji,vvoo_pdm_jk,vvoo_pdm_ki,vvoo_pdm_kj,&
                                 & ovoo,async_id,num_ids,1)

!$acc enter data copyin(ccsdpt_doubles(:,:,:,i)) async(async_id(3)) if(.not. full_no_frags)

                i_count = i_count+1

                do j = j_pos,j_pos+tile_size_tmp_j-1
         
                   if (j .gt. i) then

                      j_count = 0
                      cycle

                   endif

                   call ptr_aliasing_ijk_par(nvirt,nocc,i,j,k,i_count,j_count,k_count,&
                                     & tile_size_tmp_i,tile_size_tmp_j,tile_size_tmp_k,&
                                     & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                                     & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                     & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                     & ccsd_pdm_i,ccsd_pdm_j,ccsd_pdm_k,&
                                     & vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k,&
                                     & vvoo_pdm_ij,vvoo_pdm_ik,vvoo_pdm_ji,vvoo_pdm_jk,vvoo_pdm_ki,vvoo_pdm_kj,&
                                     & ovoo,async_id,num_ids,2)

!$acc enter data copyin(ccsdpt_doubles(:,:,:,j)) async(async_id(3)) if((.not. full_no_frags) .and. (i .gt. j))

                   j_count = j_count+1

                   do k = k_pos,k_pos+tile_size_tmp_k-1

                      if ((k .gt. j) .or. (k .gt. i)) then

                         k_count = 0
                         cycle
   
                      endif

                      call ptr_aliasing_ijk_par(nvirt,nocc,i,j,k,i_count,j_count,k_count,&
                                        & tile_size_tmp_i,tile_size_tmp_j,tile_size_tmp_k,&
                                        & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                                        & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                        & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                        & ccsd_pdm_i,ccsd_pdm_j,ccsd_pdm_k,&
                                        & vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k,&
                                        & vvoo_pdm_ij,vvoo_pdm_ik,vvoo_pdm_ji,vvoo_pdm_jk,vvoo_pdm_ki,vvoo_pdm_kj,&
                                        & ovoo,async_id,num_ids,3)

                      ! select type of tuple
                      tuple_type = -1

                      if ((i .eq. j) .and. (j .eq. k)) then
         
                         ! i == j == k
                         ! this always gives zero contribution

                         k_count = 0
                         cycle

                      endif

                      if ((i .eq. j) .and. (j .gt. k)) then
         
                         ! i == j > k
                         tuple_type = 1

!$acc enter data copyin(ccsdpt_doubles(:,:,:,k)) async(async_id(3)) if(.not. full_no_frags)
         
                      else if ((i .gt. j) .and. (j .eq. k)) then
         
                         ! i > j == k
                         tuple_type = 2

                      else
         
                         ! i > j > k 
                         tuple_type = 3

!$acc enter data copyin(ccsdpt_doubles(:,:,:,k)) async(async_id(3)) if(.not. full_no_frags)
         
                      end if

                      k_count = k_count+1

                      if(DECinfo%ccsolverskip)then
                         call random_number(trip_ampl)
                         call random_number(trip_tmp)
                      else

                      ! generate tuple(s)
                      TypeOfTuple_par_ijk: select case(tuple_type)

                      case(1)

!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_ijk_case1(i,k,nocc,nvirt,ccsd_i,ccsd_k,&
                                                 & vvvo_i,vvvo_k,&
                                                 & ovoo_ij,ovoo_ik,ovoo_ki,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

                         if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki) async(async_id(1))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_ijk_case1(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                 & vvoo_ij,vvoo_ik,vvoo_ki,&
                                                 & e4,e5,tmp_res_e5,t1_ptr(:,i),t1_ptr(:,k),&
                                                 & async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ik,vvoo_ki) async(async_id(3))

                         else

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_k) async(async_id(1))

                            call time_start_phase(PHASE_WORK, twall = time_driv )

                            call trip_denom_ijk(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))           

                            call ccsdpt_driver_ijk_case1(i,k,nocc,nvirt,vvoo_ij,vvoo_ik,vvoo_ki,&
                                                 & ovoo_ij,ovoo_ik,ovoo_ki,&
                                                 & vvvo_i,vvvo_k,&
                                                 & pt_1(:,i),pt_1(:,k),&
                                                 & pt_2(:,:,:,i),pt_2(:,:,:,k),&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

!$acc wait(async_id(5)) async(async_id(1))
!$acc exit data delete(vvvo_k,ovoo_ik,ovoo_ki) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ik,vvoo_ki) copyout(ccsdpt_doubles(:,:,:,k)) async(async_id(3))

                         endif

                      case(2)

!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_ijk_case2(i,j,nocc,nvirt,ccsd_i,ccsd_j,&
                                                 & vvvo_i,vvvo_j,&
                                                 & ovoo_ij,ovoo_ji,ovoo_jk,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

                         if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ovoo_jk) async(async_id(1))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_ijk_case2(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                 & vvoo_ij,vvoo_ji,vvoo_jk,&
                                                 & e4,e5,tmp_res_e5,t1_ptr(:,i),t1_ptr(:,j),&
                                                 & async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_jk) async(async_id(3))

                         else

                            call time_start_phase(PHASE_WORK, twall = time_driv )

                            call trip_denom_ijk(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
            
!$acc wait(async_id(3),async_id(4)) async(async_id(5))           
 
                            call ccsdpt_driver_ijk_case2(i,j,nocc,nvirt,vvoo_ij,vvoo_ji,vvoo_jk,&
                                                 & ovoo_ij,ovoo_ji,ovoo_jk,&
                                                 & vvvo_i,vvvo_j,&
                                                 & pt_1(:,i),pt_1(:,j),&
                                                 & pt_2(:,:,:,i),pt_2(:,:,:,j),&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

!$acc wait(async_id(5)) async(async_id(1))
!$acc exit data delete(ovoo_jk) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_jk) async(async_id(3))

                         endif

                      case(3)

!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_ijk_case3(i,j,k,nocc,nvirt,&
                                                 & ccsd_i,ccsd_j,ccsd_k,&
                                                 & vvvo_i,vvvo_j,vvvo_k,&
                                                 & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

                         if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki,ovoo_jk,ovoo_kj) async(async_id(1))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_ijk_case3(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                 & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                                 & e4,e5,tmp_res_e5,t1_ptr(:,i),t1_ptr(:,j),t1_ptr(:,k),&
                                                 & async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ik,vvoo_ki,vvoo_jk,vvoo_kj) async(async_id(3))

                         else

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_k) async(async_id(1))

                            call time_start_phase(PHASE_WORK, twall = time_driv ) 

                            call trip_denom_ijk(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))
            
                            call ccsdpt_driver_ijk_case3(i,j,k,nocc,nvirt,&
                                                 & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                                 & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                                 & vvvo_i,vvvo_j,vvvo_k,&
                                                 & pt_1(:,i),pt_1(:,j),pt_1(:,k),&
                                                 & pt_2(:,:,:,i),pt_2(:,:,:,j),pt_2(:,:,:,k),&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

!$acc wait(async_id(5)) async(async_id(1))
!$acc exit data delete(vvvo_k,ovoo_ik,ovoo_ki,ovoo_jk,ovoo_kj) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ik,vvoo_ki,vvoo_jk,vvoo_kj) copyout(ccsdpt_doubles(:,:,:,k)) async(async_id(3))

                         endif

                      end select TypeOfTuple_par_ijk
                      endif

                      if (k_count .eq. tile_size_tmp_k) k_count = 0

                   end do ! end k loop 

                   if (j_count .eq. tile_size_tmp_j) j_count = 0

                   if (j .eq. i) then
         
                      if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))

                      else

!$acc wait(async_id(4),async_id(5)) async(async_id(1))

                      endif

!$acc exit data delete(ovoo_ij) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ij) async(async_id(3))         

                   else ! i .gt. j

                      if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))

                      else

!$acc wait(async_id(4),async_id(5)) async(async_id(1))

                      endif

!$acc exit data delete(ccsd_j,vvvo_j,ovoo_ij,ovoo_ji) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ij,vvoo_ji) async(async_id(3)) if(full_no_frags)
!$acc exit data delete(vvoo_ij,vvoo_ji) copyout(ccsdpt_doubles(:,:,:,j)) async(async_id(3)) if(.not. full_no_frags)
         
                   endif

                end do ! end j loop
          
                if (i_count .eq. tile_size_tmp_i) i_count = 0

                if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))

                else

!$acc wait(async_id(4),async_id(5)) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data copyout(ccsdpt_doubles(:,:,:,i)) async(async_id(3))

                endif

!$acc exit data delete(ccsd_i,vvvo_i) async(async_id(1))

             end do ! end i loop

! ##########################

          needed_vvvo(k_buf_vvvo) = .false.
          needed_ccsd(k_buf_ccsd) = .false.
          needed_vvoo(ikbuf) = .false.; needed_vvoo(kibuf) = .false.
          needed_vvoo(jkbuf) = .false.; needed_vvoo(kjbuf) = .false.

          end do ! end k_tile loop

       needed_vvvo(i_buf_vvvo) = .false.; needed_vvvo(j_buf_vvvo) = .false.
       needed_ccsd(i_buf_ccsd) = .false.; needed_ccsd(j_buf_ccsd) = .false.
       needed_vvoo(ijbuf) = .false.; needed_vvoo(jibuf) = .false.

    enddo ijrun_par

    ! release ij_array
    call mem_dealloc(ij_array)

    call time_start_phase(PHASE_WORK)

!$acc wait

!$acc exit data delete(trip_tmp,trip_ampl,eivalvirt,tmp_res_e5,t1_ptr)&
!$acc& copyout(e4,e5) if(full_no_frags)
!
!$acc exit data delete(trip_tmp,trip_ampl,eivalvirt)&
!$acc& copyout(ccsdpt_singles) if(.not. full_no_frags)

    !measure idle time after loop
    call time_start_phase(PHASE_IDLE)
    call lsmpi_barrier(infpar%lg_comm)
    call time_start_phase(PHASE_WORK)

    if (alloc_in_dummy) then
       call tensor_unlock_wins(vvvo,all_nodes=.true.)
       call tensor_unlock_wins(ccsd_doubles,all_nodes=.true.)
       call tensor_unlock_wins(vvoo,all_nodes=.true.)
    endif

#ifdef VAR_CUBLAS

    ! Destroy the CUBLAS context
    stat = cublasDestroy_v2 ( cublas_handle )

#endif

    ! release async handles array
    call mem_dealloc(async_id)

    if (full_no_frags) then

       call mem_dealloc(tmp_res_e5)

       call ptr_final_ijk_par(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                      & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                      & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1_ptr)

    else

       call ptr_final_ijk_par(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                      & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                      & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj)

       call ptr_final_ijk_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    endif

    ! release triples ampl structures
#ifndef VAR_REAL_SP
    if( use_bg_buf )then
       call mem_pseudo_dealloc(trip_tmp)
       call mem_pseudo_dealloc(trip_ampl)
    else
       call mem_dealloc(trip_tmp)
       call mem_dealloc(trip_ampl)
    endif
#else
    call mem_dealloc(trip_tmp)
    call mem_dealloc(trip_ampl)
#endif

    ! release preloading stuff
    if( use_bg_buf )then
       call mem_pseudo_dealloc(vvoo_pdm_buff)
       call mem_pseudo_dealloc(ccsd_pdm_buff)
       call mem_pseudo_dealloc(vvvo_pdm_buff)
    else
       call mem_dealloc(vvoo_pdm_buff)
       call mem_dealloc(ccsd_pdm_buff)
       call mem_dealloc(vvvo_pdm_buff)
    endif
    call mem_dealloc(needed_vvvo)
    call mem_dealloc(needed_ccsd)
    call mem_dealloc(needed_vvoo)
    call mem_dealloc(req_vvvo)
    call mem_dealloc(req_ccsd)
    call mem_dealloc(req_vvoo)
    call mem_dealloc(tiles_in_buf_vvvo)
    call mem_dealloc(tiles_in_buf_ccsd)
    call mem_dealloc(tiles_in_buf_vvoo)
    call mem_dealloc(jobs)

    call time_phases_get_diff(current_wt=phase_cntrs)
    call time_start_phase( PHASE_WORK, ttot = time_pt_ijk )

    ! timings
    if (DECinfo%PL .gt. 2) then

       ! minima
       time_pt_ijk_min    = time_pt_ijk
       time_preload_min   = time_preload_tot
       time_trip_min      = time_trip_tot
       time_efull_min     = time_efull_tot
       time_driv_min      = time_driv_tot

       call lsmpi_reduce_realk_min( time_pt_ijk_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_trip_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_efull_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_driv_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_preload_min   , infpar%master, infpar%lg_comm )

       ! maxima
       time_pt_ijk_max    = time_pt_ijk
       time_preload_max   = time_preload_tot
       time_trip_max      = time_trip_tot
       time_efull_max     = time_efull_tot
       time_driv_max      = time_driv_tot

       call lsmpi_reduce_realk_max( time_pt_ijk_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_trip_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_efull_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_driv_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_preload_max   , infpar%master, infpar%lg_comm )

       ! reductions
       call lsmpi_local_reduction( time_pt_ijk   , infpar%master )
       call lsmpi_local_reduction( time_trip_tot   , infpar%master )
       call lsmpi_local_reduction( time_efull_tot   , infpar%master )
       call lsmpi_local_reduction( time_driv_tot   , infpar%master )
       call lsmpi_local_reduction( time_preload_tot   , infpar%master )

       ! unlock, waiting, and flushing
       unlock_time   = time_lsmpi_win_unlock - unlock_time
       waiting_time  = time_lsmpi_wait       - waiting_time
       flushing_time = time_lsmpi_win_flush  - flushing_time
   
       ! minima
       unlock_time_min    = unlock_time
       waiting_time_min   = waiting_time
       flushing_time_min  = flushing_time

       call lsmpi_reduce_realk_min( unlock_time_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( waiting_time_min  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( flushing_time_min , infpar%master, infpar%lg_comm )

       ! maxima
       unlock_time_max    = unlock_time
       waiting_time_max   = waiting_time
       flushing_time_max  = flushing_time
   
       call lsmpi_reduce_realk_max( unlock_time_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( waiting_time_max  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( flushing_time_max , infpar%master, infpar%lg_comm )

       ! reductions 
       call lsmpi_local_reduction( unlock_time   , infpar%master )
       call lsmpi_local_reduction( waiting_time  , infpar%master )
       call lsmpi_local_reduction( flushing_time , infpar%master )

       ! work, communication, and idle times
       ! minima
       time_w_min = phase_cntrs( PHASE_WORK_IDX )
       time_c_min = phase_cntrs( PHASE_COMM_IDX )
       time_i_min = phase_cntrs( PHASE_IDLE_IDX )

       call lsmpi_reduce_realk_min( time_w_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_c_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_i_min , infpar%master, infpar%lg_comm )

       ! maxima
       time_w_max = phase_cntrs( PHASE_WORK_IDX )
       time_c_max = phase_cntrs( PHASE_COMM_IDX )
       time_i_max = phase_cntrs( PHASE_IDLE_IDX )
   
       call lsmpi_reduce_realk_max( time_w_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_c_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_i_max , infpar%master, infpar%lg_comm )

       ! reductions 
       call lsmpi_local_reduction(phase_cntrs,nphases,infpar%master)

       if (infpar%lg_mynum .eq. infpar%master) then

          write(*,'("CCSD(T)-ijk time_trip                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_trip_max  ,  time_trip_tot   /dble(nodtotal),time_trip_min  ,  time_trip_tot   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time_preload             ",g10.3,g10.3,g10.3,g10.3)')&
             & time_preload_max  ,  time_preload_tot   /dble(nodtotal),time_preload_min  ,  time_preload_tot   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time_efull               ",g10.3,g10.3,g10.3,g10.3)')&
             & time_efull_max  ,  time_efull_tot   /dble(nodtotal),time_efull_min  ,  time_efull_tot   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time_driv                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_driv_max  ,  time_driv_tot   /dble(nodtotal),time_driv_min  ,  time_driv_tot   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time in lsmpi_win_unlock ",g10.3,g10.3,g10.3,g10.3)')&
             & unlock_time_max, unlock_time/dble(nodtotal),unlock_time_min,          unlock_time   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time in lsmpi_wait       ",g10.3,g10.3,g10.3,g10.3)')&
             & waiting_time_max,   waiting_time  /dble(nodtotal),waiting_time_min,   waiting_time  / time_pt_ijk
          write(*,'("CCSD(T)-ijk time in lsmpi_win_flush  ",g10.3,g10.3,g10.3,g10.3)')&
             & flushing_time_max,  flushing_time /dble(nodtotal),flushing_time_min,  flushing_time / time_pt_ijk
          write(*,'("CCSD(T)-ijk time WORK                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_w_max,phase_cntrs(PHASE_WORK_IDX)/dble(nodtotal),time_w_min,phase_cntrs(PHASE_WORK_IDX)/time_pt_ijk
          write(*,'("CCSD(T)-ijk time COMM                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_c_max,phase_cntrs(PHASE_COMM_IDX)/dble(nodtotal),time_c_min,phase_cntrs(PHASE_COMM_IDX)/time_pt_ijk
          write(*,'("CCSD(T)-ijk time IDLE                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_i_max,phase_cntrs(PHASE_IDLE_IDX)/dble(nodtotal),time_i_min,phase_cntrs(PHASE_IDLE_IDX)/time_pt_ijk

       endif

    endif

    if(dynamic_load)then
#ifdef VAR_HAVE_MPI3
       call lsmpi_win_unlock_all(dyn_w)
#else
       call lsquit("ERROR(ijk_loop_par): dynamic load only implemented for MPI3",-1)
#endif
       call mem_dealloc( dyn_i, dyn_c )

       call lsmpi_win_free( dyn_w )
    endif

    if (infpar%lg_mynum .eq. infpar%master) call LSTIMER('IJK_LOOP_PAR',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine ijk_loop_par
#endif


  !> \brief: main ijk-loop (serial version)
  !> \author: Janus Juul Eriksen
  !> \date: january 2014
  subroutine ijk_loop_ser(nocc,nvirt,ovoo,vvoo,vvvo,ccsd_doubles,&
                        & eivalocc,eivalvirt,ccsdpt_singles,&
                        & ccsdpt_doubles,e4,e5,t1)

    implicit none

    !> nocc,nvirt
    integer, intent(in) :: nocc,nvirt
    !> 2-el integrals
    real(realk), dimension(nocc,nvirt,nocc,nocc), target :: ovoo ! integrals (AI|JK) in the order (J,A,I,K)
    real(realk), dimension(nvirt,nvirt,nocc,nocc), target :: vvoo ! integrals (AI|BJ) in the order (A,B,I,J)
    real(realk), dimension(nvirt,nvirt,nvirt,nocc), target :: vvvo ! integrals (AI|BC) in the order (C,B,A,I)
    !> ccsd doubles amplitudes
    real(realk), dimension(nvirt,nvirt,nocc,nocc), target :: ccsd_doubles
    !> triples amplitudes and 3d work array
    real(real_pt), pointer, dimension(:,:,:) :: trip_tmp,trip_ampl
    !> ccsd(t) intermediates
    real(realk), dimension(nvirt,nocc), optional :: ccsdpt_singles
    real(realk), dimension(nvirt,nvirt,nocc,nocc), optional :: ccsdpt_doubles
    real(realk), optional :: e4,e5
    real(realk), dimension(nvirt,nocc), target, optional :: t1
    logical :: full_no_frags
    real(real_pt), pointer, dimension(:) :: tmp_res_e5
    !> pointers
    real(real_pt), pointer, dimension(:,:,:) :: ccsd_i,ccsd_j,ccsd_k
    real(real_pt), pointer, dimension(:,:,:) :: vvvo_i,vvvo_j,vvvo_k
    real(real_pt), pointer, dimension(:,:) :: vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj
    real(real_pt), pointer, dimension(:,:) :: ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj
    real(real_pt), pointer, dimension(:,:) :: t1_ptr
    !> tmp pointers
    real(real_pt), pointer :: pt_1(:,:),pt_2(:,:,:,:)
    !> orbital energiesi
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt)
    !> loop integers
    integer :: i,j,k,tuple_type
    !> async handles
    integer :: num_ids,m
    logical :: use_bg_buf
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
    integer(kind=acc_device_kind) :: acc_device_type
#else
    integer, pointer, dimension(:) :: async_id
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat
    real(realk) :: tcpu,twall,norm

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    full_no_frags = .false.

    use_bg_buf    = mem_is_background_buf_init()

    if (present(e4) .and. present(e5)) full_no_frags = .true.

    if (full_no_frags) then

       call mem_alloc(tmp_res_e5,i8*nvirt)

       call ptr_init_ijk_ser(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                      & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                      & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1,t1_ptr)

    else

       call ptr_init_ijk_ser(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                      & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                      & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj)

       call ptr_init_ijk_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    endif

#ifndef VAR_REAL_SP
    if(use_bg_buf)then
       ! init triples tuples structure
       call mem_pseudo_alloc(trip_ampl,i8*nvirt,i8*nvirt,i8*nvirt)
       ! init 3d wrk array
       call mem_pseudo_alloc(trip_tmp,i8*nvirt,i8*nvirt,i8*nvirt)
    else
       ! init triples tuples structure
       call mem_alloc(trip_ampl,nvirt,nvirt,nvirt)
       ! init 3d wrk array
       call mem_alloc(trip_tmp,nvirt,nvirt,nvirt)
    endif
#else
    ! init triples tuples structure
    call mem_alloc(trip_ampl,nvirt,nvirt,nvirt)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nvirt,nvirt,nvirt)
#endif

    ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
    ! handle 1: ccsd_doubles and vvvo / ovoo integrals
    ! handle 3: vvoo integrals and ccsdpt_doubles intermediate
    ! handle 4: triples amplitudes
    ! handle 5: energy evaluation 
    num_ids = 5
    call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

    if (DECinfo%acc_sync) then
       async_id = acc_async_sync
    else
       do m = 1,num_ids
          async_id(m) = int(m,kind=acc_handle_kind)
       enddo
    endif

#else

    if (DECinfo%acc_sync) then
       async_id = 0
    else
       do m = 1,num_ids
          async_id(m) = -m
       enddo
    endif

#endif

#ifdef VAR_CUBLAS

    ! initialize the CUBLAS context
    stat = cublasCreate_v2(cublas_handle)

#endif

!$acc enter data create(trip_tmp,trip_ampl,tmp_res_e5)&
!$acc& copyin(eivalvirt,t1_ptr,e4,e5) if(full_no_frags)
!
!$acc enter data create(trip_tmp,trip_ampl)&
!$acc& copyin(eivalvirt,ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

 irun_ser: do i=2,nocc ! i == j == k == 1 gives zero contribution

              call ptr_aliasing_ijk_ser(nvirt,nocc,i,j,k,&
                             & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                             & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                             & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                             & ccsd_doubles,vvvo,vvoo,ovoo,async_id,num_ids,1)

!$acc enter data copyin(ccsdpt_doubles(:,:,:,i)) async(async_id(3)) if(.not. full_no_frags)

    jrun_ser: do j=1,i

                 call ptr_aliasing_ijk_ser(nvirt,nocc,i,j,k,&
                                & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                                & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                & ccsd_doubles,vvvo,vvoo,ovoo,async_id,num_ids,2)

!$acc enter data copyin(ccsdpt_doubles(:,:,:,j)) async(async_id(3)) if((.not. full_no_frags) .and. (i .gt. j))

       krun_ser: do k=1,j

                    call ptr_aliasing_ijk_ser(nvirt,nocc,i,j,k,&
                                   & ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                                   & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                   & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                   & ccsd_doubles,vvvo,vvoo,ovoo,async_id,num_ids,3)

                    ! select type of tuple
                    tuple_type = -1

                    if ((i .eq. j) .and. (j .eq. k)) then

                       ! i == j == k
                       ! this always gives zero contribution
                       cycle

                    else if ((i .eq. j) .and. (j .gt. k)) then

                       ! i == j > k
                       tuple_type = 1

!$acc enter data copyin(ccsdpt_doubles(:,:,:,k)) async(async_id(3)) if(.not. full_no_frags)

                    else if ((i .gt. j) .and. (j .eq. k)) then

                       ! i > j == k
                       tuple_type = 2

                    else

                       ! i > j > k 
                       tuple_type = 3

!$acc enter data copyin(ccsdpt_doubles(:,:,:,k)) async(async_id(3)) if(.not. full_no_frags)

                    end if

                    ! generate tuple(s)
                    TypeOfTuple_ser_ijk: select case(tuple_type)

                    case(1)

!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                       call trip_generator_ijk_case1(i,k,nocc,nvirt,ccsd_i,ccsd_k,&
                                               & vvvo_i,vvvo_k,&
                                               & ovoo_ij,ovoo_ik,ovoo_ki,&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

                       if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki) async(async_id(1))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_energy_full_ijk_case1(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                               & vvoo_ij,vvoo_ik,vvoo_ki,&
                                               & e4,e5,tmp_res_e5,t1_ptr(:,i),t1_ptr(:,k),&
                                               & async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ik,vvoo_ki) async(async_id(3))
 
                       else

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_k) async(async_id(1))

                          call trip_denom_ijk(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))
 
                          call ccsdpt_driver_ijk_case1(i,k,nocc,nvirt,vvoo_ij,vvoo_ik,vvoo_ki,&
                                               & ovoo_ij,ovoo_ik,ovoo_ki,&
                                               & vvvo_i,vvvo_k,&
                                               & pt_1(:,i),pt_1(:,k),&
                                               & pt_2(:,:,:,i),pt_2(:,:,:,k),&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(1))
!$acc exit data delete(vvvo_k,ovoo_ik,ovoo_ki) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ik,vvoo_ki) copyout(ccsdpt_doubles(:,:,:,k)) async(async_id(3))

                       endif

                    case(2)

!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                       call trip_generator_ijk_case2(i,j,nocc,nvirt,ccsd_i,ccsd_j,&
                                               & vvvo_i,vvvo_j,&
                                               & ovoo_ij,ovoo_ji,ovoo_jk,&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

                       if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ovoo_jk) async(async_id(1))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_energy_full_ijk_case2(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                               & vvoo_ij,vvoo_ji,vvoo_jk,&
                                               & e4,e5,tmp_res_e5,t1_ptr(:,i),t1_ptr(:,j),&
                                               & async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_jk) async(async_id(3))

                       else

                          call trip_denom_ijk(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_driver_ijk_case2(i,j,nocc,nvirt,vvoo_ij,vvoo_ji,vvoo_jk,&
                                               & ovoo_ij,ovoo_ji,ovoo_jk,&
                                               & vvvo_i,vvvo_j,&
                                               & pt_1(:,i),pt_1(:,j),&
                                               & pt_2(:,:,:,i),pt_2(:,:,:,j),&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(1))
!$acc exit data delete(ovoo_jk) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_jk) async(async_id(3))

                       endif

                    case(3)

!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                       call trip_generator_ijk_case3(i,j,k,nocc,nvirt,&
                                               & ccsd_i,ccsd_j,ccsd_k,&
                                               & vvvo_i,vvvo_j,vvvo_k,&
                                               & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

                       if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_k,vvvo_k,ovoo_ik,ovoo_ki,ovoo_jk,ovoo_kj) async(async_id(1))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_energy_full_ijk_case3(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                               & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                               & e4,e5,tmp_res_e5,t1_ptr(:,i),t1_ptr(:,j),t1_ptr(:,k),&
                                               & async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ik,vvoo_ki,vvoo_jk,vvoo_kj) async(async_id(3))

                       else

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_k) async(async_id(1))

                          call trip_denom_ijk(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_driver_ijk_case3(i,j,k,nocc,nvirt,&
                                               & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                                               & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,&
                                               & vvvo_i,vvvo_j,vvvo_k,&
                                               & pt_1(:,i),pt_1(:,j),pt_1(:,k),&
                                               & pt_2(:,:,:,i),pt_2(:,:,:,j),pt_2(:,:,:,k),&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(1))
!$acc exit data delete(vvvo_k,ovoo_ik,ovoo_ki,ovoo_jk,ovoo_kj) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ik,vvoo_ki,vvoo_jk,vvoo_kj) copyout(ccsdpt_doubles(:,:,:,k)) async(async_id(3))

                       endif

                    end select TypeOfTuple_ser_ijk

                 end do krun_ser

                 if (j .eq. i) then

                    if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))

                    else

!$acc wait(async_id(4),async_id(5)) async(async_id(1))

                    endif

!$acc exit data delete(ovoo_ij) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ij) async(async_id(3))

                 else ! i .gt. j

                    if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))

                    else

!$acc wait(async_id(4),async_id(5)) async(async_id(1))

                    endif

!$acc exit data delete(ccsd_j,vvvo_j,ovoo_ij,ovoo_ji) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo_ij,vvoo_ji) async(async_id(3)) if(full_no_frags)
!$acc exit data delete(vvoo_ij,vvoo_ji) copyout(ccsdpt_doubles(:,:,:,j)) async(async_id(3)) if(.not. full_no_frags)

                 end if
 
              end do jrun_ser

              if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(1))

              else

!$acc wait(async_id(4),async_id(5)) async(async_id(1))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data copyout(ccsdpt_doubles(:,:,:,i)) async(async_id(3))

              endif

!$acc exit data delete(ccsd_i,vvvo_i) async(async_id(1))

           end do irun_ser

!$acc wait

!$acc exit data delete(trip_tmp,trip_ampl,eivalvirt,tmp_res_e5,t1_ptr)&
!$acc& copyout(e4,e5) if(full_no_frags)
!
!$acc exit data delete(trip_tmp,trip_ampl,eivalvirt)&
!$acc& copyout(ccsdpt_singles) if(.not. full_no_frags)

#ifdef VAR_CUBLAS

    ! Destroy the CUBLAS context
    stat = cublasDestroy_v2(cublas_handle)

#endif

    ! release async handles array
    call mem_dealloc(async_id)

    if (full_no_frags) then

       call mem_dealloc(tmp_res_e5)

       call ptr_final_ijk_ser(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                      & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                      & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj,t1_ptr)

    else

       call ptr_final_ijk_ser(nvirt,nocc,ccsd_i,ccsd_j,ccsd_k,vvvo_i,vvvo_j,vvvo_k,&
                      & vvoo_ij,vvoo_ik,vvoo_ji,vvoo_jk,vvoo_ki,vvoo_kj,&
                      & ovoo_ij,ovoo_ik,ovoo_ji,ovoo_jk,ovoo_ki,ovoo_kj)

       call ptr_final_ijk_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    endif

#ifndef VAR_REAL_SP
    if( use_bg_buf )then
       call mem_pseudo_dealloc(trip_tmp)
       call mem_pseudo_dealloc(trip_ampl)
    else
       call mem_dealloc(trip_tmp)
       call mem_dealloc(trip_ampl)
    endif
#else
    call mem_dealloc(trip_tmp)
    call mem_dealloc(trip_ampl)
#endif

    call LSTIMER('IJK_LOOP_SER',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine ijk_loop_ser


#ifdef VAR_MPI
  !> \brief: main abc-loop (mpi version)
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine abc_loop_par(nocc,nvirt,ooov,oovv,vovv,ccsd_doubles,&
                        & eivalocc,eivalvirt,nodtotal,nbuffs,tile_size,ccsdpt_singles,&
                        & ccsdpt_doubles,e4,e5,t1)

    implicit none

    !> nocc,nvirt
    integer, intent(in) :: nocc,nvirt
    !> 2-el integrals
    real(realk), dimension(nocc,nocc,nocc,nvirt), target :: ooov ! integrals (AI|JK) in the order (K,I,J,A)
    type(tensor), intent(inout)  :: oovv ! integrals (AI|BJ) in the order (I,J,A,B)
    real(realk), pointer, dimension(:) :: oovv_pdm_ab,oovv_pdm_ba,oovv_pdm_ac,oovv_pdm_ca ! o^2*tile_size tiles from oovv
    real(realk), pointer, dimension(:) :: oovv_pdm_bc,oovv_pdm_cb ! o^2*tile_size tiles from oovv
    real(realk), pointer, dimension(:,:) :: oovv_pdm_buff      ! buffers to prefetch oovv tiles
    type(tensor), intent(inout)  :: vovv ! integrals (AI|BC) in the order (B,I,A,C)
    real(realk), pointer, dimension(:) :: vovv_pdm_a,vovv_pdm_b,vovv_pdm_c ! ov^2*tile_size tiles from vovv
    real(realk), pointer, dimension(:,:) :: vovv_pdm_buff      ! buffers to prefetch vovv tiles
    integer, intent(inout) :: nodtotal, tile_size
    !> ccsd doubles amplitudes
    type(tensor), intent(inout)  :: ccsd_doubles
    real(realk), pointer, dimension(:) :: ccsd_pdm_a,ccsd_pdm_b,ccsd_pdm_c ! vo^2*tile_size tiles from ccsd_doubles
    real(realk), pointer, dimension(:,:) :: ccsd_pdm_buff ! buffers to prefetch ccsd_doubles tiles
    !> triples amplitudes and 3d work array
    real(real_pt), pointer, dimension(:,:,:) :: trip_tmp, trip_ampl
    !> ccsd(t) intermediates
    real(realk), dimension(nocc,nvirt), optional :: ccsdpt_singles
    real(realk), dimension(nocc,nocc,nvirt,nvirt), optional :: ccsdpt_doubles
    real(realk),optional :: e4,e5
    real(realk), dimension(nocc,nvirt), target, optional :: t1
    logical :: full_no_frags
    real(real_pt), pointer, dimension(:) :: tmp_res_e5
    !> pointers
    real(real_pt), pointer, dimension(:) :: ccsd_a,ccsd_b,ccsd_c
    real(real_pt), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(real_pt), pointer, dimension(:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(real_pt), pointer, dimension(:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(real_pt), pointer, dimension(:,:) :: t1_ptr
    !> tmp pointers
    real(real_pt), pointer :: pt_1(:,:),pt_2(:,:,:,:)
    !> orbital energiesi
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt) 
    !> loop integers
    integer :: b_size,njobs,ab_comp,ab_count
    integer, pointer :: ab_array(:),jobs(:)
    integer :: a,b,c,tuple_type
    integer :: a_tile,b_tile,c_tile,a_pos,b_pos,c_pos,a_count,b_count,c_count
    integer :: a_buf_vovv,a_buf_ccsd,b_buf_vovv,b_buf_ccsd,c_buf_vovv,c_buf_ccsd
    integer :: abbuf,babuf,acbuf,cabuf,bcbuf,cbbuf
    integer :: ab,ba,ac,ca,bc,cb
    integer :: total_num_tiles_1,total_num_tiles_2,total_num_tiles,dim_ts
    integer :: nelms,tile_size_tmp_a,tile_size_tmp_b,tile_size_tmp_c
    !> preloading
    integer, intent(in) :: nbuffs
    integer,pointer, dimension(:) :: tiles_in_buf_vovv,tiles_in_buf_ccsd,tiles_in_buf_oovv
    integer(kind=ls_mpik), pointer, dimension(:) :: req_vovv,req_ccsd,req_oovv
    logical,pointer,dimension(:) :: needed_vovv,needed_ccsd,needed_oovv
    !> async handles
    integer :: num_ids,m
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
    integer(kind=acc_device_kind) :: acc_device_type
#else
    integer, pointer, dimension(:) :: async_id
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat
    ! timings
    real(realk) :: tcpu,twall,time_pt_abc,time_pt_abc_min,time_pt_abc_max
    real(realk) :: time_trip,time_efull,time_driv,time_preload
    real(realk) :: time_trip_min,time_efull_min,time_driv_min,time_preload_min
    real(realk) :: time_trip_max,time_efull_max,time_driv_max,time_preload_max
    real(realk) :: time_trip_tot,time_efull_tot,time_driv_tot,time_preload_tot
    real(realk) :: phase_cntrs(nphases)
    real(realk) :: flushing_time, flushing_time_min, flushing_time_max
    real(realk) :: unlock_time,   unlock_time_min,   unlock_time_max
    real(realk) :: waiting_time,  waiting_time_min,  waiting_time_max
    real(realk) :: time_w_min, time_w_max
    real(realk) :: time_c_min, time_c_max
    real(realk) :: time_i_min, time_i_max
    logical     :: use_bg_buf

    ! init timings
    unlock_time   = time_lsmpi_win_unlock
    waiting_time  = time_lsmpi_wait
    flushing_time = time_lsmpi_win_flush
    time_trip_tot = 0.0E0_realk; time_preload_tot = 0.0E0_realk; time_efull_tot = 0.0E0_realk; time_driv_tot = 0.0E0_realk
    call time_start_phase( PHASE_WORK, twall = time_pt_abc )
    call time_phases_get_current(current_wt=phase_cntrs)
    if (infpar%lg_mynum .eq. infpar%master) call LSTIMER('START',tcpu,twall,DECinfo%output)

    full_no_frags = .false.
    use_bg_buf    = mem_is_background_buf_init()

    if (present(e4) .and. present(e5)) full_no_frags = .true.

    if (full_no_frags) then

       call mem_alloc(tmp_res_e5,i8*nocc)

       call ptr_init_abc_par(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                      & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                      & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1,t1_ptr)

    else

       call ptr_init_abc_par(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                      & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                      & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb)

       call ptr_init_abc_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    endif

    call time_start_phase(PHASE_WORK)

    ! alloc and init stuff for preloading
    if( use_bg_buf )then
       call mem_pseudo_alloc(vovv_pdm_buff,int((i8*nocc)*(i8*nvirt**2)*tile_size,kind=8),int(i8*3*nbuffs,kind=8))
       call mem_pseudo_alloc(ccsd_pdm_buff,int((i8*nvirt)*nocc**2*tile_size,kind=8),int(i8*3*nbuffs,kind=8))
       call mem_pseudo_alloc(oovv_pdm_buff,int((i8*nocc**2)*tile_size**2,kind=8),int(i8*6*nbuffs,kind=8))
    else
       call mem_alloc(vovv_pdm_buff,nocc*nvirt**2*tile_size,3*nbuffs)
       call mem_alloc(ccsd_pdm_buff,nvirt*nocc**2*tile_size,3*nbuffs)
       call mem_alloc(oovv_pdm_buff,nocc**2*tile_size**2,6*nbuffs)
    endif
    call mem_alloc(needed_vovv,3*nbuffs)
    call mem_alloc(needed_ccsd,3*nbuffs)
    call mem_alloc(needed_oovv,6*nbuffs)
    call mem_alloc(tiles_in_buf_vovv,3*nbuffs)
    call mem_alloc(tiles_in_buf_ccsd,3*nbuffs)
    call mem_alloc(tiles_in_buf_oovv,6*nbuffs)
    call mem_alloc(req_vovv,3*nbuffs)
    call mem_alloc(req_ccsd,3*nbuffs)
    call mem_alloc(req_oovv,6*nbuffs)
    if (alloc_in_dummy) then
       call tensor_lock_wins(vovv,'s',all_nodes=.true.)
       call tensor_lock_wins(ccsd_doubles,'s',all_nodes=.true.)
       call tensor_lock_wins(oovv,'s',all_nodes=.true.)
    endif
    needed_vovv       = .false.
    needed_ccsd       = .false.
    needed_oovv       = .false.
    tiles_in_buf_vovv = -1
    tiles_in_buf_ccsd = -1
    tiles_in_buf_oovv = -1

#ifndef VAR_REAL_SP
    if(use_bg_buf)then
       ! init triples tuples structure
       call mem_pseudo_alloc(trip_ampl,i8*nvirt,i8*nocc,i8*nocc)
       ! init 3d wrk array
       call mem_pseudo_alloc(trip_tmp,i8*nvirt,i8*nocc,i8*nocc)
    else
       ! init triples tuples structure
       call mem_alloc(trip_ampl,nvirt,nocc,nocc)
       ! init 3d wrk array
       call mem_alloc(trip_tmp,nvirt,nocc,nocc)
    endif
#else
    ! init triples tuples structure
    call mem_alloc(trip_ampl,nvirt,nocc,nocc)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nvirt,nocc,nocc)
#endif

    ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
    ! handle 1: ccsd_doubles
    ! handle 2: vovv and ooov integrals
    ! handle 3: oovv integrals and ccsdpt_doubles intermediate
    ! handle 4: triples amplitudes
    ! handle 5: energy evaluation
    num_ids = 5
    call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

    if (DECinfo%acc_sync) then
       async_id = acc_async_sync
    else
       do m = 1,num_ids
          async_id(m) = int(m,kind=acc_handle_kind)
       enddo
    endif

#else

    if (DECinfo%acc_sync) then
       async_id = 0
    else
       do m = 1,num_ids
          async_id(m) = -m
       enddo
    endif

#endif

#ifdef VAR_CUBLAS

    ! initialize the CUBLAS context
    stat = cublasCreate_v2(cublas_handle)

#endif

    total_num_tiles_1 = vovv%ntiles
    total_num_tiles_2 = ccsd_doubles%ntiles
    if (total_num_tiles_1 .ne. total_num_tiles_2) call lsquit('total_num_tiles_1 .ne. total_num_tiles_2 (abc)',DECinfo%output) 
    total_num_tiles = total_num_tiles_1
    dim_ts = int(nvirt / tile_size)
    if (mod(nvirt,tile_size) .gt. 0) dim_ts = dim_ts + 1 

    a_count = 0
    b_count = 0
    c_count = 0

    tile_size_tmp_a = 0
    tile_size_tmp_b = 0
    tile_size_tmp_c = 0

    ! create job distribution list
    ! first, determine common batch size from number of tasks and nodes
    ! in the ab matrix, njobs is the number of elements in the lower triangular matrix
    ! always an even number [ n(n+1) is always an even number ]
    njobs = int((dim_ts**2 + dim_ts)/2)
    b_size = int(njobs/nodtotal)

    ! ab_array stores all jobs for composite ab indices in descending order
    call mem_alloc(ab_array,njobs)
    ! init list (one more than b_size since mod(njobs,nodtotal) is not necessearily zero
    call mem_alloc(jobs,b_size + 1)

    ! create ab_array
    call create_comp_array_ccsdpt(njobs,dim_ts,ab_array)
    ! fill the list
    call job_distrib_ccsdpt(b_size,njobs,ab_array,jobs)

    ! release ab_array
    call mem_dealloc(ab_array)

    ! now follows the main loop, which is collapsed like for the ijk scheme (cf. ijk_loop_par)

!!$acc enter data create(trip_tmp,trip_ampl,&
!!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c)&
!!$acc& copyin(eivalocc,ccsdpt_singles,e4) if(full_no_frags)
!!
!!$acc enter data create(trip_tmp,trip_ampl,&
!!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c)&
!!$acc& copyin(eivalocc,ccsdpt_singles) if(.not. full_no_frags)
!
!!$acc wait

 abrun_par: do ab_count = 1,b_size + 1

          ! get value of ab from job disttribution list
          ab_comp = jobs(ab_count)

          ! no more jobs to be done? otherwise leave the loop
          if (ab_comp .lt. 0) exit

          ! calculate a and b from composite ab value
          call calc_i_leq_j(ab_comp,dim_ts,a_tile,b_tile)

          a_pos = (a_tile-1)*tile_size+1
          b_pos = (b_tile-1)*tile_size+1

          call get_tile_dim(nelms,vovv,a_tile)
          tile_size_tmp_a = int(nelms/(nocc*nvirt**2))
          call get_tile_dim(nelms,vovv,b_tile)
          tile_size_tmp_b = int(nelms/(nocc*nvirt**2))

          !FIND a and b in buffer
          call assoc_ptr_to_buf(a_tile,vovv,3*nbuffs,tiles_in_buf_vovv,needed_vovv,&
                               & vovv_pdm_a,vovv_pdm_buff,a_buf_vovv,req_vovv)
          call assoc_ptr_to_buf(b_tile,vovv,3*nbuffs,tiles_in_buf_vovv,needed_vovv,&
                               & vovv_pdm_b,vovv_pdm_buff,b_buf_vovv,req_vovv)
          call assoc_ptr_to_buf(a_tile,ccsd_doubles,3*nbuffs,tiles_in_buf_ccsd,needed_ccsd,&
                               & ccsd_pdm_a,ccsd_pdm_buff,a_buf_ccsd,req_ccsd)
          call assoc_ptr_to_buf(b_tile,ccsd_doubles,3*nbuffs,tiles_in_buf_ccsd,needed_ccsd,&
                               & ccsd_pdm_b,ccsd_pdm_buff,b_buf_ccsd,req_ccsd)
   
          !FIND ab and ba in buffer
          ab = (b_tile-1)*dim_ts+a_tile; ba = (a_tile-1)*dim_ts+b_tile
          call assoc_ptr_to_buf(ab,oovv,6*nbuffs,tiles_in_buf_oovv,needed_oovv,&
                            & oovv_pdm_ab,oovv_pdm_buff,abbuf,req_oovv)
          call assoc_ptr_to_buf(ba,oovv,6*nbuffs,tiles_in_buf_oovv,needed_oovv,&
                            & oovv_pdm_ba,oovv_pdm_buff,babuf,req_oovv)
   
          call time_start_phase(PHASE_COMM)
   
          if( alloc_in_dummy )then
   
             call lsmpi_wait(req_vovv(a_buf_vovv))
             call lsmpi_wait(req_vovv(b_buf_vovv))
             call lsmpi_wait(req_ccsd(a_buf_ccsd))
             call lsmpi_wait(req_ccsd(b_buf_ccsd))
             call lsmpi_wait(req_oovv(abbuf))
             call lsmpi_wait(req_oovv(babuf))
   
          else
   
             if(vovv%lock_set(a_tile)) call tensor_unlock_win(vovv,a_tile)
             if(vovv%lock_set(b_tile)) call tensor_unlock_win(vovv,b_tile)
             if(ccsd_doubles%lock_set(a_tile)) call tensor_unlock_win(ccsd_doubles,a_tile)
             if(ccsd_doubles%lock_set(b_tile)) call tensor_unlock_win(ccsd_doubles,b_tile)
             if(oovv%lock_set(ab)) call tensor_unlock_win(oovv,ab)
             if(oovv%lock_set(ba)) call tensor_unlock_win(oovv,ba)
   
          endif
   
          needed_vovv(a_buf_vovv) = .true.; needed_vovv(b_buf_vovv) = .true.
          needed_ccsd(a_buf_ccsd) = .true.; needed_ccsd(b_buf_ccsd) = .true.
          needed_oovv(abbuf) = .true.; needed_oovv(babuf) = .true.
        
          call time_start_phase(PHASE_WORK)

!!$acc enter data copyin(vovv_pdm_a) async(async_id(2))
!!$acc enter data copyin(vovv_pdm_b) async(async_id(2)) if(b_tile .ne. a_tile) 

          do c_tile = 1,b_tile

             c_pos = (c_tile-1)*tile_size+1

             call get_tile_dim(nelms,vovv,c_tile)
             tile_size_tmp_c = int(nelms/(nocc*nvirt**2))

             !FIND c in buffer
             call assoc_ptr_to_buf(c_tile,vovv,3*nbuffs,tiles_in_buf_vovv,needed_vovv,&
                                  & vovv_pdm_c,vovv_pdm_buff,c_buf_vovv,req_vovv)
             call assoc_ptr_to_buf(c_tile,ccsd_doubles,3*nbuffs,tiles_in_buf_ccsd,needed_ccsd,&
                                  & ccsd_pdm_c,ccsd_pdm_buff,c_buf_ccsd,req_ccsd)

             !FIND ac, ca, bc, and cb in buffer
             ac = (c_tile-1)*dim_ts+a_tile; ca = (a_tile-1)*dim_ts+c_tile
             bc = (c_tile-1)*dim_ts+b_tile; cb = (b_tile-1)*dim_ts+c_tile
             call assoc_ptr_to_buf(ac,oovv,6*nbuffs,tiles_in_buf_oovv,needed_oovv,&
                               & oovv_pdm_ac,oovv_pdm_buff,acbuf,req_oovv)
             call assoc_ptr_to_buf(ca,oovv,6*nbuffs,tiles_in_buf_oovv,needed_oovv,&
                               & oovv_pdm_ca,oovv_pdm_buff,cabuf,req_oovv)
             call assoc_ptr_to_buf(bc,oovv,6*nbuffs,tiles_in_buf_oovv,needed_oovv,&
                               & oovv_pdm_bc,oovv_pdm_buff,bcbuf,req_oovv)
             call assoc_ptr_to_buf(cb,oovv,6*nbuffs,tiles_in_buf_oovv,needed_oovv,&
                               & oovv_pdm_cb,oovv_pdm_buff,cbbuf,req_oovv)

             call time_start_phase(PHASE_COMM)
   
             if( alloc_in_dummy )then
   
                call lsmpi_wait(req_vovv(c_buf_vovv))
                call lsmpi_wait(req_ccsd(c_buf_ccsd))
                call lsmpi_wait(req_oovv(acbuf))
                call lsmpi_wait(req_oovv(cabuf))
                call lsmpi_wait(req_oovv(bcbuf))
                call lsmpi_wait(req_oovv(cbbuf))
 
             else
  
                if(vovv%lock_set(c_tile)) call tensor_unlock_win(vovv,c_tile)
                if(ccsd_doubles%lock_set(c_tile)) call tensor_unlock_win(ccsd_doubles,c_tile)   
                if(oovv%lock_set(ac)) call tensor_unlock_win(oovv,ac)
                if(oovv%lock_set(ca)) call tensor_unlock_win(oovv,ca)
                if(oovv%lock_set(bc)) call tensor_unlock_win(oovv,bc)
                if(oovv%lock_set(cb)) call tensor_unlock_win(oovv,cb)

             endif

             needed_vovv(c_buf_vovv) = .true.
             needed_ccsd(c_buf_ccsd) = .true.
             needed_oovv(acbuf) = .true.; needed_oovv(cabuf) = .true.
             needed_oovv(bcbuf) = .true.; needed_oovv(cbbuf) = .true.
 
             call time_start_phase(PHASE_WORK)

!!$acc enter data copyin(vovv_pdm_c) async(async_id(2)) if(c_tile .ne. b_tile)

             call time_start_phase(PHASE_WORK, twall = time_preload )
             call preload_tiles_in_bg_buf(vovv,jobs,b_size,nvirt,nocc,a_tile,b_tile,c_tile,ab_count,3*nbuffs,&
                                         & needed_vovv,tiles_in_buf_vovv,vovv_pdm_buff,req_vovv,&
                                         & .false.,tile_size,dim_ts,.false.)
             call preload_tiles_in_bg_buf(ccsd_doubles,jobs,b_size,nvirt,nocc,a_tile,b_tile,c_tile,ab_count,3*nbuffs,&
                                         & needed_ccsd,tiles_in_buf_ccsd,ccsd_pdm_buff,req_ccsd,&
                                         & .false.,tile_size,dim_ts,.false.)
             call preload_tiles_in_bg_buf(oovv,jobs,b_size,nvirt,nocc,a_tile,b_tile,c_tile,ab_count,6*nbuffs,&
                                         & needed_oovv,tiles_in_buf_oovv,oovv_pdm_buff,req_oovv,&
                                         & .false.,tile_size,dim_ts,.false.,vovo_array=.true.)
             call time_start_phase(PHASE_WORK, ttot = time_preload )
             time_preload_tot = time_preload_tot + time_preload

! ##########################

             do a = a_pos,a_pos+tile_size_tmp_a-1

                call ptr_aliasing_abc_par(nvirt,nocc,a,b,c,a_count,b_count,c_count,&
                              & tile_size_tmp_a,tile_size_tmp_b,tile_size_tmp_c,&
                              & ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                              & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                              & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                              & ccsd_pdm_a,ccsd_pdm_b,ccsd_pdm_c,ooov,&
                              & vovv_pdm_a,vovv_pdm_b,vovv_pdm_c,&
                              & oovv_pdm_ab,oovv_pdm_ba,oovv_pdm_ac,oovv_pdm_ca,oovv_pdm_bc,oovv_pdm_cb,async_id,num_ids,1)

                a_count = a_count+1

!!$acc enter data copyin(ccsd_doubles(:,:,:,a)) async(async_id(1))

!!$acc enter data copyin(ooov(:,:,:,a)) async(async_id(2))

                do b = b_pos,b_pos+tile_size_tmp_b-1

                   if (b .gt. a) then

                      b_count = 0
                      cycle

                   endif

                   call ptr_aliasing_abc_par(nvirt,nocc,a,b,c,a_count,b_count,c_count,&
                                 & tile_size_tmp_a,tile_size_tmp_b,tile_size_tmp_c,&
                                 & ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                                 & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                                 & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                                 & ccsd_pdm_a,ccsd_pdm_b,ccsd_pdm_c,ooov,&
                                 & vovv_pdm_a,vovv_pdm_b,vovv_pdm_c,&
                                 & oovv_pdm_ab,oovv_pdm_ba,oovv_pdm_ac,oovv_pdm_ca,oovv_pdm_bc,oovv_pdm_cb,async_id,num_ids,2)

                   b_count = b_count+1

!!$acc enter data copyin(ooov(:,:,:,b)) async(async_id(2))

!!$acc enter data copyin(oovv(:,:,a,b),oovv(:,:,b,a)) async(async_id(3)) if(full_no_frags)
!!
!!$acc enter data copyin(oovv(:,:,a,b),oovv(:,:,b,a),&
!!$acc& ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,b,a)) async(async_id(3)) if(.not. full_no_frags)

                   do c = c_pos,c_pos+tile_size_tmp_c-1

                      if ((c .gt. b) .or. (c .gt. a)) then

                         c_count = 0
                         cycle
   
                      endif

!!$acc enter data copyin(ccsd_doubles(:,:,:,c)) async(async_id(1))

!!$acc enter data copyin(ooov(:,:,:,c)) async(async_id(2))

!!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b)) async(async_id(3)) if(full_no_frags)
!!
!!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b),&
!!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a),&
!!$acc& ccsdpt_doubles(:,:,b,c),ccsdpt_doubles(:,:,c,b)) async(async_id(3)) if(.not. full_no_frags)

                      call ptr_aliasing_abc_par(nvirt,nocc,a,b,c,a_count,b_count,c_count,&
                                    & tile_size_tmp_a,tile_size_tmp_b,tile_size_tmp_c,&
                                    & ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                                    & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                                    & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                                    & ccsd_pdm_a,ccsd_pdm_b,ccsd_pdm_c,ooov,&
                                    & vovv_pdm_a,vovv_pdm_b,vovv_pdm_c,&
                                    & oovv_pdm_ab,oovv_pdm_ba,oovv_pdm_ac,oovv_pdm_ca,oovv_pdm_bc,oovv_pdm_cb,async_id,num_ids,3)

                      ! select type of tuple
                      tuple_type = -1

                      if ((a .eq. b) .and. (b .eq. c)) then
         
                         ! a == b == c
                         ! this always gives zero contribution

                         c_count = 0
                         cycle

                      endif

                      if ((a .eq. b) .and. (b .gt. c)) then
         
                         ! a == b > c
                         tuple_type = 1
         
                      else if ((a .gt. b) .and. (b .eq. c)) then
         
                         ! a > b == c
                         tuple_type = 2

                      else
         
                         ! a > b > c 
                         tuple_type = 3
         
                      end if

                      c_count = c_count+1

                      ! generate tuple(s)
                      TypeOfTuple_par_abc: select case(tuple_type)

                      case(1)

!!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_abc_case1(a,c,nocc,nvirt,ccsd_a,ccsd_c,&
                                                 & ooov_a,ooov_c,vovv_ab,vovv_ac,vovv_ca,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,c)) async(async_id(1))

                         if (full_no_frags) then

!!$acc wait(async_id(4)) async(async_id(2))
!!$acc exit data delete(ooov(:,:,:,c)) async(async_id(2))

!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_abc_case1(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                  & oovv_ab,oovv_ac,oovv_ca,&
                                                  & e4,e5,tmp_res_e5,t1_ptr(:,a),t1_ptr(:,c),&
                                                  & async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull


!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a)) async(async_id(3))
         
                         else

                            call time_start_phase(PHASE_WORK, twall = time_driv )

                            call trip_denom_abc(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))
            
                            call ccsdpt_driver_abc_case1(a,c,nocc,nvirt,&
                                                 & oovv_ab,oovv_ac,oovv_ca,&
                                                 & vovv_ab,vovv_ac,vovv_ca,&
                                                 & ooov_a,ooov_c,&
                                                 & pt_1(:,a),pt_1(:,c),&
                                                 & pt_2(:,:,:,a),pt_2(:,:,:,c),&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

!!$acc wait(async_id(5)) async(async_id(2))
!!$acc exit data delete(ooov(:,:,:,c)) async(async_id(2))

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a))&
!!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a)) async(async_id(3))

                         endif

                      case(2)

!!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_abc_case2(a,b,nocc,nvirt,ccsd_a,ccsd_b,&
                                                 & ooov_a,ooov_b,vovv_ab,vovv_ba,vovv_bc,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

                         if (full_no_frags) then

! this is different...
!!$acc wait(async_id(4)) async(async_id(2))

!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_abc_case2(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                  & oovv_ab,oovv_ba,oovv_bc,&
                                                  & e4,e5,tmp_res_e5,t1_ptr(:,a),t1_ptr(:,b),&
                                                  & async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,b,c)) async(async_id(3))

                         else

                            call time_start_phase(PHASE_WORK, twall = time_driv )

                            call trip_denom_abc(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
            
!!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))
            
                            call ccsdpt_driver_abc_case2(a,b,nocc,nvirt,&
                                                 & oovv_ab,oovv_ba,oovv_bc,&
                                                 & vovv_ab,vovv_ba,vovv_bc,&
                                                 & ooov_a,ooov_b,&
                                                 & pt_1(:,a),pt_1(:,b),&
                                                 & pt_2(:,:,:,a),pt_2(:,:,:,b),&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

! this is different...
!!$acc wait(async_id(5)) async(async_id(2))

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,b,c))&
!!$acc& copyout(ccsdpt_doubles(:,:,b,c)) async(async_id(3))

                         endif

                      case(3)

!!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_abc_case3(a,b,c,nocc,nvirt,ccsd_a,ccsd_b,ccsd_c,&
                                                 & ooov_a,ooov_b,ooov_c,&
                                                 & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,c)) async(async_id(1))

                         if (full_no_frags) then

!!$acc wait(async_id(4)) async(async_id(2))
!!$acc exit data delete(ooov(:,:,:,c)) async(async_id(2))

!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_abc_case3(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                  & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                                                  & e4,e5,tmp_res_e5,t1_ptr(:,a),t1_ptr(:,b),t1_ptr(:,c),&
                                                  & async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b)) async(async_id(3))

                         else

                            call time_start_phase(PHASE_WORK, twall = time_driv ) 

                            call trip_denom_abc(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))            
            
                            call ccsdpt_driver_abc_case3(a,b,c,nocc,nvirt,&
                                                 & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                                                 & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                                                 & ooov_a,ooov_b,ooov_c,&
                                                 & pt_1(:,a),pt_1(:,b),pt_1(:,c),&
                                                 & pt_2(:,:,:,a),pt_2(:,:,:,b),pt_2(:,:,:,c),&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

!!$acc wait(async_id(5)) async(async_id(2))
!!$acc exit data delete(ooov(:,:,:,c)) async(async_id(2))

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b))&
!!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a),&
!!$acc& ccsdpt_doubles(:,:,b,c),ccsdpt_doubles(:,:,c,b)) async(async_id(3))

                         endif

                      end select TypeOfTuple_par_abc

                      if (c_count .eq. tile_size_tmp_c) c_count = 0

                   end do ! end c loop 

                   if (b_count .eq. tile_size_tmp_b) b_count = 0

                   if (b .eq. a) then
         
                      if (full_no_frags) then

! this is different
!!$acc wait(async_id(4)) async(async_id(2))

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,b)) async(async_id(3))

                      else

! this is different
!!$acc wait(async_id(5)) async(async_id(2))

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,b))&
!!$acc& copyout(ccsdpt_doubles(:,:,a,b)) async(async_id(3))

                      endif
         
                   else ! a .gt. b

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,b)) async(async_id(1))

                      if (full_no_frags) then

!!$acc wait(async_id(4)) async(async_id(2))
!!$acc exit data delete(ooov(:,:,:,b)) async(async_id(2))

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,b),oovv(:,:,b,a)) async(async_id(3))

                      else

!!$acc wait(async_id(5)) async(async_id(2))
!!$acc exit data delete(ooov(:,:,:,b)) async(async_id(2))

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,b),oovv(:,:,b,a))&
!!$acc& ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,b,a)) async(async_id(3))

                      endif
         
                   endif

                end do ! end b loop
          
                if (a_count .eq. tile_size_tmp_a) a_count = 0

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,a)) async(async_id(1))

                if (full_no_frags) then

!!$acc wait(async_id(4)) async(async_id(2))
!!$acc exit data delete(ooov(:,:,:,a)) async(async_id(2))

                else

!!$acc wait(async_id(5)) async(async_id(2))
!!$acc exit data delete(ooov(:,:,:,a)) async(async_id(2))

                endif

             end do ! end a loop

! ##########################

!!$acc exit data delete(vovv_pdm_c) async(async_id(2))

          needed_vovv(c_buf_vovv) = .false.
          needed_ccsd(c_buf_ccsd) = .false.
          needed_oovv(acbuf) = .false.; needed_oovv(cabuf) = .false.
          needed_oovv(bcbuf) = .false.; needed_oovv(cbbuf) = .false.

          end do ! end c_tile loop

!!$acc exit data delete(vovv_pdm_b) async(async_id(2))

       needed_vovv(a_buf_vovv) = .false.; needed_vovv(b_buf_vovv) = .false.
       needed_ccsd(a_buf_ccsd) = .false.; needed_ccsd(b_buf_ccsd) = .false.
       needed_oovv(abbuf) = .false.; needed_oovv(babuf) = .false.

!!$acc exit data delete(vovv_pdm_a) async(async_id(2))

       enddo abrun_par

    call time_start_phase(PHASE_WORK)

!!$acc wait
!
!!$acc exit data delete(trip_tmp,trip_ampl,&
!!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c,&
!!$acc& eivalocc)&
!!$acc& copyout(ccsdpt_singles,e4) if(full_no_frags)
!!
!!$acc exit data delete(trip_tmp,trip_ampl,&
!!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c,&
!!$acc& eivalocc) copyout(ccsdpt_singles) if(.not. full_no_frags)
!
!!$acc wait

    if (alloc_in_dummy) then
       call tensor_unlock_wins(vovv,all_nodes=.true.)
       call tensor_unlock_wins(ccsd_doubles,all_nodes=.true.)
       call tensor_unlock_wins(oovv,all_nodes=.true.)
    endif

#ifdef VAR_CUBLAS

    ! Destroy the CUBLAS context
    stat = cublasDestroy_v2 ( cublas_handle )

#endif

    ! release async handles array
    call mem_dealloc(async_id)

    if (full_no_frags) then

       call mem_dealloc(tmp_res_e5)

       call ptr_final_abc_par(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                      & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                      & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1_ptr)

    else

       call ptr_final_abc_par(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                      & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                      & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb)

       call ptr_final_abc_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    endif

#ifndef VAR_REAL_SP
    if( use_bg_buf )then
       call mem_pseudo_dealloc(trip_tmp)
       call mem_pseudo_dealloc(trip_ampl)
    else
       call mem_dealloc(trip_tmp)
       call mem_dealloc(trip_ampl)
    endif
#else
    call mem_dealloc(trip_tmp)
    call mem_dealloc(trip_ampl)
#endif

    ! release preloading stuff
    if( use_bg_buf )then
       call mem_pseudo_dealloc(oovv_pdm_buff)
       call mem_pseudo_dealloc(ccsd_pdm_buff)
       call mem_pseudo_dealloc(vovv_pdm_buff)
    else
       call mem_dealloc(oovv_pdm_buff)
       call mem_dealloc(ccsd_pdm_buff)
       call mem_dealloc(vovv_pdm_buff)
    endif
    call mem_dealloc(needed_vovv)
    call mem_dealloc(needed_ccsd)
    call mem_dealloc(needed_oovv)
    call mem_dealloc(req_vovv)
    call mem_dealloc(req_ccsd)
    call mem_dealloc(req_oovv)
    call mem_dealloc(tiles_in_buf_vovv)
    call mem_dealloc(tiles_in_buf_ccsd)
    call mem_dealloc(tiles_in_buf_oovv)
    call mem_dealloc(jobs)
    
    call time_phases_get_diff(current_wt=phase_cntrs)
    call time_start_phase( PHASE_WORK, ttot = time_pt_abc )

    ! timings
    if (DECinfo%PL .gt. 2) then

       ! minima
       time_pt_abc_min    = time_pt_abc
       time_preload_min   = time_preload_tot
       time_trip_min      = time_trip_tot
       time_efull_min     = time_efull_tot
       time_driv_min      = time_driv_tot

       call lsmpi_reduce_realk_min( time_pt_abc_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_trip_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_efull_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_driv_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_preload_min   , infpar%master, infpar%lg_comm )

       ! maxima
       time_pt_abc_max    = time_pt_abc
       time_preload_max   = time_preload_tot
       time_trip_max      = time_trip_tot
       time_efull_max     = time_efull_tot
       time_driv_max      = time_driv_tot

       call lsmpi_reduce_realk_max( time_pt_abc_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_trip_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_efull_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_driv_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_preload_max   , infpar%master, infpar%lg_comm )

       ! reductions
       call lsmpi_local_reduction( time_pt_abc   , infpar%master )
       call lsmpi_local_reduction( time_trip_tot   , infpar%master )
       call lsmpi_local_reduction( time_efull_tot   , infpar%master )
       call lsmpi_local_reduction( time_driv_tot   , infpar%master )
       call lsmpi_local_reduction( time_preload_tot   , infpar%master )

       ! unlock, waiting, and flushing
       unlock_time   = time_lsmpi_win_unlock - unlock_time
       waiting_time  = time_lsmpi_wait       - waiting_time
       flushing_time = time_lsmpi_win_flush  - flushing_time
   
       ! minima
       unlock_time_min    = unlock_time
       waiting_time_min   = waiting_time
       flushing_time_min  = flushing_time

       call lsmpi_reduce_realk_min( unlock_time_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( waiting_time_min  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( flushing_time_min , infpar%master, infpar%lg_comm )

       ! maxima
       unlock_time_max    = unlock_time
       waiting_time_max   = waiting_time
       flushing_time_max  = flushing_time
   
       call lsmpi_reduce_realk_max( unlock_time_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( waiting_time_max  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( flushing_time_max , infpar%master, infpar%lg_comm )

       ! reductions 
       call lsmpi_local_reduction( unlock_time   , infpar%master )
       call lsmpi_local_reduction( waiting_time  , infpar%master )
       call lsmpi_local_reduction( flushing_time , infpar%master )

       ! work, communication, and idle times
       ! minima
       time_w_min = phase_cntrs( PHASE_WORK_IDX )
       time_c_min = phase_cntrs( PHASE_COMM_IDX )
       time_i_min = phase_cntrs( PHASE_IDLE_IDX )

       call lsmpi_reduce_realk_min( time_w_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_c_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_i_min , infpar%master, infpar%lg_comm )

       ! maxima
       time_w_max = phase_cntrs( PHASE_WORK_IDX )
       time_c_max = phase_cntrs( PHASE_COMM_IDX )
       time_i_max = phase_cntrs( PHASE_IDLE_IDX )
   
       call lsmpi_reduce_realk_max( time_w_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_c_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_i_max , infpar%master, infpar%lg_comm )

       ! reductions 
       call lsmpi_local_reduction(phase_cntrs,nphases,infpar%master)

       if (infpar%lg_mynum .eq. infpar%master) then

          write(*,'("CCSD(T)-abc time_trip                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_trip_max  ,  time_trip_tot   /dble(nodtotal),time_trip_min  ,  time_trip_tot   / time_pt_abc
          write(*,'("CCSD(T)-abc time_preload             ",g10.3,g10.3,g10.3,g10.3)')&
             & time_preload_max  ,  time_preload_tot   /dble(nodtotal),time_preload_min  ,  time_preload_tot   / time_pt_abc
          write(*,'("CCSD(T)-abc time_efull               ",g10.3,g10.3,g10.3,g10.3)')&
             & time_efull_max  ,  time_efull_tot   /dble(nodtotal),time_efull_min  ,  time_efull_tot   / time_pt_abc
          write(*,'("CCSD(T)-abc time_driv                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_driv_max  ,  time_driv_tot   /dble(nodtotal),time_driv_min  ,  time_driv_tot   / time_pt_abc
          write(*,'("CCSD(T)-abc time in lsmpi_win_unlock ",g10.3,g10.3,g10.3,g10.3)')&
             & unlock_time_max, unlock_time/dble(nodtotal),unlock_time_min,          unlock_time   / time_pt_abc
          write(*,'("CCSD(T)-abc time in lsmpi_wait       ",g10.3,g10.3,g10.3,g10.3)')&
             & waiting_time_max,   waiting_time  /dble(nodtotal),waiting_time_min,   waiting_time  / time_pt_abc
          write(*,'("CCSD(T)-abc time in lsmpi_win_flush  ",g10.3,g10.3,g10.3,g10.3)')&
             & flushing_time_max,  flushing_time /dble(nodtotal),flushing_time_min,  flushing_time / time_pt_abc
          write(*,'("CCSD(T)-abc time WORK                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_w_max,phase_cntrs(PHASE_WORK_IDX)/dble(nodtotal),time_w_min,phase_cntrs(PHASE_WORK_IDX)/time_pt_abc
          write(*,'("CCSD(T)-abc time COMM                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_c_max,phase_cntrs(PHASE_COMM_IDX)/dble(nodtotal),time_c_min,phase_cntrs(PHASE_COMM_IDX)/time_pt_abc
          write(*,'("CCSD(T)-abc time IDLE                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_i_max,phase_cntrs(PHASE_IDLE_IDX)/dble(nodtotal),time_i_min,phase_cntrs(PHASE_IDLE_IDX)/time_pt_abc

       endif

    endif

    if (infpar%lg_mynum .eq. infpar%master) call LSTIMER('ABC_LOOP_PAR',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine abc_loop_par
#endif


  !> \brief: main abc-loop (serial version)
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine abc_loop_ser(nocc,nvirt,ooov,oovv,vovv,ccsd_doubles,&
                        & eivalocc,eivalvirt,ccsdpt_singles,&
                        & ccsdpt_doubles,e4,e5,t1)

    implicit none

    !> nocc,nvirt
    integer, intent(in) :: nocc,nvirt
    !> 2-el integrals
    real(realk), dimension(nocc,nocc,nocc,nvirt), target, intent(inout) :: ooov ! integrals (AI|JK) in the order (K,I,J,A)
    real(realk), dimension(nocc,nocc,nvirt,nvirt), target, intent(inout) :: oovv ! integrals (AI|BJ) in the order (I,J,A,B)
    real(realk), dimension(nvirt,nocc,nvirt,nvirt), target, intent(inout) :: vovv ! integrals (AI|BC) in the order (B,I,A,C)
    !> ccsd doubles amplitudes
    real(realk), dimension(nocc,nocc,nvirt,nvirt), target, intent(inout) :: ccsd_doubles
    !> triples amplitudes and 3d work array
    real(real_pt), pointer, dimension(:,:,:) :: trip_tmp,trip_ampl
    !> ccsd(t) intermediates
    real(realk), dimension(nocc,nvirt), optional :: ccsdpt_singles
    real(realk), dimension(nocc,nocc,nvirt,nvirt), optional :: ccsdpt_doubles
    real(realk),optional :: e4,e5
    real(realk), dimension(nocc,nvirt), target, optional :: t1
    logical :: full_no_frags
    !> pointers
    real(real_pt), pointer, dimension(:) :: tmp_res_e5
    real(real_pt), pointer, dimension(:,:,:) :: ccsd_a,ccsd_b,ccsd_c
    real(real_pt), pointer, dimension(:,:,:) :: ooov_a,ooov_b,ooov_c
    real(real_pt), pointer, dimension(:,:) :: oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb
    real(real_pt), pointer, dimension(:,:) :: vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb
    real(real_pt), pointer, dimension(:,:) :: t1_ptr
    !> tmp pointers
    real(real_pt), pointer :: pt_1(:,:),pt_2(:,:,:,:)
    !> orbital energiesi
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt)
    !> loop integers
    integer :: a,b,c,tuple_type
    !> async handles
    integer :: num_ids,m
    logical :: use_bg_buf
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
    integer(kind=acc_device_kind) :: acc_device_type
#else
    integer, pointer, dimension(:) :: async_id
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat
    real(realk) :: tcpu,twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    full_no_frags = .false.

    use_bg_buf = mem_is_background_buf_init()

    if (present(e4) .and. present(e5)) full_no_frags = .true.

    if (full_no_frags) then

       call mem_alloc(tmp_res_e5,i8*nocc)

       call ptr_init_abc_ser(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                      & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                      & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1,t1_ptr)

    else

       call ptr_init_abc_ser(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                      & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                      & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb)

       call ptr_init_abc_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    endif

#ifndef VAR_REAL_SP
    if(use_bg_buf)then
       ! init triples tuples structure
       call mem_pseudo_alloc(trip_ampl,i8*nvirt,i8*nocc,i8*nocc)
       ! init 3d wrk array
       call mem_pseudo_alloc(trip_tmp,i8*nvirt,i8*nocc,i8*nocc)
    else
       ! init triples tuples structure
       call mem_alloc(trip_ampl,nvirt,nocc,nocc)
       ! init 3d wrk array
       call mem_alloc(trip_tmp,nvirt,nocc,nocc)
    endif
#else
    ! init triples tuples structure
    call mem_alloc(trip_ampl,nvirt,nocc,nocc)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nvirt,nocc,nocc)
#endif

    ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
    ! handle 1: ccsd_doubles and vovv / ooov integrals
    ! handle 3: oovv integrals and ccsdpt_doubles intermediate
    ! handle 4: triples amplitudes
    ! handle 5: energy evaluation
    num_ids = 5
    call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

    if (DECinfo%acc_sync) then
       async_id = acc_async_sync
    else
       do m = 1,num_ids
          async_id(m) = int(m,kind=acc_handle_kind)
       enddo
    endif

#else

    if (DECinfo%acc_sync) then
       async_id = 0
    else
       do m = 1,num_ids
          async_id(m) = -m
       enddo
    endif

#endif

#ifdef VAR_CUBLAS

    ! initialize the CUBLAS context
    stat = cublasCreate_v2(cublas_handle)

#endif

!!$acc enter data create(trip_tmp,trip_ampl)&
!!$acc& copyin(eivalocc,ccsdpt_singles,e4) if(full_no_frags)
!!
!!$acc enter data create(trip_tmp,trip_ampl)&
!!$acc& copyin(eivalocc,ccsdpt_singles) if(.not. full_no_frags)
!
!!$acc wait

    do a = 2,nvirt ! a == b == c == 1 gives zero contribution

       call ptr_aliasing_abc_ser(nvirt,nocc,a,b,c,&
                       & ccsd_a,ccsd_b,ccsd_c,&
                       & ooov_a,ooov_b,ooov_c,&
                       & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                       & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                       & ccsd_doubles,ooov,vovv,oovv,async_id,num_ids,1)

!!$acc enter data copyin(ccsd_doubles(:,:,:,a),&
!!$acc& ooov(:,:,:,a)) async(async_id(1))
!
!!$acc enter data copyin(ccsdpt_doubles(:,:,:,a)) async(async_id(3)) if(.not. full_no_frags)

       do b=1,a

          call ptr_aliasing_abc_ser(nvirt,nocc,a,b,c,&
                          & ccsd_a,ccsd_b,ccsd_c,&
                          & ooov_a,ooov_b,ooov_c,&
                          & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                          & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                          & ccsd_doubles,ooov,vovv,oovv,async_id,num_ids,2)

!!$acc enter data copyin(ccsd_doubles(:,:,:,b),&
!!$acc& ooov(:,:,:,b),&
!!$acc& vovv(:,:,a,b),vovv(:,:,b,a)) async(async_id(1))
!
!!$acc enter data copyin(oovv(:,:,a,b),oovv(:,:,b,a)) async(async_id(3)) if(full_no_frags)
!!
!!$acc enter data copyin(oovv(:,:,a,b),oovv(:,:,b,a),&
!!$acc& ccsdpt_doubles(:,:,:,b)) async(async_id(3)) if(.not. full_no_frags)

          do c=1,b

             ! select type of tuple
             tuple_type = -1

             if ((a .eq. b) .and. (b .eq. c)) then

                ! a == b == c
                ! this always gives zero contribution
                cycle

             else if ((a .eq. b) .and. (b .gt. c)) then

                ! a == b > c
                tuple_type = 1

!!$acc enter data copyin(ccsd_doubles(:,:,:,c),&
!!$acc& ooov(:,:,:,c),&
!!$acc& vovv(:,:,a,c),vovv(:,:,c,a)) async(async_id(1))
!
!!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a)) async(async_id(3)) if(full_no_frags)
!!
!!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),&
!!$acc& ccsdpt_doubles(:,:,:,c)) async(async_id(3)) if(.not. full_no_frags)

             else if ((a .gt. b) .and. (b .eq. c)) then

                ! a > b == c
                tuple_type = 2

!!$acc enter data copyin(vovv(:,:,b,c)) async(async_id(1))
!
!!$acc enter data copyin(oovv(:,:,b,c)) async(async_id(3))

             else

                ! a > b > c 
                tuple_type = 3

!!$acc enter data copyin(ccsd_doubles(:,:,:,c),&
!!$acc& ooov(:,:,:,c),&
!!$acc& vovv(:,:,a,c),vovv(:,:,c,a),vovv(:,:,b,c),vovv(:,:,c,b)) async(async_id(1))
!
!!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b)) &
!!$acc& async(async_id(3)) if(full_no_frags)
!!
!!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b),&
!!$acc& ccsdpt_doubles(:,:,:,c)) async(async_id(3)) if(.not. full_no_frags)

             end if

             call ptr_aliasing_abc_ser(nvirt,nocc,a,b,c,&
                             & ccsd_a,ccsd_b,ccsd_c,&
                             & ooov_a,ooov_b,ooov_c,&
                             & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                             & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                             & ccsd_doubles,ooov,vovv,oovv,async_id,num_ids,3)

             ! generate tuple(s)
             TypeOfTuple_ser_abc: select case(tuple_type)

             case(1)

!!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                call trip_generator_abc_case1(a,c,nocc,nvirt,ccsd_a,ccsd_c,&
                                        & ooov_a,ooov_c,vovv_ab,vovv_ac,vovv_ca,&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,c),&
!!$acc& ooov(:,:,:,c),&
!!$acc& vovv(:,:,a,c),vovv(:,:,c,a)) async(async_id(1)) if(full_no_frags)
!!
!!$acc exit data delete(ccsd_doubles(:,:,:,c)) async(async_id(1)) if(.not. full_no_frags)

                if (full_no_frags) then

!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_energy_full_abc_case1(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                         & oovv_ab,oovv_ac,oovv_ca,&
                                         & e4,e5,tmp_res_e5,t1_ptr(:,a),t1_ptr(:,c),&
                                         & async_id,num_ids,cublas_handle)

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a)) async(async_id(3))

                else

                   call trip_denom_abc(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_driver_abc_case1(a,c,nocc,nvirt,&
                                        & oovv_ab,oovv_ac,oovv_ca,&
                                        & vovv_ab,vovv_ac,vovv_ca,&
                                        & ooov_a,ooov_c,&
                                        & pt_1(:,a),pt_1(:,c),&
                                        & pt_2(:,:,:,a),pt_2(:,:,:,c),&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!!$acc wait(async_id(5)) async(async_id(1))
!!$acc exit data delete(ooov(:,:,:,c),&
!!$acc& vovv(:,:,a,c),vovv(:,:,c,a)) async(async_id(1))
!
!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a))&
!!$acc& copyout(ccsdpt_doubles(:,:,:,c)) async(async_id(3))

                endif

             case(2)

!!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                call trip_generator_abc_case2(a,b,nocc,nvirt,ccsd_a,ccsd_b,&
                                        & ooov_a,ooov_b,vovv_ab,vovv_ba,vovv_bc,&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

                if (full_no_frags) then

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(vovv(:,:,b,c)) async(async_id(1))
!
!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_energy_full_abc_case2(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                         & oovv_ab,oovv_ba,oovv_bc,&
                                         & e4,e5,tmp_res_e5,t1_ptr(:,a),t1_ptr(:,b),&
                                         & async_id,num_ids,cublas_handle)

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,b,c)) async(async_id(3))

                else

                   call trip_denom_abc(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_driver_abc_case2(a,b,nocc,nvirt,&
                                        & oovv_ab,oovv_ba,oovv_bc,&
                                        & vovv_ab,vovv_ba,vovv_bc,&
                                        & ooov_a,ooov_b,&
                                        & pt_1(:,a),pt_1(:,b),&
                                        & pt_2(:,:,:,a),pt_2(:,:,:,b),&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!!$acc wait(async_id(5)) async(async_id(1))
!!$acc exit data delete(vovv(:,:,b,c)) async(async_id(1))
!
!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,b,c)) async(async_id(3))

                endif

             case(3)

!!$acc wait(async_id(1),async_id(5)) async(async_id(4))

                call trip_generator_abc_case3(a,b,c,nocc,nvirt,ccsd_a,ccsd_b,ccsd_c,&
                                        & ooov_a,ooov_b,ooov_c,&
                                        & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,c),&
!!$acc& ooov(:,:,:,c),&
!!$acc& vovv(:,:,a,c),vovv(:,:,c,a),vovv(:,:,b,c),vovv(:,:,c,b)) async(async_id(1)) if(full_no_frags)
!!
!!$acc exit data delete(ccsd_doubles(:,:,:,c)) async(async_id(1)) if(.not. full_no_frags)

                if (full_no_frags) then

!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_energy_full_abc_case3(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                         & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                                         & e4,e5,tmp_res_e5,t1_ptr(:,a),t1_ptr(:,b),t1_ptr(:,c),&
                                         & async_id,num_ids,cublas_handle)

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b)) async(async_id(3))

                else

                   call trip_denom_abc(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))

!!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_driver_abc_case3(a,b,c,nocc,nvirt,&
                                        & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                                        & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,&
                                        & ooov_a,ooov_b,ooov_c,&
                                        & pt_1(:,a),pt_1(:,b),pt_1(:,c),&
                                        & pt_2(:,:,:,a),pt_2(:,:,:,b),pt_2(:,:,:,c),&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!!$acc wait(async_id(5)) async(async_id(1))
!!$acc exit data delete(ooov(:,:,:,c),&
!!$acc& vovv(:,:,a,c),vovv(:,:,c,a),vovv(:,:,b,c),vovv(:,:,c,b)) async(async_id(1))
!
!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b))&
!!$acc& copyout(ccsdpt_doubles(:,:,:,c)) async(async_id(3))

                endif

             end select TypeOfTuple_ser_abc

          end do

          if (b .eq. a) then

             if (full_no_frags) then

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(vovv(:,:,a,b)) async(async_id(1))

             else

!!$acc wait(async_id(5)) async(async_id(1))
!!$acc exit data delete(vovv(:,:,a,b)) async(async_id(1))

             endif

!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,b)) async(async_id(3))

          else ! a .gt. b

             if (full_no_frags) then

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,b),&
!!$acc& ooov(:,:,:,b),&
!!$acc& vovv(:,:,a,b),vovv(:,:,b,a)) async(async_id(1))
!
!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,b),oovv(:,:,b,a)) async(async_id(3))

             else

!!$acc wait(async_id(4),async_id(5)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,b),&
!!$acc& ooov(:,:,:,b),&
!!$acc& vovv(:,:,a,b),vovv(:,:,b,a)) async(async_id(1))
!
!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data delete(oovv(:,:,a,b),oovv(:,:,b,a))&
!!$acc& copyout(ccsdpt_doubles(:,:,:,b)) async(async_id(3))

             endif

          endif

       end do
 
       if (full_no_frags) then

!!$acc wait(async_id(4)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,a),&
!!$acc& ooov(:,:,:,a)) async(async_id(1))

       else

!!$acc wait(async_id(4),async_id(5)) async(async_id(1))
!!$acc exit data delete(ccsd_doubles(:,:,:,a),&
!!$acc& ooov(:,:,:,a)) async(async_id(1))
!
!!$acc wait(async_id(5)) async(async_id(3))
!!$acc exit data copyout(ccsdpt_doubles(:,:,:,a)) async(async_id(3))

       endif

    end do 

!!$acc wait
!
!!$acc exit data delete(trip_tmp,trip_ampl,eivalocc)&
!!$acc& copyout(ccsdpt_singles,e4) if(full_no_frags)
!!
!!$acc exit data delete(trip_tmp,trip_ampl,eivalocc)&
!!$acc& copyout(ccsdpt_singles) if(.not. full_no_frags)
!
!!$acc wait

#ifdef VAR_CUBLAS

    ! Destroy the CUBLAS context
    stat = cublasDestroy_v2 ( cublas_handle )

#endif

    ! release async handles array
    call mem_dealloc(async_id)

    if (full_no_frags) then

       call mem_dealloc(tmp_res_e5)

       call ptr_final_abc_ser(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                      & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                      & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb,t1_ptr)

    else

       call ptr_final_abc_ser(nvirt,nocc,ccsd_a,ccsd_b,ccsd_c,ooov_a,ooov_b,ooov_c,&
                      & oovv_ab,oovv_ac,oovv_ba,oovv_bc,oovv_ca,oovv_cb,&
                      & vovv_ab,vovv_ac,vovv_ba,vovv_bc,vovv_ca,vovv_cb)

       call ptr_final_abc_pt(nvirt,nocc,ccsdpt_singles,ccsdpt_doubles,pt_1,pt_2)

    endif

#ifndef VAR_REAL_SP
    if( use_bg_buf )then
       call mem_pseudo_dealloc(trip_tmp)
       call mem_pseudo_dealloc(trip_ampl)
    else
       call mem_dealloc(trip_tmp)
       call mem_dealloc(trip_ampl)
    endif
#else
    call mem_dealloc(trip_tmp)
    call mem_dealloc(trip_ampl)
#endif

    call LSTIMER('ABC_LOOP_SER',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine abc_loop_ser

!endif mod_unreleased
#endif

  subroutine dummy_ccsdpt_kernels_routine()

  end subroutine dummy_ccsdpt_kernels_routine

end module ccsdpt_kernels_module
