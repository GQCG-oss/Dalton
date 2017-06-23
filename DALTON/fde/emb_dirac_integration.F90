
module fde_dirac_matrices_integration

#ifdef PRG_DIRAC
!embed modules
   use fde_types
   use fde_cfg
   use fde_io
   use fde_data
   use fde_nadd_derv
   use fde_xcfun_interface
   use fde_max_block_length

! dirac-specific modules
   use interface_ao
   use interface_mo
   use interface_file_io
   use dft_cfg
   use interface_grid
   use ao_eval
   use density_eval

   implicit none

   public fde_dirac_emb_matrices_via_integration

   private

   type fde_omega_prefactor
      real(8), allocatable :: n(:)
      real(8), allocatable :: gn(:, :)
      real(8), allocatable :: s(:, :)
      real(8), allocatable :: gs(:, :, :)
   end type

   integer, parameter :: bllen = 1 !later more, now slowdft
   integer            :: block_length

   save

   integer, parameter :: file_unit      = 6

   real(8), allocatable :: rx(:)
   real(8), allocatable :: ry(:)
   real(8), allocatable :: rz(:)
   real(8), allocatable :: rw(:)
   real(8), allocatable :: dmat_sorted(:, :)

   real(8) :: nr_electrons_integrated
   real(8) :: fde_energy
   real(8) :: fde_mat_energy

   logical :: use_gga_qr = .false. !quaternion real
   logical :: use_gga_qi = .false. !quaternion imaginary
   logical :: use_gga    = .false.

   logical :: parallel_xc = .false.

!             number of points in current batch
   integer :: nr_points_batch

!             number of points currently looped over on this processor
   integer :: nr_points_on_this_proc

!             total nr of points (sum of nr_points_batch)
   integer :: nr_points_total

!             total nr of points that were actually used
!             we skip points with very small density
!             typically smaller than nr_points_total
   integer :: nr_points_used

   integer :: max_ao_g_order
   integer :: max_ao_m_order
   integer :: max_fun_derv

   real(8), external :: second
   real(8)           :: time_integration, time_integration_start
   real(8)           :: time_ao,          time_ao_start
   real(8)           :: time_matrix_dist, time_matrix_dist_start
   real(8)           :: time_density,     time_density_start
   real(8)           :: time_derv,        time_derv_start

   integer                      :: mat_dim
   integer                      :: nz

   real(8), pointer             :: dmat_0(:, :)
   real(8), allocatable, target :: dmat_0_container(:, :)

   real(8), allocatable         :: tmat_1(:, :)
   real(8), allocatable         :: tmat_2(:, :)

   integer                      :: nr_dmat
   real(8), pointer             :: dmat(:, :, :, :)
   real(8), allocatable, target :: dmat_container(:, :, :, :)

   integer                      :: nr_fmat
   real(8), pointer             :: fmat(:, :, :, :)
   real(8), allocatable, target :: fmat_container(:, :, :, :)

   integer, pointer             :: dmat_pg_sym(:)
   integer, allocatable, target :: dmat_pg_sym_container(:)
   integer, pointer             :: dmat_ih_sym(:)
   integer, allocatable, target :: dmat_ih_sym_container(:)
   integer, pointer             :: fmat_pg_sym(:)
   integer, allocatable, target :: fmat_pg_sym_container(:)

   integer                      :: nr_atoms
   real(8), pointer             :: property_gradient(:)
   real(8), allocatable, target :: property_gradient_container(:)

   logical :: do_potential
   integer :: response_order_mo
   integer :: response_order_ao
   logical :: do_london_rhs_direct
   logical :: do_london_rhs_ro
   logical :: do_geo
   logical :: do_geo_0
   logical :: do_geo_1
   logical :: do_geo_2
   logical :: do_london_ks
   logical :: do_london_lr
   logical :: do_overlap_diagnostic
   logical :: do_fde_update_vemb
   logical :: do_fde_response

contains


!   ----------------------------------------------------------------------------
    subroutine fde_dirac_emb_matrices_via_integration(  &
                             fde_mat_dim,               &
                             fde_nz,                    &
                             fde_dmat_0,                &
                             fde_nr_dmat,               &
                             fde_nr_fmat,               &
                             fde_dmat,                  &
                             fde_fmat,                  &
                             fde_dmat_ih_sym,           &
                             fde_dmat_pg_sym,           &
                             fde_fmat_pg_sym,           &
                             fde_do_potential,          &
                             fde_response_order_mo,     &
                             fde_response_order_ao)

! response_order_mo/_ao:  0 - do nothing (no response)
!                         1 - ks lr
!                         2 - ks qr
!                         *_mo - old code, *_ao - openrsp code
!
! do_potential            = .true.  - do energy and xc potential
!
! contributions to differentiated kohn-sham matrix:
! do_london_rhs_direct       = .true.  - xc contributions that require mag dervs of ao-s, "direct" (if london)
! do_london_rhs_ro           = .true.  - xc contributions to reotrhonormalization terms (if london)

!     --------------------------------------------------------------------------
      integer,                   intent(in)    :: fde_mat_dim
      integer,                   intent(in)    :: fde_nz
      real(8),           target, intent(in)    :: fde_dmat_0(fde_mat_dim, fde_mat_dim)
      integer,                   intent(in)    :: fde_nr_dmat
      integer,                   intent(in)    :: fde_nr_fmat
      real(8), optional, target, intent(in)    :: fde_dmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_dmat)
!     real(8), optional, target, intent(inout) :: fde_fmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_fmat)
      real(8), optional, target                :: fde_fmat(fde_mat_dim, fde_mat_dim, fde_nz, fde_nr_fmat)
      integer, optional, target, intent(in)    :: fde_dmat_ih_sym(fde_nr_dmat)
      integer, optional, target, intent(in)    :: fde_dmat_pg_sym(fde_nr_dmat)
      integer, optional, target, intent(in)    :: fde_fmat_pg_sym(fde_nr_fmat)

      logical, optional,         intent(in)    :: fde_do_potential
      integer, optional,         intent(in)    :: fde_response_order_mo
      integer, optional,         intent(in)    :: fde_response_order_ao
!     --------------------------------------------------------------------------
      real(8), external :: ddot
      integer, parameter :: max_mat_number = 10
      integer, target    :: default_pg_sym(max_mat_number)
      integer, target    :: default_ih_sym(max_mat_number)
      real(8)             :: save_ddot1, save_ddot2, h
!     --------------------------------------------------------------------------
      integer, parameter :: max_response_order_mo = 5
      
!     start timer
      time_integration_start = second()
!     reset timer
      time_ao          = 0.0d0
      time_matrix_dist = 0.0d0
      time_density     = 0.0d0
      time_derv        = 0.0d0

      nr_electrons_integrated = 0.0d0
      fde_energy    = 0.0d0

         mat_dim = fde_mat_dim
         nz      = fde_nz
         nr_dmat = fde_nr_dmat
         nr_fmat = fde_nr_fmat

!     set defaults
      do_potential          = .false.
      response_order_mo     = 0
      response_order_ao     = 0

!        change defaults
         if (present(fde_do_potential))          do_potential          = fde_do_potential
         if (present(fde_response_order_mo))     response_order_mo     = fde_response_order_mo
         if (present(fde_response_order_ao))     response_order_ao     = fde_response_order_ao

!     nullify(tmat_1)
!     nullify(tmat_2)

      nullify(dmat_0)
         dmat_0 => fde_dmat_0

      nullify(dmat)
      if (nr_dmat > 0) then
            if (present(fde_dmat)) then
               dmat => fde_dmat
            end if
      end if

      nullify(fmat)
      if (nr_fmat > 0) then
            if (present(fde_fmat)) then
               fmat => fde_fmat
            end if
      end if

      default_pg_sym = 1
      default_ih_sym = 0

      nullify(fmat_pg_sym)
      nullify(dmat_pg_sym)
      nullify(dmat_ih_sym)

         if (present(fde_fmat_pg_sym)) then
            fmat_pg_sym => fde_fmat_pg_sym
         else
            fmat_pg_sym => default_pg_sym(1:nr_fmat)
         end if
         if (present(fde_dmat_pg_sym)) then
            dmat_pg_sym => fde_dmat_pg_sym
         else
            dmat_pg_sym => default_pg_sym(1:nr_dmat)
         end if
         if (present(fde_dmat_ih_sym)) then
            dmat_ih_sym => fde_dmat_ih_sym
         else
            dmat_ih_sym => default_ih_sym(1:nr_dmat)
         end if

      if (response_order_mo > max_response_order_mo) then
          stop 'response_order_mo > max_response_order_mo'
      end if

      use_gga_qr = (fun_is_gga(fde_xcf))
      use_gga_qi = (fun_is_gga(fde_kef))

      if (response_order_mo > 0) then
         if (dft_cfg_alda_hs) use_gga_qr = .false.
         if (dft_cfg_alda_ha) use_gga_qi = .false.
      end if
      
      if (use_gga_qr .or. use_gga_qi) then
         use_gga = .true.
      end if

if (.false.) then      
      write (6,*) 'debug, in integrate_fde: use_gga ',use_gga
      write (6,*) '       fun_is_lda(fxc)       ',fun_is_lda(fde_xcf)
      write (6,*) '       fun_is_gga(fxc)       ',fun_is_gga(fde_xcf)
      write (6,*) '       fun_is_lda(fke)       ',fun_is_lda(fde_kef)
      write (6,*) '       fun_is_gga(fke)       ',fun_is_gga(fde_kef)
      write (6,*) '       is fde alda ?         ',is_fdersp_alda
endif
      
!     in qr it is not so easy to separate h+ and h-
!     contributions since they may mix
!     therefore i quit here if the user wants qr
!     and not the full alda or xalda
!     it is really unfair to quit as late as here
!     this can be checked earlier
!     but this is anyway an exotic option
      if (response_order_mo > 1) then
         if (dft_cfg_alda_hs .and. (.not. dft_cfg_alda_ha)) then
            call fde_quit('qr and partial alda is not implemented')
         end if
         if (dft_cfg_alda_ha .and. (.not. dft_cfg_alda_hs)) then
            call fde_quit('qr and partial alda is not implemented')
         end if
      end if

      max_ao_g_order = 0
      if (use_gga) then
        max_ao_g_order = 1
      end if
     
      max_ao_m_order = 0
     
      call ao_eval_init(max_ao_g_order, &
                        max_ao_m_order, &
                        .false.)

!        initial expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot1 = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
         end if

!        get average density from open shells
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) &
            call gtdoav(mat_dim, nz, 1, dmat_0, 1)

#ifdef PRG_DIRAC
         call scale_density_matrices(2.0d0)
#endif

         call insert_half_phases()

#ifdef DEBUG_XC
      write(*, *) 'debug: XC potential matrix zeroed out prior to integration'
      fmat = 0.0d0
#endif

      call loop_over_batches()

#ifdef PRG_DIRAC
         call scale_density_matrices(0.5d0)
#endif

!        get average density from open shells
         if ((nr_mo_gerade_positive_active + nr_mo_ungerade_positive_active) > 0) &
            call gtdoav(mat_dim, nz, 1, dmat_0, 1)

!        symmetrize matrices
         if (dft_cfg_blocked .or.  use_gga ) then
            if (nr_dmat > 0) then
               call gga_sym(mat_dim, nz, nr_fmat, fmat, fmat_pg_sym, dmat_ih_sym)
            else
               call gga_sym(mat_dim, nz, 1, fmat, (/1/), (/1/))
            end if
         end if

         call insert_half_phases()

!        final expectation value of the two-electron matrix
         if (do_potential) then
            save_ddot2    = ddot(mat_dim*mat_dim, dmat_0, 1, fmat, 1)
            fde_mat_energy = save_ddot2 - save_ddot1
         end if

        call report()

      call nullify_and_release()

#ifdef DEBUG_XC
      write(*, *) 'debug: real part of the XC potential matrix'
      call prqmat(fde_fmat,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  mat_dim,        &
                  1,              &
                  (/1, 2, 3, 4/), &
                  file_unit)
      call matrix_to_file('fde_fmat', mat_dim*mat_dim, fde_fmat)
#endif

   end subroutine


!   ----------------------------------------------------------------------------
   subroutine nullify_and_release()
!   ----------------------------------------------------------------------------

      if (associated(dmat)) then
         nullify(dmat)
      end if
      if (associated(fmat)) then
         nullify(fmat)
      end if
      if (associated(dmat_0)) then
         nullify(dmat_0)
      end if
      if (associated(dmat_pg_sym)) then
         nullify(dmat_pg_sym)
      end if
      if (associated(dmat_ih_sym)) then
         nullify(dmat_ih_sym)
      end if
      if (associated(fmat_pg_sym)) then
         nullify(fmat_pg_sym)
      end if

   end subroutine


!   ----------------------------------------------------------------------------
   subroutine scale_density_matrices(f)
!   ----------------------------------------------------------------------------
      real(8) :: f

      call dscal(mat_dim*mat_dim, f, dmat_0, 1)
      if (nr_dmat > 0) then
        call dscal(mat_dim*mat_dim*nz*nr_dmat, f, dmat, 1)
      end if

   end subroutine


!   ----------------------------------------------------------------------------
   subroutine insert_half_phases()
!   ----------------------------------------------------------------------------
      integer :: iz, im, iq


      if (nz < 4) then
         if (nr_dmat > 0) then
            do iz = 1, nz
               do im = 1, nr_dmat
                  iq = pq_to_uq(iz, dmat_pg_sym(im) - 1)
                  call my_q2bphase('D', iq, 1, dmat(1, 1, iz, im))
               end do
            end do
         end if
         if (nr_fmat > 0) then
            do iz = 1, nz
               do im = 1, nr_fmat
                  iq = pq_to_uq(iz, fmat_pg_sym(im) - 1)
                  call my_q2bphase('F', iq, 1, fmat(1, 1, iz, im))
               end do
            end do
         end if
      end if

   end subroutine




! ------------------------------------------------------------------------------
   subroutine integrate()
! ------------------------------------------------------------------------------
!        unperturbed density
         real(8) :: n_0(max_block_length)
!        unperturbed density gradient
         real(8) :: gn_0_pointwise(3)

         real(8) :: gnn_0(max_block_length)
         real(8) :: n_0_t(max_block_length)
         real(8) :: gnn_0_t(max_block_length)
         real(8) :: w(max_block_length)
!        --------------------------------------------------------------------------
!     geometric business
         real(8) :: r_t1(max_block_length)
         real(8) :: r_t2(max_block_length)
!     --------------------------------------------------------------------------
!     fde and fde response
!      logical :: fde_resp_verbose =.false.
         real(8) :: z0tt(3)
!     --------------------------------------------------------------------------
         integer :: id, m
         integer :: ii, ji
         real(8) :: temp
         real(8) :: block_threshold
         real(8) :: rho13, rho_outer, rho_resp, f
         integer :: ipoint
         integer :: k
         real(8) :: VT(3)
         integer              :: derv_length, im, i
         real(8), allocatable :: derv_active(:, :)
         real(8), allocatable :: derv_total(:, :)

         logical :: print_data_at_point

         real(8) :: n_b(3)
         real(8) :: gnn_b(3)
         real(8) :: gn_b(3, 3)

         logical :: un_above_threshold
         logical :: us_above_threshold
         logical :: ugn_above_threshold
         logical :: ugs_above_threshold
         logical :: utn_above_threshold
         logical :: uts_above_threshold
         logical :: us_lao_above_threshold
         logical :: ugs_lao_above_threshold

         real(8), allocatable :: ao(:, :)
         real(8), allocatable :: buffer(:, :)

         type(fde_omega_prefactor) :: u


         allocate(u%n (      nr_fmat))
         allocate(u%gn(   3, nr_fmat))
         allocate(u%s (   3, nr_fmat))
         allocate(u%gs(3, 3, nr_fmat))


!   check what functional derivatives do we need
!   ============================================

         max_fun_derv = 0

         if (do_potential) then
            max_fun_derv = 1
         end if
         if (response_order_mo > 0) then
            max_fun_derv = response_order_mo + 1
         end if

         derv_length = nr_nonzero_derv
         allocate(derv_active(max_block_length, 0:derv_length))
         allocate(derv_total(max_block_length, 0:derv_length))

!  better safe than sorry
         derv_active = 0.0d0
         derv_total  = 0.0d0

!  set it to zero, otherwise it is undefined with lda
         gnn_0 = 0.0d0

         allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
         ao = 0.0d0
         allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
         buffer = 0.0d0

         print_data_at_point = .true.

!  here starts the expensive loop over points
!  do not evaluate anything inside this loop
!  that does not change from point to point
   
         do ipoint = 1, nr_points_on_this_proc
            
            time_ao_start = second()
            call get_ao(1,          &
                  rx(ipoint), &
                  ry(ipoint), &
                  rz(ipoint), &
                  ao,         &
                  buffer)
            time_ao = time_ao + second() - time_ao_start

            time_density_start = second()

            if (need_ao_order(1, 0)) then
!        density and density gradient
               call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat_0, buffer, ao)
               gnn_0(1) = gn_0_pointwise(1)*gn_0_pointwise(1) &
                        + gn_0_pointwise(2)*gn_0_pointwise(2) &
                        + gn_0_pointwise(3)*gn_0_pointwise(3)
            else
!        density
               call get_n(n_0(1), 0, mat_dim, dmat_0, buffer, ao)
            endif

            time_density = time_density + second() - time_density_start

            if (n_0(1) > dft_cfg_tinydens) then
               nr_points_used = nr_points_used + 1

               u%n  = 0.0d0
               u%gn = 0.0d0
               u%s  = 0.0d0
               u%gs = 0.0d0

               w(1) = rw(ipoint)

               gf_active%n(ipoint)    = n_0(1)
               gf_active%gn(1,ipoint) = gn_0_pointwise(1)
               gf_active%gn(2,ipoint) = gn_0_pointwise(2)
               gf_active%gn(3,ipoint) = gn_0_pointwise(3)

               if ( gf_frozen%n(ipoint) > 0.0d0 ) then
                  n_0_t(1) = gf_active%n(ipoint)    + gf_frozen%n(ipoint)
                  z0tt(1)  = gf_active%gn(1,ipoint) + gf_frozen%gn(1,ipoint)
                  z0tt(2)  = gf_active%gn(2,ipoint) + gf_frozen%gn(2,ipoint)
                  z0tt(3)  = gf_active%gn(3,ipoint) + gf_frozen%gn(3,ipoint)

                  gnn_0_t(1) = z0tt(1)*z0tt(1) + z0tt(2)*z0tt(2) + z0tt(3)*z0tt(3)

               else
                  n_0_t(1)   = n_0(1)
                  gnn_0_t(1) = gnn_0(1)
               endif
!       number of electrons
               nr_electrons_integrated = nr_electrons_integrated + w(1)*n_0(1)

               time_derv_start = second()

! aspg: this way we are forcing the same xc functional for the active and non-additive
!       contributions. this is not desired, but we have to overcome xcfun's limitation
!       on using only up to 5 functionals...
!
               call get_xcke_fun_derv(max_fun_derv, & 
                                      bllen,        &
                                      w,            &
                                      n_0_t,        &
                                      gnn_0_t,      &
                                      derv_total)

               call get_xcke_fun_derv(max_fun_derv, & 
                                      bllen,        &
                                      w,            &
                                      n_0,          &
                                      gnn_0,        &
                                      derv_active)

               time_derv = time_derv + second() - time_derv_start


               time_matrix_dist_start = second()

!       exchange-correlation potential
!       ==============================

         if (do_potential) then

!         energy density will not be calcuated here but instead via
!         call calc_fde_interaction_energy(energy,gf_active,gf_frozen)
!         outside of the integration loop  
!           fde_energy = fde_energy + derv(1, d0000000)

           u%n(1) = (derv_total(1, d1000000) - derv_active(1, d1000000)) + w(1)*gf_frozen%elpot(ipoint) 

           if (use_gga) then
!           this factor has nothing to do
!           with symmetrization later
!           it comes from the fact that derivative of z
!           yields two identical terms (one left, one right)
              u%gn(:, 1) = 2.0d0*( derv_total(1, d0010000)*z0tt - derv_active(1, d0010000)*gn_0_pointwise)
           end if

         end if !exchange-correlation potential


!       linear response
!       ===============

        if (response_order_mo == 1) then
           call lr(bllen,          &
                   mat_dim,        &
                   nz,             &
                   nr_dmat,        &
                   dmat,           &
                   dmat_pg_sym,    &
                   dmat_ih_sym,    &
                   .false.,        &
                   ao,             &
                   n_0,            &
                   n_0_t,          &
                   gn_0_pointwise, &
                   z0tt,           &
                   gnn_0,          &
                   gnn_0_t,        &
                   w,              &
                   use_gga_qr,     &
                   use_gga_qi,     &
                   buffer,         &
                   derv_length,    &
                   derv_active,    &
                   derv_total,     &
                   u)
        endif

!         distribute over fock matrices
!         =============================

          do im = 1, nr_fmat

             un_above_threshold = .false.
             if (dabs(u%n(im)) > tiny(0.0d0)) then
                un_above_threshold = .true.
             end if

             ugn_above_threshold = .false.
             if (maxval((/dabs(u%gn(1, im)), &
                          dabs(u%gn(2, im)), &
                          dabs(u%gn(3, im))/)) > tiny(0.0d0)) then
                ugn_above_threshold = .true.
             end if

             if (ugn_above_threshold) then
!               F_pq +=   u%n    \chi_p* \chi_q
!                     + 2 u%gn_i \chi_p* \nabla_i \chi_q
                call nabla_omega_real(fmat,                &
                                      im,                  &
                                      fmat_pg_sym(im) - 1, &
                                      mat_dim,             &
                                      nz,                  &
                                      u%n(im),             &
                                      u%gn(1, im),         &
                                      ao)
             else
                if (un_above_threshold) then
!                  F_pq += u%n \chi_p* \chi_q
                   call omega_real(fmat,                &
                                   im,                  &
                                   fmat_pg_sym(im) - 1, &
                                   mat_dim,             &
                                   nz,                  &
                                   u%n(im),             &
                                   ao)
                end if
             end if

! spin-density stuff 

             us_above_threshold = .false.
             if (maxval((/dabs(u%s(1, im)), &
                          dabs(u%s(2, im)), &
                          dabs(u%s(3, im))/)) > tiny(0.0d0)) then
                us_above_threshold = .true.
             end if

             ugs_above_threshold = .false.
             if (maxval((/dabs(u%gs(1, 1, im)), &
                          dabs(u%gs(1, 2, im)), &
                          dabs(u%gs(1, 3, im)), &
                          dabs(u%gs(2, 1, im)), &
                          dabs(u%gs(2, 2, im)), &
                          dabs(u%gs(2, 3, im)), &
                          dabs(u%gs(3, 1, im)), &
                          dabs(u%gs(3, 2, im)), &
                          dabs(u%gs(3, 3, im))/)) > tiny(0.0d0)) then
                ugs_above_threshold = .true.
             end if

             if (ugs_above_threshold) then
!               F_pq +=   u%s_j   \chi_p* \Sigma_j \chi_q
!                     + 2 u%gs_ij \chi_p* \Sigma_j \nabla_i \chi_q
                call nabla_omega_imag(fmat,                &
                                      im,                  &
                                      fmat_pg_sym(im) - 1, &
                                      mat_dim,             &
                                      nz,                  &
                                      u%s(1, im),          &
                                      u%gs(1, 1, im),      &
                                      ao)
             else
                if (us_above_threshold) then
!                  F_pq += u%s_j \chi_p* \Sigma_j \chi_q
                   call omega_imag(fmat,                &
                                   im,                  &
                                   fmat_pg_sym(im) - 1, &
                                   mat_dim,             &
                                   nz,                  &
                                   u%s(1, im),          &
                                   ao)
                end if
             end if
          end do

! end spin-density stuff
          time_matrix_dist = time_matrix_dist + second() - time_matrix_dist_start

      end if !if (n_0(1) > dft_cfg_small_density_threshold) then
    end do !do ipoint = 1, nr_points_on_this_proc

    deallocate(ao)
    deallocate(buffer)
    deallocate(derv_total)
    deallocate(derv_active)

    deallocate(u%n)
    deallocate(u%gn)

   end subroutine integrate

   
!   ----------------------------------------------------------------------------
      subroutine report()
!   ----------------------------------------------------------------------------
         real(8) :: error_nr_electrons, fde_nr_electrons = 0

         error_nr_electrons = nr_electrons_integrated - fde_nr_electrons

         write(file_unit, *)
         write(file_unit, *) 'FDE integration'
         write(file_unit, *)
         write(file_unit, '(3x, a, i9)')                            &
            'number of grid points                          = ', &
            nr_points_total
         write(file_unit, '(3x, a, i9)')                            &
            'number of grid points below density threshold  = ', &
            nr_points_total - nr_points_used

      write(file_unit, '(3x, a, f26.16)')                        &
            'number of electrons from numerical integration = ', &
            nr_electrons_integrated

      if (do_potential) then
         write(file_unit, '(3x, a, f26.16)')                        &
               'FDE "energy":               = ', &
               fde_mat_energy
      end if

      write(file_unit, *)

!     timing
      time_integration = second() - time_integration_start
      call timtxt('  time spent in FDE integration                  =', &
                  time_integration, file_unit)
      write(file_unit, *)

   end subroutine report


!   ----------------------------------------------------------------------------
   real(8) function fde_get_nr_electrons_integrated()
!   ----------------------------------------------------------------------------
      fde_get_nr_electrons_integrated = nr_electrons_integrated
   end function fde_get_nr_electrons_integrated


!   ----------------------------------------------------------------------------
      subroutine loop_over_batches()
!   ----------------------------------------------------------------------------
         real(8) :: dummy
         integer :: idummy
         integer :: ibatch
         integer :: start_batch, i

         ibatch          = 0
         nr_points_total = 0
         nr_points_used  = 0
         
         start_batch = 1
         nr_points_batch = fde_grid_im%npoints

!        loop over batches
         do
            if (nr_points_batch < 0) exit

            allocate(rx(nr_points_batch))
            allocate(ry(nr_points_batch))
            allocate(rz(nr_points_batch))
            allocate(rw(nr_points_batch))

!              copy data over to rx, ry, rz, rw, fde_froz_n, fde_froz_ng
            rx(start_batch:nr_points_batch) = fde_grid_im%r(1,start_batch:nr_points_batch)
            ry(start_batch:nr_points_batch) = fde_grid_im%r(2,start_batch:nr_points_batch)
            rz(start_batch:nr_points_batch) = fde_grid_im%r(3,start_batch:nr_points_batch)
            rw(start_batch:nr_points_batch) = fde_grid_im%w(start_batch:nr_points_batch)
     
!           call get_fde_grid_coord(start_batch,nr_points_batch,rx,ry, rz)
!           call get_fde_grid_weight(start_batch,nr_points_batch,rw)
                 
            ibatch = ibatch + 1
            nr_points_total = nr_points_total + nr_points_batch
            nr_points_on_this_proc = nr_points_batch

            call integrate()

            deallocate(rx)
            deallocate(ry)
            deallocate(rz)
            deallocate(rw)

!          don't try to read another batch
           goto 999

!     end loop over batches
         end do
 999     continue
      end subroutine loop_over_batches


!   ----------------------------------------------------------------------------
      subroutine lr(block_length,       &
                    mat_dim,            &
                    nz,                 &
                    nr_mat,             &
                    dmaterturbed,       &
                    isym,               &
                    ih,                 &
                    allow_nonhermitian, &
                    ao,                 &
                    n_0,                &
                    n_0_t,              &
                    gn_0_pointwise,     &
                    gn_0_t_pointwise,   &
                    gnn_0,              &
                    gnn_0_t,            &
                    w,                  &
                    is_gga_qr,          &
                    is_gga_qi,          &
                    buffer,             &
                    n,                  &
                    derv_active,        &
                    derv_total,         &
                    u)
!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(in)    :: dmaterturbed(*)
    integer, intent(in)    :: isym(nr_mat)
    integer, intent(in)    :: ih(nr_mat)
    logical, intent(in)    :: allow_nonhermitian
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: n_0(*)
    real(8), intent(in)    :: n_0_t(*)
    real(8), intent(in)    :: gn_0_pointwise(*)
    real(8), intent(in)    :: gn_0_t_pointwise(*)
    real(8), intent(in)    :: gnn_0(*)
    real(8), intent(in)    :: gnn_0_t(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga_qr
    logical, intent(in)    :: is_gga_qi
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv_active(max_block_length, 0:n)
    real(8), intent(in)    :: derv_total(max_block_length, 0:n)
!   ----------------------------------------------------------------------------
    logical                :: calculate_qr, calculate_qi
    logical, save          :: have_triplet_contrib(100) = .false.
    logical, save          :: have_triplet = .false.
    logical, save          :: print_triplet_contrib = .true.


    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)
    real(8)                ::  t_b, ts_b(3)

    real(8)                :: z_b
    real(8)                :: z_b_t
    real(8)                :: y_b(3)
    real(8)                :: y_b_t(3)

    integer                :: k_start, k, im, im_off, j
    type(fde_omega_prefactor) :: u
!   ----------------------------------------------------------------------------

    k_start = 1

    do im = 1, nr_mat

      calculate_qr = .false.
      calculate_qi = .false.

      if (ih(im) == +1) calculate_qr = .true.
      if (ih(im) == -1) calculate_qi = .true.

      if (ih(im) == 0 .and. allow_nonhermitian) then
        calculate_qr = .true.
        calculate_qi = .true.
      end if

      if (dft_cfg_no_sdft) then
        calculate_qi = .false.
      end if

      if (calculate_qi) have_triplet_contrib(im) = .true.

      im_off = mat_dim*mat_dim*nz*(im - 1)

!     hermitian contribution (qr = quaternion real)
      if (calculate_qr) then

         if (is_gga_qr) then
            if (allow_nonhermitian) then
               call get_gn_nonhermitian(r_b,                      &
                                        gr_b,                     &
                                        isym(im) - 1,             &
                                        mat_dim,                  &
                                        dmaterturbed(1 + im_off), &
                                        buffer,                   &
                                        ao)
            else
               call get_gn(r_b,                      &
                           gr_b,                     &
                           isym(im) - 1,             &
                           mat_dim,                  &
                           dmaterturbed(1 + im_off), &
                           buffer,                   &
                           ao)
            end if
         else
            call get_n(r_b,                        &
                       isym(im) - 1,               &
                       mat_dim,                    &
                       dmaterturbed(1 + im_off), &
                       buffer,                     &
                       ao)
         end if

        u%n(im) = u%n(im) + (derv_total(1, d2000000) - derv_active(1, d2000000))*r_b

        if (is_gga_qr) then
           z_b = 2.0d0*(gn_0_pointwise(1)*gr_b(1) &
                      + gn_0_pointwise(2)*gr_b(2) &
                      + gn_0_pointwise(3)*gr_b(3))

           z_b_t = 2.0d0*(gn_0_t_pointwise(1)*gr_b(1) &
                        + gn_0_t_pointwise(2)*gr_b(2) &
                        + gn_0_t_pointwise(3)*gr_b(3))

           u%n(im) = u%n(im) + (derv_total(1, d1010000)*z_b_t - derv_active(1, d1010000)*z_b)

           u%gn(1:3, im) = u%gn(1:3, im) + 2.0d0*gn_0_t_pointwise(1:3)*r_b*derv_total(1, d1010000)   &
                                         + 2.0d0*gn_0_t_pointwise(1:3)*z_b_t*derv_total(1, d0020000) &
                                         + 2.0d0*gr_b(1:3)*derv_total(1, d0010000)                   &
                                         - 2.0d0*gn_0_pointwise(1:3)*r_b*derv_active(1, d1010000)    &
                                         - 2.0d0*gn_0_pointwise(1:3)*z_b*derv_active(1, d0020000)    &
                                         - 2.0d0*gr_b(1:3)*derv_active(1, d0010000)
        end if

      end if !end of qr contribution

!     antihermitian contribution (qi = quaternion imaginary)
      if (calculate_qi) then


        if (is_gga_qi) then
          if (allow_nonhermitian) then
            call get_gs_nonhermitian(s_b,                       &
                                     gs_b,                       &
                                     isym(im) - 1,               &
                                     mat_dim,                    &
                                     dmaterturbed(1 + im_off), &
                                     buffer,                     &
                                     ao)
          else
            call get_gs(s_b,                       &
                        gs_b,                       &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        dmaterturbed(1 + im_off), &
                        buffer,                     &
                        ao)
          end if
        else
        call get_s(s_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmaterturbed(1 + im_off), &
                   buffer,                     &
                   ao)
        end if

        do k = k_start, 3
           u%s(k, im) = u%s(k, im) + (derv_total(1, d0200000) - derv_active(1, d0200000))*s_b(k)

           if (is_gga_qi) then
              y_b(k)  = gn_0_pointwise(1)*gs_b(1, k) &
                      + gn_0_pointwise(2)*gs_b(2, k) &
                      + gn_0_pointwise(3)*gs_b(3, k)

              y_b_t(k) = gn_0_t_pointwise(1)*gs_b(1, k) &
                       + gn_0_t_pointwise(2)*gs_b(2, k) &
                       + gn_0_t_pointwise(3)*gs_b(3, k)

              u%s(k, im) = u%s(k, im) + (derv_total(1, d0101000)*y_b_t(k) - derv_active(1, d0101000)*y_b(k))

              u%gs(1:3, k, im) = u%gs(1:3, k, im) +       gn_0_t_pointwise(1:3)*s_b(k)*derv_total(1, d0101000)   &
                                                  +       gn_0_t_pointwise(1:3)*y_b_t(k)*derv_total(1, d0002000) &
                                                  + 2.0d0*gs_b(1:3, k)*derv_total(1, d0000100)                   &
                                                  -       gn_0_pointwise(1:3)*s_b(k)*derv_active(1, d0101000)    &
                                                  -       gn_0_pointwise(1:3)*y_b(k)*derv_active(1, d0002000)    &
                                                  - 2.0d0*gs_b(1:3, k)*derv_active(1, d0000100)
           end if
        end do

      end if !end of qi contribution

    end do

! now for some output
    if (print_triplet_contrib) then
       do im = 1, nr_mat
          if (have_triplet_contrib(im)) have_triplet = .true.
       end do

       if (have_triplet) then
          write(*,*) 'Triplet contributions to FDE kernel included'
       else
          write(*,*) 'Triplet contributions to FDE kernel neglected'
       endif
       print_triplet_contrib = .false.
    endif
!

  end subroutine
#else
#endif

end module fde_dirac_matrices_integration
