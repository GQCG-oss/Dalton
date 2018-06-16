!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_nadd_derv

      use fde_types

      use fde_max_block_length
      use fde_xcfun_interface

      public fde_initialize_nadd_functionals
      public fde_print_nadd_functionals

      public get_fde_nadd_fun_derv
      public fde_calc_xc_nadd_fun_derv
      public fde_calc_ke_nadd_fun_derv
      
      public fde_set_skip_nadd_xc
      public fde_set_skip_nadd_ke
      public fde_set_nadd_kef
      public fde_set_nadd_xcf
      public fde_set_nadd_all_alda

      public fde_xc_fun_is_lda
      public fde_xc_fun_is_gga
      public fde_ke_fun_is_lda
      public fde_ke_fun_is_gga
      public is_fdersp_alda

      logical,save :: skip_kin_contribution = .false.
      logical,save :: skip_xc_contribution  = .false.
      logical,save :: is_fdersp_alda        = .false.

      type(functional), save :: fde_kef,       fde_xcf
      type(functional), save :: fde_kef_alda,  fde_xcf_alda
      type(functional), save :: fde_kef_xalda, fde_xcf_xalda

! the default (orbital-free) kinetic energy functional for now is the thomas-fermi one
      character(len=80), save :: kefun_fde='kin_tf'
      character(len=80), save :: xcfun_fde='lda'
!      character(len=80), save :: kefun_fde='kin_pw91'
!      character(len=80), save :: xcfun_fde='pbe'

      character(len=80), save :: kefun_fde_alda='kin_tf'
      character(len=80), save :: xcfun_fde_alda='lda'

      logical, save :: fde_set_functionals(20) 
      logical, save :: functionals_initialized = .false.

      real(8) :: hfx_out, mu_out, beta_out

! derv_lenght here is set to the appropriate parameter from xc_derv
! looking at xcint, this would be done differently if using openrsp...
      integer :: derv_length = nr_nonzero_derv

      public fde_new_derv_result
      public fde_delete_derv_result
      
      interface fde_new_derv_result
         module procedure fde_new_derv_result_partial
      end interface

      interface fde_delete_derv_result
         module procedure fde_delete_derv_result_partial
      end interface

      contains


!   ----------------------------------------------------------------------------
      subroutine fde_new_derv_result_partial(npoints,derv)
!   ----------------------------------------------------------------------------
         real(kind=8), allocatable :: derv(:,:)
         integer :: npoints
         
         allocate(derv(npoints, 0:derv_length))
         derv = 0.0d0
      end subroutine fde_new_derv_result_partial


!   ----------------------------------------------------------------------------
      subroutine fde_del_derv_result_partial(derv)
!   ----------------------------------------------------------------------------
         real(kind=8), allocatable :: derv(:,:)
         
         deallocate(derv)
      end subroutine fde_del_derv_result_partial


!   ----------------------------------------------------------------------------
      subroutine calc_xcke_energy_density(gf,  derv)
!   ----------------------------------------------------------------------------
         type(grid_function), intent(in)  :: gf
         real(kind=8),        intent(out) :: derv(gf%npoints,0:derv_length)
         real(kind=8) :: dervtmp(max_block_length, 0:derv_length)

         integer :: i 
         integer :: order = 0
         integer :: bllen = 1
         real(kind=8) :: r_0(1)
         real(kind=8) :: z_0(1)
         real(kind=8) :: w(1)
         
         dervtmp = 0.0d0
         w(1) = 1.0d0
         
         do i = 1, gf%npoints
            r_0(1) = gf%n(i)
            z_0(1) = gf%gn(i,1)*gf%gn(i,1) &
                   + gf%gn(i,2)*gf%gn(i,2) &
                   + gf%gn(i,3)*gf%gn(i,3)
            
            call get_xcke_fun_derv(order,        &
                                   bllen,        &
                                   w,            &
                                   r_0,          &
                                   z_0,          &
                                   dervtmp)
                              
            derv(i,:) = dervtmp(1,:)
         enddo
      end subroutine calc_xcke_energy_density


!   ----------------------------------------------------------------------------
      subroutine fde_delete_derv_result_partial(derv)
!   ----------------------------------------------------------------------------
         real(kind=8), allocatable :: derv(:,:)
         
         deallocate(derv)
      end subroutine fde_delete_derv_result_partial


!   ----------------------------------------------------------------------------
      function fde_xc_fun_is_lda()
!   ----------------------------------------------------------------------------
      logical :: fde_xc_fun_is_lda

         fde_xc_fun_is_lda = (xc_get_type(fde_xcf%id) == 0)
      end function


!   ----------------------------------------------------------------------------
      function fde_xc_fun_is_gga()
!   ----------------------------------------------------------------------------
      logical :: fde_xc_fun_is_gga

         fde_xc_fun_is_gga = (xc_get_type(fde_xcf%id) == 1)
      end function
   

!   ----------------------------------------------------------------------------
      function fde_ke_fun_is_lda()
!   ----------------------------------------------------------------------------
      logical :: fde_ke_fun_is_lda

         fde_ke_fun_is_lda = (xc_get_type(fde_kef%id) == 0)
      end function


!   ----------------------------------------------------------------------------
      function fde_ke_fun_is_gga()
!   ----------------------------------------------------------------------------
      logical :: fde_ke_fun_is_gga

         fde_ke_fun_is_gga = (xc_get_type(fde_kef%id) == 1)
      end function
      
         
!   ----------------------------------------------------------------------------
         subroutine fde_set_skip_nadd_xc(value)
!   ----------------------------------------------------------------------------
            logical :: value
            skip_xc_contribution  = value
         end subroutine 
    

!   ----------------------------------------------------------------------------
         subroutine fde_set_skip_nadd_ke(value)
!   ----------------------------------------------------------------------------
            logical :: value
            skip_kin_contribution  = value
         end subroutine 


!   ----------------------------------------------------------------------------
         subroutine fde_initialize_nadd_functionals
!   ----------------------------------------------------------------------------
            integer :: i, order = 5
            logical :: parallel_fde

            parallel_fde = .false.
            
        if (.not.functionals_initialized) then
            functionals_initialized = .true.
            
            call parse_functional(kefun_fde,      fde_kef,      hfx_out, mu_out, beta_out, .false.)
            call parse_functional(xcfun_fde,      fde_xcf,      hfx_out, mu_out, beta_out, .false.)
            call parse_functional(kefun_fde_alda, fde_kef_alda, hfx_out, mu_out, beta_out, .false.)
            call parse_functional(xcfun_fde_alda, fde_xcf_alda, hfx_out, mu_out, beta_out, .false.)


!            ! we set them now so that we can verify whether or not they
            call xcfun_set_functional(fde_xcf     , order, parallel_fde) 
            call xcfun_set_functional(fde_kef,      order, parallel_fde)
            call xcfun_set_functional(fde_xcf_alda, order, parallel_fde) 
            call xcfun_set_functional(fde_kef_alda, order, parallel_fde)
        endif
        
         end subroutine fde_initialize_nadd_functionals


!   ----------------------------------------------------------------------------
         subroutine fde_set_nadd_kef(id)
!   ----------------------------------------------------------------------------
            character(len=80) :: id

            write(kefun_fde,'(A)') id
         end subroutine fde_set_nadd_kef


!   ----------------------------------------------------------------------------
         subroutine fde_set_nadd_xcf(id)
!   ----------------------------------------------------------------------------
            character(len=80) :: id

            write(xcfun_fde,'(A)') id 
         end subroutine fde_set_nadd_xcf


!   ----------------------------------------------------------------------------
         subroutine fde_set_nadd_all_alda
!   ----------------------------------------------------------------------------
             is_fdersp_alda = .true.
         end subroutine fde_set_nadd_all_alda 


!   ----------------------------------------------------------------------------
         subroutine fde_print_nadd_functionals(unit)
!   ----------------------------------------------------------------------------
            integer :: unit 

            write(unit,'(A,A20)') ' * Non-additive functional : (KE) ',kefun_fde
            write(unit,'(A,A20)') '                             (XC) ',xcfun_fde
            write(unit,'(A)') ' ' 
            write(unit,'(A,A20)') '   if .RSPLDA set will use : (XC) ',xcfun_fde_alda
            write(unit,'(A,A20)') '                             (KE) ',kefun_fde_alda
            write(unit,'(A)') ' '
         end subroutine fde_print_nadd_functionals


!   ----------------------------------------------------------------------------
         subroutine get_fde_nadd_fun_derv(order,        &
                                          block_length, &
                                          w,            &
                                          r_0,          &
                                          r_0_t,        &
                                          z_0,          &
                                          z_0_t,        &
                                          derv)
!   ----------------------------------------------------------------------------
         integer,          intent(in)  :: order
         integer,          intent(in)  :: block_length
         real(8),          intent(in)  :: w(*)
         real(8),          intent(in)  :: r_0(*)
         real(8),          intent(in)  :: z_0(*)
         real(8),          intent(in)  :: r_0_t(*)
         real(8),          intent(in)  :: z_0_t(*)
         real(8),          intent(inout) :: derv(max_block_length, 0:derv_length)

          derv      = 0.0d0

         if (.not.skip_xc_contribution) &
            call get_nadd_fun_derv(fde_xcf,       &
                                   fde_xcf_alda,  &
                                   fde_xcf_alda,  &
                                   order,         &
                                   block_length,  &
                                   w,             &
                                   r_0,           &
                                   r_0_t,         &
                                   z_0,           &
                                   z_0_t,         &
                                   derv)

         if (.not.skip_kin_contribution) &
            call get_nadd_fun_derv(fde_kef,       &
                                   fde_kef_alda,  &
                                   fde_kef_alda,  &
                                   order,         &
                                   block_length,  &
                                   w,             &
                                   r_0,           &
                                   r_0_t,         &
                                   z_0,           &
                                   z_0_t,         &
                                   derv)

         end subroutine get_fde_nadd_fun_derv


!   ----------------------------------------------------------------------------
         subroutine get_xc_nadd_fun_derv(order,        &
                                         block_length, &
                                         w,            &
                                         r_0,          &
                                         r_0_t,        &
                                         z_0,          &
                                         z_0_t,        &
                                         derv)
!   ----------------------------------------------------------------------------
            integer,          intent(in)  :: order
            integer,          intent(in)  :: block_length
            real(8),          intent(in)  :: w(*)
            real(8),          intent(in)  :: r_0(*)
            real(8),          intent(in)  :: z_0(*)
            real(8),          intent(in)  :: r_0_t(*)
            real(8),          intent(in)  :: z_0_t(*)
            real(8),          intent(inout) :: derv(max_block_length, 0:derv_length)

            call get_nadd_fun_derv(fde_xcf,       &
                                fde_xcf_alda,  &
                                fde_xcf_alda,  &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0,           &
                                r_0_t,         &
                                z_0,           &
                                z_0_t,         &
                                derv)

         end subroutine get_xc_nadd_fun_derv


!   ----------------------------------------------------------------------------
         subroutine get_ke_nadd_fun_derv(order,        &
                                      block_length, &
                                      w,            &
                                      r_0,          &
                                      r_0_t,        &
                                      z_0,          &
                                      z_0_t,        &
                                      derv)
!   ----------------------------------------------------------------------------
            integer, intent(in)  :: order
            integer, intent(in)  :: block_length
            real(8), intent(in)  :: w(*)
            real(8), intent(in)  :: r_0(*)
            real(8), intent(in)  :: z_0(*)
            real(8), intent(in)  :: r_0_t(*)
            real(8), intent(in)  :: z_0_t(*)
            real(8), intent(inout) :: derv(max_block_length, 0:derv_length)

            call get_nadd_fun_derv(fde_kef,       &
                                fde_kef_alda,  &
                                fde_kef_alda,  &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0,           &
                                r_0_t,         &
                                z_0,           &
                                z_0_t,         &
                                derv)

      end subroutine get_ke_nadd_fun_derv



!   ----------------------------------------------------------------------------
      subroutine get_nadd_fun_derv(f,            &
                                    f_alda,       &
                                    f_xalda,      &
                                    order,        &
                                    block_length, &
                                    w,            &
                                    r_0,          &
                                    r_0_t,        &
                                    z_0,          &
                                    z_0_t,        &
                                    derv)
!   ----------------------------------------------------------------------------
         type(functional)              :: f
         type(functional)              :: f_alda, f_xalda
         integer,          intent(in)  :: order
         integer,          intent(in)  :: block_length
         real(8),          intent(in)  :: w(*)
         real(8),          intent(in)  :: r_0(*)
         real(8),          intent(in)  :: z_0(*)
         real(8),          intent(in)  :: r_0_t(*)
         real(8),          intent(in)  :: z_0_t(*)
         real(8),          intent(inout) :: derv(max_block_length, 0:derv_length)
         real(8)                       :: dervtmp(max_block_length, 0:derv_length)
         integer :: i, j

         dervtmp = 0.0d0

         call get_functional_derv(f,          &
                               f_alda,        &
                               f_xalda,       &
                               order,         &
                               block_length,  &
                               w,             &
                               r_0_t,         &
                               z_0_t,         &
                               derv_length,   &
                               dervtmp,       &
                               alda_real=is_fdersp_alda, &
                               alda_imag=is_fdersp_alda)

         do i = 1, max_block_length
            do j = 0, derv_length
               derv(i,j) = derv(i,j) + dervtmp(i,j)
            enddo
         enddo

         call get_functional_derv(f,          &
                               f_alda,        &
                               f_xalda,       &
                               order,         &
                               block_length,  &
                               w,             &
                               r_0,           &
                               z_0,           &
                               derv_length,   &
                               dervtmp,       &
                               alda_real=is_fdersp_alda, &
                               alda_imag=is_fdersp_alda)

         do i = 1, max_block_length
            do j = 0, derv_length
               derv(i,j) = derv(i,j) - dervtmp(i,j)
            enddo
         enddo

      end subroutine get_nadd_fun_derv


!   ----------------------------------------------------------------------------
      subroutine get_xcke_fun_derv(order,        &
                                   block_length, &
                                   w,            &
                                   r_0,          &
                                   z_0,          &
                                   derv)
!   ----------------------------------------------------------------------------
         integer,          intent(in)  :: order
         integer,          intent(in)  :: block_length
         real(8),          intent(in)  :: w(*)
         real(8),          intent(in)  :: r_0(*)
         real(8),          intent(in)  :: z_0(*)
         real(8),        intent(inout) :: derv(max_block_length, 0:derv_length)
         integer                       :: i, j
         real(8)                       :: dervtmp(max_block_length, 0:derv_length)

         derv = 0.0d0
         dervtmp = 0.0d0

         if (.not.skip_xc_contribution) then 
         call get_functional_derv(fde_xcf,       &
                                  fde_xcf_alda,  &
                                  fde_xcf_alda,  &
                                  order,         &
                                  block_length,  &
                                  w,             &
                                  r_0,           &
                                  z_0,           &
                                  derv_length,   &
                                  dervtmp,       &
                                  alda_real=is_fdersp_alda, &
                                  alda_imag=is_fdersp_alda)

         do i = 1, max_block_length
            do j = 0, derv_length
               derv(i,j) = derv(i,j) + dervtmp(i,j)
            enddo
         enddo
         endif

         if (.not.skip_kin_contribution) then
         call get_functional_derv(fde_kef,       &
                                  fde_kef_alda,  &
                                  fde_kef_alda,  &
                                  order,         &
                                  block_length,  &
                                  w,             &
                                  r_0,           &
                                  z_0,           &
                                  derv_length,   &
                                  dervtmp,       &
                                  alda_real=is_fdersp_alda, &
                                  alda_imag=is_fdersp_alda)


         do i = 1, max_block_length
            do j = 0, derv_length
               derv(i,j) = derv(i,j) + dervtmp(i,j)
            enddo
         enddo
         endif
         
     end subroutine get_xcke_fun_derv

end module fde_nadd_derv

