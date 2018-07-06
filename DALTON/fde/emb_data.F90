!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_data
      
   use fde_types
   use fde_cfg
   use fde_io

   implicit none

   private
!
! routines
!
   public fde_cleanup_data
   public fde_test_frozen_density
   public fde_test_active_density
   public fde_test_export
   public fde_import_frozen
   public fde_import_gridout
   public fde_import_static
   public fde_initialize_import_data
!
! variables 
!
! first, those related to the static embedding potential
!  
   public fde_grid_sv
   public fde_static_vemb

   integer, save :: fde_static_vemb_grid_np = 0
   logical, save :: static_data_is_initialized = .false.

   type(fde_grid), save :: fde_grid_sv
   real(kind=8), save, allocatable, target :: fde_static_vemb(:)
!
! importing frozen density, so the active one is needed too
! the grid imported isn't necessarily that of the static 
! potential
!
   public gf_frozen
   public gf_active
   public fde_grid_im

   integer, save :: fde_import_grid_np = 0
   logical, save :: import_data_is_initialized = .false.

   type(grid_function), save :: gf_frozen
   type(grid_function), save :: gf_active
   type(fde_grid), save      :: fde_grid_im
!
! export, using a different grid
!
   public gf_export
   public fde_grid_ex

   logical, save :: export_data_is_initialized = .false.
   integer, save :: fde_export_grid_np = 0

   type(fde_grid), save      :: fde_grid_ex
   type(grid_function), save :: gf_export

   contains


!   ----------------------------------------------------------------------------
   subroutine fde_import_static
!   ----------------------------------------------------------------------------
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i
      real(kind=8), pointer :: fde_vtmp(:), gridp(:,:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)
      
      file_unit = ftmp%embpot%unit
      file_name = ftmp%embpot%name

      call fde_file_open(file_name,file_unit)
      call read_grid(file_unit,gridp,fde_vtmp)
      call fde_file_close(file_unit)

      fde_static_vemb_grid_np = size(gridp,2)

      call fde_initialize_static_data

      do i=1, fde_grid_sv%npoints
         fde_static_vemb(i) = fde_vtmp(i)
         fde_grid_sv%r(1,i) = gridp(1,i)
         fde_grid_sv%r(2,i) = gridp(2,i)
         fde_grid_sv%r(3,i) = gridp(3,i)
         fde_grid_sv%w(i)   = gridp(4,i)
      enddo

      deallocate(fde_vtmp)
      deallocate(gridp)
   end subroutine fde_import_static


   subroutine fde_import_gridout
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i
      real(kind=8), pointer :: gridp(:,:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)

      file_unit = ftmp%export%unit
      file_name = ftmp%export%name

      call fde_file_open(file_name,file_unit)
      call read_grid(file_unit,gridp)
      call fde_file_close(file_unit)

      fde_export_grid_np = size(gridp,2)
      print *," np export ",fde_export_grid_np," file "//file_name

      call fde_initialize_export_data

      do i=1, fde_grid_ex%npoints
         fde_grid_ex%r(1,i) = gridp(1,i)
         fde_grid_ex%r(2,i) = gridp(2,i)
         fde_grid_ex%r(3,i) = gridp(3,i)
         fde_grid_ex%w(i)   = gridp(4,i)
      enddo
      deallocate(gridp)
   end subroutine


!   ----------------------------------------------------------------------------
   subroutine fde_import_frozen
!   ----------------------------------------------------------------------------
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i
      real(kind=8), pointer :: fde_frozen_prp(:,:), gridp(:,:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)
      
      file_unit = ftmp%frozen%unit
      file_name = ftmp%frozen%name

      call fde_file_open(file_name,file_unit)
      call read_grid(file_unit,gridp,fde_frozen_prp)
      call fde_file_close(file_unit)

      fde_import_grid_np = size(gridp,2)

      call fde_initialize_import_data
      
      print *, "np #1: ",fde_import_grid_np,"  np #2: ",fde_grid_im%npoints
      do i = 1, fde_grid_im%npoints
         gf_frozen%elpot(i)  = fde_frozen_prp(1,i) + fde_frozen_prp(2,i)
!        gf_frozen%nucpot(i) = 0.0
         gf_frozen%n(i)      = fde_frozen_prp(3,i)
         gf_frozen%gn(1,i)   = fde_frozen_prp(4,i)
         gf_frozen%gn(2,i)   = fde_frozen_prp(5,i)
         gf_frozen%gn(3,i)   = fde_frozen_prp(6,i)
         
         fde_grid_im%r(1,i)  = gridp(1,i)
         fde_grid_im%r(2,i)  = gridp(2,i)
         fde_grid_im%r(3,i)  = gridp(3,i)
         fde_grid_im%w(i)    = gridp(4,i)
      enddo

      deallocate(fde_frozen_prp)
      deallocate(gridp)
   end subroutine fde_import_frozen


!   ----------------------------------------------------------------------------
   subroutine fde_initialize_export_data
!   ----------------------------------------------------------------------------
      if (.not.export_data_is_initialized) then
         export_data_is_initialized = .true.
         if (fde_export_grid_np > 0) then

            call new_fde_grid(fde_grid_ex,fde_export_grid_np)
            call new_grid_function(gf_export,fde_export_grid_np)
         else
            call fde_quit('there are no grid points for fde export!')
         endif
      endif

   end subroutine fde_initialize_export_data


!   ----------------------------------------------------------------------------
   subroutine fde_initialize_static_data
!   ----------------------------------------------------------------------------
      if (.not.static_data_is_initialized) then
         static_data_is_initialized = .true.

         if (fde_static_vemb_grid_np > 0) then

            call new_fde_grid(fde_grid_sv,fde_static_vemb_grid_np)

            if (.not.allocated(fde_static_vemb)) then
               allocate(fde_static_vemb(fde_static_vemb_grid_np))
               fde_static_vemb         = 0.0d0
            endif
         else
            call fde_quit('there are no grid points for vemb!')
         endif
      endif
  end subroutine fde_initialize_static_data

!   ----------------------------------------------------------------------------
   subroutine fde_initialize_import_data
!   ----------------------------------------------------------------------------
      if (.not.import_data_is_initialized) then
         import_data_is_initialized = .true.
         if (fde_import_grid_np > 0) then
            call new_fde_grid(fde_grid_im,fde_import_grid_np)
            call new_grid_function(gf_frozen,fde_import_grid_np)
            call new_grid_function(gf_active,fde_import_grid_np)
         else
            call fde_quit('there are no grid points for fde import !')
         endif
      endif
      
   end subroutine fde_initialize_import_data


!   ----------------------------------------------------------------------------
   subroutine fde_cleanup_data
!   ----------------------------------------------------------------------------
      if (static_data_is_initialized) then
         deallocate(fde_static_vemb)
         call del_fde_grid(fde_grid_sv)
      endif

      if (import_data_is_initialized) then
         call del_grid_function(gf_active)
         call del_grid_function(gf_frozen)
         call del_fde_grid(fde_grid_im)
      endif

      if (export_data_is_initialized) then
         call del_grid_function(gf_export)
         call del_fde_grid(fde_grid_ex)
      endif
   end subroutine fde_cleanup_data


!   ----------------------------------------------------------------------------
   subroutine fde_test_frozen_density
!   ----------------------------------------------------------------------------
      integer      :: i
      real(kind=8) :: frz_particle_nr = 0.0

      do i = 1, fde_grid_im%npoints
         frz_particle_nr = frz_particle_nr + fde_grid_im%w(i)*gf_frozen%n(i) 
      end do
      
      write (*,*) '  Frozen density integrates to ',frz_particle_nr,' electrons'
      write (*,*) '  Grid points processed ',fde_grid_im%npoints

      frz_particle_nr = 0.0

   end subroutine fde_test_frozen_density


!   ----------------------------------------------------------------------------
   subroutine fde_test_active_density
!   ----------------------------------------------------------------------------
      integer      :: i 
      real(kind=8) :: act_particle_nr 

      act_particle_nr = 0.0d0

      do i = 1, fde_grid_im%npoints
         act_particle_nr = act_particle_nr + fde_grid_im%w(i)*gf_active%n(i) 
      end do
      
      write (*,*) ' Active density integrates to ',act_particle_nr,' electrons'
      write (*,*) ' Grid points processed ',fde_grid_im%npoints

      act_particle_nr = 0.0

   end subroutine fde_test_active_density


!   ----------------------------------------------------------------------------
   subroutine fde_test_export
!   ----------------------------------------------------------------------------
      integer      :: i
      real(kind=8) :: exp_int_n 
      real(kind=8) :: exp_int_vc 
      real(kind=8) :: exp_int_vn 

      exp_int_n  = 0.0d0
      exp_int_vc = 0.0d0  
      exp_int_vn = 0.0d0  

      
      do i = 1, fde_grid_ex%npoints
         exp_int_n  = exp_int_n  +       fde_grid_ex%w(i)*gf_export%n(i) 
         exp_int_vc = exp_int_vc - 0.5d0*fde_grid_ex%w(i)*gf_export%n(i)*(gf_export%elpot(i)-gf_export%nucpot(i)) 
         exp_int_vn = exp_int_vn -       fde_grid_ex%w(i)*gf_export%n(i)*gf_export%nucpot(i) 
      end do
      
      write (*,*) ' '
      write (*,*) 'Exported quantities (over ',fde_grid_ex%npoints,'grid points) integrate to:'
      write (*,'(A,F16.6,A)') ' - electron density    : ',exp_int_n,' electrons'
      write (*,'(A,F16.6,A)') ' - Hartree      energy : ',exp_int_vc,' a.u.'
      write (*,'(A,F16.6,A)') ' - n.-el. attr. energy : ',exp_int_vn,' a.u.'

      exp_int_n  = 0.0d0
      exp_int_vc = 0.0d0  
      exp_int_vn = 0.0d0  

   end subroutine fde_test_export

   
end module fde_data
