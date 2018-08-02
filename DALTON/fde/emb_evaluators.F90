!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_evaluators

   use fde_types
   use fde_nadd_derv 
   use fde_evaluators_dalton

   public fde_calc_int_energy

   contains

      subroutine calc_fde_int_energy(energy,grid,gf_active,gf_frozen)
         type(grid_function), intent(in) :: gf_active
         type(grid_function), intent(in) :: gf_frozen
         type(fde_grid), intent(in)      :: grid
         real(kind=8),       intent(out) :: energy

         type(grid_function) :: gf_tot
         real(kind=8), allocatable :: derv_active(:,:)
         real(kind=8), allocatable :: derv_frozen(:,:)
         real(kind=8), allocatable :: derv_tot(:,:)
         real(kind=8) :: vna, vcb, w
         integer :: i

         energy = 0.0d0

         call new_grid_function(gf_tot,gf_active%npoints)

         gf_tot%n  = gf_active%n  + gf_frozen%n
         gf_tot%gn = gf_active%gn + gf_frozen%gn
         gf_tot%hn = gf_active%hn + gf_frozen%hn

         call fde_new_derv_result_partial(gf_active%npoints,derv_active)
         call fde_new_derv_result_partial(gf_frozen%npoints,derv_frozen)
         call fde_new_derv_result_partial(gf_tot%npoints,derv_tot)

         call calc_xcke_energy_density(gf_active,  derv_active)
         call calc_xcke_energy_density(gf_frozen,  derv_frozen)
         call calc_xcke_energy_density(gf_tot,derv_tot)

         ! get nuclear potential of active
         do i = 1, gf_active%npoints
            vna = gf_active%nucpot(i)
            vcb = gf_frozen%elpot(i)
            w   = grid%w(i)
            energy = energy                                                    &
                   + w*(vcb + derv_tot(i,0) - derv_active(i,0))*gf_active%n(i) &
                   + w*(vna + derv_tot(i,0) - derv_frozen(i,0))*gf_frozen%n(i)
         enddo

         call fde_del_derv_result_partial(derv_active)
         call fde_del_derv_result_partial(derv_frozen)
         call fde_del_derv_result_partial(derv_tot)

      end subroutine calc_fde_int_energy

end module fde_evaluators


