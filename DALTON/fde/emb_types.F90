!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_types

   public

   type fde_grid
      integer :: id
      integer :: npoints
      real(kind=8), pointer :: r(:,:)
      real(kind=8), pointer :: w(:) 
   end type

   type grid_function 
      integer :: id
      integer :: npoints
      integer :: nspin
      real(kind=8), pointer :: n(:)
      real(kind=8), pointer :: gn(:,:)
      real(kind=8), pointer :: elpot(:)
      real(kind=8), pointer :: nucpot(:)
      real(kind=8), pointer :: hn(:,:)
   end type

   public new_fde_grid
   public del_fde_grid
   public set_fde_grid_size

   public new_grid_function
   public del_grid_function
   public set_grid_function_nspin

   contains


      subroutine new_grid_function(gf,size)
         type(grid_function) :: gf
         integer :: size

         gf%npoints = size
         gf%nspin   = 1
         allocate(gf%n(size))
         allocate(gf%gn(3,size))
         allocate(gf%hn(6,size))
         allocate(gf%elpot(size))
         allocate(gf%nucpot(size))
         
         gf%n  = 0.0d0
         gf%gn = 0.0d0
         gf%hn = 0.0d0
         gf%elpot = 0.0d0
         gf%nucpot = 0.0d0
      end subroutine new_grid_function


      subroutine del_grid_function(gf)
         type(grid_function) :: gf
         integer :: size

         gf%npoints = 0
         deallocate(gf%n)
         deallocate(gf%gn)
         deallocate(gf%hn)
         deallocate(gf%elpot)
         deallocate(gf%nucpot)
      end subroutine del_grid_function


      subroutine set_grid_function_nspin(gf,nspin)
         type(grid_function) :: gf
         integer :: nspin 

         gf%nspin = nspin
      end subroutine set_grid_function_nspin
      

     subroutine new_fde_grid(g,size)
       type(fde_grid) :: g
       integer :: size

       g%npoints = size
       allocate(g%r(3,size))
       allocate(g%w(size))
     end subroutine new_fde_grid


     subroutine del_fde_grid(g)
        type(fde_grid) :: g

        g%npoints = 0
        deallocate(g%r)
        deallocate(g%w)
     end subroutine del_fde_grid

end
