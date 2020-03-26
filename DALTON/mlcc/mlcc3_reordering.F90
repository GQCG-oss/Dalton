!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
module mlcc3_reordering
!
!
!  Contains the reorder routine and various similar routines
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
   use mlcc_block_import
   use mlcc_typedef
   use mlcc_work
!
!
contains
!   
   subroutine mlcc3_reorder(original,n_original,reordered,n_reordered,list)
!
!  Copy and reorder an array
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: Read in two arrays, not necessarily the same length and a list array.
!           Copy elements from one array into the other using the order provided
!           by the list.
!
      implicit none
!      
      integer, intent(in)                             :: n_reordered, n_original
!      
      real(dp), dimension(n_original), intent(in)     :: original
      real(dp), dimension(n_reordered), intent(inout) :: reordered
!   
      integer, dimension(n_reordered), intent(in)     :: list
!      
      integer                                         :: i
!
      do i = 1,n_reordered
         reordered(i) = original(list(i))
      end do
!      
   end subroutine mlcc3_reorder
!
!
   subroutine abc_to_bac(abc,bac,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(abc) to W(bac)
!           
!
      implicit none
!      
      integer, intent(in)        :: dime
      real(dp), dimension(dime**3)  :: bac, abc
!      
      integer  :: a, b, c
      integer  :: offset_1, offset_2
!
      offset_1 = 0
!
      do c = 1,dime
         do a = 1,dime
            do b = 1, dime
!            
               offset_1 = offset_1 + 1
               offset_2 = dime**2*(c-1) + dime*(b-1) + a
!              
               bac(offset_1) = abc(offset_2)
!              
            end do
         end do
      end do
!
   end subroutine abc_to_bac
!
   subroutine abc_to_acb(abc,acb,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(abc) to W(acb)
!           
!
      implicit none
!      
      integer, intent(in)        :: dime
      real(dp), dimension(dime**3)  :: acb, abc
!      
      integer  :: a, b, c
      integer  :: offset_1, offset_2
!
      offset_1 = 0
!
      do b = 1,dime
         do c = 1,dime
            do a = 1, dime
!            
               offset_1 = offset_1 + 1
               offset_2 = dime**2*(c-1) + dime*(b-1) + a
!
               acb(offset_1) = abc(offset_2)
!
            end do
         end do
      end do
!      
   end subroutine abc_to_acb
!
   subroutine abc_to_bca(abc,bca,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(abc) to W(bca)
!           
!
      implicit none
!      
      integer, intent(in)        :: dime
      real(dp), dimension(dime**3)  :: bca, abc
!      
      integer  :: a, b, c
      integer  :: offset_1, offset_2
!
      offset_1 = 0
!
      do b = 1,dime
         do a = 1,dime
            do c = 1, dime
!            
               offset_1 = offset_1 + 1
               offset_2 = dime**2*(c-1) + dime*(b-1) + a
!              
               bca(offset_1) = abc(offset_2)
!              
            end do
         end do
      end do
!
   end subroutine abc_to_bca
!
   subroutine abc_to_cab(abc,cab,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(abc) to W(cab)
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3)  :: cab, abc
!      
      integer  :: a, b, c
      integer  :: offset_1, offset_2
!
      offset_1 = 0
!
      do a = 1,dime
         do c = 1,dime
            do b = 1, dime
!            
               offset_1 = offset_1 + 1
               offset_2 = dime**2*(c-1) + dime*(b-1) + a
!              
               cab(offset_1) = abc(offset_2)
!              
            end do
         end do
      end do
!
   end subroutine abc_to_cab
!
!
   subroutine abc_to_cba(abc,cba,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(abc) to W(cba)
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3)  :: cba, abc
!      
      integer  :: a, b, c
      integer  :: offset_1, offset_2
!
      offset_1 = 0
!
      do a = 1,dime
         do b = 1,dime
            do c = 1, dime
!            
               offset_1 = offset_1 + 1
               offset_2 = dime**2*(c-1) + dime*(b-1) + a
!              
               cba(offset_1) = abc(offset_2)
!              
            end do
         end do
      end do
!
   end subroutine abc_to_cba
!
!
!--------------------------------------------------------------
!--------------------------------------------------------------
!
   subroutine bac_to_abc(bac,abc,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(bac) to W(abc)
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3)  :: bac, abc
!      
!
      call abc_to_bac(bac,abc,dime)
!
!      
   end subroutine bac_to_abc
!
!
   subroutine acb_to_abc(acb,abc,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(acb) to W(abc)
!           
!
      implicit none
!      
      integer, intent(in)        :: dime
      real(dp), dimension(dime**3)  :: acb, abc
!      
!
      call abc_to_acb(acb,abc,dime)
!
!
   end subroutine acb_to_abc
!
!
   subroutine bca_to_abc(bca,abc,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(bca) to W(abc)
!           
!
      implicit none
!      
      integer, intent(in)        :: dime
      real(dp), dimension(dime**3)  :: bca, abc
!      
!
      call abc_to_cab(bca,abc,dime)
!
!
   end subroutine bca_to_abc
!
!
   subroutine cab_to_abc(cab,abc,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(cab) to W(abc)
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3)  :: cab, abc
!      
!
      call abc_to_bca(cab,abc,dime)
!
!
   end subroutine cab_to_abc
!
!
   subroutine cba_to_abc(cba,abc,dime)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(cba) to W(abc)
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3)  :: cba, abc
!      
!
      call abc_to_cba(cba,abc,dime)
!      
!      
   end subroutine cba_to_abc
!
!
   subroutine mlcc3_transposer(ab,ba,n_a,n_b)
!
!  W reorder
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: reorder W(cba) to W(abc)
!           
!
      implicit none
!      
      integer, intent(in)                          :: n_a, n_b
      integer                                      :: ab_off, ba_off
      integer                                      :: a, b
!      
      real(dp), dimension(n_a*n_b), intent(in)     :: ab
      real(dp), dimension(n_b*n_a), intent(inout)  :: ba
!      
!
      do a = 1,n_a
         do b = 1,n_b
!         
            ab_off   = n_a*(b-1) + a
            ba_off   = n_b*(a-1) + b
!            
            ba(ba_off) = ab(ab_off)
!            
         enddo
      enddo
!      
   end subroutine mlcc3_transposer
!
!
!----------------------------------------------------
!----------------------------------------------------
!
   subroutine u_bca_calc(t_abc,u_bca,dime)
!
!  u^bca calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate u^bca = 2t^bca - t^bac - t^cba
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_bca
!      
      integer  :: a, b, c
      integer  :: inc_off, bca_off, bac_off, cba_off
!
!
!$omp parallel do schedule(static) private(b,a,c,inc_off,bca_off,bac_off,cba_off)
      do b = 1,dime
         do a = 1,dime
            do c = 1, dime
!            
               inc_off = dime**2*(b-1) + dime*(a-1) + c
               bca_off = dime**2*(c-1) + dime*(b-1) + a
               bac_off = dime**2*(b-1) + dime*(c-1) + a
               cba_off = dime**2*(c-1) + dime*(a-1) + b
!              
               u_bca(inc_off) = 2*t_abc(bca_off) - t_abc(bac_off) - t_abc(cba_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine u_bca_calc
!
   subroutine u_bac_dax(t_abc,u_bac,dime)
!
!  u^bac calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate u^bac = 2t^bac - t^bca - t^cab
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_bac
!      
      integer  :: i,j, dime2
      integer  :: inc_off, bac_off, bca_off, cab_off
!
      dime2 = dime**2
!
      call dzero(u_bac,dime**3)
!
!$omp parallel do schedule(static) private(i,j,inc_off,bac_off,bca_off,cab_off)
      do i = 1,dime
         do j = 1,dime
!
            inc_off = dime2*(i-1) + dime*(j-1) + 1
!
            bac_off = dime2*(i-1) + j
!
            bca_off = dime*(i-1) + j
!
            cab_off = dime2*(j-1) + i
!
            call daxpy(dime,two,t_abc(bac_off),dime,u_bac(inc_off),1)
            call daxpy(dime,-one,t_abc(bca_off),dime2,u_bac(inc_off),1)
            call daxpy(dime,-one,t_abc(cab_off),dime,u_bac(inc_off),1)
!
         end do
      end do
!$omp end parallel do
!
   end subroutine u_bac_dax
!
!
   subroutine u_bac_calc(t_abc,u_bac,dime)
!
!  u^bac calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate u^bac = 2t^bac - t^bca - t^cab
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_bac
!      
      integer  :: a, b, c
      integer  :: inc_off, bac_off, bca_off, cab_off
!
!$omp parallel do schedule(static) private(c,a,b,inc_off,bac_off,bca_off,cab_off)
      do c = 1,dime
         do a = 1,dime
            do b = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(a-1) + b
               bac_off = dime**2*(c-1) + dime*(b-1) + a
               bca_off = dime**2*(b-1) + dime*(c-1) + a
               cab_off = dime**2*(a-1) + dime*(b-1) + c
!              
               u_bac(inc_off) = 2*t_abc(bac_off) - t_abc(bca_off) - t_abc(cab_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine u_bac_calc
!
!
   subroutine u_abc_calc(t_abc,u_abc,dime)
!
!  u^abc calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate u^abc = 2t^abc - t^cba - t^acb
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_abc
!      
      integer  :: a, b, c
      integer  :: inc_off, cba_off, acb_off
!
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,cba_off,acb_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               cba_off = dime**2*(a-1) + dime*(b-1) + c
               acb_off = dime**2*(b-1) + dime*(c-1) + a
!              
               u_abc(inc_off) = 2*t_abc(inc_off) - t_abc(cba_off) - t_abc(acb_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine u_abc_calc
!
!
   subroutine u_acb_calc(t_abc,u_acb,dime)
!
!  u^acb calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate u^acb = 2t^acb - t^cab - t^abc
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_acb
!      
      integer  :: a, b, c
      integer  :: inc_off, acb_off, cab_off
!
!
!$omp parallel do schedule(static) private(b,c,a,inc_off,acb_off,cab_off)
      do b = 1,dime
         do c = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(b-1) + dime*(c-1) + a
               acb_off = dime**2*(c-1) + dime*(b-1) + a
               cab_off = dime**2*(c-1) + dime*(a-1) + b
!
               u_acb(inc_off) = 2*t_abc(acb_off) - t_abc(cab_off) - t_abc(inc_off)
!
            end do
         end do
      end do
!$omp end parallel do
!      
   end subroutine u_acb_calc
!
!
   subroutine u_cba_calc(t_abc,u_cba,dime)
!
!  u^cba calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate u^cba = 2t^cba - t^abc - t^bca
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_cba
!      
      integer  :: a, b, c
      integer  :: inc_off, cba_off, bca_off
!
!
!$omp parallel do schedule(static) private(a,b,c,inc_off,cba_off,bca_off)
      do a = 1,dime
         do b = 1,dime
            do c = 1, dime
!            
               inc_off = dime**2*(a-1) + dime*(b-1) + c
               cba_off = dime**2*(c-1) + dime*(b-1) + a
               bca_off = dime**2*(c-1) + dime*(a-1) + b
!
               u_cba(inc_off) = 2*t_abc(cba_off) - t_abc(inc_off) - t_abc(bca_off)
!
            end do
         end do
      end do
!$omp end parallel do
!      
   end subroutine u_cba_calc
!
!
   subroutine u_cab_calc(t_abc,u_cab,dime)
!
!  u^cba calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate u^cab = 2t^cab - t^acb - t^bac
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_cab
!      
      integer  :: a, b, c
      integer  :: inc_off, cab_off, acb_off, bac_off
!
!
!$omp parallel do schedule(static) private(a,c,b,inc_off,cab_off,acb_off,bac_off)
      do a = 1,dime
         do c = 1,dime
            do b = 1, dime
!            
               inc_off = dime**2*(a-1) + dime*(c-1) + b
               cab_off = dime**2*(c-1) + dime*(b-1) + a
               acb_off = dime**2*(c-1) + dime*(a-1) + b
               bac_off = dime**2*(a-1) + dime*(b-1) + c
!              
               u_cab(inc_off) = 2*t_abc(cab_off) - t_abc(acb_off) - t_abc(bac_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine u_cab_calc
!
!
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
   subroutine t_abc_cba(t_abc,u_abc,dime)
!
!  abc - cba calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: calculate t^abc - t^cba
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_abc
!      
      integer  :: a, b, c
      integer  :: inc_off, cba_off
!
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,cba_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               cba_off = dime**2*(a-1) + dime*(b-1) + c
!              
               u_abc(inc_off) = t_abc(inc_off) - t_abc(cba_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
!
   end subroutine t_abc_cba
!
!
   subroutine t_acb_cab(t_abc,u_acb,dime)
!
!  acb - cab calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: calculate t^acb - t^cab
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_acb
!      
      integer  :: a, b, c
      integer  :: inc_off, acb_off, cab_off
!
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,acb_off,cab_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               acb_off = dime**2*(b-1) + dime*(c-1) + a
               cab_off = dime**2*(b-1) + dime*(a-1) + c
!              
               u_acb(inc_off) = t_abc(acb_off) - t_abc(cab_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine t_acb_cab
!
!
   subroutine t_bac_bca(t_abc,u_bac,dime)
!
!  bac - bca calclator
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: calculate t^bac - t^bca
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)  :: t_abc 
      real(dp), dimension(dime**3), intent(out) :: u_bac
!      
      integer  :: a, b, c
      integer  :: inc_off, bac_off, bca_off
!
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,bac_off,bca_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               bac_off = dime**2*(c-1) + dime*(a-1) + b
               bca_off = dime**2*(a-1) + dime*(c-1) + b
!              
               u_bac(inc_off) = t_abc(bac_off) - t_abc(bca_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine t_bac_bca
!
!
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
!
   subroutine add_bac(t_bac,t_abc,dime)
!
!  bac adder
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: add bac ordered vector to t^abc
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)     :: t_bac 
      real(dp), dimension(dime**3), intent(inout)  :: t_abc
!      
      integer  :: a, b, c
      integer  :: inc_off, bac_off
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,bac_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               bac_off = dime**2*(c-1) + dime*(a-1) + b
!              
               t_abc(inc_off) = t_abc(inc_off) + t_bac(bac_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine add_bac
!
!
   subroutine add_acb(t_acb,t_abc,dime)
!
!  acb adder
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: add acb ordered vector to t^abc
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)     :: t_acb
      real(dp), dimension(dime**3), intent(inout)  :: t_abc
!      
      integer  :: a, b, c
      integer  :: inc_off, acb_off
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,acb_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               acb_off = dime**2*(b-1) + dime*(c-1) + a
!              
               t_abc(inc_off) = t_abc(inc_off) + t_acb(acb_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine add_acb
!
!
   subroutine add_cab(t_cab,t_abc,dime)
!
!  cab adder
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: add acb ordered vector to t^abc
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)     :: t_cab
      real(dp), dimension(dime**3), intent(inout)  :: t_abc
!      
      integer  :: a, b, c
      integer  :: inc_off, cab_off
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,cab_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               cab_off = dime**2*(b-1) + dime*(a-1) + c
!              
               t_abc(inc_off) = t_abc(inc_off) + t_cab(cab_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine add_cab
!
!
   subroutine add_bca(t_bca,t_abc,dime)
!
!  bca adder
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: add bca ordered vector to t^abc
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)     :: t_bca
      real(dp), dimension(dime**3), intent(inout)  :: t_abc
!      
      integer  :: a, b, c
      integer  :: inc_off, bca_off
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,bca_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               bca_off = dime**2*(a-1) + dime*(c-1) + b
!              
               t_abc(inc_off) = t_abc(inc_off) + t_bca(bca_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine add_bca
!
!
   subroutine add_cba(t_cba,t_abc,dime)
!
!  cba adder
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: add cba ordered vector to t^abc
!           
!
      implicit none
!      
      integer, intent(in)           :: dime
      real(dp), dimension(dime**3), intent(in)     :: t_cba
      real(dp), dimension(dime**3), intent(inout)  :: t_abc
!      
      integer  :: a, b, c
      integer  :: inc_off, cba_off
!
!$omp parallel do schedule(static) private(c,b,a,inc_off,cba_off)
      do c = 1,dime
         do b = 1,dime
            do a = 1, dime
!            
               inc_off = dime**2*(c-1) + dime*(b-1) + a
               cba_off = dime**2*(a-1) + dime*(b-1) + c
!              
               t_abc(inc_off) = t_abc(inc_off) + t_cba(cba_off)
!              
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine add_cba
!
!
!-------------------------------------------------------------
!-------------------------------------------------------------
!
!  Integral reordering routines
!  ----------------------------
!
!
   subroutine bDck_order(unsorted,sorted,vir_act,vir_gen,k,n_k)
!
!  bDck order
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: order (bD|ck) integrals as D b c, k
!           
!
      implicit none
!      
      integer, intent(in)  :: vir_act, vir_gen, k, n_k
!      
      real(dp), dimension(vir_act**2*vir_gen*n_k), intent(in)  :: unsorted
      real(dp), dimension(vir_act**2*vir_gen), intent(out)     :: sorted
!      
      integer  :: c, b, D
      integer  :: int_off, sort_off
!
!$omp parallel do schedule(static) private(c,b,D,sort_off,int_off)
      do c = 1,vir_act
         do b = 1,vir_act
            do D = 1,vir_gen
!            
               sort_off = vir_act*vir_gen*(c-1) + vir_gen*(b-1) + D
!               
               int_off = vir_act*vir_gen*vir_act*(k-1) &
              &        + vir_act*vir_gen*(c-1) &
              &        + vir_act*(D-1) &
              &        + b
!           
               sorted(sort_off) = unsorted(int_off)
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine bDck_order
!
!
   subroutine Dbkc_order(unsorted,sorted,vir_act,vir_gen,k,n_k)
!
!  Dbkc order
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: order (Db|kc) integrals as b c D, k
!           
!
      implicit none
!      
      integer, intent(in)  :: vir_act, vir_gen, k, n_k
!      
      real(dp), dimension(vir_act**2*vir_gen*n_k), intent(in)  :: unsorted
      real(dp), dimension(vir_act**2*vir_gen), intent(out)     :: sorted
!      
      integer  :: c, b, D
      integer  :: int_off, sort_off
!
!$omp parallel do schedule(static) private(D,c,b,sort_off,int_off)
      do D = 1,vir_gen
         do c = 1,vir_act
            do b = 1,vir_act
!            
               sort_off = vir_act*vir_act*(D-1) + vir_act*(c-1) + b
!               
               int_off = vir_gen*vir_act*n_k*(c-1) &
              &        + vir_gen*vir_act*(k-1) &
              &        + vir_gen*(b-1) &
              &        + D
!           
               sorted(sort_off) = unsorted(int_off)
            end do
         end do
      end do
!$omp end parallel do
!
   end subroutine Dbkc_order
!
!
   subroutine Ljck_order(unsorted,sorted,vir_act,vir_gen,occ_act,occ_gen,j,k,n_k)
!
!  Ljck order
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: order (Lj|ck) integrals as L c, j k
!           
!
      implicit none
!      
      integer, intent(in)  :: occ_act, occ_gen, vir_act, vir_gen
      integer, intent(in)  :: j, k, n_k
!      
      real(dp), dimension(occ_gen*occ_act*vir_act*n_k), intent(in)   :: unsorted
      real(dp), dimension(occ_gen*vir_act), intent(out)              :: sorted
!      
      integer  :: c, L
      integer  :: int_off, sort_off
!
!$omp parallel do schedule(static) private(c,L,sort_off,int_off)
      do c = 1,vir_act
         do L = 1,occ_gen
!         
            sort_off = occ_gen*(c-1) + L
!            
            int_off = occ_gen*occ_act*vir_act*(k-1) &
           &        + occ_gen*occ_act*(c-1) &
           &        + occ_gen*(j-1) &
           &        + L
!        
            sorted(sort_off) = unsorted(int_off)
         end do
      end do
!$omp end parallel do
!
   end subroutine Ljck_order
!
!
   subroutine jLkc_order(unsorted,sorted,vir_act,vir_gen,occ_act,occ_gen,j,k,n_k)
!
!  jLkc order
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: order (jL|kc) integrals as c L, j k
!           
!
      implicit none
!      
      integer, intent(in)  :: occ_act, occ_gen, vir_act, vir_gen
      integer, intent(in)  :: j, k, n_k
!      
      real(dp), dimension(occ_gen*occ_act*vir_act*n_k), intent(in)   :: unsorted
      real(dp), dimension(occ_gen*vir_act), intent(out)              :: sorted
!      
      integer  :: c, L
      integer  :: int_off, sort_off
!
!$omp parallel do schedule(static) private(L,c,sort_off,int_off)
      do L = 1,occ_gen
         do c = 1,vir_act
!         
            sort_off = vir_act*(L-1) + c
!            
            int_off = occ_act*occ_gen*n_k*(c-1) &
           &            + occ_act*occ_gen*(k-1) &
           &            + occ_act*(L-1) &
           &            + j
!        
            sorted(sort_off) = unsorted(int_off)
         end do
      end do
!$omp end parallel do
!
   end subroutine jLkc_order
!
!
   subroutine jbkc_order(unsorted,sorted,vir_act,vir_gen,occ_act,occ_gen,j,k,n_k)
!
!  jbkc order
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: order (jb|kc) integrals as b c, j k
!           
!
      implicit none
!      
      integer, intent(in)  :: occ_act, occ_gen, vir_act, vir_gen
      integer, intent(in)  :: j, k, n_k
!      
      real(dp), dimension(vir_act**2*occ_act*n_k), intent(in)  :: unsorted
      real(dp), dimension(vir_act**2), intent(out)             :: sorted
!      
      integer  :: b, c
      integer  :: int_off, sort_off
!
!$omp parallel do schedule(static) private(c,b,sort_off,int_off)
      do c = 1,vir_act
         do b = 1,vir_act
!         
            sort_off = vir_act*(c-1) + b
!            
            int_off = occ_act*vir_act*n_k*(c-1) &
           &        + occ_act*vir_act*(k-1) &
           &        + occ_act*(b-1) &
           &        + j
!        
            sorted(sort_off) = unsorted(int_off)
         end do
      end do
!$omp end parallel do
!
   end subroutine jbkc_order
!
!
end module mlcc3_reordering
