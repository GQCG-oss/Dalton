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
module mlcc3_omega
!
!
!  Add MLCC3 contributions to omega vector
!  Authors Henrik Koch and Rolf H. Myhre
!  March 2015
!
   use mlcc_work
   use mlcc_typedef
   use mlcc_block_import
   use mlcc3_active_spaces
   use mlcc3_data
   use mlcc3_reordering
   use mlcc3_various
!
!   
!  W-intermediates
   real(dp), dimension(:), pointer, private  :: u_abc    => null()
!
!  Triples amplitudes, all permutations
   real(dp), dimension(:), pointer, private  :: t_abc    => null()
!
!  Temporary virtual Omega 2 vector ordered as a,B,i,j
   real(dp), dimension(:), pointer, private  :: omega_v  => null()
!
!  Temporary occupied Omega 2 vector ordered as a,b,I,j
   real(dp), dimension(:), pointer, private  :: omega_o  => null()
!
!  Temporary T2 and C2 array with virtual general, a,D,j,i ordering
   real(dp), dimension(:), pointer, private  :: t2_aDji  => null()
   real(dp), dimension(:), pointer, private  :: c2_aDji  => null()
!
!  Temporary T2 and C2 array with virtual general, a,b,L,i ordering
   real(dp), dimension(:), pointer, private  :: t2_abLi  => null()
   real(dp), dimension(:), pointer, private  :: c2_abLi  => null()
!
!  Option variable for response calculations
   integer  :: response_integer
!      
!  Print units
!  -----------
!
!  Ground state integrals
   integer :: bDck_g, Ljck_g, jbkc_g, Dbkc_g, jLkc_g
!
!  Response integrals
   integer :: bDck_r, Ljck_r, Dbkc_r, jLkc_r
!
!      
contains
!
   subroutine mlcc3_omega_contrib(resp_control)
!
!  contributions calculator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: loop over active occupied indices, calculate the W intermediate
!           and add to the omega vector
!
      implicit none
!      
      integer, intent(in)  :: resp_control
!      
      integer  :: i, j, k
!      
      integer  :: a, b, abij_n, abij_t
!      
!     Timing variables
      real     :: time_tot, time_w, time_t3, time_o1, time_o2
      real     :: time_add, time_t2sq
      real     :: time_start, time_1, time_allo, time_deallo, time_zero
      real     :: time_w_sum, time_t3_sum, time_o1_sum, time_o2_sum
      real     :: time_open, time_close
!      
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
      time_w_sum     = 0.0
      time_t3_sum    = 0.0
      time_o1_sum     = 0.0
      time_o2_sum     = 0.0
!
      call cpu_time(time_start)
!
!
!     Allocate arrays
!     If mlcc3_active, we need two temporary omega vectors due to
!     different lengths of n_vir and n_vir_act and n_occ and n_occ_act, 
!     so we have two smaller arrays instead. Same for T_2 amplitudes
!     C_2 amplitudes only required for [Ĥ,C_2] terms
!     -----------------------------------------------------------------
!
      call cpu_time(time_1)
!      
      if(mlcc3_active) then
!      
         call work_allocator(t2_aDji,n_vir_act*n_vir_gen*n_occ_act**2)
         call work_allocator(t2_abLi,n_vir_act**2*n_occ_gen*n_occ_act)
!      
         if(resp_control .eq. 1) then
!      
            call work_allocator(c2_aDji,n_vir_act*n_vir_gen*n_occ_act**2)
            call work_allocator(c2_abLi,n_vir_act**2*n_occ_gen*n_occ_act)
!      
         end if
!      
         call work_allocator(omega_v,n_vir_act*n_vir_gen*n_occ_act**2)
         call work_allocator(omega_o,n_vir_act**2*n_occ_gen*n_occ_act)
!      
      else
!      
         call work_allocator(t2_aDji,n_vir**2*n_occ**2)
         t2_abLi => t2_aDji
!      
         if(resp_control .eq. 1) then
!      
            call work_allocator(c2_aDji,n_vir**2*n_occ**2)
            c2_abLi => c2_aDji
!      
         end if
!      
         call work_allocator(omega_v,n_vir**2*n_occ**2)
         omega_o => omega_v
!      
      end if
!
      call work_allocator(u_abc,n_vir_act**3)
!      
      call work_allocator(t_abc,n_vir_act**3)
!      
      call cpu_time(time_allo)
      time_allo = time_allo - time_1
!
!
!     Zero out temporary omega vector(s)
!     ----------------------------------
!
      call cpu_time(time_1)
!      
      if(mlcc3_active) then
         call dzero(omega_v,n_vir_act*n_vir_gen*n_occ_act**2)
         call dzero(omega_o,n_vir_act**2*n_occ_gen*n_occ_act)
      else
         call dzero(omega_v,n_vir_act**2*n_occ**2)
      end if
!
      call cpu_time(time_zero)
      time_zero = time_zero - time_1
!
!     Square up T_2 vectors
!     ---------------------
!
      call cpu_time(time_1)
!      
      response_integer = 0
!
      call mlcc3_t2_square
!
      if(resp_control .eq. 1) then
!
!        Square up C2 amplitudes as well
!
         response_integer = 1
!
         call mlcc3_t2_square
!
      end if
!
      call cpu_time(time_t2sq)
      time_t2sq = time_t2sq - time_1
!
!
!     Open files
!     ----------
!
      call cpu_time(time_1)
!      
      call mlcc3_file_opener(resp_control)
!
      call cpu_time(time_open)
      time_open = time_open - time_1
!
!
!     Restricted loop over occupied indices i >= j >= k
!     -------------------------------------------------
!
      do i = 1,n_occ_act
         do j = 1,i
            do k = 1,j
!            
               if (i .eq. j .and. j .eq. k) then
                  cycle
               end if
!                  
!
               call cpu_time(time_1)
!               
               if(resp_control .eq. 0 .or. resp_control .eq. 2) then
!
!                 Calculate standard T_3 amplitudes
!
                  response_integer = 0
!
!                 W intermediates               
                  call mlcc3_w_calc(i,j,k)
!
               else if(resp_control .eq. 1) then
!
!                 Calculate C_3 resonse amplitudes
!
!                 [Ĥ,C_2] term in W intermdeiates
                  response_integer = 1
!
                  call mlcc3_w_calc(i,j,k)
!               
!                 [[Ĥ,C_1],T_2] term in W intermediates
                  response_integer = 2
!
                  call mlcc3_w_calc(i,j,k)
!
               else
!               
                  call quit('Error, resp_control not 0, 1 or 2')
!               
               end if
!               
               call cpu_time(time_w)
               time_w = time_w - time_1
               time_w_sum = time_w_sum + time_w
!               
!               
!              Divide by epsilon - frequency. Frequency is zero
!              for an energy calculation
               call cpu_time(time_1)
!               
               if(resp_control .eq. 0 .or. resp_control .eq. 2) then
!
                  response_integer = 0
!
                  call mlcc3_t3_calculator(i,j,k)
!               
               else if(resp_control .eq. 1) then
!
                  response_integer = 1
!
                  call mlcc3_t3_calculator(i,j,k)
!               
               end if
!               
               call cpu_time(time_t3)
               time_t3 = time_t3 - time_1
               time_t3_sum = time_t3_sum + time_t3
!               
!
               call cpu_time(time_1)
!               
!              Omega 1 and Fock contributions.
               if(resp_control .eq. 0 .or. resp_control .eq. 1) then
!
!                 Standard [Ĥ,T_3] or [Ĥ,C_3], same code
!
                  response_integer = 0
!
                  call mlcc3_omega_1(i,j,k)
!
               else if(resp_control .eq. 2) then
!
!                 Standard [[Ĥ,C_1],T_3] Fock contributions
!
                  response_integer = 1
!
                  call mlcc3_omega_1(i,j,k)
!
               end if
!               
               call cpu_time(time_o1)
               time_o1 = time_o1 - time_1
               time_o1_sum = time_o1_sum + time_o1
!               
!               
               call cpu_time(time_1)
!               
!              Omega 2 contributions
               if(resp_control .eq. 0 .or. resp_control .eq. 1) then
!
!                 [Ĥ,T_3] terms or [Ĥ,C_3] terms
!               
                  response_integer = 0
!               
                  call mlcc3_omega_2(i,j,k)
!               
               else if(resp_control .eq. 2) then
!               
!                 [[Ĥ,C_1],T_3] terms
!               
                  response_integer = 1
!               
                  call mlcc3_omega_2(i,j,k)
!               
               end if
!               
               call cpu_time(time_o2)
               time_o2 = time_o2 - time_1
               time_o2_sum = time_o2_sum + time_o2
!               
!               
!              Print timings
!              -------------
!
               if (print_mlcc3 .ge. 5) then
!
                  write(lupri,*)
                  write(lupri,*)
                  write(lupri,*) 'Timings from mlcc3_omega in loop'
                  write(lupri,*) 'i, j, k', i, j, k
                  write(lupri,9999) 'mlcc3_w_calc', time_w
                  write(lupri,9999) 'mlcc3_t3_calc', time_t3
                  write(lupri,9999) 'mlcc3_omega_1', time_o1
                  write(lupri,9999) 'mlcc3_omega_2', time_o2
                  write(lupri,*)
                  write(lupri,*)
!
               end if
!
            end do
         end do
      end do
!      
!
      if(print_mlcc3 .ge. 100) then 
!
         write(lupri,*)
         write(lupri,*) 'Omega 2 from mlcc3'
         write(lupri,*) 'resp_control',resp_control
         write(lupri,*) 'a,i,b,j ordering'
         write(lupri,*) 'xxx'
!
         do j = 1,n_occ
            do b = 1,n_vir
               do i = 1,n_occ
                  do a = 1,n_vir

                     abij_n = (n_occ*n_vir**2)*(j-1) + (n_vir**2)*(i-1) + n_vir*(b-1) + a
                     abij_t = (n_occ*n_vir**2)*(i-1) + (n_vir**2)*(j-1) + n_vir*(a-1) + b
!
                     write(lupri,'(F14.10)') (omega_o(abij_n) + omega_o(abij_t))
!
                  end do
               end do
            end do
         end do
!
         write(lupri,*) 'xxx'
         write(lupri,*)
!
      end if
!
!
!     Add temporary Omega 2 to real Omega 2
!     -------------------------------------
      call cpu_time(time_1)
!
      call mlcc3_omega_add
!
      call cpu_time(time_add)
      time_add = time_add - time_1
!
!
!     Close files
!     ----------
!
      call cpu_time(time_1)
!      
      call mlcc3_file_closer(resp_control)
!
      call cpu_time(time_close)
      time_close = time_close - time_1
!
!
!     Deallocate arrays
!     -----------------
!
      call cpu_time(time_1)
!      
      call work_deallocator(t_abc)
!
      call work_deallocator(u_abc)
!
      if(mlcc3_active) then
!
         call work_deallocator(omega_o)
         call work_deallocator(omega_v)
!
         if(resp_control .eq. 1) then
!      
            call work_deallocator(c2_abLi)
            call work_deallocator(c2_aDji)
!      
         end if
!      
         call work_deallocator(t2_abLi)
         call work_deallocator(t2_aDji)
!      
      else
!
         call work_deallocator(omega_v)
!
         if(resp_control .eq. 1) then
!      
            call work_deallocator(c2_aDji)
!      
         end if
!      
         call work_deallocator(t2_aDji)
!
      end if
!      
      call cpu_time(time_deallo)
      time_deallo = time_deallo - time_1
!
!     
!     Total timing
      call cpu_time(time_tot)
!      
      time_tot = time_tot - time_start
!      
!
!     Print timings
!     -------------
!
      if (print_mlcc3 .ge. 4) then
!
         write(lupri,*)
         write(lupri,*)
         write(lupri,*) 'Timings from mlcc3_omega after loop'
         write(lupri,9999) 'Total', time_tot
         write(lupri,9999) 'allocation', time_allo
         write(lupri,9999) 'open files', time_open
         write(lupri,9999) 'zero omega', time_zero
         write(lupri,9999) 'square T_2', time_t2sq
         write(lupri,9999) 'mlcc3_w_calc', time_w_sum
         write(lupri,9999) 'mlcc3_t3_calc', time_t3_sum
         write(lupri,9999) 'mlcc3_omega_1', time_o1_sum
         write(lupri,9999) 'mlcc3_omega_2', time_o2_sum
         write(lupri,9999) 'add to omega2', time_add
         write(lupri,9999) 'close files', time_close
         write(lupri,9999) 'deallocation', time_deallo
         write(lupri,*)
         write(lupri,*)
!
      end if
!
   end subroutine mlcc3_omega_contrib
!      
!      
   subroutine mlcc3_w_calc(i,j,k)
!
!  W intermediate calculator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate the W intermediates
!
      implicit none
!      
      integer, intent(in)  :: i, j, k
      integer              :: n_t_aD_ij, n_t_ab_iL, n_occ_read
      integer              :: rec_number, r2_off
!
      integer              :: bDck_unit, Ljck_unit, ioerror
!
      real(dp)             :: beta
!
!     File names
      character(len=12)    :: bDck_name
      character(len=12)    :: Ljck_name
!
!     (bD|ck) and (Lj|ck) integrals
!
      real(dp), dimension(:), pointer  :: bDck     => null()
      real(dp), dimension(:), pointer  :: Ljck     => null()
!      
!      
!     Timing variables
      real     :: time_tot, time_dgvir, time_dgocc
      real     :: time_start, time_ro, time_rv, time_order
      real     :: time_bac, time_acb, time_cab, time_bca, time_cba
      real     :: time_1, time_2
!      
      real(dp), dimension(:), pointer  :: r2_aDji  => null()
      real(dp), dimension(:), pointer  :: r2_abLi  => null()
!
!
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
      time_dgvir  = 0.0
      time_dgocc  = 0.0
      time_ro     = 0.0
      time_rv     = 0.0
!      
!
      if(response_integer .eq. 0) then
!
!        Calculate T_3 amplitudes
!
         bDck_unit = bDck_g
         Ljck_unit = Ljck_g
!
         beta = zero
!
         r2_aDji  => t2_aDji
         r2_abLi  => t2_abLi
!
      else if(response_integer .eq. 1) then
!
!        Calculate [Ĥ,C_2] term in C_3
!
         bDck_unit = bDck_g
         Ljck_unit = Ljck_g
!
         beta = zero
!
         r2_aDji  => c2_aDji
         r2_abLi  => c2_abLi
!
      else if(response_integer .eq. 2) then
!
!        Calculate [[Ĥ,C_1],T_2] term in C_3
!
         bDck_unit = bDck_r
         Ljck_unit = Ljck_r
!
         beta = one
!
         r2_aDji  => t2_aDji
         r2_abLi  => t2_abLi
!
      else
!
         call quit('response_integer not 0, 1 or 2')
!
      end if
!
!
      call cpu_time(time_start)
!
!     Calculate som lengths
!     ---------------------
!
      n_occ_read  = n_occ_gen*n_vir_act
      n_t_aD_ij   = n_vir_act*n_vir_gen
      n_t_ab_iL   = n_occ_gen*n_vir_act**2
!      
!     Allocate space
!     --------------
!
      call work_allocator(bDck,n_vir_aag)
      call work_allocator(Ljck,n_occ_read)
!
!
!     Open files and read k integrals. Integrals on disk 
!     are organised as D b c, k and L c, j k
!     --------------------------------------------------
!
      call cpu_time(time_1)
!
      read(bDck_unit,rec=k,iostat=ioerror) bDck
!            
      call cpu_time(time_2)
      time_rv = time_rv + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading bDck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
      rec_number = n_occ_act*(k-1) + j
!      
      call cpu_time(time_1)
!
      read(Ljck_unit,rec=rec_number,iostat=ioerror) Ljck
!            
      call cpu_time(time_2)
      time_ro = time_ro + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
!     ------------------------------------------------------------------
!     |Calculate contributions to W. Do the matrix-marix multiplication| 
!     ------------------------------------------------------------------
!
!     t^aD_ij*(bD|ck)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDck,n_vir_gen,beta,t_abc,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!               
!     t^ab_iL*(Lj|ck)
!     ---------------
!
      r2_off = (n_occ_gen*n_vir_act**2)*(i-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Ljck,n_occ_gen,one,t_abc,n_vir_act**2)
!  
      call cpu_time(time_2)
      time_dgocc = time_dgocc + time_2 - time_1
!
!      
!      
!     t^bD_ji*(aD|ck)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(i-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDck,n_vir_gen,zero,u_abc,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!     t^ba_jL*(Li|ck)
!     ---------------
      rec_number = n_occ_act*(k-1) + i
!      
      call cpu_time(time_1)
!
      read(Ljck_unit,rec=rec_number,iostat=ioerror) Ljck
!            
      call cpu_time(time_2)
      time_ro = time_ro + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
!
      r2_off = (n_occ_gen*n_vir_act**2)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Ljck,n_occ_gen,one,u_abc,n_vir_act**2)
!  
      call cpu_time(time_2)
      time_dgocc = time_dgocc + time_2 - time_1
!
!
!     Reorder and add to W
!     --------------------
!
      call cpu_time(time_1)
!
      call add_bac(u_abc,t_abc,n_vir_act)
!
      call cpu_time(time_bac)
      time_bac = time_bac - time_1
!
!
!     Read j integrals
!     ----------------
!
      call cpu_time(time_1)
!
      read(bDck_unit,rec=j,iostat=ioerror) bDck
!            
      call cpu_time(time_2)
      time_rv = time_rv + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading bDck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
!     t^aD_ik*(cD|bj)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDck,n_vir_gen,zero,u_abc,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!     t^ac_iL*(Lk|bj)
!     ---------------
      rec_number = n_occ_act*(j-1) + k
!      
      call cpu_time(time_1)
!
      read(Ljck_unit,rec=rec_number,iostat=ioerror) Ljck
!            
      call cpu_time(time_2)
      time_ro = time_ro + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
!
      r2_off = (n_occ_gen*n_vir_act**2)*(i-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Ljck,n_occ_gen,one,u_abc,n_vir_act**2)
!  
      call cpu_time(time_2)
      time_dgocc = time_dgocc + time_2 - time_1
!
!
!     Reorder and add to W
!     --------------------
!
      call cpu_time(time_1)
!
      call add_acb(u_abc,t_abc,n_vir_act)
!
      call cpu_time(time_acb)
      time_acb = time_acb - time_1
!
!      
!
!     t^cD_ki*(aD|bj)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(k-1) + (n_vir_gen*n_vir_act)*(i-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDck,n_vir_gen,zero,u_abc,n_vir_act)
!
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!     t^ca_kL*(Li|bj)
!     ---------------
      rec_number = n_occ_act*(j-1) + i
!      
      call cpu_time(time_1)
!
      read(Ljck_unit,rec=rec_number,iostat=ioerror) Ljck
!            
      call cpu_time(time_2)
      time_ro = time_ro + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
!
      r2_off = (n_occ_gen*n_vir_act**2)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Ljck,n_occ_gen,one,u_abc,n_vir_act**2)
!  
      call cpu_time(time_2)
      time_dgocc = time_dgocc + time_2 - time_1
!
!     Reorder and add to W
!     --------------------
!
      call cpu_time(time_1)
!
      call add_cab(u_abc,t_abc,n_vir_act)
!
      call cpu_time(time_cab)
      time_cab = time_cab - time_1
!
!
!      
!     Read i integrals
!     ----------------
!
      call cpu_time(time_1)
!
      read(bDck_unit,rec=i,iostat=ioerror) bDck
!            
      call cpu_time(time_2)
      time_rv = time_rv + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading bDck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
!     t^bD_jk*(cD|ai)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDck,n_vir_gen,zero,u_abc,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!               
!     t^bc_jL*(Lk|ai)
!     ---------------
      rec_number = n_occ_act*(i-1) + k
!      
      call cpu_time(time_1)
!
      read(Ljck_unit,rec=rec_number,iostat=ioerror) Ljck
!            
      call cpu_time(time_2)
      time_ro = time_ro + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
!
      r2_off = (n_occ_gen*n_vir_act**2)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Ljck,n_occ_gen,one,u_abc,n_vir_act**2)
!  
      call cpu_time(time_2)
      time_dgocc = time_dgocc + time_2 - time_1
!
!
!     Reorder and add to W
!     --------------------
!
      call cpu_time(time_1)
!
      call add_bca(u_abc,t_abc,n_vir_act)
!
      call cpu_time(time_bca)
      time_bca = time_bca - time_1
!
!
!     t^cD_kj*(bD|ai)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(k-1) + (n_vir_gen*n_vir_act)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDck,n_vir_gen,zero,u_abc,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!               
!     t^cb_kL*(Lj|ai)
!     ---------------
      rec_number = n_occ_act*(i-1) + j
!      
      call cpu_time(time_1)
!
      read(Ljck_unit,rec=rec_number,iostat=ioerror) Ljck
!            
      call cpu_time(time_2)
      time_ro = time_ro + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in w_calc')
      end if
!
!
      r2_off = (n_occ_gen*n_vir_act**2)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Ljck,n_occ_gen,one,u_abc,n_vir_act**2)
!
      call cpu_time(time_2)
      time_dgocc = time_dgocc + time_2 - time_1
!
!
!     Reorder and add to W
!     --------------------
!
      call cpu_time(time_1)
!
      call add_cba(u_abc,t_abc,n_vir_act)
!
      call cpu_time(time_cba)
      time_cba = time_cba - time_1
!
!      
!     Deallocate space
!     ----------------
!      
      call work_deallocator(Ljck)
      call work_deallocator(bDck)
!      
!     Total timing
      call cpu_time(time_tot)
!      
      time_tot = time_tot - time_start
!      
      time_order = time_bac + time_acb + time_cab + time_bca + time_cba
!
!     Print timings
!     -------------
!
      if (print_mlcc3 .ge. 5) then
!
         write(lupri,*)
         write(lupri,*)
         write(lupri,*) 'Timings from mlcc3_w_calc'
         write(lupri,9999) 'Total', time_tot
         write(lupri,9999) 'read virtual', time_rv
         write(lupri,9999) 'read occupied', time_ro
         write(lupri,9999) 'virtual MxM', time_dgvir
         write(lupri,9999) 'occupied MxM', time_dgocc
         write(lupri,9999) 'reorder bac', time_bac
         write(lupri,9999) 'reorder acb', time_acb
         write(lupri,9999) 'reorder cab', time_cab
         write(lupri,9999) 'reorder bca', time_bca
         write(lupri,9999) 'reorder cba', time_cba
         write(lupri,9999) 'total reorder', time_order
         write(lupri,*)
         write(lupri,*)
!
      end if
!
!            
   end subroutine mlcc3_w_calc
!      
!
   subroutine mlcc3_t3_calculator(i,j,k)
!
!  Divide t_abc(abc) by epsilon^abc_ijk
!  Authors: Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: Calulate the ijk triples amplitudes by dividing W by epsilon
!  Frequency is set to zero in energy calculations
!
      implicit none
!      
      integer,intent(in)   :: i, j, k
      integer              :: i_off, j_off, k_off
      integer              :: a, b, c
      real(dp)             :: epsilon_ijk, epsilon_abc
!      
!
!     Calculate occupied active indices
!     ---------------------------------
      i_off = i + n_occ_inact
      j_off = j + n_occ_inact
      k_off = k + n_occ_inact
!
      epsilon_ijk = Fock_diagonal(i_off) &
     &            + Fock_diagonal(j_off) &
     &            + Fock_diagonal(k_off) 
!
!     Subtract frequency from denominator if calculating C_3
      if(response_integer .eq. 1) then
!
         epsilon_ijk = epsilon_ijk + freq
!
      end if
!
!$omp parallel do schedule(dynamic) private(c,b,a,epsilon_abc) shared(n_vir_act,epsilon_ijk)
      do c = 1,n_vir_act
         do b = 1,n_vir_act
            do a = 1,n_vir_act
!               
               if (a .eq. b .and. b .eq. c) then
!               
                  t_abc(n_vir_act**2*(c-1) + n_vir_act*(b-1) + a) = zero
!  
               else
!  
                  epsilon_abc = -one/(Fock_diagonal(n_occ + a)  &
              &               +       Fock_diagonal(n_occ + b)  &
              &               +       Fock_diagonal(n_occ + c)  &
              &               -       epsilon_ijk)
!  
                  t_abc(n_vir_act**2*(c-1) + n_vir_act*(b-1) + a) = &
              &      t_abc(n_vir_act**2*(c-1) + n_vir_act*(b-1) + a)*epsilon_abc
!  
               endif
!  
            end do
         end do
      end do
!$omp end parallel do
!      
   end subroutine mlcc3_t3_calculator
!
!
   subroutine mlcc3_omega_1(i,j,k)
!
!  MLCC3 contributions to omega 1 and Fock contribution to Omega 2
!  Authors: Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: Calclate and add the MLCC3 contributions to the omega 1 vector and ditto for
!           the Fock contribution to the omega 2 vector
!
      implicit none
!      
      integer,intent(in)   :: i, j, k
      integer              :: i_off, j_off, k_off
      integer              :: i_off2, j_off2, k_off2
      integer              :: a, b, ab, abij
!      
      integer              :: w_off
      integer              :: omega_off, F_off
!      
      integer              :: rec_number, ioerror
!      
      real(dp)             :: alpha
!      
      real(dp), dimension(:), pointer  :: L_jbkc   => null()
      real(dp), dimension(:), pointer  :: L_kbjc   => null()
!      
      real(dp), dimension(:), pointer  :: L_ibkc   => null()
      real(dp), dimension(:), pointer  :: L_kbic   => null()
!      
      real(dp), dimension(:), pointer  :: L_ibjc   => null()
      real(dp), dimension(:), pointer  :: L_jbic   => null()
!      
      real(dp), dimension(:), pointer  :: help     => null()
!      
      real(dp), dimension(:), pointer  :: omega_ab => null()
!
      real(dp), dimension(:), pointer  :: fock_mat => null()
!      
!     Timing variables
      real     :: time_tot, time_dgmv1, time_dgmvf, time_ome, time_order
      real     :: time_start, time_read, time_trans, time_L
      real     :: time_abc, time_acb, time_bac, time_1, time_2
!      
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
!
      if(response_integer .eq. 0) then
         fock_mat => mo_fock_mat_t1
      else if(response_integer .eq. 1) then
         fock_mat => mo_fock_mat_c1
      end if
!
      time_trans = 0.0
      time_read = 0.0
      time_L = 0.0
!      
      time_abc = 0.0
      time_acb = 0.0
      time_bac = 0.0
!      
      time_dgmv1 = 0.0
      time_dgmvf = 0.0
!      
      call cpu_time(time_start)
!
!
!     Calculate occupied active indices
!     ---------------------------------
      i_off = i + n_occ_inact
      j_off = j + n_occ_inact
      k_off = k + n_occ_inact
!      
      i_off2 = i + n_occ_inact - n_gen_inact
      j_off2 = j + n_occ_inact - n_gen_inact
      k_off2 = k + n_occ_inact - n_gen_inact
!      
!     Allocate arrays
!     ---------------
!
      call work_allocator(L_jbkc,n_vir_act**2)
      call work_allocator(L_kbjc,n_vir_act**2)
!      
      call work_allocator(L_ibkc,n_vir_act**2)
      call work_allocator(L_kbic,n_vir_act**2)
!      
      call work_allocator(L_ibjc,n_vir_act**2)
      call work_allocator(L_jbic,n_vir_act**2)
!      
      call work_allocator(help,n_vir_act**2)
!      
      call work_allocator(omega_ab,n_vir_act**2)
!      
!      
!      
      if(response_integer .eq. 0) then 
!      
!     ---------------------------------------
!     |Read in (jb|kc) and (kb|jc) integrals|
!     ---------------------------------------
!
      rec_number = n_occ_act*(k-1) + j
!      
      call cpu_time(time_1)
!      
      read(jbkc_g,rec=rec_number,iostat=ioerror) L_jbkc
!            
      call cpu_time(time_2)
      time_read = time_read + time_2 - time_1
!      
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading jbkc file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in Omega1 1')
      end if
!
      call cpu_time(time_1)
!      
      if (j .ne. k) then
!         
         call mlcc3_transposer(L_jbkc,L_kbjc,n_vir_act,n_vir_act)
!
      else
!
         call dcopy(n_vir_act**2,L_jbkc,1,L_kbjc,1)
!
      end if
!            
      call cpu_time(time_2)
      time_trans = time_trans + time_2 - time_1
!      
!      
!     ---------------------------------------
!     |Read in (ib|kc) and (kb|ic) integrals|
!     ---------------------------------------
!
      if (i .ne. j) then
!
         rec_number = n_occ_act*(k-1) + i
!      
         call cpu_time(time_1)
!      
         read(jbkc_g,rec=rec_number,iostat=ioerror) L_ibkc
!            
         call cpu_time(time_2)
         time_read = time_read + time_2 - time_1
!      
         if(ioerror .ne. 0) then
            write(lupri,*) "Error reading jbkc file"
            write(lupri,*) "Error code: ", ioerror
            call quit('Error reading integral file in Omega1 3')
         end if
!
         call cpu_time(time_1)
!      
         if (i .ne. k) then
!            
            call mlcc3_transposer(L_ibkc,L_kbic,n_vir_act,n_vir_act)
!
         else
!
            call dcopy(n_vir_act**2,L_ibkc,1,L_kbic,1)
!
         end if
!
         call cpu_time(time_2)
         time_trans = time_trans + time_2 - time_1
!      
      end if
!
!     ---------------------------------------
!     |Read in (jb|ic) and (ib|jc) integrals|
!     ---------------------------------------
!
      if (i .ne. k) then
!
         rec_number = n_occ_act*(i-1) + j
!      
         call cpu_time(time_1)
!      
         read(jbkc_g,rec=rec_number,iostat=ioerror) L_jbic
!            
         call cpu_time(time_2)
         time_read = time_read + time_2 - time_1
!      
         if(ioerror .ne. 0) then
            write(lupri,*) "Error reading jbkc file"
            write(lupri,*) "Error code: ", ioerror
            call quit('Error reading integral file in Omega1 5')
         end if
!
         call cpu_time(time_1)
!      
         if (i .ne. j) then
!            
            call mlcc3_transposer(L_jbic,L_ibjc,n_vir_act,n_vir_act)
!
         else
!
            call dcopy(n_vir_act**2,L_jbic,1,L_ibjc,1)
!
         end if
!
         call cpu_time(time_2)
         time_trans = time_trans + time_2 - time_1
!      
      end if
!
!
!     -------------------------------------
!     |Calculate L_jbkc, L_jbic and L_kbic|
!     -------------------------------------
!
      call cpu_time(time_1)
!         
      if (j .ne. k) then
!
         call dcopy(n_vir_act**2,L_jbkc,1,help,1)
!
         call dscal(n_vir_act**2,two,L_jbkc,1)
!
         call daxpy(n_vir_act**2,-one,L_kbjc,1,L_jbkc,1)
!
         call dscal(n_vir_act**2,two,L_kbjc,1)
!
         call daxpy(n_vir_act**2,-one,help,1,L_kbjc,1)
!
      end if
!
      if (i .ne. j .and. i .ne. k) then
!
!
         call dcopy(n_vir_act**2,L_kbic,1,help,1)
!
         call dscal(n_vir_act**2,two,L_kbic,1)
!
         call daxpy(n_vir_act**2,-one,L_ibkc,1,L_kbic,1)
!
         call dscal(n_vir_act**2,two,L_ibkc,1)
!
         call daxpy(n_vir_act**2,-one,help,1,L_ibkc,1)
!
      end if
!
      if (i .ne. k .and. i .ne. j) then
!
!
         call dcopy(n_vir_act**2,L_jbic,1,help,1)
!
         call dscal(n_vir_act**2,two,L_jbic,1)
!
         call daxpy(n_vir_act**2,-one,L_ibjc,1,L_jbic,1)
!
         call dscal(n_vir_act**2,two,L_ibjc,1)
!
         call daxpy(n_vir_act**2,-one,help,1,L_ibjc,1)
!
      end if
!
      end if !response_integer
!
      call cpu_time(time_L)
      time_L = time_L - time_1
!      
!
!     -----------------------------------------------------------
!     |Calculate (t^abc - t^cba)*L_jbkc and 2(t^abc - t^cba)F_kc|
!     -----------------------------------------------------------
!      
!     Construct t^abc - t^cba
!     -----------------------
!
      call cpu_time(time_1)
!
      call t_abc_cba(t_abc,u_abc,n_vir_act)
!      
      call cpu_time(time_abc)
      time_abc = time_abc - time_1
!
!
      if(response_integer .eq. 0) then 
!      
!     Add to Omega 1
!     --------------
!
      omega_off   = n_vir*(i_off - 1) + 1
!
      call cpu_time(time_1)
!
      call dgemv('N',n_vir_act,n_vir_act**2,one,u_abc,n_vir_act,L_jbkc,1,one,omega1(omega_off),1)
!
      call cpu_time(time_2)
      time_dgmv1 = time_dgmv1 + time_2 - time_1
!
      end if
!
!
!     Calculate Omega 2 contribution
!     ------------------------------
!
      F_off = n_orbitals*(k_off-1) + n_occ + 1
!
      omega_off   = (n_occ_gen*n_vir_act**2)*(j-1) + (n_vir_act**2)*(i_off2-1) + 1
!
      call cpu_time(time_1)
!
      call dgemv('N',n_vir_act**2,n_vir_act,one,u_abc,n_vir_act**2,fock_mat(F_off),1, &
     &            one,omega_o(omega_off),1)
!
      call cpu_time(time_2)
      time_dgmvf = time_dgmvf + time_2 - time_1
!
!
!     -----------------------------------------------------------
!     |Calculate (t^cba - t^abc)*L_jbic and 2(t^cba - t^abc)F_ic|
!     -----------------------------------------------------------
!      
      if (i .ne. k) then
!
!
         if(response_integer .eq. 0) then 
!      
!        Add to Omega 1
!        --------------
!      
         omega_off   = n_vir*(k_off - 1) + 1
!      
         call cpu_time(time_1)
!
         call dgemv('N',n_vir_act,n_vir_act**2,-one,u_abc,n_vir_act,L_jbic,1,one,omega1(omega_off),1)
!
         call cpu_time(time_2)
         time_dgmv1 = time_dgmv1 + time_2 - time_1
!
         end if
!
!
!        Calculate Omega 2 contribution
!        ------------------------------
!
         F_off = n_orbitals*(i_off-1) + n_occ + 1
!
         call cpu_time(time_1)
!
         omega_off   = (n_occ_gen*n_vir_act**2)*(j-1) + (n_vir_act**2)*(k_off2-1) + 1
!
         call dgemv('N',n_vir_act**2,n_vir_act,-one,u_abc,n_vir_act**2,fock_mat(F_off),1, &
     &            one,omega_o(omega_off),1)
!
         call cpu_time(time_2)
         time_dgmvf = time_dgmvf + time_2 - time_1
!
!
      end if
!
!     -----------------------------------------------------------
!     |Calculate (t^acb - t^cab)*L_kbjc and 2(t^acb - t^cab)F_jc|
!     -----------------------------------------------------------
!      
      if (j .ne. k) then
!
!
!        Construct t^acb - t^cab
!        -----------------------
!
         call cpu_time(time_1)
!
         call t_acb_cab(t_abc,u_abc,n_vir_act)
!      
         call cpu_time(time_acb)
         time_acb = time_acb - time_1
!
!        Add to Omega 1
!        --------------
!      
         if(response_integer .eq. 0) then 
!      
         omega_off   = n_vir*(i_off - 1) + 1
!      
         call cpu_time(time_1)
!
         call dgemv('N',n_vir_act,n_vir_act**2,one,u_abc,n_vir_act,L_kbjc,1,one,omega1(omega_off),1)
!
         call cpu_time(time_2)
         time_dgmv1 = time_dgmv1 + time_2 - time_1
!
         end if
!
!
!        Calculate Omega 2 contribution
!        ------------------------------
!
         F_off = n_orbitals*(j_off-1) + n_occ + 1
!
         omega_off   = (n_occ_gen*n_vir_act**2)*(k-1) + (n_vir_act**2)*(i_off2-1) + 1
!
         call cpu_time(time_1)
!
         call dgemv('N',n_vir_act**2,n_vir_act,one,u_abc,n_vir_act**2,fock_mat(F_off),1, &
     &            one,omega_o(omega_off),1)
!
         call cpu_time(time_2)
         time_dgmvf = time_dgmvf + time_2 - time_1
!
!
         if (i .ne. j) then
!
!           -----------------------------------------------------------
!           |Calculate (t^cab - t^acb)*L_kbic and 2(t^cab - t^acb)F_ic|
!           -----------------------------------------------------------
!      
!           Add to Omega 1
!           --------------
!
            if(response_integer .eq. 0) then 
!      
            omega_off   = n_vir*(j_off - 1) + 1
!      
            call cpu_time(time_1)
!
            call dgemv('N',n_vir_act,n_vir_act**2,-one,u_abc,n_vir_act,L_kbic,1,one,omega1(omega_off),1)
!
            call cpu_time(time_2)
            time_dgmv1 = time_dgmv1 + time_2 - time_1
!
            end if
!
!
!           Calculate Omega 2 contribution
!           ------------------------------
!
            F_off = n_orbitals*(i_off-1) + n_occ + 1
!
            omega_off   = (n_occ_gen*n_vir_act**2)*(k-1) + (n_vir_act**2)*(j_off2-1) + 1
!
            call cpu_time(time_1)
!
            call dgemv('N',n_vir_act**2,n_vir_act,-one,u_abc,n_vir_act**2,fock_mat(F_off),1, &
     &            one,omega_o(omega_off),1)
!
            call cpu_time(time_2)
            time_dgmvf = time_dgmvf + time_2 - time_1
!
!
         end if
!         
      end if
!      
!      
!     -----------------------------------------------------------
!     |Calculate (t^bac - t^bca)*L_ibkc and 2(t^bac - t^bca)F_kc|
!     -----------------------------------------------------------
!
      if (i .ne. j) then
!
!
!        Construct t^bac - t^bca
!        -----------------------
!
         call cpu_time(time_1)
!
         call t_bac_bca(t_abc,u_abc,n_vir_act)
!      
         call cpu_time(time_bac)
         time_bac = time_bac - time_1
!
!
!        Add to Omega 1
!        --------------
!      
         if(response_integer .eq. 0) then 
!      
         omega_off   = n_vir*(j_off - 1) + 1
!      
         call cpu_time(time_1)
!
         call dgemv('N',n_vir_act,n_vir_act**2,one,u_abc,n_vir_act,L_ibkc,1,one,omega1(omega_off),1)
!
         call cpu_time(time_2)
         time_dgmv1 = time_dgmv1 + time_2 - time_1
!
         end if
!
!
!        Calculate Omega 2 contribution
!        ------------------------------
!
         F_off = n_orbitals*(k_off-1) + n_occ + 1
!
         omega_off   = (n_occ_gen*n_vir_act**2)*(i-1) + (n_vir_act**2)*(j_off2-1) + 1
!
         call cpu_time(time_1)
!
         call dgemv('N',n_vir_act**2,n_vir_act,one,u_abc,n_vir_act**2,fock_mat(F_off),1, &
        &           one,omega_o(omega_off),1)
!
         call cpu_time(time_2)
         time_dgmvf = time_dgmvf + time_2 - time_1
!
!
         if (i .ne. k) then
!
!           -----------------------------------------------------------
!           |Calculate (t^bca - t^bac)*L_jcib and 2(t^bca - t^bac)F_jc|
!           -----------------------------------------------------------
!
!      
!           Add to Omega 1
!           --------------
!      
            if(response_integer .eq. 0) then 
!      
            omega_off   = n_vir*(k_off - 1) + 1
!      
            call cpu_time(time_1)
!
            call dgemv('N',n_vir_act,n_vir_act**2,-one,u_abc,n_vir_act,L_ibjc,1,one,omega1(omega_off),1)
!      
            call cpu_time(time_2)
            time_dgmv1 = time_dgmv1 + time_2 - time_1
!
            end if
!
!
!           Calculate Omega 2 contribution
!           ------------------------------
!
            F_off = n_orbitals*(j_off-1) + n_occ + 1
!
            call cpu_time(time_1)
!
            omega_off = (n_occ_gen*n_vir_act**2)*(i-1) + (n_vir_act**2)*(k_off2-1) + 1
!
            call dgemv('N',n_vir_act**2,n_vir_act,-one,u_abc,n_vir_act**2,fock_mat(F_off),1, &
           &           one,omega_o(omega_off),1)
!
            call cpu_time(time_2)
            time_dgmvf = time_dgmvf + time_2 - time_1
!
!
         end if
!
      end if
!      
!
!     Deallocate arrays
!     -----------------
!
      call work_deallocator(omega_ab)
!      
      call work_deallocator(help)
!      
      call work_deallocator(L_jbic)
      call work_deallocator(L_ibjc)
!
      call work_deallocator(L_kbic)
      call work_deallocator(L_ibkc)
!
      call work_deallocator(L_kbjc)
      call work_deallocator(L_jbkc)
!      
!     
!     Total timing
      call cpu_time(time_tot)
!      
      time_tot = time_tot - time_start
!      
      time_order = time_abc + time_acb + time_bac
!
!
!     Print timings
!     -------------
!
      if (print_mlcc3 .ge. 5) then
!
         write(lupri,*)
         write(lupri,*)
         write(lupri,*) 'Timings from mlcc3_omega_1'
         write(lupri,9999) 'Total', time_tot
         write(lupri,9999) 'read integral', time_read
         write(lupri,9999) 'transpose', time_trans
         write(lupri,9999) 'calculate L', time_L
         write(lupri,9999) 'abc - cba', time_abc
         write(lupri,9999) 'acb - cab', time_acb
         write(lupri,9999) 'bac - bca', time_bac
         write(lupri,9999) 'total reorder', time_order
         write(lupri,9999) 'calculate omega1', time_dgmv1
         write(lupri,9999) 'calculate Fock', time_dgmvf
         write(lupri,*)
         write(lupri,*)
!
      end if
!
   end subroutine mlcc3_omega_1
!
!
!
   subroutine mlcc3_omega_2(i,j,k)
!
!  Omega 2 contributions
!  Authors: Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: Calculate and add the occupied and virtual 
!           triples contributions to the omega 2 vector
!
      implicit none
!      
      integer,intent(in)   :: i, j, k
      integer              :: i_off, j_off, k_off
      integer              :: L, t_off
!      
      integer              :: omega_off
!      
      integer              :: ioerror, rec_number
      integer              :: jLkc_print_unit, Dbkc_print_unit
!
      real(dp), dimension(:), pointer  :: jLkc     => null()
      real(dp), dimension(:), pointer  :: iLkc     => null()
      real(dp), dimension(:), pointer  :: iLjc     => null()
!
      real(dp), dimension(:), pointer  :: Dbkc     => null()
!      
      real(dp)             :: alpha
!      
!     File names
      character(len=12)    :: Dbkc_name
      character(len=12)    :: jLkc_name
!
!      
!     Timing variables
      real     :: time_tot, time_occ, time_vir, time_order
      real     :: time_start, time_occ_read, time_vir_read
      real     :: time_abc, time_bac, time_acb, time_bca, time_cba, time_cab
      real     :: time_1, time_2
!      
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
      if(response_integer .eq. 0) then
!
!        Calculate [Ĥ,T_3] or [Ĥ,C_3]
!
         Dbkc_print_unit = Dbkc_g
         jLkc_print_unit = jLkc_g
!
      else if(response_integer .eq. 1) then
!
!        Calculate [[Ĥ,C_1],T_3] 
!
         Dbkc_print_unit = Dbkc_r
         jLkc_print_unit = jLkc_r
!
      else
!
         call quit('response_integer not 0 or 1')
!
      end if
!
      time_abc = 0.0
      time_bac = 0.0
      time_acb = 0.0
      time_bca = 0.0
      time_cba = 0.0
      time_cab = 0.0
!      
      time_vir = 0.0
      time_occ = 0.0
      time_vir_read = 0.0
      time_occ_read = 0.0
!      
      call cpu_time(time_start)
!
!
!     Calculate occupied active indices
!     ---------------------------------
      i_off = i + n_occ_inact
      j_off = j + n_occ_inact
      k_off = k + n_occ_inact
!
!      
!     Allocate arrays
!     ---------------
      call work_allocator(jLkc,n_vir_act*n_occ_gen)
      call work_allocator(iLkc,n_vir_act*n_occ_gen)
      call work_allocator(iLjc,n_vir_act*n_occ_gen)
!      
      call work_allocator(Dbkc,n_vir_aag)
!      
!
!     ------------------------
!     |Read (Db|kc) integrals|
!     ------------------------
!
      call cpu_time(time_1)
!
      read(Dbkc_print_unit,rec=k,iostat=ioerror) Dbkc
!         
      call cpu_time(time_2)
      time_vir_read = time_vir_read + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading Dbkc file"
         write(lupri,*) "Record k:", k
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in Omega 2')
      end if
!
!     ------------------------
!     |Read (jL|kc) integrals|
!     ------------------------
!
      rec_number = n_occ_act*(k-1) + j
!      
      call cpu_time(time_1)
!      
      read(jLkc_print_unit,rec=rec_number,iostat=ioerror) jLkc
!         
      call cpu_time(time_2)
      time_occ_read = time_occ_read + time_2 - time_1
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error reading jLkc file"
         write(lupri,*) "Record k, j, rec_number:", k, j, rec_number
         write(lupri,*) "Error code: ", ioerror
         call quit('Error reading integral file in Omega G2')
      end if
!
      if (i .ne. j) then
!         
!        ------------------------
!        |Read (iL|kc) integrals|
!        ------------------------
!
         rec_number = n_occ_act*(k-1) + i
!      
         call cpu_time(time_1)
!      
         call cpu_time(time_2)
         time_occ_read = time_occ_read + time_2 - time_1
!
         read(jLkc_print_unit,rec=rec_number,iostat=ioerror) iLkc
!         
         if(ioerror .ne. 0) then
            write(lupri,*) "Error reading jLkc file"
            write(lupri,*) "Record k, i, rec_number:", k, i, rec_number
            write(lupri,*) "Error code: ", ioerror
            call quit('Error reading integral file in Omega G2')
         end if
!
         if (i .ne. k .and. j .ne. k) then
!         
!           ------------------------
!           |Read (iL|jc) integrals|
!           ------------------------
!
            rec_number = n_occ_act*(j-1) + i
!      
            call cpu_time(time_1)
!      
            read(jLkc_print_unit,rec=rec_number,iostat=ioerror) iLjc
!         
            call cpu_time(time_2)
            time_occ_read = time_occ_read + time_2 - time_1
!
            if(ioerror .ne. 0) then
               write(lupri,*) "Error reading jLkc file"
               write(lupri,*) "Record i, j, rec_number:", i, j, rec_number
               write(lupri,*) "Error code: ", ioerror
               call quit('Error reading integral file in Omega G2')
            end if
!
         end if
!
      end if
!
      if (i .eq. j) then
         alpha = half
      else
         alpha = one
      end if
!      
!     -------------------------------------------
!     |Construct u^abc = 2t^abc - t^acb - t^cba|
!     -------------------------------------------
!      
      call cpu_time(time_1)
!
      call u_abc_calc(t_abc,u_abc,n_vir_act)
!
      call cpu_time(time_abc)
      time_abc = time_abc - time_1
!
!
!     -----------------------------------------------
!     |Calculate u^abc*(Db|kc) and add to Omega_aDij|
!     -----------------------------------------------
!
      omega_off = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(i-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
     &           Dbkc,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
!      
      call cpu_time(time_2)
      time_vir = time_vir + time_2 - time_1
!
!      
      if (i .ne. j) then
!         
!        -----------------------------------------------
!        |Calculate u^abc*(iL|kc) and add to Omega_abLj|
!        -----------------------------------------------
!
         omega_off = (n_occ_gen*n_vir_act**2)*(j-1) + 1
!
         call cpu_time(time_1)
!
         call dgemm('N','N',n_vir_act**2,n_occ_gen,n_vir_act,-one,u_abc,n_vir_act**2, &
     &              iLkc,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
!      
         call cpu_time(time_2)
         time_occ = time_occ + time_2 - time_1
!     
!      
      end if
!            
!     ------------------------------------------
!     |Construct u^bac = 2t^bac - t^bca - t^cab|
!     ------------------------------------------
!
      call cpu_time(time_1)
!
      call u_bac_calc(t_abc,u_abc,n_vir_act)
!
      call cpu_time(time_bac)
      time_bac = time_bac - time_1
!     
!     -----------------------------------------------
!     |Calculate u^bac*(Db|kc) and add to Omega_aDji|
!     -----------------------------------------------
!
      omega_off = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
     &           Dbkc,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
!      
      call cpu_time(time_2)
      time_vir = time_vir + time_2 - time_1
!
!      
!     -----------------------------------------------
!     |Calculate u^bac*(jL|kc) and add to Omega_abLi|
!     -----------------------------------------------
!
      call cpu_time(time_1)
!
      omega_off = (n_occ_gen*n_vir_act**2)*(i-1) + 1
!
      call dgemm('N','N',n_vir_act**2,n_occ_gen,n_vir_act,-one,u_abc,n_vir_act**2, &
     &           jLkc,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
!      
      call cpu_time(time_2)
      time_occ = time_occ + time_2 - time_1
!
!      
      if(j .ne. k) then
!
!        ------------------------------------------
!        |Construct u^acb = 2t^acb - t^abc - t^cab|
!        ------------------------------------------
!
         call cpu_time(time_1)
!
         call u_acb_calc(t_abc,u_abc,n_vir_act)
!
         call cpu_time(time_acb)
         time_acb = time_acb - time_1
!     
!            
!        ------------------------
!        |Read (Db|jc) integrals|
!        ------------------------
!
         call cpu_time(time_1)
!
         read(Dbkc_print_unit,rec=j,iostat=ioerror) Dbkc
!         
         call cpu_time(time_2)
         time_vir_read = time_vir_read + time_2 - time_1
!     
         if(ioerror .ne. 0) then
            write(lupri,*) "Error reading Dbkc file"
            write(lupri,*) "Record j:", j
            write(lupri,*) "Error code: ", ioerror
            call quit('Error reading integral file in Omega H2')
         end if
!
         if (i .eq. k) then
            alpha = half
         else
            alpha = one
         end if
!      
!        ------------------------
!        |Read (kL|jc) integrals|
!        ------------------------
!
         rec_number = n_occ_act*(j-1) + k
!   
         call cpu_time(time_1)
!
         read(jLkc_print_unit,rec=rec_number,iostat=ioerror) jLkc
!      
         call cpu_time(time_2)
         time_occ_read = time_occ_read + time_2 - time_1
!     
         if(ioerror .ne. 0) then
            write(lupri,*) "Error reading jLkc file"
            write(lupri,*) "Record j, k, rec_number:", j, k, rec_number
            write(lupri,*) "Error code: ", ioerror
            call quit('Error reading integral file in Omega G2')
         end if
!   
!
!        -----------------------------------------------
!        |Calculate u^acb*(Db|jc) and add to Omega_aDik|
!        -----------------------------------------------
!
         omega_off = (n_occ_act*n_vir_gen*n_vir_act)*(k-1) + (n_vir_gen*n_vir_act)*(i-1) + 1
!
         call cpu_time(time_1)
!
         call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
        &           Dbkc,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
!      
         call cpu_time(time_2)
         time_vir = time_vir + time_2 - time_1
!     
!     
         if(i .ne. k .and. i .ne. j) then
!      
!           -----------------------------------------------
!           |Calculate u^acb*(iL|jc) and add to Omega_abLk|
!           -----------------------------------------------
!
            omega_off = (n_occ_gen*n_vir_act**2)*(k-1) + 1
!
            call cpu_time(time_1)
!
            call dgemm('N','N',n_vir_act**2,n_occ_gen,n_vir_act,-one,u_abc,n_vir_act**2, &
           &           iLjc,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
!      
            call cpu_time(time_2)
            time_occ = time_occ + time_2 - time_1
!     
!      
         end if
!              
!        ------------------------------------------
!        |Construct u^bca = 2t^bca - t^bac - t^cba|
!        ------------------------------------------
!
         call cpu_time(time_1)
!
         call u_bca_calc(t_abc,u_abc,n_vir_act)
!
         call cpu_time(time_bca)
         time_bca = time_bca - time_1
!     
!         
!        -----------------------------------------------
!        |Calculate u^bca*(Db|jc) and add to Omega_aDki|
!        -----------------------------------------------
!
         omega_off = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(k-1) + 1
!
         call cpu_time(time_1)
!
         call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
        &           Dbkc,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
!      
         call cpu_time(time_2)
         time_vir = time_vir + time_2 - time_1
!     
!     
!        -----------------------------------------------
!        |Calculate u^bca*(kL|jc) and add to Omega_abLi|
!        -----------------------------------------------
!
         omega_off = (n_occ_gen*n_vir_act**2)*(i-1) + 1
!
         call cpu_time(time_1)
!
         call dgemm('N','N',n_vir_act**2,n_occ_gen,n_vir_act,-one,u_abc,n_vir_act**2, &
     &              jLkc,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
!   
         call cpu_time(time_2)
         time_occ = time_occ + time_2 - time_1
!     
!   
      end if
!      
      if(i .ne. k .and. (j .ne. k .or. i .ne. j)) then
!
!        ------------------------------------------
!        |Construct u^cba = 2t^cba - t^bca - t^abc|
!        ------------------------------------------
!      
         call cpu_time(time_1)
!
         call u_cba_calc(t_abc,u_abc,n_vir_act)
!
         call cpu_time(time_cba)
         time_cba = time_cba - time_1
!     
!      
         if (i .ne. j) then
!         
!
!           ------------------------
!           |Read (Db|ic) integrals|
!           ------------------------
!         
            call cpu_time(time_1)
!
            read(Dbkc_print_unit,rec=i,iostat=ioerror) Dbkc
!         
            call cpu_time(time_2)
            time_vir_read = time_vir_read + time_2 - time_1
!     
            if(ioerror .ne. 0) then
               write(lupri,*) "Error reading Dbkc file"
               write(lupri,*) "Record i:", i
               write(lupri,*) "Error code: ", ioerror
               call quit('Error reading integral file in Omega H2')
            end if
!         
            if (k .eq. j) then
               alpha = half
            else
               alpha = one
            end if
!      
!           ------------------------
!           |Read (kL|ic) integrals|
!           ------------------------
!
            rec_number = n_occ_act*(i-1) + k
!      
            call cpu_time(time_1)
!
            read(jLkc_print_unit,rec=rec_number,iostat=ioerror) iLkc
!        
            call cpu_time(time_2)
            time_occ_read = time_occ_read + time_2 - time_1
!     
            if(ioerror .ne. 0) then
               write(lupri,*) "Error reading jLkc file"
               write(lupri,*) "Record i, k, rec_number:", i, k, rec_number
               write(lupri,*) "Error code: ", ioerror
               call quit('Error reading integral file in Omega G2')
            end if
!      
!
!           -----------------------------------------------
!           |Calculate u^cba*(Db|ic) and add to Omega_aDkj|
!           -----------------------------------------------
!      
            omega_off = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(k-1) + 1
!
            call cpu_time(time_1)
!
            call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
        &              Dbkc,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
!      
            call cpu_time(time_2)
            time_vir = time_vir + time_2 - time_1
!     
!     
!           -----------------------------------------------
!           |Calculate u^cba*(kL|ic) and add to Omega_abLj|
!           -----------------------------------------------
!
            omega_off = (n_occ_gen*n_vir_act**2)*(j-1) + 1
!
            call cpu_time(time_1)
!     
            call dgemm('N','N',n_vir_act**2,n_occ_gen,n_vir_act,-one,u_abc,n_vir_act**2, &
           &              iLkc,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
!      
            call cpu_time(time_2)
            time_occ = time_occ + time_2 - time_1
!     
!      
         end if
!      
!        ------------------------------------------
!        |Construct u^cab = 2t^cab - t^acb - t^bac|
!        ------------------------------------------
!      
         call cpu_time(time_1)
!
         call u_cab_calc(t_abc,u_abc,n_vir_act)
!
         call cpu_time(time_cab)
         time_cab = time_cab - time_1
!     
!      
         if (j .ne. k) then
!      
!           ------------------------
!           |Read (jL|ic) integrals|
!           ------------------------
!
            rec_number = n_occ_act*(i-1) + j
!      
            call cpu_time(time_1)
!
            read(jLkc_print_unit,rec=rec_number,iostat=ioerror) iLjc
!        
            call cpu_time(time_2)
            time_occ_read = time_occ_read + time_2 - time_1
!     
            if(ioerror .ne. 0) then
               write(lupri,*) "Error reading jLkc file"
               write(lupri,*) "Record i, j, rec_number:", i, j, rec_number
               write(lupri,*) "Error code: ", ioerror
               call quit('Error reading integral file in Omega G2')
            end if
!      
!           -----------------------------------------------
!           |Calculate u^cab*(jL|ic) and add to Omega_abLk|
!           -----------------------------------------------
!      
            omega_off = (n_occ_gen*n_vir_act**2)*(k-1) + 1
!
            call cpu_time(time_1)
!
            call dgemm('N','N',n_vir_act**2,n_occ_gen,n_vir_act,-one,u_abc,n_vir_act**2, &
        &              iLjc,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
!      
            call cpu_time(time_2)
            time_occ = time_occ + time_2 - time_1
!     
!      
         end if
!         
!         
         if (i .ne. j) then
!      
!
!           -----------------------------------------------
!           |Calculate u^cab*(Db|ic) and add to Omega_aDjk|
!           -----------------------------------------------
!      
            omega_off = (n_occ_act*n_vir_gen*n_vir_act)*(k-1) + (n_vir_gen*n_vir_act)*(j-1) + 1
!
            call cpu_time(time_1)
!
            call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
        &              Dbkc,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
!      
            call cpu_time(time_2)
            time_vir = time_vir + time_2 - time_1
!     
!     
         end if
!      
!         
      end if
!
!
!     Deallocate arrays
!     -----------------
!
      call work_deallocator(Dbkc)
!
      call work_deallocator(iLjc)
      call work_deallocator(iLkc)
      call work_deallocator(jLkc)
!      
!     
!     Total timing
      call cpu_time(time_tot)
!      
      time_tot = time_tot - time_start
!      
      time_order = time_abc + time_bac + time_acb + time_bca + time_cba + time_cab
!      
!
!     Print timings
!     -------------
!
      if (print_mlcc3 .ge. 5) then
!
         write(lupri,*)
         write(lupri,*)
         write(lupri,*) 'Timings from mlcc3_omega_2'
         write(lupri,9999) 'Total', time_tot
         write(lupri,9999) 'read vir integral', time_vir_read
         write(lupri,9999) 'read occ integral', time_occ_read
         write(lupri,9999) 'MxM vir', time_vir
         write(lupri,9999) 'MxM occ', time_occ
         write(lupri,9999) 'reorder u^abc', time_abc
         write(lupri,9999) 'reorder u^bac', time_bac
         write(lupri,9999) 'reorder u^acb', time_acb
         write(lupri,9999) 'reorder u^bca', time_bca
         write(lupri,9999) 'reorder u^cba', time_cba
         write(lupri,9999) 'reorder u^cab', time_cab
         write(lupri,9999) 'total reorder', time_order
         write(lupri,*)
         write(lupri,*)
!
      end if
!
   end subroutine mlcc3_omega_2
!
!
   subroutine mlcc3_omega_add
!
!  Omega adder
!  Authors: Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: Add temporary, unpacked, a,b,i,j ordered Omega 2 vector
!           to the real packed Omega 2 vector and symmetrise in the
!           process
!
      implicit none
!
      integer  :: i, j, a, b
      integer  :: i_off, i_off2, j_off, j_off2, ai, bj
      integer  :: a_end, abij_n, abij_t, aibj_p
!
!
      if(mlcc3_active) then
!
!        Virtual omega vector
!        --------------------
!
         do j = 1,n_occ_act
            do b = 1,n_vir_gen
               do i = 1,j
!
                  if(i .eq. j) then
                     a_end = min(b,n_vir_act)
                  else if(b .gt. n_vir_act) then
                     a_end = n_vir_act
                  else
                     a_end = n_vir_gen
                  end if
!
                  do a = 1,a_end
!
                     i_off = i + n_occ_inact
                     j_off = j + n_occ_inact
!
                     ai = n_vir*(i_off-1) + a
                     bj = n_vir*(j_off-1) + b
!
                     aibj_p = max(ai,bj)*(max(ai,bj)-3)/2 + ai + bj
!
                     if(a .le. n_vir_act) then
!
                        abij_n = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(i-1) &
                       &       + n_vir_act*(b-1) + a
!
                        omega2(aibj_p) = omega2(aibj_p) + omega_v(abij_n)
!
                     end if
!
                     if(b .le. n_vir_act) then
!
                        abij_t = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(j-1) &
                       &       + n_vir_act*(a-1) + b
!
                        omega2(aibj_p) = omega2(aibj_p) + omega_v(abij_t)
!
                     end if
!
                  end do
               end do
            end do
         end do
!
!
!        Occupied omega vector
!        ---------------------
!
         do j = 1,n_occ_act
!
            j_off = j + n_occ_inact
!
            do b = 1,n_vir_act
               do i = 1,j_off - n_gen_inact
!
                  i_off = i + n_gen_inact
!
                  if(i_off .eq. j_off) then
                     a_end = b
                  else
                     a_end = n_vir_act
                  end if
!
                  do a = 1,a_end
!
                     ai = n_vir*(i_off-1) + a
                     bj = n_vir*(j_off-1) + b
!
                     aibj_p = max(ai,bj)*(max(ai,bj)-3)/2 + ai + bj
!
                     abij_n = (n_occ_gen*n_vir_act**2)*(j-1) + (n_vir_act**2)*(i-1) &
                    &       + n_vir_act*(b-1) + a
!
                     omega2(aibj_p) = omega2(aibj_p) + omega_o(abij_n)
!
                     if(i_off .gt. n_occ_inact) then
!
                        i_off2 = i_off - n_occ_inact
                        j_off2 = j_off - n_gen_inact
!
                        abij_t = (n_occ_gen*n_vir_act**2)*(i_off2-1) + (n_vir_act**2)*(j_off2-1) &
                       &       + n_vir_act*(a-1) + b
!
                        omega2(aibj_p) = omega2(aibj_p) + omega_o(abij_t)
!
                     end if
!
!
                  end do
               end do
            end do
         end do
!
!
      else
!
!
!        Only one omega vector to symmtrise and add
!        ------------------------------------------
!
         aibj_p = 0
!
         do j = 1,n_occ
            do b = 1,n_vir
               do i = 1,j
!
                  if(i .eq. j) then
                     a_end = b
                  else
                     a_end = n_vir
                  end if
!
                  do a = 1,a_end
!
                     aibj_p = aibj_p + 1
!
                     abij_n = (n_occ*n_vir**2)*(j-1) + (n_vir**2)*(i-1) + n_vir*(b-1) + a
                     abij_t = (n_occ*n_vir**2)*(i-1) + (n_vir**2)*(j-1) + n_vir*(a-1) + b
!
                     omega2(aibj_p) = omega2(aibj_p) + (omega_v(abij_n) + omega_v(abij_t))
!
                  end do
               end do
            end do
         end do
!
      end if
!
   end subroutine mlcc3_omega_add
!
!
   subroutine mlcc3_t2_square
!
!  T2 square
!  Authors: Henrik Koch and Rolf H. Myhre
!  March 2015
!
!  Purpose: Square up T_2 amplitudes in a,b,j,i order for the
!           matrix - matrix multiplication. Two different arrays
!           in case of active spaces
!
      implicit none
!
      integer  :: i, j, a, b
      integer  :: i_off, i_off2, j_off, j_off2, ai, bj
      integer  :: a_end, abji_n, abji_t, aibj_p
!
      real(dp), dimension(:), pointer  :: r2_aDji  => null()
      real(dp), dimension(:), pointer  :: r2_abLi  => null()
!
      real(dp), dimension(:), pointer  :: r2am  => null()
!
!
      if(response_integer .eq. 0) then
!
         r2_aDji  => t2_aDji
         r2_abLi  => t2_abLi
         r2am     => t2am
!
      else if(response_integer .eq. 1) then
!
         r2_aDji  => c2_aDji
         r2_abLi  => c2_abLi
         r2am     => c2am
!
      else
!
         call quit('response_integer not 0 or 1')
!
      end if
!
      if(mlcc3_active) then
!
!
!        T^aD_ij array
!        -------------
!
         do i = 1,n_occ_act
            do j = 1,i
               do b = 1,n_vir_gen
!
                  if(i .eq. j) then
                     a_end = min(b,n_vir_act)
                  else if(b .gt. n_vir_act) then
                     a_end = n_vir_act
                  else
                     a_end = n_vir_gen
                  end if
!
                  do a = 1,a_end
!
                     i_off = i + n_occ_inact
                     j_off = j + n_occ_inact
!
                     ai = n_vir*(i_off-1) + a
                     bj = n_vir*(j_off-1) + b
!
                     aibj_p = max(ai,bj)*(max(ai,bj)-3)/2 + ai + bj
!
                     if(a .le. n_vir_act) then
!
                        abji_n = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(j-1) &
                       &       + n_vir_act*(b-1) + a
!
                        r2_aDji(abji_n) = r2am(aibj_p)
!
                     end if
!
                     if(ai .ne. bj .and. b .le. n_vir_act) then
!
                        abji_t = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(i-1) &
                       &       + n_vir_act*(a-1) + b
!
                        r2_aDji(abji_t) = r2am(aibj_p)
!
                     end if
!
                  end do
               end do
            end do
         end do
!
!
!        T^ab_iL array
!        -------------
!
         do i = 1,n_occ_act
!
            i_off = i + n_occ_inact
!
            do j = 1,i_off - n_gen_inact
!            
               j_off = j + n_gen_inact
!               
               do b = 1,n_vir_act
!
                  if(i_off .eq. j_off) then
                     a_end = b
                  else
                     a_end = n_vir_act
                  end if
!
                  do a = 1,a_end
!
                     ai = n_vir*(i_off-1) + a
                     bj = n_vir*(j_off-1) + b
!
                     aibj_p = max(ai,bj)*(max(ai,bj)-3)/2 + ai + bj
!
                     abji_n = (n_occ_gen*n_vir_act**2)*(i-1) + (n_vir_act**2)*(j-1) &
                    &       + n_vir_act*(b-1) + a
!
                     r2_abLi(abji_n) = r2am(aibj_p)
!
                     if(ai .ne. bj .and. j_off .gt. n_occ_inact) then
!
                        i_off2 = i_off - n_gen_inact
                        j_off2 = j_off - n_occ_inact
!
                        abji_t = (n_occ_gen*n_vir_act**2)*(j_off2-1) + (n_vir_act**2)*(i_off2-1) &
                       &       + n_vir_act*(a-1) + b
!
                        r2_abLi(abji_t) = r2am(aibj_p)
!
                     end if
!
                  end do
               end do
            end do
         end do
!
      else
!
!
!        Only one array
!        --------------
!
         do i = 1,n_occ
            do j = 1,i
               do b = 1,n_vir
!
                  if(i .eq. j) then
                     a_end = b
                  else
                     a_end = n_vir
                  end if
!
                  do a = 1,a_end
!
                     ai = n_vir*(i-1) + a
                     bj = n_vir*(j-1) + b
!
                     aibj_p = max(ai,bj)*(max(ai,bj)-3)/2 + ai + bj
!
                     abji_n = (n_occ*n_vir**2)*(i-1) + (n_vir**2)*(j-1) &
                    &       + n_vir*(b-1) + a
!
                     r2_aDji(abji_n) = r2am(aibj_p)
!
                     if(ai .ne. bj) then
                        abji_t = (n_occ*n_vir**2)*(j-1) + (n_vir**2)*(i-1) &
                       &       + n_vir*(a-1) + b
!
                        r2_aDji(abji_t) = r2am(aibj_p)
                     end if
!
                  end do
               end do
            end do
         end do
!
      end if
!
   end subroutine mlcc3_t2_square
!
!
   subroutine mlcc3_file_opener(io_opt)
!
!  File opener
!  Authors: Henrik Koch and Rolf H. Myhre
!  April 2015
!
!  Purpose: Open integral files required in omega
!           or rho calculations
!
      implicit none
!
      integer, intent(in)  :: io_opt
      integer              :: ioerror, n_occ_read
!
      n_occ_read  = n_occ_gen*n_vir_act
!
!
      if(io_opt .eq. 0) then
!
!        Ground state energy. Only standard integrals required for T3
!        ------------------------------------------------------------
!
         bDck_g = new_unit()
!      
         open(unit=bDck_g, file=bDck_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening bDck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
         Ljck_g = new_unit()
!         
         open(unit=Ljck_g, file=Ljck_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_occ_read, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Ljck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
!
!        Integrals for Omega1 calculations
!        ---------------------------------
!
         jbkc_g = new_unit()
!         
         open(unit=jbkc_g, file=jbkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_act**2, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
!        Integrals used in Omega2 calculation
!        ------------------------------------
!
         Dbkc_g = new_unit()
!      
         open(unit=Dbkc_g, file=Dbkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Dbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
         jLkc_g = new_unit()
!      
         open(unit=jLkc_g, file=jLkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_act*n_occ_gen, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jLkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
!
      else if(io_opt .eq. 1) then
!
!        Response call type 1. Both standard and response integrals required
!        in C3 calculation
!        -------------------------------------------------------------------
!
         bDck_g = new_unit()
!      
         open(unit=bDck_g, file=bDck_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening bDck file"
            write(lupri,*) "Error code: ", ioerror
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
         Ljck_g = new_unit()
!         
         open(unit=Ljck_g, file=Ljck_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_occ_read, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Ljck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
         bDck_r = new_unit()
!      
         open(unit=bDck_r, file=bDck_resp_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening bDck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
         Ljck_r = new_unit()
!         
         open(unit=Ljck_r, file=Ljck_resp_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_occ_read, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Ljck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
!
!        Integrals for Rho1 calculations
!        -------------------------------
!
         jbkc_g = new_unit()
!         
         open(unit=jbkc_g, file=jbkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_act**2, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
!        Integrals used in Rho2 calculation
!        ----------------------------------
!
         Dbkc_g = new_unit()
!      
         open(unit=Dbkc_g, file=Dbkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Dbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
         jLkc_g = new_unit()
!      
         open(unit=jLkc_g, file=jLkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_act*n_occ_gen, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jLkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
!
      else if(io_opt .eq. 2) then
!
!        Response call type 2. Only standard integrals in T3 calculation
!        ---------------------------------------------------------------
!
         bDck_g = new_unit()
!      
         open(unit=bDck_g, file=bDck_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening bDck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
         Ljck_g = new_unit()
!         
         open(unit=Ljck_g, file=Ljck_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_occ_read, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Ljck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
!        Integrals used in Rho2 calculation
!        ----------------------------------
!
         Dbkc_r = new_unit()
!      
         open(unit=Dbkc_r, file=Dbkc_resp_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Dbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
         jLkc_r = new_unit()
!      
         open(unit=jLkc_r, file=jLkc_resp_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_act*n_occ_gen, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jLkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
!      
      end if
!
   end subroutine mlcc3_file_opener
!
   subroutine mlcc3_file_closer(io_opt)
!
!  File closer
!  Authors: Henrik Koch and Rolf H. Myhre
!  April 2015
!
!  Purpose: Close integral files required in omega
!           or rho calculations
!
      implicit none
!
      integer, intent(in)  :: io_opt
      integer              :: ioerror
!
!
      if(io_opt .eq. 0) then
!
!      
!        Ground state energy. Only standard integrals required
!        -----------------------------------------------------
         close(bDck_g,status='keep',iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing bDck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
!      
         close(Ljck_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing Ljck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
!      
!        Integrals for Omega1 calculations
!        ---------------------------------
!
         close(jbkc_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing jbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!            
!        Integrals used in Omega2 calculation
!        ------------------------------------
!
         close(Dbkc_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing Dbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!            
         close(jLkc_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing jLkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!
!
      else if(io_opt .eq. 1) then
!
!      
!        Response call type 1. Both standard and response integrals required
!        in R3 calculation
!        -------------------------------------------------------------------
!
         close(bDck_g,status='keep',iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing bDck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
!      
         close(Ljck_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing Ljck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
         close(bDck_r,status='keep',iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing bDck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
!      
         close(Ljck_r,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing Ljck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
!        Integrals for R1 calculations
!        -----------------------------
!
         close(jbkc_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing jbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!            
!        Integrals used in Rho2 calculation
!        ----------------------------------
!
         close(Dbkc_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing Dbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!            
         close(jLkc_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing jLkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!
!      
      else if(io_opt .eq. 2) then
!
!      
!        Response call type 2. Only standard integrals in T3 calculation
!        ---------------------------------------------------------------
!
         close(bDck_g,status='keep',iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing bDck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
!      
         close(Ljck_g,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing Ljck file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
!        Integrals used in Rho2 calculation
!        ----------------------------------
!
         close(Dbkc_r,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing Dbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!            
         close(jLkc_r,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing jLkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!
      end if
!
   end subroutine mlcc3_file_closer
!
end module mlcc3_omega
