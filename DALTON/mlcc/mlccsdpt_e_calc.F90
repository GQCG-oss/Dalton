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
module mlccsdpt_energy
!
!
!  Calculate the MLCCSD(T) energy contribution
!  Authors Henrik Koch and Rolf H. Myhre
!  May 2015
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
   real(dp), dimension(:), pointer, private  :: w_abc    => null()
   real(dp), dimension(:), pointer, private  :: v_abc    => null()
!
   real(dp), dimension(:), pointer, private  :: w_cba    => null()
!
!  Temporary T2 and C2 array with virtual general, a,D,j,i ordering
   real(dp), dimension(:), pointer, private  :: t2_aDji  => null()
!
!  Temporary T2 and C2 array with virtual general, a,b,L,i ordering
   real(dp), dimension(:), pointer, private  :: t2_abLi  => null()
!
!  Integral arrays
   real(dp), dimension(:), pointer, private  :: bDci_int  => null()
   real(dp), dimension(:), pointer, private  :: bDcj_int  => null()
   real(dp), dimension(:), pointer, private  :: bDck_int  => null()
!
   real(dp), dimension(:), pointer, private  :: Lmci_int  => null()
   real(dp), dimension(:), pointer, private  :: Lmcj_int  => null()
   real(dp), dimension(:), pointer, private  :: Lmck_int  => null()
!
   real(dp), dimension(:), pointer, private  :: iajb_int  => null()
   real(dp), dimension(:), pointer, private  :: iakc_int  => null()
!
!  Batch information
   integer, private  :: batch_size, n_batch
!      
!      
!  Print units
!  -----------
!
!  Ground state integrals
   integer, private  :: bDck_u, Ljck_u, jbkc_u
!
!      
contains
!
   subroutine mlccsdpt_e_calc(energy)
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
      real(dp), intent(inout) :: energy
!
      integer  :: i, j, k
      integer  :: i_batch, j_batch, k_batch
      integer  :: i_start, j_start, k_start
      integer  :: i_end  , j_end  , k_end  
      integer  :: i_size , j_size , k_size 
!      
      integer  :: r_start,r_end, n_occ_int, n_act_int
!      
      integer  :: a, b, abij_n, abij_t, ioerror
!
      real(dp) :: ddot
      real(dp) :: w_norm, v_norm, energy_ijk
!
!     Timing variables
      real     :: time_tot, time_w, time_v, time_e
      real     :: time_add, time_t2sq
      real     :: time_start, time_1, time_allo, time_deallo
      real     :: time_w_sum, time_v_sum, time_e_sum
      real     :: time_open, time_close, time_batch
      real     :: time_read,  time_2
!      
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
      energy = zero
!      
      time_w_sum     = 0.0
      time_v_sum     = 0.0
      time_e_sum     = 0.0
      time_read      = 0.0
      time_allo      = 0.0
!
      call cpu_time(time_start)
!
!
      n_occ_int = n_vir_act*n_occ_gen*n_occ_act
      n_act_int = n_vir_act**2*n_occ_act
!
!     Allocate arrays
!     If mlcc3_active, we need two temporary amplitude vectors due to
!     different lengths of n_vir and n_vir_act and n_occ and n_occ_act, 
!     -----------------------------------------------------------------
!
      call cpu_time(time_1)
!      
      if(mlcc3_active) then
!      
         call work_allocator(t2_aDji,n_vir_act*n_vir_gen*n_occ_act**2)
         call work_allocator(t2_abLi,n_vir_act**2*n_occ_gen*n_occ_act)
!      
      else
!      
         call work_allocator(t2_aDji,n_vir**2*n_occ**2)
         t2_abLi => t2_aDji
!      
      end if
!
      call work_allocator(w_abc,n_vir_act**3)
      call work_allocator(w_cba,n_vir_act**3)
      call work_allocator(v_abc,n_vir_act**3)
!      
!      
      call cpu_time(time_2)
      time_allo = time_allo + time_2 - time_1
!
!
!     Square up T_2 vectors
!     ---------------------
!
      call cpu_time(time_1)
!      
      call mlccsdpt_t2_square
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
      call mlccsdpt_file_opener
!
      call cpu_time(time_open)
      time_open = time_open - time_1
!
!
!     Set up batching
!     ---------------
!
      call cpu_time(time_1)
!      
      call mlccsdpt_batch_setup
!
      call cpu_time(time_batch)
      time_batch = time_batch - time_1
!
!
!     Allocate integral arrays
!     ------------------------
!
      call cpu_time(time_1)
!
      call work_allocator(bDci_int,batch_size*n_vir_aag)
      call work_allocator(Lmci_int,batch_size*n_occ_int)
!
      call work_allocator(iajb_int,batch_size*n_act_int)
!
      if(n_batch .ne. 1) then
         call work_allocator(bDcj_int,batch_size*n_vir_aag)
         call work_allocator(Lmcj_int,batch_size*n_occ_int)
!
         call work_allocator(bDck_int,batch_size*n_vir_aag)
         call work_allocator(Lmck_int,batch_size*n_occ_int)
!
         call work_allocator(iakc_int,batch_size*n_act_int)
!
      else
         bDcj_int => bDci_int
         Lmcj_int => Lmci_int
!
         bDck_int => bDci_int
         Lmck_int => Lmci_int
!
         iakc_int => iajb_int
      end if
!
      call cpu_time(time_2)
      time_allo = time_allo + time_2 - time_1
!
      if(print_mlcc3 .ge. 5) then
         write(lupri,*)
         write(lupri,*) 'work_free:    ', work_free()
         write(lupri,*)
      end if
!
!     Loop over batches
!     -----------------
!
      a = 0
      i_size = batch_size
!
      do i_batch = 1,n_batch
!
!        Find i batch parameters
!        -----------------------
!
         i_start = batch_size*(i_batch-1) + 1
         if(i_batch .eq. n_batch) then
            i_size = n_occ_act - batch_size*(n_batch-1)
         end if
         i_end   = i_start + i_size - 1
!
!
!        Read i integrals
!        ----------------
!
         call cpu_time(time_1)
!
         call mlccsdpt_int_reader(bDci_int,Lmci_int,i_start,i_size)
!
         call cpu_time(time_2)
         time_read = time_read + time_2 - time_1
!
!
         j_size = batch_size
!
         do j_batch = 1,i_batch
!
!           Find j batch parameters
!           -----------------------
!
            j_start = batch_size*(j_batch-1) + 1
            if(j_batch .eq. n_batch) then
               j_size = n_occ_act - batch_size*(n_batch-1)
            end if
            j_end   = j_start + j_size - 1
!
!
!           Read j integrals
!           ----------------
!
            call cpu_time(time_1)
!
            call mlccsdpt_jbkc_reader(iajb_int,j_start,j_size)
!
            call cpu_time(time_2)
            time_read = time_read + time_2 - time_1
!
            if(j_batch .ne. i_batch) then
               call cpu_time(time_1)
!
               call mlccsdpt_int_reader(bDcj_int,Lmcj_int,j_start,j_size)
!
               call cpu_time(time_2)
               time_read = time_read + time_2 - time_1
!
            end if
!
!
            k_size = batch_size
!
            do k_batch = 1,j_batch
!
!              Find k batch parameters
!              -----------------------
!
               k_start = batch_size*(k_batch-1) + 1
               if(k_batch .eq. n_batch) then
                  k_size = n_occ_act - batch_size*(n_batch-1)
               end if
               k_end   = k_start + k_size - 1
!
!
!              Read k integrals
!              ----------------
!
               if(k_batch .ne. j_batch) then
!
                  call cpu_time(time_1)
!
                  call mlccsdpt_jbkc_reader(iakc_int,k_start,k_size)
!
                  call cpu_time(time_2)
                  time_read = time_read + time_2 - time_1
!
                  if(k_batch .ne. i_batch) then
!
                     call cpu_time(time_1)
!
                     call mlccsdpt_int_reader(bDck_int,Lmck_int,k_start,k_size)
!
                     call cpu_time(time_2)
                     time_read = time_read + time_2 - time_1
!
                  end if
               end if
!
!
!              Restricted loop over occupied indices i >= j >= k
!              -------------------------------------------------
!
               do i = i_start,i_end
                  do j = j_start,min(j_end,i)
                     do k = k_start,min(k_end,j)
!            
                        if (i .eq. j .and. j .eq. k) then
                           cycle
                        end if
!
!                       W intermediates               
!                       ---------------
!                  
                        call cpu_time(time_1)
!               
                        call mlccsdpt_w_calc(i,j,k,i_batch,j_batch,k_batch)
!               
                        call cpu_time(time_w)
                        time_w = time_w - time_1
                        time_w_sum = time_w_sum + time_w
!               
                        if(print_mlcc3 .ge. 7) then
!
                           w_norm = ddot(n_vir_act**3,w_abc,1,w_abc,1)
!
                           write(lupri,*)
                           write(lupri,*) 'w_norm: ', w_norm
                           write(lupri,*)
!
                        end if
!
!               
!                       (V^abc-V^cba) intermediates 
!                       ---------------------------
!                  
                        call cpu_time(time_1)
!               
                        call mlccsdpt_v_calc(i,j,k,i_batch,j_batch,k_batch)
!               
                        call cpu_time(time_v)
                        time_v = time_v - time_1
                        time_v_sum = time_v_sum + time_v
!               
                        if(print_mlcc3 .ge. 7) then
!
                           v_norm = ddot(n_vir_act**3,v_abc,1,v_abc,1)
!
                           write(lupri,*)
                           write(lupri,*) 'v_norm: ', v_norm
                           write(lupri,*)
!
                        end if
!
!               
!                       Energy calculation
!                       ------------------
!                  
                        call cpu_time(time_1)
!               
                        call mlccsdpt_w_prod(energy_ijk,i,j,k)
!               
                        energy = energy + energy_ijk
!               
                        call cpu_time(time_e)
                        time_e = time_e - time_1
                        time_e_sum = time_e_sum + time_e
!               
!               
!                       Print timings
!                       -------------
!
                        if (print_mlcc3 .ge. 5) then
!
                           write(lupri,*)
                           write(lupri,*)
                           write(lupri,*) 'Timings from mlccsdpt_e_calc in loop'
                           write(lupri,*) 'i, j, k', i, j, k
                           write(lupri,9999) 'mlccsdpt_w_calc', time_w
                           write(lupri,9999) 'mlccsdpt_v_calc', time_v
                           write(lupri,9999) 'mlccsdpt_w_prod', time_e
                           write(lupri,*)
                           write(lupri,*)
!
                        end if
!
                     end do
                  end do
               end do !End restricted loop
!      
            end do
         end do
      end do !End batch loop
!      
!     Divide by 3
!     -----------
!
      energy = energy/three
!
!     Close files
!     ----------
!
      call cpu_time(time_1)
!      
      call mlccsdpt_file_closer
!
      call cpu_time(time_close)
      time_close = time_close - time_1
!
!
!     Deallocate integral arrays
!     --------------------------
!
      call cpu_time(time_1)
!      
      if(n_batch .ne. 1) then
!
         call work_deallocator(iakc_int)
!
         call work_deallocator(Lmck_int)
         call work_deallocator(bDck_int)
!
         call work_deallocator(Lmcj_int)
         call work_deallocator(bDcj_int)
!
      end if
!
      call work_deallocator(iajb_int)
!
      call work_deallocator(Lmci_int)
      call work_deallocator(bDci_int)
!
!     Deallocate arrays
!     -----------------
!
      call work_deallocator(v_abc)
      call work_deallocator(w_cba)
      call work_deallocator(w_abc)
!
!
      if(mlcc3_active) then
!
         call work_deallocator(t2_abLi)
         call work_deallocator(t2_aDji)
!      
      else
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
         write(lupri,*) 'Timings from mlccsdpt_e_calc after loop'
         write(lupri,9999) 'Total', time_tot
         write(lupri,9999) 'square T_2', time_t2sq
         write(lupri,9999) 'set up batching', time_batch
         write(lupri,9999) 'allocation', time_allo
         write(lupri,9999) 'open files', time_open
         write(lupri,9999) 'read integrals', time_read
         write(lupri,9999) 'mlccsdpt_w_calc', time_w_sum
         write(lupri,9999) 'mlccsdpt_v_calc', time_v_sum
         write(lupri,9999) 'mlccsdpt_w_prod', time_e_sum
         write(lupri,9999) 'close files', time_close
         write(lupri,9999) 'deallocation', time_deallo
         write(lupri,*)
         write(lupri,*)
!
      end if
!
   end subroutine mlccsdpt_e_calc
!      
!
   subroutine mlccsdpt_t2_square
!
!  T2 square
!  Authors: Henrik Koch and Rolf H. Myhre
!  June 2015
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
!
      r2_aDji  => t2_aDji
      r2_abLi  => t2_abLi
      r2am     => t2am
!
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
   end subroutine mlccsdpt_t2_square
!
!
!      
   subroutine mlccsdpt_file_opener
!
!  File opener
!  Authors: Henrik Koch and Rolf H. Myhre
!  June 2015
!
!  Purpose: Open integral files required in MLCCSD(T)
!
      implicit none
!
      integer              :: ioerror, n_occ_read
!
      n_occ_read  = n_occ_gen*n_vir_act
!
!   
!     bDck interals
!     -------------
!
      bDck_u = new_unit()
!   
      open(unit=bDck_u, file=bDck_file_name, action='read', status='old', &
     &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening bDck file"
         write(lupri,*) "Error code:   ", ioerror
         call quit('Error in file opener')
      end if
!   
!   
!     Ljck interals
!     -------------
!
      Ljck_u = new_unit()
!      
      open(unit=Ljck_u, file=Ljck_file_name, action='read', status='old', &
     &     access='direct', form='unformatted', recl=dp*n_occ_read, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening Ljck_g file"
         write(lupri,*) "Error code:   ", ioerror
         call quit('Error in file opener')
      end if
!      
!      
!     jbkc interals
!     -------------
!
      jbkc_u = new_unit()
!      
      open(unit=jbkc_u, file=jbkc_file_name, action='read', status='old', &
     &     access='direct', form='unformatted', recl=dp*n_vir_act**2, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening jbkc file"
         write(lupri,*) "Error code:   ", ioerror
         call quit('Error in file opener')
      end if
!
!      
   end subroutine mlccsdpt_file_opener
!
   subroutine mlccsdpt_jbkc_reader(mbxc,x_start,x_size)
!
!  MLCCSD(T) Integral reader
!  Authors: Henrik Koch and Rolf H. Myhre
!  June 2015
!
!  Purpose: Read integrals in batch
!
      implicit none
!
      integer, intent(in)  :: x_start,x_size
!
      real(dp), dimension(:), pointer, intent(inout)  :: mbxc  
!
      integer  :: r_start,r_end,ioerror,x,y,n_act_int,record
!
      n_act_int = n_vir_act**2
!         
!     mbxc integrals
      do x = 1,x_size
         do y = 1,n_occ_act
!
            r_start  = n_act_int*(n_occ_act*(x-1)+y-1)+1
            r_end    = r_start + n_act_int -1
            record   = n_occ_act*(x+x_start-2) + y
!
            read(jbkc_u,rec=record,iostat=ioerror) mbxc(r_start:r_end)
!     
            if(ioerror .ne. 0) then
               write(lupri,*) "Error reading mbxc file"
               write(lupri,*) "Error code: ", ioerror
               call quit('Error reading integral file')
            end if
!     
         end do
      end do
!         
!         
   end subroutine mlccsdpt_jbkc_reader
!
   subroutine mlccsdpt_int_reader(bDcx,Lmcx,x_start,x_size)
!
!  MLCCSD(T) Integral reader
!  Authors: Henrik Koch and Rolf H. Myhre
!  June 2015
!
!  Purpose: Read integrals in batch
!
      implicit none
!
      integer, intent(in)  :: x_start,x_size
!
      real(dp), dimension(:), pointer, intent(inout)  :: bDcx
      real(dp), dimension(:), pointer, intent(inout)  :: Lmcx
!
      integer  :: r_start,r_end,ioerror,x,y,n_occ_int,record
!
      n_occ_int = n_vir_act*n_occ_gen
!
!     bDcx integrals
      do x = 1,x_size
!
         r_start  = n_vir_aag*(x-1)+1
         r_end    = r_start + n_vir_aag -1
         read(bDck_u,rec=x+x_start-1,iostat=ioerror) bDcx(r_start:r_end)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error reading bDcx_g file"
            write(lupri,*) "Error code: ", ioerror
            call quit('Error reading integral file')
         end if
!
      end do
!
!     Lmcx integrals
      do x = 1,x_size
         do y = 1,n_occ_act
!      
            r_start  = n_occ_int*(n_occ_act*(x-1)+y-1)+1
            r_end    = r_start + n_occ_int -1
            record   = n_occ_act*(x+x_start-2) + y
            read(Ljck_u,rec=record,iostat=ioerror) Lmcx(r_start:r_end)
!      
            if(ioerror .ne. 0) then
               write(lupri,*) "Error reading Lmcx_g file"
               write(lupri,*) "Error code: ", ioerror
               call quit('Error reading integral file')
            end if
!      
         end do
      end do
!         
!         
   end subroutine mlccsdpt_int_reader
!
   subroutine mlccsdpt_file_closer
!
!  File closer
!  Authors: Henrik Koch and Rolf H. Myhre
!  June 2015
!
!  Purpose: Close integral files
!
      implicit none
!
      integer              :: ioerror
!
!
!     bDck integrals
!     --------------
!
      close(bDck_u,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing bDck file"
         write(lupri,*) "Error code:   ", ioerror
         call quit('Error in file closer')
      end if
!   
!
!     Ljck integrals
!     --------------
!
      close(Ljck_u,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing Ljck file"
         write(lupri,*) "Error code:   ", ioerror
         call quit('Error in file closer')
      end if
!      
!      
!     jbkc integrals
!     --------------
!
      close(jbkc_u,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing jbkc file"
         write(lupri,*) "Error code:   ", ioerror
         call quit('Error in file closer')
      end if
!            
!
   end subroutine mlccsdpt_file_closer
!
!
   subroutine mlccsdpt_batch_setup
!
!  Batch setup
!  Authors: Henrik Koch and Rolf H. Myhre
!  June 2015
!
!  Purpose: Set up batching over integrals in the occupied loop
!
      implicit none
!
      integer  :: w_req, w_avail, n_occ_int, n_act_int
!
!
      w_avail = work_free()
!
      n_occ_int = n_vir_act*n_occ_gen*n_occ_act
      n_act_int = n_vir_act**2*n_occ_act
!
      w_req = n_vir_aag + n_occ_int + n_act_int
!
!
      if(w_avail .gt. n_occ_act*w_req) then !all integrals in memory
!
         if(print_mlcc3 .ge. 2) then
!
            write(lupri,*)
            write(lupri,*) 'All integrals in memory'
!
         end if
!
         n_batch     = 1
         batch_size  = n_occ_act
!
      else if(w_avail .lt. 3*w_req) then !not enough memory
!
         write(lupri,*)
         write(lupri,*) 'Not enough memory for a batch in batch_setup'
         write(lupri,*) 'w_avail: ', w_avail
         write(lupri,*) 'w_req:   ', w_req
         write(lupri,*)
!
         call quit('Not enough memory')
!
      else !batching required
!
         if(print_mlcc3 .ge. 2) then
!
            write(lupri,*)
            write(lupri,*) 'Batching required'
!
         end if
!
         batch_size  = w_avail/(3*w_req-n_act_int)
         n_batch     = (n_occ_act-1)/batch_size + 1
!
      end if
!
      if(print_mlcc3 .ge. 4) then
         write(lupri,*)
         write(lupri,*) 'Batch information from batch_setup'
         write(lupri,*) 'w_avail:    ', w_avail
         write(lupri,*) 'w_req:      ', w_req
         write(lupri,*) '3*w_req:    ', 3*w_req
         write(lupri,*) 'batch_size: ', batch_size
         write(lupri,*) 'n_batch:    ', n_batch
         write(lupri,*)
      end if
!
!
   end subroutine mlccsdpt_batch_setup
!
!
   subroutine mlccsdpt_w_calc(i,j,k,i_batch,j_batch,k_batch)
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
      integer, intent(in)  :: i_batch, j_batch, k_batch
!      
      integer              :: i_start, j_start, k_start
!      
      integer              :: i_int_v, j_int_v, k_int_v
!      
      integer              :: jk_int_o, ik_int_o, kj_int_o
      integer              :: ij_int_o, ki_int_o, ji_int_o
!      
      integer              :: r2_off, n_occ_int
!
      real(dp)             :: beta
!
!     (bD|ck) and (Lj|ck) integrals
      real(dp), dimension(:), pointer  :: bDci     => null()
      real(dp), dimension(:), pointer  :: bDcj     => null()
      real(dp), dimension(:), pointer  :: bDck     => null()
!
      real(dp), dimension(:), pointer  :: Ljck     => null()
      real(dp), dimension(:), pointer  :: Lick     => null()
      real(dp), dimension(:), pointer  :: Lkbj     => null()
      real(dp), dimension(:), pointer  :: Libj     => null()
      real(dp), dimension(:), pointer  :: Lkai     => null()
      real(dp), dimension(:), pointer  :: Ljai     => null()
!      
!      
!     Timing variables
      real     :: time_tot, time_dgvir, time_dgocc
      real     :: time_start, time_order
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
!
      n_occ_int = n_vir_act*n_occ_gen
!
      i_start = batch_size*(i_batch-1) + 1
      j_start = batch_size*(j_batch-1) + 1
      k_start = batch_size*(k_batch-1) + 1
!
      i_int_v = n_vir_aag*(i-i_start) + 1
      j_int_v = n_vir_aag*(j-j_start) + 1
      k_int_v = n_vir_aag*(k-k_start) + 1
!
      jk_int_o = n_occ_int*(n_occ_act*(k-k_start)+j-1) + 1
      ik_int_o = n_occ_int*(n_occ_act*(k-k_start)+i-1) + 1
      kj_int_o = n_occ_int*(n_occ_act*(j-j_start)+k-1) + 1
      ij_int_o = n_occ_int*(n_occ_act*(j-j_start)+i-1) + 1
      ki_int_o = n_occ_int*(n_occ_act*(i-i_start)+k-1) + 1
      ji_int_o = n_occ_int*(n_occ_act*(i-i_start)+j-1) + 1
!
!
      bDci  => bDci_int(i_int_v:i_int_v+n_vir_aag-1)
!
      Lkai  => Lmci_int(ki_int_o:ki_int_o+n_occ_int-1)
      Ljai  => Lmci_int(ji_int_o:ji_int_o+n_occ_int-1)
!
      if(j_batch .eq. i_batch) then
         bDcj  => bDci_int(j_int_v:j_int_v+n_vir_aag-1)
!
         Lkbj  => Lmci_int(kj_int_o:kj_int_o+n_occ_int-1)
         Libj  => Lmci_int(ij_int_o:ij_int_o+n_occ_int-1)
      else
         bDcj  => bDcj_int(j_int_v:j_int_v+n_vir_aag-1)
!
         Lkbj  => Lmcj_int(kj_int_o:kj_int_o+n_occ_int-1)
         Libj  => Lmcj_int(ij_int_o:ij_int_o+n_occ_int-1)
      end if
!
      if(k_batch .eq. i_batch) then
         bDck  => bDci_int(k_int_v:k_int_v+n_vir_aag-1)
!
         Ljck  => Lmci_int(jk_int_o:jk_int_o+n_occ_int-1)
         Lick  => Lmci_int(ik_int_o:ik_int_o+n_occ_int-1)
      else if(k_batch .eq. j_batch) then
         bDck  => bDcj_int(k_int_v:k_int_v+n_vir_aag-1)
!
         Ljck  => Lmcj_int(jk_int_o:jk_int_o+n_occ_int-1)
         Lick  => Lmcj_int(ik_int_o:ik_int_o+n_occ_int-1)
      else
         bDck  => bDck_int(k_int_v:k_int_v+n_vir_aag-1)
!
         Ljck  => Lmck_int(jk_int_o:jk_int_o+n_occ_int-1)
         Lick  => Lmck_int(ik_int_o:ik_int_o+n_occ_int-1)
      end if
!
      beta = zero
!
      r2_aDji  => t2_aDji
      r2_abLi  => t2_abLi
!
!
!
      call cpu_time(time_start)
!
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
     &           bDck,n_vir_gen,beta,w_abc,n_vir_act)
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
     &           Ljck,n_occ_gen,one,w_abc,n_vir_act**2)
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
     &           bDck,n_vir_gen,zero,w_cba,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!     t^ba_jL*(Li|ck)
!     ---------------
!
      r2_off = (n_occ_gen*n_vir_act**2)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Lick,n_occ_gen,one,w_cba,n_vir_act**2)
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
      call add_bac(w_cba,w_abc,n_vir_act)
!
      call cpu_time(time_bac)
      time_bac = time_bac - time_1
!
!
!
!     t^aD_ik*(cD|bj)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDcj,n_vir_gen,zero,w_cba,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!     t^ac_iL*(Lk|bj)
!     ---------------
!
      r2_off = (n_occ_gen*n_vir_act**2)*(i-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Lkbj,n_occ_gen,one,w_cba,n_vir_act**2)
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
      call add_acb(w_cba,w_abc,n_vir_act)
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
     &           bDcj,n_vir_gen,zero,w_cba,n_vir_act)
!
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!     t^ca_kL*(Li|bj)
!     ---------------
!
      r2_off = (n_occ_gen*n_vir_act**2)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Libj,n_occ_gen,one,w_cba,n_vir_act**2)
!  
      call cpu_time(time_2)
      time_dgocc = time_dgocc + time_2 - time_1
!
!     Reorder and add to W
!     --------------------
!
      call cpu_time(time_1)
!
      call add_cab(w_cba,w_abc,n_vir_act)
!
      call cpu_time(time_cab)
      time_cab = time_cab - time_1
!
!      
!
!     t^bD_jk*(cD|ai)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDci,n_vir_gen,zero,w_cba,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!               
!     t^bc_jL*(Lk|ai)
!     ---------------
!
      r2_off = (n_occ_gen*n_vir_act**2)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Lkai,n_occ_gen,one,w_cba,n_vir_act**2)
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
      call add_bca(w_cba,w_abc,n_vir_act)
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
     &           bDci,n_vir_gen,zero,w_cba,n_vir_act)
!  
      call cpu_time(time_2)
      time_dgvir = time_dgvir + time_2 - time_1
!
!      
!               
!     t^cb_kL*(Lj|ai)
!     ---------------
!
      r2_off = (n_occ_gen*n_vir_act**2)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Ljai,n_occ_gen,one,w_cba,n_vir_act**2)
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
      call add_cba(w_cba,w_abc,n_vir_act)
!
      call cpu_time(time_cba)
      time_cba = time_cba - time_1
!
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
         write(lupri,*) 'Timings from mlccsdpt_w_calc'
         write(lupri,9999) 'Total', time_tot
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
   end subroutine mlccsdpt_w_calc
!      
!
   subroutine mlccsdpt_v_calc(i,j,k,i_batch,j_batch,k_batch)
!
!  V intermediate calculator
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate V^abc = W^abc + t_ai(jb|kc) + t_bj(ia|kc) + t_ck(ia|jb)
!
      implicit none
!      
      integer, intent(in)  :: i, j, k
      integer, intent(in)  :: i_batch, j_batch, k_batch
!      
      integer              :: i_off, j_off, k_off
      integer              :: i_start, j_start, k_start
!      
      integer              :: n_act_int
!
      integer              :: jk_int, ik_int, ij_int
!      
      integer              :: a,b,c
      integer              :: abc, cba, bac, acb
      real(dp)             :: abc_fac, cba_fac, acb_fac, bac_fac
!
      integer              :: ai, bj, ck
      integer              :: bi, cj, ak
      integer              :: ci, aj, bk
!
      integer              :: bc, ac, ab
      integer              :: cb, ca, ba
!      
      real(dp), dimension(:), pointer  :: jbkc   => null()
      real(dp), dimension(:), pointer  :: iakc   => null()
      real(dp), dimension(:), pointer  :: iajb   => null()
!      
      real     :: time_1, time_double, time_single
!
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
!      
!     Symmetry factors, i never equal k
!     ---------------------------------
!
      cba_fac = two
      if(i .eq. j) then
         abc_fac = three
         acb_fac = zero
         bac_fac = one
      else if(j .eq. k) then
         abc_fac = three
         acb_fac = one
         bac_fac = zero
      else
         abc_fac = six
         acb_fac = two
         bac_fac = two
      end if
!      
!      
      n_act_int = n_vir_act**2
!
!     Offsets
!     -------
!
      i_off = i + n_occ_inact
      j_off = j + n_occ_inact
      k_off = k + n_occ_inact
!
      i_start = batch_size*(i_batch-1) + 1
      j_start = batch_size*(j_batch-1) + 1
      k_start = batch_size*(k_batch-1) + 1
!
      ij_int = n_act_int*(n_occ_act*(j-j_start)+i-1) + 1
      ik_int = n_act_int*(n_occ_act*(k-k_start)+i-1) + 1
      jk_int = n_act_int*(n_occ_act*(k-k_start)+j-1) + 1
!      
!     Integrals pointers
!     ------------------
!
      iajb => iajb_int(ij_int:ij_int+n_act_int-1)
!
      if(k_batch .eq. j_batch) then
!
         iakc => iajb_int(ik_int:ik_int+n_act_int-1)
         jbkc => iajb_int(jk_int:jk_int+n_act_int-1)
!
      else
!
         iakc => iakc_int(ik_int:ik_int+n_act_int-1)
         jbkc => iakc_int(jk_int:jk_int+n_act_int-1)
!
      end if
!
!     Double contributions
!     --------------------
!
      call cpu_time(time_1)
!
!$omp parallel do schedule(static) private(a,b,c,abc,cba,acb,bac)
      do c = 1,n_vir_act
         do b = 1,n_vir_act
            do a = 1, n_vir_act
            
!              
               abc = n_vir_act**2*(c-1) + n_vir_act*(b-1) + a
               cba = n_vir_act**2*(a-1) + n_vir_act*(b-1) + c
!              
               v_abc(abc) = abc_fac*w_abc(abc) - cba_fac*w_abc(cba) 
!              
               if(i .ne. j) then
!              
                  acb = n_vir_act**2*(b-1) + n_vir_act*(c-1) + a
                  v_abc(abc) = v_abc(abc) - acb_fac*w_abc(acb) 
!
               end if
!              
               if(j .ne. k) then
!              
                  bac = n_vir_act**2*(c-1) + n_vir_act*(a-1) + b
                  v_abc(abc) = v_abc(abc) - bac_fac*w_abc(bac) 
!              
               end if
!              
            end do
         end do
      end do
!$omp end parallel do
!
      call cpu_time(time_double)
      time_double = time_double - time_1
!
!
!     Single contributions
!     --------------------
!
      call cpu_time(time_1)
!
!$omp parallel do schedule(static) private(a,b,c,abc,ai,bj,ck,bi,cj,ak,ci,aj,bk,ab,ac,bc,ba,cb,ca)
      do c = 1,n_vir_act
         do b = 1,n_vir_act
            do a = 1, n_vir_act
            
               ai = n_vir*(i_off-1) + a
               bj = n_vir*(j_off-1) + b
               ck = n_vir*(k_off-1) + c
               bi = n_vir*(i_off-1) + b
               cj = n_vir*(j_off-1) + c
               ak = n_vir*(k_off-1) + a
               ci = n_vir*(i_off-1) + c
               aj = n_vir*(j_off-1) + a
               bk = n_vir*(k_off-1) + b
!              
               ab = n_vir_act*(b-1) + a
               ac = n_vir_act*(c-1) + a
               bc = n_vir_act*(c-1) + b
               ba = n_vir_act*(a-1) + b
               ca = n_vir_act*(a-1) + c
               cb = n_vir_act*(b-1) + c
!              
               abc = n_vir_act**2*(c-1) + n_vir_act*(b-1) + a
!              
               v_abc(abc) = v_abc(abc) &
              &           + abc_fac*(t1am(ai)*jbkc(bc) + t1am(bj)*iakc(ac) + t1am(ck)*iajb(ab)) &
              &           - cba_fac*(t1am(ci)*jbkc(ba) + t1am(bj)*iakc(ca) + t1am(ak)*iajb(cb))
!              
               if(i .ne. j) then
!              
                  v_abc(abc) = v_abc(abc) &
                 &           - acb_fac*(t1am(ai)*jbkc(cb) + t1am(cj)*iakc(ab) + t1am(bk)*iajb(ac))
!
               end if
!              
               if(j .ne. k) then
!              
                  v_abc(abc) = v_abc(abc) &
                 &           - bac_fac*(t1am(bi)*jbkc(ac) + t1am(aj)*iakc(bc) + t1am(ck)*iajb(ba))
!              
               end if
!              
            end do
         end do
      end do
!$omp end parallel do
!
      call cpu_time(time_single)
      time_single = time_single - time_1
!
      if (print_mlcc3 .ge. 5) then
!
         write(lupri,*)
         write(lupri,*)
         write(lupri,*) 'Timings from mlccsdpt_v_calc'
         write(lupri,9999) 'doubles', time_double
         write(lupri,9999) 'singles', time_single
         write(lupri,*)
         write(lupri,*)
!
      end if
!
   end subroutine mlccsdpt_v_calc
!
!
   subroutine mlccsdpt_w_prod(energy,i,j,k)
!
!  W dot product
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: calculate E = (4W^abc + W^bca + W^cab)(V^abc - V^cba)/epsilon^abc
!
      implicit none
!      
      integer, intent(in)     :: i, j, k
      real(dp), intent(out)   :: energy
!      
      integer  :: i_off, j_off, k_off
!      
      integer  :: a,b,c
      integer  :: abc, bca, cab, cba
!      
      real(dp) :: eps_ijk, eps_abc
!      
      real(dp) :: ddot
      real(dp) :: w_mix_norm, v_mix_norm
!      
!     Timing variables
      real     :: time_start, time_tot
      real     :: time_w_mix, time_v_mix, time_dot
      real     :: time_1, time_2
!      
!
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
!
      call cpu_time(time_start)
!
!
!     Offsets
!     -------
!
      i_off = i + n_occ_inact
      j_off = j + n_occ_inact
      k_off = k + n_occ_inact
!
      eps_ijk = Fock_diagonal(i_off) &
     &        + Fock_diagonal(j_off) &
     &        + Fock_diagonal(k_off) 
!
!
      call cpu_time(time_1)
!
!     Set up (4W^abc + W^bca + W^cab)/epsilon^abc
!     -------------------------------------------
!
!$omp parallel do schedule(static) private(a,b,c,abc,bca,cab,eps_abc) shared(n_vir_act,eps_ijk)
      do c = 1,n_vir_act
         do b = 1,n_vir_act
            do a = 1, n_vir_act
!
               abc = n_vir_act**2*(c-1) + n_vir_act*(b-1) + a
!
               eps_abc = -one/(Fock_diagonal(n_occ + a) &
              &        +       Fock_diagonal(n_occ + b) &
              &        +       Fock_diagonal(n_occ + c) - eps_ijk)
!
               bca = n_vir_act**2*(a-1) + n_vir_act*(c-1) + b
               cab = n_vir_act**2*(b-1) + n_vir_act*(a-1) + c
!
               w_cba(abc) = (four*w_abc(abc) + w_abc(bca) + w_abc(cab))*eps_abc
!
            end do
         end do
      end do
!$omp end parallel do
!
      call cpu_time(time_w_mix)
      time_w_mix = time_w_mix - time_1
!
!
      if(print_mlcc3 .ge. 7) then
!
         w_mix_norm = ddot(n_vir_act**3,w_cba,1,w_cba,1)
!
         write(lupri,*)
         write(lupri,*) 'w_mix_norm: ', w_mix_norm
         write(lupri,*)
!
      end if
!
!
!     Do the energy dot product
!     -------------------------
!
      call cpu_time(time_1)
!
      energy = ddot(n_vir_act**3,w_cba,1,v_abc,1)
!
      call cpu_time(time_dot)
      time_dot = time_dot - time_1
!
!
!     Print timings
!     -------------
!
      if (print_mlcc3 .ge. 5) then
!
         write(lupri,*)
         write(lupri,*)
         write(lupri,*) 'Timings from mlccsdpt_w_prod'
         write(lupri,9999) 'w_mix', time_w_mix
         write(lupri,9999) 'ddot', time_dot
         write(lupri,*)
         write(lupri,*)
!
      end if
!
   end subroutine mlccsdpt_w_prod
!
!
end module mlccsdpt_energy
