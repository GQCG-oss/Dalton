module mlcc3_h_omega
!
!
!  Add MLCC3 contributions to omega vector
!  High memory integral batching
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
!  Integral arrays
   real(dp), dimension(:), pointer, private  :: bDci_g   => null()
   real(dp), dimension(:), pointer, private  :: bDcj_g   => null()
   real(dp), dimension(:), pointer, private  :: bDck_g   => null()
!
   real(dp), dimension(:), pointer, private  :: bDci_r   => null()
   real(dp), dimension(:), pointer, private  :: bDcj_r   => null()
   real(dp), dimension(:), pointer, private  :: bDck_r   => null()
!
   real(dp), dimension(:), pointer, private  :: Lmci_g   => null()
   real(dp), dimension(:), pointer, private  :: Lmcj_g   => null()
   real(dp), dimension(:), pointer, private  :: Lmck_g   => null()
!
   real(dp), dimension(:), pointer, private  :: Lmci_r   => null()
   real(dp), dimension(:), pointer, private  :: Lmcj_r   => null()
   real(dp), dimension(:), pointer, private  :: Lmck_r   => null()
!
   real(dp), dimension(:), pointer, private  :: mbic     => null()
   real(dp), dimension(:), pointer, private  :: mbjc     => null()
   real(dp), dimension(:), pointer, private  :: mbkc     => null()
!
   real(dp), dimension(:), pointer, private  :: Dbic     => null()
   real(dp), dimension(:), pointer, private  :: Dbjc     => null()
   real(dp), dimension(:), pointer, private  :: Dbkc     => null()
!
   real(dp), dimension(:), pointer, private  :: mLic     => null()
   real(dp), dimension(:), pointer, private  :: mLjc     => null()
   real(dp), dimension(:), pointer, private  :: mLkc     => null()
!
   real(dp), dimension(:), pointer, private  :: help_v2  => null()
!
!  Option variable for response calculations
   integer, private  :: response_integer
!      
!  Batch information
   integer, private  :: batch_size, n_batch
!      
!      
!  Print units
!  -----------
!
!  Ground state integrals
   integer, private  :: bDck_gu, Ljck_gu, jbkc_u, Dbkc_u, jLkc_u
!
!  Response integrals
   integer, private  :: bDck_ru, Ljck_ru
!
!      
contains
!
   subroutine mlcc3_h_omega_contrib(resp_control)
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
      integer  :: i_batch, j_batch, k_batch
      integer  :: i_start, j_start, k_start
      integer  :: i_end  , j_end  , k_end  
      integer  :: i_size , j_size , k_size 
!      
      integer  :: r_start,r_end, n_occ_int, n_act_int
!      
      integer  :: a, b, abij_n, abij_t, ioerror, i_cvs_rmc
      logical  :: skip_log
!      
      real(dp) :: ddot
      real(dp) :: w_norm, t3_norm
!
!     Timing variables
      real     :: time_tot, time_w, time_t3, time_o1, time_o2
      real     :: time_add, time_t2sq
      real     :: time_start, time_1, time_allo, time_deallo, time_zero
      real     :: time_w_sum, time_t3_sum, time_o1_sum, time_o2_sum
      real     :: time_open, time_close, time_batch
      real     :: time_read, time_L, time_2
!      
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
      time_w_sum     = 0.0
      time_t3_sum    = 0.0
      time_o1_sum    = 0.0
      time_o2_sum    = 0.0
      time_read      = 0.0
      time_allo      = 0.0
      time_L         = 0.0
!
      call cpu_time(time_start)
!
!     Default no CVS
      skip_log = .false.
!
      n_occ_int = n_vir_act*n_occ_gen*n_occ_act
      n_act_int = n_vir_act**2*n_occ_act
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
!     Help array to calculate L_bjck
      if(resp_control .ne. 2) then
         call work_allocator(help_v2,n_vir_act**2)
      end if
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
      call cpu_time(time_2)
      time_allo = time_allo + time_2 - time_1
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
      call mlcc3_h_t2_square
!
      if(resp_control .eq. 1) then
!
!        Square up C2 amplitudes as well
!
         response_integer = 1
!
         call mlcc3_h_t2_square
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
      call mlcc3_h_file_opener(resp_control)
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
      call mlcc3_h_batch_setup(resp_control)
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
      call work_allocator(bDci_g,batch_size*n_vir_aag)
      call work_allocator(Lmci_g,batch_size*n_occ_int)
      call work_allocator(Dbic,batch_size*n_vir_aag)
      call work_allocator(mLic,batch_size*n_occ_int)
      if(resp_control .eq. 1) then
         call work_allocator(bDci_r,batch_size*n_vir_aag)
         call work_allocator(Lmci_r,batch_size*n_occ_int)
      end if
      if(resp_control .ne. 2) then
         call work_allocator(mbic,batch_size*n_act_int)
      end if
!
!
      if(n_batch .ne. 1) then
         call work_allocator(bDcj_g,batch_size*n_vir_aag)
         call work_allocator(Lmcj_g,batch_size*n_occ_int)
         call work_allocator(Dbjc,batch_size*n_vir_aag)
         call work_allocator(mLjc,batch_size*n_occ_int)
         if(resp_control .eq. 1) then
            call work_allocator(bDcj_r,batch_size*n_vir_aag)
            call work_allocator(Lmcj_r,batch_size*n_occ_int)
         end if
         if(resp_control .ne. 2) then
            call work_allocator(mbjc,batch_size*n_act_int)
         end if
!
         call work_allocator(bDck_g,batch_size*n_vir_aag)
         call work_allocator(Lmck_g,batch_size*n_occ_int)
         call work_allocator(Dbkc,batch_size*n_vir_aag)
         call work_allocator(mLkc,batch_size*n_occ_int)
         if(resp_control .eq. 1) then
            call work_allocator(bDck_r,batch_size*n_vir_aag)
            call work_allocator(Lmck_r,batch_size*n_occ_int)
         end if
         if(resp_control .ne. 2) then
            call work_allocator(mbkc,batch_size*n_act_int)
         end if
      else
         bDcj_g => bDci_g
         Lmcj_g => Lmci_g
         Dbjc => Dbic
         mLjc => mLic
         if(resp_control .eq. 1) then
            bDcj_r => bDci_r
            bDck_r => bDci_r
         end if
         if(resp_control .ne. 2) then
            mbjc => mbic
         end if
!
         bDck_g => bDci_g
         Lmck_g => Lmci_g
         Dbkc => Dbic
         mLkc => mLic
         if(resp_control .eq. 1) then
            Lmcj_r => Lmci_r
            Lmck_r => Lmci_r
         end if
         if(resp_control .ne. 2) then
            mbkc => mbic
         end if
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
         call mlcc3_h_integral_reader(resp_control,bDci_g,bDci_r,Lmci_g,Lmci_r,&
        &                             mbic,Dbic,mLic,i_start,i_size)
!
         call cpu_time(time_2)
         time_read = time_read + time_2 - time_1
!
         call cpu_time(time_1)
!
         if(resp_control .ne. 2) then
            call mlcc3_h_L_calc(mbic,i_size,i_start)
         end if
!
         call cpu_time(time_2)
         time_L = time_L + time_2 - time_1
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
            if(j_batch .ne. i_batch) then
               call cpu_time(time_1)
!
               call mlcc3_h_integral_reader(resp_control,bDcj_g,bDcj_r,Lmcj_g,Lmcj_r,&
              &                             mbjc,Dbjc,mLjc,j_start,j_size)
!
               call cpu_time(time_2)
               time_read = time_read + time_2 - time_1
!
               call cpu_time(time_1)
!
               if(resp_control .ne. 2) then
                  call mlcc3_h_L_calc(mbjc,j_size,j_start)
               end if
!
               call cpu_time(time_2)
               time_L = time_L + time_2 - time_1
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
               if(k_batch .ne. i_batch .and. k_batch .ne. j_batch) then
                  call cpu_time(time_1)
!
                  call mlcc3_h_integral_reader(resp_control,bDck_g,bDck_r,Lmck_g,Lmck_r,&
                 &                             mbkc,Dbkc,mLkc,k_start,k_size)
!
                  call cpu_time(time_2)
                  time_read = time_read + time_2 - time_1
!
                  call cpu_time(time_1)
!
                  if(resp_control .ne. 2) then
                     call mlcc3_h_L_calc(mbkc,k_size,k_start)
                  end if
!
                  call cpu_time(time_2)
                  time_L = time_L + time_2 - time_1
!
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
!                       Core-valence separation or remove core
                        if(resp_control .eq. 1) then
                           if(cvs_log) then
!
                              skip_log = .true.
                              do i_cvs_rmc = 1,n_cvs_rmc
                                 if(i .eq. cvs_rmc_index(i_cvs_rmc)) then
                                    skip_log = .false.
                                 else if(j .eq. cvs_rmc_index(i_cvs_rmc)) then
                                    skip_log = .false.
                                 else if(k .eq. cvs_rmc_index(i_cvs_rmc)) then
                                    skip_log = .false.
                                 end if
                              end do
!
                           else if(rmc_log) then
!
                              skip_log = .false.
                              do i_cvs_rmc = 1,n_cvs_rmc
                                 if(i .eq. cvs_rmc_index(i_cvs_rmc)) then
                                    skip_log = .true.
                                 else if(j .eq. cvs_rmc_index(i_cvs_rmc)) then
                                    skip_log = .true.
                                 else if(k .eq. cvs_rmc_index(i_cvs_rmc)) then
                                    skip_log = .true.
                                 end if
                              end do
                           end if
!
!                          Cycle if removed
                           if(skip_log) then
                              cycle
                           end if
                        end if
!                  
!                  
                        call cpu_time(time_1)
!               
                        if(resp_control .eq. 0 .or. resp_control .eq. 2) then
!
!                          Calculate standard T_3 amplitudes
!
                           response_integer = 0
!
!                          W intermediates               
                           call mlcc3_h_w_calc(i,j,k,i_batch,j_batch,k_batch)
!
                        else if(resp_control .eq. 1) then
!
!                          Calculate C_3 resonse amplitudes
!
!                          [Ĥ,C_2] term in W intermdeiates
                           response_integer = 1
!
                           call mlcc3_h_w_calc(i,j,k,i_batch,j_batch,k_batch)
!               
!                          [[Ĥ,C_1],T_2] term in W intermediates
                           response_integer = 2
!
                           call mlcc3_h_w_calc(i,j,k,i_batch,j_batch,k_batch)
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
                        if(print_mlcc3 .ge. 7) then
!
                           w_norm = ddot(n_vir_act**3,t_abc,1,t_abc,1)
!
                           write(lupri,*)
                           write(lupri,*) 'w_norm: ', w_norm
                           write(lupri,*)
!
                        end if
!
!               
!                       Divide by epsilon - frequency. Frequency is zero
!                       for an energy calculation
                        call cpu_time(time_1)
!               
                        if(resp_control .eq. 0 .or. resp_control .eq. 2) then
!
                           response_integer = 0
!
                        else if(resp_control .eq. 1) then
!
                           response_integer = 1
!
                        end if
!               
                        call mlcc3_h_t3_calculator(i,j,k)
!               
                        call cpu_time(time_t3)
                        time_t3 = time_t3 - time_1
                        time_t3_sum = time_t3_sum + time_t3
!               
                        if(print_mlcc3 .ge. 7) then
!
                           t3_norm = ddot(n_vir_act**3,t_abc,1,t_abc,1)
!
                           write(lupri,*)
                           write(lupri,*) 't3_norm: ', t3_norm
                           write(lupri,*)
!
                        end if
!
!
                        call cpu_time(time_1)
!               
!                       Omega 1 and Fock contributions.
                        if(resp_control .eq. 0 .or. resp_control .eq. 1) then
!
!                          Standard [Ĥ,T_3] or [Ĥ,C_3], same code
                           response_integer = 0
!
                        else if(resp_control .eq. 2) then
!
!                          Response [[Ĥ,C_1],T_3] Fock contributions
                           response_integer = 1
!
                        end if
!               
                        call mlcc3_h_omega_1(i,j,k,i_batch,j_batch,k_batch)
!
                        call cpu_time(time_o1)
                        time_o1 = time_o1 - time_1
                        time_o1_sum = time_o1_sum + time_o1
!               
!               
                        call cpu_time(time_1)
!               
!                       Omega 2 contributions
!               
                        call mlcc3_h_omega_2(i,j,k,i_batch,j_batch,k_batch)
!               
                        call cpu_time(time_o2)
                        time_o2 = time_o2 - time_1
                        time_o2_sum = time_o2_sum + time_o2
!               
!               
!                       Print timings
!                       -------------
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
               end do !End restricted loop
!      
            end do
         end do
      end do !End batch loop
!      
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
!
!     Add temporary Omega 2 to real Omega 2
!     -------------------------------------
      call cpu_time(time_1)
!
      call mlcc3_h_omega_add
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
      call mlcc3_h_file_closer(resp_control)
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
         if(resp_control .ne. 2) then
            call work_deallocator(mbkc)
         end if
         if(resp_control .eq. 1) then
            call work_deallocator(Lmck_r)
            call work_deallocator(bDck_r)
         end if
         call work_deallocator(mLkc)
         call work_deallocator(Dbkc)
         call work_deallocator(Lmck_g)
         call work_deallocator(bDck_g)
!
         if(resp_control .ne. 2) then
            call work_deallocator(mbjc)
         end if
         if(resp_control .eq. 1) then
            call work_deallocator(Lmcj_r)
            call work_deallocator(bDcj_r)
         end if
         call work_deallocator(mLjc)
         call work_deallocator(Dbjc)
         call work_deallocator(Lmcj_g)
         call work_deallocator(bDcj_g)
      end if
!
      if(resp_control .ne. 2) then
         call work_deallocator(mbic)
      end if
      if(resp_control .eq. 1) then
         call work_deallocator(Lmci_r)
         call work_deallocator(bDci_r)
      end if
      call work_deallocator(mLic)
      call work_deallocator(Dbic)
      call work_deallocator(Lmci_g)
      call work_deallocator(bDci_g)
!
!     Deallocate arrays
!     -----------------
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
      if(resp_control .ne. 2) then
         call work_deallocator(help_v2)
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
         write(lupri,9999) 'set up batching', time_batch
         write(lupri,9999) 'read integrals', time_read
         write(lupri,9999) 'calculate L', time_L
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
   end subroutine mlcc3_h_omega_contrib
!      
!      
   subroutine mlcc3_h_w_calc(i,j,k,i_batch,j_batch,k_batch)
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
      if(response_integer .eq. 0) then
!
!        Calculate T_3 amplitudes
!
         bDci  => bDci_g(i_int_v:i_int_v+n_vir_aag-1)
!
         Lkai  => Lmci_g(ki_int_o:ki_int_o+n_occ_int-1)
         Ljai  => Lmci_g(ji_int_o:ji_int_o+n_occ_int-1)
!
         if(j_batch .eq. i_batch) then
            bDcj  => bDci_g(j_int_v:j_int_v+n_vir_aag-1)
!
            Lkbj  => Lmci_g(kj_int_o:kj_int_o+n_occ_int-1)
            Libj  => Lmci_g(ij_int_o:ij_int_o+n_occ_int-1)
         else
            bDcj  => bDcj_g(j_int_v:j_int_v+n_vir_aag-1)
!
            Lkbj  => Lmcj_g(kj_int_o:kj_int_o+n_occ_int-1)
            Libj  => Lmcj_g(ij_int_o:ij_int_o+n_occ_int-1)
         end if
!
         if(k_batch .eq. i_batch) then
            bDck  => bDci_g(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmci_g(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmci_g(ik_int_o:ik_int_o+n_occ_int-1)
         else if(k_batch .eq. j_batch) then
            bDck  => bDcj_g(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmcj_g(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmcj_g(ik_int_o:ik_int_o+n_occ_int-1)
         else
            bDck  => bDck_g(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmck_g(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmck_g(ik_int_o:ik_int_o+n_occ_int-1)
         end if
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
         bDci  => bDci_g(i_int_v:i_int_v+n_vir_aag-1)
!
         Lkai  => Lmci_g(ki_int_o:ki_int_o+n_occ_int-1)
         Ljai  => Lmci_g(ji_int_o:ji_int_o+n_occ_int-1)
!
         if(j_batch .eq. i_batch) then
            bDcj  => bDci_g(j_int_v:j_int_v+n_vir_aag-1)
!
            Lkbj  => Lmci_g(kj_int_o:kj_int_o+n_occ_int-1)
            Libj  => Lmci_g(ij_int_o:ij_int_o+n_occ_int-1)
         else
            bDcj  => bDcj_g(j_int_v:j_int_v+n_vir_aag-1)
!
            Lkbj  => Lmcj_g(kj_int_o:kj_int_o+n_occ_int-1)
            Libj  => Lmcj_g(ij_int_o:ij_int_o+n_occ_int-1)
         end if
!
         if(k_batch .eq. i_batch) then
            bDck  => bDci_g(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmci_g(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmci_g(ik_int_o:ik_int_o+n_occ_int-1)
         else if(k_batch .eq. j_batch) then
            bDck  => bDcj_g(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmcj_g(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmcj_g(ik_int_o:ik_int_o+n_occ_int-1)
         else
            bDck  => bDck_g(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmck_g(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmck_g(ik_int_o:ik_int_o+n_occ_int-1)
         end if
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
         bDci  => bDci_r(i_int_v:i_int_v+n_vir_aag-1)
!
         Lkai  => Lmci_r(ki_int_o:ki_int_o+n_occ_int-1)
         Ljai  => Lmci_r(ji_int_o:ji_int_o+n_occ_int-1)
!
         if(j_batch .eq. i_batch) then
            bDcj  => bDci_r(j_int_v:j_int_v+n_vir_aag-1)
!
            Lkbj  => Lmci_r(kj_int_o:kj_int_o+n_occ_int-1)
            Libj  => Lmci_r(ij_int_o:ij_int_o+n_occ_int-1)
         else
            bDcj  => bDcj_r(j_int_v:j_int_v+n_vir_aag-1)
!
            Lkbj  => Lmcj_r(kj_int_o:kj_int_o+n_occ_int-1)
            Libj  => Lmcj_r(ij_int_o:ij_int_o+n_occ_int-1)
         end if
!
         if(k_batch .eq. i_batch) then
            bDck  => bDci_r(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmci_r(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmci_r(ik_int_o:ik_int_o+n_occ_int-1)
         else if(k_batch .eq. j_batch) then
            bDck  => bDcj_r(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmcj_r(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmcj_r(ik_int_o:ik_int_o+n_occ_int-1)
         else
            bDck  => bDck_r(k_int_v:k_int_v+n_vir_aag-1)
!
            Ljck  => Lmck_r(jk_int_o:jk_int_o+n_occ_int-1)
            Lick  => Lmck_r(ik_int_o:ik_int_o+n_occ_int-1)
         end if
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
!
      r2_off = (n_occ_gen*n_vir_act**2)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act**2,n_vir_act,n_occ_gen,-one,r2_abLi(r2_off),n_vir_act**2, &
     &           Lick,n_occ_gen,one,u_abc,n_vir_act**2)
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
!
!     t^aD_ik*(cD|bj)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDcj,n_vir_gen,zero,u_abc,n_vir_act)
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
     &           Lkbj,n_occ_gen,one,u_abc,n_vir_act**2)
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
     &           bDcj,n_vir_gen,zero,u_abc,n_vir_act)
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
     &           Libj,n_occ_gen,one,u_abc,n_vir_act**2)
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
!     t^bD_jk*(cD|ai)
!     ---------------
!
      r2_off = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(k-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_act**2,n_vir_gen,one,r2_aDji(r2_off),n_vir_act, &
     &           bDci,n_vir_gen,zero,u_abc,n_vir_act)
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
     &           Lkai,n_occ_gen,one,u_abc,n_vir_act**2)
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
     &           bDci,n_vir_gen,zero,u_abc,n_vir_act)
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
     &           Ljai,n_occ_gen,one,u_abc,n_vir_act**2)
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
   end subroutine mlcc3_h_w_calc
!      
!
   subroutine mlcc3_h_t3_calculator(i,j,k)
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
      integer              :: i_off, j_off, k_off, ion_int
      integer              :: a, b, c
      real(dp)             :: epsilon_ijk, epsilon_abc
      logical              :: ion_zero
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
!              Set to zero if ionising  
               if(ion_log .and. response_integer .eq. 1) then
                  ion_zero = .true.
                  do ion_int=1,n_ion
                     if(a .eq. ion_index(ion_int) .or. &
                    &   b .eq. ion_index(ion_int) .or. &
                    &   c .eq. ion_index(ion_int)) then
!  
                        ion_zero = .false.
!  
                     end if
                     if(ion_zero) then
                        t_abc(n_vir_act**2*(c-1) + n_vir_act*(b-1) + a) = 0
                     end if
                  end do
               end if
!  
!  
            end do
         end do
      end do
!$omp end parallel do
!      
   end subroutine mlcc3_h_t3_calculator
!
!
   subroutine mlcc3_h_omega_1(i,j,k,i_batch,j_batch,k_batch)
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
      integer,intent(in)   :: i_batch, j_batch, k_batch
!      
      integer              :: i_start, j_start, k_start
!
      integer              :: ji_int, ki_int, ij_int
      integer              :: kj_int, ik_int, jk_int
!      
      integer              :: i_off, j_off, k_off
      integer              :: i_off2, j_off2, k_off2
      integer              :: a, b, ab, abij
!      
      integer              :: omega_off, F_off, n_act_int
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
      real(dp), dimension(:), pointer  :: fock_mat => null()
!      
!     Timing variables
      real     :: time_tot, time_dgmv1, time_dgmvf, time_ome, time_order
      real     :: time_start, time_trans, time_L
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
!      
!      
      call cpu_time(time_1)
!
      if(response_integer .eq. 0) then 
!      
!      
!        L setup
!        -------
!
         n_act_int = n_vir_act**2
!
         i_start = batch_size*(i_batch-1) + 1
         j_start = batch_size*(j_batch-1) + 1
         k_start = batch_size*(k_batch-1) + 1
!
         ji_int = n_act_int*(n_occ_act*(i-i_start)+j-1) + 1
         ki_int = n_act_int*(n_occ_act*(i-i_start)+k-1) + 1
         ij_int = n_act_int*(n_occ_act*(j-j_start)+i-1) + 1
         kj_int = n_act_int*(n_occ_act*(j-j_start)+k-1) + 1
         ik_int = n_act_int*(n_occ_act*(k-k_start)+i-1) + 1
         jk_int = n_act_int*(n_occ_act*(k-k_start)+j-1) + 1
!
!      
         L_jbic => mbic(ji_int:ji_int+n_act_int-1)
         L_kbic => mbic(ki_int:ki_int+n_act_int-1)
!
         if(j_batch .eq. i_batch) then
            L_ibjc => mbic(ij_int:ij_int+n_act_int-1)
            L_kbjc => mbic(kj_int:kj_int+n_act_int-1)
         else
            L_ibjc => mbjc(ij_int:ij_int+n_act_int-1)
            L_kbjc => mbjc(kj_int:kj_int+n_act_int-1)
         end if
!
         if(k_batch .eq. i_batch) then
            L_ibkc => mbic(ik_int:ik_int+n_act_int-1)
            L_jbkc => mbic(jk_int:jk_int+n_act_int-1)
         else if(k_batch .eq. j_batch) then
            L_ibkc => mbjc(ik_int:ik_int+n_act_int-1)
            L_jbkc => mbjc(jk_int:jk_int+n_act_int-1)
         else
            L_ibkc => mbkc(ik_int:ik_int+n_act_int-1)
            L_jbkc => mbkc(jk_int:jk_int+n_act_int-1)
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
!        Add to Omega 1
!        --------------
!
         omega_off   = n_vir*(i_off - 1) + 1
!
         call cpu_time(time_1)
!
         call dgemv('N',n_vir_act,n_vir_act**2,one,u_abc,n_vir_act,L_jbkc,1,&
        &           one,omega1(omega_off),1)
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
!           Add to Omega 1
!           --------------
!      
            omega_off   = n_vir*(k_off - 1) + 1
!      
            call cpu_time(time_1)
!
            call dgemv('N',n_vir_act,n_vir_act**2,-one,u_abc,n_vir_act,L_jbic,1,&
           &           one,omega1(omega_off),1)
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
        &           one,omega_o(omega_off),1)
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
            call dgemv('N',n_vir_act,n_vir_act**2,one,u_abc,n_vir_act,L_kbjc,1,&
           &           one,omega1(omega_off),1)
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
               call dgemv('N',n_vir_act,n_vir_act**2,-one,u_abc,n_vir_act,L_kbic,1, &
              &           one,omega1(omega_off),1)
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
            call dgemv('N',n_vir_act,n_vir_act**2,one,u_abc,n_vir_act,L_ibkc,1,&
           &           one,omega1(omega_off),1)
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
!           |Calculate (t^bca - t^bac)*L_ibjc and 2(t^bca - t^bac)F_jc|
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
               call dgemv('N',n_vir_act,n_vir_act**2,-one,u_abc,n_vir_act,L_ibjc,1,&
              &           one,omega1(omega_off),1)
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
   end subroutine mlcc3_h_omega_1
!
!
!
   subroutine mlcc3_h_omega_2(i,j,k,i_batch,j_batch,k_batch)
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
      integer,intent(in)   :: i_batch, j_batch, k_batch
!
      integer              :: i_start, j_start, k_start
!
      integer              :: i_off, j_off, k_off
!      
      integer              :: i_int_v, j_int_v, k_int_v
!      
      integer              :: ik_int_o, jk_int_o, ij_int_o
      integer              :: kj_int_o, ki_int_o, ji_int_o
!      
      integer              :: omega_off
!
      integer              :: n_occ_int
!
!     (Db|kc) and (jL|kc) integrals
      real(dp), dimension(:), pointer  :: Dbicx    => null()
      real(dp), dimension(:), pointer  :: Dbjcx    => null()
      real(dp), dimension(:), pointer  :: Dbkcx    => null()
!
      real(dp), dimension(:), pointer  :: kLicx    => null()
      real(dp), dimension(:), pointer  :: jLicx    => null()
      real(dp), dimension(:), pointer  :: iLjcx    => null()
      real(dp), dimension(:), pointer  :: kLjcx    => null()
      real(dp), dimension(:), pointer  :: iLkcx    => null()
      real(dp), dimension(:), pointer  :: jLkcx    => null()
!
      real(dp)             :: alpha
!      
!      
!     Timing variables
      real     :: time_tot, time_occ, time_vir, time_order
      real     :: time_start
      real     :: time_abc, time_bac, time_acb, time_bca, time_cba, time_cab
      real     :: time_1, time_2
!      
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
!
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
!      
      call cpu_time(time_start)
!
!     Set up integrals
!     ----------------
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
      ki_int_o = n_occ_int*(n_occ_act*(i-i_start)+k-1) + 1
      ji_int_o = n_occ_int*(n_occ_act*(i-i_start)+j-1) + 1
      kj_int_o = n_occ_int*(n_occ_act*(j-j_start)+k-1) + 1
      ij_int_o = n_occ_int*(n_occ_act*(j-j_start)+i-1) + 1
      jk_int_o = n_occ_int*(n_occ_act*(k-k_start)+j-1) + 1
      ik_int_o = n_occ_int*(n_occ_act*(k-k_start)+i-1) + 1
!
!     i integrals
      Dbicx => Dbic(i_int_v:i_int_v+n_vir_aag-1)
!
      kLicx => mLic(ki_int_o:ki_int_o+n_occ_int-1)
      jLicx => mLic(ji_int_o:ji_int_o+n_occ_int-1)
!
!     j integrals
      if(j_batch .eq. i_batch) then
         Dbjcx => Dbic(j_int_v:j_int_v+n_vir_aag-1)
!
         iLjcx => mLic(ij_int_o:ij_int_o+n_occ_int-1)
         kLjcx => mLic(kj_int_o:kj_int_o+n_occ_int-1)
      else
         Dbjcx => Dbjc(j_int_v:j_int_v+n_vir_aag-1)
!
         iLjcx => mLjc(ij_int_o:ij_int_o+n_occ_int-1)
         kLjcx => mLjc(kj_int_o:kj_int_o+n_occ_int-1)
      end if
!
!     k integrals
      if(k_batch .eq. i_batch) then
         Dbkcx => Dbic(k_int_v:k_int_v+n_vir_aag-1)
!
         iLkcx => mLic(ik_int_o:ik_int_o+n_occ_int-1)
         jLkcx => mLic(jk_int_o:jk_int_o+n_occ_int-1)
      else if(k_batch .eq. j_batch) then
         Dbkcx => Dbjc(k_int_v:k_int_v+n_vir_aag-1)
!
         iLkcx => mLjc(ik_int_o:ik_int_o+n_occ_int-1)
         jLkcx => mLjc(jk_int_o:jk_int_o+n_occ_int-1)
      else
         Dbkcx => Dbkc(k_int_v:k_int_v+n_vir_aag-1)
!
         iLkcx => mLkc(ik_int_o:ik_int_o+n_occ_int-1)
         jLkcx => mLkc(jk_int_o:jk_int_o+n_occ_int-1)
      end if
!
!     Calculate occupied active indices
!     ---------------------------------
      i_off = i + n_occ_inact
      j_off = j + n_occ_inact
      k_off = k + n_occ_inact
!
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
      omega_off   = (n_occ_act*n_vir_gen*n_vir_act)*(j-1) + (n_vir_gen*n_vir_act)*(i-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
     &           Dbkcx,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
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
     &              iLkcx,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
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
      omega_off   = (n_occ_act*n_vir_gen*n_vir_act)*(i-1) + (n_vir_gen*n_vir_act)*(j-1) + 1
!
      call cpu_time(time_1)
!
      call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
     &           Dbkcx,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
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
     &           jLkcx,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
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
         if (i .eq. k) then
            alpha = half
         else
            alpha = one
         end if
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
        &           Dbjcx,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
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
           &           iLjcx,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
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
        &           Dbjcx,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
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
     &              kLjcx,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
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
            if (k .eq. j) then
               alpha = half
            else
               alpha = one
            end if
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
        &              Dbicx,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
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
           &              kLicx,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
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
!           -----------------------------------------------
!           |Calculate u^cab*(jL|ic) and add to Omega_abLk|
!           -----------------------------------------------
!      
            omega_off = (n_occ_gen*n_vir_act**2)*(k-1) + 1
!
            call cpu_time(time_1)
!
            call dgemm('N','N',n_vir_act**2,n_occ_gen,n_vir_act,-one,u_abc,n_vir_act**2, &
        &              jLicx,n_vir_act,one,omega_o(omega_off),n_vir_act**2)
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
!           -----------------------------------------------
!           |Calculate u^cab*(Db|ic) and add to Omega_aDjk|
!           -----------------------------------------------
!      
            omega_off = (n_occ_act*n_vir_gen*n_vir_act)*(k-1) + (n_vir_gen*n_vir_act)*(j-1) + 1
!
            call cpu_time(time_1)
!
            call dgemm('N','N',n_vir_act,n_vir_gen,n_vir_act**2,alpha,u_abc,n_vir_act, &
        &              Dbicx,n_vir_act**2,one,omega_v(omega_off),n_vir_act)
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
   end subroutine mlcc3_h_omega_2
!
!
   subroutine mlcc3_h_omega_add
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
   end subroutine mlcc3_h_omega_add
!
!
   subroutine mlcc3_h_t2_square
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
   end subroutine mlcc3_h_t2_square
!
!
   subroutine mlcc3_h_file_opener(io_opt)
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
!     Ground state integrals used for amplitudes. Always needed
!     ---------------------------------------------------------
!
      bDck_gu = new_unit()
!   
      open(unit=bDck_gu, file=bDck_file_name, action='read', status='old', &
     &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening bDck_g file"
         write(lupri,*) "Error code:   ", ioerror
         write(lupri,*) "io_opt:       ", io_opt
         call quit('Error in file opener')
      end if
!   
      Ljck_gu = new_unit()
!      
      open(unit=Ljck_gu, file=Ljck_file_name, action='read', status='old', &
     &     access='direct', form='unformatted', recl=dp*n_occ_read, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening Ljck_g file"
         write(lupri,*) "Error code:   ", ioerror
         write(lupri,*) "io_opt:       ", io_opt
         call quit('Error in file opener')
      end if
!      
!
!     Excited state integrals used for amplitudes. Needed for case 1
!     --------------------------------------------------------------
!
      if(io_opt .eq. 1) then
!      
         bDck_ru = new_unit()
!      
         open(unit=bDck_ru, file=bDck_resp_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening bDck_r file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
         Ljck_ru = new_unit()
!         
         open(unit=Ljck_ru, file=Ljck_resp_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_occ_read, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Ljck_r file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!      
      end if
!      
!      
!     jbkc integrals needed for Rho 1 calculations in case 0 and 1
!     ------------------------------------------------------------
!
      if(io_opt .eq. 0 .or. io_opt .eq. 1) then
!
         jbkc_u = new_unit()
!         
         open(unit=jbkc_u, file=jbkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_act**2, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
      end if
!      
!      
!     Omega 2 ground state integrals needed for case 0 and 1
!     ------------------------------------------------------
!
      if(io_opt .eq. 0 .or. io_opt .eq. 1) then
!
         Dbkc_u = new_unit()
!      
         open(unit=Dbkc_u, file=Dbkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Dbkc_g file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
         jLkc_u = new_unit()
!      
         open(unit=jLkc_u, file=jLkc_file_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_act*n_occ_gen, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jLkc_g file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
      else if(io_opt .eq. 2) then 
!
!     Omega 2 excited state integrals needed for case 2
!     -------------------------------------------------
!
         Dbkc_u = new_unit()
!      
         open(unit=Dbkc_u, file=Dbkc_resp_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening Dbkc_r file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
         jLkc_u = new_unit()
!      
         open(unit=jLkc_u, file=jLkc_resp_name, action='read', status='old', &
        &     access='direct', form='unformatted', recl=dp*n_vir_act*n_occ_gen, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jLkc_r file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file opener')
         end if
!
      end if
!      
!
   end subroutine mlcc3_h_file_opener
!
   subroutine mlcc3_h_file_closer(io_opt)
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
!     Ground state integrals used for amplitudes. Always needed
!     ---------------------------------------------------------
!
      close(bDck_gu,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing bDck_g file"
         write(lupri,*) "Error code:   ", ioerror
         write(lupri,*) "io_opt:       ", io_opt
         call quit('Error in file closer')
      end if
!   
      close(Ljck_gu,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing Ljck_g file"
         write(lupri,*) "Error code:   ", ioerror
         write(lupri,*) "io_opt:       ", io_opt
         call quit('Error in file closer')
      end if
!      
!
!     Excited state integrals used for amplitudes. Needed for case 1
!     --------------------------------------------------------------
!
      if(io_opt .eq. 1) then
!      
         close(bDck_ru,status='keep',iostat=ioerror)
!  
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing bDck_r file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!      
!      
         close(Ljck_ru,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing Ljck_r file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!
      end if
!      
!      
!     jbkc integrals needed for Rho 1 calculations in case 0 and 1
!     ------------------------------------------------------------
!
      if(io_opt .eq. 0 .or. io_opt .eq. 1) then
!
!        Integrals for Omega1 calculations
!        ---------------------------------
!
         close(jbkc_u,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing jbkc file"
            write(lupri,*) "Error code:   ", ioerror
            write(lupri,*) "io_opt:       ", io_opt
            call quit('Error in file closer')
         end if
!            
      end if
!      
!     Integrals used in Omega2 calculation
!     ------------------------------------
!
      close(Dbkc_u,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing Dbkc file"
         write(lupri,*) "Error code:   ", ioerror
         write(lupri,*) "io_opt:       ", io_opt
         call quit('Error in file closer')
      end if
!         
      close(jLkc_u,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing jLkc file"
         write(lupri,*) "Error code:   ", ioerror
         write(lupri,*) "io_opt:       ", io_opt
         call quit('Error in file closer')
      end if
!
!
   end subroutine mlcc3_h_file_closer
!
!
   subroutine mlcc3_h_batch_setup(resp_control)
!
!  Batch setup
!  Authors: Henrik Koch and Rolf H. Myhre
!  May 2015
!
!  Purpose: Setup batching over integrals in the occupied loop
!
      implicit none
!
      integer, intent(in)  :: resp_control
!
      integer              :: w_req, w_avail, n_occ_int, n_act_int
!
!
      w_avail = work_free()
!
      n_occ_int = n_vir_act*n_occ_gen*n_occ_act
      n_act_int = n_vir_act**2*n_occ_act
!
      if(resp_control .eq. 0) then
!
         w_req = 2*(n_vir_aag+n_occ_int) + n_act_int
!
      else if(resp_control .eq. 1) then
!
         w_req = 3*(n_vir_aag+n_occ_int) + n_act_int
!
      else if(resp_control .eq. 2) then
!
         w_req = 2*(n_vir_aag+n_occ_int)
!
      else
!
         call quit('resp_control not 0,1 or 2')
!
      end if
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
         batch_size  = w_avail/(3*w_req)
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
   end subroutine mlcc3_h_batch_setup
!
!
   subroutine mlcc3_h_integral_reader(resp_control,bDcx_g,bDcx_r,Lmcx_g,Lmcx_r,&
  &                                   mbxc,Dbxc,mLxc,x_start,x_size)
!
!  Integral reader
!  Authors: Henrik Koch and Rolf H. Myhre
!  May 2015
!
!  Purpose: Read integrals in batch
!
      implicit none
!
      integer, intent(in)  :: resp_control,x_start,x_size
!
      real(dp), dimension(:), pointer, intent(inout)  :: bDcx_g
      real(dp), dimension(:), pointer, intent(inout)  :: bDcx_r
      real(dp), dimension(:), pointer, intent(inout)  :: Lmcx_g
      real(dp), dimension(:), pointer, intent(inout)  :: Lmcx_r
      real(dp), dimension(:), pointer, intent(inout)  :: mbxc  
      real(dp), dimension(:), pointer, intent(inout)  :: Dbxc  
      real(dp), dimension(:), pointer, intent(inout)  :: mLxc  
!
      integer  :: r_start,r_end,ioerror,x,y,n_occ_int,n_act_int,record
!
      n_occ_int = n_vir_act*n_occ_gen
      n_act_int = n_vir_act**2
!
!     bDcx ground state integrals
      do x = 1,x_size
!
         r_start  = n_vir_aag*(x-1)+1
         r_end    = r_start + n_vir_aag -1
         read(bDck_gu,rec=x+x_start-1,iostat=ioerror) bDcx_g(r_start:r_end)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error reading bDcx_g file"
            write(lupri,*) "Error code: ", ioerror
            call quit('Error reading integral file')
         end if
!
      end do
!
!     Lmcx ground state integrals
      do x = 1,x_size
         do y = 1,n_occ_act
!      
            r_start  = n_occ_int*(n_occ_act*(x-1)+y-1)+1
            r_end    = r_start + n_occ_int -1
            record   = n_occ_act*(x+x_start-2) + y
            read(Ljck_gu,rec=record,iostat=ioerror) Lmcx_g(r_start:r_end)
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
!     mbxc integrals
      if(resp_control .ne. 2) then
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
      end if
!         
!     Dbxc integrals
      do x = 1,x_size
!
         r_start  = n_vir_aag*(x-1)+1
         r_end    = r_start + n_vir_aag -1
         read(Dbkc_u,rec=x+x_start-1,iostat=ioerror) Dbxc(r_start:r_end)
!         
         if(ioerror .ne. 0) then
            write(lupri,*) "Error reading Dbxc file"
            write(lupri,*) "Error code: ", ioerror
            call quit('Error reading integral file')
         end if
!
      end do
!         
!     mLxc integrals
      do x = 1,x_size
         do y = 1,n_occ_act
!      
            r_start  = n_occ_int*(n_occ_act*(x-1)+y-1)+1
            r_end    = r_start + n_occ_int -1
            record   = n_occ_act*(x+x_start-2) + y
            read(jLkc_u,rec=record,iostat=ioerror) mLxc(r_start:r_end)
!      
            if(ioerror .ne. 0) then
               write(lupri,*) "Error reading mLxc file"
               write(lupri,*) "Error code: ", ioerror
               call quit('Error reading integral file')
            end if
!      
         end do
      end do
!         
!         
      if(resp_control .eq. 1) then
!         
!        bDcx response integrals
         do x = 1,x_size
!         
            r_start  = n_vir_aag*(x-1)+1
            r_end    = r_start + n_vir_aag -1
            read(bDck_ru,rec=x+x_start-1,iostat=ioerror) bDcx_r(r_start:r_end)
!         
            if(ioerror .ne. 0) then
               write(lupri,*) "Error reading bDcx_r file"
               write(lupri,*) "Error code: ", ioerror
               call quit('Error reading integral file')
            end if
!         
         end do
!         
!        Lmcx response integrals
         do x = 1,x_size
            do y = 1,n_occ_act
!         
               r_start  = n_occ_int*(n_occ_act*(x-1)+y-1)+1
               r_end    = r_start + n_occ_int -1
               record   = n_occ_act*(x+x_start-2) + y
               read(Ljck_ru,rec=record,iostat=ioerror) Lmcx_r(r_start:r_end)
!         
               if(ioerror .ne. 0) then
                  write(lupri,*) "Error reading Lmcx_r file"
                  write(lupri,*) "Error code: ", ioerror
                  call quit('Error reading integral file')
               end if
!         
            end do
         end do
!         
      end if
!
!
   end subroutine mlcc3_h_integral_reader
!
!
   subroutine mlcc3_h_L_calc(mbxc,x_size,x_start)
!
!  L calculator
!  Authors: Henrik Koch and Rolf H. Myhre
!  May 2015
!
!  Purpose: calculate L_mbxc = 2(mb|xc) - (xb|mc)
!
      implicit none
!
      integer, intent(in)  :: x_size,x_start
!
      real(dp), dimension(:), pointer, intent(inout)  :: mbxc
!
      integer  :: n_act_int,x,m,int_start,int_end
!
      n_act_int = n_vir_act**2
!
      do x = 1,x_size
         do m = 1,n_occ_act
            if(x+x_start-1 .eq. m) then
               cycle
            else
!
               int_start = n_act_int*(n_occ_act*(x-1)+m-1)+1
               int_end   = int_start + n_act_int - 1
!
               call mlcc3_transposer(mbxc(int_start:int_end),help_v2,n_vir_act,n_vir_act)
!
               call dscal(n_act_int,two,mbxc(int_start:int_end),1)
!
               call daxpy(n_act_int,-one,help_v2,1,mbxc(int_start:int_end),1)
!
            end if
         end do
      end do
!
!
   end subroutine mlcc3_h_L_calc
!
end module mlcc3_h_omega
