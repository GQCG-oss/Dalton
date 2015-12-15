module mlcc3_intermediates
!
!
!  mlcc3 routines to calculate fock matrix and MO integrals
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
   use mlcc_typedef
   use mlcc_block_import
   use mlcc_work
   use mlcc3_active_spaces
   use mlcc3_data
   use mlcc3_reordering
   use mlcc3_various
!
!  AO integrals
   real(dp), dimension(:), pointer, private  :: ao_integrals_pack => null()
!   
!  Unsorted virtual MO integrals
   real(dp), dimension(:), pointer, private  :: b_D_c_k => null()
   real(dp), dimension(:), pointer, private  :: D_b_k_c => null()
!   
!  Unsorted occupied MO integrals
   real(dp), dimension(:), pointer, private  :: L_j_c_k => null()
   real(dp), dimension(:), pointer, private  :: j_L_k_c => null()
!
!  Unsorted active MO integrals
   real(dp), dimension(:), pointer, private  :: j_b_k_c => null()
!   
!  Sorted virtual MO integrals
   real(dp), dimension(:), pointer, private  :: integral_sort => null()
!
!  First transformed array
   real(dp), dimension(:), pointer, private  :: alpha_beta_k_pack_hole => null()
   real(dp), dimension(:), pointer, private  :: alpha_beta_k_pack_part => null()
!
!  n_basis x n_basis array, used many places
   real(dp), dimension(:), pointer, private  :: alpha_beta => null()
!
!  n_virtual_active x n_basis
   real(dp), dimension(:), pointer, private  :: b_beta => null()
!
!  n_virtual_active x n_virtual_general
   real(dp), dimension(:), pointer, private  :: b_D => null()
   real(dp), dimension(:), pointer, private  :: D_b => null()
!   
!  n_occupied_active x n_occupied_general
   real(dp), dimension(:), pointer, private  :: j_L => null()
   real(dp), dimension(:), pointer, private  :: L_j => null()
!   
!  n_occupied_active x n_virtual_active
   real(dp), dimension(:), pointer, private  :: j_b => null()
!   
contains
!      
   subroutine mlcc3_intermediates_calc(resp)
!
!  intermediates calculator
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: loop over AO's, read in AO integrals and calculate
!           Fock matrix and MO integrals
!
      implicit none
!
      logical, intent(in)  :: resp
!
      integer  :: delta
      integer  :: batch, n_batches, batch_size, memory_use, batch_size2
      logical  :: first
      integer  :: integral_size_vir, integral_size_occ
      integer  :: ioerror
      integer  :: D,b,c,k,j,L, k_off
      integer  :: batch_start, batch_end, n_occ_print, n_occ_int
      integer  :: rec_number, n_act_int
      integer  :: work_remains
!
      integer  :: bDck_print_unit, Dbkc_print_unit
      integer  :: Ljck_print_unit, jLkc_print_unit
      integer  :: jbkc_print_unit
!
      real(dp) :: ddot
      real(dp) :: bDcknorm, Dbkcnorm,Ljcknorm,jLkcnorm,jbkcnorm
!
!     File names
      character(len=12) :: bDck_name
      character(len=12) :: Dbkc_name
      character(len=12) :: Ljck_name
      character(len=12) :: jLkc_name
      character(len=12) :: jbkc_name
!
!   
!     Timing variables
      real     :: time_tot, time_start, time_1, time_2
      real     :: time_fock, time_int, time_trans, time_open
      real     :: time_read, time_write, time_close
      real     :: time_zero, time_reord
      real     :: time_bDck, time_Dbkc, time_Ljck, time_jLkc, time_jbkc
!      
      time_fock   = 0.0
      time_int    = 0.0
      time_write  = 0.0
      time_open   = 0.0
      time_read   = 0.0
      time_trans  = 0.0
      time_close  = 0.0
      time_zero   = 0.0
      time_bDck   = 0.0
      time_Dbkc   = 0.0
      time_Ljck   = 0.0
      time_jLkc   = 0.0
      time_jbkc   = 0.0
      time_reord  = 0.0
!      
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
      call cpu_time(time_start)
!
!     Energy or response integrals
!     ----------------------------
!
      bDcknorm = zero
      Dbkcnorm = zero
      Ljcknorm = zero
      jLkcnorm = zero
      jbkcnorm = zero
!
      if(resp) then
!
         bDck_name = bDck_resp_name
         Dbkc_name = Dbkc_resp_name
         Ljck_name = Ljck_resp_name
         jLkc_name = jLkc_resp_name
!
      else
!
         bDck_name = bDck_file_name
         Dbkc_name = Dbkc_file_name
         Ljck_name = Ljck_file_name
         jLkc_name = jLkc_file_name
         jbkc_name = jbkc_file_name
!
      end if
!
!     Some printing variables
!     -----------------------
      n_act_int = n_vir_act*n_vir_act*n_occ_act
      n_occ_int = n_occ_gen*n_occ_act*n_vir_act
!
      first = .true.
!      
!     Allocate space for packed AO integrals
!
      call work_allocator(ao_integrals_pack,n_ao_ints)
!      
!     Alllocate arrays independent of batch size
!
      call work_allocator(alpha_beta,n_basis_2)
      call work_allocator(b_beta,n_vir_act*n_basis)
!      
      call work_allocator(b_D,n_vir_act*n_vir_gen)
      call work_allocator(D_b,n_vir_act*n_vir_gen)
!      
      call work_allocator(L_j,n_occ_act*n_occ_gen)
      call work_allocator(j_L,n_occ_act*n_occ_gen)
!      
      call work_allocator(j_b,n_occ_act*n_vir_act)
!      
      call work_allocator(integral_sort,n_vir_aag)
!
!
!     Figure out how many k's we can batch over
!
      memory_use = 2*(n_vir_aag + n_occ_int + n_basis_2_pack) + n_act_int
!
      work_remains = work_free()
!      
      if (memory_use .gt. work_remains) then
         call quit('Not enough memory for a batch in mlcc3_intermediates')
      end if
!  
      batch_size = min(work_remains/memory_use,n_occ_act)
      n_batches  = (n_occ_act - 1)/batch_size + 1
!
      batch_size2 = batch_size !store for offsets
!
      if(print_mlcc3 .ge. 4) then
         write(lupri,*)
         write(lupri,*) 'work_remains:', work_remains
         write(lupri,*) 'batch_size:  ', batch_size 
         write(lupri,*) 'memory_use:  ', memory_use  
         write(lupri,*) 'n_batches:   ', n_batches   
         write(lupri,*)
      end if
!
!     Allocate arrays dependent on batchsize
!
!     Fully transformed virtual integrals
!
      call work_allocator(b_D_c_k,n_vir_aag*batch_size)
      call work_allocator(D_b_k_c,n_vir_aag*batch_size)
!
!     Fully transformed occupied integrals
!
      call work_allocator(L_j_c_k,n_occ_int*batch_size)
      call work_allocator(j_L_k_c,n_occ_int*batch_size)
!
!
      if(.not. resp) then
!
!        Fully transformed active integrals
!
         call work_allocator(j_b_k_c,n_act_int*batch_size)
!
      end if
!      
!     Partially transformed integrals
!
      call work_allocator(alpha_beta_k_pack_hole,n_basis_2_pack*batch_size)
      call work_allocator(alpha_beta_k_pack_part,n_basis_2_pack*batch_size)
!
!
!     -----------------------------------------------------------------------------
!     Open files for writing integrals outside batch loop, sequential direct access
!     We assume we will not reach the Fortran limit of 2GB record length
!     -----------------------------------------------------------------------------
!
      call cpu_time(time_1)
!      
      n_occ_print = n_occ_gen*n_vir_act
!
!     (bD|ck)
!
      bDck_print_unit = new_unit()
      open(unit=bDck_print_unit, file=bDck_name, action='write', status='replace', &
     &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!  
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening bDck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error with opening file in intermediates')
      end if
!
!     (Db|kc)
!
      Dbkc_print_unit = new_unit()
      open(unit=Dbkc_print_unit, file=Dbkc_name, action='write', status='replace', &
     &     access='direct', form='unformatted', recl=dp*n_vir_aag, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening Dbkc file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error with opening file in intermediates')
      end if
!
!     (Lj|ck)
!
      Ljck_print_unit = new_unit()
      open(unit=Ljck_print_unit, file=Ljck_name, action='write', status='replace', &
     &     access='direct', form='unformatted', recl=dp*n_occ_print, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error with opening file in intermediates')
      end if
!
!     (jL|kc)
!
      jLkc_print_unit = new_unit()
      open(unit=jLkc_print_unit, file=jLkc_name, action='write', status='replace', &
     &     access='direct', form='unformatted', recl=dp*n_occ_print, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening jLkc file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error with opening file in intermediates')
      end if
!
!     (jb|kc)
!
      if(.not. resp) then
         jbkc_print_unit = new_unit()
         open(unit=jbkc_print_unit, file=jbkc_name, action='write', status='replace', &
     &        access='direct', form='unformatted', recl=dp*n_vir_act**2, iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error opening jbkc file"
            write(lupri,*) "Error code: ", ioerror
            call quit('Error with opening file in intermediates')
         end if
!  
      end if
!  
      call cpu_time(time_open)
      time_open = time_open - time_1
!  
!     -------------------
!     |Loop over batches|
!     -------------------
!
      do batch = 1,n_batches
!      
!
!        Set integral arrays to zero
!        ---------------------------
!
         call cpu_time(time_1)
!
         call dzero(b_D_c_k,n_vir_aag*batch_size)
         call dzero(D_b_k_c,n_vir_aag*batch_size)
         call dzero(L_j_c_k,n_occ_int*batch_size)
         call dzero(j_L_k_c,n_occ_int*batch_size)
         if(.not. resp) then
            call dzero(j_b_k_c,n_act_int*batch_size)
         end if
!
         call cpu_time(time_2)
         time_zero = time_zero + time_2 - time_1
!
!
!        Find batch offsets
!        ------------------
!
         batch_start = batch_size*(batch-1) + 1
!
         if (batch .eq. n_batches) then
            batch_size = n_occ_act - batch_size*(n_batches - 1)
         end if
!            
         batch_end = batch_start + batch_size - 1
!         
         if(print_mlcc3 .ge. 4) then
            write(lupri,*)
            write(lupri,*) 'n_batches:   ', n_batches
            write(lupri,*) 'batch:       ', batch
            write(lupri,*) 'batch_size:  ', batch_size 
            write(lupri,*) 'batch_start: ', batch_start
            write(lupri,*) 'batch_end:   ', batch_end
            write(lupri,*)
         end if
!      
!        Loop over the delta index in (alpha beta | gamma delta)
!
         do delta = 1,n_basis
!
!           Read the packed AO integrals into the array ao_integrals_pack
!           -------------------------------------------------------------
!
            call cpu_time(time_1)
!
            call mlcc3_int_reader(delta)
!
            call cpu_time(time_2)
            time_read = time_read + time_2 - time_1
!  
!           Fock matrix calculation. Calculate both t1-transformed and regular versions
!           ---------------------------------------------------------------------------
!
            call cpu_time(time_1)
!
            if (first) then
!
               call mlcc3_fock_calc(delta,resp)
!
            end if
!
            call cpu_time(time_2)
            time_fock = time_fock + time_2 - time_1
!  
!
!           Integral transform for triples code
!           -----------------------------------
            call cpu_time(time_1)
!
            if(.not. resp) then !Energy integrals
               call mlcc3_integral_transform(batch,batch_size,batch_size2,delta)
            else !response integrals
               call mlcc3_integral_trans_resp(batch,batch_size,batch_size2,delta)
            end if
!               
            call cpu_time(time_2)
            time_int = time_int + time_2 - time_1
!  
         end do
!
!        Do not calculate Fock matrix again
         first = .false.
!
!
!        Print norms for debug
!
         if(print_mlcc3 .ge. 7) then
!
            bDcknorm = bDcknorm + ddot(n_vir_aag*batch_size,b_D_c_k,1,b_D_c_k,1)
            Dbkcnorm = Dbkcnorm + ddot(n_vir_aag*batch_size,D_b_k_c,1,D_b_k_c,1)
            Ljcknorm = Ljcknorm + ddot(n_occ_int*batch_size,L_j_c_k,1,L_j_c_k,1)
            jLkcnorm = jLkcnorm + ddot(n_occ_int*batch_size,j_L_k_c,1,j_L_k_c,1)
            if(.not. resp) then
               jbkcnorm = jbkcnorm + ddot(n_act_int*batch_size,j_b_k_c,1,j_b_k_c,1)
            end if
!
            write(lupri,*)
            write(lupri,*) 'bDcknorm',bDcknorm
            write(lupri,*) 'Dbkcnorm',Dbkcnorm
            write(lupri,*) 'Ljcknorm',Ljcknorm
            write(lupri,*) 'jLkcnorm',jLkcnorm
            if(.not. resp) then
               write(lupri,*) 'jbkcnorm',jbkcnorm
            end if
            write(lupri,*)
!
         end if
!
!
!        -----------------------------------------
!        |Reorder the integrals and print to file|
!        -----------------------------------------
!
         do k = batch_start,batch_end
!         
            k_off = k - batch_start + 1
!            
!           Reorder bDck integrals to D b c, k order
!
            call cpu_time(time_1)
!
            call bDck_order(b_D_c_k,integral_sort,n_vir_act,n_vir_gen,k_off,batch_size)
!
            call cpu_time(time_2)
            time_bDck = time_bDck + time_2 - time_1
!
!           Write the integrals to disk. 
!
            call cpu_time(time_1)
!
            write(bDck_print_unit,rec=k,iostat=ioerror) integral_sort(1:n_vir_aag)
!            
            call cpu_time(time_2)
            time_write = time_write + time_2 - time_1
!
            if(ioerror .ne. 0) then
               write(lupri,*) "Error writing to bDck file"
               write(lupri,*) "Error code: ", ioerror
               call quit('Error writing integrals to file in intermediates')
            end if
!            
!           Reorder Dbkc integrals to b c D, k order
!
            call cpu_time(time_1)
!
            call Dbkc_order(D_b_k_c,integral_sort,n_vir_act,n_vir_gen,k_off,batch_size)
!
            call cpu_time(time_2)
            time_Dbkc = time_Dbkc + time_2 - time_1
!
!           Write the integrals to disk. 
!
            call cpu_time(time_1)
!
            write(Dbkc_print_unit,rec=k,iostat=ioerror) integral_sort(1:n_vir_aag)
!            
            call cpu_time(time_2)
            time_write = time_write + time_2 - time_1
!
            if(ioerror .ne. 0) then
               write(lupri,*) "Error writing to Dbkc file"
               write(lupri,*) "Error code: ", ioerror
               call quit('Error writing integrals to file in intermediates')
            end if
!            
            do j = 1,n_occ_act
!
!
!              Find record number
!
               rec_number = n_occ_act*(k-1) + j
!               
!              Reorder Ljck integrals to L c, j k order
!               
               call cpu_time(time_1)
!
               call Ljck_order(L_j_c_k,integral_sort,n_vir_act,n_vir_gen,&
              &                n_occ_act,n_occ_gen,j,k_off,batch_size)
!
               call cpu_time(time_2)
               time_Ljck = time_Ljck + time_2 - time_1
!
!              Write the integrals to disk. 
!
               call cpu_time(time_1)
!
               write(Ljck_print_unit,rec=rec_number,iostat=ioerror) integral_sort(1:n_occ_print)
!            
               call cpu_time(time_2)
               time_write = time_write + time_2 - time_1
!
               if(ioerror .ne. 0) then
                  write(lupri,*) "Error writing to Ljck file"
                  write(lupri,*) "Error code: ", ioerror
                  call quit('Error writing integrals to file in intermediates')
               end if
!            
!              Reorder jLkc integrals to c L, j k order
!
               call cpu_time(time_1)
!
               call jLkc_order(j_L_k_c,integral_sort,n_vir_act,n_vir_gen,&
              &                n_occ_act,n_occ_gen,j,k_off,batch_size)
!
               call cpu_time(time_2)
               time_jLkc = time_jLkc + time_2 - time_1
!         
!              Write the integrals to disk. 
!
               call cpu_time(time_1)
!
               write(jLkc_print_unit,rec=rec_number,iostat=ioerror) integral_sort(1:n_occ_print)
!            
               call cpu_time(time_2)
               time_write = time_write + time_2 - time_1
!
               if(ioerror .ne. 0) then
                  write(lupri,*) "Error writing to jLkc file"
                  write(lupri,*) "Error code: ", ioerror
                  call quit('Error writing integrals to file in intermediates')
               end if
!            
!
               if(.not. resp) then
!
!                 Reorder jbkc integrals to b c, j k order
!
                  call cpu_time(time_1)
!
                  call jbkc_order(j_b_k_c,integral_sort,n_vir_act,n_vir_gen,&
                 &                n_occ_act,n_occ_gen,j,k_off,batch_size)
!
                  call cpu_time(time_2)
                  time_jbkc = time_jbkc + time_2 - time_1
!
!                 Write the integrals to disk. 
!
                  call cpu_time(time_1)
!
                  write(jbkc_print_unit,rec=rec_number,iostat=ioerror) integral_sort(1:n_vir_act**2)
!            
                  call cpu_time(time_2)
                  time_write = time_write + time_2 - time_1
!
                  if(ioerror .ne. 0) then
                     write(lupri,*) "Error writing to jbkc file"
                     write(lupri,*) "Error code: ", ioerror
                     call quit('Error writing integrals to file in intermediates')
                  end if
!            
               end if
!
            end do
!            
         end do
!            
         call cpu_time(time_write)
         time_write = time_write - time_1
!  
!
      end do
!
!     ----------------------
!     |Close integral files|
!     ----------------------
!
!  
      call cpu_time(time_1)
!
      close(bDck_print_unit,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing bDck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error closing integral file in intermediates')
      end if
!            
      close(Dbkc_print_unit,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing Dbkc file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error closing integral file in intermediates')
      end if
!            
      close(Ljck_print_unit,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error closing integral file in intermediates')
      end if
!            
      close(jLkc_print_unit,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing jLkc file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error closing integral file in intermediates')
      end if
!            
!            
      if(.not. resp) then
!            
         close(jbkc_print_unit,status='keep',iostat=ioerror)
!
         if(ioerror .ne. 0) then
            write(lupri,*) "Error closing jbkc file"
            write(lupri,*) "Error code: ", ioerror
            call quit('Error closing integral file in intermediates')
         end if
!            
      end if
!            
      call cpu_time(time_close)
      time_close = time_close - time_1
!  
!
!     -------------------------------------------------------------------------------
!     |Transform Fock matrices to MO basis. Use the alpha_beta array as a help array|
!     |Note that the Fock matrix is calculated as the transpose. Easier to use      |
!     -------------------------------------------------------------------------------
!
!
      call cpu_time(time_1)
!
      call mlcc3_fock_ao_mo(resp)
!
      call cpu_time(time_trans)
      time_trans = time_trans - time_1
!
!
!     ------------------
!     |Free memory used|
!     ------------------
!      
!     Batch dependent allocations
      call work_deallocator(alpha_beta_k_pack_part)
      call work_deallocator(alpha_beta_k_pack_hole)
!     
      if(.not. resp) then
         call work_deallocator(j_b_k_c)
      end if
!      
      call work_deallocator(j_L_k_c)
      call work_deallocator(L_j_c_k)
!      
      call work_deallocator(D_b_k_c)
      call work_deallocator(b_D_c_k)
!
!     Sorted integrals
      call work_deallocator(integral_sort)
!
!     Various help arrays
      call work_deallocator(j_b)
      call work_deallocator(j_L)
      call work_deallocator(L_j)
!      
      call work_deallocator(D_b)
      call work_deallocator(b_D)
!      
      call work_deallocator(b_beta)
      call work_deallocator(alpha_beta)
      call work_deallocator(ao_integrals_pack)
!      
!
!
!     Total timing
      call cpu_time(time_tot)
!      
      time_tot = time_tot - time_start
      time_reord  = time_bDck + time_Dbkc + time_Ljck + time_jLkc + time_jbkc
!
!     Print timings
!     -------------
!
      if (print_mlcc3 .ge. 4) then
!
         write(lupri,*)
         write(lupri,*)
         write(lupri,*) 'Timings from mlcc3_int after loop'
         write(lupri,9999) 'Total', time_tot
         write(lupri,9999) 'open files', time_open
         write(lupri,9999) 'zero arrays', time_zero
         write(lupri,9999) 'read AO ints', time_read
         write(lupri,9999) 'AO Fock mat', time_fock 
         write(lupri,9999) 'integral trans', time_int
         write(lupri,9999) 'bDck order', time_bDck
         write(lupri,9999) 'Dbkc order', time_Dbkc
         write(lupri,9999) 'Ljck order', time_Ljck
         write(lupri,9999) 'jLkc order', time_jLkc
         write(lupri,9999) 'jbkc order', time_jbkc
         write(lupri,9999) 'total order', time_reord
         write(lupri,9999) 'write ints', time_write
         write(lupri,9999) 'close files', time_close
         write(lupri,9999) 'Fock transform', time_trans
         write(lupri,*)
         write(lupri,*)
!
      end if
!
!
   end subroutine mlcc3_intermediates_calc
!
   subroutine mlcc3_int_reader(delta)
!
!  integral reader
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: read in (alpha beta | gamma delta) integrals for fixed delta
!           and alpha >= beta
!
      implicit none
!
      integer, intent(in)              :: delta
      integer                          :: dumdel         = 1
      logical                          :: direct_logic   = .false.
      real(dp), pointer, dimension(:)  :: work_point_end => null()
      integer                          :: work_remains
!
      work_point_end => work_end()
      work_remains = work_free()
!      
      call ccrdao(ao_integrals_pack,delta,dumdel,work_point_end,work_remains,dumdel,direct_logic)

   end subroutine mlcc3_int_reader
!
!
   subroutine mlcc3_fock_calc(delta,resp)
!
!  Fock matrix calculator
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: Calculate the various fock matrices using the alpha_beta array as a help array
!
!
      implicit none
!      
      integer, intent(in)  :: delta
      logical, intent(in)  :: resp 
      integer              :: gamm, gamm_off, gamm_off_pack, delta_off, gamma_delta
      real(dp)             :: ddot
!      
!      
!      
!     Calculate delta offset
      delta_off = n_basis*(delta -1) + 1
!
      do gamm = 1,n_basis
!
!        Calculate gamma offset
         gamm_off       = n_basis*(gamm-1) + 1
         gamm_off_pack  = n_basis_2_pack*(gamm-1)  + 1
         gamma_delta    = n_basis*(gamm-1) + delta
!         
!        square up integral(alpha,beta)_gamm,delta
!
         call mlcc3_square_packed(ao_integrals_pack(gamm_off_pack:n_basis_2_pack),alpha_beta,n_basis)
!         
!        Add coloumb contribution. F_gamma,delta += 2(alpa beta | gamma delta)*D_alpha,beta
!
         mo_fock_mat(gamma_delta) = mo_fock_mat(gamma_delta) + &
        &   2*ddot(n_basis_2,alpha_beta,1,ao_density,1)
!     
         mo_fock_mat_t1(gamma_delta) = mo_fock_mat_t1(gamma_delta) + &
        &   2*ddot(n_basis_2,alpha_beta,1,ao_density_t1,1)
!     
         if(resp) then
            mo_fock_mat_c1(gamma_delta) = mo_fock_mat_c1(gamma_delta) + &
           &   2*ddot(n_basis_2,alpha_beta,1,ao_density_c1,1)
         end if
!     
!        Subtract exchange contribution. F_gamma,delta -= (gamma beta | alpha delta)*D_alpha,beta
!
         call dgemv('N',n_basis,n_basis,-one,alpha_beta,n_basis, &
        &           ao_density(gamm_off),1,one,mo_fock_mat(delta_off),1)
!
         call dgemv('N',n_basis,n_basis,-one,alpha_beta,n_basis, &
        &           ao_density_t1(gamm_off),1,one,mo_fock_mat_t1(delta_off),1)
!
         if(resp) then
            call dgemv('N',n_basis,n_basis,-one,alpha_beta,n_basis, &
           &           ao_density_c1(gamm_off),1,one,mo_fock_mat_c1(delta_off),1)
         end if
!     
      end do
!
!         
   end subroutine mlcc3_fock_calc
!
!
   subroutine mlcc3_fock_ao_mo(resp)
!
!  transform AO Fock matrices
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: Transform Fock matrices from AO to MO basis
!           
!
      implicit none
!      
      logical, intent(in)  :: resp 
!
!
      if(resp) then
!
         call dgemm('T','N',n_orbitals,n_basis,n_basis,one,lambda_hole,n_basis, &
        &        mo_fock_mat_c1,n_basis,zero,alpha_beta,n_orbitals)
         call dgemm('N','N',n_orbitals,n_orbitals,n_basis,one,alpha_beta,n_orbitals, &
        &        lambda_part,n_basis,zero,mo_fock_mat_c1,n_orbitals)
!
         call dgemm('T','N',n_orbitals,n_basis,n_basis,one,lambda_hole,n_basis, &
        &        mo_fock_mat_t1,n_basis,zero,alpha_beta,n_orbitals)
         call dgemm('N','N',n_orbitals,n_orbitals,n_basis,one,alpha_beta,n_orbitals, &
        &        lambda_part_resp,n_basis,one,mo_fock_mat_c1,n_orbitals)
!
         call dgemm('T','N',n_orbitals,n_basis,n_basis,one,lambda_hole_resp,n_basis, &
        &        mo_fock_mat_t1,n_basis,zero,alpha_beta,n_orbitals)
         call dgemm('N','N',n_orbitals,n_orbitals,n_basis,one,alpha_beta,n_orbitals, &
        &        lambda_part,n_basis,one,mo_fock_mat_c1,n_orbitals)
!
      end if
!
!
      call dgemm('T','N',n_orbitals,n_basis,n_basis,one,lambda_hole,n_basis, &
     &        mo_fock_mat_t1,n_basis,zero,alpha_beta,n_orbitals)
      call dgemm('N','N',n_orbitals,n_orbitals,n_basis,one,alpha_beta,n_orbitals, &
     &        lambda_part,n_basis,zero,mo_fock_mat_t1,n_orbitals)
!
!      
      call dgemm('T','N',n_orbitals,n_basis,n_basis,one,orb_coefficients,n_basis, &
     &        mo_fock_mat,n_basis,zero,alpha_beta,n_orbitals)
      call dgemm('N','N',n_orbitals,n_orbitals,n_basis,one,alpha_beta,n_orbitals, &
     &        orb_coefficients,n_basis,zero,mo_fock_mat,n_orbitals)
!
   end subroutine mlcc3_fock_ao_mo
!
!
!      
   subroutine mlcc3_integral_transform(batch,batch_size,batch_size2,delta)
!
!  integral transformation routine
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: Transform from the ao integrals to the integrals required
!           for the MLCC3 triples and write to disk
!           (alpha beta | gamma delta) -> (bD|ck)
!           (alpha beta | gamma delta) -> (Db|kc)
!           (alpha beta | gamma delta) -> (Lj|ck)
!           (alpha beta | gamma delta) -> (jL|kc)
!           (alpha beta | gamma delta) -> (jb|kc)
!           
!
      implicit none
!      
      integer, intent(in)     :: delta, batch, batch_size, batch_size2
!      
      integer                 :: k, c
      integer                 :: lambda_vir_off, lambda_occ_off, lambda_batch_off, lambda_gen_off
      integer                 :: lambda_c_delta, int_off
      integer                 :: alpha_beta_k_off,alpha_beta_k_end
!      
!      
!      
!     Calculate lambda offsets. Inactive occupied MOs come first in ordering
!     
      lambda_occ_off       = n_basis*n_occ_inact + 1
      lambda_batch_off     = n_basis*(n_occ_inact + batch_size2*(batch-1)) + 1
      lambda_vir_off       = n_basis*n_occ + 1
      lambda_gen_off       = n_basis*n_gen_inact + 1
!      
!      
!     Start with transforming gamma to the active occupied k, hole and particle
!     -------------------------------------------------------------------------
!
      call dgemm('N','N',n_basis_2_pack,batch_size,n_basis,one,ao_integrals_pack,n_basis_2_pack, &
     &           lambda_hole(lambda_batch_off),n_basis,zero,alpha_beta_k_pack_hole,n_basis_2_pack)
!  
      call dgemm('N','N',n_basis_2_pack,batch_size,n_basis,one,ao_integrals_pack,n_basis_2_pack, &
     &           lambda_part(lambda_batch_off),n_basis,zero,alpha_beta_k_pack_part,n_basis_2_pack)
!  
!      
!     Loop over the k indices in the batch
!     ------------------------------------
!      
      do k = 1,batch_size
!
         alpha_beta_k_off = n_basis_2_pack*(k-1) + 1
         alpha_beta_k_end = n_basis_2_pack*k
!
!        ---------------------------------
!        |First set of integrals, (bD|ck)|
!        ---------------------------------
!
!        Square up hole transformed alpha beta
!        -------------------------------------
!
         call mlcc3_square_packed(alpha_beta_k_pack_hole(alpha_beta_k_off:alpha_beta_k_end), &
        &                         alpha_beta,n_basis)
!
!        Transform alpha to b for bDck
!        -----------------------------
!
         call dgemm('T','N',n_vir_act,n_basis,n_basis,one,lambda_part(lambda_vir_off),n_basis, &
        &           alpha_beta,n_basis,zero,b_beta,n_vir_act)
!
!        Transform beta to D
!        -------------------
!
         call dgemm('N','N',n_vir_act,n_vir_gen,n_basis,one,b_beta,n_vir_act, &
        &           lambda_hole(lambda_vir_off),n_basis,zero,b_D,n_vir_act)
!
!        ----------------------------------
!        |Second set of integrals, (Lj|ck)|
!        ----------------------------------
!
!        Transform beta to j for Ljck, reuse alpha_beta
!        ----------------------------------------------
!
         call dgemm('N','N',n_basis,n_occ_act,n_basis,one,alpha_beta,n_basis, &
        &           lambda_hole(lambda_occ_off),n_basis,zero,b_beta,n_basis)
!
!        Transform alpha to L
!        --------------------
!
         call dgemm('T','N', n_occ_gen,n_occ_act,n_basis,one,lambda_part(lambda_gen_off),n_basis, &
        &           b_beta,n_basis,zero,L_j,n_occ_gen)
!
!
!        ---------------------------------
!        |Third set of integrals, (Db|kc)|
!        ---------------------------------
!
!        Square up particle transformed alpha beta
!        -----------------------------------------
!
         call mlcc3_square_packed(alpha_beta_k_pack_part(alpha_beta_k_off:alpha_beta_k_end), &
        &                         alpha_beta,n_basis)
!
!        Transform beta to b for Dbkc
!        -----------------------------
!
         call dgemm('N','N',n_basis,n_vir_act,n_basis,one,alpha_beta,n_basis, &
        &           lambda_hole(lambda_vir_off),n_basis,zero,b_beta,n_basis)
!
!        Transform alpha to D
!        --------------------
!
         call dgemm('T','N', n_vir_gen,n_vir_act,n_basis,one,lambda_part(lambda_vir_off),n_basis, &
        &           b_beta,n_basis,zero,D_b,n_vir_gen)
!
!
!        ----------------------------------
!        |Fourth set of integrals, (jL|kc)|
!        ----------------------------------
!
!        Transform alpha to j for jLkc, reuse alpha_beta
!        -----------------------------------------------
!
         call dgemm('T','N',n_occ_act,n_basis,n_basis,one,lambda_part(lambda_occ_off),n_basis, &
        &           alpha_beta,n_basis,zero,b_beta,n_occ_act)
!
!        Transform beta to L
!        -------------------
!
         call dgemm('N','N',n_occ_act,n_occ_gen,n_basis,one,b_beta,n_occ_act, &
        &           lambda_hole(lambda_gen_off),n_basis,zero,j_L,n_occ_act)
!
!
!        ---------------------------------
!        |Fifth set of integrals, (jb|kc)|
!        ---------------------------------
!
!
!        Transform alpha to j for jbkc, reuse alpha_beta
!        -----------------------------------------------
!
         call dgemm('T','N',n_occ_act,n_basis,n_basis,one,lambda_part(lambda_occ_off),n_basis, &
        &           alpha_beta,n_basis,zero,b_beta,n_occ_act)
!
!        Transform beta to b for jbkc
!        -----------------------------
!
         call dgemm('N','N',n_occ_act,n_vir_act,n_basis,one,b_beta,n_occ_act, &
        &           lambda_hole(lambda_vir_off),n_basis,zero,j_b,n_occ_act)
!
!
!        Transform delta to c. Start by looping over c index
!        ---------------------------------------------------
!        
         do c = 1,n_vir_act
!
            lambda_c_delta = n_basis*(n_occ + c - 1) + delta
!            
!
!           Transform bDck
!           --------------
            int_off    = n_vir_act*n_vir_gen*n_vir_act*(k-1) + n_vir_act*n_vir_gen*(c-1) + 1
!            
            call daxpy(n_vir_act*n_vir_gen,lambda_part(lambda_c_delta),b_D,1,b_D_c_k(int_off),1)
!
!
!           Transform Ljck
!           --------------
            int_off    = n_occ_gen*n_occ_act*n_vir_act*(k-1) + n_occ_gen*n_occ_act*(c-1) + 1

            call daxpy(n_occ_act*n_occ_gen,lambda_part(lambda_c_delta),L_j,1,L_j_c_k(int_off),1)
!
!
!           Transform Dbkc
!           --------------
            int_off    = n_vir_gen*n_vir_act*batch_size*(c-1) + n_vir_gen*n_vir_act*(k-1) + 1
!            
            call daxpy(n_vir_act*n_vir_gen,lambda_hole(lambda_c_delta),D_b,1,D_b_k_c(int_off),1)
!            
!
!           Transform jLkc
!           --------------
            int_off    = n_occ_act*n_occ_gen*batch_size*(c-1) + n_occ_act*n_occ_gen*(k-1) + 1

            call daxpy(n_occ_act*n_occ_gen,lambda_hole(lambda_c_delta),j_L,1,j_L_k_c(int_off),1)
!
!
!           Transform jbkc
!           --------------
            int_off    = n_occ_act*n_vir_act*batch_size*(c-1) + n_occ_act*n_vir_act*(k-1) + 1
!
            call daxpy(n_occ_act*n_vir_act,lambda_hole(lambda_c_delta),j_b,1,j_b_k_c(int_off),1)
!            
         end do
!         
      end do
!      
!      
!
!      
   end subroutine mlcc3_integral_transform
!
!
   subroutine mlcc3_integral_trans_resp(batch,batch_size,batch_size2,delta)
!
!  integral transformation routine
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: Transform from the ao integrals to the integrals required
!           for the MLCC3 triples excitation energies and write to disk
!           (alpha beta | gamma delta) -> (b'D|ck) + (bD|c'k) + (bD|ck')
!           (alpha beta | gamma delta) -> (Db'|kc)
!           (alpha beta | gamma delta) -> (Lj'|ck) + (Lj|c'k) + (Lj|ck')
!           (alpha beta | gamma delta) -> (j'L|kc)
!           
!
      implicit none
!      
      integer, intent(in)     :: delta, batch, batch_size, batch_size2
!      
      integer                 :: k, c, i
      integer                 :: lambda_vir_off, lambda_occ_off, lambda_batch_off, lambda_gen_off
      integer                 :: lambda_c_delta, int_off
      integer                 :: alpha_beta_k_off, alpha_beta_k_end
!      
      real(dp), dimension(:), pointer  :: lambda_hole_1  => null()
!   
      real(dp), dimension(:), pointer  :: lambda_hole_2  => null()
      real(dp), dimension(:), pointer  :: lambda_part_2  => null()
!   
      real(dp), dimension(:), pointer  :: lambda_part_3  => null()
!   
!      
!      
!     Calculate lambda offsets. Inactive occupied MOs come first in ordering
!     
      lambda_occ_off       = n_basis*n_occ_inact + 1
      lambda_batch_off     = n_basis*(n_occ_inact + batch_size2*(batch-1)) + 1
      lambda_vir_off       = n_basis*n_occ + 1
      lambda_gen_off       = n_basis*n_gen_inact + 1
!      
!
      do i=1,3 !Loop to get all contributions to |ck) integrals
!      
!
!        Set up lambda matrices for different contributions
!        --------------------------------------------------
!
         if(i .eq. 1) then !b and j C_1 transformed
!
            lambda_hole_1 => lambda_hole
            lambda_hole_2 => lambda_hole_resp
            lambda_part_2 => lambda_part_resp
            lambda_part_3 => lambda_part
!
         else if(i .eq. 2) then !k C_1 transformed
!
            lambda_hole_1 => lambda_hole
            lambda_hole_2 => lambda_hole
            lambda_part_2 => lambda_part
            lambda_part_3 => lambda_part_resp
!
         else !c C_1 transformed
!
            lambda_hole_1 => lambda_hole_resp
            lambda_hole_2 => lambda_hole
            lambda_part_2 => lambda_part
            lambda_part_3 => lambda_part
!
         end if
!
!        Start with transforming gamma to the active occupied k, hole and particle
!        -------------------------------------------------------------------------
!
         if(i .ne. 2) then ! same in 1 and 2
!
            call dgemm('N','N',n_basis_2_pack,batch_size,n_basis,one,ao_integrals_pack,n_basis_2_pack, &
           &           lambda_hole_1(lambda_batch_off),n_basis,zero,alpha_beta_k_pack_hole,n_basis_2_pack)
!
         end if
!  
         if(i .eq. 1) then
!      
            call dgemm('N','N',n_basis_2_pack,batch_size,n_basis,one,ao_integrals_pack,n_basis_2_pack, &
           &           lambda_part(lambda_batch_off),n_basis,zero,alpha_beta_k_pack_part,n_basis_2_pack)
!      
         end if
!  
!      
!        Loop over the k indices in the batch
!        ------------------------------------
!      
         do k = 1,batch_size
!
            alpha_beta_k_off = n_basis_2_pack*(k-1) + 1
            alpha_beta_k_end = n_basis_2_pack*k
!
!           ---------------------------------
!           |First set of integrals, (bD|ck)|
!           ---------------------------------
!
!           Square up hole transformed alpha beta
!           -------------------------------------
!
            call mlcc3_square_packed(alpha_beta_k_pack_hole(alpha_beta_k_off:alpha_beta_k_end), &
        &                         alpha_beta,n_basis)
!
!           Transform alpha to b for bDck
!           -----------------------------
!
            call dgemm('T','N',n_vir_act,n_basis,n_basis,one,lambda_part_2(lambda_vir_off),n_basis, &
           &           alpha_beta,n_basis,zero,b_beta,n_vir_act)
!
!           Transform beta to D
!           -------------------
!
            call dgemm('N','N',n_vir_act,n_vir_gen,n_basis,one,b_beta,n_vir_act, &
           &           lambda_hole(lambda_vir_off),n_basis,zero,b_D,n_vir_act)
!
!           ----------------------------------
!           |Second set of integrals, (Lj|ck)|
!           ----------------------------------
!
!           Transform beta to j for Ljck, reuse alpha_beta
!           ----------------------------------------------
!
            call dgemm('N','N',n_basis,n_occ_act,n_basis,one,alpha_beta,n_basis, &
           &           lambda_hole_2(lambda_occ_off),n_basis,zero,b_beta,n_basis)
!
!           Transform alpha to L
!           --------------------
!
            call dgemm('T','N', n_occ_gen,n_occ_act,n_basis,one,lambda_part(lambda_gen_off),n_basis, &
           &           b_beta,n_basis,zero,L_j,n_occ_gen)
!
!
            if(i .eq. 1) then
!
!              ---------------------------------
!              |Third set of integrals, (Db|kc)|
!              ---------------------------------
!
!              Square up particle transformed alpha beta
!              -----------------------------------------
!
               call mlcc3_square_packed(alpha_beta_k_pack_part(alpha_beta_k_off:alpha_beta_k_end), &
              &                         alpha_beta,n_basis)
!
!              Transform beta to b for Dbkc
!              -----------------------------
!
               call dgemm('N','N',n_basis,n_vir_act,n_basis,one,alpha_beta,n_basis, &
              &           lambda_hole(lambda_vir_off),n_basis,zero,b_beta,n_basis)
!
!              Transform alpha to D
!              --------------------
!
               call dgemm('T','N', n_vir_gen,n_vir_act,n_basis,one,lambda_part_resp(lambda_vir_off), &
              &           n_basis,b_beta,n_basis,zero,D_b,n_vir_gen)
!
!
!              ----------------------------------
!              |Fourth set of integrals, (jL|kc)|
!              ----------------------------------
!
!              Transform alpha to j for jLkc, reuse alpha_beta
!              -----------------------------------------------
!
               call dgemm('T','N',n_occ_act,n_basis,n_basis,one,lambda_part(lambda_occ_off), &
              &           n_basis,alpha_beta,n_basis,zero,b_beta,n_occ_act)
!
!              Transform beta to L
!              -------------------
!
               call dgemm('N','N',n_occ_act,n_occ_gen,n_basis,one,b_beta,n_occ_act, &
              &           lambda_hole_resp(lambda_gen_off),n_basis,zero,j_L,n_occ_act)
!
            end if
!
!
!
!           Transform delta to c. Start by looping over c index
!           ---------------------------------------------------
!        
            do c = 1,n_vir_act
!
               lambda_c_delta = n_basis*(n_occ + c - 1) + delta
!            
!
!              Transform bDck
!              --------------
               int_off    = n_vir_act*n_vir_gen*n_vir_act*(k-1) + n_vir_act*n_vir_gen*(c-1) + 1
!            
               call daxpy(n_vir_act*n_vir_gen,lambda_part_3(lambda_c_delta),b_D,1,b_D_c_k(int_off),1)
!
!
!              Transform Ljck
!              --------------
               int_off    = n_occ_gen*n_occ_act*n_vir_act*(k-1) + n_occ_gen*n_occ_act*(c-1) + 1

               call daxpy(n_occ_act*n_occ_gen,lambda_part_3(lambda_c_delta),L_j,1,L_j_c_k(int_off),1)
!
!
               if(i .eq. 1) then
!
!                 Transform Dbkc
!                 --------------
                  int_off    = n_vir_gen*n_vir_act*batch_size*(c-1) + n_vir_gen*n_vir_act*(k-1) + 1
!            
                  call daxpy(n_vir_act*n_vir_gen,lambda_hole(lambda_c_delta),D_b,1,D_b_k_c(int_off),1)
!            
!
!                 Transform jLkc
!                 --------------
                  int_off    = n_occ_act*n_occ_gen*batch_size*(c-1) + n_occ_act*n_occ_gen*(k-1) + 1

                  call daxpy(n_occ_act*n_occ_gen,lambda_hole(lambda_c_delta),j_L,1,j_L_k_c(int_off),1)
!
               end if
!
            end do
!         
         end do
!      
      end do
!      
   end subroutine mlcc3_integral_trans_resp
!
!
end module mlcc3_intermediates
