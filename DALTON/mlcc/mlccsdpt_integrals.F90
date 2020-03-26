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
module mlccsdpt_integrals
!
!
!  mlccsdpt routines to calculate fock matrix and MO integrals
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
!  Integrals
   real(dp), dimension(:), pointer, private  :: b_D_c_k => null()
   real(dp), dimension(:), pointer, private  :: L_j_c_k => null()
   real(dp), dimension(:), pointer, private  :: j_b_k_c => null()
!   
!  Sorted virtual MO integrals
   real(dp), dimension(:), pointer, private  :: integral_sort => null()
!
!  First transformed array
   real(dp), dimension(:), pointer, private  :: alpha_beta_k_pack => null()
!
!  n_basis x n_basis array, used many places
   real(dp), dimension(:), pointer, private  :: alpha_beta => null()
!
!  n_virtual_active x n_basis
   real(dp), dimension(:), pointer, private  :: b_beta => null()
!
!  n_virtual_active x n_virtual_general
   real(dp), dimension(:), pointer, private  :: b_D => null()
!   
!  n_occupied_active x n_occupied_general
   real(dp), dimension(:), pointer, private  :: L_j => null()
!   
!  n_occupied_active x n_virtual_active
   real(dp), dimension(:), pointer, private  :: j_b => null()
!   
contains
!      
   subroutine mlccsdpt_int_calc
!
!  CCSD(T) integral calculator
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: loop over AO's, read in AO integrals and calculate
!           Fock matrix and MO integrals
!
      implicit none
!
      integer  :: delta
      integer  :: batch, n_batches, batch_size, memory_use, batch_size2
      integer  :: ioerror
      integer  :: D,b,c,k,j,L, k_off
      integer  :: batch_start, batch_end, n_occ_print, n_occ_int
      integer  :: rec_number, n_act_int
      integer  :: work_remains
!
      integer  :: bDck_print_unit
      integer  :: Ljck_print_unit
      integer  :: jbkc_print_unit
!
      real(dp) :: ddot
      real(dp) :: bDcknorm, Ljcknorm, jbkcnorm
!
!     File names
      character(len=12) :: bDck_name
      character(len=12) :: Ljck_name
      character(len=12) :: jbkc_name
!
!   
!     Timing variables
      real     :: time_tot, time_start, time_1, time_2
      real     :: time_int, time_open
      real     :: time_read, time_write, time_close
      real     :: time_zero, time_reord
      real     :: time_bDck, time_Ljck, time_jbkc
!      
      time_int    = 0.0
      time_write  = 0.0
      time_open   = 0.0
      time_read   = 0.0
      time_close  = 0.0
      time_zero   = 0.0
      time_bDck   = 0.0
      time_Ljck   = 0.0
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
      Ljcknorm = zero
      jbkcnorm = zero
!
!
      bDck_name = bDck_file_name
      Ljck_name = Ljck_file_name
      jbkc_name = jbkc_file_name
!
!
!     Some printing variables
!     -----------------------
      n_act_int = n_vir_act*n_vir_act*n_occ_act
      n_occ_int = n_occ_gen*n_occ_act*n_vir_act
!
!     Alllocate arrays independent of batch size
!      
      call work_allocator(ao_integrals_pack,n_ao_ints)
      call work_allocator(alpha_beta,n_basis_2)
      call work_allocator(b_beta,n_vir_act*n_basis)
!      
      call work_allocator(b_D,n_vir_act*n_vir_gen)
      call work_allocator(L_j,n_occ_act*n_occ_gen)
      call work_allocator(j_b,n_occ_act*n_vir_act)
!      
      call work_allocator(integral_sort,n_vir_aag)
!
!
!     Figure out how many k's we can batch over
!
      memory_use = n_basis_2_pack + n_vir_aag + n_occ_int + n_act_int
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
!
!     Fully transformed occupied integrals
!
      call work_allocator(L_j_c_k,n_occ_int*batch_size)
!
!     Fully transformed active integrals
!
      call work_allocator(j_b_k_c,n_act_int*batch_size)
!
!      
!     Partially transformed integrals
!
      call work_allocator(alpha_beta_k_pack,n_basis_2_pack*batch_size)
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
!     (jb|kc)
!
      jbkc_print_unit = new_unit()
      open(unit=jbkc_print_unit, file=jbkc_name, action='write', status='replace', &
     &     access='direct', form='unformatted', recl=dp*n_vir_act**2, iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error opening jbkc file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error with opening file in intermediates')
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
         call dzero(L_j_c_k,n_occ_int*batch_size)
         call dzero(j_b_k_c,n_act_int*batch_size)
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
         if(print_mlcc3 .ge. 5) then
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
            call mlccsdpt_int_reader(delta)
!
            call cpu_time(time_2)
            time_read = time_read + time_2 - time_1
!  
!
!           Integral transform for triples code
!           -----------------------------------
            call cpu_time(time_1)
!
            call mlccsdpt_int_transform(batch,batch_size,batch_size2,delta)
!               
            call cpu_time(time_2)
            time_int = time_int + time_2 - time_1
!  
         end do
!
!
!        Print norms for debug
!
         if(print_mlcc3 .ge. 7) then
!
            bDcknorm = bDcknorm + ddot(n_vir_aag*batch_size,b_D_c_k,1,b_D_c_k,1)
            Ljcknorm = Ljcknorm + ddot(n_occ_int*batch_size,L_j_c_k,1,L_j_c_k,1)
            jbkcnorm = jbkcnorm + ddot(n_act_int*batch_size,j_b_k_c,1,j_b_k_c,1)
!
            write(lupri,*)
            write(lupri,*) 'bDcknorm',bDcknorm
            write(lupri,*) 'Ljcknorm',Ljcknorm
            write(lupri,*) 'jbkcnorm',jbkcnorm
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
!
!              Reorder jbkc integrals to b c, j k order
!
               call cpu_time(time_1)
!
               call jbkc_order(j_b_k_c,integral_sort,n_vir_act,n_vir_gen,&
                               n_occ_act,n_occ_gen,j,k_off,batch_size)
!
               call cpu_time(time_2)
               time_jbkc = time_jbkc + time_2 - time_1
!
!              Write the integrals to disk. 
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
      close(Ljck_print_unit,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing Ljck file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error closing integral file in intermediates')
      end if
!            
!            
      close(jbkc_print_unit,status='keep',iostat=ioerror)
!
      if(ioerror .ne. 0) then
         write(lupri,*) "Error closing jbkc file"
         write(lupri,*) "Error code: ", ioerror
         call quit('Error closing integral file in intermediates')
      end if
!            
      call cpu_time(time_close)
      time_close = time_close - time_1
!  
!
!     ------------------
!     |Free memory used|
!     ------------------
!      
!     Batch dependent allocations
      call work_deallocator(alpha_beta_k_pack)
!     
      call work_deallocator(j_b_k_c)
!      
      call work_deallocator(L_j_c_k)
!      
      call work_deallocator(b_D_c_k)
!
!     Sorted integrals
      call work_deallocator(integral_sort)
!
!     Various help arrays
      call work_deallocator(j_b)
      call work_deallocator(L_j)
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
      time_reord  = time_bDck + time_Ljck + time_jbkc
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
         write(lupri,9999) 'integral trans', time_int
         write(lupri,9999) 'bDck order', time_bDck
         write(lupri,9999) 'Ljck order', time_Ljck
         write(lupri,9999) 'jbkc order', time_jbkc
         write(lupri,9999) 'total order', time_reord
         write(lupri,9999) 'write ints', time_write
         write(lupri,9999) 'close files', time_close
         write(lupri,*)
         write(lupri,*)
!
      end if
!
!
   end subroutine mlccsdpt_int_calc
!
!
   subroutine mlccsdpt_int_reader(delta)
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

   end subroutine mlccsdpt_int_reader
!
!
!      
   subroutine mlccsdpt_int_transform(batch,batch_size,batch_size2,delta)
!
!  CCSD(T) integral transformation routine
!  Authors Henrik Koch and Rolf H. Myhre
!  May 2015
!
!  Purpose: Transform from the ao integrals to the integrals required
!           for the MLCCSD(T) triples and write to disk
!           (alpha beta | gamma delta) -> (bD|ck)
!           (alpha beta | gamma delta) -> (Lj|ck)
!           (alpha beta | gamma delta) -> (jb|kc)
!           
!
      implicit none
!      
      integer, intent(in)     :: delta, batch, batch_size, batch_size2
!      
      integer                 :: k, c
      integer                 :: orb_vir_off, orb_occ_off, orb_batch_off, orb_gen_off
      integer                 :: orb_c_delta, int_off
      integer                 :: alpha_beta_k_off,alpha_beta_k_end
!      
!      
!      
!     Calculate orbital offsets. Inactive occupied MOs come first in ordering
!     
      orb_occ_off       = n_basis*n_occ_inact + 1
      orb_batch_off     = n_basis*(n_occ_inact + batch_size2*(batch-1)) + 1
      orb_vir_off       = n_basis*n_occ + 1
      orb_gen_off       = n_basis*n_gen_inact + 1
!      
!      
!     Start with transforming gamma to the active occupied k, hole and particle
!     -------------------------------------------------------------------------
!
      call dgemm('N','N',n_basis_2_pack,batch_size,n_basis,one,ao_integrals_pack,n_basis_2_pack, &
     &           orb_coefficients(orb_batch_off),n_basis,zero,alpha_beta_k_pack,n_basis_2_pack)
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
         call mlcc3_square_packed(alpha_beta_k_pack(alpha_beta_k_off:alpha_beta_k_end), &
        &                         alpha_beta,n_basis)
!
!        Transform alpha to b for bDck
!        -----------------------------
!
         call dgemm('T','N',n_vir_act,n_basis,n_basis,one,orb_coefficients(orb_vir_off),n_basis, &
        &           alpha_beta,n_basis,zero,b_beta,n_vir_act)
!
!        Transform beta to D
!        -------------------
!
         call dgemm('N','N',n_vir_act,n_vir_gen,n_basis,one,b_beta,n_vir_act, &
        &           orb_coefficients(orb_vir_off),n_basis,zero,b_D,n_vir_act)
!
!        ----------------------------------
!        |Second set of integrals, (Lj|ck)|
!        ----------------------------------
!
!        Transform beta to j for Ljck, reuse alpha_beta
!        ----------------------------------------------
!
         call dgemm('N','N',n_basis,n_occ_act,n_basis,one,alpha_beta,n_basis, &
        &           orb_coefficients(orb_occ_off),n_basis,zero,b_beta,n_basis)
!
!        Transform alpha to L
!        --------------------
!
         call dgemm('T','N', n_occ_gen,n_occ_act,n_basis,one,orb_coefficients(orb_gen_off),n_basis, &
        &           b_beta,n_basis,zero,L_j,n_occ_gen)
!
!
!        ---------------------------------
!        |Fifth set of integrals, (jb|kc)|
!        ---------------------------------
!
!        Transform alpha to j for jbkc, reuse alpha_beta
!        -----------------------------------------------
!
         call dgemm('T','N',n_occ_act,n_basis,n_basis,one,orb_coefficients(orb_occ_off),n_basis, &
        &           alpha_beta,n_basis,zero,b_beta,n_occ_act)
!
!        Transform beta to b for jbkc
!        -----------------------------
!
         call dgemm('N','N',n_occ_act,n_vir_act,n_basis,one,b_beta,n_occ_act, &
        &           orb_coefficients(orb_vir_off),n_basis,zero,j_b,n_occ_act)
!
!
!        Transform delta to c. Start by looping over c index
!        ---------------------------------------------------
!        
         do c = 1,n_vir_act
!
            orb_c_delta = n_basis*(n_occ + c - 1) + delta
!            
!
!           Transform bDck
!           --------------
            int_off    = n_vir_act*n_vir_gen*n_vir_act*(k-1) + n_vir_act*n_vir_gen*(c-1) + 1
!            
            call daxpy(n_vir_act*n_vir_gen,orb_coefficients(orb_c_delta),b_D,1,b_D_c_k(int_off),1)
!
!
!           Transform Ljck
!           --------------
            int_off    = n_occ_gen*n_occ_act*n_vir_act*(k-1) + n_occ_gen*n_occ_act*(c-1) + 1

            call daxpy(n_occ_act*n_occ_gen,orb_coefficients(orb_c_delta),L_j,1,L_j_c_k(int_off),1)
!
!
!           Transform jbkc
!           --------------
            int_off    = n_occ_act*n_vir_act*batch_size*(c-1) + n_occ_act*n_vir_act*(k-1) + 1
!
            call daxpy(n_occ_act*n_vir_act,orb_coefficients(orb_c_delta),j_b,1,j_b_k_c(int_off),1)
!            
         end do
!         
      end do
!
!      
   end subroutine mlccsdpt_int_transform
!
!
end module mlccsdpt_integrals
