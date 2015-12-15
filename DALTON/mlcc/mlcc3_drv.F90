subroutine mlcc3_drv(omega1_in,omega2_in,c1_in,c2_in,frequency,response,work,int_work,lwork)
!
!
!  mlcc3 driver
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Calculate triples contribution in MLCC3 or CC3 calculations energy or excitation energies.
!
!  omega1_in and omega2_in contains the packed omega or result vectors from DALTON CCSD for ground state
!  or excitation energies respectively
!
!  C1_in and C2_in contains the trial vectors for excitation calcuations. They are not used otherwise
!
!  Frequency contains the resulting eigenvalues. Needed in excitation calculations and only used there
!
!  Response is a logical that decides the behaviour, ground state or excitation energy
!
!  work and int_work are both work(kend) from dalton. Two needed in case of integer arrays. 
!  Only int is used in current version
!
!  Note that the diagonal of the trial vector must be scaled by two and the diagonal of the result vector
!  must be scaled by a half before entering and after leaving this rutine respectively
!
   use mlcc_work
   use mlcc_typedef
   use mlcc_block_import
   use mlcc3_init
   use mlcc3_intermediates
   use mlcc3_omega
   use mlcc3_h_omega
   use mlcc3_active_spaces
!
   implicit none
!
   integer, intent(in)                    :: lwork
   real(dp), dimension(lwork), intent(in) :: work
   integer, dimension(lwork), intent(in)  :: int_work
!
!  Energy or response calculation
   logical, intent(in)                    :: response
   real(dp), intent(in)                   :: frequency
!
!  Omega vctors from CCSD
   real(dp), dimension(*), intent(in)     :: omega1_in
   real(dp), dimension(*), intent(in)     :: omega2_in
!
!  Response trial vectors from CCSD
   real(dp), dimension(*), intent(in)     :: c1_in
   real(dp), dimension(*), intent(in)     :: c2_in
!
   integer  :: resp_control, work_required, nn, work_avail
   logical  :: resp, highmem
!
!  Timing variables
   real     :: time_tot, time_ws, time_read, time_os, time_imp
   real     :: time_start, time_1, time_int, time_ome, time_de, time_lambda
!      
   integer  :: wr, w1, w2
   real     :: wall_time
!      
9999 format(7x,'Time used in',2x,A22,2x,': ',f14.4,' seconds')
!      
!
   call system_clock(count_rate=wr)
!
   call system_clock(w1)
   call cpu_time(time_start)
!
   resp_option = response
   highmem     = .false.
!
   if(resp_option) then
      freq     = frequency
      nn       = 3
   else
      freq     = zero
      nn       = 3
   end if
!
!
!  Set up work information
!  -----------------------
   call cpu_time(time_1)
!
   call work_setup(work,int_work,lwork)
!
   call cpu_time(time_ws)
   time_ws = time_ws - time_1
!
!
!  Read information from DALTON and initialize variables.
!  ------------------------------------------------------
   call cpu_time(time_1)
!
   call sirius_reader
!
   call cpu_time(time_read)
   time_read = time_read - time_1
!
!
!  Setup pointers from Dalton
!  --------------------------
   call cpu_time(time_1)
!
   call omega_setup(omega1_in,omega2_in,c1_in,c2_in)
!   
   call cpu_time(time_os)
   time_os = time_os - time_1
!
!
!  Read AO integrals and calculate lambda matrices
!  -----------------------------------------------
   call cpu_time(time_1)
!
   call mlcc3_lambda_ao
!
   call cpu_time(time_lambda)
   time_lambda = time_lambda - time_1
!
!
!  Import active space information from cc_active
!  ----------------------------------------------
   call cpu_time(time_1)
!
   call center_h_import
!   
!
!  Import cvs input information
!  ----------------------------
!
!  Not implemented in this version
!   call cvsexci_import
!   
   call cpu_time(time_imp)
   time_imp = time_imp - time_1
!
!
!  Generate MO integrals and calculate Fock matrix
!  -----------------------------------------------
   call cpu_time(time_1)
!
!  Calculate response integrals or energy integrals. Energy integrals
!  assumed on disk in a response calculation
   call mlcc3_intermediates_calc(resp_option)
!   
   call cpu_time(time_int)
   time_int = time_int - time_1
!
!
!  Figure out which omega calculator to use
!  ----------------------------------------
!
   if(mlcc3_active) then
!
      work_required = nn*n_vir_act*n_vir_gen*n_occ_act**2 & 
     &              + nn*n_vir_act**2*n_occ_gen*n_occ_act &
     &              + 2*n_vir_act**3 + n_vir_act**2 &
     &              + 9*n_vir_aag
!
   else
!
      work_required = nn*n_vir**2*n_occ**2 & 
     &              + 2*n_vir_act**3 + n_vir_act**2 &
     &              + 9*n_vir_aag
!
   end if
!
   work_avail = work_free()
!
   if(work_avail .gt. work_required) then
!
      highmem = .true.
!
   else
!
      highmem = .false.
!
   end if
!
   if(print_mlcc3 .ge. 1) then
      if(highmem) then
!
         write(lupri,*)
         write(lupri,*) 'Batching omega chosen'
         write(lupri,*)
!
      else
!
         write(lupri,*)
         write(lupri,*) 'Non-batching omega chosen'
         write(lupri,*)
!
      end if
   end if
!
!  Calculate and add contributions to Omega Vectors
!  ------------------------------------------------
   call cpu_time(time_1)
!
!
   if(.not. resp_option) then
!
!     Energy calculation
!
      resp_control = 0
!
!     Calculate standard [Ĥ,T_3] terms in Omega_1 and Omega_2
!
      if(highmem) then
!
         call mlcc3_h_omega_contrib(resp_control)
!
      else
!
         call mlcc3_omega_contrib(resp_control)
!
      end if
!
   else
!
!     Excitation energy calculation  
!  
      resp_control = 1
!
!     Calculate [Ĥ,R_3] terms in Omega_1 and Omega_2
!
      if(highmem) then
!
         call mlcc3_h_omega_contrib(resp_control)
!
      else
!
         call mlcc3_omega_contrib(resp_control)
!
      end if
!
      resp_control = 2
!
!     Calculate [[Ĥ,R_1],T_3] term in Omega_2
!
      if(highmem) then
!
         call mlcc3_h_omega_contrib(resp_control)
!
      else
!
         call mlcc3_omega_contrib(resp_control)
!
      end if
!
   end if
!
   call cpu_time(time_ome)
   time_ome = time_ome - time_1
!
!
!  Deallocate sirius information arrays
!  ------------------------------------
   call cpu_time(time_1)
!
   call sirius_deallocator
!
   call cpu_time(time_de)
   time_de = time_de - time_1
!
!
!  Total timing
   call cpu_time(time_tot)
!   
   time_tot = time_tot - time_start
!      
   call system_clock(w2)
!
   wall_time = (w2-w1)/real(wr)
!
!  Print timings
!  -------------
!
   if (print_mlcc3 .ge. 3) then
!
      write(lupri,*)
      write(lupri,*)
      write(lupri,*) 'Timings from mlcc3_drv'
      write(lupri,9999) 'Total', time_tot
      write(lupri,9999) 'work_setup', time_ws
      write(lupri,9999) 'sirius_reader', time_read
      write(lupri,9999) 'omega_setup', time_os
      write(lupri,9999) 'mlcc3_lambda_ao', time_lambda
      write(lupri,9999) 'center_import', time_imp
      write(lupri,9999) 'mlcc3_intermediates', time_int
      write(lupri,9999) 'mlcc3_omega_contrib', time_ome
      write(lupri,9999) 'deallocation', time_de
      write(lupri,9999) 'wall time', wall_time
      write(lupri,*)
      write(lupri,*)
!
   end if
!
end subroutine mlcc3_drv
!
!
!
!
subroutine mlcc3_input(word,lunit)
!
use mlcc_typedef
use mlcc_block_import
use mlcc3_data
!
!  mlcc3 input routine
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
   implicit none
!
   integer, intent(in)  :: lunit
   character(len=7)     :: word
!   
   if (set) return
   set = .true.
!
!  initialize defaults to canonical orbitals
!
   mlcc3_active = .false. !Canonical orbitals default
   mlcc3_nrg_spa = .false. !Cholesky space default
   mlcc3_nrg_gen = .false. !Free index over all orbitals
!
   n_active  = 1 ! one active space
   n_general = 2 ! two general spaces
!
   print_mlcc3 = 1 ! print level
!
   if (word(1:7) .eq. '*MLCC3 ') then
!
      do
!
!        read new input line.
!
         read(lunit,'(a7)') word
         write(lupri,*) 'read keyword ',word
         call flshfo(lupri)
!      
         do while (word(1:1) .eq. '!' .or. word(1:1) .eq. '#' )
            read (lunit,'(a7)') word
         end do
!      
         if (word(1:1) .eq. '.') then
!
!            if (word(1:7) .eq. '.ML3ACT') then
!               
!               mlcc3_active = .true.
!               cycle
!               
!            else if(word(1:7) .eq. '.ACTSPA') then
!
!               read(lunit,*) n_active
!               read(lunit,*) n_general
!               
!               cycle
!
!            else if(word(1:7) .eq. '.ACTNRG') then
!
!               mlcc3_nrg_spa = .true.
!
!               read(lunit,*) n_occ_inp
!               read(lunit,*) n_vir_inp
!               
!               cycle
!
!            else if(word(1:7) .eq. '.ACTGEN') then
!
!               mlcc3_nrg_gen = .true.
!
!               read(lunit,*) n_gen_o_inp
!               read(lunit,*) n_gen_v_inp
!               
!               cycle
!!
!            else if(word(1:7) .eq. '.PRINT ') then
            if(word(1:7) .eq. '.PRINT ') then
!
               read(lunit,*) print_mlcc3
!               
               cycle
!
            else
!
               write (lupri,'(/5a,/)') 'prompt "',word, &
            &     '" not recognized in mlcc3_input.'
!         
               call quit('illegal prompt in mlcc3_input')
!
            end if
!            
         else if (word(1:1) .ne. '*') then
!         
            write (lupri,'(/5a,/)') 'prompt "',word, &
            '" not recognized in mlcc3_input.'
            call quit('illegal prompt in mlcc3_input')
!         
         else if (word(1:1) .eq.'*') then
!         
            backspace (lunit)
            exit
!         
         end if
!
      end do
!
   end if
!
!  Sanity checks
!
   if(mlcc3_nrg_spa .and. .not. mlcc3_active) then
      call quit('mlcc3_nrg_spa without mlcc3_active')
   end if
!
   if(mlcc3_nrg_gen .and. .not. mlcc3_nrg_spa) then
      call quit('mlcc3_nrg_gen without mlcc3_nrg_spa')
   end if
!
   if(n_active .gt. n_general) then
      call quit('n_active greater than n_general in mlccsdpt_input')
   end if
!
   if(mlcc3_nrg_gen .and. (n_gen_o_inp .gt. n_occ_inp)) then
      call quit('n_gen_o_inp greater than n_occ_inp in mlccsdpt_input')
   end if
!
   if(mlcc3_nrg_gen .and. (n_gen_v_inp .gt. n_vir_inp)) then
      call quit('n_gen_v_inp greater than n_vir_inp in mlccsdpt_input')
   end if
!
   if(print_mlcc3 .ge. 3) then
!
      write(lupri,*)
      write(lupri,*) 'Output from mlcc3_input'
      write(lupri,*) 'mlcc3_active: ', mlcc3_active
      write(lupri,*) 'n_active:     ', n_active
      write(lupri,*) 'n_general:    ', n_general
      write(lupri,*) 'print_mlcc3:  ', print_mlcc3
      write(lupri,*)
!
   end if
!
end subroutine mlcc3_input
