subroutine mlccsdpt_drv(out_energy,work,int_work,lwork)
!
!
!  mlccsdpt driver
!  Authors Henrik Koch and Rolf H. Myhre
!  May 2015
!
!  Purpose: Calculate multi-level CCSD(T) energy contribution
!
   use mlcc_work
   use mlcc_typedef
   use mlcc_block_import
   use mlcc3_init
   use mlcc3_active_spaces
   use mlccsdpt_integrals
   use mlccsdpt_energy
!
   implicit none
!
   real(dp), intent(out)                  :: out_energy
!
   integer, intent(in)                    :: lwork
   real(dp), dimension(lwork), intent(in) :: work
   integer, dimension(lwork), intent(in)  :: int_work
!
   real(dp) :: energy_pt
!
!
!  Timing variables
   real     :: time_tot, time_ws, time_read, time_imp, time_X2
   real     :: time_start, time_1, time_int, time_energy, time_de
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
   energy_pt = zero
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
!  Read in X amplitudes
!  --------------------
!
   call cpu_time(time_1)
!
   call amplitude_setup
!
   call cpu_time(time_X2)
   time_X2 = time_X2 - time_1
!
!
!  Import active space information from cc_active
!  ----------------------------------------------
   call cpu_time(time_1)
!
   call center_h_import
!   
   call cpu_time(time_imp)
   time_imp = time_imp - time_1
!
!
!  Generate MO integrals
!  ---------------------
   call cpu_time(time_1)
!
   call mlccsdpt_int_calc
!   
   call cpu_time(time_int)
   time_int = time_int - time_1
!
!
!  Calculate MLCCSD(T) energy contribution
!  ---------------------------------------
!
   call cpu_time(time_1)
!
   call mlccsdpt_e_calc(energy_pt)
!   
   out_energy = energy_pt
!   
   call cpu_time(time_energy)
   time_energy = time_energy - time_1
!
!
!  Deallocate sirius information arrays
!  ------------------------------------
   call cpu_time(time_1)
!
   call sirius_pt_deallocator
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
!
!  Print timings
!  -------------
!
   if (print_mlcc3 .ge. 3) then
!
      write(lupri,*)
      write(lupri,*)
      write(lupri,*) 'Timings from mlccsdpt_drv'
      write(lupri,9999) 'Total', time_tot
      write(lupri,9999) 'work_setup', time_ws
      write(lupri,9999) 'sirius_reader', time_read
      write(lupri,9999) 'amplitude_setup', time_X2
      write(lupri,9999) 'center_import', time_imp
      write(lupri,9999) 'mlccsdpt_integrals', time_int
      write(lupri,9999) 'mlccsdpt_energy', time_energy
      write(lupri,9999) 'deallocation', time_de
      write(lupri,9999) 'wall time', wall_time
      write(lupri,*)
      write(lupri,*)
!
   end if
!
end subroutine mlccsdpt_drv
!
!
subroutine mlccsdpt_input(word,lunit)
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
   mlcc3_active = .false. !Full spce default
   mlcc3_nrg_spa = .false. !Cholesky space default
   mlcc3_nrg_gen = .false. !Free index over all orbitals
!
   n_active  = 1 ! one active space
   n_general = 2 ! two general spaces
!
   print_mlcc3 = 1 ! print level
!
   if (word(1:7) .eq. '*MLCCPT') then
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
!!               
!               mlcc3_active = .true.
!               cycle
!!               
!            else if(word(1:7) .eq. '.ACTSPA') then
!!
!               read(lunit,*) n_active
!               read(lunit,*) n_general
!!               
!               cycle
!!
!            else if(word(1:7) .eq. '.ACTNRG') then
!!
!               mlcc3_nrg_spa = .true.
!!
!               read(lunit,*) n_occ_inp
!               read(lunit,*) n_vir_inp
!!               
!               cycle
!!
!            else if(word(1:7) .eq. '.ACTGEN') then
!!
!               mlcc3_nrg_gen = .true.
!!
!               read(lunit,*) n_gen_o_inp
!               read(lunit,*) n_gen_v_inp
!!               
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
            &     '" not recognized in mlccsdpt_input.'
!         
               call quit('illegal prompt in mlccsdpt_input')
!
            end if
!            
         else if (word(1:1) .ne. '*') then
!         
            write (lupri,'(/5a,/)') 'prompt "',word, &
            '" not recognized in mlccsdpt_input.'
            call quit('illegal prompt in mlccsdpt_input')
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
      write(lupri,*) 'Output from mlccsdpt_input'
      write(lupri,*) 'mlcc3_active: ', mlcc3_active
      write(lupri,*) 'n_active:     ', n_active
      write(lupri,*) 'n_general:    ', n_general
      write(lupri,*) 'print_mlcc3:  ', print_mlcc3
      write(lupri,*)
!
   end if
!
end subroutine mlccsdpt_input
