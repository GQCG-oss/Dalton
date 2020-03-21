  subroutine read_input_sections(word, kw_section)

    use input_reader

    implicit none

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------

    select case (kw_section)

      case ('*LUCITA')
        call read_input_lucita(word, kw_section)

      case ('*ENSEMB')
        call read_input_ensemble_dft(word, kw_section)
      
      case ('*EXTPCM')
#ifdef HAS_PCMSOLVER
        call read_input_pcm(word, kw_section)
#else
        call quit('you need to recompile with ENABLE_PCMSOLVER to use *EXTPCM')
#endif

      case default
!       activate as soon as everything is merged to here
!       call quit('section '//kw_section//' not recognized')

    end select

  end subroutine

  subroutine read_input_lucita(word, kw_section)

    use gasci_input_cfg
    use input_reader

    implicit none

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    character(len=400)               :: input_line
    integer                          :: i, j
    integer                          :: ios, islash
!   ----------------------------------------------------------------------------

#include "priunit.h"
#include "inforb.h"

    call reset_available_kw_list()

    if (kw_matches(word, '.TITLE ')) then
      call kw_read(word, gasci_input_run_title)
    end if

    if (kw_matches(word, '.INIWFC')) then
      call kw_read(word, gasci_input_ini_wavef)
    end if

    if (kw_matches(word, '.CITYPE')) then
      call kw_read(word, gasci_input_ci_type)
    end if

    if (kw_matches(word, '.NROOTS')) then
      call kw_read(word, gasci_input_nr_roots)
    end if

    if (kw_matches(word, '.SYMMET')) then
      call kw_read(word, gasci_input_ptg_symmetry)
    end if

    if (kw_matches(word, '.NACTEL')) then
      call kw_read(word, gasci_input_nr_active_e)
    end if

    if (kw_matches(word, '.MULTIP')) then
      call kw_read(word, gasci_input_is_spin_multiplett)
    end if

    if (kw_matches(word, '.PRINTG')) then
      call kw_read(word, gasci_input_global_print_lvl)
    end if

    if (kw_matches(word, '.PRINTL')) then
      call kw_read(word, gasci_input_local_print_lvl)
    end if

    if (kw_matches(word, '.INACTI')) then

      gasci_input_inactive_shell_set = .true.
      read(get_file_unit(), *) (nish_gasci_input(i), i=1,nsym)

    end if

    if (kw_matches(word, '.GAS SH')) then

!     set GAS control variables
      gasci_input_init_wave_f_type   = 1
      gasci_input_nr_gas_space_set   = .true.
      gasci_input_minmax_occ_gas_set = .true.

!     read the # of GAS spaces
      call kw_read(word, gasci_input_nr_gas_spaces)

      if(gasci_input_nr_gas_spaces > 6) call quit('# gas spaces too high. max = 6.')

!     process the min max occupation / per orbital occupation in each GAS shell
      do i = 1, gasci_input_nr_gas_spaces

        read(get_file_unit(),'(a)') input_line
        call upcase(input_line)

        islash = index(input_line,'/')
        if(islash <= 1)then
          write(lupri,'(/a,i2,a,i2/a/2a)') 'ERROR for *LUCITA .GAS SHELL shell no.',                                   &
                i,'/',gasci_input_nr_gas_spaces,'- the defining input line does not contain a "/":',            &
                '- the bad line : ',input_line
          write(lupri,'(a)') '- Or the specification of the number of GAS-shells might be wrong! Check the input.'
          call quit('Input error for .GAS SHELL under *LUCITA')
        end if

!       min max occupation
        read(input_line(1:islash-1),*,iostat=ios) ngso_gasci_input(i,1), ngso_gasci_input(i,2)

        if(ios /= 0)then
          write(lupri,'(/a,i2,a,i2/a/2a)') 'ERROR for *LUCITA .GAS SHELL shell no.',                                   &
                i,'/',gasci_input_nr_gas_spaces, '- the input line does not contain correct min max electrons', &
                '- the bad line : ',input_line
          call quit('Input error for .GAS SHELL under *LUCITA')
        end if

!       GAS shell orbital occupation
        read(input_line(islash+1:),*,iostat=ios) (ngsh_gasci_input(i,j), j=1,nsym)

        if(ios /= 0)then
          write(lupri,'(/a,i2,a,i2/a,i2,a/2a)') 'ERROR for *LUCITA .GAS SHELL shell no.',                              &
                i,'/',gasci_input_nr_gas_spaces,'- the input line does not contain',nsym,' occupations',        &
                '- the bad line : ',input_line
          call quit('Input error for .GAS SHELL under *LUCITA')
        end if
      end do

!     set # of active electrons ==> max occupation in last GAS space
      gasci_input_nr_active_e = ngso_gasci_input(gasci_input_nr_gas_spaces,2)
    end if

    if (kw_matches(word, '.RAS1  ')) then

      gasci_input_ras1_set         = .true.
      gasci_input_init_wave_f_type = 2

      read(get_file_unit(), *) (nas1_gasci_input(i), i=1,nsym)
!     call kw_read(word, gasci_input_max_holes_ras1)

    end if

    if (kw_matches(word, '.RAS1 E')) then
      read(get_file_unit(), *) gasci_input_min_e_ras1, gasci_input_max_e_ras1
    end if

    if (kw_matches(word, '.RAS2  ')) then
      gasci_input_ras2_set = .true.
      read(get_file_unit(), *) (nas2_gasci_input(i), i=1,nsym)
    end if

    if (kw_matches(word, '.RAS3  ')) then

      gasci_input_ras3_set = .true.
      read(get_file_unit(), *) (nas3_gasci_input(i), i=1,nsym)
!     call kw_read(word, gasci_input_max_e_ras3)

    end if

    if (kw_matches(word, '.RAS3 E')) then
      read(get_file_unit(), *) gasci_input_min_e_ras3, gasci_input_max_e_ras3
    end if

    if (kw_matches(word, '.DENSI ')) then
      call kw_read(word, gasci_input_density_calc_lvl)
    end if

    if (kw_matches(word, '.SPDENS')) then
      call kw_read(word, gasci_input_spindensity_calc_lvl)
    end if

    if (kw_matches(word, '.RSTART')) then
      call kw_read(word, gasci_input_restart_ci)
    end if

    if (kw_matches(word, '.MXCIVE')) then
      call kw_read(word, gasci_input_max_dav_subspace_dim)
    end if

    if (kw_matches(word, '.FCIDUM')) then
      gasci_input_fci_dump = .true.
    end if

    if (kw_matches(word, '.MAXITR')) then
      call kw_read(word, gasci_input_max_nr_dav_ci_iter)
    end if

    if (kw_matches(word, '.LBLKSZ')) then
      call kw_read(word, gasci_input_max_batch_size)
    end if

    if (kw_matches(word, '.DISTRT')) then
      call kw_read(word, gasci_input_ttss_dist_strategy)
    end if

    if (kw_matches(word, '.TRUNCF')) then
      call kw_read(word, gasci_input_accepted_truncation)
    end if

    if (kw_matches(word, '.MEMFAC')) then
      call kw_read(word, gasci_input_mem_reduction_multp)
    end if

    if (kw_matches(word, '.ANALYZ')) then
      gasci_input_analyze_cvec = .true.
    end if

    if (kw_matches(word, '.TIMING')) then
      gasci_input_timing_par   = .true.
    end if

    if (kw_matches(word, '.THRESH')) then
      call kw_read(word, gasci_input_convergence_thr)
    end if

    if (kw_matches(word, '.AUXTHR')) then
      call kw_read(word, gasci_input_aux_convergence_thr)
    end if

    if (kw_matches(word, '.NATORB')) then
      gasci_input_natural_orb_occ_nr = .true.
      if(gasci_input_density_calc_lvl < 2) gasci_input_density_calc_lvl = 1
    end if

    if (kw_matches(word, '.SKIP4I')) then
      gasci_input_skip_4index_trafo = .true.
    end if
! Mickael: I add Charge Tranfer control for GASCI 
!    if (kw_matches(word, '.CT    ')) then
!      gasci_input_charge_transfer = .true.
!    end if
!    if (kw_matches(word, '.CTAB  ')) then
!      gasci_input_occlse_manual_control = .true.
!    end if
!    if (kw_matches(word, '.GAS_CT')) then
!      gasci_input_gas_conf_for_ct = .true.
!    end if
! 

    if (kw_matches(word, '.PLUS C')) then
      gasci_input_plus_combi = .true.
    end if

    call check_whether_kw_found(word, kw_section)

  end subroutine

  subroutine read_input_ensemble_dft(word, kw_section)

    use lucita_mcscf_srdftci_cfg 
    use input_reader

    implicit none

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    integer                          :: i
!   ----------------------------------------------------------------------------

    call reset_available_kw_list()

    !> set logical to enable self-consistent ensemble-DFT calculation
    if (kw_matches(word, '.SC-ENS')) then
      do_sc_ensemble_dft = .true.
    end if

    if (kw_matches(word, '.WEIGHT')) then

      call kw_read(word, nr_of_weights)
      read(get_file_unit(),*) (weights(i), i=1,nr_of_weights)
      
    end if

    call check_whether_kw_found(word, kw_section)

  end subroutine
