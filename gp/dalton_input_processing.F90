!dalton_copyright_start
!
!
!dalton_copyright_end

module dalton_input_processing

! stefan: - this module reads the input section of a specific module
!           within Dalton.
!
!radovan: please feel free to copy/use/improve this inside DIRAC/Dalton
!         without asking

  use keyword
  use character_processing

  implicit none

  public read_dalton_input

#include "inforb.h"
#include "priunit.h"

contains

  subroutine read_dalton_input(word, kw_section)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------

    select case (kw_section)
         
      case ('*LUCITA')
        call read_input_lucita(word, kw_section)
         
      case default
!       activate as soon as everything is merged to here
!       call quit('section '//kw_section//' not recognized')
         
    end select

  end subroutine read_dalton_input

  subroutine read_input_lucita(word, kw_section)

    use lucita_cfg

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    character(len=400)               :: input_line
    integer                          :: i, j
    integer                          :: ios, islash
!   ----------------------------------------------------------------------------

    call reset_available_kw_list()

    if (kw_matches(word, '.TITLE ')) then
      call kw_read(word, lucita_cfg_run_title)
    end if

    if (kw_matches(word, '.INIWFC')) then
      call kw_read(word, lucita_cfg_ini_wavef)
    end if

    if (kw_matches(word, '.CITYPE')) then
      call kw_read(word, lucita_cfg_ci_type)
    end if

    if (kw_matches(word, '.NROOTS')) then
      call kw_read(word, lucita_cfg_nr_roots)
    end if

    if (kw_matches(word, '.SYMMET')) then
      call kw_read(word, lucita_cfg_ptg_symmetry)
    end if

    if (kw_matches(word, '.NACTEL')) then
      call kw_read(word, lucita_cfg_nr_active_e)
    end if

    if (kw_matches(word, '.MULTIP')) then
      call kw_read(word, lucita_cfg_is_spin_multiplett)
    end if

    if (kw_matches(word, '.PRINTG')) then
      call kw_read(word, lucita_cfg_global_print_lvl)
    end if

    if (kw_matches(word, '.PRINTL')) then
      call kw_read(word, lucita_cfg_local_print_lvl)
    end if

    if (kw_matches(word, '.INACTI')) then

      lucita_cfg_inactive_shell_set = .true.
      read(unit_in, *) (nish_lucita(i), i=1,nsym)

    end if

    if (kw_matches(word, '.GAS SH')) then

!     set GAS control variables
      lucita_cfg_init_wave_f_type   = 1
      lucita_cfg_nr_gas_space_set   = .true.
      lucita_cfg_minmax_occ_gas_set = .true.

!     read the # of GAS spaces
      call kw_read(word, lucita_cfg_nr_gas_spaces)

      if(lucita_cfg_nr_gas_spaces > 6) call quit('# gas spaces too high. max = 6.')

!     process the min max occupation / per orbital occupation in each GAS shell
      do i = 1, lucita_cfg_nr_gas_spaces

        read(unit_in,'(a)') input_line
        call upcase(input_line)
 
        islash = index(input_line,'/')
        if(islash <= 1)then
          write(lupri,'(/a,i2,a,i2/a/2a)') 'ERROR for *LUCITA .GAS SHELL shell no.',                           &
                i,'/',lucita_cfg_nr_gas_spaces,'- the defining input line does not contain a "/":',            &
                '- the bad line : ',input_line
          write(lupri,'(a)') '- Or the specification of the number of GAS-shells might be wrong! Check the input.'
          call quit('Input error for .GAS SHELL under *LUCITA')
        end if

!       min max occupation
        read(input_line(1:islash-1),*,iostat=ios) ngso_lucita(i,1), ngso_lucita(i,2)

        if(ios /= 0)then
          write(lupri,'(/a,i2,a,i2/a/2a)') 'ERROR for *LUCITA .GAS SHELL shell no.',                           &
                i,'/',lucita_cfg_nr_gas_spaces, '- the input line does not contain correct min max electrons', &
                '- the bad line : ',input_line
          call quit('Input error for .GAS SHELL under *LUCITA')
        end if

!       GAS shell orbital occupation
        read(input_line(islash+1:),*,iostat=ios) (ngsh_lucita(i,j), j=1,nsym)

        if(ios /= 0)then
          write(lupri,'(/a,i2,a,i2/a,i2,a/2a)') 'ERROR for *LUCITA .GAS SHELL shell no.',                      &
                i,'/',lucita_cfg_nr_gas_spaces,'- the input line does not contain',nsym,' occupations',        &
                '- the bad line : ',input_line
          call quit('Input error for .GAS SHELL under *LUCITA')
        end if
      end do

!     set # of active electrons ==> max occupation in last GAS space
      lucita_cfg_nr_active_e = ngso_lucita(lucita_cfg_nr_gas_spaces,2)
    end if

    if (kw_matches(word, '.RAS1  ')) then

      lucita_cfg_ras1_set         = .true.
      lucita_cfg_init_wave_f_type = 2

      read(unit_in, *) (nas1_lucita(i), i=1,nsym)
      call kw_read(word, lucita_cfg_max_holes_ras1)

    end if

    if (kw_matches(word, '.RAS2  ')) then
      lucita_cfg_ras2_set = .true.
      read(unit_in, *) (nas2_lucita(i), i=1,nsym)
    end if

    if (kw_matches(word, '.RAS3  ')) then

      lucita_cfg_ras3_set = .true.
      read(unit_in, *) (nas3_lucita(i), i=1,nsym)
      call kw_read(word, lucita_cfg_max_e_ras3)

    end if

    if (kw_matches(word, '.DENSI ')) then
      call kw_read(word, lucita_cfg_density_calc_lvl)
    end if

    if (kw_matches(word, '.RSTART')) then
      call kw_read(word, lucita_cfg_restart_ci)
    end if

    if (kw_matches(word, '.MXCIVE')) then
      call kw_read(word, lucita_cfg_max_dav_subspace_dim)
    end if

    if (kw_matches(word, '.MAXITR')) then
      call kw_read(word, lucita_cfg_max_nr_dav_ci_iter)
    end if

    if (kw_matches(word, '.LBLKSZ')) then
      call kw_read(word, lucita_cfg_max_batch_size)
    end if

    if (kw_matches(word, '.DISTRT')) then
      call kw_read(word, lucipar_cfg_ttss_dist_strategy)
    end if

    if (kw_matches(word, '.TRUNCF')) then
      call kw_read(word, lucita_cfg_accepted_truncation)
    end if

    if (kw_matches(word, '.MEMFAC')) then
      call kw_read(word, lucipar_cfg_mem_reduction_multp)
    end if

    if (kw_matches(word, '.ANALYZ')) then
      lucita_cfg_analyze_cvec = .true.
    end if

    if (kw_matches(word, '.TIMING')) then
      lucita_cfg_timing_par   = .true.
    end if

    if (kw_matches(word, '.NATORB')) then
      lucita_cfg_natural_orb_occ_nr = .true.
      if(lucita_cfg_density_calc_lvl < 2) lucita_cfg_density_calc_lvl = 1
    end if

    call check_whether_kw_found(word, kw_section)

  end subroutine

end module
