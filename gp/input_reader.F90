!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module input_reader

! radovan: - this module should not initialize anything
!            it is a reader
!          - when we have only modules and no common blocks the initialization
!            problem goes away since variables can be initialized at definition

!radovan: please feel free to copy/use/improve this inside DIRAC/Dalton
!         without asking

  use keyword
  use character_processing
#ifdef PRG_DIRAC
  use dirac_cfg
#endif

  implicit none

  public read_input
  public move_to_next_star

  private

#ifdef PRG_DIRAC
#include "dgroup.h"
#else
! Dalton
#include "inforb.h"
#endif

contains

#ifdef PRG_DIRAC
  subroutine set_lscale_dft_gaunt()

#include "dcbham.h"

    lscale_dft_gaunt = .true.

  end subroutine
#endif

  subroutine read_input()

!   radovan: - this routine should not include any common blocks
!              this should be done in called subroutines
!            - keep ifndef traps inside subroutines for better readability
!              of this subroutine

!   ----------------------------------------------------------------------------
    character(kw_length) :: word
    character(kw_length) :: kw_section
!   ----------------------------------------------------------------------------

!   this is to catch keywords that appear before some section starts
    kw_section = '       '

#ifdef PRG_DIRAC
    open(unit_in, file='DIRAC.INP')
    rewind unit_in

    do while (.true.)

      read(unit_in, '(a7)', end=1) word

      if (word == '       ') then
!         blank line
      else
         select case (word(1:1))
         
           case ('!', '#')
!            comment
         
           case ('*')
!            section
             kw_section = uppercase(word)
         
           case default
!            keyword
             select case (kw_section)
         
!              case ('**METHO', '**WAVE ')
!                call read_input_method(word, kw_section)
         
               case ('*DFT   ')
                 call read_input_dft(word, kw_section)
         
!              *VISUAL is for backward compatibility
!              *VISUAL used to be under **ANALYZE
               case ('**VISUA', '*VISUAL')
                 call read_input_visual(word, kw_section)
         
               case ('*OPENRS')
                 call read_input_openrsp(word, kw_section)
         
               case ('**RELCC')
                 call read_input_relcc(word, kw_section)
         
               case ('*END OF')
                 go to 1
         
!              case ('       ')
!                call quit('keyword '//word//' without section')
         
               case default
!                activate as soon as everything is merged to here
!                call quit('section '//kw_section//' not recognized')
         
             end select
         
         end select
      end if

    end do

1   close(unit_in)

#else /* PRG_DIRAC */

    open(unit_in, file='DALTON.INP')
    rewind unit_in

    do while (.true.)

      read(unit_in, '(a7)', end=1) word

      if (word == '       ') then
!         blank line
      else
         select case (word(1:1))
         
           case ('!', '#')
!            comment
         
           case ('*')
!            section
             kw_section = uppercase(word)
         
           case default
!            keyword
             select case (kw_section)
         
!              case ('**METHO', '**WAVE ')
!                call read_input_method(word, kw_section)
         
               case ('*LUCITA')
                 call read_input_lucita(word, kw_section)
         
!              case ('*OPENRS')
!                call read_input_openrsp(word, kw_section)
         
               case ('**END OF')
                 go to 1
         
!              case ('       ')
!                call quit('keyword '//word//' without section')
         
               case default
!                activate as soon as everything is merged to here
!                call quit('section '//kw_section//' not recognized')
         
             end select
         
         end select
      end if

    end do

1   close(unit_in,status='keep')
#endif

  end subroutine

  subroutine move_to_next_star(word_io)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(inout) :: word_io
!   ----------------------------------------------------------------------------
    character(kw_length)                :: word
!   ----------------------------------------------------------------------------

    do while (.true.)
      read(unit_in, '(a7)') word
      if (word(1:1) == '*') then
        word_io = word
        return
      end if
    end do

  end subroutine

#ifdef PRG_DIRAC
  subroutine read_input_dft(word, kw_section)

    use dft_cfg
    use dirac_cfg

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    character(80)                    :: line
    character(4)                     :: first_4
    integer                          :: i, n
!   ----------------------------------------------------------------------------

    if (.not. dirac_cfg_dft_calculation) return

    call reset_available_kw_list()

    if (kw_matches(word, '.SMALLD')) then
      call kw_read(word, dft_cfg_dfthr0)
    end if

    if (kw_matches(word, '.RADINT')) then
      call kw_read(word, dft_cfg_radint)
    end if

    if (kw_matches(word, '.ANGINT')) then
      call kw_read(word, dft_cfg_angint)
    end if

    if (kw_matches(word, '.NOPRUN')) then
!     turn off pruning of angular grid
      dft_cfg_no_pruning = .true.
    end if

    if (kw_matches(word, '.ANGMIN')) then
      call kw_read(word, dft_cfg_angmin)
    end if

    if (kw_matches(word, '.ATSIZE')) then
!     estimate relative atomic sizes for use in the Becke
!     partitioning scheme from atomic contributions
      dft_cfg_estimate_radii = .true.
    end if

    if (kw_matches(word, '.ERRELS')) then
!     threshold for accepted error in integrated number of electrons
      call kw_read(word, dft_cfg_accepted_error)
    end if

    if (kw_matches(word, '.EXPORT')) then
      dft_cfg_export_grid = .true.
      call kw_read(word, dft_cfg_gridfile)
    end if

    if (kw_matches(word, '.IMPORT')) then
      dft_cfg_import_grid = .true.
      call kw_read(word, dft_cfg_gridfile)
    end if

!radovan: not documented on wiki, i plan to remove this option
    if (kw_matches(word, '.MULTIG')) then
      dft_cfg_multigrid = .true.
      call kw_read(word, n)
      dft_cfg_multigrid_max = n
      allocate(dft_cfg_multigrid_file(0:n))
      allocate(dft_cfg_multigrid_file_unit(0:n))
      do i = 0, n
        call kw_read(word, dft_cfg_multigrid_file(i))
      end do
    end if

!radovan: not documented on wiki
    if (kw_matches(word, '.ZIPGRI')) then
!     is true by default, keep for backward compatibimity
      dft_cfg_zipgrid = .true.
    end if

    if (kw_matches(word, '.NOZIP ')) then
      dft_cfg_zipgrid = .false.
    end if

    if (kw_matches(word, '.4CGRID')) then
!     include small component in the grid generation
!     also for 1- and 2-c dft calculations
      dft_cfg_force_4c_grid = .true.
    end if

    if (kw_matches(word, '.DEBUG ')) then
      dft_cfg_radint = 1.0d-3
      dft_cfg_angint = 10
    end if

    if (kw_matches(word, '.COARSE')) then
      dft_cfg_radint = 1.0d-11
      dft_cfg_angint = 35
    end if

    if (kw_matches(word, '.ULTRAF')) then
      dft_cfg_radint = 2.0d-15
      dft_cfg_angint = 64
    end if

    if (kw_matches(word, '.INTCHK')) then
      call kw_read(word, dft_cfg_integration_check_level)
    end if

    if (kw_matches(word, '.GRAC  ')) then
!     gradient regulated asymptotic correction
      dft_cfg_grac = .true.
      read(unit_in, '(a4)') first_4
      if      (lowercase(first_4) == 'lb94') then
        dft_cfg_asymptote_is_lb94 = .true.
      else if (lowercase(first_4) == 'lbal') then
        dft_cfg_asymptote_is_lbalpha = .true.
      else
        call quit('keyword following .GRAC not recognized')
      end if
      call kw_read(word,               &
                   dft_cfg_grac_alpha, &
                   dft_cfg_grac_beta,  &
                   dft_cfg_ac_ip,      &
                   dft_cfg_ac_threshold)
    end if

    if (kw_matches(word, '.SAOP  ')) then
!     statistical averaging of (model) orbital potentials
      dft_cfg_saop = .true.
      read(unit_in, '(a4)') first_4
      if      (lowercase(first_4) == 'lb94') then
        dft_cfg_asymptote_is_lb94 = .true.
      else if (lowercase(first_4) == 'lbal') then
        dft_cfg_asymptote_is_lbalpha = .true.
      else
        call quit('keyword following .SAOP not recognized')
      end if
    end if

    if (kw_matches(word, '.SAOP! ')) then
!     "original" SAOP as in JCP 112, 1344 (2000).
      dft_cfg_saop                    = .true.
      dft_cfg_saop_with_response_part = .true.
      dft_cfg_asymptote_is_lbalpha    = .true.
      dft_cfg_alda_hs                 = .true.
      dft_cfg_alda_ha                 = .true.
      call setaldahs(.true.) !activate alda for the h+ part
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.NOSDFT')) then
      dft_cfg_no_sdft = .true.
    end if

    if (kw_matches(word, '.COLLIN')) then
      dft_cfg_sdft_collinear = .true.
    end if

    if (kw_matches(word, '.BETASI')) then
      dft_cfg_betasigma = .true.
    end if

    if (kw_matches(word, '.ALDA  ')) then
      dft_cfg_alda_hs = .true.
      dft_cfg_alda_ha = .true.
      call setaldahs(.true.) !activate alda for the h+ part
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.XALDA ')) then
!     use        (1-f)*S + VWN + f*HFx
!     instead of       S + VWN
      dft_cfg_xalda   = .true.
      dft_cfg_alda_hs = .true.
      dft_cfg_alda_ha = .true.
      call setxalda(.true.)
      call setaldahs(.true.) !activate alda for the h+ part
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.ALDA+ ')) then
!     use alda only for the      hermitian part
      dft_cfg_alda_hs = .true.
      call setaldahs(.true.) !activate alda for the h+ part
    end if

    if (kw_matches(word, '.ALDA- ')) then
!     use alda only for the anti-hermitian part
      dft_cfg_alda_ha = .true.
      call setaldaha(.true.) !activate alda for the h- part
    end if

    if (kw_matches(word, '.XALDA+')) then
!     use xalda only for the      hermitian part
      dft_cfg_xalda   = .true.
      dft_cfg_alda_hs = .true.
      call setxalda(.true.)
      call setaldahs(.true.) !activate alda for the h+ part
    end if

    if (kw_matches(word, '.XALDA-')) then
!     use xalda only for the anti-hermitian part
      dft_cfg_xalda   = .true.
      dft_cfg_alda_ha = .true.
      call setxalda(.true.)
      call setaldaha(.true.) !activate alda for the h- part
    end if

!radovan: not documented on wiki, code not perfectly tested
    if (kw_matches(word, '.ALDAOR')) then
      call kw_read(word, dft_cfg_alda_order)
    end if

    if (kw_matches(word, '.PRINT ')) then
      call kw_read(word, dft_cfg_print_level)
    end if

    if (kw_matches(word, '.GAUTSC')) then
      call set_lscale_dft_gaunt()
    end if

!radovan: not documented on wiki
    if (kw_matches(word, '.FROMVX')) then
      call quit('explicit xc potential (.FROMVX) not available')
    end if

!radovan: not documented on wiki, code not perfectly tested
    if (kw_matches(word, '.VRG   ')) then
      dft_cfg_cdft_vrg = .true.
    end if

!radovan: not documented on wiki, code not perfectly tested
    if (kw_matches(word, '.SCREEN')) then
      dft_cfg_screening = .true.
    end if

!radovan: not documented on wiki, code not perfectly tested
    if (kw_matches(word, '.BLOCKE')) then
      dft_cfg_blocked = .true.
    end if

!radovan: not documented on wiki, code not perfectly tested
    if (kw_matches(word, '.LAMBDA')) then
      dirac_cfg_overlap_diagnostic = .true.
    end if

    call check_whether_kw_found(word, kw_section)

  end subroutine

  subroutine read_input_visual(word, kw_section)

    use visual_cfg

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    integer                          :: i
    character(6)                     :: reply
!   ----------------------------------------------------------------------------

    dirac_cfg_visual = .true.

    call reset_available_kw_list()

    if (kw_matches(word, '.LIST  ')) then
      visual_cfg_list = .true.
      call kw_read(word, visual_cfg_nr_points_in_list)
      do i = 1, visual_cfg_nr_points_in_list
        call kw_read(word, visual_cfg_xyz_list(i, 1), &
                           visual_cfg_xyz_list(i, 2), &
                           visual_cfg_xyz_list(i, 3))
      end do
    end if

    if (kw_matches(word, '.LINE  ')) then
      visual_cfg_line = .true.
      call kw_read(word, visual_cfg_line_nr_steps)
      call kw_read(word, visual_cfg_line_from(1), &
                         visual_cfg_line_from(2), &
                         visual_cfg_line_from(3))
      call kw_read(word, visual_cfg_line_to(1), &
                         visual_cfg_line_to(2), &
                         visual_cfg_line_to(3))
    end if

    if (kw_matches(word, '.2D    ')) then
       visual_cfg_2d = .true.
       read(unit_in, *) visual_cfg_2d_p_origin(1:3)
       read(unit_in, *) visual_cfg_2d_p_right(1:3)
       read(unit_in, *) visual_cfg_2d_nr_right
       read(unit_in, *) visual_cfg_2d_p_top(1:3)
       read(unit_in, *) visual_cfg_2d_nr_top
    end if

    if (kw_matches(word, '.2D_INT')) then
       visual_cfg_2d_integration = .true.
       read(unit_in, *) visual_cfg_2d_integration_p_origin(1:3)
       read(unit_in, *) visual_cfg_2d_integration_p_right(1:3)
       read(unit_in, *) visual_cfg_2d_integration_nr_right
       read(unit_in, *) visual_cfg_2d_integration_p_top(1:3)
       read(unit_in, *) visual_cfg_2d_integration_nr_top
       read(unit_in, *) visual_cfg_2d_integration_order
    end if

    if (kw_matches(word, '.3D    ')) then
      visual_cfg_3d = .true.
    end if

    if (kw_matches(word, '.QUANTI')) then
      call kw_read(word, visual_cfg_nr_xvectors)
      do i = 1, visual_cfg_nr_xvectors
        read(unit_in, *) visual_cfg_xvector_file(i),   &
                         visual_cfg_iwhich_xvector(i), &
                         visual_cfg_iq_xvector(i),     &
                         visual_cfg_ioption_xvector(i)
      end do
    end if

    if (kw_matches(word, '.DENSIT')) then
      visual_cfg_nr_xvectors        = 1
      visual_cfg_xvector_file(1)    = 'ZORDER'
      visual_cfg_iwhich_xvector(1)  = 1
      visual_cfg_iq_xvector(1)      = 1
      visual_cfg_ioption_xvector(1) = 0
    end if

    if (kw_matches(word, '.GAMMA5')) then
      visual_cfg_nr_xvectors        = 1
      visual_cfg_xvector_file(1)    = 'ZORDER'
      visual_cfg_iwhich_xvector(1)  = 1
      visual_cfg_iq_xvector(1)      = 15
      visual_cfg_ioption_xvector(1) = 0
    end if

    if (kw_matches(word, '.ELF   ')) then
      visual_cfg_nr_xvectors        = 1
      visual_cfg_xvector_file(1)    = 'ZORDER'
      visual_cfg_iwhich_xvector(1)  = 1
      visual_cfg_iq_xvector(1)      = 12
      visual_cfg_ioption_xvector(1) = 0
    end if

    if (kw_matches(word, '.INTEGR')) then
      visual_cfg_integrate = .true.
    end if

    if (kw_matches(word, '.RADIAL')) then
      visual_cfg_radial = .true.
    end if

    if (kw_matches(word, '.SMALLA')) then
      visual_cfg_force_small_ao = .true.
    end if

    if (kw_matches(word, '.ORBSTR')) then
      visual_cfg_use_orbital_string = .true.
      call kw_read(word, visual_cfg_orbital_string(1))
      if (nfsym > 1) then
        call kw_read(word, visual_cfg_orbital_string(2))
      end if
    end if

    if (kw_matches(word, '.COMMON')) then
      call kw_read(word, visual_cfg_common_factor)
    end if

    if (kw_matches(word, '.3DADD ')) then
      call kw_read(word, visual_cfg_add_3d)
    end if

    if (kw_matches(word, '.CUBE  ')) then
      call kw_read(word, visual_cfg_ncube(1), &
                         visual_cfg_ncube(2), &
                         visual_cfg_ncube(3))
    end if

    if (kw_matches(word, '.GAUGE ')) then
      call kw_read(word, visual_cfg_gauge_origin(1), &
                         visual_cfg_gauge_origin(2), &
                         visual_cfg_gauge_origin(3))
    end if

#ifdef MOD_UNRELEASED
    if (kw_matches(word, '.LONDON')) then
      visual_cfg_london = .true.
      do i = 1, visual_cfg_nr_xvectors
        read(unit_in, *) visual_cfg_conmat_file(i),      &
                         visual_cfg_lao_chcomp(i),       &
                         visual_cfg_connection_key(i),   &
                         visual_cfg_lao_contributions(i)
!     gosia todo: the last one should be made optional
      end do
    end if
#endif

    call check_whether_kw_found(word, kw_section)

  end subroutine

  subroutine read_input_openrsp(word, kw_section)

#ifdef MOD_OPENRSP
    use openrsp_cfg
#endif

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    integer                          :: i
!   ----------------------------------------------------------------------------

#ifndef MOD_OPENRSP
    call quit(kw_section//' not available in this version')
#else

    call reset_available_kw_list()

    if (kw_matches(word, '.THRESH')) then
      call kw_read(word, openrsp_cfg_threshold_response, openrsp_cfg_threshold_rhs)
    end if

    if (kw_matches(word, '.FREQ  ')) then
       call kw_read(word, openrsp_cfg_nr_freq)
       allocate(openrsp_cfg_freq(openrsp_cfg_nr_freq))
       do i = 1, openrsp_cfg_nr_freq
          call kw_read(word, openrsp_cfg_freq(i))
       end do
    end if

    if (kw_matches(word, '.FREQ_P')) then
!     give frequency for a process
      call kw_read(word, openrsp_cfg_freq_process)
    end if

    if (kw_matches(word, '.SKIPEE')) then
      openrsp_cfg_skip_pp = .true.
    end if

    if (kw_matches(word, '.SKIPEP')) then
      openrsp_cfg_skip_pn = .true.
    end if

    if (kw_matches(word, '.SKIP1E')) then
      openrsp_cfg_skip_1el = .true.
    end if

    if (kw_matches(word, '.ALPHA ')) then
!     frequency-dependent \alpha
      openrsp_cfg_alpha = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.0ALPHA')) then
!     static \alpha
      openrsp_cfg_alpha_lin_static = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.BETA  ')) then
!     frequency-dependent \beta by the n+1 rule
      openrsp_cfg_beta = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.BETA2N')) then
!     frequency-dependent \beta by the 2n+1 rule
      openrsp_cfg_beta_2np1 = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.0BETA ')) then
!     static \beta
      openrsp_cfg_beta_lin_static = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.GAMMA ')) then
!     frequency-dependent \gamma
      openrsp_cfg_gamma = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.0GAMMA')) then
!     static \gamma
      openrsp_cfg_gamma_lin_static = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.DC-SHG')) then
      openrsp_cfg_gamma_lin_dc_shg = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.THG   ')) then
      openrsp_cfg_gamma_lin_thg = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.IDRI  ')) then
      openrsp_cfg_gamma_lin_idri = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.DC-KER')) then
      openrsp_cfg_gamma_lin_dc_kerr = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.MAGNET')) then
      openrsp_cfg_magnetizability = .true.
      need_1el_o = .true.
      need_1el_b = .true.
    end if

    if (kw_matches(word, '.QUADRU')) then
      openrsp_cfg_quadrupole = .true.
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.ROA   ')) then
      openrsp_cfg_roa = .true.
      need_1el_f = .true.
      need_1el_q = .true.
      need_1el_o = .true.
      need_1el_g = .true.
    end if

    if (kw_matches(word, '.CARS  ')) then
      openrsp_cfg_cars = .true.
      need_1el_f = .true.
      need_1el_q = .true.
      need_1el_o = .true.
      need_1el_g = .true.
    end if

    if (kw_matches(word, '.PVBETA')) then
       openrsp_cfg_pv_beta = .true.
       need_1el_f = .true.
       need_1el_g = .true.
    end if

    if (kw_matches(word, '.PVGAMM')) then
       openrsp_cfg_pv_gamma = .true.
       need_1el_f = .true.
       need_1el_g = .true.
    end if

    if (kw_matches(word, '.NOLAO ')) then
      openrsp_cfg_nolao = .true.
    end if

    if (kw_matches(word, '.EFGB  ')) then
      openrsp_cfg_efgb = .true.
      need_1el_f = .true.
      need_1el_q = .true.
      need_1el_v = .true.
      need_1el_o = .true.
      need_1el_b = .true.
    end if

    if (kw_matches(word, '.JONES ')) then
      openrsp_cfg_jones = .true.
      need_1el_f = .true.
      need_1el_q = .true.
      need_1el_o = .true.
      need_1el_b = .true.
    end if

    if (kw_matches(word, '.CME   ')) then
      openrsp_cfg_cme = .true.
      need_1el_f = .true.
      need_1el_q = .true.
      need_1el_o = .true.
      need_1el_b = .true.
    end if

    if (kw_matches(word, '.MOLGRD')) then
      openrsp_cfg_dg_energy = .true.
      need_1el_f = .true.
      need_1el_g = .true.
    end if

    if (kw_matches(word, '.INTEGR')) then
      openrsp_cfg_integrate = .true.
    end if

    if (kw_matches(word, '.VISUAL')) then
      openrsp_cfg_visual = .true.
    end if

    if (kw_matches(word, '.ZZZ   ')) then
!     calculate static derivatives E_{zzz...}
!     of given order
      openrsp_cfg_zzz = .true.
      call kw_read(word, openrsp_cfg_zzz_order)
      need_1el_f = .true.
    end if

    if (kw_matches(word, '.DELTA ')) then
      openrsp_cfg_test_delta = .true.
      need_1el_f = .true.
    end if

    call check_whether_kw_found(word, kw_section)
#endif /* ifndef MOD_OPENRSP */

  end subroutine

  subroutine read_input_relcc(word, kw_section)

    use relcc_cfg

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    character(80)                    :: line
    character(4)                     :: first_4
!   ----------------------------------------------------------------------------

    call reset_available_kw_list()

    if (kw_matches(word, '.DEBUG ')) then
      relcc_debug = .true.
    end if

    if (kw_matches(word, '.TIMING')) then
      relcc_timing = .true.
    end if

    if (kw_matches(word, '.NOSORT')) then
      relcc_do_sort = .false.
    end if

    if (kw_matches(word, '.ENERGY')) then
      relcc_do_energy = .true.
    end if

    if (kw_matches(word, '.GRADIE')) then
      relcc_do_gradient = .true.
    end if

    if (kw_matches(word, '.FOCKSP')) then
      relcc_do_fspc     = .true.
    end if

    if (kw_matches(word, '.PRINT ')) then
      call kw_read(word, relcc_print)
    end if

    call check_whether_kw_found(word, kw_section)

  end subroutine

#else /* ifdef PRG_DIRAC */

  subroutine read_input_lucita(word, kw_section)

    use lucita_cfg

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    integer                          :: i, j
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
      allocate(nish_lucita(max_number_of_ptg_irreps))
      read(unit_in, *) (nish_lucita(i), i=1,nsym)

    end if

    if (kw_matches(word, '.GASSHE')) then

      lucita_cfg_init_wave_f_type = 1
      lucita_cfg_nr_gas_space_set = .true.

      call kw_read(word, lucita_cfg_nr_gas_spaces)

      if(lucita_cfg_nr_gas_spaces <= max_number_of_gas_spaces)then

        allocate(ngsh_lucita(max_number_of_gas_spaces,max_number_of_ptg_irreps))
        do i = 1, lucita_cfg_nr_gas_spaces
          read(unit_in, *) (ngsh_lucita(i,j), j=1,nsym)
        end do

      else
        call quit('# gas spaces too high. max = 6.')
      end if

    end if

    if (kw_matches(word, '.GASSPC')) then

      lucita_cfg_minmax_occ_gas_set = .true.
      call kw_read(word, lucita_cfg_nr_calc_sequences)

      if(.not.lucita_cfg_nr_gas_space_set)then 
        call quit(' error in input reading: .GASSHE has to be specified before .GASSPC.')
      else
       
        allocate(ngso_lucita(max_number_of_gas_spaces,2))
        do i = 1, lucita_cfg_nr_gas_spaces
          read(unit_in, *) (ngso_lucita(i,j), j=1,2)
        end do
      end if

    end if

    if (kw_matches(word, '.RAS1  ')) then

      lucita_cfg_ras1_set         = .true.
      lucita_cfg_init_wave_f_type = 2

      allocate(nas1_lucita(max_number_of_ptg_irreps))
      read(unit_in, *) (nas1_lucita(i), i=1,nsym)
      call kw_read(word, lucita_cfg_max_holes_ras1)

    end if

    if (kw_matches(word, '.RAS2  ')) then
      lucita_cfg_ras2_set = .true.
      allocate(nas2_lucita(max_number_of_ptg_irreps))
      read(unit_in, *) (nas2_lucita(i), i=1,nsym)
    end if

    if (kw_matches(word, '.RAS3  ')) then

      lucita_cfg_ras3_set = .true.
      allocate(nas3_lucita(max_number_of_ptg_irreps))
      read(unit_in, *) (nas3_lucita(i), i=1,nsym)
      call kw_read(word, lucita_cfg_max_e_ras3)

    end if

    if (kw_matches(word, '.DENSI ')) then
      call kw_read(word, lucita_cfg_density_calc_lvl)
    end if

    if (kw_matches(word, '.RSTRCI')) then
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

    call check_whether_kw_found(word, kw_section)

  end subroutine

#endif
 end module
