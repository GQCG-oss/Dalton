!dalton_copyright_start
!
!
!dalton_copyright_end

module mcscf_or_gasci_2_define_cfg

! stefan: - this module provides an interface between the MCSCF input / GASCI input 
!           in order to define the generic LUCITA cfg variables / LUCITA orbital spaces.
!
!           output: variables in lucita_cfg
!
!           ---------------------------------------
!           called for MCSCF from SETCI  (getcix.F)
!           called for GASCI from SIRCTL (sirctl.F)
!           ---------------------------------------


! output:
  use lucita_cfg
  use lucita_mcscf_ci_cfg
  use lucita_orbital_spaces

  implicit none

  public define_lucita_cfg_static
  public define_lucita_cfg_dynamic


contains

  subroutine define_lucita_cfg_static(run_type)

! GASCI input:
  use gasci_input_cfg
! MCSCF/general input:
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"


#include "priunit.h"
!   ----------------------------------------------------------------------------
    integer             , intent(in) :: run_type
!   ----------------------------------------------------------------------------
    integer                          :: i, j
    integer                          :: ngas_opt   ! local variable: make it global for MCSCF GAS input
!   ----------------------------------------------------------------------------

!     global settings
      lucita_cfg_nr_ptg_irreps = nsym ! sync with Dalton nsym (inforb.h)
      naos_lucita(1:8)         = nbas(1:8)
      nmos_lucita(1:8)         = norb(1:8)
      lucita_cfg_core_energy   = potnuc

      select case(run_type)

        case(1) ! MCSCF

!         character block
          lucita_cfg_run_title(1:36)              = 'SIRIUS-MCSCF run supported by LUCITA'
          if(mod(mctype,2) /= 0) then
            lucita_cfg_ci_type(1:6)   = 'GASCI '
          else
            lucita_cfg_ci_type(1:6)   = 'RASCI '
          end if

!         logical block
          lucita_cfg_inactive_shell_set = .true.
          lucita_cfg_minmax_occ_gas_set = .true.
          lucita_cfg_timing_par         = .false.
          lucita_cfg_plus_combi         = flag(58) ! plus combination of degenerate start vectors

          select case(mctype) 
            case(1,3) ! CASSCF/GASSCF
              lucita_cfg_ras1_set = .false.
              lucita_cfg_ras2_set = .false.
              lucita_cfg_ras3_set = .false.
            case(2)   ! RASSCF
              lucita_cfg_ras1_set = .true.
              lucita_cfg_ras2_set = .true.
              lucita_cfg_ras3_set = .true.
          end select
        
!         integer block
          lucita_cfg_init_wave_f_type              = mod(mctype+1,2) + 1
          lucita_cfg_nr_active_e                   = nactel
          ngas_opt                                 =  1
          if(mctype == 2) ngas_opt                 = -1
          lucita_cfg_nr_gas_spaces                 = ngas_opt
!         lucita_cfg_restart_ci                    =  0         ! use automatic detection
          lucita_cfg_max_dav_subspace_dim          =  0         ! default
          lucita_cfg_max_batch_size                =  64000000  ! default
          lucipar_cfg_ttss_dist_strategy           =  2         ! default
          lucipar_cfg_mem_reduction_multp          =  3         ! default
          lucita_cfg_spindensity_calc_lvl          = spindens_lvl


          select case(mctype)
            case(1,3) ! CASSCF/GASSCF
              if(mctype == 3)then
                call quit('*** GASSCF input for MCSCF runs has not been written yet... FIXME. ***')
              else
                ngso_lucita(1,1) = lucita_cfg_nr_active_e
                ngso_lucita(1,2) = lucita_cfg_nr_active_e
                do i = 1, lucita_cfg_nr_gas_spaces
                  do j = 1, nsym
                    ngsh_lucita(i,j) = nas2(j)
                  end do
                end do
              end if
            case(2)   ! RASSCF
              j = 0
              do i = 1, nsym
                j = j + nas1(i)
              end do
              if(j > 0)then
                lucita_cfg_min_e_ras1       = NELMN1
                lucita_cfg_max_e_ras1       = NELMX1
              end if
!             write(lupri,*) ' min/max e ras1 ==> ',lucita_cfg_min_e_ras1,lucita_cfg_max_e_ras1
              lucita_cfg_min_e_ras3           =       NELMN3
              lucita_cfg_max_e_ras3           =       NELMX3
!             write(lupri,*) ' min e-    ras3 ==> ',lucita_cfg_min_e_ras3
!             write(lupri,*) ' max e-    ras3 ==> ',lucita_cfg_max_e_ras3
              call icopy(max_number_of_ptg_irreps,nas1,1,nas1_lucita,1)
              call icopy(max_number_of_ptg_irreps,nas2,1,nas2_lucita,1)
              call icopy(max_number_of_ptg_irreps,nas3,1,nas3_lucita,1)
          end select
          
          call icopy(max_number_of_ptg_irreps,nish,1,nish_lucita,1)

        case(2) ! GASCI

!         character block
          lucita_cfg_run_title(1:72)      = gasci_input_run_title(1:72)
          lucita_cfg_ci_type(1:72)        = gasci_input_ci_type(1:72)

!         logical block
          lucita_cfg_inactive_shell_set   = gasci_input_inactive_shell_set
          lucita_cfg_minmax_occ_gas_set   = gasci_input_minmax_occ_gas_set
          lucita_cfg_ras1_set             = gasci_input_ras1_set
          lucita_cfg_ras2_set             = gasci_input_ras2_set
          lucita_cfg_ras3_set             = gasci_input_ras3_set
          lucita_cfg_timing_par           = gasci_input_timing_par
          lucita_cfg_fci_dump             = gasci_input_fci_dump
          lucita_cfg_plus_combi           = gasci_input_plus_combi! plus combination of degenerate start vectors

!         integer block
          lucita_cfg_init_wave_f_type     = gasci_input_init_wave_f_type
          lucita_cfg_nr_active_e          = gasci_input_nr_active_e
          lucita_cfg_nr_gas_spaces        = gasci_input_nr_gas_spaces
          lucita_cfg_min_e_ras1           = gasci_input_min_e_ras1
          lucita_cfg_max_e_ras1           = gasci_input_max_e_ras1
          lucita_cfg_min_e_ras3           = gasci_input_min_e_ras3
          lucita_cfg_max_e_ras3           = gasci_input_max_e_ras3
          lucita_cfg_restart_ci           = gasci_input_restart_ci
          lucita_cfg_max_dav_subspace_dim = gasci_input_max_dav_subspace_dim
          lucita_cfg_max_batch_size       = gasci_input_max_batch_size
          lucipar_cfg_ttss_dist_strategy  = gasci_input_ttss_dist_strategy
          lucipar_cfg_mem_reduction_multp = gasci_input_mem_reduction_multp
          lucita_cfg_spindensity_calc_lvl = gasci_input_spindensity_calc_lvl

          call icopy(max_number_of_gas_spaces*2,ngso_gasci_input,1,ngso_lucita,1)
          call icopy(max_number_of_gas_spaces*max_number_of_ptg_irreps,ngsh_gasci_input,1,ngsh_lucita,1)
          call icopy(max_number_of_ptg_irreps,nas1_gasci_input,1,nas1_lucita,1)
          call icopy(max_number_of_ptg_irreps,nas2_gasci_input,1,nas2_lucita,1)
          call icopy(max_number_of_ptg_irreps,nas3_gasci_input,1,nas3_lucita,1)
          call icopy(max_number_of_ptg_irreps,nish_gasci_input,1,nish_lucita,1)

        case default 
          write(lupri,*) ' define_lucita_cfg interface: unknown run type ==> ',run_type
          call quit('*** error in define_lucita_cfg:  unknown run type. ***')
      end select

!     define LUCITA orbital spaces
      call define_lucita_orb_spaces(lucita_cfg_init_wave_f_type)

  end subroutine define_lucita_cfg_static
!*******************************************************************************

  subroutine define_lucita_cfg_dynamic(c_sym,                    &
                                       hc_sym,                   &
                                       state_of_interest,        &
                                       spin1_2el_op,             &
                                       spin2_2el_op,             &
                                       nr_states,                &
                                       spin_multiplicity,        &
                                       nr_davidson_iterations,   &
                                       density_matrix_level,     &
                                       print_lvl,                &
                                       truncation_thr,           &
                                       convergence_thr,          &
                                       aux_convergence_thr,      &
                                       calc_natorb,              &
                                       perform_c_vec_analysis,   &
                                       no_2el_part,              &
                                       mcscf_provides_ints,      &
                                       two_el_dens_matrix,       &
                                       calc_transition_densm,    &
                                       xctype1,                  &
                                       xctype2,                  &
                                       vec_xc_parallel,          &
                                       io2io_vector_exchange)

#ifdef MOD_SRDFT
  use lucita_mcscf_srdftci_cfg
#endif

#include "priunit.h"
!   ----------------------------------------------------------------------------
    integer             , intent(in) :: c_sym
    integer             , intent(in) :: hc_sym
    integer             , intent(in) :: state_of_interest
    integer             , intent(in) :: spin1_2el_op
    integer             , intent(in) :: spin2_2el_op
    integer             , intent(in) :: nr_states
    integer             , intent(in) :: spin_multiplicity
    integer             , intent(in) :: nr_davidson_iterations
    integer             , intent(in) :: density_matrix_level
    integer             , intent(in) :: xctype1
    integer             , intent(in) :: xctype2
    integer             , intent(in) :: print_lvl
    real(8)             , intent(in) :: truncation_thr
    real(8)             , intent(in) :: convergence_thr
    real(8)             , intent(in) :: aux_convergence_thr
    logical             , intent(in) :: calc_natorb
    logical             , intent(in) :: perform_c_vec_analysis
    logical             , intent(in) :: no_2el_part
    logical             , intent(in) :: mcscf_provides_ints
    logical             , intent(in) :: two_el_dens_matrix
    logical             , intent(in) :: calc_transition_densm
    logical             , intent(in) :: io2io_vector_exchange
    logical             , intent(in) :: vec_xc_parallel
!   ----------------------------------------------------------------------------

!     set 'dynamic' variables
!     -----------------------

!     logical block
      lucita_cfg_natural_orb_occ_nr                = calc_natorb
      lucita_cfg_analyze_cvec                      = perform_c_vec_analysis
      lucita_cfg_transition_densm                  = calc_transition_densm
!     logical block in module lucita_mcscf_ci_cfg
      integrals_from_mcscf_env                     = mcscf_provides_ints
      io2io_vector_exchange_mc2lu_lu2mc            = io2io_vector_exchange
      vector_exchange_mc2lu_lu2mc_active           = vec_xc_parallel

!     real(8) block
      lucita_cfg_accepted_truncation               = truncation_thr
      lucita_cfg_convergence_c                     = convergence_thr
      lucita_cfg_auxiliary_conv_c                  = aux_convergence_thr

!     integer block
      lucita_cfg_csym                              = c_sym
      lucita_cfg_hcsym                             = hc_sym
      lucita_cfg_ptg_symmetry                      = lucita_cfg_csym
      lucita_cfg_icstate                           = state_of_interest
      lucita_cfg_spin1_2el_op                      = spin1_2el_op
      lucita_cfg_spin2_2el_op                      = spin2_2el_op
      lucita_cfg_el_operator_level                 = 2
      if(no_2el_part) lucita_cfg_el_operator_level = 1
      lucita_cfg_nr_roots                          = nr_states
      lucita_cfg_is_spin_multiplett                = spin_multiplicity
      lucita_cfg_global_print_lvl                  = min(print_lvl,4)
      lucita_cfg_local_print_lvl                   = print_lvl
      lucita_cfg_max_nr_dav_ci_iter                = nr_davidson_iterations
      lucita_cfg_density_calc_lvl                  = density_matrix_level
      if(lucita_cfg_natural_orb_occ_nr)            &
      lucita_cfg_density_calc_lvl                  = max(lucita_cfg_density_calc_lvl,1)
      if(two_el_dens_matrix)                       &
      lucita_cfg_density_calc_lvl                  = max(lucita_cfg_density_calc_lvl,2)

      if(lucita_cfg_spindensity_calc_lvl > lucita_cfg_density_calc_lvl)then
        lucita_cfg_spindensity_calc_lvl  = lucita_cfg_density_calc_lvl
      end if

!     set for mcscf (if not set otherwise)
      if(lucita_cfg_max_dav_subspace_dim == 0) lucita_cfg_max_dav_subspace_dim = 16 * lucita_cfg_nr_roots

!     integer/logical block in module lucita_mcscf_ci_cfg
      vector_exchange_type1 = xctype1
      vector_exchange_type2 = xctype2
      if(vector_exchange_type1 > 1) &
      vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type1) = .true.
      if(vector_exchange_type2 > 1) & 
      vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type2) = .true.

!     set restart option if MOD_SRDFT and CIsrDFT
#ifdef MOD_SRDFT
      if(srdft_ci_with_lucita .and. srdft_restart_ci == 1)then
        lucita_cfg_restart_ci = 1
      end if
#endif

  end subroutine define_lucita_cfg_dynamic

end module
