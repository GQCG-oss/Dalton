!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_orbital_spaces

! stefan: - this module provides all necessary functionality
!           to setup the orbital spaces in LUCITA mcscf/ci calculations.
!
!           Written by sknecht for DALTON, december 2010.
!
!     nish_lucita(:)  : total number of inactive (doubly occupied) shells per irrep
!     nash_lucita(:)  : total number of active shells (GAS and RAS case) per irrep
!     nas1_lucita(:)  : total number of active shells in RAS1 (RAS case only) per irrep
!     nas2_lucita(:)  : total number of active shells in RAS2 (RAS case only) per irrep
!     nas3_lucita(:)  : total number of active shells in RAS3 (RAS case only) per irrep
!     nocc_lucita(:)  : total number of occupied shells (nish_lucita + nash_lucita) per irrep
!     ngsh_lucita(:,:): total number of active shells per gas space and irrep
!     ngso_lucita(:,2): min and max # accumulated electrons per gas space

  use lucita_cfg

  implicit none

  public define_lucita_orb_spaces
  public set_inforb_lucita

  private

contains 
 
  subroutine define_lucita_orb_spaces(number_of_ptg_irreps,        &
                                      number_of_gas_spaces,        &
                                      init_input_type,             &
                                      init_wave_f_type)
!*******************************************************************************
!
!    purpose:  define the orbital spaces for LUCITA ci/mcscf calculations.
!
!*******************************************************************************
    integer, intent(in) :: number_of_ptg_irreps
    integer, intent(in) :: number_of_gas_spaces
    integer, intent(in) :: init_input_type
    integer, intent(in) :: init_wave_f_type
!-------------------------------------------------------------------------------

!     (reset and) initialize lucita orbital space arrays
      nfro_lucita = 0
      nash_lucita = 0
      nocc_lucita = 0

!     transfer information from temporary arrays provided in: 
!     case 1: lucita input 
!     case 2:  mcscf input
!     ------------------------------------------
      select case(init_input_type)

        case(1)

          call fill_lucita_orb_spaces_ci(number_of_ptg_irreps,        &
                                         number_of_gas_spaces,        &
                                         max_number_of_ptg_irreps,    &
                                         max_number_of_gas_spaces,    &
                                         nish_lucita,                 &
                                         nash_lucita,                 &
                                         nocc_lucita,                 &
                                         nfro_lucita,                 &
                                         nas1_lucita,                 &
                                         nas2_lucita,                 &
                                         nas3_lucita,                 &
                                         ngsh_lucita,                 &
                                         init_wave_f_type)

        case(2)

          call quit(' *** error in lucita_orb_spaces_init: mcscf orbital' & 
                    //' input conversion not available yet. ***')
          
!         call fill_lucita_orb_spaces_mcscf

      end select

  end subroutine define_lucita_orb_spaces
!*******************************************************************************

  subroutine fill_lucita_orb_spaces_ci(number_of_ptg_irreps,        &
                                       number_of_gas_spaces,        &
                                       mx_number_of_ptg_irreps,     &
                                       mx_number_of_gas_spaces,     &
                                       is_nish_lucita,              &
                                       is_nash_lucita,              &
                                       is_nocc_lucita,              &
                                       is_nfro_lucita,              &
                                       is_nas1_lucita,              &
                                       is_nas2_lucita,              &
                                       is_nas3_lucita,              &
                                       is_ngsh_lucita,              &
                                       init_wave_f_type)
!*******************************************************************************
!
!    purpose:  fill orbital spaces for LUCITA from ci input.
!
!*******************************************************************************
    integer, intent(in)    :: number_of_ptg_irreps
    integer, intent(in)    :: number_of_gas_spaces
    integer, intent(in)    :: mx_number_of_ptg_irreps
    integer, intent(in)    :: mx_number_of_gas_spaces
    integer, intent(in)    :: init_wave_f_type
    integer, intent(in)    :: is_nish_lucita(mx_number_of_ptg_irreps)
    integer, intent(in)    :: is_nas1_lucita(mx_number_of_ptg_irreps)
    integer, intent(in)    :: is_nas2_lucita(mx_number_of_ptg_irreps)
    integer, intent(in)    :: is_nas3_lucita(mx_number_of_ptg_irreps)
    integer, intent(out)   :: is_nash_lucita(mx_number_of_ptg_irreps)
    integer, intent(out)   :: is_nocc_lucita(mx_number_of_ptg_irreps)
    integer, intent(out)   :: is_nfro_lucita(mx_number_of_ptg_irreps)
    integer, intent(inout) :: is_ngsh_lucita(mx_number_of_gas_spaces,mx_number_of_ptg_irreps)
!-------------------------------------------------------------------------------
    integer             :: tmp_nr_gas_spaces
    integer             :: i
    integer             :: j
!-------------------------------------------------------------------------------

    select case(init_wave_f_type)

      case(1) ! GAS
       
!       a. nothing else to do than setting the # of GA spaces
        tmp_nr_gas_spaces = number_of_gas_spaces

      case(2) ! RAS

        tmp_nr_gas_spaces = 3

!       a. active shells for each point group irrep in RAS1, RAS2 and RAS3
        do i = 1, tmp_nr_gas_spaces
          select case(i)
            case(1) 
              do j = 1, number_of_ptg_irreps
                is_ngsh_lucita(i,j) = is_nas1_lucita(j)
              end do
            case(2) 
              do j = 1, number_of_ptg_irreps
                is_ngsh_lucita(i,j) = is_nas2_lucita(j)
              end do
            case(3) 
              do j = 1, number_of_ptg_irreps
                is_ngsh_lucita(i,j) = is_nas3_lucita(j)
              end do
          end select
        end do
 
    end select

!   b. total number of active shells for each point group irrep
    do i = 1, tmp_nr_gas_spaces
      do j = 1, number_of_ptg_irreps
        is_nash_lucita(j) = is_nash_lucita(j) + is_ngsh_lucita(i,j)
      end do
    end do

!   c. total number of occupied shells for each point group irrep
    do i = 1, number_of_ptg_irreps
      is_nocc_lucita(i) = is_nash_lucita(i) + is_nish_lucita(i)
    end do

!   d. total number of frozen shells for each point group irrep
    do i = 1, number_of_ptg_irreps
      is_nfro_lucita(i) = 0
    end do

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
    print *, ' debug print of orbital spaces:'
    print *, ' nfro_lucita        : ',(is_nfro_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nish_lucita        : ',(is_nish_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nash_lucita        : ',(is_nash_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nocc_lucita        : ',(is_nocc_lucita(i),i=1,number_of_ptg_irreps)
    do j = 1, tmp_nr_gas_spaces
      print *, ' ngsh_lucita per gas: ',(is_ngsh_lucita(j,i),i=1,number_of_ptg_irreps)
    end do
#endif
!#undef LUCI_DEBUG

  end subroutine fill_lucita_orb_spaces_ci
!*******************************************************************************
 
  subroutine set_inforb_lucita(keyword)
!*******************************************************************************
!
!    purpose: transfer LUCITA orbital information to common blocks in inforb.h
!             (inforb.h defines the current orbital spaces in SIRIUS, including
!              the integral transformation module)
!
!*******************************************************************************
!-------------------------------------------------------------------------------
      character, intent(in) :: keyword*(*)
!
#include "priunit.h"
#include "inforb.h"
!
      integer, save         :: NFRO_SAVE(8), NISH_SAVE(8), NASH_SAVE(8)
!*******************************************************************************

      select case(keyword)

        case('LUCITA')

          NFRO_SAVE(1:8) = NFRO(1:8)
          NISH_SAVE(1:8) = NISH(1:8)
          NASH_SAVE(1:8) = NASH(1:8)
          NFRO(1:8) = 0
          NISH(1:8) = nish_lucita(1:8)
          NASH(1:8) = nash_lucita(1:8)

          CALL SETORB ! set derived orbital spaces in inforb.h
          
        case('RESET')

          NFRO(1:8) = NFRO_SAVE(1:8)
          NISH(1:8) = NISH_SAVE(1:8)
          NASH(1:8) = NASH_SAVE(1:8)

          CALL SETORB ! reset derived orbital spaces in inforb.h

        case DEFAULT

          WRITE(LUPRI,*) 'Keyword "',keyword,'" not understood in set_inforb_lucita'
          CALL QUIT('FATAL ERROR: Unrecognized keyword in set_inforb_lucita')

      end select

  end subroutine set_inforb_lucita
  
end module
