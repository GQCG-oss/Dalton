!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_orbital_spaces

! stefan: - this module provides all necessary functionality
!           to setup the orbital spaces in LUCITA mcscf/ci 
!           calculations.
!
!           written by sknecht for DALTON, december 2010.
!
!     nish_lucita(:)  : total number of inactive (doubly occupied) shells per irrep
!     nash_lucita(:)  : total number of active shells (GAS and RAS case) per irrep
!     nas1_lucita(:)  : total number of active shells in RAS1 (RAS case only) per irrep
!     nas2_lucita(:)  : total number of active shells in RAS2 (RAS case only) per irrep
!     nas3_lucita(:)  : total number of active shells in RAS3 (RAS case only) per irrep
!     nocc_lucita(:)  : total number of occupied shells (nish_lucita + nash_lucita) per irrep
!     nssh_lucita(:)  : total number of secondary (virtual) shells per irrep
!     ngsh_lucita(:,:): total number of active shells per gas space and irrep

  implicit none

  public define_lucita_orb_spaces
  public lucita_orb_spaces_setup_delete

  private

  integer, allocatable, public :: nish_lucita(:)
  integer, allocatable, public :: nash_lucita(:)
  integer, allocatable, public :: nas1_lucita(:)
  integer, allocatable, public :: nas2_lucita(:)
  integer, allocatable, public :: nas3_lucita(:)
  integer, allocatable, public :: nocc_lucita(:)
  integer, allocatable, public :: nssh_lucita(:)
  integer, allocatable, public :: ngsh_lucita(:,:)

  integer, private             :: my_dummy

contains 
 
  subroutine define_lucita_orb_spaces(number_of_ptg_irreps,        &
                                      number_of_gas_spaces,        &
                                      max_number_of_ptg_irreps,    &
                                      max_number_of_gas_spaces,    &
                                      tmp_inactive_orb_space,      &
                                      tmp_active_orb_space,        &
                                      nr_orb_tot,                  &
                                      init_input_type,             &
                                      init_wave_f_type)
!*******************************************************************************
!
!    purpose:  define the orbital spaces for LUCITA ci/mcscf calculations.
!
!*******************************************************************************
    integer, intent(in) :: number_of_ptg_irreps
    integer, intent(in) :: number_of_gas_spaces
    integer, intent(in) :: max_number_of_ptg_irreps
    integer, intent(in) :: max_number_of_gas_spaces
    integer, intent(in) :: tmp_inactive_orb_space(max_number_of_ptg_irreps)
    integer, intent(in) :: tmp_active_orb_space(max_number_of_gas_spaces,max_number_of_ptg_irreps)
    integer, intent(in) :: nr_orb_tot(max_number_of_ptg_irreps)
    integer, intent(in) :: init_input_type
    integer, intent(in) :: init_wave_f_type
!-------------------------------------------------------------------------------

!     (reset and) initialize lucita orbital space arrays
      call lucita_orb_spaces_init(max_number_of_ptg_irreps,max_number_of_gas_spaces)

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
                                         tmp_inactive_orb_space,      &
                                         tmp_active_orb_space,        &
                                         nr_orb_tot,                  &
                                         nish_lucita,                 &
                                         nash_lucita,                 &
                                         nas1_lucita,                 &
                                         nas2_lucita,                 &
                                         nas3_lucita,                 &
                                         nocc_lucita,                 &
                                         nssh_lucita,                 &
                                         ngsh_lucita,                 &
                                         init_wave_f_type)

        case(2)

          call quit(' *** error in lucita_orb_spaces_init: mcscf orbital& 
 input conversion not available yet. ***')
          
!         call fill_lucita_orb_spaces_mcscf

      end select

  end subroutine define_lucita_orb_spaces
!*******************************************************************************

  subroutine fill_lucita_orb_spaces_ci(number_of_ptg_irreps,        &
                                       number_of_gas_spaces,        &
                                       max_number_of_ptg_irreps,    &
                                       max_number_of_gas_spaces,    &
                                       tmp_inactive_orb_space,      &
                                       tmp_active_orb_space,        &
                                       nr_orb_tot,                  &
                                       nish_lucita,                 &
                                       nash_lucita,                 &
                                       nas1_lucita,                 &
                                       nas2_lucita,                 &
                                       nas3_lucita,                 &
                                       nocc_lucita,                 &
                                       nssh_lucita,                 &
                                       ngsh_lucita,                 &
                                       init_wave_f_type)
!*******************************************************************************
!
!    purpose:  fill orbital spaces for LUCITA from ci input.
!
!*******************************************************************************
    integer, intent(in) :: number_of_ptg_irreps
    integer, intent(in) :: number_of_gas_spaces
    integer, intent(in) :: max_number_of_ptg_irreps
    integer, intent(in) :: max_number_of_gas_spaces
    integer, intent(in) :: tmp_inactive_orb_space(max_number_of_ptg_irreps)
    integer, intent(in) :: tmp_active_orb_space(max_number_of_gas_spaces,max_number_of_ptg_irreps)
    integer, intent(in) :: nr_orb_tot(max_number_of_ptg_irreps)
    integer, intent(in) :: init_wave_f_type
    integer, intent(out):: nish_lucita(max_number_of_ptg_irreps)
    integer, intent(out):: nash_lucita(max_number_of_ptg_irreps)
    integer, intent(out):: nas1_lucita(max_number_of_ptg_irreps)
    integer, intent(out):: nas2_lucita(max_number_of_ptg_irreps)
    integer, intent(out):: nas3_lucita(max_number_of_ptg_irreps)
    integer, intent(out):: nocc_lucita(max_number_of_ptg_irreps)
    integer, intent(out):: nssh_lucita(max_number_of_ptg_irreps)
    integer, intent(out):: ngsh_lucita(max_number_of_gas_spaces,max_number_of_ptg_irreps)
!-------------------------------------------------------------------------------
    integer             :: tmp_nr_gas_spaces
    integer             :: i
    integer             :: j
!-------------------------------------------------------------------------------

!   a. inactive shells for each point group irrep
    call icopy(number_of_ptg_irreps,tmp_inactive_orb_space,1,nish_lucita,1)

    select case(init_wave_f_type)

      case(1) ! GAS
       
        tmp_nr_gas_spaces = number_of_gas_spaces

!       b. active shells for each point group irrep and gas space
        call icopy(tmp_nr_gas_spaces*number_of_ptg_irreps,      &
                   tmp_active_orb_space,1,ngsh_lucita,1)
      
      case(2) ! RAS

        tmp_nr_gas_spaces = 3

!       b1. active shells for each point group irrep and ras space
        call icopy(tmp_nr_gas_spaces*number_of_ptg_irreps,      &
                   tmp_active_orb_space,1,ngsh_lucita,1)

!       b2. active shells for each point group irrep in RAS1, RAS2 and RAS3
        do i = 1, tmp_nr_gas_spaces
          select case(i)
            case(1) 
              call icopy(number_of_ptg_irreps,                  &
                         ngsh_lucita(i,1),1,nas1_lucita,1)
            case(2) 
              call icopy(number_of_ptg_irreps,                  &
                         ngsh_lucita(i,1),1,nas2_lucita,1)
            case(3) 
              call icopy(number_of_ptg_irreps,                  &
                         ngsh_lucita(i,1),1,nas3_lucita,1)
          end select
        end do
 
    end select

!   c. total number of active shells for each point group irrep
    do i = 1, tmp_nr_gas_spaces
      do j = 1, number_of_ptg_irreps
        nash_lucita(j) = nash_lucita(j) + ngsh_lucita(i,j)
      end do
    end do

!   d. total number of occupied shells for each point group irrep
    do i = 1, number_of_ptg_irreps
      nocc_lucita(i) = nash_lucita(i) + nish_lucita(i)
    end do

!   e. total number of secondary shells for each point group irrep
    do i = 1, number_of_ptg_irreps
      nssh_lucita(i) = nr_orb_tot(i) - nocc_lucita(i)
    end do

    print *, ' debug print of orbital spaces:'
    print *, ' nish_lucita        : ',(nish_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nash_lucita        : ',(nash_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nocc_lucita        : ',(nocc_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nssh_lucita        : ',(nssh_lucita(i),i=1,number_of_ptg_irreps)
    do j = 1, tmp_nr_gas_spaces
      print *, ' ngsh_lucita per gas: ',(ngsh_lucita(j,i),i=1,number_of_ptg_irreps)
    end do

  end subroutine fill_lucita_orb_spaces_ci
!*******************************************************************************
 
  subroutine lucita_orb_spaces_init(max_number_of_ptg_irreps,    &
                                    max_number_of_gas_spaces)
!*******************************************************************************
!
!    purpose:  initialize (if necessary) the orbital spaces for LUCITA
!              calculations.
!
!*******************************************************************************
      integer, intent(in) :: max_number_of_ptg_irreps
      integer, intent(in) :: max_number_of_gas_spaces
!-------------------------------------------------------------------------------

!     reset arrays (if necessary)
      call lucita_orb_spaces_setup_delete

!     allocate CI orbital space arrays
      allocate(nish_lucita(max_number_of_ptg_irreps))
      allocate(nash_lucita(max_number_of_ptg_irreps))
      allocate(nas1_lucita(max_number_of_ptg_irreps))
      allocate(nas2_lucita(max_number_of_ptg_irreps))
      allocate(nas3_lucita(max_number_of_ptg_irreps))
      allocate(nocc_lucita(max_number_of_ptg_irreps))
      allocate(nssh_lucita(max_number_of_ptg_irreps))
      allocate(ngsh_lucita(max_number_of_gas_spaces,max_number_of_ptg_irreps))

!     initialize CI orbital spaces
      nish_lucita = 0
      nash_lucita = 0
      nas1_lucita = 0
      nas2_lucita = 0
      nas3_lucita = 0
      nocc_lucita = 0
      nssh_lucita = 0
      ngsh_lucita = 0

  end subroutine lucita_orb_spaces_init
!*******************************************************************************

  subroutine lucita_orb_spaces_setup_delete
!*******************************************************************************
!
!    purpose: deallocate all CI orbital space arrays for LUCITA.
!
!*******************************************************************************
!-------------------------------------------------------------------------------

     !if(allocated(nish_lucita))deallocate(nish_lucita)
     !if(allocated(nash_lucita))deallocate(nash_lucita)
     !if(allocated(nas1_lucita))deallocate(nas1_lucita)
     !if(allocated(nas2_lucita))deallocate(nas2_lucita)
     !if(allocated(nas3_lucita))deallocate(nas3_lucita)
     !if(allocated(nocc_lucita))deallocate(nocc_lucita)
     !if(allocated(nssh_lucita))deallocate(nssh_lucita)
     !if(allocated(ngsh_lucita))deallocate(ngsh_lucita)
      if(allocated(ngsh_lucita)) then 
        print *, ' ngsh_lucita is allocated'
        print *, ' shells per gas',ngsh_lucita(1,1),ngsh_lucita(1,2),ngsh_lucita(1,3),& 
        ngsh_lucita(1,4),ngsh_lucita(1,5),ngsh_lucita(1,6),ngsh_lucita(1,7), ngsh_lucita(1,8)
      end if


  end subroutine lucita_orb_spaces_setup_delete
!*******************************************************************************
  
end module
