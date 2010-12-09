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
!     ngso_lucita(:,2): min and max # accumulated electrons per gas space

  use lucita_cfg

  implicit none

  public define_lucita_orb_spaces

  private

contains 
 
  subroutine define_lucita_orb_spaces(number_of_ptg_irreps,        &
                                      number_of_gas_spaces,        &
                                      mx_number_of_ptg_irreps,     &
                                      mx_number_of_gas_spaces,     &
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
    integer, intent(in) :: mx_number_of_ptg_irreps
    integer, intent(in) :: mx_number_of_gas_spaces
    integer, intent(in) :: nr_orb_tot(mx_number_of_ptg_irreps)
    integer, intent(in) :: init_input_type
    integer, intent(in) :: init_wave_f_type
!-------------------------------------------------------------------------------

!     (reset and) initialize lucita orbital space arrays
      call lucita_orb_spaces_init(mx_number_of_ptg_irreps,mx_number_of_gas_spaces)

!     transfer information from temporary arrays provided in: 
!     case 1: lucita input 
!     case 2:  mcscf input
!     ------------------------------------------
      select case(init_input_type)

        case(1)

          if(init_wave_f_type == 2)then
            allocate(ngsh_lucita(mx_number_of_gas_spaces,mx_number_of_ptg_irreps))
          end if

          call fill_lucita_orb_spaces_ci(number_of_ptg_irreps,        &
                                         number_of_gas_spaces,        &
                                         mx_number_of_ptg_irreps,     &
                                         mx_number_of_gas_spaces,     &
                                         nr_orb_tot,                  &
                                         nish_lucita,                 &
                                         nash_lucita,                 &
                                         nocc_lucita,                 &
                                         nssh_lucita,                 &
                                         nas1_lucita,                 &
                                         nas2_lucita,                 &
                                         nas3_lucita,                 &
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
                                       mx_number_of_ptg_irreps,     &
                                       mx_number_of_gas_spaces,     &
                                       nr_orb_tot,                  &
                                       is_nish_lucita,              &
                                       is_nash_lucita,              &
                                       is_nocc_lucita,              &
                                       is_nssh_lucita,              &
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
    integer, intent(in)    :: nr_orb_tot(mx_number_of_ptg_irreps)
    integer, intent(in)    :: init_wave_f_type
    integer, intent(in)    :: is_nish_lucita(mx_number_of_ptg_irreps)
    integer, intent(in)    :: is_nas1_lucita(mx_number_of_ptg_irreps)
    integer, intent(in)    :: is_nas2_lucita(mx_number_of_ptg_irreps)
    integer, intent(in)    :: is_nas3_lucita(mx_number_of_ptg_irreps)
    integer, intent(out)   :: is_nash_lucita(mx_number_of_ptg_irreps)
    integer, intent(out)   :: is_nocc_lucita(mx_number_of_ptg_irreps)
    integer, intent(out)   :: is_nssh_lucita(mx_number_of_ptg_irreps)
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
              call icopy(mx_number_of_ptg_irreps,                &
                         is_nas1_lucita,1,is_ngsh_lucita(i,1),1)
            case(2) 
              call icopy(mx_number_of_ptg_irreps,                &
                         is_nas2_lucita,1,is_ngsh_lucita(i,1),1)
            case(3) 
              call icopy(mx_number_of_ptg_irreps,                &
                         is_nas3_lucita,1,is_ngsh_lucita(i,1),1)
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

!   d. total number of secondary shells for each point group irrep
    do i = 1, number_of_ptg_irreps
      is_nssh_lucita(i) = nr_orb_tot(i) - is_nocc_lucita(i)
    end do

    print *, ' debug print of orbital spaces:'
    print *, ' nish_lucita        : ',(is_nish_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nash_lucita        : ',(is_nash_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nocc_lucita        : ',(is_nocc_lucita(i),i=1,number_of_ptg_irreps)
    print *, ' nssh_lucita        : ',(is_nssh_lucita(i),i=1,number_of_ptg_irreps)
    do j = 1, tmp_nr_gas_spaces
      print *, ' ngsh_lucita per gas: ',(is_ngsh_lucita(j,i),i=1,number_of_ptg_irreps)
    end do

  end subroutine fill_lucita_orb_spaces_ci
!*******************************************************************************
 
  subroutine lucita_orb_spaces_init(mx_number_of_ptg_irreps,    &
                                    mx_number_of_gas_spaces)
!*******************************************************************************
!
!    purpose:  initialize (if necessary) the orbital spaces for LUCITA
!              calculations.
!
!*******************************************************************************
      integer, intent(in) :: mx_number_of_ptg_irreps
      integer, intent(in) :: mx_number_of_gas_spaces
!-------------------------------------------------------------------------------

!     reset arrays (if necessary)
      call lucita_orb_spaces_setup_delete_subset()

!     allocate CI orbital space arrays
      allocate(nash_lucita(mx_number_of_ptg_irreps))
      allocate(nocc_lucita(mx_number_of_ptg_irreps))
      allocate(nssh_lucita(mx_number_of_ptg_irreps))

!     initialize CI orbital spaces
      nash_lucita = 0
      nocc_lucita = 0
      nssh_lucita = 0

  end subroutine lucita_orb_spaces_init
!*******************************************************************************

  subroutine lucita_orb_spaces_setup_delete_subset
!*******************************************************************************
!
!    purpose: deallocate selected CI orbital space arrays for LUCITA, 
!             i.e. those which are not allocated in the input reader.
!
!*******************************************************************************
!-------------------------------------------------------------------------------

      if(allocated(nash_lucita))deallocate(nash_lucita)
      if(allocated(nocc_lucita))deallocate(nocc_lucita)
      if(allocated(nssh_lucita))deallocate(nssh_lucita)

  end subroutine lucita_orb_spaces_setup_delete_subset
!*******************************************************************************
  
end module
