!dalton_copyright_start
!
!
!dalton_copyright_end

! RMA windows object (used for MPI-one-sided communication).
!
! written by sknecht for Dalton - linkoeping jan 2013
!
!
module rma_windows

  implicit none

! rma-windows definition
! ----------------------------------------------------------------------------
  type rma_win_dalton

    integer ::                   &
      dmat_win,                  &             ! density matrix window (hermit)
      fmat_win                                 ! fock matrix window (hermit)

    logical ::                   &
      rma_win_dalton_init = .false.            ! status of the rma_win_dalton type

  end type rma_win_dalton

! rma-windows object
  type(rma_win_dalton), public, save :: rma_win_dalton_info
! ----------------------------------------------------------------------------

end module rma_windows
