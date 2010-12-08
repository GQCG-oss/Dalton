!     Common block for LUCITA orbital spaces
!     initialized in LUCITAINP
!     information is required for e.g. integral transformation.
!
!     nish_lucita(8): number of inactive (doubly occupied) shells
!     nash_lucita(8): total number of active shells (GAS and RAS case)
!     nas1_lucita(8): total number of active shells in RAS1 (RAS case only)
!     nas2_lucita(8): total number of active shells in RAS2 (RAS case only)
!     nas3_lucita(8): total number of active shells in RAS3 (RAS case only)
!     nocc_lucita(8): total number of occupied shells (nish_lucita + nash_lucita)
!     nssh_lucita(8): total number of secondary (virtual) shells
!
      INTEGER nish_lucita, nash_lucita, nas1_lucita, nas2_lucita,       &
     &        nas3_lucita, nocc_lucita, nssh_lucita
      COMMON/luciorbsp/nish_lucita(8), nash_lucita(8), nas1_lucita(8),  &
     &                 nas2_lucita(8), nas3_lucita(8), nocc_lucita(8),  &
     &                 nssh_lucita(8)
