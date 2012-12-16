      integer :: IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS, NROOT, IDC
      integer :: INTSEL,IALIGN, c_state_of_interest
      integer :: root_converged, nfinal_vec
      integer :: irefsm_c, irefsm_hc
      real(8) :: PSSIGN,PLSIGN, ETHRES,CTHRES,EROOT,root_residual

      integer, parameter :: maxrts_lu = 150
     
      COMMON/CSTATE/IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,               &
     &              NROOT,IDC,INTSEL,IALIGN,c_state_of_interest,        &
     &              root_converged(maxrts_lu),nfinal_vec,               &
     &              irefsm_c, irefsm_hc,                                &
     &              PSSIGN,PLSIGN,ETHRES,CTHRES,EROOT(maxrts_lu),       &
     &              root_residual(maxrts_lu)
