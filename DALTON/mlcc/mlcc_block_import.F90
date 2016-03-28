module mlcc_block_import
!
!
!  mlcc types
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: import print units form Dalton
!
   integer, parameter   :: mxcorb = 2400
!  priunit, contains variables used for I/O
#if defined(__CVERSION__)
struct common_priunit {
  int lucmd, lupri, luerr, lustat, luw4, lupot, ninfo, nwarn, iprstatr;
};
extern struct common_priunit priunit_;
#else
!     FILE: priunit.h
      COMMON /PRIUNIT/ LUCMD,                                           &
     &        LUPRI, LUERR, LUSTAT, LUW4, LUPOT, NINFO, NWARN, IPRSTAT, &
     &        LQM3PCM, LPCMQM3
! ---  end of priunit.h ---
#endif
!
!     cc_tap contains information about the integral files
!
       COMMON /CC_TAP/ LUIAJB, LUBFDN, LURES, LUCSOL, LUPRPC,           &
     &                 LUMMMO, LUCCEF, LUCCPO
!
end module mlcc_block_import

