#if defined(__CVERSION__)
struct common_priunit {
  int lucmd, lupri, luerr, lustat, luw4, lupot, ninfo, nwarn, iprstatr;
};
extern struct common_priunit priunit_;
#else
!     FILE: priunit.h
      CHARACTER*80 SEPARATOR
      PARAMETER (SEPARATOR = '----------------------------------------' &
     &                     //'----------------------------------------')
      INTEGER LUCMD
      INTEGER LUPRI, LUERR, LUSTAT, LUW4, LUPOT, NINFO, NWARN, IPRSTAT, &
     &        LQM3PCM, LPCMQM3
      COMMON /PRIUNIT/ LUCMD,                                           &
     &        LUPRI, LUERR, LUSTAT, LUW4, LUPOT, NINFO, NWARN, IPRSTAT, &
     &        LQM3PCM, LPCMQM3
! ---  end of priunit.h ---
#endif
