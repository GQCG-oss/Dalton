#if defined(__CVERSION__)
struct common_priunt {
  int lucmd, lupri, luerr, luw4, ninfo, nwarn, iprerr;
};
extern struct common_priunt priunt_;
#else
C     FILE: priunit.h
      INTEGER LUCMD, LUPRI, LUERR , LUW4, NINFO, NWARN, IPRERR
      COMMON /PRIUNT/ LUCMD, LUPRI, LUERR, LUW4, NINFO, NWARN, IPRERR
#endif

