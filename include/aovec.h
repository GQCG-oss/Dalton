C
C$Id: aovec.h,v 1.1.1.1 2001-02-08 13:33:26 hjj Exp $
C
#if defined (VAR_TESTIBM)
      PARAMETER (MXAOVC = 10, MXCONT = MXAOVC)
#else
      PARAMETER (MXAOVC = 25, MXCONT = MXAOVC)
#endif
