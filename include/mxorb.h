C
C$Id: mxorb.h,v 1.1.1.1 2001-02-08 13:33:29 hjj Exp $
C
#if defined (VAR_TESTIBM)
      PARAMETER (MXSHEL =  15, MXPRIM =  60, MXCORB =  30,
     *           MXORBT = MXCORB*(MXCORB + 1)/2)
#else
#if defined (VAR_SIRBIG)
      PARAMETER (MXSHEL = 750, MXPRIM = 4000, MXCORB = 1800,
     *           MXORBT = MXCORB*(MXCORB + 1)/2)
#else
      PARAMETER (MXSHEL = 200, MXPRIM = 800, MXCORB = 400,
     *           MXORBT = MXCORB*(MXCORB + 1)/2)
#endif
#endif
