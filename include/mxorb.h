C
C$Id: mxorb.h,v 1.3 2001-12-13 08:54:39 ruden Exp $
C
C Note: MXCORB should be equal to MAXORB in maxorb.h
#if defined (VAR_ABASMALL)
      PARAMETER (MXSHEL = 200, MXPRIM = 800, MXCORB = 400,
     *           MXORBT = MXCORB*(MXCORB + 1)/2)
#else
      PARAMETER (MXSHEL = 750, MXPRIM = 4000, MXCORB = 1200,
     *           MXORBT = MXCORB*(MXCORB + 1)/2)
#endif
