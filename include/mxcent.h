C
C$Id: mxcent.h,v 1.1.1.1 2001-02-08 13:33:26 hjj Exp $
C
#if defined (VAR_TESTIBM)
      PARAMETER (MXCENT =  8, MXCOOR = 3*MXCENT)
#else
#if defined (VAR_ABASMALL)
      PARAMETER (MXCENT = 20, MXCOOR = 3*MXCENT)
#else
      PARAMETER (MXCENT = 80, MXCOOR = 3*MXCENT)
#endif
#endif


