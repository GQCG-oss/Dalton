C
C$Id: mxcent.h,v 1.3 2001-12-13 08:54:39 ruden Exp $
C
#if defined (VAR_ABASMALL)
      PARAMETER (MXCENT = 20, MXCOOR = 3*MXCENT)
#else
      PARAMETER (MXCENT = 80, MXCOOR = 3*MXCENT)
#endif
