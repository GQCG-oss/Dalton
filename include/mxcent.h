C
C$Id: mxcent.h,v 1.2 2001-05-09 13:08:08 hjj Exp $
C
#if defined (VAR_ABASMALL)
      PARAMETER (MXCENT = 20, MXCOOR = 3*MXCENT)
#else
      PARAMETER (MXCENT = 80, MXCOOR = 3*MXCENT)
#endif
