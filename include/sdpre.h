C
C$Id: sdpre.h,v 1.1.1.1 2001-02-08 13:33:23 hjj Exp $
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      SDPRE(I) = FLOAT(I)
#else
      SDPRE(I) = DFLOAT(I)
#endif
