C
C$Id: scalar.h,v 1.1.1.1 2001-02-08 13:33:28 hjj Exp $
C
#if defined (SYS_ALLIANT)
CVD$ NOVECTOR
#endif
#if defined (SYS_CONVEX) || defined (SYS_HAL)
C$DIR SCALAR
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
CDIR$ NEXTSCALAR
#endif
#if defined (SYS_NEC)
*VDIR NOVECTOR
#endif
