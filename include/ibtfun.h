C
C$Id: ibtfun.h,v 1.1.1.1 2001-02-08 13:33:26 hjj Exp $
C
#if defined (SYS_VAX) || defined (SYS_ALLIANT) || defined (SYS_IBM) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_HPUX) || defined (SYS_SUN) || defined (SYS_NEC) || defined (SYS_HAL)
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = IEOR(I,J)
#endif
#if defined (SYS_LINUX)
#define IBTAND(I,J) IAND(I,J)
#define IBTOR(I,J)  IOR(I,J)
#define IBTSHL(I,J) ISHFT(I,J)
#define IBTSHR(I,J) ISHFT(I,-J)
#define IBTXOR(I,J) IEOR(I,J)
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = SHIFTL(I,J)
      IBTSHR(I,J) = SHIFTR(I,J)
      IBTXOR(I,J) = XOR(I,J)
#endif
#if !defined (SYS_VAX) && !defined (SYS_ALLIANT) && !defined (SYS_IBM) && !defined (SYS_CONVEX) && !defined (SYS_CRAY) && !defined (SYS_AIX) && !defined (SYS_PARAGON) && !defined (SYS_DEC) && !defined (SYS_IRIX) && !defined (SYS_HPUX) && !defined (SYS_SUN) && !defined (SYS_T3D) && !defined (SYS_LINUX) && !defined (SYS_NEC) && !defined (SYS_HAL) && !defined (SYS_T90)
      You must define IBTFUN in comdeck file for this computer.
#endif
