C
C$Id: maxorb.h,v 1.1.1.1 2001-02-08 13:33:30 hjj Exp $
C
C     MAXORB = maximum number of orbitals
C     MAXOCC = maximum number of occupied orbitals
#if defined (VAR_SIRBIG)
      PARAMETER ( MAXORB = 1200, MAXOCC = 400 )
#else
      PARAMETER ( MAXORB = 400, MAXOCC = 120 )
#endif
