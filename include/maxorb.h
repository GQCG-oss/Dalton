C
C$Id: maxorb.h,v 1.2 2001-05-09 13:08:07 hjj Exp $
C
C     MAXORB = maximum number of orbitals
C             (Note: should be equal to MXCORB in mxorb.h)
C     MAXOCC = maximum number of occupied orbitals
#if defined (VAR_ABASMALL)
      PARAMETER ( MAXORB = 400, MAXOCC = 120 )
#else
      PARAMETER ( MAXORB = 1200, MAXOCC = 400 )
#endif
