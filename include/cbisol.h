C
C$Id: cbisol.h,v 1.1.1.1 2001-02-08 13:33:25 hjj Exp $
C
      LOGICAL SOLVNT
      COMMON /CBISOL/ RCAV(3),EPDIEL,
     &                LCAVMX,LMTOT,LMNTOT,NCNTCV,SOLVNT
C     ... LCAVMX = max l in moment expansion for solvent cavity
C         LMTOT  = number of spherical components for R(l,m)
C         LMNTOT = number of cartesian components for RC(L,M,N)
