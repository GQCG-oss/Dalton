      LOGICAL SOLVNT
      INTEGER LCAVMX,LMTOT,LMNTOT,NCNTCV
      REAL*8 EPDIEL,RCAV
      COMMON /CBISOL/ EPDIEL,RCAV(3),LCAVMX,LMTOT,LMNTOT,NCNTCV,SOLVNT
C     ... LCAVMX = max l in moment expansion for solvent cavity
C         LMTOT  = number of spherical components for R(l,m)
C         LMNTOT = number of cartesian components for RC(L,M,N)
