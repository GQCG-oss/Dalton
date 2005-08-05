#if defined (VAR_STAR2)
      INTEGER*2 JWOP
#else
      INTEGER JWOP
#endif
      INTEGER NCONF, NWOPT, NVAR, JWOPSY, NWOP, NWOPH, NVARH, MAXWOP, 
     *        NCDETS
C     MAXWOP = maximum number of rotations (dimension of JWOP)
      PARAMETER ( MAXWOP = 110000 )
      COMMON /INFVAR/ JWOP(2,MAXWOP),
     *                NCONF,NWOPT,NVAR,JWOPSY,NWOP(8),NWOPH,NVARH,NCDETS
