#if defined (VAR_STAR2)
      INTEGER*2 JWOP
#else
      INTEGER JWOP
#endif
      INTEGER NCONF, NWOPT, NVAR, JWOPSY, NWOP, NWOPH, NVARH, MAXWOP, 
     *        NCDETS
#if defined (VAR_STAR2)
      INTEGER*2 JWOP
#endif
#if defined (VAR_STAR4)
      INTEGER*4 JWOP
#endif
C     MAXWOP = maximum number of rotations (dimension of JWOP)
      PARAMETER ( MAXWOP = 110000 )
      COMMON /INFVAR/ JWOP(2,MAXWOP),
     *                NCONF,NWOPT,NVAR,JWOPSY,NWOP(8),NWOPH,NVARH,NCDETS
