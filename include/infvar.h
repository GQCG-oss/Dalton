#if defined (VAR_STAR2)
      INTEGER*2 JWOP
#endif
#if defined (VAR_STAR4)
      INTEGER*4 JWOP
#endif
C     MAXWOP = maximum number of rotations (dimension of JWOP)
#if defined (VAR_ABASMALL)
      PARAMETER ( MAXWOP = 5 000 )
#else
      PARAMETER ( MAXWOP = 110 000 )
#endif
      COMMON /INFVAR/ JWOP(2,MAXWOP),
     *                NCONF,NWOPT,NVAR,JWOPSY,NWOP(8),NWOPH,NVARH,NCDETS
