C
C file: infvar.h (originally from SIRIUS)
C       contains INFormation on number and symmetry of mcscf VARiables
C
      INTEGER MAXWOP, JWOP,
     &        NCONF, NWOPT, NVAR, JWOPSY, NWOPH, NVARH, NCDETS
C
C     MAXWOP = maximum number of orbital rotations (dimension of JWOP)
C     JWOP(1:2,*) = /P, Q/, rotation from orbital P to orbital Q
C     NCONF  = number of configurations (CSFs or determinants)
C     NCDETS = number of determinants
C     NWOPT  = number of orbital rotations = number of non-redundant elements in JWOP
C     NVAR   = NCONF + NWOPT = total number of variables
C     JWOPSY = symmetry of orbital rotations
C     NWOPH  = NWOPT + no. of active-active rotations needed for Hessian
C     NVARH  = NCONF + NWOPH
C
      PARAMETER ( MAXWOP = 200 000 )
      COMMON /INFVAR/ NCONF,NWOPT,NVAR,JWOPSY,NWOPH,NVARH,NCDETS,
     &                JWOP(2,MAXWOP)
C --- end of infvar.h ---
