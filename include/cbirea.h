C
C    cbirea.h - common block for reading of MOLECULE.INP (herrdn.F)
C
C    MAXPRD = default value for MAXPRI
#if defined (SYS_T3D)
      PARAMETER ( MAXPRD = 14 )
#else
      PARAMETER ( MAXPRD = 25 )
#endif
C
      LOGICAL BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  ATOMBA,
     &        UNCONT, OLDNRM, WRTLIN, ANGS
      COMMON /CBIREA/ ZCMVAL,
     &                IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND,
     &                BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  ATOMBA,
     &                UNCONT, OLDNRM, WRTLIN, ANGS
      CHARACTER*80 BASNAM
      COMMON /CBIRET/ BASNAM
