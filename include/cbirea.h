C
C     cbirea.h - common block for reading of MOLECULE.INP (herrdn.F)
C
C     MAXPRD = default value for MAXPRI
      PARAMETER ( MAXPRD = 35 )
C
      LOGICAL         BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  ATOMBA,
     &                UNCONT, OLDNRM, WRTLIN, ANGS,   ATOMDF, NODDYDF
      COMMON /CBIREA/ ZCMVAL, TOL_SYMADD,
     &                IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND,
     &                BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  ATOMBA,
     &                UNCONT, OLDNRM, WRTLIN, ANGS,   ATOMDF, NODDYDF
      CHARACTER*80    BASNAM,AUXNAM
      COMMON /CBIRET/ BASNAM,AUXNAM
