!
!     File: cbirea.h - common block for reading of MOLECULE.INP (herrdn.F)
!     Last revision: Apr 2010 /hjaaj
!
!     MAXPRD = default value for MAXPRI
      INTEGER MAXPRD
      INTEGER CBIREAlast, CBIREA_Clast, CMMBASlast
      PARAMETER ( MAXPRD = 35 )
!     MAXFAMEXP = maximum number of exponents in family basis sets
      INTEGER MXFAMEXP
      PARAMETER (MXFAMEXP = 100)

      REAL*8  ZCMVAL, TOL_SYMADD, FAMEXP(MXFAMEXP, 2), FAMPAR(4)
      INTEGER IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND, NFAMEXP(2)
      LOGICAL BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  PRIBAS, ATOMBA,   &
     &        UNCONT, WRTLIN, ANGS,   ATOMDF, CNTMAT, NODDYDF,LCNTNUUM
      COMMON /CBIREA/ ZCMVAL, TOL_SYMADD,     FAMEXP, FAMPAR,           & ! real*8
     &        IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND, NFAMEXP,  & ! integer
     &        BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  PRIBAS, ATOMBA,   & ! logical
     &        UNCONT, WRTLIN, ANGS,   ATOMDF, CNTMAT, NODDYDF,LCNTNUUM, &
     &   CBIREAlast !  Very important:
      !  Always keep CBIREAlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
!
!     Info for multiple basis sets (p.t. used with r12int.h in Dalton)

      INTEGER    MXMULB
      PARAMETER (MXMULB = 2)
      CHARACTER*80 MULNAM(MXMULB)
      COMMON /CBIREA_C/ MULNAM,                                         &
     &   CBIREA_Clast !  Very important:
      !  Always keep CBIREA_Clast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

      INTEGER NMULBS
      LOGICAL LMULBS
      COMMON /CMMBAS/ NMULBS, LMULBS,                                   &
     &   CMMBASlast !  Very important:
      !  Always keep CMBASlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
! -- end of cbirea.h --
