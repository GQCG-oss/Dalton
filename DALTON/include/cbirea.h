!
!     File: cbirea.h - common block for reading of MOLECULE.INP (herrdn.F)
!     Last revision: Apr 2010 /hjaaj
!
!     MAXPRD = default value for MAXPRI
      INTEGER MAXPRD
      INTEGER CBIREALAST, CBIREA_CLAST, CMMBASLAST
!
      PARAMETER ( MAXPRD = 35 )
!     MAXFAMEXP = maximum number of exponents in family basis sets
      INTEGER MXFAMEXP
      PARAMETER (MXFAMEXP = 100)

      REAL*8  ZCMVAL, TOL_SYMADD, FAMEXP(MXFAMEXP, 2), FAMPAR(4)
      INTEGER IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND, NFAMEXP(2)
      LOGICAL BIGVC,                  DIRAC,  BASIS,          ATOMBA,   &
     &        UNCONT, WRTLIN, ANGS,   ATOMDF, CNTMAT, NODDYDF
      LOGICAL LCNTNUUM
      COMMON /CBIREA/ ZCMVAL, TOL_SYMADD,     FAMEXP, FAMPAR,           &
     &        IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND, NFAMEXP,  &
     &        BIGVC,                  DIRAC,  BASIS,          ATOMBA,   &
     &        UNCONT, WRTLIN, ANGS,   ATOMDF, CNTMAT, NODDYDF,          &
     &        LCNTNUUM
!
      COMMON /CBIREA/ CBIREALAST
      !  Very important !!!
      !  Always keep CBIREALAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
!
!
!     Info for multiple basis sets (p.t. used with r12int.h in Dalton)

      INTEGER    MXMULB
      PARAMETER (MXMULB = 2)
      CHARACTER*80 MULNAM(MXMULB)
      COMMON /CBIREA_C/ MULNAM
!
      COMMON /CBIREA_C/ CBIREA_CLAST
      !  Very important !!!
      !  Always keep CBIREA_CLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
!
      INTEGER NMULBS
      LOGICAL LMULBS
      COMMON /CMMBAS/ NMULBS, LMULBS
!
      COMMON /CMMBAS/ CMMBASLAST
      !  Very important !!!
      !  Always keep CMMBASLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
! -- end of cbirea.h --
