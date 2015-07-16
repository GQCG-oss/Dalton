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
      LOGICAL BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  PRIBAS, ATOMBA,   &
     &        UNCONT, WRTLIN, ANGS,   ATOMDF, CNTMAT, NODDYDF
      LOGICAL LCNTNUUM
      COMMON /CBIREA/ ZCMVAL, TOL_SYMADD,     FAMEXP, FAMPAR,           &
     &        IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND, NFAMEXP,  &
     &        BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  PRIBAS, ATOMBA,   &
     &        UNCONT, WRTLIN, ANGS,   ATOMDF, CNTMAT, NODDYDF,          &
     &        LCNTNUUM
!
      COMMON /CBIREA/ CBIREALAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
!
!
!     Info for multiple basis sets (p.t. used with r12int.h in Dalton)

      INTEGER    MXMULB
      PARAMETER (MXMULB = 2)
      CHARACTER*80 MULNAM(MXMULB)
      COMMON /CBIREA_C/ MULNAM
!
      COMMON /CBIREA_C/ CBIREA_CLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
!
      INTEGER NMULBS
      LOGICAL LMULBS
      COMMON /CMMBAS/ NMULBS, LMULBS
!
      COMMON /CMMBAS/ CMMBASLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
! -- end of cbirea.h --
