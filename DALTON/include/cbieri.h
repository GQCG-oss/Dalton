C File : cbieri.h
C
C     Parameters NCBI? must be updated after changes (for parallelization)
C
C     NOTE: New logicals should appear after the last logical (NCLERI)
C           New integers should appear after the last integer (IFITDM)
C           Reals should appear at the end.
C
      INTEGER NCBII,  NCBIL
C     
C     IANGMO takes up space for 4 ints. That is not included in this counter
      PARAMETER (NCBII = 15,
     &           NCBIL = 31)

      INTEGER IPRERI, IAOBCH, IPRNT1, IPRNT2, IPROD1, IPROD2

      INTEGER LBFINP, MAXDST, NDMAT,  IANGMO, NSPMAX, MAXDSD,
     &        MXBCH,  MXBCH0, IFITDM, CBIERILAST

      LOGICAL RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI, OFFCNT,
     &        DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT, EXTPRI,
     &        INTPRI, DODIST, INTSKP, DISTST, WRTINT, FCKINT, PROFIL,
     &        NOLOCS, GRDZER, OLDDER, EXPERI, DOERIP, ERITWO, CCRUN,
     &        COMPRS, GENCON_ERI, NCLERI, DGABAB

      DIMENSION IANGMO(4)

      COMMON /CBIERI/ RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI,
     &                OFFCNT, DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT,
     &                EXTPRI, INTPRI, DODIST, INTSKP, DISTST, WRTINT,
     &                FCKINT, PROFIL, NOLOCS, GRDZER, OLDDER, EXPERI,
     &                DOERIP, ERITWO, CCRUN, COMPRS, GENCON_ERI, NCLERI,
     &                DGABAB, NDMAT,  IPROD1, IPROD2, IAOBCH, IPRNT1,
     &                IPRNT2, IPRERI, LBFINP, MAXDST, IANGMO, NSPMAX, 
     &                MAXDSD, MXBCH,  MXBCH0, IFITDM

      COMMON /CBIERI/ CBIERILAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
C -- end of cbieri.h --
