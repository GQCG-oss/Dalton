      LOGICAL SKIP,   DIRECT, DIRGRD, CCRSTR,
     &        FROIMP, FROEXP, NOCCIT,
     &        CCSAVE, STOLD,  JACEXP, LHTR,
     &        DEBUG, CCSTST, ANAAOD,
     &        HERDIR, FREEZE, KEEPAOIN, NOEONL, NOSORT,
     &        HURWITZ_CHECK,
     &        SIRSOP, LVVVV, ONLYMO, AOSOPPA,
     &        CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &        MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &        CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &        rCCD, drCCD, SOSEX, rTCCD,
     &        CCSDT,CCR12, MTRIP, CHOPT, DCPT2,
     &        ETACCPT, DIRKAPB

      INTEGER MXDIIS, MXLRV,
     &        IT2UPD, IT2START,
     &        ICHANG, IPRINT, KEEPAOTWO, CCSDGNINPLAST

      COMMON /CCSDGNINP/ SKIP, DIRECT, DIRGRD, CCRSTR,
     &                   FROIMP,  FROEXP, NOCCIT,
     &                   CCSAVE,  STOLD,  JACEXP, LHTR,
     &                   DEBUG,   CCSTST, ANAAOD,
     &                   MXDIIS,  MXLRV,   IT2UPD, IT2START,
     &                   ICHANG,  IPRINT,  KEEPAOTWO, HERDIR,
     &                   ETACCPT, DIRKAPB,
     &                   FREEZE,  KEEPAOIN, NOEONL, NOSORT,
     &                   HURWITZ_CHECK,
     &                   SIRSOP,  LVVVV, ONLYMO, AOSOPPA


      COMMON /CCSDGNINP/ CCSDGNINPLAST 
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)


      INTEGER CCMODELSLAST
      COMMON /CCMODELS/ CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &                  MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &                  CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &                  rCCD, drCCD, SOSEX, rTCCD,
     &                  CCSDT,CCR12,MTRIP, CHOPT, DCPT2
      COMMON /CCMODELS/ CCMODELSLAST 
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)





C For the RPA models we can calcualte several energies using the
C same set of amplitudes -- this temporary variable is used to store
C one of the energies, the usual variables are reserved for use with
c the choice that could be used in subsequent response or gradient
C calculations
C
C
C
C
      INTEGER ETMPLAST
      DOUBLE PRECISION  ETMP
      COMMON/ETMP/ ETMP
C
      COMMON /ETMP/ ETMPLAST 
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)


