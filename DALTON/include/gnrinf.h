!
!     -*- mode: fortran; fortran-continuation-string: "&" -*-
!     File: gnrinf.h -- general information for DALTON
!
!
!     EMBEDDING : QM part is embedded in environment (solvent or e.g. protein)
!                 May 2011/hjaaj: EMBEDDING = FLAG(16) .or. PCM .or. QM3 .or. QMMM
!                 (For now, EMBEDDING is defined in sirius/sirinp.F because this is
!                 the first instance where all of FLAG(16), PCM, QM3, QMMM are set)
!
      LOGICAL TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,           &
     &        HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,           &
     &        WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,           &
     &        DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,           &
     &        TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,           &
     &        ERFEXP, DOSRIN, SRINTS, CHI1ST, DKHINT,                   &
     &        EMBEDDING, QM3, QMMM,   QMNPMM, QFIT,   USE_LSLIB
      REAL*8  GRADML, PANAS,  CHIVAL, THR_REDFAC
      INTEGER KCHARG, ITERNR, ITERMX, IPRUSR, LENBAS
      INTEGER GNRINFLAST
!
      COMMON /GNRINF/                                                   &
              ! real*8ZZ
     &        GRADML, PANAS,  CHIVAL, THR_REDFAC,                       &
              ! integer
     &        KCHARG, ITERNR, ITERMX, IPRUSR, LENBAS,                   &
              ! logical
     &        TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,           &
     &        HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,           &
     &        WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,           &
     &        DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,           &
     &        TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,           &
     &        ERFEXP, DOSRIN, SRINTS, CHI1ST, DKHINT,                   &
     &        EMBEDDING, QM3, QMMM,   QMNPMM, QFIT,   USE_LSLIB
!
      COMMON /GNRINF/ GNRINFLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, this variable has no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)

      INTEGER LBASDIR
      PARAMETER (LBASDIR = 600)
      CHARACTER*(LBASDIR) BASDIR
      CHARACTER*12        WFTYPE
      COMMON /GNRCHR/ BASDIR, WFTYPE
