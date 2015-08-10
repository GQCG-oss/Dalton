      INTEGER  MAXSHL, NLRGBL, NSMLBL, NHKTSH, KHKTSH, KCKTSH,
     &         ISTBSH, NUCOSH, NORBSH, NSTRSH, NCNTSH, NSETSH,
     &         IORBSB, NRCSH,  LCLASH, NO2INT, NLRBL,  ISYMBL, NSYMBL,
     &         MBSISH, BLOCKSLAST

#if defined (SYS_CRAY)
      REAL CENTSH
#else
      DOUBLE PRECISION CENTSH
#endif

      LOGICAL BIGVEC, SEGMEN, SEGMSH, SPHRSH

      COMMON /BLOCKS/ CENTSH(MXSHEL,3),
     &                MAXSHL, BIGVEC, SEGMEN,NLRGBL,NSMLBL,
     &                NHKTSH(MXSHEL), KHKTSH(MXSHEL), KCKTSH(MXSHEL),
     &                ISTBSH(MXSHEL), NUCOSH(MXSHEL), NORBSH(MXSHEL),
     &                NSTRSH(MXSHEL), NCNTSH(MXSHEL), NSETSH(MXSHEL,2),
     &                IORBSB(0:MXCORB-1), NRCSH(MXSHEL), SEGMSH(MXSHEL),
     &                LCLASH(MXSHEL), SPHRSH(MXCORB),
     &                NLRBL,ISYMBL(MXSHEL,8),NSYMBL,
     &                MBSISH(MXSHEL)
C     MBSISH has been added for multiple basis sets (WK/UniKA/31-10-2002).

      COMMON /BLOCKS/ BLOCKSLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
