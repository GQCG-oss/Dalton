!
!------------------------------------------------------------
! information for the packing of the AO integrals on the
! presorted CCAOIN_* files
!------------------------------------------------------------
!
      LOGICAL LPACKINT
      INTEGER IPCKTABINT(0:255)
      INTEGER IOFFINT(MXCORB)
      INTEGER NPCKINT(MXCORB)
      INTEGER NTOTINT, NTOTPCK
      INTEGER CCPACKLAST
!
!
!
#if defined (SYS_CRAY)
      REAL THRPCKINT, PCKRATIO, PCKTIME
#else
      DOUBLE PRECISION THRPCKINT, PCKRATIO, PCKTIME
#endif
      COMMON /CCPACK/ THRPCKINT, PCKRATIO, PCKTIME,                     &
     &                IPCKTABINT, IOFFINT, NPCKINT,                     &
     &                NTOTINT, NTOTPCK, LPACKINT
!
      COMMON /CCPACK/ CCPACKLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
