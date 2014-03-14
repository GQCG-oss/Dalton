      INTEGER LBUF, NBUFS, NIBUF, NBITS, NBUFX, LUAORC, ERIBUFLAST 

      LOGICAL NEWDIS




      COMMON /ERIBUF/ LBUF, NBUFS, NIBUF, NBITS, NBUFX(0:MXCOOR),
     &                NEWDIS

      COMMON /ERIBUF/ ERIBUFLAST 
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)


      COMMON /ERITAP/ LUAORC(0:MXCOOR)
