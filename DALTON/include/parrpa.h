C     PARRPA = Parallel Random Phase Approximation
C
      integer :: NUMDIS, KODCL1, KODCL2, KODBC1, 
     &           KODBC2, KRDBC1, KRDBC2, KODPP1,
     &           KODPP2, KRDPP1, KRDPP2, KCCFB1,
     &           KINDXB, NTOT,
     &           KCMO  , KEND1 , LWORK1,
     &           LTR1E , LTR1D , LRES1E, LRES1D,
     &           LFOCK , LDENS , LBTR1E, LBTR1D, 
     &           LBTJ1E, LBTJ1D, KTR1E , KTR1D ,
     &           KRES1E, KRES1D, KFOCK , KDENS ,
     &           KBTR1E, KBTR1D, KBTJ1E, KBTJ1D,
     &           KEND2,  LWORK2, KENDSV, LWORKSV,
     &           ISYMPAR, parrpalast
      logical    forceupdate
C
C
      common /parrpa/
     & NTOT  , NUMDIS, KODCL1, KODCL2, KODBC1, LWORKSV,
     & KODBC2, KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1, 
     & KRDPP2, KCCFB1, KINDXB, KCMO  , KEND1 , LWORK1,
     & LTR1E , LTR1D , LRES1E, LRES1D, LFOCK , LDENS , 
     & LBTR1E, LBTR1D, LBTJ1E, LBTJ1D, KTR1E , KTR1D ,
     & KRES1E, KRES1D, KFOCK , KDENS , KBTR1E, KBTR1D, 
     & KBTJ1E, KBTJ1D, KEND2 , LWORK2, KENDSV, ISYMPAR,
     & forceupdate
C
C
C
      common /parrpa/ parrpalast
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
C
      save /parrpa/
