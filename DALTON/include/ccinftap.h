!      INTEGER         LUIAJB, LUCKJD, LU3VI2, LUDKBC, LU3VI, LUDELD,
!     &                LUTOC, LU3SRT, LUBFDN, LURES,
!     &                LUCSOL, LUCTLM, LU3VI3, LU3VI4, LU3SRT2, LU3SRT3,
!     &                LUCKJD2, LUCKJD3, LUDELD2, LUDELD3, LUDELD4,
!     &                LUDELD5, LUDKBC2, LUDKBC3, LUDKBC4, LUDKBC5
!      CHARACTER       FNTOC*8, FN3VI*6, FN3VI2*8, FN3SRT*8,
!     &                FNDELD*6, FNDKBC*4, FNCKJD*6,
!     &                FN3VI3*8, FN3VI4*8, FN3SRT2*9,
!     &                FN3SRT3*9, FNCKJD2*7, FNCKJD3*7,
!     &                FNDELD2*7, FNDELD3*7, FNDELD4*7, FNDELD5*7,
!     &                FNDKBC2*5, FNDKBC3*5, FNDKBC4*5, FNDKBC5*5
!      COMMON /CC_TAP/ LUIAJB, LUCKJD, LU3VI2, LUDKBC, LU3VI, LUDELD,
!     &                LUTOC, LU3SRT, LUBFDN, LURES, LUCSOL,
!     &                LUCTLM, LU3VI3, LU3VI4, LU3SRT2, LU3SRT3,
!     &                LUCKJD2, LUCKJD3, LUDELD2, LUDELD3, LUDELD4,
!     &                LUDELD5, LUDKBC2, LUDKBC3, LUDKBC4, LUDKBC5,
!     &                FNTOC, FN3VI,
!     &                FN3VI2, FN3SRT, FNDELD, FNDKBC, FNCKJD,
!     &                FN3VI3, FN3VI4, FN3SRT2, FN3SRT3, FNCKJD2,
!     &                FNCKJD3, FNDELD2, FNDELD3, FNDELD4, FNDELD5,
!     &                FNDKBC2, FNDKBC3, FNDKBC4, FNDKBC5

       INTEGER         LUIAJB, LUBFDN, LURES, LUCSOL, LUPRPC,           &
     &                 LUMMMO, LUCCEF, LUCCPO, CC_TAPLAST
       COMMON /CC_TAP/ LUIAJB, LUBFDN, LURES, LUCSOL, LUPRPC,           &
     &                 LUMMMO, LUCCEF, LUCCPO
!
       COMMON /CC_TAP/ CC_TAPLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
