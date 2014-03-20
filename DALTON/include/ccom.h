!     File : ccom.h -- permanent info about basis functions
!         (l = J + 1, i.e. e.g. p-orbital is J = 2)
!         NHTYP - maximum angular quantum number + 1 for ALL orbitals
!         KHK(J) - number of spherical (cartesian) components for given J
!         KCK(J) - number of Cartesian components for given J
!         NHKOFF(J) - offset for components in list of l-functions
!         GTOTYP(index(K,J)) - label for K component of J-type orbital
!         DOCART - Cartesian basis fu., if false:  spherical or own def.
!         SPH(J) - true if cartesian to spherical/own basis needed
!         SPHNRM - true if all basis functions are normalized
!                  (only true for spherical, not for cartesian or own)
!
!     NHTYP  is set in BASINP(herrdn.F),
!     KHK, KCK, SPH are set in SPHINP(herrdn.F), based on user input
!     NHKOFF is set in BASPAR(herrdn.F)
!     DOCART are set in line 4 from .mol READI1(herrdn.F)
!     SPHNRM true if spherical basis funcitons (herrdn.F)
!     GTOTYP is set in CARLAB, SPHLAB or input (herrdn.F)


      REAL*8  THRS
      INTEGER NHTYP,  KHK(MXQN), KCK(MXQN), NHKOFF(MXQN), CCOMLAST
      LOGICAL DOCART, SPH(MXQN), SPHNRM
      COMMON /CCOM/ THRS,                                               &
     &              NHTYP, KHK, KCK, NHKOFF,                            &
     &              DOCART, SPH, SPHNRM
      COMMON /CCOM/ CCOMLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
      INTEGER CCOMCLAST
      CHARACTER*4 GTOTYP(MXQN*(MXQN+1)*(MXQN+2)/6)
      COMMON /CCOMC/ GTOTYP
      COMMON /CCOMC/ CCOMCLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
!     --- end of ccom.h ---
