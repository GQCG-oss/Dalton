!
!     File: infind.h - orbital index vectors for sirius module + abacus and response;
!                      vectors generally defined in sirius/sirset.F
!
!     IROW      : for indexing of triangular packed matrices
!     ISMO/ISAO : symmetry of each MO/AO
!     ISW       : reorders to inactive-active-secondary ordering from symmetry ordering
!     ISX       : reorders back
!     ICH       : if neg, inactive no.; if pos., active no.; if zero, empty orbital
!     IOBTYP    : define type of each MO with JTFRO, JTINAC, JTACT, or JTSEC
!     NSM       : symmetry of active orbitals
!     IACTYP    : RAS block for each active orbital
!     ISSMO     : super symmetry of each MO (defined in sirius/sirave.F)
!     ISSORD    : reorder array to super symmetry ordering (defined in sirius/sirave.F)
!
      INTEGER LIROW, JTFRO, JTINAC, JTACT, JTSEC, JTFRFR, JTINFR,       &
     &        JTACFR, JTSEFR, JTININ, JTACIN, JTACAC, JTSEIN, JTSEAC,   &
     &        JTSESE

      PARAMETER (LIROW  = 3*MXCORB + 3)
!     Definition of readable codes for orbital types and 2-index orbital types
      PARAMETER (JTFRO  = 1, JTINAC = 2, JTACT  = 3, JTSEC  = 4)
      PARAMETER (JTFRFR =-1, JTINFR =-2, JTACFR =-3, JTSEFR =-4,        &
     &           JTININ = 1, JTACIN = 2, JTACAC = 3, JTSEIN = 4,        &
     &           JTSEAC = 5, JTSESE = 6)
!

      INTEGER IROW, ISMO, ISAO, ISW, ISX, ICH, IOBTYP,                  &
     &        NSM, IACTYP, ISSMO, ISSORD, INFINDLAST

      COMMON /INFIND/ IROW(LIROW),  ISMO(MXCORB),ISAO(MXCORB),          &
     &                ISW(MXCORB),  ISX(MXCORB),                        &
     &                ICH(MXCORB),  IOBTYP(MXCORB),                     &
     &                NSM(MAXASH),  IACTYP(MAXASH),                     &
     &                ISSMO(MXCORB),ISSORD(MXCORB)

      COMMON /INFIND/ INFINDLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
