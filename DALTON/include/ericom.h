      LOGICAL DOPTH1, DOPTH2, GDER,   BDER,  UNDIFF,                    &
     &        FIRST,  LAST,   LSTCLS, DIACLS,                           &
     &        SPHRA,  SPHRB,  SPHRC,  SPHRD, SPHRAB, SPHRCD,            &
     &        GCONA,  GCONB,  GCONC,  GCOND, GCONAB, GCONCD,            &
     &        DIAGAB, DIAGCD, DIAGPQ, TKMPAB, TKMPCD,TCMPAB, TCMPCD

      INTEGER ERICOMLAST
!
      COMMON /ERICOM/ SCRMAB, SCRMCD,                                   &
     &                MAXDER, IPATH, NPDIMA, NPDIMB, NCORS, NPERTS,     &
     &                NHKTA,  NHKTB,  NHKTC,  NHKTD,                    &
     &                JMAXA,  JMAXB,  JMAXC,  JMAXD,                    &
     &                JMAXAB, JMAXCD, JMAX,                             &
     &                NTUV,   NTUVAB, NTUVCD,                           &
     &                KHKTA,  KHKTB,  KHKTC,  KHKTD, KHKTAB, KHKTCD,    &
     &                KCKTA,  KCKTB,  KCKTC,  KCKTD, KCKTAB, KCKTCD,    &
     &                ISTBLA, ISTBLB, ISTBLC, ISTBLD,                   &
     &                ISTBLR, ISTBLS, ISTBLT,                           &
     &                MLTPA,  MLTPB,  MLTPC,  MLTPD, MLTPX,             &
     &                MLTPR,  MLTPS,  MLTPT,                            &
     &                NREDZ,  MLTPZ,                                    &
     &                NGTOAB, NGTOCD,                                   &
     &                NHKMAX, KCKMAX, KC2MAX, KCREC1, NRDER,            &
     &                NPRFA,  NPRFB,  NPRFC,  NPRFD,                    &
     &                NPRFAB, NPRFCD, NPRFPQ,                           &
     &                NCTFA,  NCTFB,  NCTFC,  NCTFD,                    &
     &                NCTFAB, NCTFCD, NCTFPQ, NCNTAB, NCNTCD,           &
     &                ICMATA, ICMATB, ICMATC, ICMATD,                   &
     &                NCSQ1,  NCSQ2,                                    &
     &                NODSAB, NODSCD, NODSPQ,                           &
     &                NODTAB, NODTCD, KODSAB, KODSCD, MEMBCH, MAXBCH,   &
     &                NBCHES, NBTPAS, NPQBCX, NPPBCX, NPQBCS, NPPBCS,   &
     &                NPPX,   NPCX,   NCPX,   NCCX, NCCT, NAOINT,       &
     &                NPPS,   NPCS,   NCPS,   NCCS,                     &
     &                NPPAB,  NPPCD,  NITPQ,                            &
     &                IPQ0(3),IPQXYZ, IAB0(3), IABXYZ, ICD0(3), ICDXYZ, &
     &                IHHXYZ, IHCXYZ, ICHXYZ, ICCXYZ,                   &
     &                DOPTH1, DOPTH2, GDER,   BDER,  UNDIFF,            &
     &                FIRST,  LAST,   LSTCLS, DIACLS,                   &
     &                SPHRA,  SPHRB,  SPHRC,  SPHRD, SPHRAB, SPHRCD,    &
     &                GCONA,  GCONB,  GCONC,  GCOND, GCONAB, GCONCD,    &
     &                DIAGAB, DIAGCD, DIAGPQ,                           &
     &                TKMPAB, TKMPCD,TCMPAB, TCMPCD



      COMMON /ERICOM/ ERICOMLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
