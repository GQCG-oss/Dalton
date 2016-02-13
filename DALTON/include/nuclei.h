!
! FILE: include/nuclei.h
!
      REAL*8  CHARGE, CORD, GNUEXP

      INTEGER NUCPRE, NUCNUM, NUCDEG, ISTBNU, NCTOT,                    &
     &        NUCIND, NUCDEP, NTRACO, ITRACO, NATOMS, NFLOAT,           &
     &        NBASIS, NLARGE, NSMALL, NPBAS,  NPLRG,  NPSML,            &
     &        NCHTOT, INCENT, INUNIQ, NDEGNM, ISOTOP, IZATOM,           &
     &        NBASISAUX, NPBASAUX,    NAUX,   NPAUX,  MULBSI,           &
     &        NUCLEILAST, NUCLECLAST
!     MULBSI has been added for multiple basis sets (WK/UniKA/31-10-2002).
      LOGICAL NOORBT, GAUNUC
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(3,MXCENT), GNUEXP(MXCENT),   &
     &                NUCPRE(MXCENT), NUCNUM(MXCENT,8), NUCDEG(MXCENT), &
     &                ISTBNU(MXCENT), NDEGNM(MXCENT), NUCIND, NUCDEP,   &
     &                NTRACO, ITRACO(3),                                &
     &                NATOMS, NFLOAT, NBASIS, NLARGE, NSMALL, NPBAS,    &
     &                NPLRG, NPSML, NCHTOT, INCENT(MXCENT), NCTOT,      &
     &                INUNIQ(MXCENT), ISOTOP(MXCENT),IZATOM(MXCENT),    &
     &                NBASISAUX, NPBASAUX, NAUX, NPAUX,                 &
     &                MULBSI(MXCENT),                                   &
     &                GAUNUC, NOORBT(MXCENT)
      COMMON /NUCLEI/ NUCLEILAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. (see below)
!
      CHARACTER       NAMN*4, NAMEX*6, NAMDEP*6, NAMDPX*8
      COMMON /NUCLEC/ NAMN(MXCENT),   NAMEX(MXCOOR),                    &
     &                NAMDEP(MXCENT), NAMDPX(MXCOOR)
      COMMON /NUCLEC/ NUCLECLAST
      !   Very important !!!
      !   Always keep this variable as the last variable in the common block. 
      !   If you add more variables to the block add them before <name>last.
      !   This variable is used to synchronize slaves for parallel
      !   calculations. Other than acting as a target to calculate the size of a common
      !   block, they have no use.
      !   Use CALL GETBYTESPAN(firstvar, <name>last, SizeInBytes) from all processes 
      !   to get the number of bytes needed to transfer the common block.
      !   Then transfer the block with mpi_bcast(firstvar, SizeInBytes, mpi_byte, 0, mpi_comm_world, ierr)
!  -- end of nuclei.h --
