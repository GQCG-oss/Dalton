!     eribuf.h
      INTEGER         LBUF, NBUFS, NIBUF, NBITS, IBIT1, IBIT2
      LOGICAL         NEWDIS
      INTEGER         ERIBUFLAST

      COMMON /ERIBUF/ LBUF, NBUFS, NIBUF, NBITS, IBIT1, IBIT2,          &
     &                NEWDIS,                                           &
     &                ERIBUFLAST
      !  Very important !!!
      !  Always keep ERIBUFLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

! -- end of eribuf.h --
