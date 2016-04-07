!     eribuf.h
      INTEGER         LBUF, NBUFS, NIBUF, NBITS, IBIT1, IBIT2
      LOGICAL         NEWDIS
      INTEGER         ERIBUFlast

      COMMON /ERIBUF/ LBUF, NBUFS, NIBUF, NBITS, IBIT1, IBIT2,
     &                NEWDIS,
     &   ERIBUFlast !  Very important:
      !  Always keep ERIBUFlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

! -- end of eribuf.h --
