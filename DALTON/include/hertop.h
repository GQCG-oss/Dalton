      INTEGER HERTOPLAST
!
      COMMON /HERTOP/ JTOP, NRTOP
!
!
      COMMON /HERTOP/ HERTOPLAST
      !  Very important !!!
      !  Always keep HERTOPLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
