      INTEGER JTOP, NRTOP, HERTOPlast
C
      COMMON /HERTOP/ JTOP, NRTOP,
     &   HERTOPlast !  Very important:
      !  Always keep HERTOPlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
