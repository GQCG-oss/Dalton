!     eritap.h
      INTEGER         LUAORC,           NBUFX
      INTEGER         ERITAPlast

      COMMON /ERITAP/ LUAORC(0:MXCOOR), NBUFX(0:MXCOOR),
     &   ERITAPlast !  Very important:
      !  Always keep ERITAPlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

! -- end of eritap.h --
