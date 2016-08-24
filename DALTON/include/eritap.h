!     eritap.h
      INTEGER         LUAORC,           NBUFX
      INTEGER         ERITAPLAST

      COMMON /ERITAP/ LUAORC(0:MXCOOR), NBUFX(0:MXCOOR),                &
     &                ERITAPLAST

      !  Very important !!!
      !  Always keep ERITAPLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

! -- end of eritap.h --
