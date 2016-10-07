! erimem.h
      LOGICAL MEMOK
      INTEGER ERIMEMlast

      COMMON /ERIMEM/ MEMADD, MODAB, MODCD, MEMOK

      COMMON /ERIMEM/ ERIMEMLAST
      !  Very important !!!
      !  Always keep ERIMEMlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
