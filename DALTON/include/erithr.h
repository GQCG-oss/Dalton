      INTEGER ERITHRlast

      COMMON /ERITHR/ THRSH

      COMMON /ERITHR/ ERITHRLAST
      !  Very important !!!
      !  Always keep ERITHRLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
