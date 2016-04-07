      INTEGER ERITHRlast

      COMMON /ERITHR/ THRSH,
     &   ERITHRlast !  Very important:
      !  Always keep ERITHRlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
