      INTEGER ERISELlast
      COMMON /ERISEL/ NSELCT(4), NACTAO(MXPRIM,4),
     &   ERISELlast !  Very important:
      !  Always keep ERISELlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
