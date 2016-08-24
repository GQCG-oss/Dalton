      INTEGER ERISELLAST
      COMMON /ERISEL/ NSELCT(4), NACTAO(MXPRIM,4)
!
!
      COMMON /ERISEL/ ERISELLAST
      !  Very important !!!
      !  Always keep ERISELLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
