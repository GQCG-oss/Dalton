      INTEGER I2BST,IAODPK
      INTEGER SYMSQlast
      COMMON /SYMSQ/ I2BST(8),IAODPK(8,8),                              &
     &   SYMSQlast !  Very important:
      !  Always keep SYMSQlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

