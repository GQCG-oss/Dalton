      INTEGER I2BST,IAODPK
      INTEGER SYMSQLAST
      COMMON /SYMSQ/ I2BST(8),IAODPK(8,8)
!      
      COMMON /SYMSQ/ SYMSQLAST
      !  Very important !!!
      !  Always keep SYMSQLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

