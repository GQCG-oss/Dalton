      INTEGER ISAO
      INTEGER CCISAOLAST
!
      COMMON /CCISAO/ ISAO(MXCORB)
!
      COMMON /CCISAO/ CCISAOLAST
      !  Very important !!!
      !  Always keep CCISAOLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
