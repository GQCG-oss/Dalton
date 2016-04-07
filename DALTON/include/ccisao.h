      INTEGER ISAO
      INTEGER CCISAOlast
C
      COMMON /CCISAO/ ISAO(MXCORB),                                     &
     &   CCISAOlast !  Very important:
      !  Always keep CCISAOlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
! -- end of ccisao.h
