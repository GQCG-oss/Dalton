      INTEGER MAXSIM

      PARAMETER(MAXSIM = 540)

      LOGICAL DUMPCD

      INTEGER IT2DEL, IT2DLR, CCSDIOLAST


      COMMON /CCSDIO/ IT2DEL(MXCORB),                                   &
     &                DUMPCD,IT2DLR(MXCORB,MAXSIM)


      COMMON /CCSDIO/ CCSDIOLAST
      !  Very important !!!
      !  Always keep CCSDIOLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
