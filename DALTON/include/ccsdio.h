      INTEGER MAXSIM

      PARAMETER(MAXSIM = 540)

      LOGICAL DUMPCD

      INTEGER IT2DEL, IT2DLR, CCSDIOlast


      COMMON /CCSDIO/ IT2DEL(MXCORB),
     &                DUMPCD,IT2DLR(MXCORB,MAXSIM),
     &   CCSDIOlast !  Very important:
      !  Always keep CCSDIOlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
