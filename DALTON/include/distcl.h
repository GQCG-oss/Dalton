      INTEGER MXCALL, ICLFRS, ICLLST, ICLBCH, ICLDST,                   &
     &        DISTCLLAST                                                


      COMMON /DISTCL/ MXCALL,                                           &
     &                ICLFRS(MXSHEL), ICLLST(MXSHEL),                   &
     &                ICLBCH(MXSHEL), ICLDST(MXSHEL)                    


      COMMON /DISTCL/ DISTCLLAST
      !  Very important !!!
      !  Always keep DISTCLLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
