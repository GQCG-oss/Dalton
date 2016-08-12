      INTEGER MXCALL, ICLFRS, ICLLST, ICLBCH, ICLDST, 
     &        DISTCLlast


      COMMON /DISTCL/ MXCALL,
     &                ICLFRS(MXSHEL), ICLLST(MXSHEL),
     &                ICLBCH(MXSHEL), ICLDST(MXSHEL),
     &   DISTCLlast !  Very important:
      !  Always keep DISTCLlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
