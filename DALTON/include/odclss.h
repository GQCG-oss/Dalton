      LOGICAL ODTRI1, ODTRI2, ODTR12


      INTEGER ODCLSSLAST

      COMMON /ODCLSS/ NITCL,  NODCL1, NODCL2,                           &
     &                NITBC,  NRTBC,  NODBC1, NODBC2,                   &
     &                NITPP,  NRTPP,  NODPP1, NODPP2,                   &
     &                ODTRI1, ODTRI2, ODTR12,                           &
     &                IODAB,  IODCD,                                    &
     &                NITTR

      COMMON /ODCLSS/ ODCLSSLAST
      !  Very important !!!
      !  Always keep ODCLSSLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
