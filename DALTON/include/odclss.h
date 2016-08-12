      LOGICAL ODTRI1, ODTRI2, ODTR12

      INTEGER ODCLSSlast

      COMMON /ODCLSS/ NITCL,  NODCL1, NODCL2,
     &                NITBC,  NRTBC,  NODBC1, NODBC2,
     &                NITPP,  NRTPP,  NODPP1, NODPP2,
     &                ODTRI1, ODTRI2, ODTR12,
     &                IODAB,  IODCD, NITTR,
     &   ODCLSSlast !  Very important:
      !  Always keep ODCLSSlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
