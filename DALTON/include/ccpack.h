C
C------------------------------------------------------------
C information for the packing of the AO integrals on the
C presorted CCAOIN_* files
C------------------------------------------------------------
C
      REAL*8  THRPCKINT, PCKRATIO, PCKTIME
      INTEGER IPCKTABINT(0:255)
      INTEGER IOFFINT(MXCORB)
      INTEGER NPCKINT(MXCORB)
      INTEGER NTOTINT, NTOTPCK
      LOGICAL LPACKINT
      INTEGER CCPACKlast

      COMMON /CCPACK/ THRPCKINT, PCKRATIO, PCKTIME,
     &                IPCKTABINT, IOFFINT, NPCKINT,
     &                NTOTINT, NTOTPCK, LPACKINT,
     &   CCPACKlast !  Very important:
      !  Always keep CCPACKlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
! -- end of ccpack.h
