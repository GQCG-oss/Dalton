!
!------------------------------------------------------------
! information for the packing of the AO integrals on the
! presorted CCAOIN_* files
!------------------------------------------------------------
!
      LOGICAL LPACKINT
      INTEGER IPCKTABINT(0:255)
      INTEGER IOFFINT(MXCORB)
      INTEGER NPCKINT(MXCORB)
      INTEGER NTOTINT, NTOTPCK
      INTEGER CCPACKLAST
!
!
!
#if defined (SYS_CRAY)
      REAL THRPCKINT, PCKRATIO, PCKTIME
#else
      DOUBLE PRECISION THRPCKINT, PCKRATIO, PCKTIME
#endif
      COMMON /CCPACK/ THRPCKINT, PCKRATIO, PCKTIME,                     &
     &                IPCKTABINT, IOFFINT, NPCKINT,                     &
     &                NTOTINT, NTOTPCK, LPACKINT
!
      COMMON /CCPACK/ CCPACKLAST
      !  Very important !!!
      !  Always keep CCPACKLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
