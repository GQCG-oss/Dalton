      PARAMETER (MXDIST = 250)

      INTEGER ERIDSTlast

      COMMON /ERIDST/ NDISTR, MLTDST, KHKDST, NACDST, IRPDST,
     &                MAXCML, NTCLAS, NDST,   LPNBCH,
     &                INDDST(MXDIST), IACDST(MXDIST),
     &                INDXDS(MXSHEL), NCLASS(MXSHEL), MCLASS(MXSHEL),
     &                KCLASS(MXSHEL), KLAOBT(MXSHEL), LUINTD(MXCOOR),
     &                NCCDST(MXDIST), NCCFST(MXDIST),
     &   ERIDSTlast !  Very important:
      !  Always keep ERIDSTlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
