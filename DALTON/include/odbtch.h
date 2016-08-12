      PARAMETER (NPSORT = 9)
      INTEGER NODBCH, NODCLS, NODPPR, ODBTCHlast
      LOGICAL DOSORT

C     NPSORT has been increased from 8 to 9 for basis-set identifiers (WK/UniKA/31-10-2002).
      COMMON /ODBTCH/ NODBCH, NODCLS, NODPPR, DOSORT(NPSORT),
     &   ODBTCHlast !  Very important:
      !  Always keep ODBTCHlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
