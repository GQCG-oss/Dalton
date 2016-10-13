      PARAMETER (NPSORT = 9)
      LOGICAL DOSORT
      INTEGER ODBTCHLAST

      COMMON /ODBTCH/ NODBCH, NODCLS, NODPPR, DOSORT(NPSORT)
!     NPSORT has been increased from 8 to 9 for basis-set identifiers
!     (WK/UniKA/31-10-2002).



      COMMON /ODBTCH/ ODBTCHLAST
      !  Very important !!!
      !  Always keep ODBTCHLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
