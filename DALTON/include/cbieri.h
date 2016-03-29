C File : cbieri.h
C
C     Parameters NCBI? must be updated after changes (for parallelization)

      INTEGER NCBII,  NCBIL
      PARAMETER (NCBII = 18, NCBIL = 31)

C     NOTE: First integer MUST be NDMAT,  number of integers is NCBII
C           First logical MUST be RUNERI, number of logicals is NCBIL
C

      INTEGER NDMAT,  IPROD1, IPROD2, IAOBCH, IPRNT1, IPRNT2,
     &        IPRERI, LBFINP, MAXDST, IANGMO, NSPMAX, 
     &        MAXDSD, MXBCH,  MXBCH0, IFITDM

      LOGICAL RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI, OFFCNT,
     &        DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT, EXTPRI,
     &        INTPRI, DODIST, INTSKP, DISTST, WRTINT, FCKINT, PROFIL,
     &        NOLOCS, GRDZER, OLDDER, EXPERI, DOERIP, ERITWO, CCRUN,
     &        COMPRS, GENCON_ERI, NCLERI, DGABAB

      INTEGER CBIERIlast

      COMMON /CBIERI/ NDMAT,  IPROD1, IPROD2, IAOBCH, IPRNT1, IPRNT2,             ! integers
     &                IPRERI, LBFINP, MAXDST, IANGMO(4), NSPMAX, 
     &                MAXDSD, MXBCH,  MXBCH0, IFITDM,
     &                RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI,             ! logicals
     &                OFFCNT, DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT,
     &                EXTPRI, INTPRI, DODIST, INTSKP, DISTST, WRTINT,
     &                FCKINT, PROFIL, NOLOCS, GRDZER, OLDDER, EXPERI,
     &                DOERIP, ERITWO, CCRUN, COMPRS, GENCON_ERI, NCLERI,
     &                DGABAB,
     &   CBIERIlast !  Very important:
      !  Always keep CBIERIlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
C -- end of cbieri.h --
