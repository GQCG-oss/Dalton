C
C     Parameters NCBI? must be updated after changes (for parallelization)
C
C     NOTE: New logicals should appear after PROFIL
C           New integers should appear after IANGMO.
C           Reals should appear at the end.
C
      INTEGER NCBII,  NCBIL
      PARAMETER (NCBII = 13,
     &           NCBIL = 22)
      INTEGER IPRERI, IAOBCH, IPRNT1, IPRNT2, IPROD1, IPROD2
      INTEGER LBFINP, MAXDST, NDMAT,  IANGMO, NSPMAX, MAXDSD,
     &        MXBCH,  MXBCH0
      LOGICAL RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI, OFFCNT,
     &        DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT, EXTPRI,
     &        INTPRI, DODIST, INTSKP, DISTST, WRTINT, FCKINT, PROFIL,
     &        NOLOCS, GRDZER, OLDDER, EXPERI, DOERIP, ERITWO, CCRUN,
     &        COMPRS, GENCON, NCLERI
      DIMENSION IANGMO(4)
      COMMON /CBIERI/ RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI,
     &                OFFCNT, DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT,
     &                EXTPRI, INTPRI, DODIST, INTSKP, DISTST,
     &                WRTINT, FCKINT, PROFIL, EXPERI, DOERIP, ERITWO,
     &                NDMAT,  IPROD1, IPROD2, IAOBCH, IPRNT1, IPRNT2,
     &                IPRERI, LBFINP, MAXDST, IANGMO, NSPMAX, NOLOCS,
     &                GRDZER, OLDDER, CCRUN,  COMPRS, MAXDSD, MXBCH,
     &                MXBCH0, GENCON, NCLERI
