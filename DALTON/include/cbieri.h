! File : cbieri.h
!
!     Parameters NCBI? must be updated after changes (for parallelization)
!
!     NOTE: New logicals should appear after the last logical (NCLERI)
!           New integers should appear after the last integer (IFITDM)
!           Reals should appear at the end.
!
      INTEGER NCBII,  NCBIL
!     
!     IANGMO takes up space for 4 ints. That is not included in this counter
      PARAMETER (NCBII = 15,                                            &
     &           NCBIL = 31)

      INTEGER IPRERI, IAOBCH, IPRNT1, IPRNT2, IPROD1, IPROD2

      INTEGER LBFINP, MAXDST, NDMAT,  IANGMO, NSPMAX, MAXDSD,           &
     &        MXBCH,  MXBCH0, IFITDM, CBIERILAST

      LOGICAL RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI, OFFCNT,   &
     &        DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT, EXTPRI,           &
     &        INTPRI, DODIST, INTSKP, DISTST, WRTINT, FCKINT, PROFIL,   &
     &        NOLOCS, GRDZER, OLDDER, EXPERI, DOERIP, ERITWO, CCRUN,    &
     &        COMPRS, GENCON_ERI, NCLERI, DGABAB  

      DIMENSION IANGMO(4)

      COMMON /CBIERI/ RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI,   &
     &                OFFCNT, DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT,   &
     &                EXTPRI, INTPRI, DODIST, INTSKP, DISTST, WRTINT,   &
     &                FCKINT, PROFIL, NOLOCS, GRDZER, OLDDER, EXPERI,   &
     &                DOERIP, ERITWO, CCRUN, COMPRS, GENCON_ERI, NCLERI,&
     &                DGABAB, NDMAT,  IPROD1, IPROD2, IAOBCH, IPRNT1,   &
     &                IPRNT2, IPRERI, LBFINP, MAXDST, IANGMO, NSPMAX,   &
     &                MAXDSD, MXBCH,  MXBCH0, IFITDM

      COMMON /CBIERI/ CBIERILAST
      !  Very important !!!
      !  Always keep CBIERILAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
! -- end of cbieri.h --
