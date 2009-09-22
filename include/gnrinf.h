C File : gnrinf.h
C
c     -*- mode: fortran; fortran-continuation-string: "&" -*-
c     File: gnrinf.h -- general information for DALTON
c
      LOGICAL TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,
     &        HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,
     &        WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,
     &        DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,
     &        TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,
     &        ERFEXP, DOSRIN, SRINTS, CHI1ST, QM3, QMMM
      COMMON /GNRINF/ GRADML, PANAS,  CHIVAL, TESTIN, OPTWLK, RNHERM,
     &                RNSIRI, RNABAC, GEOCNV, HRINPC, SRINPC, RDINPC, 
     &                RDMLIN, PARCAL, DIRCAL, KCHARG, WRINDX, 
     &                WLKREJ, WALKIN, RNRESP, USRIPR, ITERNR,
     &                ITERMX, IPRUSR, SEGBAS, DOCCSD, OPTNEW, NEWSYM,
     &                LENBAS, NEWBAS, NEWPRP, RELCAL, TOTSYM, NMWALK,
     &                DKTRAN, GEOALL, WESTA,  SEGAUX,
     &                ERFEXP, DOSRIN, SRINTS, CHI1ST, QM3, QMMM

      INTEGER LBASDIR
      PARAMETER (LBASDIR = 600)
      CHARACTER*(LBASDIR) BASDIR
      COMMON /GNRCHR/ BASDIR
