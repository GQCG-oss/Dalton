c     -*- mode: fortran; fortran-continuation-string: "&" -*-
c     File: gnrinf.h -- general information for DALTON
c
      LOGICAL TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,
     &        HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,
     &        WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,
     &        DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,
     &        TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX
      COMMON /GNRINF/ GRADML, PANAS,  TESTIN, OPTWLK, RNHERM,
     &                RNSIRI, RNABAC, GEOCNV, HRINPC, SRINPC, RDINPC, 
     &                RDMLIN, PARCAL, DIRCAL, KCHARG, WRINDX, 
     &                WLKREJ, WALKIN, RNRESP, USRIPR, ITERNR,
     &                ITERMX, IPRUSR, SEGBAS, DOCCSD, OPTNEW, NEWSYM,
     &                LENBAS, NEWBAS, NEWPRP, RELCAL, TOTSYM, NMWALK,
     &                DKTRAN, GEOALL, WESTA,  SEGAUX
      CHARACTER*70 BASDIR
      COMMON /GNRCHR/ BASDIR
