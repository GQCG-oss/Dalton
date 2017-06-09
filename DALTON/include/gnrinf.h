!
!     -*- mode: fortran; fortran-continuation-string: "&" -*-
!     File: gnrinf.h -- general information for DALTON
!
!
!     EMBEDDING : QM part is embedded in environment (solvent or e.g. protein)
!                 May 2011/hjaaj: EMBEDDING = FLAG(16) .or. PCM .or. QM3 .or. QMMM
!                 (For now, EMBEDDING is defined in sirius/sirinp.F because this is
!                 the first instance where all of FLAG(16), PCM, QM3, QMMM are set)
!
      LOGICAL TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,           &
     &        HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,           &
     &        WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,           &
     &        DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,           &
     &        TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,           &
     &        ERFEXP, DOSRIN, SRINTS, CHI1ST, DKHINT,                   &
     &        EMBEDDING, QM3, QMMM,   QMNPMM, QFIT,   USE_LSLIB,        &
     &        USE_OPENRSP, SIR_INPPRC,NEWGEO, DOFDE
      REAL*8  GRADML, PANAS,  CHIVAL, THR_REDFAC
      INTEGER KCHARG, ITERNR, ITERMX, IPRUSR, LENBAS
      INTEGER GNRINFlast
!
      COMMON /GNRINF/                                                   &
              ! real*8ZZ
     &        GRADML, PANAS,  CHIVAL, THR_REDFAC,                       &
              ! integer
     &        KCHARG, ITERNR, ITERMX, IPRUSR, LENBAS,                   &
              ! logical
     &        TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,           &
     &        HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,           &
     &        WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,           &
     &        DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,           &
     &        TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,           &
     &        ERFEXP, DOSRIN, SRINTS, CHI1ST, DKHINT,                   &
     &        EMBEDDING, QM3, QMMM,   QMNPMM, QFIT,   USE_LSLIB,        &
     &        USE_OPENRSP, SIR_INPPRC,NEWGEO, DOFDE
      COMMON /GNRINF/ GNRINFlast
      ! Very important: Always keep GNRINFlast as the last variable in the common block. 
      ! See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

      INTEGER LBASDIR
      PARAMETER (LBASDIR = 600)
      CHARACTER*(LBASDIR) BASDIR
      CHARACTER*12        WFTYPE
      COMMON /GNRCHR/ BASDIR, WFTYPE

! -- end of gnrinf.h
