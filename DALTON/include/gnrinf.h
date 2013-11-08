C File : gnrinf.h
C
c     -*- mode: fortran; fortran-continuation-string: "&" -*-
c     File: gnrinf.h -- general information for DALTON
c
C
C     EMBEDDING : QM part is embedded in environment (solvent or e.g. protein)
C                 May 2011/hjaaj: EMBEDDING = FLAG(16) .or. PCM .or. QM3 .or. QMMM .or. PEQM
C                 (For now, EMBEDDING is defined in sirius/sirinp.F because this is
C                 the first instance where all of FLAG(16), PCM, QM3, QMMM, PEQM are set)
C
      LOGICAL TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,
     &        HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,
     &        WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,
     &        DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,
     &        TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,
     &        ERFEXP, DOSRIN, SRINTS, CHI1ST, DKHINT,
     &        EMBEDDING, QM3, QMMM,   QMNPMM, PEQM
      REAL*8  GRADML, PANAS,  CHIVAL, THR_REDFAC
      INTEGER KCHARG, ITERNR, ITERMX, IPRUSR, LENBAS
C
      COMMON /GNRINF/
     &        ! real*8ZZ
     &        GRADML, PANAS,  CHIVAL, THR_REDFAC,
     &        ! integer:
     &        KCHARG, ITERNR, ITERMX, IPRUSR, LENBAS,
     &        ! logical:
     &        TESTIN, OPTWLK, RNHERM, RNSIRI, RNABAC, GEOCNV,
     &        HRINPC, SRINPC, RDINPC, RDMLIN, PARCAL, DIRCAL,
     &        WRINDX, WLKREJ, WALKIN, RNRESP, USRIPR, SEGBAS,
     &        DOCCSD, OPTNEW, NEWSYM, NEWBAS, NEWPRP, RELCAL,
     &        TOTSYM, NMWALK, DKTRAN, GEOALL, WESTA,  SEGAUX,
     &        ERFEXP, DOSRIN, SRINTS, CHI1ST, DKHINT,
     &        EMBEDDING, QM3, QMMM,   QMNPMM, PEQM

      INTEGER LBASDIR
      PARAMETER (LBASDIR = 600)
      CHARACTER*(LBASDIR) BASDIR
      CHARACTER*12        WFTYPE
      COMMON /GNRCHR/ BASDIR, WFTYPE
C --- end of gnrinf.h ---
