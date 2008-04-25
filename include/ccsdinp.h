      LOGICAL SKIP,   DIRECT, DIRGRD, CCRSTR, 
     &        FROIMP, FROEXP, NOCCIT,
     &        CCSAVE, STOLD,  JACEXP, LHTR,
     &        DEBUG,  CCSTST, ANAAOD,
     &        HERDIR, FREEZE, KEEPAOIN, NOEONL,
     &        SIRSOP, LVVVV, ONLYMO

      LOGICAL ETACCPT, DIRKAPB

      INTEGER MXDIIS, MXLRV,
     &        ICHANG, IPRINT, KEEPAOTWO

      LOGICAL CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &        MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &        CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &        CCSDT,CCR12, MTRIP

      COMMON /CCSDGNINP/ SKIP,   DIRECT, DIRGRD, CCRSTR,
     &                   FROIMP, FROEXP, NOCCIT,
     &                   CCSAVE, STOLD,  JACEXP, LHTR,
     &                   DEBUG,  CCSTST, ANAAOD,
     &                   MXDIIS, MXLRV,
     &                   ICHANG, IPRINT, KEEPAOTWO, HERDIR,
     &                   ETACCPT, DIRKAPB,
     &                   FREEZE, KEEPAOIN, NOEONL,
     &                   SIRSOP, LVVVV, ONLYMO

      COMMON /CCMODELS/ CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &                  MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &                  CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &                  CCSDT,CCR12,MTRIP

