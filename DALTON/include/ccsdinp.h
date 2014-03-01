      LOGICAL SKIP,   DIRECT, DIRGRD, CCRSTR,
     &        FROIMP, FROEXP, NOCCIT,
     &        CCSAVE, STOLD,  JACEXP, LHTR,
     &        DEBUG, CCSTST, ANAAOD,
     &        HERDIR, FREEZE, KEEPAOIN, NOEONL, NOSORT,
     &        HURWITZ_CHECK,
CSPAS 15.11.2009 implementing AO-SOPPA
C    &        SIRSOP, LVVVV, ONLYMO
     &        SIRSOP, LVVVV, ONLYMO, AOSOPPA
CKeinSPASmehr
C


      LOGICAL ETACCPT, DIRKAPB

      INTEGER MXDIIS, MXLRV,
     &        IT2UPD, IT2START,
     &        ICHANG, IPRINT, KEEPAOTWO

      LOGICAL CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &        MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &        CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &        rCCD, drCCD, SOSEX, rTCCD,
     &        CCSDT,CCR12, MTRIP, CHOPT, DCPT2

      COMMON /CCSDGNINP/ SKIP,   DIRECT, DIRGRD, CCRSTR,
     &                   FROIMP, FROEXP, NOCCIT,
     &                   CCSAVE, STOLD,  JACEXP, LHTR,
     &                   DEBUG,  CCSTST, ANAAOD,
     &                   MXDIIS, MXLRV,IT2UPD, IT2START,
     &                   ICHANG, IPRINT, KEEPAOTWO, HERDIR,
     &                   ETACCPT, DIRKAPB,
     &                   FREEZE, KEEPAOIN, NOEONL, NOSORT,
     &                   HURWITZ_CHECK,
CSPAS 15.11.2009 implementing AO-SOPPA
C    &                   SIRSOP, LVVVV, ONLYMO
     &                   SIRSOP, LVVVV, ONLYMO, AOSOPPA
CKeinSPASmehr

      COMMON /CCMODELS/ CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &                  MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &                  CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &                  rCCD, drCCD, SOSEX, rTCCD,
     &                  CCSDT,CCR12,MTRIP, CHOPT, DCPT2

C For the RPA models we can calcualte several energies using the
C same set of amplitudes -- this temporary variable is used to store
C one of the energies, the usual variables are reserved for use with
C the choice that could be used in subsequent response or gradient
C calculations
      DOUBLE PRECISION  ETMP
      COMMON/ETMP/ETMP

