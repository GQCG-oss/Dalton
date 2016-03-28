! FILE: ccsdinp.h
      LOGICAL SKIP,   DIRECT, DIRGRD, CCRSTR,
     &        FROIMP, FROEXP, NOCCIT,
     &        CCSAVE, STOLD,  JACEXP, LHTR,
     &        DEBUG, CCSTST, ANAAOD,
     &        HERDIR, FREEZE, KEEPAOIN, NOEONL, NOSORT,
     &        HURWITZ_CHECK,
     &        SIRSOP, LVVVV, ONLYMO, AOSOPPA,
     &        CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &        MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &        CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &        rCCD, drCCD, SOSEX, rTCCD,
     &        CCSDT,CCR12, MTRIP, CHOPT, DCPT2,
     &        ETACCPT, DIRKAPB, MLCC3, MLCCSDPT

      INTEGER MXDIIS, MXLRV,
     &        IT2UPD, IT2START,
     &        ICHANG, IPRINT, KEEPAOTWO, CCSDGNINPlast

      COMMON /CCSDGNINP/ SKIP, DIRECT, DIRGRD, CCRSTR,
     &                   FROIMP,  FROEXP, NOCCIT,
     &                   CCSAVE,  STOLD,  JACEXP, LHTR,
     &                   DEBUG,   CCSTST, ANAAOD,
     &                   MXDIIS,  MXLRV,   IT2UPD, IT2START,
     &                   ICHANG,  IPRINT,  KEEPAOTWO, HERDIR,
     &                   ETACCPT, DIRKAPB,
     &                   FREEZE,  KEEPAOIN, NOEONL, NOSORT,
     &                   HURWITZ_CHECK, MLCC3, MLCCSDPT,
     &                   SIRSOP,  LVVVV, ONLYMO, AOSOPPA,
     &   CCSDGNINPlast !  Very important:
      !  Always keep CCSDGNINPlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

      INTEGER CCMODELSlast
      COMMON /CCMODELS/ CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &                  MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &                  CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &                  rCCD, drCCD, SOSEX, rTCCD,
     &                  CCSDT,CCR12,MTRIP, CHOPT, DCPT2,
     &   CCMODELSlast !  Very important:
      !  Always keep CCMODELSlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

C For the RPA models we can calculate several energies using the
C same set of amplitudes -- this temporary variable is used to store
C one of the energies, the usual variables are reserved for use with
c the choice that could be used in subsequent response or gradient
C calculations
C
      INTEGER ETMPlast
      DOUBLE PRECISION  ETMP
      COMMON/ETMP/ ETMP,
     &   ETMPlast !  Very important:
      !  Always keep ETMPlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

! -- end of ccsdinp.h
