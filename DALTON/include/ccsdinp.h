      LOGICAL SKIP,   DIRECT, DIRGRD, CCRSTR,                           &
     &        FROIMP, FROEXP, NOCCIT,                                   &
     &        CCSAVE, STOLD,  JACEXP, LHTR,                             &
     &        DEBUG, CCSTST, ANAAOD,                                    &
     &        HERDIR, FREEZE, KEEPAOIN, NOEONL, NOSORT,                 &
     &        HURWITZ_CHECK,                                            &
     &        SIRSOP, LVVVV, ONLYMO, AOSOPPA,                           &
     &        CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,                     &
     &        MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,                     &
     &        CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,                      &
     &        rCCD, drCCD, SOSEX, rTCCD,                                &
     &        CCSDT,CCR12, MTRIP, CHOPT, DCPT2,                         &
     &        ETACCPT, DIRKAPB, MLCC3, MLCCSDPT

      INTEGER MXDIIS, MXLRV,                                            &
     &        IT2UPD, IT2START,                                         &
     &        ICHANG, IPRINT, KEEPAOTWO, CCSDGNINPLAST

      COMMON /CCSDGNINP/ SKIP, DIRECT, DIRGRD, CCRSTR,                  &
     &                   FROIMP,  FROEXP, NOCCIT,                       &
     &                   CCSAVE,  STOLD,  JACEXP, LHTR,                 &
     &                   DEBUG,   CCSTST, ANAAOD,                       &
     &                   MXDIIS,  MXLRV,   IT2UPD, IT2START,            &
     &                   ICHANG,  IPRINT,  KEEPAOTWO, HERDIR,           &
     &                   ETACCPT, DIRKAPB,                              &
     &                   FREEZE,  KEEPAOIN, NOEONL, NOSORT,             &
     &                   HURWITZ_CHECK, MLCC3, MLCCSDPT,                &
     &                   SIRSOP,  LVVVV, ONLYMO, AOSOPPA


      COMMON /CCSDGNINP/ CCSDGNINPLAST 
      !  Very important !!!
      !  Always keep CCSDGNINPLSDT as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.


      INTEGER CCMODELSLAST
      COMMON /CCMODELS/ CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,           &
     &                  MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,           &
     &                  CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,            &
     &                  rCCD, drCCD, SOSEX, rTCCD,                      &
     &                  CCSDT,CCR12,MTRIP, CHOPT, DCPT2                 
      COMMON /CCMODELS/ CCMODELSLAST 
      !  Very important !!!
      !  Always keep CCMODELSLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.





! For the RPA models we can calcualte several energies using the
! same set of amplitudes -- this temporary variable is used to store
! one of the energies, the usual variables are reserved for use with
! the choice that could be used in subsequent response or gradient
! calculations
!
!
!
!
      INTEGER ETMPLAST
      DOUBLE PRECISION  ETMP
      COMMON/ETMP/ ETMP
!
      COMMON /ETMP/ ETMPLAST 
      !  Very important !!!
      !  Always keep ETMPLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
