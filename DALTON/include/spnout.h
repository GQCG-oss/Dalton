! -- FILE: spnout.h --
      LOGICAL DOSD, DODSO, DOFC, DOSDFC, DOPSO, DOSELE, ANISON,
     &        FCFIN, SPNISO, NCSPNI,
     &        SOS, SOSSPN, SOSOCC, SOSOCS
      REAL*8 ABUND
      INTEGER ISPPRI, ISOTPS, NSTATS, NSTATT, NSTATI, NSTATF,
     &        NITRST, NUCSPI
      COMMON /SPNOUT/ ABUND,                                        ! real*8
     &        ISPPRI, ISOTPS(MXCENT),                               ! integer
     &        NSTATS, NSTATT, NSTATI, NSTATF, NITRST, NUCSPI,
     &        DOSD, DODSO, DOFC, DOSDFC, DOPSO, DOSELE, ANISON,     ! logical
     &        FCFIN, SPNISO, NCSPNI(MXCENT),
     &        SOS, SOSSPN, SOSOCC, SOSOCS
! -- end of spnout.h --
