C
C     FILE: cbiher.h
C
C     Purpose: define control variables for Hermit
C
      REAL*8  EXPKR,  THRESH
      INTEGER IPRDEF, IORCAR, IORSPH, NPQUAD, NPATOM, IPATOM
      LOGICAL SKIP,   TSTINP,
     &        HAMILT, SUPMAT, SPNORB, DAR2EL, ORBORB,
     &        ONEPRP, PROPRI, ALLATM, TRIANG, NOTWO,  NO2SO,  SOTEST,
     &        DIRAC_HER,      TEST_GEN1INT,   TWOBYTEPACKING

C     Selection of specific one-electron (property) integrals:
      LOGICAL KINENE, KINADI, NUCPOT,
     &        DIPLEN, DIPVEL, QUADRU, THETA,  SECMOM, THRMOM,
     &        CARMOM, SPHMOM,
     &        FERMI,  PSO, SPIDIP, DSO, NMRISS, SDFC, HDO,
     &        S1MAG,  S2MAG,  ANGMOM, ANGLON, LONMOM, MAGMOM,
     &        S1MAGT, MGMOMT, S2MAGT, DSUSNL, DSUSLL, DSUSLH,
     &        DIASUS, DSUTST, NUCSNL, NUCSLO, NUCSHI, NSNLTS, NSLTST,
     &        NELFLD, NSTTST, EFGCAR, EFGSPH, S1MAGL, S1MAGR,
     &        HDOBR,  S1MLT,  S1MRT,  HDOBRT, NPOTST,MGMO2T,HBDO,SUSCGO,
     &        NSTCGO, EXPIKR, DARWIN, MASSVL, CM1,    CM2,
     &        SQHDOL, SQHDOR, SQHD2O, S1ELE,  S1ELB,  ONEELD,
     &        DERHAM, DEROVL, SOMM,   SOFLD,  ROTSTR, OCTGRA,
     &        QUAGRA, DPLGRA, ELGDIL, ELGDIA, MNF_SO, DPTOVL, DPTPOT,
     &        XDDXR3, PVIOLA, PVPINT, POTENE, QDBINT, QDBTST,
     &        RANGMO, RPSO, PXPINT, OZKE, PSOKE, DNSKE, SDKE, FCKE,
     &        DSOKE,  PSOOZ,  EFBDER, EFB2DR, MAGQDP, MQDPTS,
     &        DERAM,  DIPANH, RMAOTWO,LFDIPLN,
     &        S2MBRA, S2MKET, S2MMIX, DOLRINTS
C
      COMMON /CBIHER/ EXPKR(3), THRESH,                                         ! real*8
     &        IPRDEF, IORCAR, IORSPH, NPQUAD, NPATOM, IPATOM(MXCENT),           ! integer
     &        SKIP,   TSTINP,                                                   ! logical control variables
     &        HAMILT, SUPMAT, SPNORB, DAR2EL, ORBORB,
     &        ONEPRP, PROPRI, ALLATM, TRIANG, NOTWO,  NO2SO,  SOTEST,
     &        DIRAC_HER,      TEST_GEN1INT,   TWOBYTEPACKING,
     &        KINENE, KINADI, NUCPOT,                                           ! logicals for specific one-electron integrals
     &        DIPLEN, DIPVEL, QUADRU, THETA,  SECMOM, THRMOM,
     &        CARMOM, SPHMOM,
     &        FERMI,  PSO, SPIDIP, DSO, NMRISS, SDFC, HDO,
     &        S1MAG,  S2MAG,  ANGMOM, ANGLON, LONMOM, MAGMOM,
     &        S1MAGT, MGMOMT, S2MAGT, DSUSNL, DSUSLL, DSUSLH,
     &        DIASUS, DSUTST, NUCSNL, NUCSLO, NUCSHI, NSNLTS, NSLTST,
     &        NELFLD, NSTTST, EFGCAR, EFGSPH, S1MAGL, S1MAGR,
     &        HDOBR,  S1MLT,  S1MRT,  HDOBRT, NPOTST,MGMO2T,HBDO,SUSCGO,
     &        NSTCGO, EXPIKR, DARWIN, MASSVL, CM1,    CM2,
     &        SQHDOL, SQHDOR, SQHD2O, S1ELE,  S1ELB,  ONEELD,
     &        DERHAM, DEROVL, SOMM,   SOFLD,  ROTSTR, OCTGRA,
     &        QUAGRA, DPLGRA, ELGDIL, ELGDIA, MNF_SO, DPTOVL, DPTPOT,
     &        XDDXR3, PVIOLA, PVPINT, POTENE, QDBINT, QDBTST,
     &        RANGMO, RPSO, PXPINT, OZKE, PSOKE, DNSKE, SDKE, FCKE,
     &        DSOKE,  PSOOZ,  EFBDER, EFB2DR, MAGQDP, MQDPTS,
     &        DERAM,  DIPANH, RMAOTWO,LFDIPLN,
     &        S2MBRA, S2MKET, S2MMIX, DOLRINTS
C --- end of cbiher.h ---
