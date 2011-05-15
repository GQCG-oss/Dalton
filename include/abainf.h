!
!     File: abainf.h
!     Purpose: Control of what to do in ABACUS module
!
!     NOTE: DOSYM(NSYML) is sent to slaves in a parallel calculation and
!           should not be moved from this common block !
!
CSPAS:20/3-2011: allowing for perturbations in only one irrep
C     INTEGER NSYML, IPRDEF, NWNABA
      INTEGER NSYML, IPRDEF, NWNABA, IRVIBG
CKeinSPASmehr
      PARAMETER (NSYML = 8)
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,
     &        VIB,    RESTAR, DOWALK, GDALL,  CCSD,
     &        H2MO,   DOSYM(NSYML),   DOLRES, DOEXCI, SHIELD,
     &        SPNSPN, MAGSUS, VCD,    NACME,  AAT,
     &        NOLOND, FCKDDR, ECD,    NODIFC, DODRCT,
     &        SUPMAT, MOLGFA, OPTROT, SPINRO, MASSVE,
     &        DARWIN, ABALNR, VROA,   NOCMC,  EXPFCK,
     &        RAMAN,  QUADRU, NQCC,   HYPER,  VERDET,
     &        MCD,    HELFEY, LINCPL, ABASOP, SKIPAB,
     &        ABA_ALPHA,      EXPGRD, CTOCD,  NUMHES, DOD2DQ2,
     &        OECD,   MVEOR,  QPGRAD, SECNDM, THIRDM, VIB_G,
CSPAS:20/3-2011: allowing for perturbations in only one irrep
C    &        GSDIP,  GSQUAD, GSOCT,  GSDIDI, GSDIQU, GSQUQU
     &        GSDIP,  GSQUAD, GSOCT,  GSDIDI, GSDIQU, GSQUQU,
     &        VIBGIR
CKeinSPASmehr
      COMMON /ABAINF/ IPRDEF, NWNABA, 
     &        MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,
     &        VIB,    RESTAR, DOWALK, GDALL,  CCSD,
     &        H2MO,   DOSYM,          DOLRES, DOEXCI, SHIELD,
     &        SPNSPN, MAGSUS, VCD,    NACME,  AAT,
     &        NOLOND, FCKDDR, ECD,    NODIFC, DODRCT,
     &        SUPMAT, MOLGFA, OPTROT, SPINRO, MASSVE,
     &        DARWIN, ABALNR, VROA,   NOCMC,  EXPFCK,
     &        RAMAN,  QUADRU, NQCC,   HYPER,  VERDET,
     &        MCD,    HELFEY, LINCPL, ABASOP, SKIPAB,
     &        ABA_ALPHA,      EXPGRD, CTOCD,  NUMHES, DOD2DQ2,
     &        OECD,   MVEOR,  QPGRAD, SECNDM, THIRDM, VIB_G,
CSPAS:20/3-2011: allowing for perturbations in only one irrep
C    &        GSDIP,  GSQUAD, GSOCT,  GSDIDI, GSDIQU, GSQUQU
     &        GSDIP,  GSQUAD, GSOCT,  GSDIDI, GSDIQU, GSQUQU,
     &        VIBGIR, IRVIBG
CKeinSPASmehr
! -- end of abainf.h --
