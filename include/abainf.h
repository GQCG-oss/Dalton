C
C     File: abainf.h
C     Purpose: Control of what to do in ABACUS module
C
C     NOTE: DOSYM(NSYML) is sent to slaves in a parallel calculation and
C           should not be move from this common block
C     Parameter NSYML must be updated after changes (for parallelization)
C
      INTEGER NSYML, IPRDEF, NWNABA
      PARAMETER (NSYML = 8)
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,
     &        VIB,    RESTAR, DOWALK, GDALL,  CCSD,
     &        H2MO,   DOSYM,  DOLRES, DOEXCI, SHIELD,
     &        SPNSPN, MAGSUS, VCD,    NACME,  AAT,
     &        NOLOND, FCKDDR, ECD,    NODIFC, DODRCT,
     &        SUPMAT, MOLGFA, OPTROT, SPINRO, MASSVE,
     &        DARWIN, ABALNR, VROA,   NOCMC,  EXPFCK,
     &        RAMAN,  QUADRU, NQCC,   HYPER,  VERDET,
     &        MCD,    HELFEY, LINCPL, ABASOP, SKIPAB,
     &        ABA_ALPHA, EXPGRD, CTOCD, NUMHES, DOD2DQ2,
     &        OECD, MVEOR
      COMMON /ABAINF/ IPRDEF, NWNABA, 
     &                MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,
     &                VIB,    RESTAR, DOWALK, GDALL,  CCSD,
     &                H2MO,   DOSYM(NSYML),   DOLRES, DOEXCI, SHIELD,
     &                SPNSPN, MAGSUS, VCD,    NACME,  AAT,
     &                NOLOND, FCKDDR, ECD,    NODIFC, DODRCT,
     &                SUPMAT, MOLGFA, OPTROT, SPINRO, MASSVE,
     &                DARWIN, ABALNR, VROA,   NOCMC,  EXPFCK,
     &                RAMAN,  QUADRU, NQCC,   HYPER,  VERDET,
     &                MCD,    HELFEY, LINCPL, ABASOP, SKIPAB,
     &                ABA_ALPHA,      EXPGRD, CTOCD,  NUMHES, DOD2DQ2,
     &                OECD, MVEOR
