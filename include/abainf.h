C
C     Parameters NSYML must be updated after changes (for parallelization)
C     NOTE: DOSYM are sent to slaves in a parallel calculation and
C           should not be move from this subroutine
C
      INTEGER NSYML, IPRDEF, NWNABA, NTCSY, NSCSY
      PARAMETER (NSYML = 8)
cLig <> added CTOCD
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,
     &        VIB,    RESTAR, DOWALK, GDALL,  CCSD,
     &        H2MO,   DOSYM,  DOLRES, DOEXCI, SHIELD,
     &        SPNSPN, MAGSUS, VCD,    NACME,  AAT,
     &        NOLOND, FCKDDR, ECD,    NODIFC, DODRCT,
     &        SUPMAT, MOLGFA, OPTROT, SPINRO, MASSVE,
     &        DARWIN, ABALNR, VROA,   NOCMC,  EXPFCK,
     &        RAMAN,  QUADRU, NQCC,   HYPER,  VERDET,
     &        MCD,    HELFEY, LINCPL, ABASOP, SKIPAB,
     &        ABA_ALPHA, EXPGRD, CTOCD
      COMMON /ABAINF/ NTCSY(8), NSCSY(8), IPRDEF, NWNABA, 
     &                MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,
     &                VIB,    RESTAR, DOWALK, GDALL,  CCSD,
     &                H2MO,   DOSYM(8),DOLRES,DOEXCI, SHIELD,
     &                SPNSPN, MAGSUS, VCD,    NACME,  AAT,
     &                NOLOND, FCKDDR, ECD,    NODIFC, DODRCT,
     &                SUPMAT, MOLGFA, OPTROT, SPINRO, MASSVE,
     &                DARWIN, ABALNR, VROA,   NOCMC,  EXPFCK,
     &                RAMAN,  QUADRU, NQCC,   HYPER,  VERDET,
     &                MCD,    HELFEY, LINCPL, ABASOP, SKIPAB,
     &                ABA_ALPHA, EXPGRD, CTOCD
