C
C$Id: abainf.h,v 1.1.1.1 2001-02-08 13:33:28 hjj Exp $
C
C     Parameters NSYML must be updated after changes (for parallelization)
C     NOTE: DOSYM are sent to slaves in a parallel calculation and
C           should not be move from this subroutine
C
      PARAMETER (NSYML = 8)
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,
     &        VIB,    RESTAR, DOWALK, GDALL,  CCSD,
     &        H2MO,   DOSYM,  DOLRES, DOEXCI, SHIELD,
     &        SPNSPN, MAGSUS, VCD,    NACME,  AAT,
     &        NOLOND, FCKDDR, ECD, NODIFC, DODRCT,
     &        SUPMAT, MOLGFA, OPTROT,
     &        SPINRO, MASSVE, DARWIN, ABALNR, VROA, NOCMC,
     &        EXPFCK, RAMAN,  QUADRU, NQCC, HYPER, VERDET, MCD,
     &        HELFEY, LINCPL, ABASOP, SKIPAB
      COMMON /ABAINF/ IPRDEF, NWNABA,
     &                MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,
     &                VIB,    RESTAR, DOWALK, GDALL,  CCSD,
     &                H2MO,   DOSYM(8), DOLRES, DOEXCI,
     &                NTCSY(8), NSCSY(8), SHIELD, SPNSPN,
     &                MAGSUS, VCD, NACME, AAT, NOLOND, FCKDDR, ECD,
     &                DODRCT, SUPMAT, MOLGFA, OPTROT, SPINRO,
     &                MASSVE, DARWIN, ABALNR, VROA,
     &                NODIFC, ISOTOP(MXCENT), NOCMC, EXPFCK,
     &                RAMAN, QUADRU, NQCC, HYPER, VERDET, MCD, HELFEY,
     &                LINCPL, ABASOP, SKIPAB

