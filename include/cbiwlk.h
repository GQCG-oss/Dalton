C
C$Id: cbiwlk.h,v 1.5 2001-12-13 08:54:39 ruden Exp $
C
      PARAMETER (MAXTMP = 20)
      LOGICAL WFPRED, REJECT, KEEPSY, VIBCNV,
     *        START, DOREPW, IMAGE, STRICT, NOORTH, NATCON, V3CAL,
     *        VIBAVE, NMODIF, ECKART, DOTEMP, DOCENT, REUSED
      COMMON /CBIWLK/ TOLST, TRUSTR, TRUSTI, TRUSTD,
     &                REJMIN, REJMAX, RTRMIN, RTRGOD,
     &                XMXNUC, ZERGRD, DISPLC,
     &                STRMOM(3*MXCENT), ECKGEO(3,MXCOOR),
     &                SCALCO(3,MXCENT), THRLDP, ANHFAC, TRUMAX,
     &                ISTMOM(3*MXCENT), NSTMOM,
     &                IPRWLK, IWKTYP, IWKIND, IMODE, ISCTYP,
     &                IPART(MXCENT), DOREPW(0:7),
     &                NZEROG, IZEROG(MXCOOR), START,
     &                WFPRED, REJECT, KEEPSY, VIBCNV,
     &                IMAGE, ISOTPS(MXCENT), STRICT,
     &                NOORTH, NATCON, V3CAL, VIBAVE, NMODIF,
     &                ECKART, DOTEMP, DOCENT, TEMP(MAXTMP), NTEMP,
     &                REUSED
