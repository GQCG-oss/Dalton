C
C$Id: cbiwlk.h,v 1.2 2001-02-23 15:34:25 ruud Exp $
C
      PARAMETER (MAXTMP = 10)
      LOGICAL WFPRED, REJECT, KEEPSY, VIBCNV,
     *        START, DOREPW, IMAGE, STRICT, NOORTH, NATCON, V3CAL,
     *        VIBAVE, NMODIF, ECKART, DOTEMP, DOCENT
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
     &                ECKART, DOTEMP, DOCENT, TEMP(MAXTMP), NTEMP
