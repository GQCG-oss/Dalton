C
C$Id: cbiwlk.h,v 1.1.1.1 2001-02-08 13:33:27 hjj Exp $
C
      LOGICAL WFPRED, REJECT, KEEPSY, VIBCNV,
     *        START, DOREPW, IMAGE, STRICT, NOORTH, NATCON, V3CAL,
     *        VIBAVE, NMODIF, ECKART
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
     &                ECKART
