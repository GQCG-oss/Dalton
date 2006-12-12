c  -*- mode:fortran; fortran-continuation-string: "&" -*-
C     this common block contains general Sirius input read in sirius/sirinp.F
C     (specified under **WAVE FUNCTION).  /hjaaj Oct 2003
C     - SCF specific input is in scbrhf.h
C     - orbital specifications are in inforb.h
      INTEGER NFLAG, MAXRTS, MXFELT
      PARAMETER (NFLAG = 80, MAXRTS = 20, MXFELT = 20)
      INTEGER         NFIELD, ISPIN,ISTATE,LSYM,NACTEL, MCTYPE,
     *                LSOLMX,NLMSOL,NELMN1,NELMX1,NELMN3,NELMX3,
     *                LROOTS,NROOTS,IROOT,
     *                NOROT        ,IMOORD        ,
     *                IORTO,ICI0,KDEL,ICHECK,NTIT,
     *                MAXMAC,MAXMIC,MAXJT,MAXCIT,MAXUIT,MAXAPM,MAXABS,
     *                ITRLVL,ITRFIN,JCHSYM,JCHORB,
     *                NROOCI,ISTACI, MXCIMA, ICICNO,IMCCNO
      COMMON /INTINP/ NFIELD, ISPIN,ISTATE,LSYM,NACTEL, MCTYPE,
     *                LSOLMX,NLMSOL,NELMN1,NELMX1,NELMN3,NELMX3,
     *                LROOTS,NROOTS,IROOT(MAXRTS),
     *                NOROT(MXCORB),IMOORD(MXCORB),
     *                IORTO,ICI0,KDEL,ICHECK,NTIT,
     *                MAXMAC,MAXMIC,MAXJT,MAXCIT,MAXUIT,MAXAPM,MAXABS,
     *                ITRLVL,ITRFIN,JCHSYM,JCHORB,
     *                NROOCI,ISTACI, MXCIMA, ICICNO,IMCCNO
      LOGICAL FLAG,  DORHF, DOMP2, DOCINO,DOCI,  DOMC,  DORSP,FCVORB,
     &        LNOROT,LMOORD,DIRFCK,CORHOL,CORRLX,RESPHP,JOLSEN,ABAIPH,
     &        INERSI,INERSF,SUPSYM,DODFT, DONEVPT,HSROHF ,BOYORB,PIPORB,
     &        ADDMP2
      COMMON /LOGINP/ FLAG(NFLAG), DORHF,DOMP2,DOCINO,DOCI,DOMC,DORSP,
     &                FCVORB,LNOROT,LMOORD,DIRFCK,CORHOL,CORRLX,RESPHP,
     &                JOLSEN,ABAIPH,INERSI,INERSF,DODFT,DONEVPT,HSROHF,
     &                BOYORB,PIPORB,ADDMP2
      EQUIVALENCE (SUPSYM,FLAG(17))
#if defined (SYS_CRAY)
      REAL             SPIN,POTNUC, EPSOL,EPSTAT,EPPN,RSOL,
#else
      DOUBLE PRECISION SPIN,POTNUC, EPSOL,EPSTAT,EPPN,RSOL,
#endif
     &                THRGRD, THRPWF, THRRHF, THRCI, THRMC, THRCGR,
     &                EFIELD, TITMOL, CMAXMO, THROVL,
     &                THRSSY, DEFLVL
      COMMON /RELINP/ SPIN,POTNUC, EPSOL,EPSTAT,EPPN,RSOL(3),
     &                THRGRD, THRPWF, THRCI, THRMC, THRCGR,
     &                EFIELD(MXFELT), TITMOL(12,2), CMAXMO, THROVL,
     &                THRSSY, DEFLVL
      CHARACTER*60 TITLE
      CHARACTER*4  CENT,   TYPE
      CHARACTER*8  LFIELD
      COMMON /CHRINP/ TITLE(6), CENT(MXCORB), TYPE(MXCORB),
     &                LFIELD(MXFELT)
