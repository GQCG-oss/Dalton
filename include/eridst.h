C
C$Id: eridst.h,v 1.1.1.1 2001-02-08 13:33:26 hjj Exp $
C
      PARAMETER (MXDIST = 100)
      COMMON /ERIDST/ NDISTR, MLTDST, KHKDST, NACDST, IRPDST,
     &                MAXCML, NTCLAS,
     &                INDDST(MXDIST), IACDST(MXDIST),
     &                INDXDS(MXSHEL), NCLASS(MXSHEL), MCLASS(MXSHEL),
     &                KCLASS(MXSHEL), KLAOBT(MXSHEL)
