C
C$Id: shells.h,v 1.1.1.1 2001-02-08 13:33:28 hjj Exp $
C
      LOGICAL SHARE, SEGM, SPHR
      COMMON /SHELLS/ CENT(MXSHEL,3,2), NHKT(MXSHEL),   KHKT(MXSHEL),
     &                KCKT(MXSHEL),     ISTBAO(MXSHEL), NUCO(MXSHEL),
     &                JSTRT(MXSHEL),    NSTRT(MXSHEL),  MST(MXSHEL),
     &                NCENT(MXSHEL),    SHARE(MXSHEL),  NRCO(MXSHEL),
     &                NUMCF(MXSHEL),    NBCH(MXSHEL),   KSTRT(MXSHEL),
     &                SEGM(MXSHEL),     LCLASS(MXSHEL),
     &                IPTSHL(MXSHEL),   NUMCFT(MXSHEL), SPHR(MXSHEL),
     &                KMAX, NLRGSH, NSMLSH, NORBS
