c     -*- mode: fortran -*-
c     see abacus/herrdn.F for description.
      LOGICAL SHARE, SEGM, SPHR
      COMMON /SHELLS/ CENT(MXSHEL,3,2), NHKT(MXSHEL),   KHKT(MXSHEL),
     &                KCKT(MXSHEL),     ISTBAO(MXSHEL), NUCO(MXSHEL),
     &                JSTRT(MXSHEL),    NSTRT(MXSHEL),  MST(MXSHEL),
     &                NCENT(MXSHEL),    SHARE(MXSHEL),  NRCO(MXSHEL),
     &                NUMCF(MXSHEL),    NBCH(MXSHEL),   KSTRT(MXSHEL),
     &                SEGM(MXSHEL),     LCLASS(MXSHEL),
     &                IPTSHL(MXSHEL),   NUMCFT(MXSHEL), SPHR(MXSHEL),
     &                MBSID(MXSHEL), KMAX, NLRGSH, NSMLSH, NORBS,
     &                MBIDBT(MXSHEL), 
     &                NAUXSH, KMAXAUX, NORBAUX, KMAXTOT, IFSTSHAUX
C     MBSID and MBIDBT have been added for multiple basis sets (WK/UniKA/31-10-2002).
