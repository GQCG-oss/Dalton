* ==================================================================== *
* blocks_dec.h
* -------------------------------------------------------------------- *
* explicite variable declaration for blocks.h
* -------------------------------------------------------------------- *
       INTEGER  MAXSHL, NLRGBL, NSMLBL, NHKTSH, KHKTSH, KCKTSH,
     &          ISTBSH, NUCOSH, NORBSH, NSTRSH, NCNTSH, NSETSH,
     &          IORBSB, NRCSH,  LCLASH, NO2INT, NLRBL,  ISYMBL, NSYMBL,
     &          MBSISH
C     MBSISH has been added for multiple basis sets (WK/UniKA/31-10-2002).
#if defined (SYS_CRAY)
      REAL CENTSH
#else
      DOUBLE PRECISION CENTSH
#endif
* -------------------------------------------------------------------- *
