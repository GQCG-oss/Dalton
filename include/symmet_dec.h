* ==================================================================== *
* SYMMET_DEC.H
* -------------------------------------------------------------------- *
* explicite variable declaration for SYMMET.H
* -------------------------------------------------------------------- *
      INTEGER MAXREP, MAXOPR, MULT, ISYMAX, ISYMAO, NPARSU,
     &        NAOS, NPARNU, IPTSYM, IPTCNT, NCRREP, IPTCOR,
     &        NAXREP, IPTAX, IPTXYZ, IPTNUC, ISOP, 
     &        NROTS, NINVC, NREFL, IXVAL, NCOS, ICLASS, ICNTAO

#if defined (SYS_CRAY)
       REAL FMULT, PT
#else
       DOUBLE PRECISION FMULT, PT
#endif
* -------------------------------------------------------------------- *
