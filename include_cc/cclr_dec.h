* ==================================================================== *
* CCLR_DEC.H
* -------------------------------------------------------------------- *
* explicite variable declaration for CCLR.H
* -------------------------------------------------------------------- *
      INTEGER MAXRED, MAXITE
      INTEGER ISYMTR, NSIMLE, NCCVAR

#if defined (SYS_CRAY)
      REAL ECURR, THREXC, THRLEQ, THRENR, THRLDPHF
#else
      DOUBLE PRECISION ECURR, THREXC, THRLEQ, THRENR, THRLDPHF
#endif
* -------------------------------------------------------------------- *
