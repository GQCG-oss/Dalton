* ==================================================================== *
* CCFIELD_DEC.H
* -------------------------------------------------------------------- *
* explicite variable declaration for CCFIELD.H
* -------------------------------------------------------------------- *
      INTEGER MXFELT, NFIELD
      LOGICAL NHFFIELD

#if defined (SYS_CRAY)
      REAL EFIELD
#else
      DOUBLE PRECISION EFIELD
#endif
* -------------------------------------------------------------------- *
