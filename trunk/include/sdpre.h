#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      SDPRE(I) = FLOAT(I)
#else
      SDPRE(I) = DFLOAT(I)
#endif
