      LOGICAL RESTPP,ABCHK,ABSYM,RESTLR,RESTC6,DETERM
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      REAL
#else
      DOUBLE PRECISION
#endif
     * AFREQ
      INTEGER
     *                KSYMOP, KZVAR,  KZYVAR, KZCONF, KZYCON, KZWOPT, 
     *                KZYWOP, KZRED , KZYRED, KEXCNV, KEXSIM, KEXSTV,
     *                KLRSTV, JEXSIM, KOFFTY, KCONV,  KSYMST
      COMMON /WRKRSP/ AFREQ,  
     *                KSYMOP, KZVAR,  KZYVAR, KZCONF, KZYCON, KZWOPT, 
     *                KZYWOP, KZRED , KZYRED, KEXCNV, KEXSIM, KEXSTV,
     *                KLRSTV, JEXSIM, KOFFTY, KCONV,  KSYMST, 
     *                RESTPP, ABCHK,  ABSYM,  RESTLR, RESTC6, DETERM
