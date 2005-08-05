C     PARAMETER (MXTSK=(MXSHEL*(MXSHEL+1)/2), IEIR=0, MMWORK=0,
C    &     MXSHL_KEEP=MAX(MXTSK,IEIR), MNSHL_KEEP=MIN(MXTSK,IEIR) )
Chj 11-May-2005: changes needed for compilation with g77:
C     MAX and MIN not allowed by "g77"; zero size arrays not allowed
      PARAMETER (MXTSK=(MXSHEL*(MXSHEL + 1)/2), IEIR=3000000,
     &           MMWORK=20000000,
     &          MXSHL_KEEP=MXTSK, MNSHL_KEEP=IEIR)
      DOUBLE PRECISION CORE
      LOGICAL AOSAVE, LINTSV, INITX, LINTMP, MSAVE
      COMMON /INCORE/ CORE(MMWORK), ISCORE, MMCORE, I_SHL, N_SHL,
     &     INDX_SHL(MXSHL_KEEP), INDX_SHL1, INDX_SHL2, INDX_SHL3, 
     &     INDX_SHL4, AOSAVE, LINTSV, INITX, LMCORE, JSCORE,
     &     INDX_C(3,MNSHL_KEEP), ITOTNT, LINTMP, MSAVE,
     &     IPREVA, IPREVB, IPREVC, IPREVD
#if 0
#if defined (VAR_MPI)
      LOGICAL INDMPI
      COMMON /MPISAV/ INDMPI(MXTSK)
#endif

/* MMWORK and IEIR must be set to nonzero values here if the 
  AOSAVE-feature of DALTON is to be used. Suggested values: 
  MMWORK=80000000 and IEIR=3000000, 
  if your hardware can take it.*/
#endif
