      PARAMETER (MXTSK=(MXSHEL*(MXSHEL+1)/2), IEIR=0, 
     &     MXSHL_KEEP=MAX(MXTSK,IEIR), MMWORK=0)
      LOGICAL AOSAVE, LINTSV, INITX, LINTMP, MSAVE
      COMMON /INCORE/ CORE(MMWORK), ISCORE, MMCORE, I_SHL, N_SHL,
     &     INDX_SHL(MXSHL_KEEP), INDX_SHL1, INDX_SHL2, INDX_SHL3, 
     &     INDX_SHL4, AOSAVE, LINTSV, INITX, LMCORE, JSCORE,
     &     INDX_C(3,MIN(MXTSK,IEIR)), ITOTNT, LINTMP, MSAVE,
     & 	   IPREVA, IPREVB, IPREVC, IPREVD
      /*#if defined (VAR_MPI)
      LOGICAL INDMPI
      COMMON /MPISAV/ INDMPI(MXTSK)
      #endif*/
/* MMWORK and IEIR must be set to nonzero values here if the 
  AOSAVE-feature of DALTON is to be used. Suggested values: 
  MMWORK=80000000 and IEIR=3000000, 
  if your hardware can take it.*/
