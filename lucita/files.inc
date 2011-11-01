!     Files to be used in LUCITA
      integer :: NCCFIL, MSLVOUT
      PARAMETER (NCCFIL=11) ! Number of word addressable files
!     Define the unit numbers
      PARAMETER (MSLVOUT = 12) ! Output from slaves
!     The standard input and output units may be redefined in the code
      COMMON /FILES/ IN,IW
      COMMON/PARNAM/ LUCITABASF, LUCIFILN, LUDIA_MPIN,NLUDIA_MPI,       &
     &               LUC_MPIN, NLUC_MPI, LUHC_MPIN, NLUHC_MPI,          &
     &               LUSC1_MPIN, LUSC2_MPIN, LUSC3_MPIN, NLUSC1_MPI,    &
     &               NLUSC2_MPI, NLUSC3_MPI, LUSC34_MPIN, LUSC35_MPIN,  &
     &               LUSC36_MPIN, LUSC37_MPIN, LUSC38_MPIN,             &
     &               LUSC40_MPIN, NLUSC34_MPI, NLUSC35_MPI, NLUSC36_MPI,&
     &               NLUSC37_MPI, NLUSC38_MPI, NLUSC40_MPI,             &
     &               NLUSC41_MPI,NDUMPF
      CHARACTER*10 LUCITABASF  ! generic name for LUCITA output file
      CHARACTER*14 LUCIFILN    ! node-unique LUCITA output file name
      CHARACTER*5  LUSC1_MPIN, LUSC2_MPIN, LUSC3_MPIN
      CHARACTER*9  NLUSC1_MPI, NLUSC2_MPI, NLUSC3_MPI
      CHARACTER*12 LUC_MPIN, LUHC_MPIN, LUDIA_MPIN
      CHARACTER*16 NLUC_MPI, NLUHC_MPI, NLUDIA_MPI
      CHARACTER*6  LUSC34_MPIN,LUSC35_MPIN,LUSC36_MPIN,LUSC37_MPIN
      CHARACTER*6  LUSC38_MPIN, LUSC40_MPIN
      CHARACTER*6  LUSC41_MPIN, NDUMPF
      CHARACTER*10 NLUSC34_MPI, NLUSC35_MPI, NLUSC36_MPI,NLUSC37_MPI
      CHARACTER*10 NLUSC38_MPI, NLUSC40_MPI, NLUSC41_MPI
