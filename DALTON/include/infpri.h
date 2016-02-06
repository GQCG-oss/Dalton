! ---- include/infpri.h ---
!      INFo on PRInt in SIRIUS module (primarily)
      LOGICAL P4FLAG,P6FLAG, CMOPRI
      integer, parameter :: NPFLAG = 60
      integer         IPRI4, IPRI6, IPRSIR,                             &
     &                IPRCNO,IPRDIA,IPRSIG,IPRDNS,IPRSOL,               &
     &                IPRKAP,IPRMUL,IPRCIX,IPRRHF,IPRAVE,IPRFCK,        &
     &                MPRI4, MPRI6, MPRSIR,LIM_POPPRI
      COMMON /INFPRI/ IPRI4, IPRI6, IPRSIR,                             &
     &                IPRCNO,IPRDIA,IPRSIG,IPRDNS,IPRSOL,               &
     &                IPRKAP,IPRMUL,IPRCIX,IPRRHF,IPRAVE,IPRFCK,        &
     &                MPRI4, MPRI6, MPRSIR,LIM_POPPRI,                  &
     &                P4FLAG(NPFLAG),P6FLAG(NPFLAG), CMOPRI
! ---- end of include/infpri.h

