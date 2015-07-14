      PARAMETER ( MXSMOP = MAXLBL , MBSMFR = MAXLBL)
      LOGICAL SOMOM,TWOPHO,ASMOP,BSMOP, MCDCAL,LTPCD
      LOGICAL MCDPRJ
      CHARACTER*8 ASMLB,BSMLB
      COMMON /INFSMO/ BSMFR(MBSMFR), NASMOP(8),                         &
     &                NBSMOP(8), NSMCNV(8),                             &
     &                ASMOP(MXSMOP), BSMOP(MXSMOP),                     &
     &                SOMOM,TWOPHO,MCDCAL,IPRSMO,LTPCD,                 &
     &                MCDPRJ,                                           &
     &                NBSMFR
      COMMON /CHRSMO/ ASMLB(8,MXSMOP), BSMLB(8,MXSMOP)
