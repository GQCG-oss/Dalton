!
! Logical variables for pcm calculations
!
! PCM:    True if we perform a PCM calculation
! OUTFLD: Local field correction for pure liquids
      LOGICAL PCM,OUTFLD,NEWMAT,LOCFLD,NONEQ,NEQRSP,                           &
     &        NPCMIN,OLDCEN,NEWQR
      COMMON /PCMLOG/ PCM,OUTFLD,NEWMAT,LOCFLD,NONEQ,NEQRSP,                   &
     &                NPCMIN,OLDCEN,NEWQR
