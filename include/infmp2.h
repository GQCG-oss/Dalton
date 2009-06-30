C  -- infmp2.h -- Do not change without checking with Hans Joergen Aa. Jensen  /l.r. Sep 2007, hjaaj/
C
C     Information for ".MP2" (sirius/sirmp2.F) and ".SOPPA" in RESPONS (rsp/rspsoppa.F)
C     Also used for ".SOPPA" in abacus/abadrv.F.
C
C     Off-set to (bj) block of kappa(ai,bj) is IADR1(B,J), which is calculated when needed.
C
      PARAMETER (MAXVIR = MXCORB)
      LOGICAL MP2FRO, SC_SRMP2, SOSMP2, SCSMP2, TDAMP2, MP2_NO_OV
      COMMON /INFMP2/ IPRMP2, NFRMP2(8), NPHSYM(8), NPHMAX, NP, NH,
     &        IPHORD(MXCORB), INDXPH(MXCORB), IFRMP2(MXCORB), NPHTOT(8),
     &        MP2FRO, SC_SRMP2, SOSMP2, SCSMP2, TDAMP2, MP2_NO_OV
C  -- end of infmp2.h --
