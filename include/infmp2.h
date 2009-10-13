C  -- infmp2.h -- Do not change without checking with Hans Joergen Aa. Jensen  /l.r. Oct 2009, hjaaj/
C
C     Information for ".MP2" (sirius/sirmp2.F) and ".SOPPA" in RESPONS (rsp/rspsoppa.F)
C     Also used for ".SOPPA" in abacus/abadrv.F.
C
C     Off-set to (bj) block of kappa(ai,bj) is IADR1(B,J), which is calculated when needed.
C
      PARAMETER (MAXVIR = MXCORB)
      LOGICAL MP2FRO, SOS_MP2, SCS_MP2, TDA_MP2, MP2_NO_OV,
     &        SC_SRMP2, SRINTS_SRMP2
      COMMON /INFMP2/ IPRMP2, NFRMP2(8), NPHSYM(8), NPHMAX, NP, NH,
     &        IPHORD(MXCORB), INDXPH(MXCORB), IFRMP2(MXCORB), NPHTOT(8),
     &        MP2FRO, SOS_MP2, SCS_MP2, TDA_MP2, MP2_NO_OV,
     &        SC_SRMP2, SRINTS_SRMP2
C  -- end of infmp2.h --
