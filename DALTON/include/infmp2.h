C  -- infmp2.h -- Do not change without checking with Hans Joergen Aa. Jensen  /l.r. Nov 2009, hjaaj/
C
C     Information for ".MP2" (sirius/sirmp2.F) and ".SOPPA" in RESPONS (rsp/rspsoppa.F)
C     Also used for ".SOPPA" in abacus/abadrv.F.
C
C     Off-set to (bj) block of kappa(ai,bj) is IADR1(B,J), which is calculated when needed.
C
C     MP2_LSHIFT      = 0.0d0 ! level shift for 2p-2h orbital energy difference
C     MP2_SCALEFAC(1) = 1.0d0 ! scale factor on (ij|ab)
C     MP2_SCALEFAC(2) = 1.0d0 ! scale factor on <ij|ab> = (ia|jb)
C     MP2_SCALED    = .FALSE. ! MP2_SCALEFAC .ne. 1, 1
C     MP2_SCS       = .FALSE. ! Grimme's SCS scale factors
C     MP2_SOS       = .FALSE. ! Head-Gordon's SOS scale factors
C     MP2_TDA       = .FALSE. ! save info for SOPPA(TDA)
C     MP2_NO_OCCVIR = .FALSE. ! ignore occ-vir block of MP2 density matrix
C     SRMP2_SELFCONSISTENT = .FALSE. ! self-consistent iterations on occ-vir block of srMP2 density matrix
C     SRMP2_SRINTS  = .FALSE. ! calculate both sr-sr and lr-lr MP2 energies - advanced option for study of MP2 energy contributions
C     LAMSR                   ! fraction of srHFX and lr-sr MP2 coupling in the RSDHf energy 
C
C     DCPT2       : degeneracy corrected PT2
C     SAVE_MP2WF1 : calculate and save MP2 coefficients on SIRIFC
C
      REAL*8  MP2_SCALEFAC, MP2_LSHIFT, LAMSR
      LOGICAL MP2FRO, MP2_SOS, MP2_SCS, MP2_TDA, MP2_NO_OCCVIR,
     &        MP2_SCALED, SRMP2_SELFCONSISTENT, SRMP2_SRINTS,
     &        DCPT2, DO_RSDHF, SAVE_MP2WF1
      COMMON /INFMP2/ MP2_SCALEFAC(2), MP2_LSHIFT, LAMSR,
     &        IPRMP2, NFRMP2(8), NPHSYM(8), NPHMAX, NP, NH,
     &        IPHORD(MXCORB), INDXPH(MXCORB), IFRMP2(MXCORB), NPHTOT(8),
     &        MP2FRO, MP2_SOS, MP2_SCS, MP2_TDA, MP2_NO_OCCVIR,
     &        MP2_SCALED, SRMP2_SELFCONSISTENT, SRMP2_SRINTS,
     &        DCPT2, DO_RSDHF, SAVE_MP2WF1
C  -- end of infmp2.h --
