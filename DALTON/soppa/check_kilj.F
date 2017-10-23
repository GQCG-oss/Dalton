C   /* Deck check_kilj */
      SUBROUTINE CHECK_KILJ(XINT,LXINT,NBAS,NRHF,NVIR,CMO,
     *      NLAMDT,T2AMSQ,LT2SQ,KILJ,WORK,LWORK,DELTA,ISYMD,LUPRI)
C
C   PORPOSE: 1)Set up memory allocation need here and in CHECK_AGBD1
C            2)Square up ao-integrals and move to 3dim array (INTG)
C            3)Move MO coef. to 2dim array (C) 
C   List of Variables:
C       XINT:      AO integrals (alpha gamma | beta delta)
C       LXINT:     The lenth of XINT
C       NBAS:      Number of Basis functions
C       NRHF:      Number of Restricted Hartree-Fock orbitals.
C       NVIR:      Number of Virtual orbitals.
C       CMO:       MO coeffecients.
C       NLAMDT:    Length of CMO array.
C       T2AMSQ:    Squared T2 amplitudes (first order doubles)
C       LT2SQ:     Lenght of LT2SQ
C       DELTA:     Outer loop index.
C       ISYMD:     Symmetry of delta. (The outer loop)
C       LUPRI:     Output file
C
C   Written by Lilli Irene Ør Kristensen, Summer 2017.
C
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      REAL*8 INTG(NBAS,NBAS,NBAS)
      REAL*8 XINT(LXINT),CMO(NLAMDT),T2AMSQ(LT2SQ),WORK(LWORK)
      REAL*8 XOUT(NBAS,NBAS,NBAS),KILJ(NRHF,NRHF,NRHF)
      INTEGER DELTA,GAM,LXINT,NBAS,NRHF,NVIR,NLAMDT,LT2SQ
      INTEGER NLTME
C
C
      CALL DZERO(WORK,LWORK)
C
C--------------------------------------------------------------------
C       Dynamic allocation
C--------------------------------------------------------------------
      LINTG = (NBAS*NBAS)*NBAS
      LKC = (NBAS*(NRHF+NVIR))
C
      KINTG  = 1
      KC     = KINTG  + LINTG
      KXOUT  = KC     + LKC
      KABGJ  = KXOUT  + LINTG
      KKBGJ  = KABGJ  + NBAS*NBAS*NBAS*NRHF
      KKIGJ  = KKBGJ  + NRHF*NBAS*NBAS*NRHF
      KKIJLD = KKIGJ  + NRHF*NRHF*NBAS*NRHF
      KEND1  = KKIJLD + NRHF*NRHF*NRHF*NRHF
C
      LWRK1 = LWORK - KEND1
C

      IF (LWRK1 .LT. 0) THEN
        WRITE(LUPRI,*) 'NEED:', KEND1,'Available:',LWORK
        CALL QUIT('Insufficitent space in CHECK_AGBD')
      ENDIF
C
C------------------------------------------------------------------
C       Square up AO integrals and move to 3dim array.
C------------------------------------------------------------------
      CALL DZERO(WORK(KXOUT),LINTG)
C
      NLTME = (((NBAS+1)*NBAS)/TWO) 
C               Number of Lower Triangular NBAS*NBAS-Matrix Elements
C
      DO 111 GAM = 1,NBAS !Loop over Gamma

          KOFF1 = 1 +(GAM-1)*NLTME
          KOFF2 = KXOUT+((NBAS*NBAS)*(GAM-1))
C
C          WRITE(LUPRI,*) 'KOFF1',KOFF1,'KOFF2',KOFF2,'NLTME',NLTME
C
          CALL CCSD_SYMSQ(XINT(KOFF1),ISYMD,WORK(KOFF2)) !square up XINT
C
  111 CONTINUE
C
       CALL DCOPY(LINTG,WORK(KXOUT),1,WORK(KINTG),1) 
C           Move integrals to WORK(KINTG)
C
C-------------------------------------------------------------------
C       Move MO-coeffecients to 2dim array.
C-------------------------------------------------------------------
C
      CALL DCOPY(NBAS*(NVIR+NRHF),CMO,1,WORK(KC),1)
C
C------------------------------------------------------------------
C       Call routines to calculate contribution
C---------------------------------------------------------------------
C
C  Transform the a and b index to gamma and delta. Then multiply on XINT
      CALL CHECK_KILJ1(WORK(KINTG),NBAS,WORK(KC),NVIR,NRHF,T2AMSQ,
     *            LT2SQ,KILJ,WORK(KABGJ),WORK(KKBGJ),WORK(KKIGJ),
     *             WORK(KKIJLD),WORK(KEND1),LWRK1,DELTA,LUPRI)
C
      RETURN
      END

