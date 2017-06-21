C   /* Deck check_agbd1 */
      SUBROUTINE CHECK_AGBD1(XINT,NBAS,CMO,NVIR,NRHF,T2AMSQ,LT2SQ,
     *                   PQRDEL,T2IJAD,T2IJGD,WORK,LWORK,DELTA,LUPRI)
C
C   PORPOSE: Calculate Sum_cd T^(c delta)_(ij) (alpha gamma | beta delta)
C           1) Transform the c index to delta: T_(ij)^(c delta).
C           2) Contract the integral with the amplitude.
C
C   List of Variables:
C       XINT:    Squared up AO integrals (alpha gamma | beta delta)
C       NBAS:    Number of basis functions
C       CMO:     MO coeffecients in (AO,MO)-array.
C       NVIR:    Number of virtual orbitals
C       NRHF:    Number of Restricted Hartree-Fock arbitals.
C       T2AMSQ:  Squared T2 amplitudes (first order doubles)
C       LT2SQ:   Lenth of T2AMSQ
C       PQRDEL:  Array for (alpha beta|gamma delta) T_(ij)^(gamma delta)
C       T2IJAD:  Array for T_(ij)^(a delta)
C       T2IJGD:  Array for T_(ij)^(gamma delta)
C       DELTA:   A loop counter, looping over the ao index, delta.
C       LUPRI:   Name of output file.
C
C   Written by Lilli Irene Ør Kristensen, Summer 2017.
C
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      REAL*8 XINT(NBAS,NBAS,NBAS),PQRDEL(NBAS,NBAS,NBAS)
      REAL*8 CMO(NBAS,NRHF+NVIR),T2AMSQ(NRHF*NVIR,NRHF*NVIR)
      REAL*8 T2IJAD(NRHF+NVIR,NRHF+NVIR),T2IJGD(NRHF+NVIR,NRHF+NVIR)
      REAL*8 WORK(LWORK)
      INTEGER A,I,B,J,C,D,AI,BJ,DELTA,GAM,BETA,ALPHA,GA,LXINT,NRHF,NVIR
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------------------------------------------------------------
C      1) Transfrom occupied index in T2 ampltude to Delta.
C------------------------------------------------------------------------
C
      ICOUNT = 1
      ICOUN2 = 1
C
      LXINT = NBAS*NBAS*NBAS
C
      DO 110 A = NRHF+1,NRHF+NVIR !Loop over virual orbitals
         DO 120 I = 1,NRHF !Loop over occupied orbitals
            AI = ICOUNT
            ICOUNT = ICOUNT+1
            DO 130 B = NRHF+1,NRHF+NVIR !loop over virtual orbitals
               DO 140 J = 1,NRHF !Loop over occupied orbitals
                  BJ = ICOUN2 
                  ICOUN2 = ICOUN2+1
                  DO 150 C = NRHF+1,NRHF+VIR
                     DO 160 D = NRHF+1,NRHF+VIR
                        DO 170 GAM = 1,NBAS
                           DO 180 BETA = 1,NBAS
                              DO 190 ALPHA = 1,NBAS
C   Transforming the b index of the amplitude to delta.
                              T2IJAD(AI,BJ) = T2AMSQ(AI,BJ)*CMO(DELTA,B)
C   Transforming the a index of the amplitude to gamma
                              T2IJGD(AI,BJ) = T2IJAD(AI,BJ)*CMO(GAM,A)
C Computing the the AGBD-intermediate: T_(ij)^(gamma delta) (alpha beta|gamma delta)
                              PQRDEL(ALPHA,BETA,GAM) = T2IJGD(AI,BJ) *
     *                                  XINT(ALPHA,BETA,GAM)
C Tranforming the AGBD-intermediate to MO-basis.
                               WORK  =
     *                         PQRDEL(ALPHA,BETA,GAM)*CMO(ALPHA,A)
     *                          *CMO(BETA,D)*CMO(GAM,B)*CMO(DELTA,D)
C                               WRITE(LUPRI,*) 'WORK(LILLI)',WORK

  190                   CONTINUE
  180                 CONTINUE
  170                CONTINUE
  160              CONTINUE
  150            CONTINUE
  140           CONTINUE
  130       CONTINUE
  120   CONTINUE
  110 CONTINUE
C
      RETURN
      END

