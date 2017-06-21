C   /* Deck check_KILJ1 */
      SUBROUTINE CHECK_KILJ1(XINT,NBAS,CMO,NVIR,NRHF,T2AMSQ,LT2SQ,
     *                   KILJ,WORK,LWORK,DELTA,LUPRI)
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
C       KILJ:  Array for (k i|l j)
C       DELTA:   A loop counter, looping over the ao index, delta.
C       LUPRI:   Name of output file.
C
C   Written by Lilli Irene Ør Kristensen, Summer 2017.
C
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      REAL*8 XINT(NBAS,NBAS,NBAS),KILJ(NBAS,NBAS,NBAS)
      REAL*8 CMO(NBAS,NRHF+NVIR),T2AMSQ(NRHF*NVIR,NRHF*NVIR)
      REAL*8 WORK(LWORK)
      INTEGER A,B,C,D,DELTA,GAM,BETA,ALPHA,LXINT
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------------------------------------------------------------
C      1) Transfrom occupied index in T2 ampltude to Delta.
C------------------------------------------------------------------------
C
      ICOUNT = 1
C
      LXINT = NBAS*NBAS*NBAS
C
      DO 110 I  = 1,NRHF !Loop over occupied orbitals
        DO 120 J = 1,NRHF !Loop over occupied orbitals
C            AI = ICOUNT
C            ICOUNT = ICOUNT+1
            DO 130 K = 1,NRHF !loop over occupied orbitals
                DO 140 L = 1,NRHF !Loop over occupied orbitals
C                     BJ = J
                     DO 150 GAM = 1,NBAS
                      DO 160 BETA = 1,NBAS
                            DO 170 ALPHA = 1,NBAS
C
C                            WRITE(LUPRI,*) 'CMO(A,alpha)',CMO(K,ALPHA)
C                            WRITE(LUPRI,*) 'CMO(B,beta)',CMO(I,BETA)
C                            WRITE(LUPRI,*) 'CMO(C,gamma)',CMO(L,GAM)
C                            WRITE(LUPRI,*) 'CMO(A,delta)',CMO(J,DELTA)
C   Computing the GAMMA-intermediate - (ki|lj)
C                            WRITE(LUPRI,*) 'CHECK-(K beta|gamma delta)'
C                            KILJ(ALPHA,BETA,GAM) = XINT(ALPHA,BETA,GAM)*
C     *                                         CMO(K,ALPHA)
C                            WRITE(LUPRI,*) KILJ(ALPHA,BETA,GAM)
C                            WRITE(LUPRI,*) 'CHECK-(ki | gamma delta)'
C                            KILJ(ALPHA,BETA,GAM) = XINT(ALPHA,BETA,GAM)*
C     *                                         CMO(K,ALPHA)*CMO(I,BETA)
C                            WRITE(LUPRI,*) KILJ(ALPHA,BETA,GAM)
C                            WRITE(LUPRI,*) 'CHECK-(K i| l delta)'
C                            KILJ(ALPHA,BETA,GAM) = XINT(ALPHA,BETA,GAM)*
C     *                                         CMO(K,ALPHA)*CMO(I,BETA)*
C     *                                         CMO(L,GAM)
C                            WRITE(LUPRI,*) KILJ(ALPHA,BETA,GAM)
C                            WRITE(LUPRI,*) 'CHECK-(K i| l j)'
                            KILJ(ALPHA,BETA,GAM) = XINT(ALPHA,BETA,GAM)*
     *                                         CMO(K,ALPHA)*CMO(I,BETA)*
     *                                         CMO(L,GAM)*CMO(J,DELTA)
C                            WRITE(LUPRI,*) KILJ(ALPHA,BETA,GAM)
  170                   CONTINUE
  160                 CONTINUE
  150                CONTINUE
  140           CONTINUE
  130       CONTINUE
  120   CONTINUE
  110 CONTINUE
C
      RETURN
      END
