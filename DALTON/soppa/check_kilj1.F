C   /* Deck check_KILJ1 */
      SUBROUTINE CHECK_KILJ1(XINT,NBAS,CMO,NVIR,NRHF,T2AMSQ,LT2SQ,
     *                   KILJ,ABGJ,KBGJ,KIGJ,KIJLD,
     *                   WORK,LWORK,DELTA,LUPRI)
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
      REAL*8 XINT(NBAS,NBAS,NBAS),KILJ(NRHF,NRHF)
      REAL*8 CMO(NBAS,NRHF+NVIR),T2AMSQ(NRHF*NVIR,NRHF*NVIR)
      REAL*8 ABGJ(NBAS,NBAS,NRHF*NBAS),KBGJ(NRHF,NBAS,NRHF*NBAS)
      REAL*8 KIGJ(NRHF,NRHF,NBAS*NRHF),KIJLD(NRHF,NRHF,NRHF)
      REAL*8 WORK(LWORK)
      INTEGER A,B,C,D,DELTA,GAM,BETA,ALPHA,DELJ,LXINT
C
C     INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      CALL DZERO(KBGJ(1,1,1),NRHF*NBAS*NBAS*NRHF)
      CALL DZERO(KIJLD,NRHF*NRHF*NRHF*NRHF)
C

C------------------------------------------------------------------------
C      1) Transfrom occupied index in T2 ampltude to Delta.
C------------------------------------------------------------------------
C
      ICOUNT = 1
C
      LXINT = NBAS*NBAS*NBAS
      LXINT2 = LXINT*NRHF
C
      CALL DZERO(ABGJ,LXINT2)
c      DO 110 L = 1,NRHF !loop over occupied orbitals
      DO 120 J = 1,NRHF !Loop over occupied orbitals
c      DO 130 I = 1,NRHF
c      DO 140 K = 1,NRHF !Loop over occupied orbitals   DELJ = J + (DELTA-1)*NRHF
C         DELJ = J + (DELTA-1)*NRHF
      DO 150 GAM = 1,NBAS
      DO 160 BETA = 1,NBAS
      DO 170 ALPHA = 1,NBAS
C
      ABGJ(ALPHA,BETA,J) = ABGJ(ALPHA,BETA,J)
     *                      +  XINT(ALPHA,BETA,GAM) * CMO(GAM,J)
  170 CONTINUE
  160 CONTINUE
  150 CONTINUE
c  140 CONTINUE
c  130 CONTINUE
  120 CONTINUE
c  110 CONTINUE
C
C      DO 137 J =1,NRHF
C      WRITE(LUPRI,*) 'ABGJ'
C      CALL OUTPUT(ABGJ(1,1,J),1,NBAS,1,NBAS,NBAS,NBAS,1,LUPRI)
C 137  CONTINUE
C
C
C      DO 210 L = 1,NRHF !loop over occupied orbitals
      DO 220 J = 1,NRHF !Loop over occupied orbitals
C      DO 230 I = 1,NRHF
C      DO 240 K = 1,NRHF !Loop over occupied orbitals
C         LJ = L +(J-1)*NRHF
C      DO 250 GAM = 1,NBAS
      DO 260 BETA = 1,NBAS
      DO 240 K = 1,NRHF !Loop over occupied orbitals
      DO 270 ALPHA = 1,NBAS
C
      KBGJ(K,BETA,J) = KBGJ(K,BETA,J) 
     *   +  ABGJ(ALPHA,BETA,J) *  CMO(ALPHA,K)

c      IF ((K .EQ. 1) .AND. (BETA .EQ. 1) .AND. (J .EQ. 1)) THEN
C        WRITE(LUPRI,*) '(1,1,1)'
C        WRITE(LUPRI,*) 'ABGJ', ABGJ(ALPHA,BETA,J)
C        WRITE(LUPRI,*) 'CMO ', CMO(ALPHA,K)
C        WRITE(LUPRI,*) 'KBGJ', KBGJ(K,BETA,J)
c      END IF

c      IF ((K .EQ. 1) .AND. (BETA .EQ. 2) .AND. (J .EQ. 1)) THEN
C        WRITE(LUPRI,*) '(1,2,1)'
C        WRITE(LUPRI,*) 'ABGJ', ABGJ(ALPHA,BETA,J)
C        WRITE(LUPRI,*) 'CMO ', CMO(ALPHA,K)
C        WRITE(LUPRI,*) 'KBGJ', KBGJ(K,BETA,J)
C      END IF
  270 CONTINUE
  240 CONTINUE
  260 CONTINUE
C  250 CONTINUE
C  240 CONTINUE
C  230 CONTINUE
  220 CONTINUE
C  210 CONTINUE


C      DO 237 J =1,NRHF
C      WRITE(LUPRI,*) 'KGBJ'
C      CALL OUTPUT(KBGJ(1,1,J),1,NRHF,1,NBAS,NRHF,NBAS,1,LUPRI)
c 237  CONTINUE


      DO 320 J = 1,NRHF !Loop over occupied orbitals
      DO 330 I = 1,NRHF
      DO 340 K = 1,NRHF !Loop over occupied orbitals
C         LJ = L +(J-1)*NRHF
      DO 360 BETA = 1,NBAS
C
      KIGJ(K,I,J) = KIGJ(K,I,J) +  KBGJ(K,BETA,J)*
     *                                         CMO(BETA,I)
C
  360 CONTINUE
C        WRITE(LUPRI,*) 'KIGJ', KIGJ(K,I,J)
  340 CONTINUE
  330 CONTINUE
  320 CONTINUE
C      WRITE(LUPRI,*) 'KIGJ'
C      DO 337 J =1,NRHF
C      CALL OUTPUT(KIGJ(1,1,J),1,NRHF,1,NRHF,NRHF,NRHF,1,LUPRI)
C 337  CONTINUE

      DO 410 J = 1,NRHF
      DO 420 L = 1,NRHF !Loop over occupied orbitals
        LJ = J + (L-1)*NRHF
      DO 430 I = 1,NRHF
      DO 440 K = 1,NRHF !Loop over occupied orbitals
      
C
      KIJLD(K,I,J) = KIJLD(K,I,J) + KIGJ(K,I,J) * CMO(DELTA,L)

  440 CONTINUE
  430 CONTINUE
  420 CONTINUE
  410 continue

      WRITE(LUPRI,*) 'KIJLD'
      DO 437 L =1,NRHF
      CALL OUTPUT(KIJLD(1,1,L),1,NRHF,1,NRHF,NRHF,NRHF,1,LUPRI)
 437  CONTINUE

C----------------------------------------------------------------
C           HER TIL VIRKER KODEN!!!!
C-------------------------------------------------------------------

      DO 520 LJ = 1,NRHF !Loop over occupied orbitals
      DO 530 I = 1,NRHF !Loop over occupied orbitals
      DO 540 K = 1,NRHF !Loop over occupied orbitals
        KI = I

      KILJ(KI,LJ)= KILJ(KI,LJ)+ KIJLD(K,I,LJ) 

      WRITE(LUPRI,*) KI,LJ,KILJ(KI,LJ)
  540 CONTINUE
  530 CONTINUE
  520 CONTINUE

      WRITE(LUPRI,*) 'KILJ'
      CALL OUTPUT(KILJ(1,1),1,25,1,1,1,1,1,LUPRI)






      RETURN
      END
