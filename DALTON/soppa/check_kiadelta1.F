C   /* Deck check_kiadelta1 */
      SUBROUTINE CHECK_KIADELTA1(XINT,NBAS,CMO,NVIR,NRHF,T2AMSQ,LT2SQ,
     *               KGBD,KIBD,KIAD,KIAC,KIACT2,WORK,LWORK,DELTA,LUPRI)
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
      REAL*8 XINT(NBAS,NBAS,NBAS),KIAC(NRHF,NRHF,NVIR)
      REAL*8 KIAD(NRHF,NRHF,NVIR),KIBD(NRHF,NRHF,NBAS)
      REAL*8 KGBD(NRHF,NBAS,NBAS),KIACT2(NRHF,NRHF,NVIR)
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
      DO 110 GAM = 1,NBAS
C        DO 115 I = 1,NRHF
          DO 120 BETA = 1,NBAS
            DO 130 K = 1,NRHF !loop over occupied orbitals
               DO 140 ALPHA = 1,NBAS
C
                  KGBD(K,BETA,GAM) =KGBD(K,BETA,GAM)+
     *                          XINT(ALPHA,BETA,GAM)*CMO(ALPHA,K)
C
  140          CONTINUE
  130       CONTINUE
  120    CONTINUE
  110 CONTINUE
C
      DO 210 GAM = 1,NBAS
         DO 220 I = 1,NRHF
            DO 230 K = 1,NRHF
               DO 240 BETA = 1,NBAS  !loop over occupied orbitals
C
                  KIBD(K,I,GAM) = KIBD(K,I,GAM) + KGBD(K,BETA,GAM)
     *                                              *CMO(BETA,I)
  240           CONTINUE
  230        CONTINUE 
  220     CONTINUE
  210  CONTINUE
C
      DO 310 A = NRHF+1,NRHF+NVIR
         DO 320 I = 1,NRHF
            DO 330 K = 1,NRHF !loop over occupied orbitals
               DO 340 GAM = 1,NBAS
C
C                              WRITE(LUPRI,*) 'CMO(GAM,A)',CMO(GAM,A)
C                              WRITE(LUPRI,*) 'A',A,'I',I,'K',K,'GAM',GAM
C                              WRITE(LUPRI,*) 'XINT',KIBD(K,I,GAM)
                  KIAD(K,I,A) = KIAD(K,I,A) + KIBD(K,I,GAM)
     *                                              *CMO(GAM,A)
C                  WRITE(LUPRI,*) 'K',K,'I',I,'A',A,'GAM',GAM
C                  WRITE(LUPRI,*) 'KIAD(K,I,A)',KIAD(K,I,A)
  340          CONTINUE 
  330       CONTINUE
  320    CONTINUE
  310 CONTINUE
C
      DO 410 A = NRHF+1,NRHF+NVIR
         DO 420 I = 1,NRHF
            DO 430 K = 1,NRHF !loop over occupied orbitals
               DO 440 C = NRHF+1,NRHF+NVIR
C
C                              WRITE(LUPRI,*) 'CMO(GAM,A)',CMO(GAM,A)
C                              WRITE(LUPRI,*) 'A',A,'I',I,'K',K,'GAM',GAM
C                              WRITE(LUPRI,*) 'XINT',KIBD(K,I,GAM)
                  KIAC(K,I,A) = KIAC(K,I,A) + KIAD(K,I,A)
     *                                              *CMO(DELTA,C)
C                  WRITE(LUPRI,*) 'K',K,'I',I,'A',A,'DELTA',DELTA,'C',C
C                  WRITE(LUPRI,*) 'KIAC(K,I,A)',KIAC(K,I,A)
  440          CONTINUE 
  430       CONTINUE
  420    CONTINUE
  410 CONTINUE
C
      ICOUNT = 1
      ICOUN2 = 1
C
      DO 510 A = NRHF+1,NRHF+NVIR !Loop over virual orbitals
         DO 520 I = 1,NRHF !Loop over occupied orbitals
            AI = ICOUNT
            ICOUNT = ICOUNT+1
            DO 530 B = NRHF+1,NRHF+NVIR !loop over virtual orbitals
               DO 540 J = 1,NRHF !Loop over occupied orbitals
                  BJ = ICOUN2 
                  ICOUN2 = ICOUN2+1
                  DO 550 K = 1,NRHF
                     DO 560 C = NRHF+1,NRHF+NVIR
C
                        KIACT2(K,I,A) = KIACT2(K,I,A)
     *                          + KIAC(K,I,A) * T2AMSQ(AI,BJ)

  560              CONTINUE
  550            CONTINUE
  540           CONTINUE
  530       CONTINUE
  520   CONTINUE
  510 CONTINUE
      RETURN
      END
