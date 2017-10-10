C   /* Deck check_kiadelta1 */
      SUBROUTINE CHECK_KIADELTA1(XINT,NBAS,CMO,NVIR,NRHF,T2AMSQ,LT2SQ,
     *               KGBD,KIBD,KIAD,WORK,LWORK,DELTA,LUPRI)
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
      REAL*8 XINT(NBAS,NBAS,NBAS)
      REAL*8 KIAD(NRHF,NRHF,NVIR*NBAS),KIBD(NRHF,NRHF,NBAS)
      REAL*8 KGBD(NRHF,NBAS,NBAS)
      REAL*8 CMO(NBAS,NRHF+NVIR),T2AMSQ(NRHF*NVIR,NRHF*NVIR)
      REAL*8 WORK(LWORK)
      INTEGER A,B,C,D,DELTA,GAM,BETA,ALPHA,LXINT,ADEL
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------------------------------------------------------------
C      1) Transfrom occupied index in T2 ampltude to Delta.
C------------------------------------------------------------------------
C
      ICOUNT = 1
C
      LXINT = NBAS*NBAS*NBAS
C
      DO 140 GAM = 1,NBAS
      GAMDEL = GAM + (DELTA-1)*NBAS
      DO 110 BETA = 1,NBAS
      DO 130 K = 1,NRHF !loop over occupied orbitals
      DO 120 ALPHA = 1,NBAS
              KGBD(K,BETA,GAM) =KGBD(K,BETA,GAM)+
     *                          XINT(ALPHA,BETA,GAM)*CMO(ALPHA,K)
C
  120          CONTINUE
C              WRITE(LUPRI,*) 'KGBD',KGBD(K,BETA,GAM)
  130       CONTINUE
  110    CONTINUE
  140 CONTINUE
C
      DO 210 GAM = 1,NBAS
         GAMDEL = GAM + (DELTA-1)*NBAS
         DO 220 I = 1,NRHF
            DO 230 K = 1,NRHF
               DO 240 BETA = 1,NBAS  !loop over occupied orbitals
C
                  KIBD(K,I,GAM) = KIBD(K,I,GAM) 
     *                              + KGBD(K,BETA,GAM)
     *                                              *CMO(BETA,I)
  240           CONTINUE
C              WRITE(LUPRI,*) 'KIBD',KIBD(K,I,GAM)
  230        CONTINUE 
  220     CONTINUE
  210  CONTINUE
C
C        WRITE(LUPRI,*) KIAD(1,1,1) 
C
      DO 320 I = 1,NRHF
      DO 310 A = NRHF+1,NRHF+NVIR
      ADEL = (A-NRHF) + (DELTA-1)*(A-NRHF)
      DO 330 K = 1,NRHF !loop over occupied orbitals
      DO 340 GAM = 1,NBAS

C
             KIAD(K,I,ADEL) = KIAD(K,I,ADEL) + KIBD(K,I,GAM)
     *                                              *CMO(GAM,A)
C
C            IF ((K .EQ. 1) .AND. (I .EQ. 1) .AND. (ADEL .EQ. 1)) THEN
C                WRITE(LUPRI,*) 'CMO ', CMO(GAM,A)
C                WRITE(LUPRI,*) 'KIBD', KIBD(K,I,GAM)
C                WRITE(LUPRI,*) 'KIAD', KIAD(K,I,ADEL)
C            ENDIF

  340          CONTINUE 
C                  WRITE(LUPRI,*) 'KIAD',ADEL,KIAD(K,I,ADEL)
  330       CONTINUE
  310    CONTINUE
  320 CONTINUE
C
C      DO 420 I = 1,NRHF
C      DO 410 A = NRHF+1,NRHF+NVIR
C         DO 420 I = 1,NRHF
c            DO 430 K = 1,NRHF !loop over occupied orbitals
C               DO 440 C = NRHF+1,NRHF+NVIR
C
C                              WRITE(LUPRI,*) 'CMO(GAM,A)',CMO(GAM,A)
C                              WRITE(LUPRI,*) 'A',A,'I',I,'K',K,'GAM',GAM
C                              WRITE(LUPRI,*) 'XINT',KIBD(K,I,GAM)
C                  KIAC(K,I,A) = KIAC(K,I,A) + KIAD(K,I,A)
C     *                                              *CMO(DELTA,C)
CC                  WRITE(LUPRI,*) 'K',K,'I',I,'A',A,'DELTA',DELTA,'C',C
CC                  WRITE(LUPRI,*) 'KIAC(K,I,A)',KIAC(K,I,A)
C  440          CONTINUE 
CC                WRITE(LUPRI,*) 'KIAC',KIAC(K,I,A)
C  430       CONTINUE
CC  420    CONTINUE
C  410 CONTINUE
C  420    CONTINUE
CC
C      ICOUNT = 1
C      ICOUN2 = 1
CC
C      DO 510 A = NRHF+1,NRHF+NVIR !Loop over virual orbitals
C         DO 520 I = 1,NRHF !Loop over occupied orbitals
C            AI = ICOUNT
C            ICOUNT = ICOUNT+1
C            DO 530 B = NRHF+1,NRHF+NVIR !loop over virtual orbitals
C               DO 540 J = 1,NRHF !Loop over occupied orbitals
C                  BJ = ICOUN2 
C                  ICOUN2 = ICOUN2+1
C                  DO 550 K = 1,NRHF
C                     DO 560 C = NRHF+1,NRHF+NVIR
CC
C                        KIACT2(K,I,A) = KIACT2(K,I,A)
C     *                          + KIAC(K,I,A) * T2AMSQ(AI,BJ)
C  560              CONTINUE
C  550            CONTINUE
C  540           CONTINUE
C  530       CONTINUE
C  520   CONTINUE
C  510 CONTINUE

      RETURN
      END
