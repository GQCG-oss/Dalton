C   /* Deck check_aikdelta1 */
      SUBROUTINE CHECK_AIKDELTA1(XINT,NBAS,CMO,NVIR,NRHF,T2AMSQ,LT2SQ,
     *              ALGKD,ALIKD,AIKD,AIKC,AIKCT2,WORK,LWORK,DELTA,LUPRI)
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
      REAL*8 XINT(NBAS,NBAS,NBAS),AIKC(NVIR,NRHF,NRHF)
      REAL*8 AIKD(NVIR,NRHF,NRHF),ALIKD(NBAS,NRHF,NRHF)
      REAL*8 ALGKD(NBAS,NBAS,NRHF),AIKCT2(NVIR,NRHF,NRHF)
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
      DO 130 K = 1,NRHF
         DO 110 BETA = 1,NBAS
            DO 120 ALPHA = 1,NBAS
               DO 140 GAM = 1,NBAS
C
                   ALGKD(ALPHA,BETA,K)=ALGKD(ALPHA,BETA,K)+
     *                          XINT(ALPHA,BETA,GAM)*CMO(GAM,K)
C
  140          CONTINUE
  120       CONTINUE
  110    CONTINUE
  130 CONTINUE
C
      DO 137 K = 1,NRHF
        WRITE(LUPRI,*) '(ALPHA BETA| K DELTA)'
        CALL OUTPUT(ALGKD(1,1,K),1,NBAS,1,NBAS,NBAS,NBAS,1,LUPRI)
  137 CONTINUE

      DO 210 K = 1,NRHF
         DO 220 ALPHA = 1,NBAS
            DO 230 I = 1,NRHF
               DO 240 BETA = 1,NBAS 
C
                 ALIKD(ALPHA,I,K) = ALIKD(ALPHA,I,K)+ALGKD(ALPHA,BETA,K)
     *                                              *CMO(BETA,I)
  240           CONTINUE
  230        CONTINUE 
  220     CONTINUE
  210  CONTINUE
C
      DO 237 K = 1,NRHF
        WRITE(LUPRI,*) '(ALPHA I| K DELTA)'
        CALL OUTPUT(ALIKD(1,1,K),1,NBAS,1,NRHF,NBAS,NRHF,1,LUPRI)
  237 CONTINUE
C
C
      DO 310 K = 1,NRHF
         DO 320 I = 1,NRHF
            DO 330 A = NRHF+1,NRHF+NVIR !loop over occupied orbitals
               DO 340 ALPHA = 1,NBAS
C
                  AIKD(A,I,K) = AIKD(A,I,K) + ALIKD(ALPHA,I,K)
     *                                              *CMO(ALPHA,A)
C                  WRITE(LUPRI,*) 'K',K,'I',I,'A',A,'ALPHA',ALPHA
C                  WRITE(LUPRI,*)'CMO',CMO(ALPHA,A),'(AL I|K D)'
C     *    ,ALIKD(ALPHA,I,K) 
  340          CONTINUE 
               WRITE(LUPRI,*) 'AIKD',AIKD(A,I,K)
  330       CONTINUE
  320    CONTINUE
  310 CONTINUE

C      DO 337 K = 1,NRHF
C        WRITE(LUPRI,*) '(AI|K DELTA)'
C        CALL OUTPUT(AIKD(,I,K),1,NVIR,1,NRHF,NVIR,NRHF,1,LUPRI)
C  337 CONTINUE
C      DO 137 K = 1,NRHF
C        WRITE(LUPRI,*) '(AIK DELTA)'
C        CALL OUTPUT(AIKD(1,1,K),1,NVIR,1,NRHF,NVIR,NRHF,1,LUPRI)
C  137 CONTINUE


C
C      DO 410 K = 1,NRHF
C         DO 420 I = 1,NRHF
C            DO 430 A = NRHF+1,NRHF+NVIR !loop over occupied orbitals
C               DO 440 C = NRHF+1,NRHF+NVIR
CC
CC                              WRITE(LUPRI,*) 'CMO(GAM,A)',CMO(GAM,A)
CC                              WRITE(LUPRI,*) 'A',A,'I',I,'K',K,'GAM',GAM
CC                              WRITE(LUPRI,*) 'XINT',KIBD(K,I,GAM)
C                  AIKC(A,I,K) = AIKC(A,I,K) + AIKD(A,I,K)
C     *                                              *CMO(DELTA,C)
C                  WRITE(LUPRI,*) 'K',K,'I',I,'A',A,'DELTA',DELTA,'C',C
C                  WRITE(LUPRI,*) 'AIKC(A,I,K)',AIKC(A,I,K)
C  440          CONTINUE 
C  430       CONTINUE
C  420    CONTINUE
C  410 CONTINUE
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
C                        AIKCT2(A,I,K) = AIKCT2(A,I,K)
C     *                          + AIKC(A,I,K) * T2AMSQ(AI,BJ)
C
C  560              CONTINUE
C  550            CONTINUE
C  540           CONTINUE
C  530       CONTINUE
C  520   CONTINUE
C  510 CONTINUE
      RETURN
      END
