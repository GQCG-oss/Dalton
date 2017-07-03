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
      REAL*8 T2IJAD(NRHF*NVIR,NRHF*NBAS),T2IJGD(NRHF*NBAS,NRHF*NBAS)
      REAL*8 WORK(LWORK)
      INTEGER A,I,B,J,C,D,AI,BJ,DELTA,GAM,BETA,ALPHA,GA,LXINT,NRHF,NVIR
      INTEGER DELJ,DEL,GAMI
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------------------------------------------------------------
C      1) Transfrom occupied index in T2 ampltude to Delta.
C------------------------------------------------------------------------
C
      ICOUNT = 1
      ICOUN2 = 1
      WRITE(LUPRI,*) 'T2SQ'
      CALL OUTPUT(T2AMSQ,1,5,1,5,5,5,1,LUPRI)
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
                  IF (ICOUN2 .GT. NRHF*NVIR) ICOUN2=1
                  DELJ = DELTA + J - 1
C
                  T2IJAD(AI,DELJ) = T2IJAD(AI,DELJ)
     *                                      + T2AMSQ(AI,BJ)*CMO(DELTA,B)
C                  WRITE(LUPRI,*) 'T_(ij)^(c DELTA)',T2IJAD(AI,DELJ)
  140           CONTINUE
  130       CONTINUE
  120   CONTINUE
  110 CONTINUE

      ICOUN3 = 1
      ICOUN4 = 1
      DO 210 A = NRHF+1,NRHF+NVIR !Loop over virual orbitals
         DO 220 I = 1,NRHF !Loop over occupied orbitals
            AI = ICOUN3
            ICOUN3 = ICOUN3+1
            DO 230 B = NRHF+1,NRHF+NVIR !loop over virtual orbitals
               DO 240 J = 1,NRHF !Loop over occupied orbitals
                  DELJ =  J 
                     DO 250 GAM = 1,NBAS
                        GAMI = GAM + (I-1)*NBAS
                        T2IJGD(GAMI,DELJ) = T2IJGD(GAMI,DELJ)
     *                                    + T2IJAD(AI,DELJ)*CMO(GAM,A)
C
                              WRITE(LUPRI,*) 'T_(ij)^(gamma delta)'
                              WRITE(LUPRI,*) 'GAMI',GAMI,'DELJ'
     *                                      ,DELJ,T2IJGD(GAMI,DELJ)
  250              CONTINUE
  240           CONTINUE
  230       CONTINUE
  220   CONTINUE
  210 CONTINUE




C      DO 110 A = NRHF+1,NRHF+NVIR !Loop over virual orbitals
C         DO 120 I = 1,NRHF !Loop over occupied orbitals
C            AI = ICOUNT
C            ICOUNT = ICOUNT+1
C            DO 130 B = NRHF+1,NRHF+NVIR !loop over virtual orbitals
C               DO 140 J = 1,NRHF !Loop over occupied orbitals
C                  BJ = ICOUN2 
C                  ICOUN2 = ICOUN2+1
C                  IF (ICOUN2 .GT. NRHF*NVIR)  ICOUN2 = 1
C                  DO 150 C = NRHF+1,NRHF+NVIR
C                     DO 160 D = NRHF+1,NRHF+NVIR
C                        DO 170 GAM = 1,NBAS
C                           DO 180 BETA = 1,NBAS
C                              DO 190 ALPHA = 1,NBAS
CC   Transforming the b index of the amplitude to delta.
C                              T2IJAD(AI,BJ) = T2IJAD(AI,BJ)
C     *                                      + T2AMSQ(AI,BJ)*CMO(DELTA,B)
CC   Transforming the a index of the amplitude to gamma
C                              T2IJGD(AI,BJ) = T2IJGD(AI,BJ)
C     *                                      + T2IJAD(AI,BJ)*CMO(GAM,A)
CC                              WRITE(LUPRI,*) 'T_(ij)^(gamma delta)'
CC                              WRITE(LUPRI,*) 'A',A,'GAM',GAM,'AI',AI,
CC     *                                       'BJ',BJ,T2IJGD(AI,BJ)
CC Computing the the AGBD-intermediate: T_(ij)^(gamma delta) (alpha beta|gamma delta)
CC                              PQRDEL(ALPHA,BETA,GAM) =
CC     *                        PQRDEL(ALPHA,BETA,GAM) + T2IJGD(AI,BJ) *
CC     *                                              XINT(ALPHA,BETA,GAM)
CC Tranforming the AGBD-intermediate to MO-basis.
CC                               WORK  =
CC     *                         PQRDEL(ALPHA,BETA,GAM)*CMO(ALPHA,A)
CC     *                          *CMO(BETA,D)*CMO(GAM,B)*CMO(DELTA,D)
CC                               WRITE(LUPRI,*) 'WORK(LILLI)',WORK
C
C  190                   CONTINUE
C  180                 CONTINUE
C  170                CONTINUE
C  160              CONTINUE
C  150            CONTINUE
C                              WRITE(LUPRI,*) 'T_(ij)^(gamma delta)'
C                              WRITE(LUPRI,*) 'AI',AI,
C     *                                       'BJ',BJ,T2IJGD(AI,BJ)
C  140           CONTINUE
C  130       CONTINUE
C  120   CONTINUE
C  110 CONTINUE
CC
      FLUSH(LUPRI)
      RETURN
      END

