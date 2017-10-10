C   /* Deck check_agbd1 */
      SUBROUTINE CHECK_AGBD1(XINT,NBAS,CMO,NVIR,NRHF,T2AMSQ,LT2SQ,
     *                       T2IJAD,T2IJGD,AGBD,
     *                       WORK,LWORK,DELTA,LUPRI)
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
C       LUPRI:   Logical unit number of output file.
C
C   Written by Lilli Irene Ør Kristensen, Summer 2017.
C
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      REAL*8 XINT(NBAS,NBAS,NBAS),AGBD(NBAS,NRHF,NBAS*NRHF)
      REAL*8 CMO(NBAS,NRHF+NVIR),T2AMSQ(NRHF*NVIR,NRHF*NVIR)
      REAL*8 T2IJAD(NRHF*NVIR,NRHF*NBAS),T2IJGD(NRHF*NBAS,NRHF*NBAS)
      REAL*8 WORK(LWORK)!,NAIBEJ(NVIR,NRHF,NBAS*NRHF)
C      REAL*8 NAIBJ(NVIR*NRHF,NVIR*NRHF)
      INTEGER A,I,B,J,C,D,AI,BJ,DELTA,GAM,BETA,ALPHA,GA,LXINT,NRHF,NVIR
      INTEGER DELJ,DEL,GAMI,BETJ
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
      CALL DZERO(T2IJAD,NRHF*NRHF*NVIR*NVIR)
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
                  DELJ =  J 
C
                  T2IJAD(AI,DELJ) = T2IJAD(AI,DELJ)
     *                                      + T2AMSQ(AI,BJ)*CMO(DELTA,B)
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
C                         WRITE(LUPRI,*) 'T_(ij)^(gamma delta)'
c                         WRITE(LUPRI,*) 'GAMI',GAMI,'DELJ'
c     *                                      ,DELJ,T2IJGD(GAMI,DELJ)
c
c
  250              CONTINUE
  240           CONTINUE
  230       CONTINUE
  220   CONTINUE
  210 CONTINUE
C
      DO 310 BETA = 1,NBAS
         DO 320 J = 1,NRHF
            DELJ =  J
            BETJ = J + (BETA-1)*NRHF 
            DO 330 I  = 1,NRHF
               DO 340 GAM =1,NBAS
                  GAMI= GAM + (I-1)*NBAS
C                  IGAM= I + (GAM-1)*NRHF
                  DO 350 ALPHA = 1,NBAS
C
                     AGBD(ALPHA,I,BETJ) = AGBD(ALPHA,I,BETJ)
     *                          +XINT(ALPHA,GAM,BETA)*T2IJGD(GAMI,DELJ)
C
 350              CONTINUE
 340           CONTINUE
 330        CONTINUE
 320     CONTINUE
 310  CONTINUE
C
C      WRITE(LUPRI,*) 'AGBD'
C      DO 370 DELJ = 1,NRHF*NBAS
C        CALL OUTPUT(AGBD(1,1,DELJ),1,NBAS,1,NRHF,NBAS,NRHF,1,LUPRI)
C 370  CONTINUE
C
C
      RETURN
      END

