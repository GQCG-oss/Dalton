C   /* Deck check_afterloop */
      SUBROUTINE CHECK_AFTERLOOP(AGBD,AIBD,AIBJ,CMO,NBAS,NVIR,NRHF,
     *                   T2AMSQ,LT2SQ,
     *                     WORK,LWORK,DELTA,LUPRI)
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
      REAL*8 CMO(NBAS,NRHF+NVIR),T2AMSQ(NRHF*NVIR,NRHF*NVIR)
      REAL*8 AGBD(NBAS,NRHF,NRHF*NBAS),AIBD(NBAS,NRHF,NBAS*NVIR)
      REAL*8 AIBJ(NRHF*NVIR,NRHF*NVIR),NAIBJ(NRHF*NVIR,NRHF*NVIR)
      REAL*8 WORK(LWORK)
      INTEGER A,I,B,J,C,D,AI,BJ,DELTA,GAM,BETA,ALPHA,GA,LXINT,NRHF,NVIR
      INTEGER DELJ,DEL,GAMI,BETJ
C
      DO 370 DELJ = 1,NRHF*NBAS
      CALL OUTPUT(AGBD(1,1,DELJ),1,NBAS,1,NRHF,NBAS,NRHF,1,LUPRI)
 370  CONTINUE
       CALL DZERO(AIBD(1,1,1),NBAS*NRHF*NRHF*NVIR)
      DO 111 J = 1,NRHF
      DO 110 I = 1,NRHF
      DO 112 B = NRHF + 1, NRHF+NVIR 
      DO 119 BETA = 1,NBAS
      DO 118 ALPHA = 1,NBAS
      BETJ = J +(BETA-1)*NRHF
      BJ = J + (B-NRHF-1)*NRHF

      AIBD(ALPHA,I,BJ) = AIBD(ALPHA,I,BJ)
     *     + AGBD(ALPHA,I,BETJ)*CMO(BETA,B)

 118  CONTINUE
 119  CONTINUE
C      WRITE(LUPRI,*) 'AIBD','BJ',BJ, (AIBD(ALPHA,I,BJ),ALPHA=1,NBAS)
 112  CONTINUE
 110  CONTINUE
 111  CONTINUE

      DO 211 J = 1,NRHF
      DO 210 I = 1,NRHF
      DO 212 B = NRHF + 1, NRHF+NVIR 
      DO 217 A = NRHF +1, NRHF+NVIR
      AI = I + (A-NRHF-1)*NRHF
      BJ = J + (B-NRHF-1)*NRHF
      DO 218 ALPHA = 1,NBAS

      AIBJ(AI,BJ) = AIBJ(AI,BJ)
     *     + AIBD(ALPHA,I,BJ)*CMO(ALPHA,A)
C
 218  CONTINUE
C      WRITE(LUPRI,*) 'AIBJ',AI,BJ,AIBJ(AI,BJ)
 217  CONTINUE
 212  CONTINUE
 210  CONTINUE
 211  CONTINUE



      DO 317 A = NRHF +1, NRHF+NVIR
      DO 310 I = 1,NRHF
      DO 312 B = NRHF + 1, NRHF+NVIR 
      DO 311 J = 1,NRHF
      AI = I + (A-NRHF-1)*NRHF
      BJ = J + (B-NRHF-1)*NRHF
      NAI=AI
      NBJ=BJ
      IF (NAI .LT. NBJ) GO TO 311 

      NAIBJ(NAI,NBJ)= NAIBJ(NAI,NBJ) +AIBJ(AI,BJ) 
C
      WRITE(LUPRI,*) 'NAIBJ',NAI,NBJ,NAIBJ(NAI,NBJ)
 311  CONTINUE
 312  CONTINUE
 310  CONTINUE
 317  CONTINUE


      FLUSH(LUPRI)
      RETURN
      END
