C   /* Deck check_afterloop */
      SUBROUTINE CHECK_AFTERLOOP(AGBD,AIBD,AIBJ,KIAD,KIAC,
     *                     KTCKBJ,CMO,NBAS,NVIR,NRHF,
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
      REAL*8 KIAD(NRHF,NRHF,NVIR*NBAS),KIAC(NRHF,NRHF,NVIR*NVIR)
      REAL*8 TCKBJ(NRHF*NVIR,NRHF*NVIR)
      REAL*8 WORK(LWORK)
      INTEGER A,I,B,J,C,D,AI,BJ,DELTA,GAM,BETA,ALPHA,GA,LXINT,NRHF,NVIR
      INTEGER DELJ,DEL,GAMI,BETJ,AC,ADEL,CK
C
      CALL DZERO(TCKBJ(1,1),NRHF*NRHF*NVIR*NVIR) 
      CALL DZERO(AIBD(1,1,1),NBAS*NRHF*NRHF*NVIR)
C
C---------------------------------------------------------------
C   Transform the ALPHA and BETA index to A and B in 
C    sum_gamma( t_(ij)^(gamma delta)(alpha gamma| beta delta) )
C---------------------------------------------------------------
C
      DO 111 J = 1,NRHF
         DO 110 I = 1,NRHF
            DO 112 B = NRHF + 1, NRHF+NVIR 
               DO 119 BETA = 1,NBAS
                  DO 118 ALPHA = 1,NBAS
                     BETJ = J +(BETA-1)*NRHF
                     BJ = J + (B-NRHF-1)*NRHF

                     AIBD(ALPHA,I,BJ) = AIBD(ALPHA,I,BJ)
     *                   + AGBD(ALPHA,I,BETJ)*CMO(BETA,B)

 118              CONTINUE
 119           CONTINUE
C      WRITE(LUPRI,*) 'AIBD','BJ',BJ, (AIBD(ALPHA,I,BJ),ALPHA=1,NBAS)
 112        CONTINUE
 110     CONTINUE
 111  CONTINUE

      DO 211 J = 1,NRHF
         DO 210 I = 1,NRHF
            DO 212 B = NRHF + 1, NRHF+NVIR 
               DO 217 A = NRHF +1, NRHF+NVIR
                  AI = I + (A-NRHF-1)*NRHF
                  BJ = J + (B-NRHF-1)*NRHF
                  DO 218 ALPHA = 1,NBAS

                     AIBJ(AI,BJ) = AIBJ(AI,BJ)
     *                 + AIBD(ALPHA,I,BJ)*CMO(ALPHA,A)
C
 218              CONTINUE
C      WRITE(LUPRI,*) 'AIBJ',AI,BJ,AIBJ(AI,BJ)
 217           CONTINUE
 212        CONTINUE
 210     CONTINUE
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
C      WRITE(LUPRI,*) 'NAIBJ',NAI,NBJ,NAIBJ(NAI,NBJ)
 311           CONTINUE
 312        CONTINUE
 310     CONTINUE
 317  CONTINUE


C-----------------------------------------------------------------------
C   KIAC BIDRAGET
C----------------------------------------------------------------------

  
      DO 410 A = NRHF+1,NRHF+NVIR
      do 450 I = 1,NRHF
      DO 430 K = 1,NRHF
      DO 420 C = NRHF+1,NRHF+NVIR
      DO 440 DELTA = 1,NBAS !loop over occupied orbitals
               ADEL = (A-NRHF)+(DELTA-1)*NVIR
               AC = (A-NRHF) +(C-NRHF)*NVIR
C
        KIAC(K,I,AC) = KIAC(K,I,AC) + KIAD(K,I,ADEL) *CMO(DELTA,C)
C                  WRITE(LUPRI,*) 'K',K,'I',I,'A',A,'DELTA',DELTA,'C',C
C                  WRITE(LUPRI,*) K,I,A,DELTA
C                  WRITE(LUPRI,*) 'KIAD ', K,I,ADEL,KIAD(K,I,ADEL)
c                  WRITE(LUPRI,*) 'CMO  ', DELTA,C ,CMO(DELTA,C)
  440          CONTINUE 
                WRITE(LUPRI,*) 'KIAC',KIAC(K,I,AC)
  420       CONTINUE
  430 CONTINUE
  450 CONTINUE
  410    CONTINUE
C
C
      DO 510 A = NRHF+1,NRHF+NVIR !Loop over virual orbitals
      DO 520 I = 1,NRHF !Loop over occupied orbitals
      DO 530 B = NRHF+1,NRHF+NVIR !loop over virtual orbitals
      DO 540 J = 1,NRHF !Loop over occupied orbitals
      DO 550 K = 1,NRHF
      DO 560 C = NRHF+1,NRHF+NVIR
C 
      AI = (A-NRHF)+(I-1)*NVIR
      CK = (C-NRHF)+(K-1)*NVIR
      BJ = (B-NRHF)+(J-1)*NVIR
      AC = (A-NRHF) +(C-NRHF)*NVIR
C      KAIBJ = AI+(BJ-1)*NRHF*NVIR
C      KOFF = (C-NRHF)+(K-1)*NVIR+(A-NRHF-1)*NVIR*NRHF
C     *              +(I-1)*NVIR*NRHF*NVIR
C
      TCKBJ(BJ,AI) = TCKBJ(BJ,AI)
     *             + KIAC(K,I,AC) * T2AMSQ(CK,BJ)
C      IF ((BJ .EQ. 3) .AND. (AI .EQ. 1)) THEN
C      WRITE(LUPRI,*) ' CK ',CK,' BJ ',BJ,
c     *              ' K ',K,' I ',I,' AC ',AC
c      WRITE(LUPRI,*) 'KIAC(K,I,AC) ',KIAC(K,I,AC)
C      WRITE(LUPRI,*) 'T2AMSQ(CK,BJ)',T2AMSQ(CK,BJ)
C      WRITE(LUPRI,*) 'TCKBJ(AI,BJ) ',TCKBJ(BJ,AI)
C      END IF
  560              CONTINUE
  550            CONTINUE
  540           CONTINUE
  530       CONTINUE
  520   CONTINUE
  510 CONTINUE

C      DO 511 A = NRHF+1,NRHF+NVIR !Loop over virual orbitals
C      DO 521 I = 1,NRHF !Loop over occupied orbitals
C      DO 531 B = NRHF+1,NRHF+NVIR !loop over virtual orbitals
C      DO 541 J = 1,NRHF !Loop over occupied orbitals
C
C      AI = (A-NRHF)+(I-1)*NVIR
C      BJ = (B-NRHF)+(J-1)*NVIR
C
C      WRITE(LUPRI,*) 'TCKBJ(BJ,AI) ',BJ,AI,TCKBJ(BJ,AI)
C  541           CONTINUE
C  531       CONTINUE
C  521   CONTINUE
C  511 CONTINUE
C
C
C-----------------------------------------------------------------------
C   KIAC BIDRAGET
C----------------------------------------------------------------------










      FLUSH(LUPRI)
      RETURN
      END
