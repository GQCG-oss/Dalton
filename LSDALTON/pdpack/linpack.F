C  /* Deck dgelinpack */
#if !defined (VAR_NOLINPACK)
      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(*)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL A(LDA,*),Z(*)
      REAL RCOND
#else
      DOUBLE PRECISION A(LDA,*),Z(*)
      DOUBLE PRECISION RCOND
#endif
C
C     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,DSIGN
C
C     INTERNAL VARIABLES
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL DDOT,EK,T,WK,WKM
      REAL ANORM,S,DASUM,SM,YNORM
#else
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
#endif
      INTEGER INFO,J,K,KB,KP1,L
C
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 30
            S = ABS(A(K,K))/ABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 150
            S = ABS(A(K,K))/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(*),INFO
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL A(LDA,*)
#else
      DOUBLE PRECISION A(LDA,*)
#endif
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL T
#else
      DOUBLE PRECISION T
#endif
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(*),JOB
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL A(LDA,*),B(*)
#else
      DOUBLE PRECISION A(LDA,*),B(*)
#endif
C
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
C        OR DGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C
C     INTERNAL VARIABLES
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL DDOT,T
#else
      DOUBLE PRECISION DDOT,T
#endif
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(*),JOB
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL A(LDA,*),DET(2),WORK(*)
#else
      DOUBLE PRECISION A(LDA,*),DET(2),WORK(*)
#endif
C
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     DOUBLE PRECISION(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,DSWAP
C     FORTRAN DABS,MOD
C
C     INTERNAL VARIABLES
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
      REAL T
      REAL TEN
#else
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
#endif
      INTEGER I,J,K,KB,KP1,L,NM1
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
#endif

c$$$C  /* Deck dsplinpck */
c$$$#if !defined (VAR_NOLINPACK)
c$$$      SUBROUTINE DSPCO(AP,N,KPVT,RCOND,Z)
c$$$      INTEGER N,KPVT(*)
c$$$#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
c$$$      REAL AP(*),Z(*)
c$$$      REAL RCOND
c$$$#else
c$$$      DOUBLE PRECISION AP(*),Z(*)
c$$$      DOUBLE PRECISION RCOND
c$$$#endif
c$$$C
c$$$C     DSPCO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN
c$$$C     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES
c$$$C     THE CONDITION OF THE MATRIX.
c$$$C
c$$$C     IF  RCOND  IS NOT NEEDED, DSPFA IS SLIGHTLY FASTER.
c$$$C     TO SOLVE  A*X = B , FOLLOW DSPCO BY DSPSL.
c$$$C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPCO BY DSPSL.
c$$$C     TO COMPUTE  INVERSE(A) , FOLLOW DSPCO BY DSPDI.
c$$$C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPCO BY DSPDI.
c$$$C     TO COMPUTE  INERTIA(A), FOLLOW DSPCO BY DSPDI.
c$$$C
c$$$C     ON ENTRY
c$$$C
c$$$C        AP      DOUBLE PRECISION (N*(N+1)/2)
c$$$C                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
c$$$C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
c$$$C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
c$$$C                SEE COMMENTS BELOW FOR DETAILS.
c$$$C
c$$$C        N       INTEGER
c$$$C                THE ORDER OF THE MATRIX  A .
c$$$C
c$$$C     OUTPUT
c$$$C
c$$$C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
c$$$C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
c$$$C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
c$$$C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
c$$$C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
c$$$C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
c$$$C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
c$$$C
c$$$C        KPVT    INTEGER(N)
c$$$C                AN INTEGER VECTOR OF PIVOT INDICES.
c$$$C
c$$$C        RCOND   DOUBLE PRECISION
c$$$C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
c$$$C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
c$$$C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
c$$$C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
c$$$C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
c$$$C                           1.0 + RCOND .EQ. 1.0
c$$$C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
c$$$C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
c$$$C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
c$$$C                UNDERFLOWS.
c$$$C
c$$$C        Z       DOUBLE PRECISION(N)
c$$$C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
c$$$C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
c$$$C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
c$$$C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
c$$$C
c$$$C     PACKED STORAGE
c$$$C
c$$$C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
c$$$C          TRIANGLE OF A SYMMETRIC MATRIX.
c$$$C
c$$$C                K = 0
c$$$C                DO 20 J = 1, N
c$$$C                   DO 10 I = 1, J
c$$$C                      K = K + 1
c$$$C                      AP(K) = A(I,J)
c$$$C             10    CONTINUE
c$$$C             20 CONTINUE
c$$$C
c$$$C     LINPACK. THIS VERSION DATED 08/14/78 .
c$$$C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
c$$$C
c$$$C     SUBROUTINES AND FUNCTIONS
c$$$C
c$$$C     LINPACK DSPFA
c$$$C     BLAS DAXPY,DDOT,DSCAL,DASUM
c$$$C     FORTRAN ABS,DMAX1,IABS,DSIGN
c$$$C
c$$$C     INTERNAL VARIABLES
c$$$C
c$$$#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
c$$$      REAL AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T
c$$$      REAL ANORM,S,DASUM,YNORM
c$$$#else
c$$$      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T
c$$$      DOUBLE PRECISION ANORM,S,DASUM,YNORM
c$$$#endif
c$$$      INTEGER I,IJ,IK,IKM1,IKP1,INFO,J,JM1,J1
c$$$      INTEGER K,KK,KM1K,KM1KM1,KP,KPS,KS
c$$$C
c$$$C
c$$$C     FIND NORM OF A USING ONLY UPPER HALF
c$$$C
c$$$      J1 = 1
c$$$      DO 30 J = 1, N
c$$$         Z(J) = DASUM(J,AP(J1),1)
c$$$         IJ = J1
c$$$         J1 = J1 + J
c$$$         JM1 = J - 1
c$$$         IF (JM1 .LT. 1) GO TO 20
c$$$         DO 10 I = 1, JM1
c$$$            Z(I) = Z(I) + ABS(AP(IJ))
c$$$            IJ = IJ + 1
c$$$   10    CONTINUE
c$$$   20    CONTINUE
c$$$   30 CONTINUE
c$$$      ANORM = 0.0D0
c$$$      DO 40 J = 1, N
c$$$         ANORM = DMAX1(ANORM,Z(J))
c$$$   40 CONTINUE
c$$$C
c$$$C     FACTOR
c$$$C
c$$$      CALL DSPFA(AP,N,KPVT,INFO)
c$$$C
c$$$C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
c$$$C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
c$$$C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
c$$$C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
c$$$C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
c$$$C
c$$$C     SOLVE U*D*W = E
c$$$C
c$$$      EK = 1.0D0
c$$$      DO 50 J = 1, N
c$$$         Z(J) = 0.0D0
c$$$   50 CONTINUE
c$$$      K = N
c$$$      IK = (N*(N - 1))/2
c$$$   60 IF (K .EQ. 0) GO TO 120
c$$$         KK = IK + K
c$$$         IKM1 = IK - (K - 1)
c$$$         KS = 1
c$$$         IF (KPVT(K) .LT. 0) KS = 2
c$$$         KP = IABS(KPVT(K))
c$$$         KPS = K + 1 - KS
c$$$         IF (KP .EQ. KPS) GO TO 70
c$$$            T = Z(KPS)
c$$$            Z(KPS) = Z(KP)
c$$$            Z(KP) = T
c$$$   70    CONTINUE
c$$$         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))
c$$$         Z(K) = Z(K) + EK
c$$$         CALL DAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
c$$$         IF (KS .EQ. 1) GO TO 80
c$$$            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))
c$$$            Z(K-1) = Z(K-1) + EK
c$$$            CALL DAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
c$$$   80    CONTINUE
c$$$         IF (KS .EQ. 2) GO TO 100
c$$$            IF (ABS(Z(K)) .LE. ABS(AP(KK))) GO TO 90
c$$$               S = ABS(AP(KK))/ABS(Z(K))
c$$$               CALL DSCAL(N,S,Z,1)
c$$$               EK = S*EK
c$$$   90       CONTINUE
c$$$            IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)
c$$$            IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0
c$$$         GO TO 110
c$$$  100    CONTINUE
c$$$            KM1K = IK + K - 1
c$$$            KM1KM1 = IKM1 + K - 1
c$$$            AK = AP(KK)/AP(KM1K)
c$$$            AKM1 = AP(KM1KM1)/AP(KM1K)
c$$$            BK = Z(K)/AP(KM1K)
c$$$            BKM1 = Z(K-1)/AP(KM1K)
c$$$            DENOM = AK*AKM1 - 1.0D0
c$$$            Z(K) = (AKM1*BK - BKM1)/DENOM
c$$$            Z(K-1) = (AK*BKM1 - BK)/DENOM
c$$$  110    CONTINUE
c$$$         K = K - KS
c$$$         IK = IK - K
c$$$         IF (KS .EQ. 2) IK = IK - (K + 1)
c$$$      GO TO 60
c$$$  120 CONTINUE
c$$$      S = 1.0D0/DASUM(N,Z,1)
c$$$      CALL DSCAL(N,S,Z,1)
c$$$C
c$$$C     SOLVE TRANS(U)*Y = W
c$$$C
c$$$      K = 1
c$$$      IK = 0
c$$$  130 IF (K .GT. N) GO TO 160
c$$$         KS = 1
c$$$         IF (KPVT(K) .LT. 0) KS = 2
c$$$         IF (K .EQ. 1) GO TO 150
c$$$            Z(K) = Z(K) + DDOT(K-1,AP(IK+1),1,Z(1),1)
c$$$            IKP1 = IK + K
c$$$            IF (KS .EQ. 2)
c$$$     *         Z(K+1) = Z(K+1) + DDOT(K-1,AP(IKP1+1),1,Z(1),1)
c$$$            KP = IABS(KPVT(K))
c$$$            IF (KP .EQ. K) GO TO 140
c$$$               T = Z(K)
c$$$               Z(K) = Z(KP)
c$$$               Z(KP) = T
c$$$  140       CONTINUE
c$$$  150    CONTINUE
c$$$         IK = IK + K
c$$$         IF (KS .EQ. 2) IK = IK + (K + 1)
c$$$         K = K + KS
c$$$      GO TO 130
c$$$  160 CONTINUE
c$$$      S = 1.0D0/DASUM(N,Z,1)
c$$$      CALL DSCAL(N,S,Z,1)
c$$$C
c$$$      YNORM = 1.0D0
c$$$C
c$$$C     SOLVE U*D*V = Y
c$$$C
c$$$      K = N
c$$$      IK = N*(N - 1)/2
c$$$  170 IF (K .EQ. 0) GO TO 230
c$$$         KK = IK + K
c$$$         IKM1 = IK - (K - 1)
c$$$         KS = 1
c$$$         IF (KPVT(K) .LT. 0) KS = 2
c$$$         IF (K .EQ. KS) GO TO 190
c$$$            KP = IABS(KPVT(K))
c$$$            KPS = K + 1 - KS
c$$$            IF (KP .EQ. KPS) GO TO 180
c$$$               T = Z(KPS)
c$$$               Z(KPS) = Z(KP)
c$$$               Z(KP) = T
c$$$  180       CONTINUE
c$$$            CALL DAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
c$$$            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
c$$$  190    CONTINUE
c$$$         IF (KS .EQ. 2) GO TO 210
c$$$            IF (ABS(Z(K)) .LE. ABS(AP(KK))) GO TO 200
c$$$               S = ABS(AP(KK))/ABS(Z(K))
c$$$               CALL DSCAL(N,S,Z,1)
c$$$               YNORM = S*YNORM
c$$$  200       CONTINUE
c$$$            IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)
c$$$            IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0
c$$$         GO TO 220
c$$$  210    CONTINUE
c$$$            KM1K = IK + K - 1
c$$$            KM1KM1 = IKM1 + K - 1
c$$$            AK = AP(KK)/AP(KM1K)
c$$$            AKM1 = AP(KM1KM1)/AP(KM1K)
c$$$            BK = Z(K)/AP(KM1K)
c$$$            BKM1 = Z(K-1)/AP(KM1K)
c$$$            DENOM = AK*AKM1 - 1.0D0
c$$$            Z(K) = (AKM1*BK - BKM1)/DENOM
c$$$            Z(K-1) = (AK*BKM1 - BK)/DENOM
c$$$  220    CONTINUE
c$$$         K = K - KS
c$$$         IK = IK - K
c$$$         IF (KS .EQ. 2) IK = IK - (K + 1)
c$$$      GO TO 170
c$$$  230 CONTINUE
c$$$      S = 1.0D0/DASUM(N,Z,1)
c$$$      CALL DSCAL(N,S,Z,1)
c$$$      YNORM = S*YNORM
c$$$C
c$$$C     SOLVE TRANS(U)*Z = V
c$$$C
c$$$      K = 1
c$$$      IK = 0
c$$$  240 IF (K .GT. N) GO TO 270
c$$$         KS = 1
c$$$         IF (KPVT(K) .LT. 0) KS = 2
c$$$         IF (K .EQ. 1) GO TO 260
c$$$            Z(K) = Z(K) + DDOT(K-1,AP(IK+1),1,Z(1),1)
c$$$            IKP1 = IK + K
c$$$            IF (KS .EQ. 2)
c$$$     *         Z(K+1) = Z(K+1) + DDOT(K-1,AP(IKP1+1),1,Z(1),1)
c$$$            KP = IABS(KPVT(K))
c$$$            IF (KP .EQ. K) GO TO 250
c$$$               T = Z(K)
c$$$               Z(K) = Z(KP)
c$$$               Z(KP) = T
c$$$  250       CONTINUE
c$$$  260    CONTINUE
c$$$         IK = IK + K
c$$$         IF (KS .EQ. 2) IK = IK + (K + 1)
c$$$         K = K + KS
c$$$      GO TO 240
c$$$  270 CONTINUE
c$$$C     MAKE ZNORM = 1.0
c$$$      S = 1.0D0/DASUM(N,Z,1)
c$$$      CALL DSCAL(N,S,Z,1)
c$$$      YNORM = S*YNORM
c$$$C
c$$$      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
c$$$      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE DSPFA(AP,N,KPVT,INFO)
c$$$      INTEGER N,KPVT(*),INFO
c$$$#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
c$$$      REAL AP(*)
c$$$#else
c$$$      DOUBLE PRECISION AP(*)
c$$$#endif
c$$$C
c$$$C     DSPFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN
c$$$C     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.
c$$$C
c$$$C     TO SOLVE  A*X = B , FOLLOW DSPFA BY DSPSL.
c$$$C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPFA BY DSPSL.
c$$$C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPFA BY DSPDI.
c$$$C     TO COMPUTE  INERTIA(A) , FOLLOW DSPFA BY DSPDI.
c$$$C     TO COMPUTE  INVERSE(A) , FOLLOW DSPFA BY DSPDI.
c$$$C
c$$$C     ON ENTRY
c$$$C
c$$$C        AP      DOUBLE PRECISION (N*(N+1)/2)
c$$$C                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
c$$$C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
c$$$C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
c$$$C                SEE COMMENTS BELOW FOR DETAILS.
c$$$C
c$$$C        N       INTEGER
c$$$C                THE ORDER OF THE MATRIX  A .
c$$$C
c$$$C     OUTPUT
c$$$C
c$$$C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
c$$$C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
c$$$C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
c$$$C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
c$$$C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
c$$$C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
c$$$C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
c$$$C
c$$$C        KPVT    INTEGER(N)
c$$$C                AN INTEGER VECTOR OF PIVOT INDICES.
c$$$C
c$$$C        INFO    INTEGER
c$$$C                = 0  NORMAL VALUE.
c$$$C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
c$$$C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
c$$$C                     BUT IT DOES INDICATE THAT DSPSL OR DSPDI MAY
c$$$C                     DIVIDE BY ZERO IF CALLED.
c$$$C
c$$$C     PACKED STORAGE
c$$$C
c$$$C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
c$$$C          TRIANGLE OF A SYMMETRIC MATRIX.
c$$$C
c$$$C                K = 0
c$$$C                DO 20 J = 1, N
c$$$C                   DO 10 I = 1, J
c$$$C                      K = K + 1
c$$$C                      AP(K)  = A(I,J)
c$$$C             10    CONTINUE
c$$$C             20 CONTINUE
c$$$C
c$$$C     LINPACK. THIS VERSION DATED 08/14/78 .
c$$$C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
c$$$C
c$$$C     SUBROUTINES AND FUNCTIONS
c$$$C
c$$$C     BLAS DAXPY,DSWAP,IDAMAX
c$$$C     FORTRAN DABS,DMAX1,DSQRT
c$$$C
c$$$C     INTERNAL VARIABLES
c$$$C
c$$$#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
c$$$      REAL AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
c$$$      REAL ABSAKK,ALPHA,COLMAX,ROWMAX
c$$$#else
c$$$      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
c$$$      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
c$$$#endif
c$$$      INTEGER IDAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
c$$$      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
c$$$      LOGICAL SWAP
c$$$C
c$$$C
c$$$C     INITIALIZE
c$$$C
c$$$C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
c$$$      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
c$$$C
c$$$      INFO = 0
c$$$C
c$$$C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
c$$$C
c$$$      K = N
c$$$      IK = (N*(N - 1))/2
c$$$   10 CONTINUE
c$$$C
c$$$C        LEAVE THE LOOP IF K=0 OR K=1.
c$$$C
c$$$C     ...EXIT
c$$$         IF (K .EQ. 0) GO TO 200
c$$$         IF (K .GT. 1) GO TO 20
c$$$            KPVT(1) = 1
c$$$            IF (AP(1) .EQ. 0.0D0) INFO = 1
c$$$C     ......EXIT
c$$$            GO TO 200
c$$$   20    CONTINUE
c$$$C
c$$$C        THIS SECTION OF CODE DETERMINES THE KIND OF
c$$$C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
c$$$C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
c$$$C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
c$$$C        REQUIRED.
c$$$C
c$$$         KM1 = K - 1
c$$$         KK = IK + K
c$$$         ABSAKK = ABS(AP(KK))
c$$$C
c$$$C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
c$$$C        COLUMN K.
c$$$C
c$$$         IMAX = IDAMAX(K-1,AP(IK+1),1)
c$$$         IMK = IK + IMAX
c$$$         COLMAX = ABS(AP(IMK))
c$$$         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
c$$$            KSTEP = 1
c$$$            SWAP = .FALSE.
c$$$         GO TO 90
c$$$   30    CONTINUE
c$$$C
c$$$C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
c$$$C           ROW IMAX.
c$$$C
c$$$            ROWMAX = 0.0D0
c$$$            IMAXP1 = IMAX + 1
c$$$            IM = IMAX*(IMAX - 1)/2
c$$$            IMJ = IM + 2*IMAX
c$$$            DO 40 J = IMAXP1, K
c$$$               ROWMAX = DMAX1(ROWMAX,ABS(AP(IMJ)))
c$$$               IMJ = IMJ + J
c$$$   40       CONTINUE
c$$$            IF (IMAX .EQ. 1) GO TO 50
c$$$               JMAX = IDAMAX(IMAX-1,AP(IM+1),1)
c$$$               JMIM = JMAX + IM
c$$$               ROWMAX = DMAX1(ROWMAX,ABS(AP(JMIM)))
c$$$   50       CONTINUE
c$$$            IMIM = IMAX + IM
c$$$            IF (ABS(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
c$$$               KSTEP = 1
c$$$               SWAP = .TRUE.
c$$$            GO TO 80
c$$$   60       CONTINUE
c$$$            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
c$$$               KSTEP = 1
c$$$               SWAP = .FALSE.
c$$$            GO TO 80
c$$$   70       CONTINUE
c$$$               KSTEP = 2
c$$$               SWAP = IMAX .NE. KM1
c$$$   80       CONTINUE
c$$$   90    CONTINUE
c$$$         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
c$$$C
c$$$C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
c$$$C
c$$$            KPVT(K) = K
c$$$            INFO = K
c$$$         GO TO 190
c$$$  100    CONTINUE
c$$$         IF (KSTEP .EQ. 2) GO TO 140
c$$$C
c$$$C           1 X 1 PIVOT BLOCK.
c$$$C
c$$$            IF (.NOT.SWAP) GO TO 120
c$$$C
c$$$C              PERFORM AN INTERCHANGE.
c$$$C
c$$$               CALL DSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
c$$$               IMJ = IK + IMAX
c$$$               DO 110 JJ = IMAX, K
c$$$                  J = K + IMAX - JJ
c$$$                  JK = IK + J
c$$$                  T = AP(JK)
c$$$                  AP(JK) = AP(IMJ)
c$$$                  AP(IMJ) = T
c$$$                  IMJ = IMJ - (J - 1)
c$$$  110          CONTINUE
c$$$  120       CONTINUE
c$$$C
c$$$C           PERFORM THE ELIMINATION.
c$$$C
c$$$            IJ = IK - (K - 1)
c$$$            DO 130 JJ = 1, KM1
c$$$               J = K - JJ
c$$$               JK = IK + J
c$$$               MULK = -AP(JK)/AP(KK)
c$$$               T = MULK
c$$$               CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
c$$$               IJJ = IJ + J
c$$$               AP(JK) = MULK
c$$$               IJ = IJ - (J - 1)
c$$$  130       CONTINUE
c$$$C
c$$$C           SET THE PIVOT ARRAY.
c$$$C
c$$$            KPVT(K) = K
c$$$            IF (SWAP) KPVT(K) = IMAX
c$$$         GO TO 190
c$$$  140    CONTINUE
c$$$C
c$$$C           2 X 2 PIVOT BLOCK.
c$$$C
c$$$            KM1K = IK + K - 1
c$$$            IKM1 = IK - (K - 1)
c$$$            IF (.NOT.SWAP) GO TO 160
c$$$C
c$$$C              PERFORM AN INTERCHANGE.
c$$$C
c$$$               CALL DSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
c$$$               IMJ = IKM1 + IMAX
c$$$               DO 150 JJ = IMAX, KM1
c$$$                  J = KM1 + IMAX - JJ
c$$$                  JKM1 = IKM1 + J
c$$$                  T = AP(JKM1)
c$$$                  AP(JKM1) = AP(IMJ)
c$$$                  AP(IMJ) = T
c$$$                  IMJ = IMJ - (J - 1)
c$$$  150          CONTINUE
c$$$               T = AP(KM1K)
c$$$               AP(KM1K) = AP(IMK)
c$$$               AP(IMK) = T
c$$$  160       CONTINUE
c$$$C
c$$$C           PERFORM THE ELIMINATION.
c$$$C
c$$$            KM2 = K - 2
c$$$            IF (KM2 .EQ. 0) GO TO 180
c$$$               AK = AP(KK)/AP(KM1K)
c$$$               KM1KM1 = IKM1 + K - 1
c$$$               AKM1 = AP(KM1KM1)/AP(KM1K)
c$$$               DENOM = 1.0D0 - AK*AKM1
c$$$               IJ = IK - (K - 1) - (K - 2)
c$$$               DO 170 JJ = 1, KM2
c$$$                  J = KM1 - JJ
c$$$                  JK = IK + J
c$$$                  BK = AP(JK)/AP(KM1K)
c$$$                  JKM1 = IKM1 + J
c$$$                  BKM1 = AP(JKM1)/AP(KM1K)
c$$$                  MULK = (AKM1*BK - BKM1)/DENOM
c$$$                  MULKM1 = (AK*BKM1 - BK)/DENOM
c$$$                  T = MULK
c$$$                  CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
c$$$                  T = MULKM1
c$$$                  CALL DAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
c$$$                  AP(JK) = MULK
c$$$                  AP(JKM1) = MULKM1
c$$$                  IJJ = IJ + J
c$$$                  IJ = IJ - (J - 1)
c$$$  170          CONTINUE
c$$$  180       CONTINUE
c$$$C
c$$$C           SET THE PIVOT ARRAY.
c$$$C
c$$$            KPVT(K) = 1 - K
c$$$            IF (SWAP) KPVT(K) = -IMAX
c$$$            KPVT(K-1) = KPVT(K)
c$$$  190    CONTINUE
c$$$         IK = IK - (K - 1)
c$$$         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
c$$$         K = K - KSTEP
c$$$      GO TO 10
c$$$  200 CONTINUE
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE DSPSL(AP,N,KPVT,B)
c$$$      INTEGER N,KPVT(*)
c$$$#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
c$$$      REAL AP(*),B(*)
c$$$#else
c$$$      DOUBLE PRECISION AP(*),B(*)
c$$$#endif
c$$$C
c$$$C     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM
c$$$C     A * X = B
c$$$C     USING THE FACTORS COMPUTED BY DSPFA.
c$$$C
c$$$C     ON ENTRY
c$$$C
c$$$C        AP      DOUBLE PRECISION(N*(N+1)/2)
c$$$C                THE OUTPUT FROM DSPFA.
c$$$C
c$$$C        N       INTEGER
c$$$C                THE ORDER OF THE MATRIX  A .
c$$$C
c$$$C        KPVT    INTEGER(N)
c$$$C                THE PIVOT VECTOR FROM DSPFA.
c$$$C
c$$$C        B       DOUBLE PRECISION(N)
c$$$C                THE RIGHT HAND SIDE VECTOR.
c$$$C
c$$$C     ON RETURN
c$$$C
c$$$C        B       THE SOLUTION VECTOR  X .
c$$$C
c$$$C     ERROR CONDITION
c$$$C
c$$$C        A DIVISION BY ZERO MAY OCCUR IF  DSPCO  HAS SET RCOND .EQ. 0.0
c$$$C        OR  DSPFA  HAS SET INFO .NE. 0  .
c$$$C
c$$$C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
c$$$C     WITH  P  COLUMNS
c$$$C           CALL DSPFA(AP,N,KPVT,INFO)
c$$$C           IF (INFO .NE. 0) GO TO ...
c$$$C           DO 10 J = 1, P
c$$$C              CALL DSPSL(AP,N,KPVT,C(1,J))
c$$$C        10 CONTINUE
c$$$C
c$$$C     LINPACK. THIS VERSION DATED 08/14/78 .
c$$$C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
c$$$C
c$$$C     SUBROUTINES AND FUNCTIONS
c$$$C
c$$$C     BLAS DAXPY,DDOT
c$$$C     FORTRAN IABS
c$$$C
c$$$C     INTERNAL VARIABLES.
c$$$C
c$$$#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
c$$$      REAL AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP
c$$$#else
c$$$      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP
c$$$#endif
c$$$      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
c$$$C
c$$$C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
c$$$C     D INVERSE TO B.
c$$$C
c$$$      K = N
c$$$      IK = (N*(N - 1))/2
c$$$   10 IF (K .EQ. 0) GO TO 80
c$$$         KK = IK + K
c$$$         IF (KPVT(K) .LT. 0) GO TO 40
c$$$C
c$$$C           1 X 1 PIVOT BLOCK.
c$$$C
c$$$            IF (K .EQ. 1) GO TO 30
c$$$               KP = KPVT(K)
c$$$               IF (KP .EQ. K) GO TO 20
c$$$C
c$$$C                 INTERCHANGE.
c$$$C
c$$$                  TEMP = B(K)
c$$$                  B(K) = B(KP)
c$$$                  B(KP) = TEMP
c$$$   20          CONTINUE
c$$$C
c$$$C              APPLY THE TRANSFORMATION.
c$$$C
c$$$               CALL DAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
c$$$   30       CONTINUE
c$$$C
c$$$C           APPLY D INVERSE.
c$$$C
c$$$            B(K) = B(K)/AP(KK)
c$$$            K = K - 1
c$$$            IK = IK - K
c$$$         GO TO 70
c$$$   40    CONTINUE
c$$$C
c$$$C           2 X 2 PIVOT BLOCK.
c$$$C
c$$$            IKM1 = IK - (K - 1)
c$$$            IF (K .EQ. 2) GO TO 60
c$$$               KP = IABS(KPVT(K))
c$$$               IF (KP .EQ. K - 1) GO TO 50
c$$$C
c$$$C                 INTERCHANGE.
c$$$C
c$$$                  TEMP = B(K-1)
c$$$                  B(K-1) = B(KP)
c$$$                  B(KP) = TEMP
c$$$   50          CONTINUE
c$$$C
c$$$C              APPLY THE TRANSFORMATION.
c$$$C
c$$$               CALL DAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
c$$$               CALL DAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
c$$$   60       CONTINUE
c$$$C
c$$$C           APPLY D INVERSE.
c$$$C
c$$$            KM1K = IK + K - 1
c$$$            KK = IK + K
c$$$            AK = AP(KK)/AP(KM1K)
c$$$            KM1KM1 = IKM1 + K - 1
c$$$            AKM1 = AP(KM1KM1)/AP(KM1K)
c$$$            BK = B(K)/AP(KM1K)
c$$$            BKM1 = B(K-1)/AP(KM1K)
c$$$            DENOM = AK*AKM1 - 1.0D0
c$$$            B(K) = (AKM1*BK - BKM1)/DENOM
c$$$            B(K-1) = (AK*BKM1 - BK)/DENOM
c$$$            K = K - 2
c$$$            IK = IK - (K + 1) - K
c$$$   70    CONTINUE
c$$$      GO TO 10
c$$$   80 CONTINUE
c$$$C
c$$$C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
c$$$C
c$$$      K = 1
c$$$      IK = 0
c$$$   90 IF (K .GT. N) GO TO 160
c$$$         IF (KPVT(K) .LT. 0) GO TO 120
c$$$C
c$$$C           1 X 1 PIVOT BLOCK.
c$$$C
c$$$            IF (K .EQ. 1) GO TO 110
c$$$C
c$$$C              APPLY THE TRANSFORMATION.
c$$$C
c$$$               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
c$$$               KP = KPVT(K)
c$$$               IF (KP .EQ. K) GO TO 100
c$$$C
c$$$C                 INTERCHANGE.
c$$$C
c$$$                  TEMP = B(K)
c$$$                  B(K) = B(KP)
c$$$                  B(KP) = TEMP
c$$$  100          CONTINUE
c$$$  110       CONTINUE
c$$$            IK = IK + K
c$$$            K = K + 1
c$$$         GO TO 150
c$$$  120    CONTINUE
c$$$C
c$$$C           2 X 2 PIVOT BLOCK.
c$$$C
c$$$            IF (K .EQ. 1) GO TO 140
c$$$C
c$$$C              APPLY THE TRANSFORMATION.
c$$$C
c$$$               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
c$$$               IKP1 = IK + K
c$$$               B(K+1) = B(K+1) + DDOT(K-1,AP(IKP1+1),1,B(1),1)
c$$$               KP = IABS(KPVT(K))
c$$$               IF (KP .EQ. K) GO TO 130
c$$$C
c$$$C                 INTERCHANGE.
c$$$C
c$$$                  TEMP = B(K)
c$$$                  B(K) = B(KP)
c$$$                  B(KP) = TEMP
c$$$  130          CONTINUE
c$$$  140       CONTINUE
c$$$            IK = IK + K + K + 1
c$$$            K = K + 2
c$$$  150    CONTINUE
c$$$      GO TO 90
c$$$  160 CONTINUE
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE DSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
c$$$      INTEGER N,JOB
c$$$#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
c$$$      REAL AP(*),WORK(*)
c$$$      REAL DET(2)
c$$$#else
c$$$      DOUBLE PRECISION AP(*),WORK(*)
c$$$      DOUBLE PRECISION DET(2)
c$$$#endif
c$$$      INTEGER KPVT(*),INERT(3)
c$$$C
c$$$C     DSPDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
c$$$C     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM
c$$$C     DSPFA, WHERE THE MATRIX IS STORED IN PACKED FORM.
c$$$C
c$$$C     ON ENTRY
c$$$C
c$$$C        AP      DOUBLE PRECISION (N*(N+1)/2)
c$$$C                THE OUTPUT FROM DSPFA.
c$$$C
c$$$C        N       INTEGER
c$$$C                THE ORDER OF THE MATRIX A.
c$$$C
c$$$C        KPVT    INTEGER(N)
c$$$C                THE PIVOT VECTOR FROM DSPFA.
c$$$C
c$$$C        WORK    DOUBLE PRECISION(N)
c$$$C                WORK VECTOR.  CONTENTS IGNORED.
c$$$C
c$$$C        JOB     INTEGER
c$$$C                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE
c$$$C                   IF  C .NE. 0, THE INVERSE IS COMPUTED,
c$$$C                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,
c$$$C                   IF  A .NE. 0, THE INERTIA IS COMPUTED.
c$$$C
c$$$C                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.
c$$$C
c$$$C     ON RETURN
c$$$C
c$$$C        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
c$$$C
c$$$C        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
c$$$C               THE ORIGINAL MATRIX, STORED IN PACKED FORM.
c$$$C               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED
c$$$C               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY.
c$$$C
c$$$C        DET    DOUBLE PRECISION(2)
c$$$C               DETERMINANT OF ORIGINAL MATRIX.
c$$$C               DETERMINANT = DET(1) * 10.0**DET(2)
c$$$C               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0
c$$$C               OR DET(1) = 0.0.
c$$$C
c$$$C        INERT  INTEGER(3)
c$$$C               THE INERTIA OF THE ORIGINAL MATRIX.
c$$$C               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
c$$$C               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
c$$$C               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
c$$$C
c$$$C     ERROR CONDITION
c$$$C
c$$$C        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED
c$$$C        AND  DSPCO  HAS SET RCOND .EQ. 0.0
c$$$C        OR  DSPFA  HAS SET  INFO .NE. 0 .
c$$$C
c$$$C     LINPACK. THIS VERSION DATED 08/14/78 .
c$$$C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
c$$$C
c$$$C     SUBROUTINES AND FUNCTIONS
c$$$C
c$$$C     BLAS DAXPY,DCOPY,DDOT,DSWAP
c$$$C     FORTRAN DABS,IABS,MOD
c$$$C
c$$$C     INTERNAL VARIABLES.
c$$$C
c$$$#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_REAL)
c$$$      REAL AKKP1,DDOT,TEMP
c$$$      REAL TEN,D,T,AK,AKP1
c$$$#else
c$$$      DOUBLE PRECISION AKKP1,DDOT,TEMP
c$$$      DOUBLE PRECISION TEN,D,T,AK,AKP1
c$$$#endif
c$$$      INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
c$$$      INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
c$$$      LOGICAL NOINV,NODET,NOERT
c$$$C
c$$$      NOINV = MOD(JOB,10) .EQ. 0
c$$$      NODET = MOD(JOB,100)/10 .EQ. 0
c$$$      NOERT = MOD(JOB,1000)/100 .EQ. 0
c$$$C
c$$$      IF (NODET .AND. NOERT) GO TO 140
c$$$         IF (NOERT) GO TO 10
c$$$            INERT(1) = 0
c$$$            INERT(2) = 0
c$$$            INERT(3) = 0
c$$$   10    CONTINUE
c$$$         IF (NODET) GO TO 20
c$$$            DET(1) = 1.0D0
c$$$            DET(2) = 0.0D0
c$$$            TEN = 10.0D0
c$$$   20    CONTINUE
c$$$         T = 0.0D0
c$$$         IK = 0
c$$$         DO 130 K = 1, N
c$$$            KK = IK + K
c$$$            D = AP(KK)
c$$$C
c$$$C           CHECK IF 1 BY 1
c$$$C
c$$$            IF (KPVT(K) .GT. 0) GO TO 50
c$$$C
c$$$C              2 BY 2 BLOCK
c$$$C              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)
c$$$C                      (S  C)
c$$$C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
c$$$C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
c$$$C
c$$$               IF (T .NE. 0.0D0) GO TO 30
c$$$                  IKP1 = IK + K
c$$$                  KKP1 = IKP1 + K
c$$$                  T = ABS(AP(KKP1))
c$$$                  D = (D/T)*AP(KKP1+1) - T
c$$$               GO TO 40
c$$$   30          CONTINUE
c$$$                  D = T
c$$$                  T = 0.0D0
c$$$   40          CONTINUE
c$$$   50       CONTINUE
c$$$C
c$$$            IF (NOERT) GO TO 60
c$$$               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1
c$$$               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1
c$$$               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1
c$$$   60       CONTINUE
c$$$C
c$$$            IF (NODET) GO TO 120
c$$$               DET(1) = D*DET(1)
c$$$               IF (DET(1) .EQ. 0.0D0) GO TO 110
c$$$   70             IF (ABS(DET(1)) .GE. 1.0D0) GO TO 80
c$$$                     DET(1) = TEN*DET(1)
c$$$                     DET(2) = DET(2) - 1.0D0
c$$$                  GO TO 70
c$$$   80             CONTINUE
c$$$   90             IF (ABS(DET(1)) .LT. TEN) GO TO 100
c$$$                     DET(1) = DET(1)/TEN
c$$$                     DET(2) = DET(2) + 1.0D0
c$$$                  GO TO 90
c$$$  100             CONTINUE
c$$$  110          CONTINUE
c$$$  120       CONTINUE
c$$$            IK = IK + K
c$$$  130    CONTINUE
c$$$  140 CONTINUE
c$$$C
c$$$C     COMPUTE INVERSE(A)
c$$$C
c$$$      IF (NOINV) GO TO 270
c$$$         K = 1
c$$$         IK = 0
c$$$  150    IF (K .GT. N) GO TO 260
c$$$            KM1 = K - 1
c$$$            KK = IK + K
c$$$            IKP1 = IK + K
c$$$            KKP1 = IKP1 + K
c$$$            IF (KPVT(K) .LT. 0) GO TO 180
c$$$C
c$$$C              1 BY 1
c$$$C
c$$$               AP(KK) = 1.0D0/AP(KK)
c$$$               IF (KM1 .LT. 1) GO TO 170
c$$$                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)
c$$$                  IJ = 0
c$$$                  DO 160 J = 1, KM1
c$$$                     JK = IK + J
c$$$                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
c$$$                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
c$$$                     IJ = IJ + J
c$$$  160             CONTINUE
c$$$                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
c$$$  170          CONTINUE
c$$$               KSTEP = 1
c$$$            GO TO 220
c$$$  180       CONTINUE
c$$$C
c$$$C              2 BY 2
c$$$C
c$$$               T = ABS(AP(KKP1))
c$$$               AK = AP(KK)/T
c$$$               AKP1 = AP(KKP1+1)/T
c$$$               AKKP1 = AP(KKP1)/T
c$$$               D = T*(AK*AKP1 - 1.0D0)
c$$$               AP(KK) = AKP1/D
c$$$               AP(KKP1+1) = AK/D
c$$$               AP(KKP1) = -AKKP1/D
c$$$               IF (KM1 .LT. 1) GO TO 210
c$$$                  CALL DCOPY(KM1,AP(IKP1+1),1,WORK,1)
c$$$                  IJ = 0
c$$$                  DO 190 J = 1, KM1
c$$$                     JKP1 = IKP1 + J
c$$$                     AP(JKP1) = DDOT(J,AP(IJ+1),1,WORK,1)
c$$$                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
c$$$                     IJ = IJ + J
c$$$  190             CONTINUE
c$$$                  AP(KKP1+1) = AP(KKP1+1)
c$$$     *                         + DDOT(KM1,WORK,1,AP(IKP1+1),1)
c$$$                  AP(KKP1) = AP(KKP1)
c$$$     *                       + DDOT(KM1,AP(IK+1),1,AP(IKP1+1),1)
c$$$                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)
c$$$                  IJ = 0
c$$$                  DO 200 J = 1, KM1
c$$$                     JK = IK + J
c$$$                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
c$$$                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
c$$$                     IJ = IJ + J
c$$$  200             CONTINUE
c$$$                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
c$$$  210          CONTINUE
c$$$               KSTEP = 2
c$$$  220       CONTINUE
c$$$C
c$$$C           SWAP
c$$$C
c$$$            KS = IABS(KPVT(K))
c$$$            IF (KS .EQ. K) GO TO 250
c$$$               IKS = (KS*(KS - 1))/2
c$$$               CALL DSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
c$$$               KSJ = IK + KS
c$$$               DO 230 JB = KS, K
c$$$                  J = K + KS - JB
c$$$                  JK = IK + J
c$$$                  TEMP = AP(JK)
c$$$                  AP(JK) = AP(KSJ)
c$$$                  AP(KSJ) = TEMP
c$$$                  KSJ = KSJ - (J - 1)
c$$$  230          CONTINUE
c$$$               IF (KSTEP .EQ. 1) GO TO 240
c$$$                  KSKP1 = IKP1 + KS
c$$$                  TEMP = AP(KSKP1)
c$$$                  AP(KSKP1) = AP(KKP1)
c$$$                  AP(KKP1) = TEMP
c$$$  240          CONTINUE
c$$$  250       CONTINUE
c$$$            IK = IK + K
c$$$            IF (KSTEP .EQ. 2) IK = IK + K + 1
c$$$            K = K + KSTEP
c$$$         GO TO 150
c$$$  260    CONTINUE
c$$$  270 CONTINUE
c$$$      RETURN
c$$$      END
c$$$#endif
