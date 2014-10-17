c CHOLESKY DECOMPOSITION
      subroutine dchdc(a,lda,p,work,jpvt,job,info)
      integer lda,p,jpvt(1),job,info
      double precision a(lda,p),work(p)
c
c     dchdc computes the cholesky decomposition of a positive definite
c     matrix.  a pivoting option allows the user to estimate the
c     condition of a positive definite matrix or determine the rank
c     of a positive semidefinite matrix.
c
c     on entry
c
c         a      double precision(lda,p).
c                a contains the matrix whose decomposition is to
c                be computed.  onlt the upper half of a need be stored.
c                the lower part of the array a is not referenced.
c
c         lda    integer.
c                lda is the leading dimension of the array a.
c
c         p      integer.
c                p is the order of the matrix.
c
c         work   double precision.
c                work is a work array.
c
c         jpvt   integer(p).
c                jpvt contains integers that control the selection
c                of the pivot elements, if pivoting has been requested.
c                each diagonal element a(k,k)
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      element.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free element.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final element.
c
c                before the decomposition is computed, initial elements
c                are moved by symmetric row and column interchanges to
c                the beginning of the array a and final
c                elements to the end.  both initial and final elements
c                are frozen in place during the computation and only
c                free elements are moved.  at the k-th stage of the
c                reduction, if a(k,k) is occupied by a free element
c                it is interchanged with the largest free element
c                a(l,l) with l .ge. k.  jpvt is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c         a      a contains in its upper half the cholesky factor
c                of the matrix a as it has been permuted by pivoting.
c
c         jpvt   jpvt(j) contains the index of the diagonal element
c                of a that was moved into the j-th position,
c                provided pivoting was requested.
c
c         info   contains the index of the last positive diagonal
c                element of the cholesky factor.
c
c     for positive definite matrices info = p is the normal return.
c     for pivoting with positive semidefinite matrices info will
c     in general be less than p.  however, info may be greater than
c     the rank of a, since rounding error can cause an otherwise zero
c     element to be positive. indefinite systems will always cause
c     info to be less than p.
c
c     linpack. this version dated 08/14/78 .
c     j.j. dongarra and g.w. stewart, argonne national laboratory and
c     university of maryland.
c
c
c     blas daxpy,dswap
c     fortran dsqrt
c
c     internal variables
c
      integer pu,pl,plp1,i,j,jp,jt,k,kb,km1,kp1,l,maxl
      double precision temp
      double precision maxdia
      logical swapk,negk
c
      pl = 1
      pu = 0
      info = p
      if (job .eq. 0) go to 160
c
c        pivoting has been requested. rearrange the
c        the elements according to jpvt.
c
         do 70 k = 1, p
            swapk = jpvt(k) .gt. 0
            negk = jpvt(k) .lt. 0
            jpvt(k) = k
            if (negk) jpvt(k) = -jpvt(k)
            if (.not.swapk) go to 60
               if (k .eq. pl) go to 50
                  call dswap(pl-1,a(1,k),1,a(1,pl),1)
                  temp = a(k,k)
                  a(k,k) = a(pl,pl)
                  a(pl,pl) = temp
                  plp1 = pl + 1
                  if (p .lt. plp1) go to 40
                  do 30 j = plp1, p
                     if (j .ge. k) go to 10
                        temp = a(pl,j)
                        a(pl,j) = a(j,k)
                        a(j,k) = temp
                     go to 20
   10                continue
                     if (j .eq. k) go to 20
                        temp = a(k,j)
                        a(k,j) = a(pl,j)
                        a(pl,j) = temp
   20                continue
   30             continue
   40             continue
                  jpvt(k) = jpvt(pl)
                  jpvt(pl) = k
   50          continue
               pl = pl + 1
   60       continue
   70    continue
         pu = p
         if (p .lt. pl) go to 150
         do 140 kb = pl, p
            k = p - kb + pl
            if (jpvt(k) .ge. 0) go to 130
               jpvt(k) = -jpvt(k)
               if (pu .eq. k) go to 120
                  call dswap(k-1,a(1,k),1,a(1,pu),1)
                  temp = a(k,k)
                  a(k,k) = a(pu,pu)
                  a(pu,pu) = temp
                  kp1 = k + 1
                  if (p .lt. kp1) go to 110
                  do 100 j = kp1, p
                     if (j .ge. pu) go to 80
                        temp = a(k,j)
                        a(k,j) = a(j,pu)
                        a(j,pu) = temp
                     go to 90
   80                continue
                     if (j .eq. pu) go to 90
                        temp = a(k,j)
                        a(k,j) = a(pu,j)
                        a(pu,j) = temp
   90                continue
  100             continue
  110             continue
                  jt = jpvt(k)
                  jpvt(k) = jpvt(pu)
                  jpvt(pu) = jt
  120          continue
               pu = pu - 1
  130       continue
  140    continue
  150    continue
  160 continue
      do 270 k = 1, p
c
c        reduction loop.
c
         maxdia = a(k,k)
         kp1 = k + 1
         maxl = k
c
c        determine the pivot element.
c
         if (k .lt. pl .or. k .ge. pu) go to 190
            do 180 l = kp1, pu
               if (a(l,l) .le. maxdia) go to 170
                  maxdia = a(l,l)
                  maxl = l
  170          continue
  180       continue
  190    continue
c
c        quit if the pivot element is not positive.
c
         if (maxdia .gt. 0.0d0) go to 200
            info = k - 1
c     ......exit
            go to 280
  200    continue
         if (k .eq. maxl) go to 210
c
c           start the pivoting and update jpvt.
c
            km1 = k - 1
            call dswap(km1,a(1,k),1,a(1,maxl),1)
            a(maxl,maxl) = a(k,k)
            a(k,k) = maxdia
            jp = jpvt(maxl)
            jpvt(maxl) = jpvt(k)
            jpvt(k) = jp
  210    continue
c
c        reduction step. pivoting is contained across the rows.
c
         work(k) = dsqrt(a(k,k))
         a(k,k) = work(k)
         if (p .lt. kp1) go to 260
         do 250 j = kp1, p
            if (k .eq. maxl) go to 240
               if (j .ge. maxl) go to 220
                  temp = a(k,j)
                  a(k,j) = a(j,maxl)
                  a(j,maxl) = temp
               go to 230
  220          continue
               if (j .eq. maxl) go to 230
                  temp = a(k,j)
                  a(k,j) = a(maxl,j)
                  a(maxl,j) = temp
  230          continue
  240       continue
            a(k,j) = a(k,j)/work(k)
            work(j) = a(k,j)
            temp = -a(k,j)
            call daxpy(j-k,temp,work(kp1),1,a(kp1,j),1)
  250    continue
  260    continue
  270 continue
  280 continue
      return
      end
c
c  S^1/2 - S^-1/2 DECOMPOSITION (from Jeppe)
      SUBROUTINE LS_SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)
*
* Calculate square root of positive definite symmetric matrix A
* if(ITASK .EQ. 2 ) Inverted square root matrix is also calculated
* In case of singularities in A A -1/2 is defined to have the same
* singularity
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
*
      DIMENSION A(NDIM,NDIM)
      DIMENSION ASQRT(NDIM,NDIM),AMSQRT(NDIM,NDIM)
      DIMENSION SCR(*)
* Length of SCR should at least be 2 * NDIM ** 2 + NDIM*(NDIM+1)/2
      KLFREE = 1
*
      KLASYM = KLFREE
      KLAVAL = KLASYM
      KLFREE = KLASYM + NDIM*(NDIM+1)/2
*
      KLAVEC = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      NTEST = 0
*
      CALL LS_TRIPAK(A,SCR(KLASYM),1,NDIM,NDIM)
      CALL LS_EIGEN(SCR(KLASYM),SCR(KLAVEC),NDIM,0,1)
      CALL LS_COPDIA(SCR(KLASYM),SCR(KLAVAL),NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        WRITE(6,*) ' Wrong use of LS_WRTMAT - Do not print : '
!        CALL( LS_WRTMATSCR(KLAVAL),NDIM,1,NDIM,1)
      END IF
*. Check for negative eigenvalues
      DO I = 1, NDIM
       IF(SCR(KLAVAL-1+I).LT.0.0D0) THEN
         WRITE(6,*) ' LS_SQRTMT : Negative eigenvalue ', SCR(KLAVAL-1+I)
         WRITE(6,*) ' LS_SQRTMT : I will STOP '
         STOP       ' LS_SQRTMT : Negative eigenvalue '
       END IF
      END DO
*
      DO 100 I = 1,NDIM
        SCR(KLAVAL-1+I) = SQRT(SCR(KLAVAL-1+I))
  100 CONTINUE
      CALL LS_XDIAXT(ASQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
*
      IF(ITASK .EQ. 2 ) THEN
        DO 200 I = 1,NDIM
          IF(SCR(KLAVAL-1+I) .GT. 1.0D-13 ) then
            SCR(KLAVAL-1+I) = 1.0D0/SCR(KLAVAL-1+I)
          ELSE
            SCR(KLAVAL-1+I) = SCR(KLAVAL-1+I)
          END IF
  200   CONTINUE
        CALL LS_XDIAXT(AMSQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
      END IF
*
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Info from LS_SQRTMT '
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' Input matrix to LS_SQRTMT '
        CALL LS_WRTMAT(A,NDIM,NDIM,NDIM,NDIM,1)
        WRITE(6,*) ' Square root of matrix '
        CALL LS_WRTMAT(ASQRT,NDIM,NDIM,NDIM,NDIM,1)
        IF(ITASK .EQ. 2 ) THEN
          WRITE(6,*) ' Inverse square root of matrix '
          CALL LS_WRTMAT(AMSQRT,NDIM,NDIM,NDIM,NDIM,1)
        END IF
      END IF
*
      RETURN
      END
C  /* Deck tripak */
      SUBROUTINE LS_TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
C
C ( NOT A SIMPLIFIED VERSION OF TETRAPAK )
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C   IWAY = 2 : PACKED TO FULL FORM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO 100 I = 1,NDIM
          DO 50  J = 1, I
           APAK(IJ+J) = AUTPAK(J,I)
   50     CONTINUE
          IJ = IJ + I
  100   CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO 200 I = 1,NDIM
          DO 150  J = 1, I
           AUTPAK(I,J) = APAK(IJ+J)
           AUTPAK(J,I) = APAK(IJ+J)
  150     CONTINUE
          IJ = IJ + I
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(*,*) ' AUTPAK AND APAK FROM LS_TRIPAK '
        CALL LS_WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM,0)
        CALL LS_PRSYM(APAK,NDIM)
      END IF
C
      RETURN
      END
C  /* Deck eigen */
      SUBROUTINE LS_EIGEN(A,R,N,MV,MFKR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N*N),R(N*N)
      DATA TESTIT/1.D-20/
      DATA TESTX/1.D-26/
      DATA TESTY/1.D-18/
C
C        PURPOSE
C           COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C           MATRIX
C
C        USAGE
C           CALL EIGEN(A,R,N,MV,MFKR)
C
C        DESCRIPTION OF PARAMETERS
C           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF
C               MATRIX A IN ASSCENDING ORDER.
C           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,
C               IN SAME SEQUENCE AS EIGENVALUES)
C           N - ORDER OF MATRICES A AND R
C           MV- INPUT CODE
C   0   COMPUTE EIGENVALUES AND EIGENVECTORS
C   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE
C       DIMENSIONED BUT MUST STILL APPEAR IN CALLING
C       SEQUENCE)
C           MFKR=0 NO SORT
C               =1 SORT
C
C        REMARKS
C           ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1)
C           MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED
C           BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN ?MATHEMATICAL
C           METHODS FOR DIGITAL COMPUTERS?, EDITED BY A. RALSTON AND
C           H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7
C
C     ..................................................................
C
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
C     DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
C    1 COSX2,SINCS,RANGE
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
C        40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT
C        62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD
C        BE CHANGED TO 1.0D-12.
C        900111-hjaaj: use generic SQRT and ABS
C
C        ...............................................................
C
C        GENERATE IDENTITY MATRIX
C
      RANGE=1.0D-12
      IF(MV-1) 10,25,10
   10 IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0D+00
      IF(I-J) 20,15,20
   15 R(IJ)=1.0D+00
   20 CONTINUE
C
C        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)
C
   25 ANORM=0.0D+00
      DO 35 I=1,N
      DO 35 J=I,N
      IF(I-J) 30,35,30
   30 IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF(ANORM) 165,165,40
   40 ANORM=1.414D+00*SQRT(ANORM)
      ANRMX=ANORM*RANGE/DFLOAT(N)
C
C        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
C
      IND=0
      THR=ANORM
   45 THR=THR/DFLOAT(N)
      IF(THR.LT.TESTY)THR=0.D0
   50 L=1
   55 M=L+1
C
C        COMPUTE SIN AND COS
C
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
      IF(ABS(A(LM)).LT.TESTY)A(LM)=0.D0
      IF(ABS(A(LM)).EQ.0.D0.AND.THR.EQ.0.D0)GO TO 130
      IF(ABS(A(LM))-THR) 130,65,65
   65 IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5D+00*(A(LL)-A(MM))
      AJUK=(A(LM)*A(LM)+X*X)
      AJUK=SQRT(AJUK)
      IF(ABS(AJUK).LT.TESTIT)WRITE(*,3000)TESTIT,AJUK,A(LM)
 3000 FORMAT(1H0,'***DENOMINATOR LT ',D14.6,'. VALUE=',D16.8,
     &'. NUMERATOR=',D16.8)
      Y=0.D0
      IF(ABS(AJUK).LT.TESTIT)GO TO 67
      Y=-A(LM)/AJUK
   67 CONTINUE
C  68 Y=-A(LM)/ SQRT(A(LM)*A(LM)+X*X)
      IF(X) 70,75,75
   70 Y=-Y
   75 AJUK=(1.D0-Y*Y)
      IF(AJUK.LT.0.D0)WRITE(*,3001) AJUK
 3001 FORMAT(1H0,'***SQRT OF ',D16.8)
      IF(AJUK.LT.0.D0)AJUK=0.D0
      AJUK=SQRT(AJUK)
      AJUK=2.D0*(1.D0+AJUK)
      AJUK=SQRT(AJUK)
      SINX=Y/AJUK
C     SINX=Y/ SQRT(2.0D+00*(1.0D+00+( SQRT(1.0D+00-Y*Y))))
      SINX2=SINX*SINX
C     COSX= SQRT(1.0D+00-SINX2)
      AJUK=1.D0-SINX2
      IF(AJUK.LT.TESTX)AJUK=0.D0
      COSX=SQRT(AJUK)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C
C        ROTATE L AND M COLUMNS
C
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      IF(I-L) 80,115,80
   80 IF(I-M) 85,115,90
   85 IM=I+MQ
      GO TO 95
   90 IM=M+IQ
   95 IF(I-L) 100,105,105
  100 IL=I+LQ
      GO TO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 IF(MV-1) 120,125,120
  120 ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0D+00*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C        TESTS FOR COMPLETION
C
C        TEST FOR M = LAST COLUMN
C
  130 IF(M-N) 135,140,135
  135 M=M+1
      GO TO 60
C
C        TEST FOR L = SECOND FROM LAST COLUMN
C
  140 IF(L-(N-1)) 145,150,145
  145 L=L+1
      GO TO 55
  150 IF(IND-1) 160,155,160
  155 IND=0
      GO TO 50
C
C        COMPARE THRESHOLD WITH FINAL NORM
C
  160 IF(THR-ANRMX) 165,165,45
C
C        SORT EIGENVALUES AND EIGENVECTORS
C
  165 IQ=-N
      IF(MFKR.EQ.0)GO TO 186
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF(A(MM)-A(LL)) 170,185,185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF(MV-1) 175,185,175
  175 DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
186   CONTINUE
      RETURN
      END
C  /* Deck wrtmat */
      SUBROUTINE LS_WRTMAT(AMATRX,NRDIM,NCDIM,NRMAX,NCMAX,ITRANS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AMATRX(NRMAX,NCMAX)
C
      IF ( NRDIM .EQ. 1 ) THEN
         WRITE(*,1011) (AMATRX(1,J),J=1,NCDIM)
      ELSE IF ( ITRANS.EQ. 0 ) THEN
         DO 100 I = 1,NRDIM
            WRITE(*,1010) I,(AMATRX(I,J),J=1,NCDIM)
  100    CONTINUE
      ELSE
         DO 101 I = 1, NCDIM
            WRITE(*,1010) I,(AMATRX(J,I),J=1,NRDIM)
  101    CONTINUE
      END IF
C
 1010 FORMAT(/,I6,1P,4E16.8,/,(6X,1P,4E16.8) )
 1011 FORMAT(/,(6X,1P,4E16.8) )
      RETURN
      END
C  /* Deck prsym */
      SUBROUTINE LS_PRSYM(A,MATDIM)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)
      JSTART=1
      JSTOP=0
      DO 100 I=1,MATDIM
        JSTART=JSTART+I-1
        JSTOP=JSTOP +I
        WRITE(*,1010) I,(A(J),J=JSTART,JSTOP)
  100 CONTINUE
      RETURN
 1010 FORMAT(1H0,2X,I3,5(1X,E14.7),/,(1H ,5X,5(1X,E14.7)))
      END

