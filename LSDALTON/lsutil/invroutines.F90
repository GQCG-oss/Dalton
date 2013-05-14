module inv_mod
! Routines for inversion                                                
! From Jeppe - late november 03                                         
contains
      SUBROUTINE INVERT_BY_DIAG2 (A, B, SCR, VEC, NDIM) 
!                                                                       
! Invert a symmetric  - hopefully nonsingular - matrix A                
! by diagonalization.                                                   
!                                                                       
! Jeppe Olsen, Oct 97 to check INVMAT                                   
!              March 00 : Scale initial matrix to obtain unit diagonal  
!                                                                       
!. The input matrix is delivered in A in complete form, and             
!. the output matrix is returned in A, also in complete form            
      IMPLICIT REAL (8)(A - H, O - Z) 
!. Input and output matrix                                              
      DIMENSION A (NDIM * NDIM) 
!. Scratch matrices and vector                                          
      DIMENSION B (NDIM * NDIM), SCR (NDIM * NDIM), VEC (NDIM) 
!                                                                       
      NTEST = 00 
!. Reform a to symmetric packed form                                    
      CALL LS_TRIPAK (A, SCR, 1, NDIM, NDIM) 
!. Extract diagonal                                                     
      CALL LS_COPDIA (SCR, VEC, NDIM, 1) 
!                                                                       
!.scale                                                                 
!                                                                       
      DO I = 1, NDIM 
!. Scaling vector                                                       
      IF (VEC (I) .EQ.0.0D0) THEN 
      VEC (I) = 1.0D0 
      ELSE 
      VEC (I) = 1.0D0 / SQRT (ABS (VEC (I) ) ) 
      ENDIF 
      ENDDO 
!. Scale matrix                                                         
      IJ = 0 
      DO I = 1, NDIM 
      DO J = 1, I 
      IJ = IJ + 1 
      SCR (IJ) = SCR (IJ) * VEC (I) * VEC (J) 
      ENDDO 
      ENDDO 
!     DO I = 1, NDIM                                                    
!       VEC(I) = 1.0D0/VEC(I)                                           
!     END DO                                                            
!. Diagonalize                                                          
      CALL LS_EIGEN (SCR, B, NDIM, 0, 1) 
!. Scale eigenvectors                                                   
      DO IVEC = 1, NDIM 
      IOFF = 1 + (IVEC - 1) * NDIM 
      CALL LS_VVTOV (B (IOFF), VEC, B (IOFF), NDIM) 
      ENDDO 
!.                                                                      
      CALL LS_COPDIA (SCR, VEC, NDIM, 1) 
      IF (NTEST.GE.1) THEN 
      WRITE (6, * ) ' Eigenvalues of scaled matrix : ' 
      CALL LS_WRTMAT (VEC, NDIM, 1, NDIM, 1,0) 
!      call lsquit('make proper code! wrong number of arguments',-1)
      ENDIF 
!. Invert diagonal elements                                             
      DO I = 1, NDIM 
      IF (ABS (VEC (I) ) .GT.1.0D-15) THEN 
      VEC (I) = 1.0D0 / VEC (I) 
      ELSE 
      VEC (I) = 0.0D0 
      WRITE (6, * ) ' Singular mode activated ' 
      ENDIF 
      ENDDO 
!. and obtain inverse matrix by transformation                          
!     XDIAXT(XDX,X,DIA,NDIM,SCR)                                        
      CALL XDIAXT (A, B, VEC, NDIM, SCR) 
!                                                                       
      IF (NTEST.GE.100) THEN 
      WRITE (6, * ) ' Inverse matrix from INVERSE_BY_DIAG' 
      CALL LS_WRTMAT (A, NDIM, NDIM, NDIM, NDIM,0) 
!      call lsquit('make proper code! wrong number of arguments',-1)
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE INVERT_BY_DIAG2                
      SUBROUTINE XDIAXT (XDX, X, DIA, NDIM, SCR) 
!                                                                       
! Obtain XDX = X * DIA * X(Transposed)                                  
! where DIA is an diagonal matrix stored as a vector                    
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION XDX (NDIM, NDIM) 
      DIMENSION X (NDIM, NDIM), DIA (NDIM) 
      DIMENSION SCR (NDIM, NDIM) 
!                                                                       
! DIA * X(transposed)                                                   
      DO 100 I = 1, NDIM 
      CALL LS_COPVEC (X (1, I), SCR (1, I), NDIM) 
      CALL LS_SCALVE (SCR (1, I), DIA (I), NDIM) 
  100 END DO 
! X * DIA * X(Transposed)                                               
      CALL LS_MATML4 (XDX, X, SCR, NDIM, NDIM, NDIM, NDIM, NDIM, NDIM,  &
      2)                                                                
!                                                                       
      RETURN 
      END SUBROUTINE XDIAXT                         
                                                                        
      SUBROUTINE LS_COPDIA (A, VEC, NDIM, IPACK) 
!                                                                       
! Copy diagonal of matrix A into vector VEC                             
!                                                                       
!   IPACK = 0 : Full matrix                                             
!   IPACK = 1 : Lower triangular matrix                                 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ), VEC ( * ) 
!                                                                       
      IF (IPACK.EQ.0) THEN 
      DO 100 I = 1, NDIM 
      VEC (I) = A ( (I - 1) * NDIM + I) 
  100 END DO 
      ELSE 
      DO 200 I = 1, NDIM 
      VEC (I) = A (I * (I + 1) / 2) 
  200 END DO 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE LS_COPDIA                      
!  /* Deck vvtov */                                                     
      SUBROUTINE LS_VVTOV (VECIN1, VECIN2, VECUT, NDIM) 
!                                                                       
! VECUT(I) = VECIN1(I) * VECIN2(I)                                      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION VECIN1 (NDIM), VECIN2 (NDIM), VECUT (NDIM) 
!                                                                       
      DO 100 I = 1, NDIM 
      VECUT (I) = VECIN1 (I) * VECIN2 (I) 
  100 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE LS_VVTOV                       
!  /* Deck scalve */                                                    
      SUBROUTINE LS_SCALVE (VECTOR, FACTOR, NDIM) 
!                                                                       
! CALCULATE SCALAR(FACTOR) TIMES VECTOR                                 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION VECTOR ( * ) 
!                                                                       
      DO 100 I = 1, NDIM 
      VECTOR (I) = VECTOR (I) * FACTOR 
  100 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE LS_SCALVE                      
!  /* Deck matml4 */                                                    
      SUBROUTINE LS_MATML4 (C, A, B, NCROW, NCCOL, NAROW, NACOL, NBROW, &
      NBCOL, ITRNSP)                                                    
!                                                                       
! MULTIPLY A AND B TO GIVE C                                            
!                                                                       
!     C = A * B             FOR ITRNSP = 0                              
!                                                                       
!     C = A(TRANSPOSED) * B FOR ITRNSP = 1                              
!                                                                       
!     C = A * B(TRANSPOSED) FOR ITRNSP = 2                              
!                                                                       
!... JEPPE OLSEN, LAST REVISION JULY 24 1987                            
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (NAROW, NACOL), B (NBROW, NBCOL) 
      DIMENSION C (NCROW, NCCOL) 
!                                                                       
      NTEST = 0 
      IF (NTEST.NE.0) THEN 
      WRITE ( *, * ) 
      WRITE ( * , * ) ' A AND B MATRIX FROM LS_MATML4, ITRNSP =',       &
      ITRNSP                                                            
      WRITE ( *, * ) 
      CALL LS_WRTMAT (A, NAROW, NACOL, NAROW, NACOL, 0) 
      CALL LS_WRTMAT (B, NBROW, NBCOL, NBROW, NBCOL, 0) 
      ENDIF 
!                                                                       
      CALL LS_SETVEC (C, 0.0D0, NCROW * NCCOL) 
!                                                                       
      IF (ITRNSP.NE.0) GOTO 001 
      DO 50 J = 1, NCCOL 
      DO 40 K = 1, NBROW 
      BKJ = B (K, J) 
      DO 30 I = 1, NCROW 
      C (I, J) = C (I, J) + A (I, K) * BKJ 
   30 END DO 
   40 END DO 
   50 END DO 
!                                                                       
!                                                                       
    1 CONTINUE 
!                                                                       
      IF (ITRNSP.NE.1) GOTO 101 
!... C = A(T) * B                                                       
      DO 150 J = 1, NCCOL 
      DO 140 K = 1, NBROW 
      BKJ = B (K, J) 
      DO 130 I = 1, NCROW 
      C (I, J) = C (I, J) + A (K, I) * BKJ 
  130 END DO 
  140 END DO 
  150 END DO 
!                                                                       
  101 CONTINUE 
!                                                                       
      IF (ITRNSP.NE.2) GOTO 201 
!... C = A*B(T)                                                         
      DO 250 J = 1, NCCOL 
      DO 240 K = 1, NBCOL 
      BJK = B (J, K) 
      DO 230 I = 1, NCROW 
      C (I, J) = C (I, J) + A (I, K) * BJK 
  230 END DO 
  240 END DO 
  250 END DO 
!                                                                       
!                                                                       
  201 CONTINUE 
!                                                                       
      IF (NTEST.NE.0) THEN 
      WRITE ( *, * ) 
      WRITE ( * , * ) ' C MATRIX FROM LS_MATML4, ITRNSP =', ITRNSP 
      WRITE ( *, * ) 
      CALL LS_WRTMAT (C, NCROW, NCCOL, NCROW, NCCOL, 0) 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE LS_MATML4                      
!  /* Deck copvec */                                                    
      SUBROUTINE LS_COPVEC (VECIN, VECOUT, NDIM) 
!                                                                       
! 880717 - HJAaJ - written based on a qualified guess                   
!                  about Jeppe's original                               
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION VECIN (NDIM), VECOUT (NDIM) 
      DO 100 I = 1, NDIM 
      VECOUT (I) = VECIN (I) 
  100 END DO 
      RETURN 
      END SUBROUTINE LS_COPVEC                      
!  /* Deck setvec */                                                    
      SUBROUTINE LS_SETVEC (VECTOR, VALUE, NDIM) 
!                                                                       
! VECTOR(*) = VALUE                                                     
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION VECTOR (NDIM) 
!                                                                       
      DO 100 I = 1, NDIM 
      VECTOR (I) = VALUE 
  100 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE LS_SETVEC                      

!  /* Deck wrtmat */                                                    
      SUBROUTINE LS_WRTMAT (AMATRX, NRDIM, NCDIM, NRMAX, NCMAX, ITRANS) 
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION AMATRX (NRMAX, NCMAX) 
!                                                                       
      IF (NRDIM.EQ.1) THEN 
      WRITE ( *, 1011) (AMATRX (1, J), J = 1, NCDIM) 
      ELSEIF (ITRANS.EQ.0) THEN 
      DO 100 I = 1, NRDIM 
      WRITE ( *, 1010) I, (AMATRX (I, J), J = 1, NCDIM) 
  100 END DO 
      ELSE 
      DO 101 I = 1, NCDIM 
      WRITE ( *, 1010) I, (AMATRX (J, I), J = 1, NRDIM) 
  101 END DO 
      ENDIF 
!                                                                       
 1010 FORMAT(/,I6,1P,4E16.8,/,(6X,1P,4E16.8) ) 
 1011 FORMAT(/,(6X,1P,4E16.8) ) 
      RETURN 
      END SUBROUTINE LS_WRTMAT                      
!  /* Deck prsym */                                                     
      SUBROUTINE LS_PRSYM (A, MATDIM) 
! PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.           
! THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.           
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A ( * ) 
      JSTART = 1 
      JSTOP = 0 
      DO 100 I = 1, MATDIM 
      JSTART = JSTART + I - 1 
      JSTOP = JSTOP + I 
      WRITE ( *, 1010) I, (A (J), J = JSTART, JSTOP) 
  100 END DO 
      RETURN 
 1010 FORMAT(1H0,2X,I3,5(1X,E13.7),/,(1H ,5X,5(1X,E13.7))) 
      END SUBROUTINE LS_PRSYM                       
   
!  /* Deck tripak */                                                    
      SUBROUTINE LS_TRIPAK (AUTPAK, APAK, IWAY, MATDIM, NDIM) 
!                                                                       
! ( NOT A SIMPLIFIED VERSION OF TETRAPAK )                              
!                                                                       
!.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING                        
!   AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX                         
!                                                                       
!   IWAY = 1 : FULL TO PACKED                                           
!   IWAY = 2 : PACKED TO FULL FORM                                      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION AUTPAK (MATDIM, MATDIM), APAK ( * ) 
!                                                                       
      IF (IWAY.EQ.1) THEN 
      IJ = 0 
      DO 100 I = 1, NDIM 
      DO 50 J = 1, I 
      APAK (IJ + J) = AUTPAK (J, I) 
   50 END DO 
      IJ = IJ + I 
  100 END DO 
      ENDIF 
!                                                                       
      IF (IWAY.EQ.2) THEN 
      IJ = 0 
      DO 200 I = 1, NDIM 
      DO 150 J = 1, I 
      AUTPAK (I, J) = APAK (IJ + J) 
      AUTPAK (J, I) = APAK (IJ + J) 
  150 END DO 
      IJ = IJ + I 
  200 END DO 
      ENDIF 
!                                                                       
      NTEST = 0 
      IF (NTEST.NE.0) THEN 
      WRITE ( * , * ) ' AUTPAK AND APAK FROM LS_TRIPAK ' 
      CALL LS_WRTMAT (AUTPAK, NDIM, MATDIM, NDIM, MATDIM, 0) 
      CALL LS_PRSYM (APAK, NDIM) 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE LS_TRIPAK                      
!  /* Deck eigen */                                                     
      SUBROUTINE LS_EIGEN (A, R, N, MV, MFKR) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (N * N), R (N * N) 
      DATA TESTIT / 1.D-20 / 
      DATA TESTX / 1.D-26 / 
      DATA TESTY / 1.D-18 / 
!                                                                       
!        PURPOSE                                                        
!           COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC    
!           MATRIX                                                      
!                                                                       
!        USAGE                                                          
!           CALL EIGEN(A,R,N,MV,MFKR)                                   
!                                                                       
!        DESCRIPTION OF PARAMETERS                                      
!           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.  
!               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF      
!               MATRIX A IN ASSCENDING ORDER.                           
!           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,    
!               IN SAME SEQUENCE AS EIGENVALUES)                        
!           N - ORDER OF MATRICES A AND R                               
!           MV- INPUT CODE                                              
!   0   COMPUTE EIGENVALUES AND EIGENVECTORS                            
!   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE                         
!       DIMENSIONED BUT MUST STILL APPEAR IN CALLING                    
!       SEQUENCE)                                                       
!           MFKR=0 NO SORT                                              
!               =1 SORT                                                 
!                                                                       
!        REMARKS                                                        
!           ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1)   
!           MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R         
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED     
!           BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN ?MATHEMATICAL
!           METHODS FOR DIGITAL COMPUTERS?, EDITED BY A. RALSTON AND    
!           H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7   
!                                                                       
!     ..................................................................
!                                                                       
!                                                                       
!        ...............................................................
!                                                                       
!        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  
!        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      
!        STATEMENT WHICH FOLLOWS.                                       
!                                                                       
!     DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,         
!    1 COSX2,SINCS,RANGE                                                
!                                                                       
!        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    
!        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      
!        ROUTINE.                                                       
!                                                                       
!        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      
!        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
!        40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT 
!        62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD 
!        BE CHANGED TO 1.0D-12.                                         
!        900111-hjaaj: use generic SQRT and ABS                         
!                                                                       
!        ...............................................................
!                                                                       
!        GENERATE IDENTITY MATRIX                                       
!                                                                       
      RANGE = 1.0D-12 
      IF (MV - 1) 10, 25, 10 
   10 IQ = - N 
      DO 20 J = 1, N 
      IQ = IQ + N 
      DO 20 I = 1, N 
      IJ = IQ + I 
      R (IJ) = 0.0D+00 
      IF (I - J) 20, 15, 20 
   15 R (IJ) = 1.0D+00 
   20 CONTINUE 
!                                                                       
!        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)             
!                                                                       
   25 ANORM = 0.0D+00 
      DO 35 I = 1, N 
      DO 35 J = I, N 
      IF (I - J) 30, 35, 30 
   30 IA = I + (J * J - J) / 2 
      ANORM = ANORM + A (IA) * A (IA) 
   35 CONTINUE 
      IF (ANORM) 165, 165, 40 
   40 ANORM = 1.414D+00 * SQRT (ANORM) 
      ANRMX = ANORM * RANGE / DFLOAT (N) 
!                                                                       
!        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR               
!                                                                       
      IND = 0 
      THR = ANORM 
   45 THR = THR / DFLOAT (N) 
      IF (THR.LT.TESTY) THR = 0.D0 
   50 L = 1 
   55 M = L + 1 
!                                                                       
!        COMPUTE SIN AND COS                                            
!                                                                       
   60 MQ = (M * M - M) / 2 
      LQ = (L * L - L) / 2 
      LM = L + MQ 
      IF (ABS (A (LM) ) .LT.TESTY) A (LM) = 0.D0 
      IF (ABS (A (LM) ) .EQ.0.D0.AND.THR.EQ.0.D0) GOTO 130 
      IF (ABS (A (LM) ) - THR) 130, 65, 65 
   65 IND = 1 
      LL = L + LQ 
      MM = M + MQ 
      X = 0.5D+00 * (A (LL) - A (MM) ) 
      AJUK = (A (LM) * A (LM) + X * X) 
      AJUK = SQRT (AJUK) 
      IF (ABS (AJUK) .LT.TESTIT) WRITE ( *, 3000) TESTIT, AJUK, A (LM) 
 3000 FORMAT(1H0,'***DENOMINATOR LT ',D12.6,'. VALUE=',D14.8,           &
     &'. NUMERATOR=',D14.8)                                             
      Y = 0.D0 
      IF (ABS (AJUK) .LT.TESTIT) GOTO 67 
      Y = - A (LM) / AJUK 
   67 CONTINUE 
!  68 Y=-A(LM)/ SQRT(A(LM)*A(LM)+X*X)                                   
      IF (X) 70, 75, 75 
   70 Y = - Y 
   75 AJUK = (1.D0 - Y * Y) 
      IF (AJUK.LT.0.D0) WRITE ( *, 3001) AJUK 
 3001 FORMAT(1H0,'***SQRT OF ',D14.8) 
      IF (AJUK.LT.0.D0) AJUK = 0.D0 
      AJUK = SQRT (AJUK) 
      AJUK = 2.D0 * (1.D0 + AJUK) 
      AJUK = SQRT (AJUK) 
      SINX = Y / AJUK 
!     SINX=Y/ SQRT(2.0D+00*(1.0D+00+( SQRT(1.0D+00-Y*Y))))              
      SINX2 = SINX * SINX 
!     COSX= SQRT(1.0D+00-SINX2)                                         
      AJUK = 1.D0 - SINX2 
      IF (AJUK.LT.TESTX) AJUK = 0.D0 
      COSX = SQRT (AJUK) 
      COSX2 = COSX * COSX 
      SINCS = SINX * COSX 
!                                                                       
!        ROTATE L AND M COLUMNS                                         
!                                                                       
      ILQ = N * (L - 1) 
      IMQ = N * (M - 1) 
      DO 125 I = 1, N 
      IQ = (I * I - I) / 2 
      IF (I - L) 80, 115, 80 
   80 IF (I - M) 85, 115, 90 
   85 IM = I + MQ 
      GOTO 95 
   90 IM = M + IQ 
   95 IF (I - L) 100, 105, 105 
  100 IL = I + LQ 
      GOTO 110 
  105 IL = L + IQ 
  110 X = A (IL) * COSX - A (IM) * SINX 
      A (IM) = A (IL) * SINX + A (IM) * COSX 
      A (IL) = X 
  115 IF (MV - 1) 120, 125, 120 
  120 ILR = ILQ + I 
      IMR = IMQ + I 
      X = R (ILR) * COSX - R (IMR) * SINX 
      R (IMR) = R (ILR) * SINX + R (IMR) * COSX 
      R (ILR) = X 
  125 END DO 
      X = 2.0D+00 * A (LM) * SINCS 
      Y = A (LL) * COSX2 + A (MM) * SINX2 - X 
      X = A (LL) * SINX2 + A (MM) * COSX2 + X 
      A (LM) = (A (LL) - A (MM) ) * SINCS + A (LM) * (COSX2 - SINX2) 
      A (LL) = Y 
      A (MM) = X 
!                                                                       
!        TESTS FOR COMPLETION                                           
!                                                                       
!        TEST FOR M = LAST COLUMN                                       
!                                                                       
  130 IF (M - N) 135, 140, 135 
  135 M = M + 1 
      GOTO 60 
!                                                                       
!        TEST FOR L = SECOND FROM LAST COLUMN                           
!                                                                       
  140 IF (L - (N - 1) ) 145, 150, 145 
  145 L = L + 1 
      GOTO 55 
  150 IF (IND-1) 160, 155, 160 
  155 IND = 0 
      GOTO 50 
!                                                                       
!        COMPARE THRESHOLD WITH FINAL NORM                              
!                                                                       
  160 IF (THR - ANRMX) 165, 165, 45 
!                                                                       
!        SORT EIGENVALUES AND EIGENVECTORS                              
!                                                                       
  165 IQ = - N 
      IF (MFKR.EQ.0) GOTO 186 
      DO 185 I = 1, N 
      IQ = IQ + N 
      LL = I + (I * I - I) / 2 
      JQ = N * (I - 2) 
      DO 185 J = I, N 
      JQ = JQ + N 
      MM = J + (J * J - J) / 2 
      IF (A (MM) - A (LL) ) 170, 185, 185 
  170 X = A (LL) 
      A (LL) = A (MM) 
      A (MM) = X 
      IF (MV - 1) 175, 185, 175 
  175 DO 180 K = 1, N 
      ILR = IQ + K 
      IMR = JQ + K 
      X = R (ILR) 
      R (ILR) = R (IMR) 
  180 R (IMR) = X 
  185 CONTINUE 
  186 CONTINUE 
      RETURN 
      END SUBROUTINE LS_EIGEN                       
   
end module
