module chol_decomp_mod
  use inv_mod
! CHOLESKY DECOMPOSITION                                                
  contains
      SUBROUTINE dchdc (a, lda, p, work, jpvt, job, info) 
      INTEGER lda, p, jpvt (p), job, info 
      DOUBLEPRECISION a (lda, p), work (p) 
!                                                                       
!     dchdc computes the cholesky decomposition of a positive definite  
!     matrix.  a pivoting option allows the user to estimate the        
!     condition of a positive definite matrix or determine the rank     
!     of a positive semidefinite matrix.                                
!                                                                       
!     on entry                                                          
!                                                                       
!         a      double precision(lda,p).                               
!                a contains the matrix whose decomposition is to        
!                be computed.  onlt the upper half of a need be stored. 
!                the lower part of the array a is not referenced.       
!                                                                       
!         lda    integer.                                               
!                lda is the leading dimension of the array a.           
!                                                                       
!         p      integer.                                               
!                p is the order of the matrix.                          
!                                                                       
!         work   double precision.                                      
!                work is a work array.                                  
!                                                                       
!         jpvt   integer(p).                                            
!                jpvt contains integers that control the selection      
!                of the pivot elements, if pivoting has been requested. 
!                each diagonal element a(k,k)                           
!                is placed in one of three classes according to the     
!                value of jpvt(k).                                      
!                                                                       
!                   if jpvt(k) .gt. 0, then x(k) is an initial          
!                                      element.                         
!                                                                       
!                   if jpvt(k) .eq. 0, then x(k) is a free element.     
!                                                                       
!                   if jpvt(k) .lt. 0, then x(k) is a final element.    
!                                                                       
!                before the decomposition is computed, initial elements 
!                are moved by symmetric row and column interchanges to  
!                the beginning of the array a and final                 
!                elements to the end.  both initial and final elements  
!                are frozen in place during the computation and only    
!                free elements are moved.  at the k-th stage of the     
!                reduction, if a(k,k) is occupied by a free element     
!                it is interchanged with the largest free element       
!                a(l,l) with l .ge. k.  jpvt is not referenced if       
!                job .eq. 0.                                            
!                                                                       
!        job     integer.                                               
!                job is an integer that initiates column pivoting.      
!                if job .eq. 0, no pivoting is done.                    
!                if job .ne. 0, pivoting is done.                       
!                                                                       
!     on return                                                         
!                                                                       
!         a      a contains in its upper half the cholesky factor       
!                of the matrix a as it has been permuted by pivoting.   
!                                                                       
!         jpvt   jpvt(j) contains the index of the diagonal element     
!                of a that was moved into the j-th position,            
!                provided pivoting was requested.                       
!                                                                       
!         info   contains the index of the last positive diagonal       
!                element of the cholesky factor.                        
!                                                                       
!     for positive definite matrices info = p is the normal return.     
!     for pivoting with positive semidefinite matrices info will        
!     in general be less than p.  however, info may be greater than     
!     the rank of a, since rounding error can cause an otherwise zero   
!     element to be positive. indefinite systems will always cause      
!     info to be less than p.                                           
!                                                                       
!     linpack. this version dated 08/14/78 .                            
!     j.j. dongarra and g.w. stewart, argonne national laboratory and   
!     university of maryland.                                           
!                                                                       
!                                                                       
!     blas daxpy,dswap                                                  
!     fortran dsqrt                                                     
!                                                                       
!     internal variables                                                
!                                                                       
      INTEGER pu, pl, plp1, i, j, jp, jt, k, kb, km1, kp1, l, maxl 
      DOUBLEPRECISION temp 
      DOUBLEPRECISION maxdia 
      LOGICAL swapk, negk 
!                                                                       
      pl = 1 
      pu = 0 
      info = p 
      IF (job.eq.0) goto 160 
!                                                                       
!        pivoting has been requested. rearrange the                     
!        the elements according to jpvt.                                
!                                                                       
      DO 70 k = 1, p 
      swapk = jpvt (k) .gt.0 
      negk = jpvt (k) .lt.0 
      jpvt (k) = k 
      IF (negk) jpvt (k) = - jpvt (k) 
      IF (.not.swapk) goto 60 
      IF (k.eq.pl) goto 50 
      CALL dswap (pl - 1, a (1, k), 1, a (1, pl), 1) 
      temp = a (k, k) 
      a (k, k) = a (pl, pl) 
      a (pl, pl) = temp 
      plp1 = pl + 1 
      IF (p.lt.plp1) goto 40 
      DO 30 j = plp1, p 
      IF (j.ge.k) goto 10 
      temp = a (pl, j) 
      a (pl, j) = a (j, k) 
      a (j, k) = temp 
      GOTO 20 
   10 CONTINUE 
      IF (j.eq.k) goto 20 
      temp = a (k, j) 
      a (k, j) = a (pl, j) 
      a (pl, j) = temp 
   20 CONTINUE 
   30 END DO 
   40 CONTINUE 
      jpvt (k) = jpvt (pl) 
      jpvt (pl) = k 
   50 CONTINUE 
      pl = pl + 1 
   60 CONTINUE 
   70 END DO 
      pu = p 
      IF (p.lt.pl) goto 150 
      DO 140 kb = pl, p 
      k = p - kb + pl 
      IF (jpvt (k) .ge.0) goto 130 
      jpvt (k) = - jpvt (k) 
      IF (pu.eq.k) goto 120 
      CALL dswap (k - 1, a (1, k), 1, a (1, pu), 1) 
      temp = a (k, k) 
      a (k, k) = a (pu, pu) 
      a (pu, pu) = temp 
      kp1 = k + 1 
      IF (p.lt.kp1) goto 110 
      DO 100 j = kp1, p 
      IF (j.ge.pu) goto 80 
      temp = a (k, j) 
      a (k, j) = a (j, pu) 
      a (j, pu) = temp 
      GOTO 90 
   80 CONTINUE 
      IF (j.eq.pu) goto 90 
      temp = a (k, j) 
      a (k, j) = a (pu, j) 
      a (pu, j) = temp 
   90 CONTINUE 
  100 END DO 
  110 CONTINUE 
      jt = jpvt (k) 
      jpvt (k) = jpvt (pu) 
      jpvt (pu) = jt 
  120 CONTINUE 
      pu = pu - 1 
  130 CONTINUE 
  140 END DO 
  150 CONTINUE 
  160 CONTINUE 
      DO 270 k = 1, p 
!                                                                       
!        reduction loop.                                                
!                                                                       
      maxdia = a (k, k) 
      kp1 = k + 1 
      maxl = k 
!                                                                       
!        determine the pivot element.                                   
!                                                                       
      IF (k.lt.pl.or.k.ge.pu) goto 190 
      DO 180 l = kp1, pu 
      IF (a (l, l) .le.maxdia) goto 170 
      maxdia = a (l, l) 
      maxl = l 
  170 CONTINUE 
  180 END DO 
  190 CONTINUE 
!                                                                       
!        quit if the pivot element is not positive.                     
!                                                                       
      IF (maxdia.gt.0.0d0) goto 200 
      info = k - 1 
!     ......exit                                                        
      GOTO 280 
  200 CONTINUE 
      IF (k.eq.maxl) goto 210 
!                                                                       
!           start the pivoting and update jpvt.                         
!                                                                       
      km1 = k - 1 
      CALL dswap (km1, a (1, k), 1, a (1, maxl), 1) 
      a (maxl, maxl) = a (k, k) 
      a (k, k) = maxdia 
      jp = jpvt (maxl) 
      jpvt (maxl) = jpvt (k) 
      jpvt (k) = jp 
  210 CONTINUE 
!                                                                       
!        reduction step. pivoting is contained across the rows.         
!                                                                       
      work (k) = dsqrt (a (k, k) ) 
      a (k, k) = work (k) 
      IF (p.lt.kp1) goto 260 
      DO 250 j = kp1, p 
      IF (k.eq.maxl) goto 240 
      IF (j.ge.maxl) goto 220 
      temp = a (k, j) 
      a (k, j) = a (j, maxl) 
      a (j, maxl) = temp 
      GOTO 230 
  220 CONTINUE 
      IF (j.eq.maxl) goto 230 
      temp = a (k, j) 
      a (k, j) = a (maxl, j) 
      a (maxl, j) = temp 
  230 CONTINUE 
  240 CONTINUE 
      a (k, j) = a (k, j) / work (k) 
      work (j) = a (k, j) 
      temp = - a (k, j) 
      CALL daxpy (j - k, temp, work (kp1), 1, a (kp1, j), 1) 
  250 END DO 
  260 CONTINUE 
  270 END DO 
  280 CONTINUE 
      RETURN 
      END SUBROUTINE dchdc                          
!                                                                       
!  S^1/2 - S^-1/2 DECOMPOSITION (from Jeppe)                            
      SUBROUTINE SQRTMT (A, NDIM, ITASK, ASQRT, AMSQRT, SCR) 
!                                                                       
! Calculate square root of positive definite symmetric matrix A         
! if(ITASK .EQ. 2 ) Inverted square root matrix is also calculated      
! In case of singularities in A A -1/2 is defined to have the same      
! singularity                                                           
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DIMENSION A (NDIM, NDIM) 
      DIMENSION ASQRT (NDIM, NDIM), AMSQRT (NDIM, NDIM) 
      DIMENSION SCR ( * ) 
! Length of SCR should at least be 2 * NDIM ** 2 + NDIM*(NDIM+1)/2      
      KLFREE = 1 
!                                                                       
      KLASYM = KLFREE 
      KLAVAL = KLASYM 
      KLFREE = KLASYM + NDIM * (NDIM + 1) / 2 
!                                                                       
      KLAVEC = KLFREE 
      KLFREE = KLFREE+NDIM**2 
!                                                                       
      NTEST = 0 
!                                                                       
      CALL LS_TRIPAK (A, SCR (KLASYM), 1, NDIM, NDIM) 
      CALL LS_EIGEN (SCR (KLASYM), SCR (KLAVEC), NDIM, 0, 1) 
      CALL LS_COPDIA (SCR (KLASYM), SCR (KLAVAL), NDIM, 1) 
      IF (NTEST.GE.1) THEN 
      WRITE (6, * ) ' Eigenvalues of matrix : ' 
      WRITE (6, * ) ' Wrong use of LS_WRTMAT - Do not print : ' 
!        CALL( LS_WRTMATSCR(KLAVAL),NDIM,1,NDIM,1)                      
      ENDIF 
!. Check for negative eigenvalues                                       
      DO I = 1, NDIM 
      IF (SCR (KLAVAL - 1 + I) .LT.0.0D0) THEN 
      WRITE (6, * ) ' SQRTMT : Negative eigenvalue ', SCR (KLAVAL - 1 + &
      I)                                                                
      WRITE (6, * ) ' SQRTMT : I will STOP ' 
      STOP ' SQRTMT : Negative eigenvalue ' 
      ENDIF 
      ENDDO 
!                                                                       
      DO 100 I = 1, NDIM 
      SCR (KLAVAL - 1 + I) = SQRT (SCR (KLAVAL - 1 + I) ) 
  100 END DO 
      CALL XDIAXT (ASQRT, SCR (KLAVEC), SCR (KLAVAL), NDIM, SCR (KLFREE)&
      )                                                                 
!                                                                       
      IF (ITASK.EQ.2) THEN 
      DO 200 I = 1, NDIM 
      IF (SCR (KLAVAL - 1 + I) .GT.1.0D-13) then 
      SCR (KLAVAL - 1 + I) = 1.0D0 / SCR (KLAVAL - 1 + I) 
      ELSE 
      SCR (KLAVAL - 1 + I) = SCR (KLAVAL - 1 + I) 
      ENDIF 
  200 END DO 
      CALL XDIAXT (AMSQRT, SCR (KLAVEC), SCR (KLAVAL), NDIM, SCR (      &
      KLFREE) )                                                         
      ENDIF 
!                                                                       
      IF (NTEST.GE.1) THEN 
      WRITE (6, * ) ' Info from SQRTMT ' 
      WRITE (6, * ) ' =================' 
      WRITE (6, * ) ' Input matrix to SQRTMT ' 
      CALL LS_WRTMAT (A, NDIM, NDIM, NDIM, NDIM, 1) 
      WRITE (6, * ) ' Square root of matrix ' 
      CALL LS_WRTMAT (ASQRT, NDIM, NDIM, NDIM, NDIM, 1) 
      IF (ITASK.EQ.2) THEN 
      WRITE (6, * ) ' Inverse square root of matrix ' 
      CALL LS_WRTMAT (AMSQRT, NDIM, NDIM, NDIM, NDIM, 1) 
      ENDIF 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE SQRTMT                         
   end module
