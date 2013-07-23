c     @file blas_wrapper.F This file wraps around blas routines allowing
c     to call native, hardware optimized routines and still use
c     nonnative 64-bit integers throughout the program.  The entire
c     remaining code of dalton uses -Defines to transform relevant blas
c     calls into names prefixed with 6. This code carefully undefs so
c     that proper calls to native BLAS can be generated.  For time
c     being, we wrap around only BLAS-3 routines, they are the ones
c     where there gain is most noticeable.
#undef DROT
      SUBROUTINE  DROT6(n,dx,incx,dy,incy,c,s)
      double precision dx(*),dy(*),c,s
      INTEGER N, INCX,INCY
      INTEGER*4 N4, INCX4, INCY4
      N4 = N
      INCX4 = INCX
      INCY4 = INCY
      CALL DROT(N4,DX,INCX4,DY,INCY4,C,S)
      END
#undef DASUM
      DOUBLE PRECISION FUNCTION DASUM6(N,DX,INCX)
      DOUBLE PRECISION DX(1),DASUM
      INTEGER*4 N4,INCX
      N4 = N
      INCX4 = INCX
      DASUM6 = DASUM(N4,DX,INCX4)
      END
#undef DDOT
      FUNCTION DDOT6(N,DX,INCX,DY,INCY)
      DOUBLE PRECISION DDOT, DDOT6,DX(1),DY(1)
      INTEGER INCX,INCY,N
      INTEGER*4 N4, INCX4, INCY4
      N4 = N
      INCX4 = INCX
      INCY4 = INCY
      DDOT6 = DDOT(N4,DX,INCX4,DY,INCY4)
      END
#undef DCOPY
      SUBROUTINE DCOPY6(N,DX,INCX,DY,INCY)
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER*4 N4,INCX4,INCY4
      N4 = N
      INCX4 = INCX
      INCY4 = INCY
      CALL DCOPY(N4,DX,INCX4,DY,INCY4)
      END
#undef DSWAP
      SUBROUTINE  DSWAP6(N,DX,INCX,DY,INCY)
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER INCX,INCY,N
      INTEGER*4 N4, INCX4, INCY4
      N4 = N
      INCX4 = INCX
      INCY4 = INCY
      CALL DSWAP(N4,DX,INCX4,DY,INCY4)
      END
#undef DSCAL
      SUBROUTINE DSCAL6(N,DA,DX,INCX)
      DOUBLE PRECISION DA,DX(*)
      INTEGER INCX,N
      INTEGER N4, INCX4
      N4 = N
      INCX4 = INCX
      CALL DSCAL(N4,DA,DX,INCX4)
      END
#undef DAXPY
      SUBROUTINE DAXPY6(N,DA,DX,INCX,DY,INCY)
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER INCX,INCY,N
      INTEGER*4 INCX4, INCY4,N4
      N4 = N
      INCX4 = INCX
      INCY4 = INCY
      CALL DAXPY(N4,DA,DX,INCX4,DY,INCY4)
      END
#undef DSYMV
      SUBROUTINE DSYMV6( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
      INTEGER*4 N4, LDA4, INCX4, INCY4
      N4 = N
      LDA4 = LDA
      INCX4 = INCX
      INCY4 = INCY
      CALL DSYMV(UPLO,N4,ALPHA,A,LDA4,X,INCX4,BETA,Y,INCY4)
      END
#undef DGEMV
      SUBROUTINE DGEMV6( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
      INTEGER*4 M4, N4, LDA4, INCX4,INCY4
      M4 = M
      N4 = N
      LDA4 = LDA
      INCX4 = INCX
      INCY4 = INCY
      CALL DGEMV(TRANS,M4,N4,ALPHA,A,LDA4,X,INCX4,BETA,Y,INCY4)
      END
#undef DTRMV
      SUBROUTINE DTRMV6( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
      DOUBLE PRECISION   A( LDA, * ), X( * )
      INTEGER*4 N4, INCX4, LDA4
      N4 = N
      LDA4 = LDA
      INCX4 = INCX
      CALL DTRMV(UPLO,TRANS,DIAG,N4,A,LDA4,X,INCX4)
      END
#undef DGEMM
      subroutine dgemm6(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BET,C,LDC)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
      INTEGER*4 K4,LDA4,LDB4,LDC4,M4,N4
      K4   = K
      LDA4 = LDA
      LDB4 = LDB
      LDC4 = LDC
      M4   = M
      N4   = N
c      write (2,'(3I5,A,3I5,A,3I5)') M4,N4,K4,' orig MNK ', M,N,K,' d ',
c     & M4-M,N4-N,K4-K
c      if(M4.GT.100000) THEN
c         CALL QUIT("dgemm6: I am afraid I can't do that")
c      end if
      call DGEMM(TRANSA,TRANSB,M4,N4,K4,ALPHA,A,LDA4,B,LDB4,BET,C,LDC4)
c      write(2,*) 'dgemm completed'
      end
#undef DSPEV
      SUBROUTINE DSPEV6( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDZ, N
      DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
      INTEGER*4 N4, LDZ4, INFO4
      N4 = N
      LDZ4 = LDZ
      CALL DSPEV(JOBZ,UPLO,N4,AP,W,Z,LDZ4,WORK,INFO4)
      INFO = INFO4
      END
#undef DSYMM
      SUBROUTINE DSYMM6( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
      INTEGER*4 M4, N4, LDA4, LDB4, LDC4
      M4 = M
      N4 = N
      LDA4 = LDA
      LDB4 = LDB
      LDC4 = LDC
      CALL DSYMM(SIDE, UPLO, M4, N4, ALPHA, A, LDA4, B, LDB4,
     $                   BETA, C, LDC4 )
      END
#undef DSYR2K
      SUBROUTINE DSYR2K6( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
      INTEGER*4 N4, K4, LDA4, LDB4, LDC4
c     
      N4 = N
      K4 = K
      LDA4 = LDA
      LDB4 = LDB
      LDC4 = LDC
      call DSYR2K( UPLO, TRANS, N4, K4, ALPHA, A, LDA4, B, LDB4,
     $             BETA, C, LDC4 )
      End
#undef DSYRK
      SUBROUTINE DSYRK6( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      DOUBLE PRECISION   ALPHA, BETA
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
      INTEGER*4 N4, K4, LDA4, LDC4
c
      N4 = N
      K4 = K
      LDA4 = LDA
      LDC4 = LDC
      call DSYRK( UPLO, TRANS, N4, K4, ALPHA, A, LDA4,
     $                   BETA, C, LDC4 )
      end
#undef DTRMM
      SUBROUTINE DTRMM6( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      INTEGER*4 M4, N4, LDA4, LDB4
c
      M4 = M
      N4 = N
      LDA4 = LDA
      LDB4 = LDB
      call DTRMM( SIDE, UPLO, TRANSA, DIAG, M4, N4, ALPHA, A, LDA4,
     $                   B, LDB4 )
      end
#undef DTRSM
      SUBROUTINE DTRSM6( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c
      INTEGER*4 M4, N4, LDA4, LDB4
      M4 = M
      N4 = N
      LDA4 = LDA
      LDB4 = LDB
      call DTRSM( SIDE, UPLO, TRANSA, DIAG, M4, N4, ALPHA, A, LDA4,
     $            B, LDB4 )
      end
#undef DPOTRF
      SUBROUTINE DPOTRF6( UPLO, N, A, LDA, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      INTEGER*4 N4, LDA4, INFO4
      N4 = N
      LDA4 = LDA
      CALL DPOTRF(UPLO,N4,A,LDA4,INFO4)
      INFO = INFO4
      END
#undef DTRTRI
      SUBROUTINE DTRTRI6( UPLO, DIAG, N, A, LDA, INFO )
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      INTEGER*4  N4, LDA4, INFO4
      N4 = N
      LDA4 = LDA
      CALL DTRTRI(UPLO,DIAG,N4,A,LDA4,INFO4)
      INFO = INFO4
      END
#undef DGESV
      SUBROUTINE DGESV6( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      INTEGER*4 IPIV4(N),N4,NRHS4,LDA4,LDB4,INFO4,I
      N4 = N
      NRHS4 = NRHS
      LDA4 = LDA
      LDB4 = LDB
      CALL DGESV(N4,NRHS4,A,LDA4,IPIV4,B,LDB4,INFO4)
      DO I = 1, N4
         IPIV(I) = IPIV4(I)
      END DO
      INFO = INFO4
      END
#undef DSYEV
      SUBROUTINE DSYEV6( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
      INTEGER*4 N4, LDA4,LWORK4,INFO4
      N4 = N
      LDA4 = LDA
      LWORK4 = LWORK
      CALL DSYEV(JOBZ,UPLO,N4,A,LDA4,W,WORK,LWORK4,INFO4)
      INFO = INFO4
      END
#undef DSYEVX
      SUBROUTINE DSYEVX6( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
     $                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     $                   IFAIL, INFO )
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
      INTEGER*4 N4, LDA4, IL4, IU4, M4,LDZ4,LWORK4
      INTEGER*4 IFAIL4(N), INFO4,I
      N4 = N
      LDA4 = LDA
      IL4 = IL
      IU4 = IU
      M4   = M
      LDZ4 = LDZ
      LWORK4 = LWORK
      CALL DSYEVX(JOBZ,RANGE,UPLO,N4,A,LDA4,VL,VU,IL4,IU4,ABSTOL,
     $     M4,W,Z,LDZ4,WORK,LWORK4,IWORK,IFAIL4,INFO4)
      DO I=1, N4
         IFAIL(I) = IFAIL4(I)
      END DO
      INFO = INFO4
      END
#undef DSYGV
      SUBROUTINE DSYGV6( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     $                  LWORK, INFO )
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
      INTEGER*4 ITYPE4, N4, LDA4, LDB4,LWORK4,INFO4
      ITYPE4 = ITYPE
      N4 = N
      LDA4 = LDA
      LDB4 = LDB
      LWORK4 = LWORK
      CALL DSYGV(ITYPE,JOBZ,UPLO,N4,A,LDA4,B,LDB4,W,WORK,LWORK4,INFO4)
      INFO = INFO4
      END
