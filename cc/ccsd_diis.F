!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
#if defined (SYS_CRAY)
      SUBROUTINE CCSD_DIIS(S,T,NDIM,NITER)
C
C.... SUBROUTINE TO GET FROM S TO EXTRAPOLATED VECTOR ON T
C.... USES S AS PARAMETER VECTOR AND S-T AS ERROR VECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include <priunit.h>
      PARAMETER (MAXVEC=8)
      COMMON/DICED/IPOSP(MAXVEC+1),IPOSE(MAXVEC+1),A(MAXVEC+1,MAXVEC+1)
      DIMENSION B(MAXVEC+1),C(MAXVEC+1),WK1(MAXVEC+1),WK2(MAXVEC+1),
     +             AA(MAXVEC+1,MAXVEC+1)
C
      DIMENSION S(NDIM)
      DIMENSION T(NDIM)
C
      LEN=NDIM
CRDa      LENGF=LENSEC(LEN)
      IFILS=18
      IFAIL=0
C
      DO 1 I=1,LEN
1     T(I)= S(I) - T(I)
C
      IVEC=NITER
      IF(NITER.GT.MAXVEC)IVEC=NITER- ((NITER-1)/MAXVEC) * MAXVEC
C
      IF(IVEC.EQ.1)IPOSP(IVEC)=1
CRDa      CALL LOKFOR(IPOSP(IVEC),IFILS)
CRDa      CALL WRT3S(S,LEN,IFILS)
CRDA      IPOSE(IVEC)=IPOSP(IVEC)+LENGF
CRDA  CALL WRT3S(T,LEN,IFILS)
      call putwa(ifils,s,iposp(ivec),len)
      ipose(ivec)=iposp(ivec)+len
      call putwa(ifils,t,ipose(ivec),len)
      iposp(ivec+1)=ipose(ivec)+len
C     IPOSP(IVEC+1)=IPOSE(IVEC)+LENGF
C
      NLOOP=NITER
      IF(NITER.GT.MAXVEC)NLOOP=MAXVEC
C
      DO 2 I=1,NLOOP
CRDA      CALL READ3(S,LEN,IPOSE(I),IFILS)
      call getwa(ifils,s,ipose(i),len)
CRDA      ERR=SUMVEC(T,S,LEN)
      ERR=SDOT(LEN,T,1,S,1)
      A(I,IVEC)=ERR
      A(IVEC,I)=ERR
      A(I,NLOOP+1)=-1D0
      A(NLOOP+1,I)=-1D0
2     B(I)=0.D0
      A(NLOOP+1,NLOOP+1)=0.D0
      B(NLOOP+1)=-1D0
      NSIZ=NLOOP+1
C
      CALL F04ATF(A,MAXVEC+1,B,NSIZ,C,AA,MAXVEC+1,WK1,WK2,IFAIL)
C
      IF(IFAIL.NE.0)WRITE(LUPRI,*)  'IFAIL ',IFAIL
C
C
      SUM=0.D0
      DO 22 I=1,NLOOP
22    SUM=SUM+C(I)
C
      CALL DZERO(T,LEN)
      DO 4 I=1,NLOOP
CRDA      CALL READ3(S,LEN,IPOSP(I),IFILS)
      call getwa(ifils,s,iposp(i),len)
      DO 5 J=1,LEN
5     T(J)=T(J)+C(I)*S(J)
4     CONTINUE
      RETURN
      END
      SUBROUTINE F04ATF(A,IA,B,N,C,AA,IAA,WKS1,IPVT,IFAIL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(IA,N),B(N),C(N),AA(IAA,N),WKS1(N),IPVT(N)
C
C     THIS IS NOT THE NAG ROUTINE F04ATF, THOUGH THE
C     ARGUMENTS ARE THE SAME
C     IT IS A (HOPEFULLY) LESS MACHINE DEPENDENT REPLACEMENT
C     BASED ON VECTOR ALGORITHM OF DONGARRA USING
C     LU DECOMPOSITION
C     SEE ACM TRANSACTIONS OF MATHEMATICAL SOFTWARE
C     VOL 10 SEPT 1984
C
      DO 10 I=1,N
      DO 10 J=1,N
10    AA(J,I)=A(J,I)
      CALL LU(AA,IAA,N,IPVT,INFO)
      CALL LUS(AA,IAA,N,IPVT,C,B)
      IFAIL=0
      RETURN
      END
C
C
C
      SUBROUTINE LUS(A,LDA,N,IPVT,X,B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include <priunit.h>
      DIMENSION A(LDA,1),X(1),B(1),IPVT(1)
C
C     SIMULTANEOUS EQUATION AX=B GIVEN PRIOR LU
C     DECOMPOSITION OF A
C     TAKEN FROM DONGARRA ACM TRANS MATH SOFTWARE VOL 10
C
      DO 10 K=1,N
10    X(K)=B(K)
      DO 20 K=1,N
      L=IPVT(K)
      XK=X(L)
      X(L)=X(K)
20    X(K)=XK
      DO 40 K=1,N
      XK=X(K)*A(K,K)
      IF(XK.NE.0.0D0)THEN
      DO 30 I=K+1,N
30    X(I)=X(I)-A(I,K)*XK
      ENDIF
40    X(K)=XK
      DO 60 K=N,1,-1
      XK=X(K)
      IF(XK.NE.0.0D0)THEN
      DO 50 I=1,K-1
50    X(I)=X(I)+A(I,K)*XK
      ENDIF
60    CONTINUE
      RETURN
      END
      SUBROUTINE LU(A,LDA,N,IPVT,INFO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LDA,*),IPVT(*)
C
C     LU DECOMPOSITION OF A
C
      INFO=0
      DO 40 J=1,N
      CALL SMXPY(N-J+1,A(J,J),J-1,LDA,A(1,J),A(J,1))
C
C     FIND PIVOT
C     NOTE : THIS IS ISAMAX OPERATION
C     (USE CRAY ROUTINE OR VECTOR MERGE)
C
      T=DABS(A(J,J))
      K=J
      DO 10 I=J+1,N
      IF(DABS(A(I,J)).GT.T)THEN
      T=DABS(A(I,J))
      K=I
      ENDIF
10    CONTINUE
      IPVT(J)=K
C
      IF(T.EQ.0.0D0)THEN
      WRITE(LUPRI,*)'ZERO PIVOT - MATRIX SINGULAR'
      CALL QUIT(' ')
      ENDIF
C
C     SWOP ROWS
C     NOTE : THIS IS SSWAP OPERATION ON CRAY
C
      DO 20 I=1,N
      T=A(J,I)
      A(J,I)=A(K,I)
      A(K,I)=T
20    CONTINUE
C
      A(J,J)=1.0D0/A(J,J)
      CALL SXMPY(N-J,LDA,A(J,J+1),J-1,LDA,A(J,1),LDA,A(1,J+1))
      T=-A(J,J)
      DO 30 I=J+1,N
30    A(J,I)=T*A(J,I)
40    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE SMXPY(N1,Y,N2,LDA,X,A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),X(*),A(LDA,*)
C
C     Y = Y +A X
C
      J=MOD(N2,2)
      IF(J.GE.1)THEN
      DO 10 I=1,N1
10    Y(I)=(Y(I))+X(J)*A(I,J)
      ENDIF
      J=MOD(N2,4)
      IF(J.GE.2)THEN
      DO 20 I=1,N1
20    Y(I)=((Y(I))+X(J-1)*A(I,J-1))+X(J)*A(I,J)
      ENDIF
      JMIN=J+4
      DO 40 J=JMIN,N2,4
      DO 30 I=1,N1
      Y(I)=((((Y(I))+X(J-3)*A(I,J-3))+X(J-2)*A(I,J-2))
     1    +X(J-1)*A(I,J-1))+X(J)*A(I,J)
30    CONTINUE
40    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE SXMPY(N1,LDY,Y,N2,LDX,X,LDA,A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(LDY,*),X(LDX,*),A(LDA,*)
C
C
      J=MOD(N2,2)
      IF(J.GE.1)THEN
      DO 10 I=1,N1
10    Y(1,I)=(Y(1,I))+X(1,J)*A(J,I)
      ENDIF
      J=MOD(N2,4)
      IF(J.GE.2)THEN
      DO 20 I=1,N1
20    Y(1,I)=((Y(1,I))+X(1,J-1)*A(J-1,I))+X(1,J)*A(J,I)
      ENDIF
      JMIN=J+4
      DO 40 J=JMIN,N2,4
      DO 30 I=1,N1
      Y(1,I)=((((Y(1,I))+X(1,J-3)*A(J-3,I))+X(1,J-2)*A(J-2,I))
     1    +X(1,J-1)*A(J-1,I))+X(1,J)*A(J,I)
30    CONTINUE
40    CONTINUE
      RETURN
      END
#endif
#if !defined (SYS_CRAY)
      SUBROUTINE CCSD_DIIS(S,T,NDIM,NITER)
C
C
C     Written by Poul Joegensen and Henrik Koch 3-Aug-1995
C
C     The S vector contain the perturbation estimate and T the previous
C     amplitudes. On exit the T contain the diis estimate of the new
C     amplitudes.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
      PARAMETER (MAXVEC=80)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C
      LOGICAL FIRST
C
      DIMENSION S(NDIM)
      DIMENSION T(NDIM)
C
      COMMON/LUDIIS/ LUTDIS, LUSDIS
      COMMON/DICED/A(MAXVEC+1,MAXVEC+1)
      DIMENSION B(MAXVEC+1),C(MAXVEC+1),WK1(MAXVEC+1),WK2(MAXVEC+1),
     +          AA(MAXVEC+1,MAXVEC+1)
cKeld      DIMENSION B(MAXVEC+1),C(MAXVEC+1),WK1(MAXVEC+1),WK2(MAXVEC+1),
cKeld     +          A(MAXVEC+1,MAXVEC+1),AA(MAXVEC+1,MAXVEC+1)
C
#include <ccsdinp.h>
C
      SAVE FIRST
      DATA FIRST /.TRUE./
C
      IF (FIRST) THEN
C
         FIRST = .FALSE.
         LUTDIS = -9003
         LUSDIS = -9004
C
         CALL GPOPEN(LUTDIS,'DIIST','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         CALL GPOPEN(LUSDIS,'DIISS','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
C
      ENDIF
C
      DO 100 I=1,NDIM
         T(I)= S(I) - T(I)
  100 CONTINUE
C
      IVEC=NITER
C
      IF (NITER .GT. MXDIIS) THEN
         IVEC = NITER - ((NITER-1)/MXDIIS)*MXDIIS
      ENDIF
C
      IF (IVEC.EQ.1) THEN
         IF (LUTDIS .LE. 0) CALL GPOPEN(LUTDIS,'DIIST','UNKNOWN',' ',
     *               'UNFORMATTED',IDUMMY,.FALSE.)
         IF (LUSDIS .LE. 0) CALL GPOPEN(LUSDIS,'DIISS','UNKNOWN',' ',
     *               'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND(LUTDIS)
         REWIND(LUSDIS)
      ENDIF
C
      WRITE(LUTDIS) S
      WRITE(LUSDIS) T
C
      NLOOP=NITER
      IF (NITER .GT. MXDIIS) THEN
         NLOOP = IVEC
      ENDIF
C
      REWIND(LUSDIS)
C
      DO 200 I = 1,NLOOP
C
         READ(LUSDIS) S
C
         ERR=DDOT(NDIM,T,1,S,1)
C
         A(I,IVEC)    =  ERR
         A(IVEC,I)    =  ERR
         A(I,NLOOP+1) = -ONE
         A(NLOOP+1,I) = -ONE
         B(I)         = ZERO
C
  200 CONTINUE
C
      A(NLOOP+1,NLOOP+1) =  ZERO
      B(NLOOP+1)         = -ONE
C
      NSIZ=NLOOP+1
C
      IFAIL = 0
      CALL F04ATF(A,MAXVEC+1,B,NSIZ,C,AA,MAXVEC+1,WK1,WK2,IFAIL)
C
      IF (IFAIL .NE. 0) THEN
        WRITE(LUPRI,*)  'IFAIL ',IFAIL
      ENDIF
C
      CALL DZERO(T,NDIM)
C
      REWIND(LUTDIS)
C
      DO 300 I = 1,NLOOP
C
         READ(LUTDIS) S
C
         CALL DAXPY(NDIM,C(I),S,1,T,1)
C
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE LUS(A,LDA,N,IPVT,X,B)
#include <implicit.h>
      DIMENSION A(LDA,1),X(1),B(1),IPVT(1)
C
C     SIMULTANEOUS EQUATION AX=B GIVEN PRIOR LU
C     DECOMPOSITION OF A
C     TAKEN FROM DONGARRA ACM TRANS MATH SOFTWARE VOL 10
C
      DO 10 K=1,N
10    X(K)=B(K)
      DO 20 K=1,N
      L=IPVT(K)
      XK=X(L)
      X(L)=X(K)
20    X(K)=XK
      DO 40 K=1,N
      XK=X(K)*A(K,K)
      IF(XK.NE.0.D00)THEN
      DO 30 I=K+1,N
30    X(I)=X(I)-A(I,K)*XK
      ENDIF
40    X(K)=XK
      DO 60 K=N,1,-1
      XK=X(K)
      IF(XK.NE.0.D00)THEN
      DO 50 I=1,K-1
50    X(I)=X(I)+A(I,K)*XK
      ENDIF
60    CONTINUE
      RETURN
      END
      SUBROUTINE LU(A,LDA,N,IPVT,INFO)
#include <implicit.h>
#include <priunit.h>
      DIMENSION A(LDA,*),IPVT(*)
C
C     LU DECOMPOSITION OF A
C
      INFO=0
      DO 40 J=1,N
      CALL SMXPY(N-J+1,A(J,J),J-1,LDA,A(1,J),A(J,1))
C
C     FIND PIVOT
C     NOTE : THIS IS ISAMAX OPERATION
C     (USE CRAY ROUTINE OR VECTOR MERGE)
C
      T=DABS(A(J,J))
      K=J
      DO 10 I=J+1,N
      IF(DABS(A(I,J)).GT.T)THEN
      T=DABS(A(I,J))
      K=I
      ENDIF
10    CONTINUE
      IPVT(J)=K
C
      IF(T.EQ.0.0D0)THEN
      WRITE(LUPRI,*)'ZERO PIVOT - MATRIX SINGULAR'
      CALL QUIT(' ')
      ENDIF
C
C     SWOP ROWS
C     NOTE : THIS IS SSWAP OPERATION ON CRAY
C
      DO 20 I=1,N
      T=A(J,I)
      A(J,I)=A(K,I)
      A(K,I)=T
20    CONTINUE
C
      A(J,J)=1.0D0/A(J,J)
      CALL SXMPY(N-J,LDA,A(J,J+1),J-1,LDA,A(J,1),LDA,A(1,J+1))
      T=-A(J,J)
      DO 30 I=J+1,N
30    A(J,I)=T*A(J,I)
40    CONTINUE
      RETURN
      END
      SUBROUTINE SMXPY(N1,Y,N2,LDA,X,A)
#include <implicit.h>
      DIMENSION Y(*),X(*),A(LDA,*)
C
C     Y = Y +A X
C
      J=MOD(N2,2)
      IF(J.GE.1)THEN
      DO 10 I=1,N1
10    Y(I)=(Y(I))+X(J)*A(I,J)
      ENDIF
      J=MOD(N2,4)
      IF(J.GE.2)THEN
      DO 20 I=1,N1
20    Y(I)=((Y(I))+X(J-1)*A(I,J-1))+X(J)*A(I,J)
      ENDIF
      JMIN=J+4
      DO 40 J=JMIN,N2,4
      DO 30 I=1,N1
      Y(I)=((((Y(I))+X(J-3)*A(I,J-3))+X(J-2)*A(I,J-2))
     1    +X(J-1)*A(I,J-1))+X(J)*A(I,J)
30    CONTINUE
40    CONTINUE
      RETURN
      END
      SUBROUTINE SXMPY(N1,LDY,Y,N2,LDX,X,LDA,A)
#include <implicit.h>
      DIMENSION Y(LDY,*),X(LDX,*),A(LDA,*)
C
C
      J=MOD(N2,2)
      IF(J.GE.1)THEN
      DO 10 I=1,N1
10    Y(1,I)=(Y(1,I))+X(1,J)*A(J,I)
      ENDIF
      J=MOD(N2,4)
      IF(J.GE.2)THEN
      DO 20 I=1,N1
20    Y(1,I)=((Y(1,I))+X(1,J-1)*A(J-1,I))+X(1,J)*A(J,I)
      ENDIF
      JMIN=J+4
      DO 40 J=JMIN,N2,4
      DO 30 I=1,N1
      Y(1,I)=((((Y(1,I))+X(1,J-3)*A(J-3,I))+X(1,J-2)*A(J-2,I))
     1    +X(1,J-1)*A(J-1,I))+X(1,J)*A(J,I)
30    CONTINUE
40    CONTINUE
      RETURN
      END
      SUBROUTINE F04ATF(A,IA,B,N,C,AA,IAA,WKS1,IPVT,IFAIL)
#include <implicit.h>
      DIMENSION A(IA,N),B(N),C(N),AA(IAA,N),WKS1(N),IPVT(N)
C
C     THIS IS NOT THE NAG ROUTINE F04ATF, THOUGH THE
C     ARGUMENTS ARE THE SAME
C     IT IS A (HOPEFULLY) LESS MACHINE DEPENDENT REPLACEMENT
C     BASED ON VECTOR ALGORITHM OF DONGARRA USING
C     LU DECOMPOSITION
C     SEE ACM TRANSACTIONS OF MATHEMATICAL SOFTWARE
C     VOL 10 SEPT 1984
C
      DO 10 I=1,N
      DO 10 J=1,N
10    AA(J,I)=A(J,I)
      CALL LU(AA,IAA,N,IPVT,INFO)
      CALL LUS(AA,IAA,N,IPVT,C,B)
      IFAIL=0
      RETURN
      END
#endif
