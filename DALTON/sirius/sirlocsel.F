C#######################################################################
      SUBROUTINE LOCCTL(CMO,WORK,LWORK)
#include "implicit.h"
      REAL*8 CMO(*)
      REAL*8 WORK(LWORK)
!
!  Localize orbitals 
!
#include "priunit.h"
#include "infloc.h"
#include "inforb.h"
!
      CALL QENTER('LOCCTL')
      CALL TITLER('Orbital localization','*',105)
      KFREE = 1
      LFREE = LWORK
!
! R2 localization
!
C     IF (R2LOC) THEN
C        CALL HEADER('Orbital localization for center',6)
C        WRITE(*,'(5X3F10.6)') LOCCEN
C        CALL LOCR2(CMO,WORK,LWORK)
C     END IF
C
      IF (BOYSEL) THEN
         CALL HEADER('Boys orbital localization ',LUPRI)
         DO I=1,NSYM
            IF (NBOYS(I).GT.0) THEN
               CALL MEMGET('REAL',KCMO,NBOYS(I)*NBAS(I),
     &            WORK,KFREE,LFREE)
               CALL MEMGET('REAL',KCMOLOC,NBOYS(I)*NBAS(I),
     &            WORK,KFREE,LFREE)
               KCMO1=KCMO
               DO J=1,NBOYS(I)
                  IJ=(BOYSORB(J,I)-1)*NBAS(I)
                  CALL DCOPY(NBAS(I),
     &               CMO(ICMO(I)+IJ+1),1,
     &               WORK(KCMO1),1
     &               )
                  KCMO1=KCMO1+NBAS(I)
               END DO
               WRITE(LUPRI,'(A,I5)')
     &            '   Localization of orbitals, symmetry',I
               WRITE(LUPRI,'(/A,10I4)')
     &            '   Input orbitals',(BOYSORB(J,I),J=1,NBOYS(I))
               CALL OUTPUT(WORK(KCMO),1,NBAS(I),1,NBOYS(I),
     &            NBAS(I),NBOYS(I),1,LUPRI)
               CALL BOYS(NBAS(I),NBOYS(I),WORK(KCMO),WORK(KCMOLOC),
     &            WORK(KFREE),LFREE)
               WRITE(LUPRI,'(/A)')
     &            '   Output orbitals'
               CALL OUTPUT(WORK(KCMOLOC),1,NBAS(I),1,NBOYS(I),
     &            NBAS(I),NBOYS(I),1,LUPRI)
C
C Copy back loclized orbitals to CMO
C
               KCMO1=KCMOLOC
               DO J=1,NBOYS(I)
                  IJ=(BOYSORB(J,I)-1)*NBAS(I)
                  CALL DCOPY(NBAS(I),
     &               WORK(KCMO1),1,
     &               CMO(ICMO(I)+IJ+1),1
     &               )
                  KCMO1=KCMO1+NBAS(I)
               END DO
               CALL MEMREL('LOCCTL',WORK,KCMO,KCMO,KFREE,LFREE)
            END IF
         END DO
      END IF
C
C Save localized orbitals
C
      CALL NEWORB('LOCCTL',CMO,.FALSE.)
      CALL QEXIT('LOCCTL')
      END
C#######################################################################
C /* Deck locr2 */
C      SUBROUTINE LOCR2(CMO,WORK,LWORK)
C#include "implicit.h"
C      INTEGER LWORK
C      REAL*8  WORK(LWORK)
C      REAL*8  CMO(*)
CC
CC  Localize orbitals of a given core hole by minimizing average distance
CC  to the center.
CC
C#include "inforb.h"
C#include "inftap.h"
C#include "dummy.h"
C#include "orgcom.h"
C#include "infloc.h"
C#include "wrkrsp.h"
CC
CC  External
CC
C      LOGICAL FNDLB2
C      EXTERNAL FNDLB2
CC
CC  Local
CC
C      CHARACTER*7 RTNLBL(2)
C      INTEGER IPRINT
C      data iprint/11/
CC
C      CALL QENTER('LOCR2')
C      KFREE = 1
C      LFREE = LWORK
C      CALL MEMGET('REAL',KR2,N2ORBX,WORK,KFREE,LFREE)
CC
CC   Calculate quadrupole integrals, if necessary
CC
C      GDIFF=(DIPORG(1)-LOCCEN(1))**2
C     &     +(DIPORG(2)-LOCCEN(2))**2
C     &     +(DIPORG(3)-LOCCEN(3))**2
C      IF (.NOT.FNDLB2('XXQUADRU',RTNLBL,LUPROP)
C     &   .OR. GDIFF .GT.1E-5) THEN
C         CALL DCOPY(3,LOCCEN,1,DIPORG,1)
C         IF (IPRINT.GT.10) THEN
C            CALL HEADER(
C     &         'Recalculate quadrupole integrals with center',
C     &         6
C     &         )
C            PRINT '(3F10.5)',DIPORG
C         END IF
C         BACKSPACE(LUPROP)
C         CALL NEWLB2('EOFLABEL',RTNLBL,LUPROP,LUERR)
C         CALL PR1INT('QUADRUP',WORK(KFREE),LFREE,IDUMMY,
C     &      IDYMMY,.TRUE.,IPRINT.GT.10,0)
C      END IF
CC
CC  Get x2 integrals and copy to r2
CC
C      KSYMOP = 1
C      CALL MEMGET('REAL',KX2,N2ORBX,WORK,KFREE,LFREE)
C      CALL PRPGET('XXQUADRU',CMO,WORK(KR2),WORK(KFREE),LFREE)
C      IF (IPRINT.GT.10) THEN
C         CALL OUTFUL(WORK(KR2),NORBT,NORBT,'XXQUADRU mo basis')
C      END IF
CC
CC  Get y2 integrals and copy to r2
CC
C      D1 = 1.0D0
C      CALL PRPGET('YYQUADRU',CMO,WORK(KX2),WORK(KFREE),LFREE)
C      IF (IPRINT.GT.10) THEN
C         CALL OUTFUL(WORK(KX2),NORBT,NORBT,'YYQUADRU mo basis')
C      END IF
C      CALL DAXPY(N2ORBX,D1,WORK(KX2),1,WORK(KR2),1)
CC
CC  Get z2 integrals and copy to r2
CC
C      CALL PRPGET('ZZQUADRU',CMO,WORK(KX2),WORK(KFREE),LFREE)
C      IF (IPRINT.GT.10) THEN
C         CALL OUTFUL(WORK(KX2),NORBT,NORBT,'ZZQUADRU mo basis')
C      END IF
C      CALL DAXPY(N2ORBX,D1,WORK(KX2),1,WORK(KR2),1)
C      IF (IPRINT.GT.10) THEN
C         CALL OUTFUL(WORK(KR2),NORBT,NORBT,'LOCR2: R2 matrix MO basis')
C      END IF
C      CALL MEMREL('LOCR2: after reading X2...',WORK,KR2,KX2,KFREE,LFREE)
CC
CC  Diagonalize
CC
C      NISH1=NISH(1)
C      NBAS1=NBAS(1)
C      CALL MEMGET('REAL',KR2II,NISH1*NISH1,WORK,KFREE,LFREE)
C      CALL MCOPY(NISH1,NISH1,WORK(KR2),NORBT,WORK(KR2II),NISH1)
C      CALL MEMGET('REAL',KW,NISH1,WORK,KFREE,LFREE)
C      CALL DSYEV(
C     &   'V','U',NISH1,WORK(KR2II),NISH1,WORK(KW),
C     &   WORK(KFREE),LFREE,INFO
C     &   )
C      IF (INFO.NE.0) THEN
C         CALL QUIT('LOCR2:SSYEV diagonalization failed')
C      END IF
C      IF (IPRINT.GT.10) THEN
C         CALL OUTFUL(WORK(KR2II),NISH1,NISH1,'LOCR2:eigenvectors')
C         CALL OUTFUL(WORK(KW),1,NISH1,'LOCR2:eigenvalues ')
C      END IF
C
CC
CC  Use the eigenvectors of R2 to transform occupied molecular orbitals
CC  of symmetry one
CC
C      KTMP=KR2
C      CALL MPAB(
C     &   CMO,        NBAS1,NISH1,NBAS1,NISH1,
C     &   WORK(KR2II),NISH1,NISH1,NISH1,NISH1,
C     &   WORK(KTMP),NBAS1,NISH1
C     &   )
C      IOFF = NBAS1*NISH1
C      CALL DCOPY(
C     &   NCMOT-IOFF,
C     &   CMO(1+IOFF),1,
C     &   WORK(KTMP+IOFF),1
C     &   )
C
C      IF (IPRINT.GT.1) THEN
C         CALL HEADER('LOCR2:Input molecular orbitals',6)
C         CALL PRORB(CMO,IPRINT.LT.10,6)
C         CALL HEADER('LOCR2:Output molecular orbitals',6)
C         CALL PRORB(WORK(KTMP),IPRINT.LT.10,6)
C      END IF
C      CALL DCOPY(NCMOT,WORK(KTMP),1,CMO,1)
C      CALL MEMREL('LOCR2',WORK,1,1,KFREE,LFREE)
C      CALL QEXIT('LOCR2')
C      END
C
C /* Deck boys */
      SUBROUTINE BOYS(NAO,NMO,CMOIN,CMOOUT,WORK,LWORK)
#ifdef __MAIN__
      IMPLICIT NONE
#else
#include "implicit.h"
#endif
      INTEGER NAO,NMO,LWORK
      REAL*8 CMOIN(NAO,NMO), CMOOUT(NAO,NMO), WORK(LWORK)
      REAL*8 THRESB
      REAL*8 THRESNR
      DATA THRESB /1D-8/, THRESNR/1.0D0/
      INTEGER KFREE,LFREE,KDAO,KDMO,KIJ,KDELTA,KKAPPA,KC,KD,KTMP,LTMP,
     &   NZ,KG,KH
#ifdef __MAIN__
      INTEGER LUPRI
      LUPRI=6
#else
#include "priunit.h"
#endif
C
C Input CMOIN
C Output CMOOUT: transformed to minimize  sum(i=1,nmo) (x_ii)^2
C
      CALL QENTER('BOYS')
      CALL DCOPY(NAO*NMO,CMOIN,1,CMOOUT,1)
C
      KFREE=1
      LFREE=LWORK
      CALL MEMGET('REAL',KDAO,NAO*NAO,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMO,3*NAO*NAO,WORK,KFREE,LFREE)
      NZ=NMO*(NMO-1)/2
      CALL MEMGET('INTE',KIJ,NZ*2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDELTA,NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KKAPPA,NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KC,NMO*NMO,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KD,NMO*NMO,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KG,NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH,NZ*NZ,WORK,KFREE,LFREE)
C
      CALL BOYS1(
     &   NAO,NMO,NZ,WORK(KIJ),
     &   CMOOUT,WORK(KDAO),WORK(KDMO),
     &   WORK(KDELTA),WORK(KKAPPA),WORK(KC),WORK(KD),
     &   WORK,KFREE,LFREE,
     &   THRESB,THRESNR,
     &   WORK(KG),WORK(KH),
     &   LUPRI
     &   )
      CALL MEMREL('BOYS',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('BOYS')
      END
C
      SUBROUTINE BOYS1(
     &   NAO,NMO,NZ,IJ,
     &   CMO,DAO,X,DELTA,KAPPA,C,D,WORK,KFREE,LFREE,
     &   THRESB,THRESNR,
     &   G0,H0,
     &   LUPR
     &   )
      IMPLICIT NONE
C
C Input
C
      INTEGER NAO,NMO,NZ,IJ(NZ,2)
      REAL*8 THRESB, THRESNR
      INTEGER LUPR
C
C Input/output
C
      REAL*8 CMO(NAO,NMO)
C
C Scratch
C
      REAL*8 DAO(NAO,NAO),X(NMO,NMO,3),DELTA(NZ),KAPPA(NZ),C(NMO,NMO),
     &   D(NMO,NMO), G0(NZ),H0(NZ,NZ), WORK(*)
C
C Local
C
      CHARACTER*8 DLAB(3)
      DATA DLAB /'XDIPLEN ','YDIPLEN ','ZDIPLEN '/
      REAL*8 F0, D1, D0, GNORM, DELTA0, DELTA1, KAPPA0, KAPPA1,
     &   SIN4K, COS4K, TAN4K, XKL, DKL, THRESROT , KAPPAM, PI
      INTEGER IMAX, LU, NAOT, I, J, ITER, K, L, KTMP,LTMP,KFREE,LFREE
      INTEGER HSIG(3), MAXIT
      DATA THRESROT /1D-10/ , MAXIT /30/
      PARAMETER ( D1=1.0D0,  D0=0.0D0)
      PARAMETER (PI = 3.14159 26535 89793 23846 D00)
      LOGICAL ISMAX
C
C External funcitons (blas)
C
      INTEGER IDAMAX 
      EXTERNAL IDAMAX
      CALL QENTER('BOYS1')
C
C Set up index matrix
C
      I=0
      DO K=2,NMO
         DO L=1,K-1
            I=I+1
            IJ(I,1)=K
            IJ(I,2)=L
         END DO
      END DO
      CALL MEMCHK('IJ',WORK,1)
C
C Read dipole matrix and transform to selected MO:s
C
      LU=-1
      CALL GPOPEN(LU,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',0,
     &   .FALSE.)
      NAOT=NAO*(NAO+1)/2
      LTMP=MAX(NAOT,NMO*NAO)
      CALL MEMGET('REAL',KTMP,LTMP,WORK,KFREE,LFREE)
      DO I=1,3
         REWIND LU
         CALL MOLLAB(DLAB(I),LU,LUPR)
         CALL READT(LU,NAOT,WORK(KTMP))
C        CALL OUTPAK(WORK(KTMP),NAO,1,LUPR)
         CALL DSPTSI(NAO,WORK(KTMP),DAO)
C        WRITE (LUPR,'(A)')DLAB(I)//'(AO)'
C        CALL OUTPUT(DAO,1,NAO,1,NAO,NAO,NAO,1,LUPR)
         CALL DGEMM('T','N',NMO,NAO,NAO,
     &      D1,CMO,NAO,
     &         DAO,NAO,
     &      D0,WORK(KTMP),NMO
     &      )
         CALL DGEMM('N','N',NMO,NMO,NAO,
     &      D1,WORK(KTMP),NMO,
     &         CMO,NAO,
     &      D0,X(1,1,I),NMO
     &      )
C        WRITE (LUPR,'(A)')DLAB(I)//'(AO)'
C        CALL OUTPUT(X(1,1,I),1,NMO,1,NMO,NMO,NMO,1,LUPR)
      END DO
      CALL MEMCHK('CXC(T)',WORK,1)
      CALL GPCLOSE(LU,'KEEP')
C
C Initialize for localization iterations
C
      CALL DUNIT(D,NMO)
      ITER=0
      CALL LOCSUM(NMO,NZ,IJ,D,X,F0,G0,H0,HSIG,ISMAX,GNORM,
     &   WORK,KFREE,LFREE)
      WRITE (LUPR,'(/3X,A5,A,3A14)') 'Iter',':','Boys sum',
     &   'Gradient norm','   Hessian signature'
      WRITE (LUPR,'(3X,I5,A,2F14.8,A5,3I3,A)')
     &   ITER,':',F0,GNORM,'(',HSIG,')'
      DO WHILE ((.NOT. ISMAX.OR.ISMAX.AND.GNORM.GT.THRESB) 
     &   .AND. ITER.LT.MAXIT)
         ITER=ITER+1
C
C Close to convergence to some Newton-Raphson steps
C
         IF (ISMAX.AND.GNORM .LT. THRESNR) THEN
             call memchk('before nr',work,1)
             CALL NR(NMO,NZ,G0,H0,C,IJ,WORK,KFREE,LFREE)
             call memchk('after nr',work,1)
         ELSE
C        
C Loop over individual pairs and find largest rotation
C
            CALL DZERO(DELTA,NZ)
            DO I=1,NZ
               K=IJ(I,1)
               L=IJ(I,2)
               SIN4K=0D0
               COS4K=0D0
               DO J=1,3
                  XKL=X(K,L,J)
                  DKL=X(K,K,J)-X(L,L,J)
                  SIN4K=SIN4K + 4*XKL*DKL
                  COS4K=COS4K + DKL**2 - 4*XKL**2
               END DO
               IF (ABS(COS4K).LT.THRESROT) THEN
                  DELTA(I)=D0
                  KAPPA(I)=D0
               ELSE
                  TAN4K=SIN4K/COS4K
C
C Two solutions to the equation for the max effect of a k-l rotation
C
                  KAPPA0=ATAN(TAN4K)/4
                  KAPPA1=(PI+ATAN(+TAN4K))/4
C
C The actual change for these two rotations
C
                  DELTA0=D0
                  DELTA1=D0
                  DO J=1,3
                     XKL=X(K,L,J)
                     DKL=X(K,K,J)-X(L,L,J)
                     DELTA0 = DELTA0 + SIN(4*KAPPA0)*XKL*DKL
     &                  - (1-COS(4*KAPPA0))/4*(DKL**2-4*XKL**2)
                     DELTA1 = DELTA1 + SIN(4*KAPPA1)*XKL*DKL
     &                  - (1-COS(4*KAPPA1))/4*(DKL**2-4*XKL**2)
                  END DO
                  IF (DELTA1 .GT. DELTA0) THEN
                     DELTA(I)=DELTA1
                     KAPPA(I)=KAPPA1
                  ELSE
                     DELTA(I)=DELTA0
                     KAPPA(I)=KAPPA0
                  END IF
               END IF
            END DO
C
C Done with all pair rotation, carry out the transformation for max increase
C in the localisation sum
C
            IMAX=IDAMAX(NZ,DELTA,1)
            K=IJ(IMAX,1)
            L=IJ(IMAX,2)
            KAPPAM=KAPPA(IMAX)
            CALL DUNIT(C,NMO)
            C(K,K)=COS(KAPPAM)
            C(L,L)=C(K,K)
            C(K,L)=SIN(KAPPAM)
            C(L,K)=-C(K,L)
         END IF
C
C Transform the dipoles and MOs with C obtained from either the NR or the pair 
C rotation 
C
         DO I=1,3
            ! X = C*X*C(T)
            CALL DGEMM('N','N',NMO,NMO,NMO,
     &         D1,C,NMO,
     &            X(1,1,I),NMO,
     &         D0,WORK(KTMP),NMO
     &         )
            CALL DGEMM('N','T',NMO,NMO,NMO,
     &         D1,WORK(KTMP),NMO,
     &            C,NMO,
     &         D0,X(1,1,I),NMO
     &            )
         END DO
         CALL DGEMM('N','T',NAO,NMO,NMO,
     &      D1,CMO,NAO,
     &         C,NMO,
     &      D0,WORK(KTMP),NAO
     &      )
         CALL DCOPY(NAO*NMO,WORK(KTMP),1,CMO,1)
         CALL LOCSUM(NMO,NZ,IJ,
     &      D,X,F0,G0,H0,HSIG,ISMAX,GNORM,WORK,KFREE,LFREE)
         WRITE (LUPR,'(3X,I5,A,2F14.8,A5,3I3,A)')
     &      ITER,':',F0,GNORM,'(',HSIG,')'
      END DO !WHILE
      CALL MEMREL('BOYS1:TMP',WORK,KTMP,KTMP,KFREE,LFREE)
      CALL QEXIT('BOYS1')
      END
C
      SUBROUTINE LOCSUM(NMO,NZ,IJ,D,X,F0,G0,H0,HSIG,ISMAX,GNORM,W,KW,LW)
      IMPLICIT NONE
C
C Input
C
      INTEGER NMO, NZ, IJ(NZ,2)
      REAL*8 X(NMO,NMO,3),D(NMO,NMO)
C
C Output
C
      INTEGER HSIG(3)
      REAL*8 F0, G0(NZ), H0(NZ,NZ)
      LOGICAL ISMAX
      REAL*8 GNORM
C
C I/O KW,LW
C
      INTEGER KW,LW
C
C Scratch
C
      REAL*8 W(*)
C
C External
C
      REAL*8 DNRM2
      EXTERNAL DNRM2
C
C Local
C
      INTEGER I,J,K,L,M,N,IZ,JZ,KH,INFO,KHCOPY
      REAL*8 D0, D1, D2
      PARAMETER (D0=0.0D0, D1=0.0D0, D2=2.0D0)
      REAL*8 H2
      EXTERNAL H2
C
      CALL QENTER('LOCSUM')
      call memchk('in locsum',W,1)
      F0=D0
      CALL DZERO(G0,NZ)
      CALL DZERO(H0,NZ*NZ)
      call memchk('dzero',W,1)
C
C Localization sum
C
      DO J=1,NMO
         DO I=1,3
            F0=F0 + X(J,J,I)**2
         END DO
      END DO
C
C Localization gradient
C
      DO IZ=1,NZ
         K=IJ(IZ,1)
         L=IJ(IZ,2)
         DO I=1,3
            G0(IZ)=G0(IZ) + 2*(X(K,K,I)-X(L,L,I))*X(K,L,I)
         END DO
      END DO
      GNORM=DNRM2(NZ,G0,1)
C
C Localization hessian
C
      CALL DGEMM('N','T',NZ,NZ,1,
     &   D2,G0,NZ,
     &      G0,NZ,
     &   D0,H0,NZ
     &   )
      DO IZ=1,NZ
         K=IJ(IZ,1)
         L=IJ(IZ,2)
         DO JZ=1,NZ
            M=IJ(JZ,1)
            N=IJ(JZ,2)
            DO I=1,3
               H0(IZ,JZ) = H0(IZ,JZ) 
     &         + H2(NMO,X(1,1,I),D,K,L,M,N) - H2(NMO,X(1,1,I),D,K,L,N,M)
            END DO
         END DO
      END DO
C
C Check if hessian eigenvalues negative, set ismax if yes
C
      CALL MEMCHK('locsum:before dsyev',w,1)
      CALL MEMGET('REAL',KH,NZ,W,KW,LW)
      CALL MEMGET('REAL',KHCOPY,NZ*NZ,W,KW,LW)
      CALL DCOPY(NZ*NZ,H0,1,W(KHCOPY),1)
      CALL DSYEV('N','U',NZ,W(KHCOPY),NZ,W(KH),W(KW),LW,INFO)
      IF (INFO.NE.0) STOP 'DSYEV FAILED'
      ISMAX=.TRUE.
      HSIG(1)=0
      HSIG(2)=0
      HSIG(3)=0
      DO I=0,NZ-1
         ISMAX=ISMAX.AND.W(KH+I).LT.D0
         IF (W(KH+I).GT.0 ) THEN
            HSIG(1)=HSIG(1)+1
         ELSE IF (W(KH+I).LT.0 ) THEN
            HSIG(3)=HSIG(3)+1
         ELSE
            HSIG(2)=HSIG(2)+1
         END IF
      END DO
      CALL MEMREL('LOCSUM',W,KH,KH,KW,LW)
      CALL QEXIT('LOCSUM')
      END
C
      REAL*8 FUNCTION H2(NMO,X,D,R,S,T,U)
      IMPLICIT NONE
      INTEGER NMO
      REAL*8 X(NMO,NMO),D(NMO,NMO)
      INTEGER R,S,T,U
      H2=(X(R,U)*D(T,S) - X(T,S)*D(R,U))
     &      *(X(R,R)+X(U,U)-X(T,T)-X(S,S))
     &   + 2*X(R,S)*X(T,U)*(D(R,T)-D(R,U)-D(S,T)+D(S,U))
      END
C
      SUBROUTINE NR (NMO,NZ,G,H,C,IJ,W,KW,LW)
      IMPLICIT NONE
      INTEGER NMO,NZ,KW,LW
      REAL*8 G(NZ), H(NZ,NZ), C(NMO,NMO), W(*)
      INTEGER IJ(NZ,2)
C
      REAL*8 D1
      PARAMETER (D1=1.0D0)
      INTEGER I, K, L, INFO, IPIV
      integer kg,kh
C
      CALL QENTER('NR')
C
C Solve   Hx=-g
C
      CALL DSCAL(NZ,-D1,G,1)
      CALL MEMGET('INTE',IPIV,NZ,W,KW,LW)
      call memget('REAL',kg,nz,w,kw,lw)
      call memget('REAL',kh,nz*nz,w,kw,lw)
      call dcopy(nz,g,1,w(kg),1)
      call dcopy(nz*nz,h,1,w(kh),1)
      CALL DGESV(NZ,1,H,NZ,W(IPIV),G,NZ,INFO)
      IF (INFO.NE.0) STOP  'NR FAILED'
      CALL MEMREL('NR:dgesv',W,IPIV,IPIV,KW,LW)

C
C G now contains the solution vector, unpack to C
C
      DO I=1,NZ
         K=IJ(I,1)
         L=IJ(I,2)
         C(K,L)=G(I)
         C(L,K)=-G(I)
      END DO
C
C Form the exponential operator
C
      CALL MEXP(NZ,NMO,C,W,KW,LW)
      CALL QEXIT('NR')
      END 
C
      SUBROUTINE MEXP(NZ,NMO,C,W,KW,LW)
      IMPLICIT NONE
C
C Input NZ,NMO
C
C Input/output C, on exit EXP(C)
C
C
C Scratch complex W
C
      INTEGER NZ,NMO,KW,LW
      REAL*8 C(NMO,NMO)
      COMPLEX*16 W(*)
C
C Local
C
      INTEGER INFO
      DOUBLE PRECISION D1,D0
      PARAMETER (D1=1.0D0, D0=0.0D0)
      DOUBLE COMPLEX C1,C0,J, Z
      PARAMETER (C1=(1D0,0D0), C0=(0D0,0D0), J=(0D0,1D0))
C
      INTEGER IVR,IC,IW,I,RW,EW,K,L
C
      CALL QENTER('MEXP')
C     CALL MEMGET('COMP',IVR,NMO*NMO,W,KW,LW)
C     CALL MEMGET('COMP',IC ,NMO*NMO,W,KW,LW)
      CALL MEMGET('REAL',IVR,2*NMO*NMO,W,KW,LW)
      CALL MEMGET('REAL',IC ,2*NMO*NMO,W,KW,LW)
      CALL MEMGET('REAL',IW ,2*NMO*NMO,W,KW,LW)
C
C Hermition (pure imaginary) copy of C
C
      DO I=1,NMO*NMO
         W(IC+I-1)=J*C(I,1)
      END DO
C
C Eigenvalues of C (now hermitian)
C
      CALL MEMGET('REAL',RW,3*NMO-2,W,KW,LW)
      CALL MEMGET('REAL',EW,NMO*NMO,W,KW,LW)
      call memchk('before zheev',w,1)
      CALL ZHEEV('V','U',NMO,W(IC),NMO,C,W(KW),LW,W(RW),INFO)
      call memchk('after zheev',w,1)
      CALL DZERO(W(EW),2*NMO*NMO)
      DO I=1,NMO
         Z=-J*C(I,1)
         W(EW+I-1+(I-1)*NMO) = EXP(Z)
      END DO
      call memchk('after ivr..',w,1)
C
C Transform with eigenvectors for final matrix T exp(c) T+
C
      CALL ZGEMM('N','N',NMO,NMO,NMO,
     &   C1,W(IC),NMO,
     &      W(EW),NMO,
     &   C0,W(IW),NMO
     &   )
      CALL ZGEMM('N','C',NMO,NMO,NMO,
     &   C1,W(IW),NMO,
     &      W(IC),NMO,
     &   C0,W(IVR),NMO
     &   )
      call memchk('after zgemm',w,1)
C
C Resulting matrix is real
C
      DO I=1,NMO*NMO
         C(I,1) = DBLE(W(IVR + I -1))
      END DO
      CALL MEMREL('NR',W,IVR,IVR,KW,LW)
      CALL QEXIT('MEXP')
      END
C
C
C
#ifdef __MAIN__
      PROGRAM BOYSTEST
      IMPLICIT NONE
      INTEGER NAO,NMO,LUSIR,LUPRI,LWORK,NSYM,NORB(8),NBAS(8),KW,LW,
     &   KC1,KC2
      PARAMETER (LWORK=1000)
      REAL*8 WORK(LWORK)
      LOGICAL FNDLAB
      EXTERNAL FNDLAB
      LUSIR=1
      LUPRI=6
      OPEN(LUSIR,FILE='SIRIUS.RST',STATUS='OLD',FORM='UNFORMATTED')
      REWIND LUSIR
      CALL MOLLAB('BASINFO ',LUSIR,LUPRI)
      READ(LUSIR) NSYM,NBAS,NORB
      NAO=NBAS(1)
      NMO=NORB(1)
      KW = 1
      LW = LWORK
      CALL MEMGET('REAL',KC1,NAO*NMO,WORK,KW,LW)
      CALL MEMGET('REAL',KC2,NAO*NMO,WORK,KW,LW)
      CALL MOLLAB('NEWORB  ',LUSIR,LUPRI)
      CALL READT(LUSIR,NAO*NMO,WORK(KC1))
      CLOSE(LUSIR)
      WRITE(LUPRI,'(/A)') 'Input orbitals'
      CALL OUTPUT(WORK(KC1),1,NAO,1,NMO,NAO,NMO,1,6)
      CALL GPIOIN
      CALL BOYS(NAO,NMO,WORK(KC1),WORK(KC2),WORK(KW),LW)
      WRITE(LUPRI,'(/A)') 'Localized orbitals'
      CALL OUTPUT(WORK(KC2),1,NAO,1,NMO,NAO,NMO,1,6)
      CALL MEMREL('BOYSTEST',WORK,1,1,KW,LW)
      END
#endif
C#######################################################################
      BLOCK DATA INFLOC
#include "infloc.h"
      DATA BOYSEL /.FALSE./
      END
