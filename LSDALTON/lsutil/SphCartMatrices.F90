!> @file
!> Module that contains information about angular indeces
!> The Ordering of Cartesian functions are within a shell
!> using d functions as an example
!> 200  D_XX
!> 110  D_XY
!> 101  D_XZ
!> 020  D_YY
!> 011  D_YZ
!> 002  D_ZZ
!> so always sorting after first maximum X then Y then Z
MODULE SphCart_Matrices
use precision
use math_fun
use OverlapType
use memory_handling
CONTAINS
!> \brief calculate spherical transformation matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param SPH_MAT the spherical matrices 
!> \param nSPHMAT the number of spherical matrices
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE setPrecalculatedSphmat(SPH_MAT,nSPHMAT,LUPRI,IPRINT)
implicit none
TYPE(SPHMAT),pointer     :: SPH_MAT(:)
INTEGER                  :: nSPHMAT,LUPRI,IPRINT
!
INTEGER                  :: I,L,nLM,nXYZ

NULLIFY(SPH_MAT)
ALLOCATE(SPH_MAT(nSPHMAT))
DO I=1,nSPHMAT
  L = I-1
  nLM  = 2*L+1
  nXYZ = (L+1)*(L+2)/2
  NULLIFY(SPH_MAT(I)%elms)
  call mem_alloc(SPH_MAT(I)%elms,nLM*nXYZ)
  CALL Sph_to_Cart_matrix(L,SPH_MAT(I)%elms,nLM,nXYZ,LUPRI,IPRINT)
ENDDO

END SUBROUTINE setPrecalculatedSphmat

!> \brief free the calculated spherical transformation matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param SPH_MAT the spherical matrices 
!> \param nSPHMAT the number of spherical matrices
SUBROUTINE freePrecalculatedSphmat(SPH_MAT,nSPHMAT)
implicit none
TYPE(SPHMAT),pointer     :: SPH_MAT(:)
INTEGER                  :: nSPHMAT
!
INTEGER                  :: I

DO I=1,nSPHMAT
  CALL MEM_DEALLOC(SPH_MAT(I)%elms)
  NULLIFY(SPH_MAT(I)%elms)
ENDDO
DEALLOCATE(SPH_MAT)
NULLIFY(SPH_MAT)

END SUBROUTINE freePrecalculatedSphmat

!> \brief calculate spherical to cartesian matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param L angular moment
!> \param SCMAT the transformation matrix
!> \param nLM the number of spherical components
!> \param nXYZ the number of cartesian components
!> \param lupri the logical unit number for the output file
!> \param iprint the printlevel, determining how much output should be generated
SUBROUTINE Sph_to_Cart_matrix(L,SCMAT,nLM,nXYZ,LUPRI,IPRINT)
implicit none
INTEGER,intent(in) :: L,nLM,nXYZ,LUPRI,IPRINT
REAL(REALK),intent(inout) :: SCMAT(nLM,nXYZ)
!
Real(realk), parameter :: DM1 = -1.0E0_realk, D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER     :: M1,MADR,MABS,V0
REAL(realk) :: FACNRM,FAC3,FACTOR
INTEGER     :: T,U,V,A,B,C,P,Q,R,X,Y,Z,M
INTEGER     :: iLM,iXYZ

IF(L .EQ. 0)THEN ! S
 SCMAT(1,1)=1E0_realk
ELSEIF(L .EQ. 1)THEN ! P
 CALL LS_DZERO(SCMAT,9)
 SCMAT(1,1)=1E0_realk
 SCMAT(2,2)=1E0_realk
 SCMAT(3,3)=1E0_realk
ELSEIF(L .GT. 1)THEN
 CALL LS_DZERO(SCMAT,nLM*nXYZ)
 DO M1 = 0, 2*L 
  M = M1 - L
  IF (L.EQ. 1) THEN
    IF (M .EQ. -1) MADR =  0  
    IF (M .EQ.  0) MADR =  1 
    IF (M .EQ.  1) MADR = -1 
  ELSE
    MADR = M
  END IF
  MABS = ABS(M)
  V0 = 0
  IF (M .LT. 0) V0 = 1 
  FACNRM = D1
  IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(LUPRI,L+MABS)*&
                         &FACULT(LUPRI,L-MABS))/(FACULT(LUPRI,L)*(D2**MABS))
  FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
  FACNRM = FACNRM/SQRT(FACUL2(LUPRI,2*L-1))
  DO T = 0, L - MABS, 2
   DO U = 0, T, 2
    DO V = V0, MABS, 2
          !        almost 6.4.48 in the book
     FAC3 = FACNRM*BINOM(LUPRI,L,T/2)*BINOM(LUPRI,L-T/2,MABS+T/2)&
          &                    *BINOM(LUPRI,T/2,U/2)*BINOM(LUPRI,MABS,V)
     DO A = 0, MIN(0,T+MABS-U-V) 
      DO B = 0, MIN(0,U+V)
       DO C = 0, MIN(0,L-T-MABS)
             !           6.4.47 in the book
        DO P = 0, - A, 2
         DO Q = 0, - B, 2
          DO R = 0, - C, 2
           FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
                &   D2**(-A-B-C-P-Q-R-T)*FAC3
           X = T+MABS-U-V-2*A-P
           Y = U+V-2*B-Q
           Z = L-T-MABS-2*C-R
           iLM = 1 + L + MADR
           iXYZ = NCRT(X,Y,Z)
           SCMAT(iLM,iXYZ) = SCMAT(iLM,iXYZ) + FACTOR 
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDIF

END SUBROUTINE Sph_to_Cart_matrix

!!$!> \brief calculate cartesian to spherical matrix
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param L angular moment
!!$!> \param CSMAT the transformation matrix
!!$!> \param nXYZ the number of cartesian components
!!$!> \param nLM the number of spherical components
!!$!> \param lupri the logical unit number for the output file
!!$!> \param iprint the printlevel, determining how much output should be generated
!!$SUBROUTINE Cart_to_Sph_matrix(L,CSMAT,nXYZ,nLM,LUPRI,IPRINT)
!!$implicit none
!!$REAL(REALK) :: CSMAT(nXYZ,nLM)
!!$INTEGER     :: L,nXYZ,nLM,LUPRI,IPRINT
!!$!
!!$Real(realk), parameter :: DM1 = -1.0E0_realk, D1 = 1.0E0_realk, D2 = 2.0E0_realk
!!$INTEGER     :: M1,MADR,MABS,V0
!!$REAL(realk) :: FACNRM,FAC3,FACTOR
!!$INTEGER     :: T,U,V,A,B,C,P,Q,R,X,Y,Z,M
!!$INTEGER     :: iXYZ,iLM
!!$
!!$IF(L .EQ. 0)THEN ! S
!!$  CSMAT(1,1)=1E0_realk
!!$ELSEIF(L .EQ. 1)THEN ! P
!!$  CALL LS_DZERO(CSMAT,9)
!!$  CSMAT(1,1)=1E0_realk
!!$  CSMAT(2,2)=1E0_realk
!!$  CSMAT(3,3)=1E0_realk
!!$ELSEIF(L .GT. 1)THEN
!!$  CALL LS_DZERO(CSMAT,nXYZ*nLM)
!!$  DO M1 = 0, 2*L 
!!$    M = M1 - L
!!$    IF (L.EQ. 1) THEN
!!$      IF (M .EQ. -1) MADR =  0  
!!$      IF (M .EQ.  0) MADR =  1 
!!$      IF (M .EQ.  1) MADR = -1 
!!$    ELSE
!!$      MADR = M
!!$    END IF
!!$    MABS = ABS(M)
!!$    V0 = 0
!!$    IF (M .LT. 0) V0 = 1 
!!$    FACNRM = D1
!!$    IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(LUPRI,L+MABS)*&
!!$                           &FACULT(LUPRI,L-MABS))/(FACULT(LUPRI,L)*(D2**MABS))
!!$    FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
!!$    FACNRM = FACNRM/SQRT(FACUL2(LUPRI,2*L-1))
!!$    DO T = 0, L - MABS, 2
!!$    DO U = 0, T, 2
!!$    DO V = V0, MABS, 2
!!$            !        almost 6.4.48 in the book
!!$      FAC3 = FACNRM*BINOM(LUPRI,L,T/2)*BINOM(LUPRI,L-T/2,MABS+T/2)&
!!$            &                    *BINOM(LUPRI,T/2,U/2)*BINOM(LUPRI,MABS,V)
!!$      DO A = 0, MIN(0,T+MABS-U-V) 
!!$      DO B = 0, MIN(0,U+V)
!!$      DO C = 0, MIN(0,L-T-MABS)
!!$               !           6.4.47 in the book
!!$        DO P = 0, - A, 2
!!$        DO Q = 0, - B, 2
!!$        DO R = 0, - C, 2
!!$          FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
!!$                  &   D2**(-A-B-C-P-Q-R-T)*FAC3
!!$          X = T+MABS-U-V-2*A-P
!!$          Y = U+V-2*B-Q
!!$          Z = L-T-MABS-2*C-R
!!$	  iXYZ = NCRT(X,Y,Z)
!!$	  iLM = 1 + L + MADR
!!$          CSMAT(iXYZ,iLM) = CSMAT(iXYZ,iLM) + FACTOR 
!!$        ENDDO
!!$        ENDDO
!!$        ENDDO
!!$      ENDDO
!!$      ENDDO
!!$      ENDDO
!!$    ENDDO
!!$    ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDIF
!!$
!!$END SUBROUTINE Cart_to_Sph_matrix

!> \brief Does the SPHERICAL TRANSFORMATION of a 5 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TEMP the input 5 dim array i cartesian 
!> \param SPHINT the input 5 dim array i spherical
!> \param ndim1 the size of dimension 1
!> \param ndim2 the size of dimension 2
!> \param nMAT the number of input matrices
!> \param nsphmat the number og output matrices
!> \param nangmom the number of angular moments
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param LUPRI the logical unit number for the output file
SUBROUTINE SPHERICAL_TRANSFORMATION(TEMP,SPHINT,ndim1,ndim2,nMAT,nSPHMAT,nANGMOM,iprint,lupri)
use precision
IMPLICIT NONE
INTEGER     :: ANGMOM,nANGMOM,nLM,nXYZ,lupri,ELM,nMAT,nSPHMAT,iprint,ndim1,ndim2
REAL(REALK) :: TEMP(ndim1,ndim2,1,1,nMAT),SPHINT(ndim1,ndim2,1,1,nSPHMAT),COEF
REAL(REALK),pointer :: TRANSMAT(:)
REAL(REALK),PARAMETER :: D0 = 0E0_realk, D1 = 1E0_realk, D2 = 2E0_realk, D3 = 3.0E0_realk
INTEGER     :: STARTI,STARTJ,I,J,K1,K2,TSIZE
nLM  = 2*nANGMOM + 1
nXYZ = (nANGMOM + 1)*(nANGMOM + 2)/2
!call mem_alloc(TRANSMAT,nXYZ,nLM)
TSIZE=nXYZ*nLM
call mem_alloc(TRANSMAT,nXYZ*nLM)

STARTJ=0
STARTI=0
DO ANGMOM=0,nANGMOM
   nLM  = 2*ANGMOM + 1
   nXYZ = (ANGMOM + 1)*(ANGMOM + 2)/2
   CALL BUILD_CART_TO_SPH_MAT(ANGMOM,TRANSMAT,TSIZE,nLM,nXYZ,lupri,iprint)
   DO I = 1, nLM
      DO J = 1, nXYZ
         COEF = TRANSMAT(J+(I-1)*nXYZ)
         IF (ABS(COEF) .GT. D0) THEN
            DO K2 = 1, ndim2
               DO K1 = 1, ndim1
                  SPHINT(K1,K2,1,1,STARTI+I) = SPHINT(K1,K2,1,1,STARTI+I) + COEF*TEMP(K1,K2,1,1,STARTJ+J)
               ENDDO
            ENDDO
         END IF
      ENDDO
   ENDDO
   STARTI=STARTI+nLM
   STARTJ=STARTJ+nXYZ
ENDDO
!call mem_dealloc(TRANSMAT)
CALL MEM_DEALLOC(TRANSMAT)

END SUBROUTINE SPHERICAL_TRANSFORMATION

!!$!> \brief Does the SPHERICAL TRANSFORMATION of a 2 dim array
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param TEMP the input 2 dim array i cartesian 
!!$!> \param SPHINT the input 2 dim array i spherical
!!$!> \param ELM the size of dimension 1
!!$!> \param nMAT the size of input dimension 2 (number of cartesian comp) 
!!$!> \param nSPHMAT the size of output dimension 2 (number of spherical comp) 
!!$!> \param nangmom the number of angular moments
!!$!> \param IPRINT the printlevel, determining how much output should be generated
!!$!> \param LUPRI the logical unit number for the output file
!!$SUBROUTINE SPHERICAL_TRANSFORMATION2(TEMP,SPHINT,ELM,nMAT,nSPHMAT,nANGMOM,iprint,lupri)
!!$  use precision
!!$IMPLICIT NONE
!!$INTEGER     :: ANGMOM,nANGMOM,nLM,nXYZ,lupri,ELM,nMAT,nSPHMAT,iprint,TSIZE
!!$REAL(REALK) :: TEMP(ELM,nMAT),SPHINT(ELM,nSPHMAT),COEF
!!$REAL(REALK),pointer :: TRANSMAT(:)
!!$REAL(REALK),PARAMETER :: D0 = 0E0_realk, D1 = 1E0_realk, D2 = 2E0_realk, D3 = 3.0E0_realk
!!$INTEGER     :: STARTI,STARTJ,I,J,K
!!$nLM  = 2*nANGMOM + 1
!!$nXYZ = (nANGMOM + 1)*(nANGMOM + 2)/2
!!$!CALL MEM_ALLOC(TRANSMAT,nXYZ*nLM)
!!$TSIZE=nXYZ*nLM
!!$call mem_alloc(TRANSMAT,TSIZE)
!!$STARTJ=0
!!$STARTI=0
!!$DO ANGMOM=0,nANGMOM
!!$   nLM  = 2*ANGMOM + 1
!!$   nXYZ = (ANGMOM + 1)*(ANGMOM + 2)/2
!!$   CALL BUILD_CART_TO_SPH_MAT(ANGMOM,TRANSMAT,TSIZE,nLM,nXYZ,lupri,iprint)
!!$   DO I = 1, nLM
!!$      DO J = 1, nXYZ
!!$         COEF = TRANSMAT(J+(I-1)*nXYZ)
!!$         IF (ABS(COEF) .GT. D0) THEN
!!$            DO K = 1, ELM
!!$               SPHINT(K,STARTI+I) = SPHINT(K,STARTI+I) + COEF*TEMP(K,STARTJ+J)
!!$            ENDDO
!!$         END IF
!!$      ENDDO
!!$   ENDDO
!!$   STARTI=STARTI+nLM
!!$   STARTJ=STARTJ+nXYZ
!!$ENDDO
!!$CALL MEM_DEALLOC(TRANSMAT)
!!$END SUBROUTINE SPHERICAL_TRANSFORMATION2

!> \brief Does the SPHERICAL TRANSFORMATION of a 2 dim array
!> \author T. Kjaergaard
!> \date 2010
!> \param TEMP the input 2 dim array i cartesian 
!> \param SPHINT the input 2 dim array i spherical
!> \param ELM the size of dimension 1
!> \param nMAT the size of input dimension 2 (number of cartesian comp) 
!> \param nSPHMAT the size of output dimension 2 (number of spherical comp) 
!> \param nangmom the number of angular moments
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param LUPRI the logical unit number for the output file
SUBROUTINE SPHERICAL_TRANSFORMATION3(TEMP,SPHINT,ELM,nMAT,nSPHMAT,nANGMOM,iprint,lupri)
  use precision
IMPLICIT NONE
INTEGER     :: ANGMOM,nANGMOM,nLM,nXYZ,lupri,ELM,nMAT,nSPHMAT,iprint,TSIZE
REAL(REALK) :: TEMP(nMAT,ELM),SPHINT(nSPHMAT,ELM),COEF
REAL(REALK),pointer :: TRANSMAT(:)
REAL(REALK),PARAMETER :: D0 = 0E0_realk, D1 = 1E0_realk, D2 = 2E0_realk, D3 = 3.0E0_realk
INTEGER     :: STARTI,STARTJ,I,J,K
nLM  = 2*nANGMOM + 1
nXYZ = (nANGMOM + 1)*(nANGMOM + 2)/2
!CALL MEM_ALLOC(TRANSMAT,nXYZ*nLM)
TSIZE=nXYZ*nLM
call mem_alloc(TRANSMAT,TSIZE)
STARTJ=0
STARTI=0
DO ANGMOM=0,nANGMOM
   nLM  = 2*ANGMOM + 1
   nXYZ = (ANGMOM + 1)*(ANGMOM + 2)/2
   CALL BUILD_CART_TO_SPH_MAT(ANGMOM,TRANSMAT,TSIZE,nLM,nXYZ,lupri,iprint)
   DO I = 1, nLM
      DO J = 1, nXYZ
         COEF = TRANSMAT(J+(I-1)*nXYZ)
         IF (ABS(COEF) .GT. D0) THEN
            DO K = 1, ELM
               SPHINT(STARTI+I,K) = SPHINT(STARTI+I,K) + COEF*TEMP(STARTJ+J,K)
            ENDDO
         END IF
      ENDDO
   ENDDO
   STARTI=STARTI+nLM
   STARTJ=STARTJ+nXYZ
ENDDO
CALL MEM_DEALLOC(TRANSMAT)
END SUBROUTINE SPHERICAL_TRANSFORMATION3

!> \brief build cartesian to spherical matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param ANGMOM the angular momentum
!> \param TRAMAT the output matrix
!> \param TSIZE the dim of TRAMAT
!> \param nLM the number of spherical components
!> \param nXYZ the number of cartesian components
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BUILD_CART_TO_SPH_MAT(ANGMOM,TRAMAT,TSIZE,nLM,nXYZ,lupri,iprint)
use precision
use math_fun
IMPLICIT NONE
INTEGER               :: LVAL1,ANGMOM,nLM,nXYZ,IPRINT,lupri,TSIZE
REAL(REALK),PARAMETER :: D0 = 0E0_realk, D1 = 1E0_realk, D2 = 2E0_realk, D3 = 3.0E0_realk
REAL(REALK)           :: TRAMAT(TSIZE), COSSIN(0:ANGMOM,0:ANGMOM), PL(0:ANGMOM)
INTEGER               :: K,M,IX,I,J,L,ILM,N,IXYZ
REAL(REALK)           :: CM,CMK,CMKIJ,CMKI,FAC
LVAL1 = ANGMOM+1
CALL LS_DZERO(PL(0),LVAL1)
DO K = 0, ANGMOM/2
   FAC = (-1)**K
   PL(ANGMOM-2*K) = (FAC/(2**ANGMOM))&
   &  *BINOM(lupri,ANGMOM,K)*BINOM(lupri,2*(ANGMOM-K),ANGMOM)
ENDDO
!IF(IPRINT .GT. 10)THEN
!   CALL LSHEADER(lupri,'Legendre polynomial')
!   CALL LS_OUTPUT(PL(0:LVAL1),1,1,1,LVAL1,1,LVAL1,1,LUPRI)
!ENDIF
CALL LS_DZERO(COSSIN(0,0),LVAL1*LVAL1)
DO M = 0, ANGMOM
   COSSIN(M,0) = D1
   DO K = 1, M
      COSSIN(M,K) = COSSIN(M-1,K-1)*((-1)**(K-1))
      IF (M .GT. K) COSSIN(M,K) = COSSIN(M,K) + COSSIN(M-1,K)
   ENDDO
ENDDO
!IF(IPRINT .GT. 10)THEN
!   CALL LSHEADER(lupri,'Cosine and sine factors')
!   CALL LS_OUTPUT(COSSIN(0:LVAL1,0:LVAL1),1,LVAL1,1,LVAL1,LVAL1,LVAL1,1,LUPRI)
!ENDIF

! Transformation coefficients
! ---------------------------
!

CALL LS_DZERO(TRAMAT,nXYZ*nLM)
DO M = 0, ANGMOM
   CM = SQRT(D2*FACULT(LUPRI,ANGMOM-M)/FACULT(LUPRI,ANGMOM+M))
   IF (M .EQ. 0) CM = D1
!   IF (MINTEG.EQ. 2) CM = CM/SQRT(FACUL2(2*ANGMOM-1))
   DO K = MOD(ANGMOM - M,2), ANGMOM - M, 2
      IF (M .GT. 0) PL(K) = (K+1)*PL(K+1)
      CMK = CM*PL(K)
      DO I = 0, (ANGMOM - K - M)/2
         CMKI = CMK*BINOM(lupri,(ANGMOM - K - M)/2,I)
         DO J = 0, I
            CMKIJ = CMKI*BINOM(lupri,I,J)
            DO N = 0, M
               IX = ANGMOM - 2*J - M + N
               IX = IX*(IX + 1)/2 + LVAL1 - M - 2*I
!               IF (MORDER .EQ. 0) THEN
                  ILM = MAX(1,2*M + MOD(N,2))
!               ELSE
!                  IF (MOD(N,2) .EQ. 1) THEN
!                     ILM = 1 + ANGMOM - M
!                  ELSE
!                     ILM = 1 + ANGMOM + M
!                  END IF
!               END IF
               TRAMAT(IX+(ILM-1)*nXYZ) = TRAMAT(IX+(ILM-1)*nXYZ) + CMKIJ*COSSIN(M,N) 
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

!IF (IPRINT .GT. 4) THEN
!   CALL LSHEADER(lupri,'Cartesian to spherical transformation matrix')
!   WRITE (LUPRI,'(29X,A,I2)') ' Moment order:',ANGMOM
!   IXYZ = (ANGMOM+1)*(ANGMOM+2)/2
!   ILM  = 2*ANGMOM + 1
!   CALL LS_OUTPUT(TRAMAT,1,IXYZ,1,ILM,NXYZ,NLM,1,LUPRI)
!END IF

END SUBROUTINE BUILD_CART_TO_SPH_MAT

END MODULE SphCart_Matrices
