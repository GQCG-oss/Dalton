!> @file
!> Module containing main exchange-correlation integral driver, and routines to evaluate AOs and electron-densities
!> \brief 
!> \author T. Kjaergaard
!> \date 2009 
MODULE THC_UTIL
use THCgridgenerationmodule
use LSparameters
use memory_handling
use precision
use TYPEDEF
use pseudoinverseMod
use BUILDAOBATCH
private
public :: Get_THC_AO_grid_ngrid, Get_THC_AO_grid_X, get_thc_grid_overlap, &
     & get_thc_grid_overlap_inv, get_thc_grid_overlap2
CONTAINS
!> \brief 
!> \author T. Kjaergaard
!> \date 2015
SUBROUTINE Get_THC_AO_grid_ngrid(LUPRI,IPRINT,SETTING,NBAST,NGRID,&
     & THCradint,THCangint,THCHRDNES,THCNOPRUN,THCTURBO,THCRADIALGRID,&
     & THCZdependenMaxAng,THCPARTITIONING,THC_MIN_RAD_PT)
IMPLICIT NONE
REAL(REALK),intent(in) :: THCradint
INTEGER,intent(in)     :: LUPRI,IPRINT,NBAST,THCangint,THCHRDNES
INTEGER,intent(in)     :: THCTURBO,THCRADIALGRID,THCPARTITIONING,THC_MIN_RAD_PT
logical,intent(in)     :: THCNOPRUN,THCZdependenMaxAng
INTEGER,intent(inout)  :: NGRID
TYPE(LSSETTING) :: SETTING
!
TYPE(BASINF)  :: BAS
INTEGER :: THCGRIDDONE,IPRUNE
THCGRIDDONE = 0
CALL BUILD_BASINF(LUPRI,IPRINT,BAS,SETTING,THCGRIDDONE,.FALSE.)

IPRUNE = 1 !pruning: per default on
IF (THCNOPRUN) IPRUNE = 0

CALL THCGenerateGrid(NBAST,THCradint,THCangint,THCHRDNES,&
     & iprune,BAS%natoms,BAS%X,BAS%Y,BAS%Z,BAS%Charge,&
     & BAS%SHELL2ATOM,BAS%SHELLANGMOM,BAS%SHELLNPRIM,BAS%MAXANGMOM,&
     & BAS%MAXNSHELL,BAS%MXPRIM,BAS%PRIEXP,BAS%PRIEXPSTART,BAS%RSHEL,&
     & THCTURBO,THCRADIALGRID,THCZdependenMaxAng,THCPARTITIONING,&
     & LUPRI,IPRINT,THC_MIN_RAD_PT,NGRID)

CALL FREE_BASINF(BAS)
END SUBROUTINE Get_THC_AO_grid_ngrid

!> \brief 
!> \author T. Kjaergaard
!> \date 2015
SUBROUTINE Get_THC_AO_grid_X(LUPRI,IPRINT,SETTING,NBAST,NGRID,X)
use BUILDAOBATCH
IMPLICIT NONE
INTEGER,intent(in)     :: LUPRI,IPRINT,NBAST,NGRID
TYPE(LSSETTING)        :: SETTING
REAL(REALK),intent(inout):: X(NGRID,NBAST)
!
TYPE(BASINF)  :: BAS
INTEGER :: THCGRIDDONE,spSIZE,L,spSIZE2,IJ,J,I,KCKTA
INTEGER :: MMM,SHELLANGMOMA
REAL(REALK),pointer :: SPHMAT(:)
INTEGER,pointer     :: SPINDEX(:)
integer,pointer :: LVALUE(:,:),MVALUE(:,:),NVALUE(:,:)
call ls_dzero(X,ngrid*nbast)
THCGRIDDONE = 0
CALL BUILD_BASINF(LUPRI,IPRINT,BAS,SETTING,THCGRIDDONE,.FALSE.)

! Setting spherical transformation quantities
IF(BAS%MAXANGMOM .GE. 3)THEN
   spSIZE=0
   DO L=2,BAS%MAXANGMOM-1
      spSIZE=spSIZE+(2*L+1)*(L+1)*(L+2)/2
   ENDDO
   call mem_alloc(SPHMAT,spSIZE)
   CALL LS_DZERO(SPHMAT,spSIZE)
   call mem_alloc(SPINDEX,BAS%MAXANGMOM)
   spsize2 = BAS%maxangmom
   CALL THC_Build_PRECALCULATED_SPHMAT(LUPRI,BAS%MAXANGMOM-1,spSIZE,SPHMAT,SPINDEX)
ELSE
   spSIZE = 9
   call mem_alloc(SPHMAT,spSIZE)
   call mem_alloc(SPINDEX,2)
   SPINDEX(1)=1
   SPINDEX(2)=1
   spsize2 = 2
ENDIF
MMM = BAS%MAXANGMOM*(BAS%MAXANGMOM+1)/2*BAS%MAXANGMOM*(BAS%MAXANGMOM+1)/2
call mem_alloc(LVALUE,MMM,BAS%MAXANGMOM)
call mem_alloc(MVALUE,MMM,BAS%MAXANGMOM)
call mem_alloc(NVALUE,MMM,BAS%MAXANGMOM)
LVALUE=0
MVALUE=0
NVALUE=0
DO SHELLANGMOMA=1,BAS%MAXANGMOM
   KCKTA = SHELLANGMOMA*(SHELLANGMOMA+1)/2
   IJ=0
   DO I = 1,KCKTA
      DO J = 1,I
         IJ=IJ+1
         LVALUE(IJ,SHELLANGMOMA)=SHELLANGMOMA-I
         MVALUE(IJ,SHELLANGMOMA)=I-J
         NVALUE(IJ,SHELLANGMOMA)=J-1
      ENDDO
   ENDDO
ENDDO
CALL THC_AO_grid_X(X,ngrid,nbast,BAS%maxnshell,BAS%CC,BAS%ushells,&
     & BAS%CCSTART,BAS%CCINDEX,BAS%CENT,BAS%PRIEXPSTART,lupri,&
     & BAS%mxprim,BAS%shellangmom,MMM,BAS%maxangmom,BAS%nstart,BAS%priexp,&
     & spsize,sphmat,spsize2,spindex,LVALUE,MVALUE,NVALUE)
call mem_dealloc(SPHMAT)
call mem_dealloc(SPINDEX)
call mem_dealloc(LVALUE)
call mem_dealloc(MVALUE)
call mem_dealloc(NVALUE)
CALL FREE_BASINF(BAS)
call mem_dealloc(THCCOOR)
call mem_dealloc(THCWG)

END SUBROUTINE GET_THC_AO_GRID_X

!X = GAO
SUBROUTINE THC_AO_grid_X(GAO,ngrid,nbast,maxnshell,CC,ushells,CCSTART,CCINDEX,&
     & CENT,PRIEXPSTART,lupri,mxprim,shellangmom,MMM,maxangmom,nstart,priexp,&
     & spsize,sphmat,spsize2,spindex,LVALUE,MVALUE,NVALUE)
implicit none
integer,intent(in) :: LUPRI,MAXNSHELL,nbast,MXPRIM,NGRID
REAL(REALK),intent(inout):: GAO(NGRID,NBAST)
INTEGER,intent(in)  :: SHELLANGMOM(MAXNSHELL)
INTEGER,intent(in)  :: MAXANGMOM,MMM,ushells
TYPE(LSMATRIX),intent(in):: CC(ushells)
INTEGER,intent(in) :: CCSTART(ushells),CCINDEX(MAXNSHELL)
REAL(REALK),intent(in):: CENT(3,MAXNSHELL)
INTEGER,intent(in) :: PRIEXPSTART(MAXNSHELL),NSTART(MAXNSHELL)
REAL(REALK),intent(in) :: PRIEXP(MXPRIM)
integer,intent(in) :: spsize,spsize2
real(realk),intent(in) :: sphmat(spsize)
integer,intent(in) :: spindex(spsize2)
integer,intent(in) :: LVALUE(MMM,MAXANGMOM),MVALUE(MMM,MAXANGMOM),NVALUE(MMM,MAXANGMOM)
CALL II_BLGETSOSTHC(LUPRI,GAO,NBAST,MMM,MAXANGMOM,spSIZE,SPHMAT,SPINDEX,MAXNSHELL,NSTART,SHELLANGMOM,&
     & CENT,PRIEXPSTART,CC,ushells,CCSTART,CCINDEX,MXPRIM,PRIEXP,LVALUE,MVALUE,NVALUE,NGRID)
END SUBROUTINE THC_AO_GRID_X

SUBROUTINE II_BLGETSOSTHC(LUPRI,GAO,NBAST,MMM,MAXANGMOM,SIZE,SPHMAT,SPINDEX,MAXNSHELL,NSTART,&
     & SHELLANGMOM,CENT,PRIEXPSTART,CC,ushells,CCSTART,CCINDEX,MXPRIM,PRIEXP,LVALUE,MVALUE,&
     & NVALUE,NGRID)
IMPLICIT NONE
INTEGER,intent(in) :: LUPRI,SIZE,MAXANGMOM,MAXNSHELL,NGRID,NBAST
INTEGER,intent(in) :: SPINDEX(MAXANGMOM),NSTART(MAXNSHELL),SHELLANGMOM(MAXNSHELL),ushells
INTEGER,intent(in) :: PRIEXPSTART(MAXNSHELL),MXPRIM,CCSTART(ushells),MMM
integer,intent(in) :: LVALUE(MMM,MAXANGMOM),MVALUE(MMM,MAXANGMOM)
integer,intent(in) :: NVALUE(MMM,MAXANGMOM)
INTEGER,intent(in)     :: CCINDEX(MAXNSHELL)
REAL(REALK),intent(in) :: PRIEXP(MXPRIM)
real(realk),intent(in) :: CENT(3,MAXNSHELL)
REAL(REALK),intent(inout) :: GAO(NGRID,NBAST)
REAL(REALK),intent(in) :: SPHMAT(SIZE)
TYPE(LSMATRIX),intent(in) :: CC(ushells)
!
REAL(REALK),pointer :: PA(:,:), PA2(:)
REAL(REALK) :: DFTHRI,CENX,CENY,CENZ
INTEGER     :: I,IADR,JSTA,ISHELA,KHKTA,KCKTA,SHELLANGMOMA,IBL,J,oooo,IJ,SPVAL,K,IBL2
call mem_alloc(PA,3,NGRID)
call mem_alloc(PA2,NGRID)
SPVAL=1
IADR = 0
DO ISHELA = 1,MAXNSHELL
   SHELLANGMOMA = SHELLANGMOM(ISHELA)
   KHKTA  = 2*(SHELLANGMOMA-1)+1  
   KCKTA  = SHELLANGMOMA*(SHELLANGMOMA+1)/2  
   SPVAL= KCKTA*KCKTA
   JSTA = PRIEXPSTART(ISHELA)
   DO I=1,NGRID
      PA(1,i) = THCCOOR(1,i)-CENT(1,ISHELA)
      PA(2,i) = THCCOOR(2,i)-CENT(2,ISHELA)
      PA(3,i) = THCCOOR(3,i)-CENT(3,ISHELA)
      PA2(i) = PA(1,i)*PA(1,i)+PA(2,i)*PA(2,i)+PA(3,i)*PA(3,i)
   END DO
!   print*,'THC: SHELLANGMOMA',SHELLANGMOMA,KHKTA,KCKTA,SPVAL,JSTA,ISHELA
!   print*,'THC: THCCOOR(2,1)',THCCOOR(2,1),'CENT(2,ISHELA)',CENT(2,ISHELA)
!   print*,'THC: THCCOOR(2,10)',THCCOOR(2,10),'CENT(2,ISHELA)',CENT(2,ISHELA)
!   print*,'THC: THCCOOR(2,13)',THCCOOR(2,13),'CENT(2,ISHELA)',CENT(2,ISHELA)
!   print*,'THC: THCCOOR(2,19)',THCCOOR(2,19),'CENT(2,ISHELA)',CENT(2,ISHELA)
!   print*,'THC: PA',PA(2,:)
   CALL II_BLGETGAOTHC(LUPRI,NGRID,NBAST,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
        & SPHMAT(SPINDEX(SHELLANGMOMA):SPINDEX(SHELLANGMOMA)+KHKTA*KCKTA-1),&
        & GAO,IADR,PA,PA2,CC(CCINDEX(ISHELA)),&
        & CCSTART(CCINDEX(ISHELA)),PRIEXP,LVALUE(1:SPVAL,SHELLANGMOMA),&
        & MVALUE(1:SPVAL,SHELLANGMOMA),NVALUE(1:SPVAL,SHELLANGMOMA),SPVAL)
   IADR = IADR + KHKTA !UPDATE ACTIVEORBITALINDEX
END DO

!write (lupri,*) 'GAO(NGRID,NBAST)'
!call ls_output(gao,1,ngrid,1,NBAST,ngrid,NBAST,1,lupri)

!APPLY WEIGHT ACCORDING TO EQ. 29 of JCP 137, 224106
!To obtain collocation matrix X(NGRID,NBAST)  THCWG(:)
DO I=1,NGRID
   PA2(I) = THCWG(I)**0.25E0_realk
ENDDO
DO J=1,NBAST
   DO I=1,NGRID
      GAO(I,J) = GAO(I,J)*PA2(I)
   END DO
ENDDO
call mem_dealloc(PA)
call mem_dealloc(PA2)

END SUBROUTINE II_BLGETSOSTHC

!> \brief evaluates the pure GAO(gaussian atomic orbitals) , so no derivatives
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGAOTHC(LUPRI,NGRID,NBAST,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
                     & CSP,GAO,IADR,PA,PA2,CC,CCSTART,PRIEXP,LVALUE,MVALUE,&
                     & NVALUE,SPVAL)
IMPLICIT NONE
INTEGER,intent(in)        :: NGRID,NBAST,MXPRIM,KHKTA,JSTA,KCKTA,LUPRI,SHELLANGMOMA
INTEGER,intent(in)        :: SPVAL,CCSTART,IADR,LVALUE(SPVAL),MVALUE(SPVAL)
INTEGER,intent(in)        :: NVALUE(SPVAL)
REAL(REALK),PARAMETER     :: D0 = 0.0E0_realk, D1 = 1.0E0_realk
REAL(REALK),intent(inout) :: GAO(NGRID,NBAST)
REAL(REALK),intent(in)    :: PA(3,NGRID),PA2(NGRID)
REAL(REALK),intent(in)    :: CSP(KHKTA,KCKTA),PRIEXP(MXPRIM)
TYPE(LSMATRIX),intent(in) :: CC
!
REAL(REALK)    :: GAZ,GAXX,GAX,GAY,GAYY,GMAX,PA_1,PA_2,PA_3
REAL(REALK)    :: PRICCFVAL,PRIEXPVAL
INTEGER        :: I,J,TEMPI,LVALJ, MVALJ, NVALJ, LVALI, MVALI, NVALI,K
INTEGER        :: ISTART1,ISTART2,ISTART3,ISTART4,ISTART5,ISTART6,ISTART
REAL(REALK),pointer  :: GA(:),CINT(:),CINT2(:,:)

call mem_alloc(GA,NGRID)
!loop over primitives
ISTART=IADR
ISTART1=ISTART+1
ISTART2=ISTART+2
ISTART3=ISTART+3
ISTART4=ISTART+4
ISTART5=ISTART+5
ISTART6=ISTART+6
I=1
J = JSTA+CCSTART-1+I
PRICCFVAL = CC%elms(I)
PRIEXPVAL = -PRIEXP(J)
DO K = 1, NGRID
   GA(K) = PRICCFVAL*EXP(PRIEXPVAL*PA2(K))
END DO
DO I = 2,CC%nrow
   J = JSTA+CCSTART-1+I
   PRICCFVAL = CC%elms(I)
   PRIEXPVAL = -PRIEXP(J)
   DO K = 1, NGRID
      GA(K) = GA(K) + PRICCFVAL*EXP(PRIEXPVAL*PA2(K))
   END DO
ENDDO
!     
!     contracted orbitals
!
IF (SHELLANGMOMA .EQ. 1) THEN      !s orbitals
   DO K = 1, NGRID
      GAO(K,ISTART1) = GA(K)
   END DO
ELSEIF (SHELLANGMOMA .EQ. 2) THEN !p orbitals
   DO K = 1, NGRID
      GAO(K,ISTART1) = PA(1,K)*GA(K)
      GAO(K,ISTART2) = PA(2,K)*GA(K)
      GAO(K,ISTART3) = PA(3,K)*GA(K)
   END DO
ELSEIF (SHELLANGMOMA .EQ. 3) THEN !d orbitals
   DO K = 1, NGRID
      GAO(K,ISTART1) = PA(2,K)*PA(1,K)*GA(K)
      GAO(K,ISTART2) = PA(2,K)*PA(3,K)*GA(K)
      GAO(K,ISTART3) = -0.288675134594813E0_realk*(PA(1,K)*PA(1,K) &
           & + PA(2,K)*PA(2,K))*GA(K)&
           & + 0.577350269189626E0_realk*PA(3,K)*PA(3,K)*GA(K)
      GAO(K,ISTART4) = PA(1,K)*PA(3,K)*GA(K)
      GAO(K,ISTART5) = 0.5E0_realk*(PA(1,K)*PA(1,K) - PA(2,K)*PA(2,K))*GA(K)
   END DO
ELSEIF (SHELLANGMOMA .EQ. 4) THEN !f orbitals
   !KCKTA=10
   call mem_alloc(CINT2,NGRID,KCKTA)
   DO J = 1, KCKTA
      LVALJ = LVALUE(J)
      MVALJ = MVALUE(J)
      NVALJ = NVALUE(J)
      DO K = 1, NGRID
         CINT2(K,J) = (PA(1,K)**LVALJ)*(PA(2,K)**MVALJ)&
              &                 *(PA(3,K)**NVALJ)*GA(K)
      END DO
   ENDDO
   DO K = 1, NGRID
      GAO(K,ISTART+1) = &
           & + 0.612372435695794E0_realk*CINT2(K,2) &
           & - 0.204124145231932E0_realk*CINT2(K,7) 
   ENDDO
   DO K = 1, NGRID
      GAO(K,ISTART+2) = CINT2(K,5)
   ENDDO
   DO K = 1, NGRID
      GAO(K,ISTART+3) = &
           & -0.158113883008419E0_realk*(CINT2(K,2)+CINT2(K,7))&
           & +0.632455532033676E0_realk*CINT2(K,9)
   ENDDO
   DO K = 1, NGRID
      GAO(K,ISTART+4) = &
           & -0.387298334620742E0_realk*(CINT2(K,3)+CINT2(K,8))&
           & +0.258198889747161E0_realk*CINT2(K,10)
   ENDDO
   DO K = 1, NGRID
      GAO(K,ISTART+5) = &
           & -0.158113883008419E0_realk*(CINT2(K,1)+CINT2(K,4))&
           & +0.632455532033676E0_realk*CINT2(K,6)
   END DO
   DO K = 1, NGRID
      GAO(K,ISTART+6) = &
           & 0.5E0_realk*(CINT2(K,3)-CINT2(K,8))
   END DO
   DO K = 1, NGRID
      GAO(K,ISTART+7) = &
           & +0.204124145231932E0_realk*CINT2(K,1)&
           & -0.612372435695794E0_realk*CINT2(K,4)
   END DO
   call mem_dealloc(CINT2)
ELSE !higher than f orbitals
   call mem_alloc(CINT,NGRID)
   DO I = 1, KHKTA
      CALL LS_DZERO(GAO(1,ISTART+I),NGRID)
   END DO
   DO J = 1, KCKTA
      LVALJ = LVALUE(J)
      MVALJ = MVALUE(J)
      NVALJ = NVALUE(J)
      DO K = 1, NGRID
         CINT(K) = (PA(1,K)**LVALJ)*(PA(2,K)**MVALJ)&
              &                 *(PA(3,K)**NVALJ)*GA(K)
      END DO
      DO I = 1, KHKTA
         TEMPI=ISTART+I
         DO K = 1, NGRID
            GAO(K,TEMPI) = GAO(K,TEMPI) + CSP(I,J)*CINT(K)
         END DO
      END DO
   END DO
   call mem_dealloc(CINT)
END IF
call mem_dealloc(GA)

END SUBROUTINE II_BLGETGAOTHC

!> \brief build spherical transformation matrices
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE THC_Build_PRECALCULATED_SPHMAT(LUPRI,MAXANGMOM,SIZE,SPHMAT,SPINDEX) 
use math_fun
IMPLICIT NONE
INTEGER,intent(in)        :: MAXANGMOM,SIZE,LUPRI
INTEGER,intent(inout)     :: SPINDEX(MAXANGMOM+1)
REAL(REALK),intent(inout) :: SPHMAT(SIZE)
!
INTEGER          :: L,I
Real(realk), parameter :: DM1 = -1.0E0_realk, DO = 0.0E0_realk, D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER  :: M1,MADR,MABS,V0, NDER, IOFF,nAngmom
INTEGER  :: EE,FF,GG,BX,BY,BZ,II,JJ,KK,PX,PY,PZ
REAL(realk)  :: FACNRM,FAC1,FAC2, FAC3,FAC4, FACTOR
INTEGER  :: T,U,V,A,B,C,P,Q,R,X,Y,Z,TOT,IADR,AX0,AY0,AZ0,NSIZE
INTEGER  :: M,Ncol,Nrow,INDEX,NDIM,STARTINDEX

IF(MAXANGMOM .LE. 1)CALL LSQUIT('ERROR IN Build_PRECALCULATED_SPHMAT',lupri)
! CALL LS_DZERO(SPHMAT,SIZE) SHOULD BE DONE OUTSIDE

nANGMOM=MAXANGMOM+1 
NSIZE=0
STARTINDEX=1
SPINDEX(1)=1
SPINDEX(2)=1
DO I=3,nANGMOM
   SPINDEX(I)=STARTINDEX
   L = I-1 !angmom
   NRow = 2*L+1
   NCol = (L+1)*(L+2)/2
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
            TOT = X + Y + Z
            IADR = 1 + (2*L+1)*(NCRT(X,Y,Z)-1) + L + MADR+NSIZE
            SPHMAT(IADR) = SPHMAT(IADR) + FACTOR 
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
   NSIZE= NSIZE+Nrow*Ncol 
   STARTINDEX=STARTINDEX+Nrow*Ncol
ENDDO
END SUBROUTINE THC_BUILD_PRECALCULATED_SPHMAT

subroutine Get_THC_grid_overlap(S,X,nbasis,ngrid)
  implicit none
  integer,intent(in) :: nbasis,ngrid
  real(realk),intent(inout) :: S(ngrid,ngrid)
  real(realk),intent(in) :: X(ngrid,nbasis)
  !
  integer :: A,B,P,Q
  real(realk) :: tmpA,tmpB
  DO Q=1,ngrid
     DO P=1,ngrid
        tmpB = 0.0E0_realk
        DO B=1,nbasis
           tmpB = tmpB + X(P,B)*X(Q,B)
        ENDDO
        tmpA = 0.0E0_realk
        DO A=1,nbasis
           tmpA = tmpA + X(P,A)*X(Q,A)
        ENDDO
        S(P,Q) = tmpA*tmpB
     enddo
  enddo
end subroutine Get_THC_grid_overlap

subroutine Get_THC_grid_overlap2(S,XO,XV,nocc,nvirt,ngrid)
  implicit none
  integer,intent(in) :: nocc,nvirt,ngrid
  real(realk),intent(inout) :: S(ngrid,ngrid)
  real(realk),intent(in) :: XO(ngrid,nocc)
  real(realk),intent(in) :: XV(ngrid,nvirt)
  !
  integer :: A,B,P,Q
  real(realk) :: tmpA,tmpB
  DO Q=1,ngrid
     DO P=1,ngrid
        tmpB = 0.0E0_realk
        DO B=1,nocc
           tmpB = tmpB + XO(P,B)*XO(Q,B)
        ENDDO
        tmpA = 0.0E0_realk
        DO A=1,nvirt
           tmpA = tmpA + XV(P,A)*XV(Q,A)
        ENDDO
        S(P,Q) = tmpA*tmpB
     enddo
  enddo
end subroutine Get_THC_grid_overlap2

subroutine Get_THC_grid_overlap_inv(S,S_inv,ngrid,epsilon)
  implicit none
  integer,intent(in) :: ngrid
  real(realk),intent(inout) :: S_inv(ngrid,ngrid)
  real(realk),intent(in)    :: S(ngrid,ngrid)
  real(realk),intent(in)    :: epsilon
  !
  !Make an SVD look at S^-1/2 ? 
  CALL DCOPY(ngrid*ngrid,S,1,S_inv,1)
  CALL PSEUDOINVERSE(ngrid,ngrid,S_inv,epsilon)

end subroutine Get_THC_grid_overlap_inv

END MODULE THC_UTIL
