!> @file
!> Contains the main MBIE integral drivers

!> \brief Main MBIE drivers for the calculation of absolute spherical multipoles used in MBIE. (JCP 123,184101) 
!> \author T. Kjaergaard
!> \date 2011 
MODULE MBIEintegraldriver
  use integraldriver
  use TYPEDEF
  use ODbatches
  use OD_type
  use precision
  use lstiming
  use thermite_OD
  use Thermite_integrals
  SAVE
CONTAINS
!> \brief Main driver for the calculation of absolute spherical multipoles used in MBIE
!> \author T. Kjaergaard
!> \date 2011
!>
!>  This is the Main driver for the calculation of absolute spherical multipoles used in 
!>  Multipole Based Integral Estimate (MBIE) Screening.
!>  It sets op the overlap distribution (call ODbatch) from the Atomic orbitals given in 
!>  input INPUT%AO.
!>  
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param INTOUT the integral output specifications, determines how the output should be given
SUBROUTINE MBIE_INTEGRAL_DRIVER(LUPRI,IPRINT,INPUT,INTOUT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: INTOUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD
Real(realk)          :: TS,TE

CALL Create_ODbatches(OD,INPUT,'LHS',LUPRI)
IF(OD%nbatches .NE. 0)THEN
   CALL PRINT_OD(OD,LUPRI,IPRINT)
   CALL LSTIMER('START ',TS,TE,LUPRI)
   CALL MBIE_INT(OD,INPUT,INTOUT,LUPRI,IPRINT)
   CALL LSTIMER('MBIE IntDriver',TS,TE,LUPRI)
   CALL FREE_ODitem(OD)
ENDIF

END SUBROUTINE MBIE_INTEGRAL_DRIVER

!> \brief Absolute spherical multipole integral driver, used for MBIE
!> \author T. Kjaergaard
!> \date 2010
!>
!>  \param OD the ODbatch
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param OUTPUT the integral output specifications, determines how the output should be given
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE MBIE_INT(OD,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD
!
Integer               :: ILHS,nPrimP
TYPE(Overlap)         :: P
TYPE(Allocitem)       :: Allocations
TYPE(Integralitem)    :: Integral
TYPE(TUVitem),target  :: SharedTUV
integer :: nthreads,tid
#ifdef VAR_OMP
integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif

IF (IPRINT.GT. 5) THEN
  CALL LSHEADER(LUPRI,'MBIE_INT')
ENDIF
CALL initTUVitem(sharedTUV,Input,OD,OD,LUPRI,IPRINT)
Integral%TUV => sharedTUV
!CALL initAlloc(Allocations,LUPRI,IPRINT,'LHS')
call allocitem_zero(Allocations,'Both')
CALL SET_ALLOC(Allocations,Input,OD,'LHS',IPRINT,LUPRI)

IF(.NOT.INPUT%noOMP)call mem_TurnONThread_Memory()
!$OMP PARALLEL IF(.NOT.INPUT%noOMP) DEFAULT(none) PRIVATE(integral,&
!$OMP P,ILHS,nPrimP,tid,nthreads) SHARED(INPUT,Allocations,OD,lupri,iprint,output,sharedTUV)
IF(.NOT.INPUT%noOMP)call init_threadmemvar()

#ifdef VAR_OMP
nthreads=OMP_GET_NUM_THREADS()
tid = omp_get_thread_num()
#else
nthreads=1
tid=0
#endif

call INIT_BUFCOUNTERS(1)
call MEM_INIT_OVERLAP(Allocations,0,Input,1,IPRINT,LUPRI)
call ALLOC_ODLHS_BUFFERS
integral%TUV => sharedTUV 
CALL INIT_OVERLAP(P,Allocations,0,Input,1,IPRINT,LUPRI)

!!$OMP DO SCHEDULE(DYNAMIC,1)
!DO ILHS=1,OD%nbatches
DO ILHS=1+tid,OD%nbatches,nthreads
   CALL SET_Overlap(P,Input,SharedTUV,Integral,OD%BATCH(ILHS),1,LUPRI,IPRINT,.FALSE.)
   CALL BUILD_MBIE_MM(P,OUTPUT%ScreenTensor,lupri,iprint,ILHS,OD%BATCH(ILHS)%IA,OD%BATCH(ILHS)%IB,Input%sameLHSaos)
ENDDO
!!$OMP END DO 

CALL FREE_OVERLAP(P)
call DEALLOC_ODLHS_BUFFERS
IF(.NOT.INPUT%noOMP)call collect_thread_memory()
!$OMP END PARALLEL
IF(.NOT.INPUT%noOMP)call mem_TurnOffThread_Memory()

CALL freeTUVitem(sharedTUV,Input)

END SUBROUTINE MBIE_INT

SUBROUTINE BUILD_MBIE_MM(P,OUT,lupri,iprint,ILHS,IA,IB,sameLHSaos)
implicit none
logical :: sameLHSaos
type(overlap),intent(in)      :: P
TYPE(lstensor),intent(inout)  :: OUT
Integer,intent(in) :: LUPRI,iprint,ILHS,IA,IB
!
real(realk) :: GAMMA(0:15),JVAL,e1,e2,exponent,VAL,m,X,Y,Z,Kab
integer :: nAngA,nAngB,iAngA,iAngB,l1,l2,nP1,nP2,nC1,nC2,jmax,i12
integer :: i1,i2,ic1,ic2,J,J1,Xorder,nPrim,i,t,e,e12,ediff
integer :: AtomA,AtomB,batchA,BatchB,Gindex,MBIELEVEL
REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
REAL(REALK), PARAMETER :: PI = 3.1415926535897932384626E0_realk
REAL(REALK), PARAMETER :: D4 = 4E0_realk,D2=2E0_realk,D1=1E0_realk,D6=6E0_realk,D3=3E0_realk
real(realk),pointer :: DIST(:)!,pexp2(:)
real(realk),pointer :: TMP1(:,:),MX(:,:,:),ETIJ2(:,:,:,:,:)
real(realk),pointer :: TMP2(:,:),M000(:),ETIJ(:,:,:,:,:),ABSM(:,:,:)
Real(realk),pointer      :: CC(:,:),VAL2(:)
real(realk) :: MBIEMOM(0:1),FAC

MBIELEVEL = 1
! Definition of a helping array
! argument to gammafunction is ((angA+angB+n+3)/2) argument to GAMMA array(angA+angB+n)
! SO if argument to GAMMA is even then the argument to the function is odd abd we use the formula 
! Gamma(m+1/2)=1*3*5**(2m-1)*sqrtpi/(2**m)   (1)
! if odd argument to GAMMA, even argument to function and we can use the formula
! Gamma(n+1) = n!                            (2)
 GAMMA(0)=0.5E0_realk*sqrtpi            !m=1  in Eq. 1
 GAMMA(1)=1                       !n=1  in Eq. 2
 GAMMA(2)=0.75E0_realk*sqrtpi           !m=2  in Eq. 1
 GAMMA(3)=2                       !n=2  in Eq. 2 
 GAMMA(4)=1.875*sqrtpi            !m=3  in Eq. 1 
 GAMMA(5)=6                       !n=3  in Eq. 2 
 GAMMA(6)=6.5625*sqrtpi           !m=4  in Eq. 1 
 GAMMA(7)=24                      !n=4  in Eq. 2 
 GAMMA(8)=29.53125*sqrtpi         !m=5  in Eq. 1  
 GAMMA(9)=120                     !n=5  in Eq. 2 
 GAMMA(10)=162.421875*sqrtpi      !m=6  in Eq. 1  
 GAMMA(11)=720                    !n=6  in Eq. 2 
 GAMMA(12)=1055.7421875*sqrtpi    !m=7  in Eq. 1  
 GAMMA(13)=5040                   !n=7  in Eq. 2 
 GAMMA(14)=7918.06640625E0_realk*sqrtpi !m=8  in Eq. 1
 GAMMA(15)=40320                 

 nAngA = P%orbital1%nAngmom
 nAngB = P%orbital2%nAngmom
 MBIEMOM(0:1)=0E0_realk
 DO IangA=1,nAngA
  DO IangB=1,nAngB
   l1=P%orbital1%ANGMOM(IangA)
   l2=P%orbital2%ANGMOM(IangB)
   nPrim=P%nPrimitives
   JMAX = l1+l2
   i12=0
!   call mem_alloc(pexp2,nPrim)
   call mem_alloc(DIST,nPrim)
   DO i12=1,nPrim
!     pexp2(i12) = 1/(D2*P%exponents(i12))
     X = P%center(1+(i12-1)*3) - P%ODcenter(1)
     Y = P%center(2+(i12-1)*3) - P%ODcenter(2)
     Z = P%center(3+(i12-1)*3) - P%ODcenter(3)
     DIST(i12) = SQRT(X*X+Y*Y+Z*Z)
  ENDDO
  !Calculate the Ecoefficients
  call mem_alloc(ETIJ,nprim,l1+l2,l1,l2,3,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
  CALL GET_ECOEFF(ETIJ,nprim,l1+l2,l1,l2,P%nPrimitives,1,P,LUPRI,iprint)
  call mem_alloc(ETIJ2,nprim,l1+l2,l1,l2,1,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
  DO j=0,l2
   DO i=0,l1
    DO t=0,i+j
     DO i12=1,nprim
      ETIJ2(i12,t,i,j,1)=MAX(ABS(ETIJ(i12,t,i,j,1)),ABS(ETIJ(i12,t,i,j,2)),ABS(ETIJ(i12,t,i,j,3)))
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  call mem_dealloc(ETIJ)

  !Calculate the Coefficients comming from the multipole moments
  !At the moment we only do MBIE-1 as suggested by the paper (JCP 123,184101)
  !but this have been test up to 3
  call mem_alloc(MX,nprim,Mbielevel,Mbielevel,.FALSE.,.TRUE.,.TRUE.)
  CALL LS_DZERO(MX,nPrim*(Mbielevel+1)*(Mbielevel+1))
!  IF(Mbielevel .GE. 1)THEN
!   IF(Mbielevel .GT. 1)THEN
!    CALL STANDARDLOOPX(P%preExpFac,MX,DIST,Mbielevel,nPrim,Pexp2,Mbielevel,lupri)
!   ELSEIF(Mbielevel .EQ. 1)THEN
    CALL ORDER_EQ_ONE_LOOP(P%preExpFac,MX,DIST,Mbielevel,nPrim,Mbielevel,lupri)
!   ENDIF
!  ELSE
!   CALL DCOPY(nPrim,P%preExpFac,1,MX(:,0,0),1)
!  ENDIF
!  call mem_dealloc(pexp2)
  call mem_dealloc(DIST)

  !Calculate the actual integral \int | r^{n}_{P} \Omega(r_{P}) | r^{2}_{P} dr_{P} 
  call mem_alloc(ABSM,nprim,jmax,Mbielevel,.FALSE.,.TRUE.,.TRUE.)
  DO t=0,jmax
   DO e = 0,Mbielevel
    JVAL = e+t+3
    J1 = e+t
    FAC = D2*PI*GAMMA(J1)
    DO i12=1,nPrim
     ABSM(i12,t,e) = FAC/(P%exponents(i12)**(JVAL/D2))
    ENDDO
   ENDDO
  enddo  

  !Combine
  call mem_alloc(VAL2,nPrim)
  call mem_alloc(TMP1,nPrim,MBIELEVEL,.FALSE.,.TRUE.)
  DO e = 0,Mbielevel
     DO i12=1,nPrim
        VAL2(i12) = 0E0_realk
     ENDDO
     DO e12=0,e
        ediff=e-e12
        DO t=0,l1+l2
           DO i12=1,nPrim
              VAL2(i12) = VAL2(i12) + MX(i12,e,ediff)*ETIJ2(i12,t,l1,l2,1)*ABSM(i12,t,ediff)
           ENDDO
        ENDDO
     ENDDO
     DO i12=1,nPrim
        TMP1(i12,e)=VAL2(i12)
     ENDDO
  ENDDO
  call mem_dealloc(VAL2)
  call mem_dealloc(MX)
  call mem_dealloc(ETIJ2)
  call mem_dealloc(ABSM)

  nC1 = P%orbital1%nContracted(iAngA)
  nC2 = P%orbital2%nContracted(iAngB)
  call mem_alloc(TMP2,nC2*nC1,MBIELEVEL,.FALSE.,.TRUE.)
  call mem_alloc(CC,nC2*nC1,nPrim)
  CALL ConstructContraction_PA(CC,P,nPrim,nC1,nC2,iAngA,iAngB,LUPRI,IPRINT)
  DO J=0,MBIELEVEL   
     DO iC1 = 1,nC1*nC2
        VAL = CC(iC1,1)*TMP1(1,J)
        DO i12 = 2,nPrim
           VAL = VAL+CC(iC1,i12)*TMP1(i12,J)
        ENDDO
        TMP2(iC1,J) = VAL
     ENDDO
  ENDDO
  call mem_dealloc(CC)
  call mem_dealloc(TMP1)

  DO J=0,MBIELEVEL   
     DO iC1=1,nC1*nC2
        MBIEMOM(J)=MAX(MBIEMOM(J),ABS(TMP2(iC1,J)))
     ENDDO
  ENDDO
  call mem_dealloc(TMP2)
 ENDDO
ENDDO
IF(IPRINT.GT. 100)THEN
   WRITE(LUPRI,'(1X,A20)') 'Final MBIE (Absolute spherical moments)'
   WRITE(LUPRI,'(5F16.12/,(11X,5F16.12))') (MBIEMOM(J),J=0,MBIELEVEL)
ENDIF

IF(MBIELEVEL.GT.1)call lsquit('MBIE2 where to put it. TK',-1)

OUT%MBIE(1,IA,IB) = MBIEMOM(0)       !MBIE0
OUT%MBIE(2,IA,IB) = MBIEMOM(1)       !MBIE1
IF(sameLHSaos)THEN
   OUT%MBIE(1,IB,IA) = MBIEMOM(0)    !MBIE0
   OUT%MBIE(2,IB,IA) = MBIEMOM(1)    !MBIE1
ENDIF

END SUBROUTINE BUILD_MBIE_MM

END MODULE MBIEintegraldriver
