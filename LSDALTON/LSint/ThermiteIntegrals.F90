!> @file 
!> Contains routines to evaluate integrals of primitive hermite gaussian 
!> Thermite integral module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE Thermite_integrals
use typedef
use thermite_OD
use Integralparameters
use OverlapType
!****** INTEGRAND
TYPE Integrand
TYPE(OverlapPointer) :: P,Q
Integer              :: nAngmom
Integer              :: nPrimitives
Integer              :: maxContracted
Integer              :: minAngmom
Integer              :: maxAngmom
Integer              :: startAngmom
Integer              :: endAngmom
!Dimensions according to nPrimitives
!Integer,pointer      :: iprimP(:)
!Integer,pointer      :: iprimQ(:)
Real(realk),pointer  :: distance(:,:)
Real(realk),pointer  :: squaredDistance(:)
Real(realk),pointer  :: exponents(:)
Real(realk),pointer  :: reducedExponents(:)
Real(realk),pointer  :: integralPrefactor(:)
!One dimension turning into six: nOrbA,nOrbB,nOrbC,nOrbD,nPassP,nPassQ)
Real(realk),pointer  :: theta(:) 
!One dimension turning into seven: nOrbA,nOrbB,nOrbC,nOrbD,nPassP,nPassQ,iDeriv)
Real(realk),pointer  :: integral(:) 
Real(realk)          :: MOM_CENTER(3)
LOGICAL              :: samePQ
LOGICAL              :: kinetic
LOGICAL              :: reverseOrder
INTEGER              :: Operator
END TYPE Integrand

CONTAINS
!> \brief wrapper to recurrence relation for hermite integrals Eq. (9.9.18-9.9.20) in the book \f[  R^{n}_{t+1,u,v}(p,\textbf{R}_{PC}) = t R^{n+1}_{t-1,u,v}(p,\textbf{R}_{PC}) + X_{pC} R^{n+1}_{t,u,v}(p,\textbf{R}_{PC})  \f]
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param WTUV the output of the recurrence 
!> \param WJ000 Eq. (9.9.14) in the book \f$ R^{n}_{000}(p,\textbf{R}_{PC}) = (-2p)^{n}F_{n}(pR_{PC}^{2}) \f$
!> \param sharedTUV a TUVitem which contains tabulated boys function and TUVindexing
!> \param Rpq distance between overlap distribution P and Q (X,Y and Z distance) 
!> \param jmin minimum 
!> \param jmax maximum angular momentum ( N where t+u+v =< N )
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param ntuv number of TUV components
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE wtuvRecurrence(WTUV,WJ000,sharedTUV,Rpq,jmin,jmax,&
     &                    nPrim,ntuv,lupri,iprint)
  use memory_handling
  implicit none
  Integer,intent(in)           :: jmin,jmax,nPrim,ntuv,lupri,iprint
  TYPE(TUVitem),intent(in)     :: SharedTUV
  Real(realk),intent(in)       :: WJ000(0:jmax,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,ntuv)
  !
  !INTEGER,PARAMETER :: MAXJ = 50
  Real(realk) :: Xtemp,Ytemp,Ztemp,WTUV2(nPrim,ntuv)
  Integer     :: n,j,t,u,v,ntuvfull,ituv,dir
  Real(realk),pointer :: CUR(:,:),OLD(:,:),TEMP(:,:)
  Real(realk) :: W100,W010,W001,W200,W110,W101,W020,W011,W002,WJ,WJ1,WJ2,WJ3
  integer :: tm1,um1,vm1,tm2,um2,vm2,m1,ituvm,ituvm2,ioffp
  integer :: idp1,idir,k
  IF(jmin.EQ.0.AND.jmax.LT.9)then
     IF (jmax.EQ. 0) THEN
        call wtuvRecurrenceJMIN0JMAX0(WTUV,WJ000,nPrim)
     ELSEIF (jmax.EQ. 1) THEN
        call wtuvRecurrenceJMIN0JMAX1(WTUV,WJ000,Rpq,nPrim)
     ELSEIF (jmax.EQ. 2) THEN
        call wtuvRecurrenceJMIN0JMAX2(WTUV,WJ000,Rpq,nPrim)
     ELSEIF (jmax.EQ. 3) THEN
        call wtuvRecurrenceJMIN0JMAX3(WTUV,WJ000,Rpq,nPrim)
     ELSEIF (jmax.EQ. 4) THEN
        call wtuvRecurrenceJMIN0JMAX4(WTUV,WJ000,Rpq,nPrim)
     ELSEIF (jmax.EQ. 5) THEN
        call wtuvRecurrenceJMIN0JMAX5(WTUV,WJ000,Rpq,nPrim)
     ELSEIF (jmax.EQ. 6) THEN
        call wtuvRecurrenceJMIN0JMAX6(WTUV,WJ000,Rpq,nPrim)
     ELSEIF (jmax.EQ. 7) THEN
        call wtuvRecurrenceJMIN0JMAX7(WTUV,WJ000,Rpq,nPrim)
     ELSEIF (jmax.EQ. 8) THEN
        call wtuvRecurrenceJMIN0JMAX8(WTUV,WJ000,Rpq,nPrim)         
     ENDIF
  ELSE
     IF (jmax.EQ. 0) THEN
        CALL DCOPY(nPrim,WJ000,1,WTUV,1)
     ELSE IF (jmax.EQ. 1) THEN
      IF (jmin.EQ. 1) THEN
       DO idir=1,3
        DO n=1,nPrim
         WTUV(n,idir) = Rpq(n,idir)*WJ000(1,n)
        ENDDO
       ENDDO
      ELSE
       DO n=1,nPrim
        WTUV(n,1) = WJ000(0,n)
       ENDDO
       DO idir=1,3
        idp1 = idir+1
        DO n=1,nPrim
         WTUV(n,idp1) = Rpq(n,idir)*WJ000(1,n)
        ENDDO
       ENDDO
      ENDIF
     ELSE IF (jmax.EQ. 2) THEN
      IF (jmin.EQ. 2) THEN
       DO n=1,nPrim
        Xtemp = Rpq(n,1)
        Ytemp = Rpq(n,2)
        Ztemp = Rpq(n,3)
        WJ1 = WJ000(1,n)
        WJ2 = WJ000(2,n)
        WTUV(n,1) = WJ1 + Xtemp*Xtemp*WJ2
        WTUV(n,2) =       Xtemp*Ytemp*WJ2
        WTUV(n,3) =       Xtemp*Ztemp*WJ2
        WTUV(n,4) = WJ1 + Ytemp*Ytemp*WJ2
        WTUV(n,5) =       Ytemp*Ztemp*WJ2
        WTUV(n,6) = WJ1 + Ztemp*Ztemp*WJ2
       ENDDO
      ELSE IF (jmin.EQ. 1) THEN
       DO n=1,nPrim
        Xtemp = Rpq(n,1)
        Ytemp = Rpq(n,2)
        Ztemp = Rpq(n,3)
        WJ1 = WJ000(1,n)
        WJ2 = WJ000(2,n)
        WTUV(n,1) = Xtemp*WJ1
        WTUV(n,2) = Ytemp*WJ1
        WTUV(n,3) = Ztemp*WJ1
        WTUV(n,4) = WJ1 + Xtemp*Xtemp*WJ2
        WTUV(n,5) =       Xtemp*Ytemp*WJ2
        WTUV(n,6) =       Xtemp*Ztemp*WJ2
        WTUV(n,7) = WJ1 + Ytemp*Ytemp*WJ2
        WTUV(n,8) =       Ytemp*Ztemp*WJ2
        WTUV(n,9) = WJ1 + Ztemp*Ztemp*WJ2
       ENDDO
      ELSE
       DO n=1,nPrim
        Xtemp = Rpq(n,1)
        Ytemp = Rpq(n,2)
        Ztemp = Rpq(n,3)
        WJ  = WJ000(0,n)
        WJ1 = WJ000(1,n)
        WJ2 = WJ000(2,n)
        WTUV(n,1)  = WJ
        WTUV(n,2)  = Xtemp*WJ1
        WTUV(n,3)  = Ytemp*WJ1
        WTUV(n,4)  = Ztemp*WJ1
        WTUV(n,5)  = WJ1 + Xtemp*Xtemp*WJ2
        WTUV(n,6)  =       Xtemp*Ytemp*WJ2
        WTUV(n,7)  =       Xtemp*Ztemp*WJ2
        WTUV(n,8)  = WJ1 + Ytemp*Ytemp*WJ2
        WTUV(n,9)  =       Ytemp*Ztemp*WJ2
        WTUV(n,10) = WJ1 + Ztemp*Ztemp*WJ2
       ENDDO
      ENDIF
     ELSE  ! J > 2
      ntuvfull = (jmax+1)*(jmax+2)*(jmax+3)/6
      ioffp    = jmin*(jmin+1)*(jmin+2)/6+1
      call mem_alloc(CUR,nPrim,ntuvfull)
      call mem_alloc(OLD,nPrim,ntuvfull)
      
      DO j=jmax-3,0,-1
       TEMP => CUR
       CUR  => OLD
       OLD  => TEMP
       DO n=1,nPrim
        Xtemp = Rpq(n,1)
        Ytemp = Rpq(n,2)
        Ztemp = Rpq(n,3)
        WJ   = WJ000(j,n)
        WJ1  = WJ000(j+1,n)
        WJ2  = WJ000(j+2,n)
        WJ3  = WJ000(j+3,n)
        W100 = Xtemp*WJ2
        W010 = Ytemp*WJ2
        W001 = Ztemp*WJ2
        W200 = WJ2 + Xtemp*Xtemp*WJ3
        W110 =       Xtemp*Ytemp*WJ3
        W101 =       Xtemp*Ztemp*WJ3
        W020 = WJ2 + Ytemp*Ytemp*WJ3
        W011 =       Ytemp*Ztemp*WJ3
        W002 = WJ2 + Ztemp*Ztemp*WJ3
        CUR(n,1)  = WJ                           !000
        CUR(n,2)  = Xtemp*WJ1                    !100
        CUR(n,3)  = Ytemp*WJ1                    !010
        CUR(n,4)  = Ztemp*WJ1                    !001
        CUR(n,5)  = WJ1 + Xtemp*W100             !200
        CUR(n,6)  =       Ytemp*W100             !110
        CUR(n,7)  =       Ztemp*W100             !101
        CUR(n,8)  = WJ1 + Ytemp*W010             !020
        CUR(n,9)  =       Ztemp*W010             !011
        CUR(n,10) = WJ1 + Ztemp*W001             !002
        CUR(n,11) = 2.0E0_realk*W100 + Xtemp*W200      !300
        CUR(n,12) =       W010 + Xtemp*W110      !210
        CUR(n,13) =       W001 + Xtemp*W101      !201
        CUR(n,14) =       W100 + Ytemp*W110      !120
        CUR(n,15) =              Xtemp*W011      !111
        CUR(n,16) =       W100 + Ztemp*W101      !102
        CUR(n,17) = 2.0E0_realk*W010 + Ytemp*W020      !030
        CUR(n,18) =       W001 + Ytemp*W011      !021
        CUR(n,19) =       W010 + Ztemp*W011      !012
        CUR(n,20) = 2.0E0_realk*W001 + Ztemp*W002      !003
       ENDDO
       ituv = 20
       DO k=4,jmax-j
        DO t=k,0,-1
         DO u=k-t,0,-1
          v=k-t-u
          IF ((t.ge.u).AND.(t.ge.v)) THEN
           tm1 = t-1
           um1 = u
           vm1 = v
           tm2 = t-2
           um2 = u
           vm2 = v
           m1  = tm1
           dir = 1
          ELSEIF (u.ge.v) THEN
           tm1 = t
           um1 = u-1
           vm1 = v
           tm2 = t
           um2 = u-2
           vm2 = v
           m1  = um1
           dir = 2
          ELSE
           tm1 = t
           um1 = u
           vm1 = v-1
           tm2 = t
           um2 = u
           vm2 = v-2
           m1  = vm1
           dir = 3
          ENDIF
          ituv   = ituv + 1
          ituvm  = sharedTUV%tuvIndex(tm1,um1,vm1)
          ituvm2 = sharedTUV%tuvIndex(tm2,um2,vm2)
          DO n=1,nPrim
           CUR(n,ituv) = m1*OLD(n,ituvm2) + Rpq(n,dir)*OLD(n,ituvm)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      do t=1,ntuv
         DO n=1,nPrim
            WTUV(n,t)=CUR(n,ioffp+t-1)
         enddo
      ENDDO
!      CALL DCOPY(ntuv*nPrim,CUR(1:,ioffp:),1,WTUV,1)
      call mem_dealloc(CUR)
      call mem_dealloc(OLD)
   ENDIF
ENDIF
END SUBROUTINE wtuvRecurrence

SUBROUTINE wtuvRecurrenceJMIN0JMAX0(WTUV,WJ000,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:0,nPrim)
  Real(realk),intent(inout)    :: WTUV(nPrim,1)
  !
  integer :: n
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
END SUBROUTINE WTUVRECURRENCEJMIN0JMAX0

SUBROUTINE wtuvRecurrenceJMIN0JMAX1(WTUV,WJ000,Rpq,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:1,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,4)
  !
  integer :: n
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
  DO n=1,nPrim
     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
  ENDDO
END SUBROUTINE WTUVRECURRENCEJMIN0JMAX1

SUBROUTINE wtuvRecurrenceJMIN0JMAX2(WTUV,WJ000,Rpq,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:2,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,10)
  !
  integer :: n
  !000
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*Rpq(n,1)*WJ000(2,n)
     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*Rpq(n,2)*WJ000(2,n)
     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*Rpq(n,3)*WJ000(2,n)
     WTUV(n,6) = Rpq(n,1)*Rpq(n,2)*WJ000(2,n)
     WTUV(n,7) = Rpq(n,1)*Rpq(n,3)*WJ000(2,n)
     WTUV(n,9) = Rpq(n,2)*Rpq(n,3)*WJ000(2,n)
  ENDDO
end SUBROUTINE wtuvRecurrenceJMIN0JMAX2

SUBROUTINE wtuvRecurrenceJMIN0JMAX3(WTUV,WJ000,Rpq,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:3,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,20)
  !
  Real(realk) :: WTUVTEMP(2:4,2)
  Real(realk) :: WTUVTEMP25,WTUVTEMP28
  Real(realk) :: WTUVTEMP210,WTUVTEMP29
  integer :: n
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(n,2)     = Rpq(n,1)*WJ000(1,n)
     WTUVTEMP(2,1) = Rpq(n,1)*WJ000(2,n)
     WTUVTEMP(2,2) = Rpq(n,1)*WJ000(3,n)
     WTUV(n,3)     = Rpq(n,2)*WJ000(1,n)
     WTUVTEMP(3,1) = Rpq(n,2)*WJ000(2,n)
     WTUVTEMP(3,2) = Rpq(n,2)*WJ000(3,n)
     WTUV(n,4)     = Rpq(n,3)*WJ000(1,n)
     WTUVTEMP(4,1) = Rpq(n,3)*WJ000(2,n)
     WTUVTEMP(4,2) = Rpq(n,3)*WJ000(3,n)
     WTUV(n,6) = Rpq(n,1)*WTUVTEMP(3,1)
     WTUV(n,7) = Rpq(n,1)*WTUVTEMP(4,1)
     WTUV(n,9)  = Rpq(n,2)*WTUVTEMP(4,1)
     WTUVTEMP29 = Rpq(n,2)*WTUVTEMP(4,2)
     WTUV(n,5)  = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
     WTUV(n,8)  = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
     WTUV(n,10) = WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
     WTUVTEMP25 = WJ000(2,n) + Rpq(n,1)*WTUVTEMP(2,2)
     WTUVTEMP28 = WJ000(2,n) + Rpq(n,2)*WTUVTEMP(3,2)
     WTUVTEMP210 = WJ000(2,n) + Rpq(n,3)*WTUVTEMP(4,2)
     WTUV(n,14)= Rpq(n,1)*WTUVTEMP28
     WTUV(n,15)= Rpq(n,1)*WTUVTEMP29
     WTUV(n,16)= Rpq(n,1)*WTUVTEMP210
     WTUV(n,12)= Rpq(n,2)*WTUVTEMP25
     WTUV(n,19)= Rpq(n,2)*WTUVTEMP210
     WTUV(n,13)= Rpq(n,3)*WTUVTEMP25
     WTUV(n,18)= Rpq(n,3)*WTUVTEMP28
     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP25
     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP28
     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP210
  ENDDO
END SUBROUTINE WTUVRECURRENCEJMIN0JMAX3

SUBROUTINE wtuvRecurrenceJMIN0JMAX4(WTUV,WJ000,Rpq,nPrim)
  use memory_handling
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:4,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,35)
  !
  Integer     :: n,M
  Real(realk) :: WTUVTEMP(2:4,3)
  Real(realk) :: WTUVTEMP2(5:10,2)
  Real(realk) :: WTUVTEMP3(11:20)
  !000
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
     WTUVTEMP(2,1) = Rpq(n,1)*WJ000(2,n)
     WTUVTEMP(3,1) = Rpq(n,2)*WJ000(2,n)
     WTUVTEMP(4,1) = Rpq(n,3)*WJ000(2,n)     
     WTUVTEMP(2,2) = Rpq(n,1)*WJ000(3,n)
     WTUVTEMP(3,2) = Rpq(n,2)*WJ000(3,n)
     WTUVTEMP(4,2) = Rpq(n,3)*WJ000(3,n)     
     WTUVTEMP(2,3) = Rpq(n,1)*WJ000(4,n)
     WTUVTEMP(3,3) = Rpq(n,2)*WJ000(4,n)
     WTUVTEMP(4,3) = Rpq(n,3)*WJ000(4,n)
     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
     WTUVTEMP2(5,1) = WJ000(2,n) + Rpq(n,1)*WTUVTEMP(2,2)
     WTUVTEMP2(6,1) =              Rpq(n,1)*WTUVTEMP(3,2)
     WTUVTEMP2(7,1) =              Rpq(n,1)*WTUVTEMP(4,2)
     WTUVTEMP2(8,1) = WJ000(2,n) + Rpq(n,2)*WTUVTEMP(3,2)
     WTUVTEMP2(9,1) =              Rpq(n,2)*WTUVTEMP(4,2)
     WTUVTEMP2(10,1)= WJ000(2,n) + Rpq(n,3)*WTUVTEMP(4,2)     
     WTUVTEMP2(5,2) = WJ000(3,n) + Rpq(n,1)*WTUVTEMP(2,3)
!     WTUVTEMP2(6,2) =              Rpq(n,1)*WTUVTEMP(3,3)
!     WTUVTEMP2(7,2) =              Rpq(n,1)*WTUVTEMP(4,3)
     WTUVTEMP2(8,2) = WJ000(3,n) + Rpq(n,2)*WTUVTEMP(3,3)
     WTUVTEMP2(9,2) =              Rpq(n,2)*WTUVTEMP(4,3)
     WTUVTEMP2(10,2)= WJ000(3,n) + Rpq(n,3)*WTUVTEMP(4,3)
     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
     WTUVTEMP3(11)= 2*WTUVTEMP(2,2) + Rpq(n,1)*WTUVTEMP2(5,2)
    ! WTUVTEMP3(12)=  Rpq(n,2)*WTUVTEMP2(5,2)
     WTUVTEMP3(13)= Rpq(n,3)*WTUVTEMP2(5,2)
     WTUVTEMP3(14)= Rpq(n,1)*WTUVTEMP2(8,2)
    ! WTUVTEMP3(15)= Rpq(n,1)*WTUVTEMP2(9,2)
     WTUVTEMP3(16)= Rpq(n,1)*WTUVTEMP2(10,2)
     WTUVTEMP3(17)= 2*WTUVTEMP(3,2) + Rpq(n,2)*WTUVTEMP2(8,2)
     WTUVTEMP3(18)= Rpq(n,3)*WTUVTEMP2(8,2)
     WTUVTEMP3(19)= Rpq(n,2)*WTUVTEMP2(10,2)
     WTUVTEMP3(20)= 2*WTUVTEMP(4,2) + Rpq(n,3)*WTUVTEMP2(10,2)
     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11)
     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11)
     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11)
     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14)
     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13)
     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16)
     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17)
     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18)
     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19)
     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20)
     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17) 
     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17)
     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19)
     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20)
     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20)
  ENDDO
END SUBROUTINE wtuvRecurrenceJMIN0JMAX4

SUBROUTINE wtuvRecurrenceJMIN0JMAX5(WTUV,WJ000,Rpq,nPrim)
  use memory_handling
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:5,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,56)
  !
  Integer     :: n,M
  Real(realk) :: WTUVTEMP(2:4,4)
  Real(realk) :: WTUVTEMP2(5:10,3)
  Real(realk) :: WTUVTEMP3(11:20,2)
  Real(realk) :: WTUVTEMP4(21:35)
  !000
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
     DO M=1,4
        WTUVTEMP(2,M) = Rpq(n,1)*WJ000(M+1,n)
        WTUVTEMP(3,M) = Rpq(n,2)*WJ000(M+1,n)
        WTUVTEMP(4,M) = Rpq(n,3)*WJ000(M+1,n)     
     ENDDO
     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
     DO M=1,3
        WTUVTEMP2(5,M) = WJ000(M+1,n) + Rpq(n,1)*WTUVTEMP(2,M+1)
!        WTUVTEMP2(6,M) =                Rpq(n,1)*WTUVTEMP(3,M+1)
!        WTUVTEMP2(7,M) =                Rpq(n,1)*WTUVTEMP(4,M+1)
        WTUVTEMP2(8,M) = WJ000(M+1,n) + Rpq(n,2)*WTUVTEMP(3,M+1)
        WTUVTEMP2(9,M) =                Rpq(n,2)*WTUVTEMP(4,M+1)
        WTUVTEMP2(10,M)= WJ000(M+1,n) + Rpq(n,3)*WTUVTEMP(4,M+1)     
     ENDDO
     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
     DO M=1,2
        WTUVTEMP3(11,M)= 2*WTUVTEMP(2,M+1) + Rpq(n,1)*WTUVTEMP2(5,M+1)
!        WTUVTEMP3(12,M)=  Rpq(n,2)*WTUVTEMP2(5,M+1)
        WTUVTEMP3(13,M)= Rpq(n,3)*WTUVTEMP2(5,M+1)
        WTUVTEMP3(14,M)= Rpq(n,1)*WTUVTEMP2(8,M+1)
!        WTUVTEMP3(15,M)= Rpq(n,1)*WTUVTEMP2(9,M+1)
        WTUVTEMP3(16,M)= Rpq(n,1)*WTUVTEMP2(10,M+1)
        WTUVTEMP3(17,M)= 2*WTUVTEMP(3,M+1) + Rpq(n,2)*WTUVTEMP2(8,M+1)
        WTUVTEMP3(18,M)= Rpq(n,3)*WTUVTEMP2(8,M+1)
        WTUVTEMP3(19,M)= Rpq(n,2)*WTUVTEMP2(10,M+1)
        WTUVTEMP3(20,M)= 2*WTUVTEMP(4,M+1) + Rpq(n,3)*WTUVTEMP2(10,M+1)
     ENDDO

     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11,1)
     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11,1)
     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11,1)
     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14,1)
     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13,1)
     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16,1)
     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17,1)
     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18,1)
     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19,1)
     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20,1)
     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17,1) 
     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17,1)
     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19,1)
     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20,1)
     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20,1)

     WTUVTEMP4(21)= 3*WTUVTEMP2(5,2) + Rpq(n,1)*WTUVTEMP3(11,2)
!     WTUVTEMP4(22)= Rpq(n,2)*WTUVTEMP3(11,2)
     WTUVTEMP4(23)= Rpq(n,3)*WTUVTEMP3(11,2)
     WTUVTEMP4(24)= WTUVTEMP2(8,2) + Rpq(n,1)*WTUVTEMP3(14,2)
!     WTUVTEMP4(25)= Rpq(n,2)*WTUVTEMP3(13,2)
     WTUVTEMP4(26)= WTUVTEMP2(10,2) + Rpq(n,1)*WTUVTEMP3(16,2)
     WTUVTEMP4(27)= Rpq(n,1)*WTUVTEMP3(17,2)
!     WTUVTEMP4(28)= Rpq(n,1)*WTUVTEMP3(18,2)
!     WTUVTEMP4(29)= Rpq(n,1)*WTUVTEMP3(19,2)
     WTUVTEMP4(30)= Rpq(n,1)*WTUVTEMP3(20,2)
     WTUVTEMP4(31)= 3*WTUVTEMP2(8,2) + Rpq(n,2)*WTUVTEMP3(17,2) 
     WTUVTEMP4(32)= Rpq(n,3)*WTUVTEMP3(17,2)
     WTUVTEMP4(33)= WTUVTEMP2(10,2) + Rpq(n,2)*WTUVTEMP3(19,2)
     WTUVTEMP4(34)= Rpq(n,2)*WTUVTEMP3(20,2)
     WTUVTEMP4(35)= 3*WTUVTEMP2(10,2) + Rpq(n,3)*WTUVTEMP3(20,2)

     WTUV(n,  36) = 4*WTUVTEMP3(11,1) + Rpq(n,1)*WTUVTEMP4(21)
     WTUV(n,  37) = Rpq(n,2)*WTUVTEMP4(21)
     WTUV(n,  38) = Rpq(n,3)*WTUVTEMP4(21)
     WTUV(n,  39) = 2*WTUVTEMP3(14,1) + Rpq(n,1)*WTUVTEMP4(24)
     WTUV(n,  40) = Rpq(n,2)*WTUVTEMP4(23)
     WTUV(n,  41) = 2*WTUVTEMP3(16,1) + Rpq(n,1)*WTUVTEMP4(26)
     WTUV(n,  42) = WTUVTEMP3(17,1) + Rpq(n,1)*WTUVTEMP4(27)
     WTUV(n,  43) = Rpq(n,3)*WTUVTEMP4(24)
     WTUV(n,  44) = Rpq(n,2)*WTUVTEMP4(26)
     WTUV(n,  45) = WTUVTEMP3(20,1) + Rpq(n,1)*WTUVTEMP4(30)
     WTUV(n,  46) = Rpq(n,1)*WTUVTEMP4(31)
     WTUV(n,  47) = Rpq(n,1)*WTUVTEMP4(32)
     WTUV(n,  48) = Rpq(n,1)*WTUVTEMP4(33)
     WTUV(n,  49) = Rpq(n,1)*WTUVTEMP4(34)
     WTUV(n,  50) = Rpq(n,1)*WTUVTEMP4(35)
     WTUV(n,  51) = 4*WTUVTEMP3(17,1) + Rpq(n,2)*WTUVTEMP4(31)
     WTUV(n,  52) = Rpq(n,3)*WTUVTEMP4(31)
     WTUV(n,  53) = 2*WTUVTEMP3(19,1) + Rpq(n,2)*WTUVTEMP4(33)
     WTUV(n,  54) = WTUVTEMP3(20,1) + Rpq(n,2)*WTUVTEMP4(34)
     WTUV(n,  55) = Rpq(n,2)*WTUVTEMP4(35)
     WTUV(n,  56) = 4*WTUVTEMP3(20,1) + Rpq(n,3)*WTUVTEMP4(35)
  ENDDO
END SUBROUTINE wtuvRecurrenceJMIN0JMAX5

SUBROUTINE wtuvRecurrenceJMIN0JMAX6(WTUV,WJ000,Rpq,nPrim)
  use memory_handling
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:6,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,84)
  !
  Integer     :: n,M
  Real(realk) :: WTUVTEMP(2:4,5)
  Real(realk) :: WTUVTEMP2(5:10,4)
  Real(realk) :: WTUVTEMP3(11:20,3)
  Real(realk) :: WTUVTEMP4(21:35,2)
  Real(realk) :: WTUVTEMP5(36:56,1)
  !000
  !ERROR IN THIS SUBROUTINE SOMEWHERE
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
     DO M=1,5
        WTUVTEMP(2,M) = Rpq(n,1)*WJ000(M+1,n)
        WTUVTEMP(3,M) = Rpq(n,2)*WJ000(M+1,n)
        WTUVTEMP(4,M) = Rpq(n,3)*WJ000(M+1,n)     
     ENDDO
     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
     DO M=1,4
        WTUVTEMP2(5,M) = WJ000(M+1,n) + Rpq(n,1)*WTUVTEMP(2,M+1)
!        WTUVTEMP2(6,M) =                Rpq(n,1)*WTUVTEMP(3,M+1)
!        WTUVTEMP2(7,M) =                Rpq(n,1)*WTUVTEMP(4,M+1)
        WTUVTEMP2(8,M) = WJ000(M+1,n) + Rpq(n,2)*WTUVTEMP(3,M+1)
        WTUVTEMP2(9,M) =                Rpq(n,2)*WTUVTEMP(4,M+1)
        WTUVTEMP2(10,M)= WJ000(M+1,n) + Rpq(n,3)*WTUVTEMP(4,M+1)     
     ENDDO
     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
     DO M=1,3
        WTUVTEMP3(11,M)= 2*WTUVTEMP(2,M+1) + Rpq(n,1)*WTUVTEMP2(5,M+1)
!        WTUVTEMP3(12,M)=  Rpq(n,2)*WTUVTEMP2(5,M+1)
        WTUVTEMP3(13,M)= Rpq(n,3)*WTUVTEMP2(5,M+1)
        WTUVTEMP3(14,M)= Rpq(n,1)*WTUVTEMP2(8,M+1)
!        WTUVTEMP3(15,M)= Rpq(n,1)*WTUVTEMP2(9,M+1)
        WTUVTEMP3(16,M)= Rpq(n,1)*WTUVTEMP2(10,M+1)
        WTUVTEMP3(17,M)= 2*WTUVTEMP(3,M+1) + Rpq(n,2)*WTUVTEMP2(8,M+1)
        WTUVTEMP3(18,M)= Rpq(n,3)*WTUVTEMP2(8,M+1)
        WTUVTEMP3(19,M)= Rpq(n,2)*WTUVTEMP2(10,M+1)
        WTUVTEMP3(20,M)= 2*WTUVTEMP(4,M+1) + Rpq(n,3)*WTUVTEMP2(10,M+1)
     ENDDO

     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11,1)
     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11,1)
     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11,1)
     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14,1)
     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13,1)
     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16,1)
     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17,1)
     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18,1)
     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19,1)
     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20,1)
     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17,1) 
     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17,1)
     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19,1)
     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20,1)
     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20,1)

     DO M=1,2
        WTUVTEMP4(21,M)= 3*WTUVTEMP2(5,M+1) + Rpq(n,1)*WTUVTEMP3(11,M+1)
!        WTUVTEMP4(22,M)= Rpq(n,2)*WTUVTEMP3(11,M+1) !not used
        WTUVTEMP4(23,M)= Rpq(n,3)*WTUVTEMP3(11,M+1)
        WTUVTEMP4(24,M)= WTUVTEMP2(8,M+1) + Rpq(n,1)*WTUVTEMP3(14,M+1)
!        WTUVTEMP4(25,M)= Rpq(n,2)*WTUVTEMP3(13,M+1) !not used
        WTUVTEMP4(26,M)= WTUVTEMP2(10,M+1) + Rpq(n,1)*WTUVTEMP3(16,M+1)
        WTUVTEMP4(27,M)= Rpq(n,1)*WTUVTEMP3(17,M+1)
!        WTUVTEMP4(28,M)= Rpq(n,1)*WTUVTEMP3(18,M+1) !not used
!        WTUVTEMP4(29,M)= Rpq(n,1)*WTUVTEMP3(19,M+1) !not used
        WTUVTEMP4(30,M)= Rpq(n,1)*WTUVTEMP3(20,M+1) 
        WTUVTEMP4(31,M)= 3*WTUVTEMP2(8,M+1) + Rpq(n,2)*WTUVTEMP3(17,M+1) 
        WTUVTEMP4(32,M)= Rpq(n,3)*WTUVTEMP3(17,M+1)
        WTUVTEMP4(33,M)= WTUVTEMP2(10,M+1) + Rpq(n,2)*WTUVTEMP3(19,M+1)
        WTUVTEMP4(34,M)= Rpq(n,2)*WTUVTEMP3(20,M+1)
        WTUVTEMP4(35,M)= 3*WTUVTEMP2(10,M+1) + Rpq(n,3)*WTUVTEMP3(20,M+1)
     ENDDO

     WTUV(n,  36) = 4*WTUVTEMP3(11,1) + Rpq(n,1)*WTUVTEMP4(21,1)
     WTUV(n,  37) = Rpq(n,2)*WTUVTEMP4(21,1)
     WTUV(n,  38) = Rpq(n,3)*WTUVTEMP4(21,1)
     WTUV(n,  39) = 2*WTUVTEMP3(14,1) + Rpq(n,1)*WTUVTEMP4(24,1)
     WTUV(n,  40) = Rpq(n,2)*WTUVTEMP4(23,1)
     WTUV(n,  41) = 2*WTUVTEMP3(16,1) + Rpq(n,1)*WTUVTEMP4(26,1)
     WTUV(n,  42) = WTUVTEMP3(17,1) + Rpq(n,1)*WTUVTEMP4(27,1)
     WTUV(n,  43) = Rpq(n,3)*WTUVTEMP4(24,1)
     WTUV(n,  44) = Rpq(n,2)*WTUVTEMP4(26,1)
     WTUV(n,  45) = WTUVTEMP3(20,1) + Rpq(n,1)*WTUVTEMP4(30,1)
     WTUV(n,  46) = Rpq(n,1)*WTUVTEMP4(31,1)
     WTUV(n,  47) = Rpq(n,1)*WTUVTEMP4(32,1)
     WTUV(n,  48) = Rpq(n,1)*WTUVTEMP4(33,1)
     WTUV(n,  49) = Rpq(n,1)*WTUVTEMP4(34,1)
     WTUV(n,  50) = Rpq(n,1)*WTUVTEMP4(35,1)
     WTUV(n,  51) = 4*WTUVTEMP3(17,1) + Rpq(n,2)*WTUVTEMP4(31,1)
     WTUV(n,  52) = Rpq(n,3)*WTUVTEMP4(31,1)
     WTUV(n,  53) = 2*WTUVTEMP3(19,1) + Rpq(n,2)*WTUVTEMP4(33,1)
     WTUV(n,  54) = WTUVTEMP3(20,1) + Rpq(n,2)*WTUVTEMP4(34,1)
     WTUV(n,  55) = Rpq(n,2)*WTUVTEMP4(35,1)
     WTUV(n,  56) = 4*WTUVTEMP3(20,1) + Rpq(n,3)*WTUVTEMP4(35,1)

     WTUVTEMP5(36,1) = 4*WTUVTEMP3(11,2) + Rpq(n,1)*WTUVTEMP4(21,2)
!     WTUVTEMP5(37,1) = Rpq(n,2)*WTUVTEMP4(21,2) !not used
     WTUVTEMP5(38,1) = Rpq(n,3)*WTUVTEMP4(21,2)
     WTUVTEMP5(39,1) = 2*WTUVTEMP3(14,2) + Rpq(n,1)*WTUVTEMP4(24,2)
!     WTUVTEMP5(40,1) = Rpq(n,2)*WTUVTEMP4(23,2)!not used
     WTUVTEMP5(41,1) = 2*WTUVTEMP3(16,2) + Rpq(n,1)*WTUVTEMP4(26,2)
     WTUVTEMP5(42,1) = WTUVTEMP3(17,2) + Rpq(n,1)*WTUVTEMP4(27,2)
!     WTUVTEMP5(43,1) = Rpq(n,3)*WTUVTEMP4(24,2)!not used
!     WTUVTEMP5(44,1) = Rpq(n,2)*WTUVTEMP4(26,2)!not used
     WTUVTEMP5(45,1) = WTUVTEMP3(20,2) + Rpq(n,1)*WTUVTEMP4(30,2)
     WTUVTEMP5(46,1) = Rpq(n,1)*WTUVTEMP4(31,2)
!     WTUVTEMP5(47,1) = Rpq(n,1)*WTUVTEMP4(32,2)!not used
     WTUVTEMP5(48,1) = Rpq(n,1)*WTUVTEMP4(33,2)
!     WTUVTEMP5(49,1) = Rpq(n,1)*WTUVTEMP4(34,2)!not used
     WTUVTEMP5(50,1) = Rpq(n,1)*WTUVTEMP4(35,2)
     WTUVTEMP5(51,1) = 4*WTUVTEMP3(17,2) + Rpq(n,2)*WTUVTEMP4(31,2)
     WTUVTEMP5(52,1) = Rpq(n,3)*WTUVTEMP4(31,2)
     WTUVTEMP5(53,1) = 2*WTUVTEMP3(19,2) + Rpq(n,2)*WTUVTEMP4(33,2)
     WTUVTEMP5(54,1) = WTUVTEMP3(20,2) + Rpq(n,2)*WTUVTEMP4(34,2)
     WTUVTEMP5(55,1) = Rpq(n,2)*WTUVTEMP4(35,2)
     WTUVTEMP5(56,1) = 4*WTUVTEMP3(20,2) + Rpq(n,3)*WTUVTEMP4(35,2)

     WTUV(n,57) = 5*WTUVTEMP4(21,1) + Rpq(n,1)*WTUVTEMP5(36,1)
     WTUV(n,58) = Rpq(n,2)*WTUVTEMP5(36,1)
     WTUV(n,59) = Rpq(n,3)*WTUVTEMP5(36,1)
     WTUV(n,60) = 3*WTUVTEMP4(24,1) + Rpq(n,1)*WTUVTEMP5(39,1)
     WTUV(n,61) = Rpq(n,2)*WTUVTEMP5(38,1)
     WTUV(n,62) = 3*WTUVTEMP4(26,1) + Rpq(n,1)*WTUVTEMP5(41,1)
     WTUV(n,63) = 2*WTUVTEMP4(27,1) + Rpq(n,1)*WTUVTEMP5(42,1)
     WTUV(n,64) = Rpq(n,3)*WTUVTEMP5(39,1)
     WTUV(n,65) = Rpq(n,2)*WTUVTEMP5(41,1)
     WTUV(n,66) = 2*WTUVTEMP4(30,1) + Rpq(n,1)*WTUVTEMP5(45,1)
     WTUV(n,67) = WTUVTEMP4(31,1) + Rpq(n,1)*WTUVTEMP5(46,1)
     WTUV(n,68) = Rpq(n,3)*WTUVTEMP5(42,1)
     WTUV(n,69) = WTUVTEMP4(33,1) + Rpq(n,1)*WTUVTEMP5(48,1)
     WTUV(n,70) = Rpq(n,2)*WTUVTEMP5(45,1)
     WTUV(n,71) = WTUVTEMP4(35,1) + Rpq(n,1)*WTUVTEMP5(50,1)
     WTUV(n,72) = Rpq(n,1)*WTUVTEMP5(51,1)
     WTUV(n,73) = Rpq(n,1)*WTUVTEMP5(52,1)
     WTUV(n,74) = Rpq(n,1)*WTUVTEMP5(53,1)
     WTUV(n,75) = Rpq(n,1)*WTUVTEMP5(54,1)
     WTUV(n,76) = Rpq(n,1)*WTUVTEMP5(55,1)
     WTUV(n,77) = Rpq(n,1)*WTUVTEMP5(56,1)
     WTUV(n,78) = 5*WTUVTEMP4(31,1) + Rpq(n,2)*WTUVTEMP5(51,1)
     WTUV(n,79) = Rpq(n,3)*WTUVTEMP5(51,1)
     WTUV(n,80) = 3*WTUVTEMP4(33,1) + Rpq(n,2)*WTUVTEMP5(53,1)
     WTUV(n,81) = 2*WTUVTEMP4(34,1) + Rpq(n,2)*WTUVTEMP5(54,1)
     WTUV(n,82) = WTUVTEMP4(35,1) + Rpq(n,2)*WTUVTEMP5(55,1)
     WTUV(n,83) = Rpq(n,2)*WTUVTEMP5(56,1)
     WTUV(n,84) = 5*WTUVTEMP4(35,1) + Rpq(n,3)*WTUVTEMP5(56,1)
  ENDDO
END SUBROUTINE wtuvRecurrenceJMIN0JMAX6

SUBROUTINE wtuvRecurrenceJMIN0JMAX7(WTUV,WJ000,Rpq,nPrim)
  use memory_handling
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:7,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,120)
  !
  Integer     :: n,M
  Real(realk) :: WTUVTEMP(2:4,6)
  Real(realk) :: WTUVTEMP2(5:10,5)
  Real(realk) :: WTUVTEMP3(11:20,4)
  Real(realk) :: WTUVTEMP4(21:35,3)
  Real(realk) :: WTUVTEMP5(36:56,2)
  Real(realk) :: WTUVTEMP6(57:84,1)
  !000
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
     DO M=1,6
        WTUVTEMP(2,M) = Rpq(n,1)*WJ000(M+1,n)
        WTUVTEMP(3,M) = Rpq(n,2)*WJ000(M+1,n)
        WTUVTEMP(4,M) = Rpq(n,3)*WJ000(M+1,n)     
     ENDDO
     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
     DO M=1,5
        WTUVTEMP2(5,M) = WJ000(M+1,n) + Rpq(n,1)*WTUVTEMP(2,M+1)
!        WTUVTEMP2(6,M) =                Rpq(n,1)*WTUVTEMP(3,M+1)
!        WTUVTEMP2(7,M) =                Rpq(n,1)*WTUVTEMP(4,M+1)
        WTUVTEMP2(8,M) = WJ000(M+1,n) + Rpq(n,2)*WTUVTEMP(3,M+1)
        WTUVTEMP2(9,M) =                Rpq(n,2)*WTUVTEMP(4,M+1)
        WTUVTEMP2(10,M)= WJ000(M+1,n) + Rpq(n,3)*WTUVTEMP(4,M+1)     
     ENDDO
     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
     DO M=1,4
        WTUVTEMP3(11,M)= 2*WTUVTEMP(2,M+1) + Rpq(n,1)*WTUVTEMP2(5,M+1)
!        WTUVTEMP3(12,M)=  Rpq(n,2)*WTUVTEMP2(5,M+1)
        WTUVTEMP3(13,M)= Rpq(n,3)*WTUVTEMP2(5,M+1)
        WTUVTEMP3(14,M)= Rpq(n,1)*WTUVTEMP2(8,M+1)
!        WTUVTEMP3(15,M)= Rpq(n,1)*WTUVTEMP2(9,M+1)
        WTUVTEMP3(16,M)= Rpq(n,1)*WTUVTEMP2(10,M+1)
        WTUVTEMP3(17,M)= 2*WTUVTEMP(3,M+1) + Rpq(n,2)*WTUVTEMP2(8,M+1)
        WTUVTEMP3(18,M)= Rpq(n,3)*WTUVTEMP2(8,M+1)
        WTUVTEMP3(19,M)= Rpq(n,2)*WTUVTEMP2(10,M+1)
        WTUVTEMP3(20,M)= 2*WTUVTEMP(4,M+1) + Rpq(n,3)*WTUVTEMP2(10,M+1)
     ENDDO

     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11,1)
     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11,1)
     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11,1)
     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14,1)
     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13,1)
     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16,1)
     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17,1)
     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18,1)
     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19,1)
     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20,1)
     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17,1) 
     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17,1)
     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19,1)
     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20,1)
     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20,1)

     DO M=1,3
        WTUVTEMP4(21,M)= 3*WTUVTEMP2(5,M+1) + Rpq(n,1)*WTUVTEMP3(11,M+1)
        !     WTUVTEMP4(22,M)= Rpq(n,2)*WTUVTEMP3(11,M+1)
        WTUVTEMP4(23,M)= Rpq(n,3)*WTUVTEMP3(11,M+1)
        WTUVTEMP4(24,M)= WTUVTEMP2(8,M+1) + Rpq(n,1)*WTUVTEMP3(14,M+1)
        !     WTUVTEMP4(25,M)= Rpq(n,2)*WTUVTEMP3(13,M+1)
        WTUVTEMP4(26,M)= WTUVTEMP2(10,M+1) + Rpq(n,1)*WTUVTEMP3(16,M+1)
        WTUVTEMP4(27,M)= Rpq(n,1)*WTUVTEMP3(17,M+1)
        !     WTUVTEMP4(28,M)= Rpq(n,1)*WTUVTEMP3(18,M+1)
        !     WTUVTEMP4(29,M)= Rpq(n,1)*WTUVTEMP3(19,M+1)
        WTUVTEMP4(30,M)= Rpq(n,1)*WTUVTEMP3(20,M+1)
        WTUVTEMP4(31,M)= 3*WTUVTEMP2(8,M+1) + Rpq(n,2)*WTUVTEMP3(17,M+1) 
        WTUVTEMP4(32,M)= Rpq(n,3)*WTUVTEMP3(17,M+1)
        WTUVTEMP4(33,M)= WTUVTEMP2(10,M+1) + Rpq(n,2)*WTUVTEMP3(19,M+1)
        WTUVTEMP4(34,M)= Rpq(n,2)*WTUVTEMP3(20,M+1)
        WTUVTEMP4(35,M)= 3*WTUVTEMP2(10,M+1) + Rpq(n,3)*WTUVTEMP3(20,M+1)
     ENDDO

     WTUV(n,  36) = 4*WTUVTEMP3(11,1) + Rpq(n,1)*WTUVTEMP4(21,1)
     WTUV(n,  37) = Rpq(n,2)*WTUVTEMP4(21,1)
     WTUV(n,  38) = Rpq(n,3)*WTUVTEMP4(21,1)
     WTUV(n,  39) = 2*WTUVTEMP3(14,1) + Rpq(n,1)*WTUVTEMP4(24,1)
     WTUV(n,  40) = Rpq(n,2)*WTUVTEMP4(23,1)
     WTUV(n,  41) = 2*WTUVTEMP3(16,1) + Rpq(n,1)*WTUVTEMP4(26,1)
     WTUV(n,  42) = WTUVTEMP3(17,1) + Rpq(n,1)*WTUVTEMP4(27,1)
     WTUV(n,  43) = Rpq(n,3)*WTUVTEMP4(24,1)
     WTUV(n,  44) = Rpq(n,2)*WTUVTEMP4(26,1)
     WTUV(n,  45) = WTUVTEMP3(20,1) + Rpq(n,1)*WTUVTEMP4(30,1)
     WTUV(n,  46) = Rpq(n,1)*WTUVTEMP4(31,1)
     WTUV(n,  47) = Rpq(n,1)*WTUVTEMP4(32,1)
     WTUV(n,  48) = Rpq(n,1)*WTUVTEMP4(33,1)
     WTUV(n,  49) = Rpq(n,1)*WTUVTEMP4(34,1)
     WTUV(n,  50) = Rpq(n,1)*WTUVTEMP4(35,1)
     WTUV(n,  51) = 4*WTUVTEMP3(17,1) + Rpq(n,2)*WTUVTEMP4(31,1)
     WTUV(n,  52) = Rpq(n,3)*WTUVTEMP4(31,1)
     WTUV(n,  53) = 2*WTUVTEMP3(19,1) + Rpq(n,2)*WTUVTEMP4(33,1)
     WTUV(n,  54) = WTUVTEMP3(20,1) + Rpq(n,2)*WTUVTEMP4(34,1)
     WTUV(n,  55) = Rpq(n,2)*WTUVTEMP4(35,1)
     WTUV(n,  56) = 4*WTUVTEMP3(20,1) + Rpq(n,3)*WTUVTEMP4(35,1)

     DO M=1,2
        WTUVTEMP5(36,M) = 4*WTUVTEMP3(11,M+1) + Rpq(n,1)*WTUVTEMP4(21,M+1)
        !     WTUVTEMP5(37,M) = Rpq(n,2)*WTUVTEMP4(21,M+1)
        WTUVTEMP5(38,M) = Rpq(n,3)*WTUVTEMP4(21,M+1)
        WTUVTEMP5(39,M) = 2*WTUVTEMP3(14,M+1) + Rpq(n,1)*WTUVTEMP4(24,M+1)
        !     WTUVTEMP5(40,M) = Rpq(n,2)*WTUVTEMP4(23,M+1)
        WTUVTEMP5(41,M) = 2*WTUVTEMP3(16,M+1) + Rpq(n,1)*WTUVTEMP4(26,M+1)
        WTUVTEMP5(42,M) = WTUVTEMP3(17,M+1) + Rpq(n,1)*WTUVTEMP4(27,M+1)
        !     WTUVTEMP5(43,M) = Rpq(n,3)*WTUVTEMP4(24,M+1)
        !     WTUVTEMP5(44,M) = Rpq(n,2)*WTUVTEMP4(26,M+1)
        WTUVTEMP5(45,M) = WTUVTEMP3(20,M+1) + Rpq(n,1)*WTUVTEMP4(30,M+1)
        WTUVTEMP5(46,M) = Rpq(n,1)*WTUVTEMP4(31,M+1)
        !     WTUVTEMP5(47,M) = Rpq(n,1)*WTUVTEMP4(32,M+1)
        WTUVTEMP5(48,M) = Rpq(n,1)*WTUVTEMP4(33,M+1)
        !     WTUVTEMP5(49,M) = Rpq(n,1)*WTUVTEMP4(34,M+1)
        WTUVTEMP5(50,M) = Rpq(n,1)*WTUVTEMP4(35,M+1)
        WTUVTEMP5(51,M) = 4*WTUVTEMP3(17,M+1) + Rpq(n,2)*WTUVTEMP4(31,M+1)
        WTUVTEMP5(52,M) = Rpq(n,3)*WTUVTEMP4(31,M+1)
        WTUVTEMP5(53,M) = 2*WTUVTEMP3(19,M+1) + Rpq(n,2)*WTUVTEMP4(33,M+1)
        WTUVTEMP5(54,M) = WTUVTEMP3(20,M+1) + Rpq(n,2)*WTUVTEMP4(34,M+1)
        WTUVTEMP5(55,M) = Rpq(n,2)*WTUVTEMP4(35,M+1)
        WTUVTEMP5(56,M) = 4*WTUVTEMP3(20,M+1) + Rpq(n,3)*WTUVTEMP4(35,M+1)
     ENDDO
     WTUV(n,57) = 5*WTUVTEMP4(21,1) + Rpq(n,1)*WTUVTEMP5(36,1)
     WTUV(n,58) = Rpq(n,2)*WTUVTEMP5(36,1)
     WTUV(n,59) = Rpq(n,3)*WTUVTEMP5(36,1)
     WTUV(n,60) = 3*WTUVTEMP4(24,1) + Rpq(n,1)*WTUVTEMP5(39,1)
     WTUV(n,61) = Rpq(n,2)*WTUVTEMP5(38,1)
     WTUV(n,62) = 3*WTUVTEMP4(26,1) + Rpq(n,1)*WTUVTEMP5(41,1)
     WTUV(n,63) = 2*WTUVTEMP4(27,1) + Rpq(n,1)*WTUVTEMP5(42,1)
     WTUV(n,64) = Rpq(n,3)*WTUVTEMP5(39,1)
     WTUV(n,65) = Rpq(n,2)*WTUVTEMP5(41,1)
     WTUV(n,66) = 2*WTUVTEMP4(30,1) + Rpq(n,1)*WTUVTEMP5(45,1)
     WTUV(n,67) = WTUVTEMP4(31,1) + Rpq(n,1)*WTUVTEMP5(46,1)
     WTUV(n,68) = Rpq(n,3)*WTUVTEMP5(42,1)
     WTUV(n,69) = WTUVTEMP4(33,1) + Rpq(n,1)*WTUVTEMP5(48,1)
     WTUV(n,70) = Rpq(n,2)*WTUVTEMP5(45,1)
     WTUV(n,71) = WTUVTEMP4(35,1) + Rpq(n,1)*WTUVTEMP5(50,1)
     WTUV(n,72) = Rpq(n,1)*WTUVTEMP5(51,1)
     WTUV(n,73) = Rpq(n,1)*WTUVTEMP5(52,1)
     WTUV(n,74) = Rpq(n,1)*WTUVTEMP5(53,1)
     WTUV(n,75) = Rpq(n,1)*WTUVTEMP5(54,1)
     WTUV(n,76) = Rpq(n,1)*WTUVTEMP5(55,1)
     WTUV(n,77) = Rpq(n,1)*WTUVTEMP5(56,1)
     WTUV(n,78) = 5*WTUVTEMP4(31,1) + Rpq(n,2)*WTUVTEMP5(51,1)
     WTUV(n,79) = Rpq(n,3)*WTUVTEMP5(51,1)
     WTUV(n,80) = 3*WTUVTEMP4(33,1) + Rpq(n,2)*WTUVTEMP5(53,1)
     WTUV(n,81) = 2*WTUVTEMP4(34,1) + Rpq(n,2)*WTUVTEMP5(54,1)
     WTUV(n,82) = WTUVTEMP4(35,1) + Rpq(n,2)*WTUVTEMP5(55,1)
     WTUV(n,83) = Rpq(n,2)*WTUVTEMP5(56,1)
     WTUV(n,84) = 5*WTUVTEMP4(35,1) + Rpq(n,3)*WTUVTEMP5(56,1)
     
     WTUVTEMP6(57,1) = 5*WTUVTEMP4(21,2) + Rpq(n,1)*WTUVTEMP5(36,2)
!     WTUVTEMP6(58,1) = Rpq(n,2)*WTUVTEMP5(36,2)
     WTUVTEMP6(59,1) = Rpq(n,3)*WTUVTEMP5(36,2)
     WTUVTEMP6(60,1) = 3*WTUVTEMP4(24,2) + Rpq(n,1)*WTUVTEMP5(39,2)
!     WTUVTEMP6(61,1) = Rpq(n,2)*WTUVTEMP5(38,2)
     WTUVTEMP6(62,1) = 3*WTUVTEMP4(26,2) + Rpq(n,1)*WTUVTEMP5(41,2)
     WTUVTEMP6(63,1) = 2*WTUVTEMP4(27,2) + Rpq(n,1)*WTUVTEMP5(42,2)
!     WTUVTEMP6(64,1) = Rpq(n,3)*WTUVTEMP5(39,2)
!     WTUVTEMP6(65,1) = Rpq(n,2)*WTUVTEMP5(41,2)
     WTUVTEMP6(66,1) = 2*WTUVTEMP4(30,2) + Rpq(n,1)*WTUVTEMP5(45,2)
     WTUVTEMP6(67,1) = WTUVTEMP4(31,2) + Rpq(n,1)*WTUVTEMP5(46,2)
!     WTUVTEMP6(68,1) = Rpq(n,3)*WTUVTEMP5(42,2)
     WTUVTEMP6(69,1) = WTUVTEMP4(33,2) + Rpq(n,1)*WTUVTEMP5(48,2)
!     WTUVTEMP6(70,1) = Rpq(n,2)*WTUVTEMP5(45,2)
     WTUVTEMP6(71,1) = WTUVTEMP4(35,2) + Rpq(n,1)*WTUVTEMP5(50,2)
     WTUVTEMP6(72,1) = Rpq(n,1)*WTUVTEMP5(51,2)
!     WTUVTEMP6(73,1) = Rpq(n,1)*WTUVTEMP5(52,2)
     WTUVTEMP6(74,1) = Rpq(n,1)*WTUVTEMP5(53,2)
     WTUVTEMP6(75,1) = Rpq(n,1)*WTUVTEMP5(54,2)
!     WTUVTEMP6(76,1) = Rpq(n,1)*WTUVTEMP5(55,2)
     WTUVTEMP6(77,1) = Rpq(n,1)*WTUVTEMP5(56,2)
     WTUVTEMP6(78,1) = 5*WTUVTEMP4(31,2) + Rpq(n,2)*WTUVTEMP5(51,2)
     WTUVTEMP6(79,1) = Rpq(n,3)*WTUVTEMP5(51,2)
     WTUVTEMP6(80,1) = 3*WTUVTEMP4(33,2) + Rpq(n,2)*WTUVTEMP5(53,2)
     WTUVTEMP6(81,1) = 2*WTUVTEMP4(34,2) + Rpq(n,2)*WTUVTEMP5(54,2)
     WTUVTEMP6(82,1) = WTUVTEMP4(35,2) + Rpq(n,2)*WTUVTEMP5(55,2)
     WTUVTEMP6(83,1) = Rpq(n,2)*WTUVTEMP5(56,2)
     WTUVTEMP6(84,1) = 5*WTUVTEMP4(35,2) + Rpq(n,3)*WTUVTEMP5(56,2)

     WTUV(n, 85) = 6*WTUVTEMP5(36,1) + Rpq(n,1)*WTUVTEMP6(57,1)
     WTUV(n, 86) = Rpq(n,2)*WTUVTEMP6(57,1)
     WTUV(n, 87) = Rpq(n,3)*WTUVTEMP6(57,1)
     WTUV(n, 88) = 4*WTUVTEMP5(39,1) + Rpq(n,1)*WTUVTEMP6(60,1)
     WTUV(n, 89) = Rpq(n,2)*WTUVTEMP6(59,1)
     WTUV(n, 90) = 4*WTUVTEMP5(41,1) + Rpq(n,1)*WTUVTEMP6(62,1)
     WTUV(n, 91) = 3*WTUVTEMP5(42,1) + Rpq(n,1)*WTUVTEMP6(63,1)
     WTUV(n, 92) = Rpq(n,3)*WTUVTEMP6(60,1)
     WTUV(n, 93) = Rpq(n,2)*WTUVTEMP6(62,1)
     WTUV(n, 94) = 3*WTUVTEMP5(45,1) + Rpq(n,1)*WTUVTEMP6(66,1)
     WTUV(n, 95) = 2*WTUVTEMP5(46,1) + Rpq(n,1)*WTUVTEMP6(67,1)
     WTUV(n, 96) = Rpq(n,3)*WTUVTEMP6(63,1)
     WTUV(n, 97) = 2*WTUVTEMP5(48,1) + Rpq(n,1)*WTUVTEMP6(69,1)
     WTUV(n, 98) = Rpq(n,2)*WTUVTEMP6(66,1)
     WTUV(n, 99) = 2*WTUVTEMP5(50,1) + Rpq(n,1)*WTUVTEMP6(71,1)
     WTUV(n,100) = WTUVTEMP5(51,1) + Rpq(n,1)*WTUVTEMP6(72,1)
     WTUV(n,101) = Rpq(n,3)*WTUVTEMP6(67,1)
     WTUV(n,102) = WTUVTEMP5(53,1) + Rpq(n,1)*WTUVTEMP6(74,1)
     WTUV(n,103) = WTUVTEMP5(54,1) + Rpq(n,1)*WTUVTEMP6(75,1)
     WTUV(n,104) = Rpq(n,2)*WTUVTEMP6(71,1)
     WTUV(n,105) = WTUVTEMP5(56,1) + Rpq(n,1)*WTUVTEMP6(77,1)
     WTUV(n,106) = Rpq(n,1)*WTUVTEMP6(78,1)
     WTUV(n,107) = Rpq(n,1)*WTUVTEMP6(79,1)
     WTUV(n,108) = Rpq(n,1)*WTUVTEMP6(80,1)
     WTUV(n,109) = Rpq(n,1)*WTUVTEMP6(81,1)
     WTUV(n,110) = Rpq(n,1)*WTUVTEMP6(82,1)
     WTUV(n,111) = Rpq(n,1)*WTUVTEMP6(83,1)
     WTUV(n,112) = Rpq(n,1)*WTUVTEMP6(84,1)
     WTUV(n,113) = 6*WTUVTEMP5(51,1) + Rpq(n,2)*WTUVTEMP6(78,1)
     WTUV(n,114) = Rpq(n,3)*WTUVTEMP6(78,1)
     WTUV(n,115) = 4*WTUVTEMP5(53,1) + Rpq(n,2)*WTUVTEMP6(80,1)
     WTUV(n,116) = 3*WTUVTEMP5(54,1) + Rpq(n,2)*WTUVTEMP6(81,1)
     WTUV(n,117) = 2*WTUVTEMP5(55,1) + Rpq(n,2)*WTUVTEMP6(82,1)
     WTUV(n,118) = WTUVTEMP5(56,1) + Rpq(n,2)*WTUVTEMP6(83,1)
     WTUV(n,119) = Rpq(n,2)*WTUVTEMP6(84,1)
     WTUV(n,120) = 6*WTUVTEMP5(56,1) + Rpq(n,3)*WTUVTEMP6(84,1)
  ENDDO
END SUBROUTINE wtuvRecurrenceJMIN0JMAX7

SUBROUTINE wtuvRecurrenceJMIN0JMAX8(WTUV,WJ000,Rpq,nPrim)
  use memory_handling
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:8,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,165)
  !
  Integer     :: n,M
  Real(realk) :: WTUVTEMP(2:4,7)
  Real(realk) :: WTUVTEMP2(5:10,6)
  Real(realk) :: WTUVTEMP3(11:20,5)
  Real(realk) :: WTUVTEMP4(21:35,4)
  Real(realk) :: WTUVTEMP5(36:56,3)
  Real(realk) :: WTUVTEMP6(57:84,2)
  Real(realk) :: WTUVTEMP7(85:120,1)
  !000
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
     DO M=1,7
        WTUVTEMP(2,M) = Rpq(n,1)*WJ000(M+1,n)
        WTUVTEMP(3,M) = Rpq(n,2)*WJ000(M+1,n)
        WTUVTEMP(4,M) = Rpq(n,3)*WJ000(M+1,n)     
     ENDDO
     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
     DO M=1,6
        WTUVTEMP2(5,M) = WJ000(M+1,n) + Rpq(n,1)*WTUVTEMP(2,M+1)
!        WTUVTEMP2(6,M) =                Rpq(n,1)*WTUVTEMP(3,M+1)
!        WTUVTEMP2(7,M) =                Rpq(n,1)*WTUVTEMP(4,M+1)
        WTUVTEMP2(8,M) = WJ000(M+1,n) + Rpq(n,2)*WTUVTEMP(3,M+1)
        WTUVTEMP2(9,M) =                Rpq(n,2)*WTUVTEMP(4,M+1)
        WTUVTEMP2(10,M)= WJ000(M+1,n) + Rpq(n,3)*WTUVTEMP(4,M+1)     
     ENDDO
     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
     DO M=1,5
        WTUVTEMP3(11,M)= 2*WTUVTEMP(2,M+1) + Rpq(n,1)*WTUVTEMP2(5,M+1)
!        WTUVTEMP3(12,M)=  Rpq(n,2)*WTUVTEMP2(5,M+1)
        WTUVTEMP3(13,M)= Rpq(n,3)*WTUVTEMP2(5,M+1)
        WTUVTEMP3(14,M)= Rpq(n,1)*WTUVTEMP2(8,M+1)
!        WTUVTEMP3(15,M)= Rpq(n,1)*WTUVTEMP2(9,M+1)
        WTUVTEMP3(16,M)= Rpq(n,1)*WTUVTEMP2(10,M+1)
        WTUVTEMP3(17,M)= 2*WTUVTEMP(3,M+1) + Rpq(n,2)*WTUVTEMP2(8,M+1)
        WTUVTEMP3(18,M)= Rpq(n,3)*WTUVTEMP2(8,M+1)
        WTUVTEMP3(19,M)= Rpq(n,2)*WTUVTEMP2(10,M+1)
        WTUVTEMP3(20,M)= 2*WTUVTEMP(4,M+1) + Rpq(n,3)*WTUVTEMP2(10,M+1)
     ENDDO

     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11,1)
     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11,1)
     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11,1)
     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14,1)
     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13,1)
     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16,1)
     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17,1)
     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18,1)
     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19,1)
     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20,1)
     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17,1) 
     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17,1)
     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19,1)
     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20,1)
     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20,1)

     DO M=1,4
        WTUVTEMP4(21,M)= 3*WTUVTEMP2(5,M+1) + Rpq(n,1)*WTUVTEMP3(11,M+1)
        !     WTUVTEMP4(22,M)= Rpq(n,2)*WTUVTEMP3(11,M+1)
        WTUVTEMP4(23,M)= Rpq(n,3)*WTUVTEMP3(11,M+1)
        WTUVTEMP4(24,M)= WTUVTEMP2(8,M+1) + Rpq(n,1)*WTUVTEMP3(14,M+1)
        !     WTUVTEMP4(25,M)= Rpq(n,2)*WTUVTEMP3(13,M+1)
        WTUVTEMP4(26,M)= WTUVTEMP2(10,M+1) + Rpq(n,1)*WTUVTEMP3(16,M+1)
        WTUVTEMP4(27,M)= Rpq(n,1)*WTUVTEMP3(17,M+1)
        !     WTUVTEMP4(28,M)= Rpq(n,1)*WTUVTEMP3(18,M+1)
        !     WTUVTEMP4(29,M)= Rpq(n,1)*WTUVTEMP3(19,M+1)
        WTUVTEMP4(30,M)= Rpq(n,1)*WTUVTEMP3(20,M+1)
        WTUVTEMP4(31,M)= 3*WTUVTEMP2(8,M+1) + Rpq(n,2)*WTUVTEMP3(17,M+1) 
        WTUVTEMP4(32,M)= Rpq(n,3)*WTUVTEMP3(17,M+1)
        WTUVTEMP4(33,M)= WTUVTEMP2(10,M+1) + Rpq(n,2)*WTUVTEMP3(19,M+1)
        WTUVTEMP4(34,M)= Rpq(n,2)*WTUVTEMP3(20,M+1)
        WTUVTEMP4(35,M)= 3*WTUVTEMP2(10,M+1) + Rpq(n,3)*WTUVTEMP3(20,M+1)
     ENDDO

     WTUV(n,  36) = 4*WTUVTEMP3(11,1) + Rpq(n,1)*WTUVTEMP4(21,1)
     WTUV(n,  37) = Rpq(n,2)*WTUVTEMP4(21,1)
     WTUV(n,  38) = Rpq(n,3)*WTUVTEMP4(21,1)
     WTUV(n,  39) = 2*WTUVTEMP3(14,1) + Rpq(n,1)*WTUVTEMP4(24,1)
     WTUV(n,  40) = Rpq(n,2)*WTUVTEMP4(23,1)
     WTUV(n,  41) = 2*WTUVTEMP3(16,1) + Rpq(n,1)*WTUVTEMP4(26,1)
     WTUV(n,  42) = WTUVTEMP3(17,1) + Rpq(n,1)*WTUVTEMP4(27,1)
     WTUV(n,  43) = Rpq(n,3)*WTUVTEMP4(24,1)
     WTUV(n,  44) = Rpq(n,2)*WTUVTEMP4(26,1)
     WTUV(n,  45) = WTUVTEMP3(20,1) + Rpq(n,1)*WTUVTEMP4(30,1)
     WTUV(n,  46) = Rpq(n,1)*WTUVTEMP4(31,1)
     WTUV(n,  47) = Rpq(n,1)*WTUVTEMP4(32,1)
     WTUV(n,  48) = Rpq(n,1)*WTUVTEMP4(33,1)
     WTUV(n,  49) = Rpq(n,1)*WTUVTEMP4(34,1)
     WTUV(n,  50) = Rpq(n,1)*WTUVTEMP4(35,1)
     WTUV(n,  51) = 4*WTUVTEMP3(17,1) + Rpq(n,2)*WTUVTEMP4(31,1)
     WTUV(n,  52) = Rpq(n,3)*WTUVTEMP4(31,1)
     WTUV(n,  53) = 2*WTUVTEMP3(19,1) + Rpq(n,2)*WTUVTEMP4(33,1)
     WTUV(n,  54) = WTUVTEMP3(20,1) + Rpq(n,2)*WTUVTEMP4(34,1)
     WTUV(n,  55) = Rpq(n,2)*WTUVTEMP4(35,1)
     WTUV(n,  56) = 4*WTUVTEMP3(20,1) + Rpq(n,3)*WTUVTEMP4(35,1)

     DO M=1,3
        WTUVTEMP5(36,M) = 4*WTUVTEMP3(11,M+1) + Rpq(n,1)*WTUVTEMP4(21,M+1)
        !     WTUVTEMP5(37,M) = Rpq(n,2)*WTUVTEMP4(21,M+1)
        WTUVTEMP5(38,M) = Rpq(n,3)*WTUVTEMP4(21,M+1)
        WTUVTEMP5(39,M) = 2*WTUVTEMP3(14,M+1) + Rpq(n,1)*WTUVTEMP4(24,M+1)
        !     WTUVTEMP5(40,M) = Rpq(n,2)*WTUVTEMP4(23,M+1)
        WTUVTEMP5(41,M) = 2*WTUVTEMP3(16,M+1) + Rpq(n,1)*WTUVTEMP4(26,M+1)
        WTUVTEMP5(42,M) = WTUVTEMP3(17,M+1) + Rpq(n,1)*WTUVTEMP4(27,M+1)
        !     WTUVTEMP5(43,M) = Rpq(n,3)*WTUVTEMP4(24,M+1)
        !     WTUVTEMP5(44,M) = Rpq(n,2)*WTUVTEMP4(26,M+1)
        WTUVTEMP5(45,M) = WTUVTEMP3(20,M+1) + Rpq(n,1)*WTUVTEMP4(30,M+1)
        WTUVTEMP5(46,M) = Rpq(n,1)*WTUVTEMP4(31,M+1)
        !     WTUVTEMP5(47,M) = Rpq(n,1)*WTUVTEMP4(32,M+1)
        WTUVTEMP5(48,M) = Rpq(n,1)*WTUVTEMP4(33,M+1)
        !     WTUVTEMP5(49,M) = Rpq(n,1)*WTUVTEMP4(34,M+1)
        WTUVTEMP5(50,M) = Rpq(n,1)*WTUVTEMP4(35,M+1)
        WTUVTEMP5(51,M) = 4*WTUVTEMP3(17,M+1) + Rpq(n,2)*WTUVTEMP4(31,M+1)
        WTUVTEMP5(52,M) = Rpq(n,3)*WTUVTEMP4(31,M+1)
        WTUVTEMP5(53,M) = 2*WTUVTEMP3(19,M+1) + Rpq(n,2)*WTUVTEMP4(33,M+1)
        WTUVTEMP5(54,M) = WTUVTEMP3(20,M+1) + Rpq(n,2)*WTUVTEMP4(34,M+1)
        WTUVTEMP5(55,M) = Rpq(n,2)*WTUVTEMP4(35,M+1)
        WTUVTEMP5(56,M) = 4*WTUVTEMP3(20,M+1) + Rpq(n,3)*WTUVTEMP4(35,M+1)
     ENDDO
     WTUV(n,57) = 5*WTUVTEMP4(21,1) + Rpq(n,1)*WTUVTEMP5(36,1)
     WTUV(n,58) = Rpq(n,2)*WTUVTEMP5(36,1)
     WTUV(n,59) = Rpq(n,3)*WTUVTEMP5(36,1)
     WTUV(n,60) = 3*WTUVTEMP4(24,1) + Rpq(n,1)*WTUVTEMP5(39,1)
     WTUV(n,61) = Rpq(n,2)*WTUVTEMP5(38,1)
     WTUV(n,62) = 3*WTUVTEMP4(26,1) + Rpq(n,1)*WTUVTEMP5(41,1)
     WTUV(n,63) = 2*WTUVTEMP4(27,1) + Rpq(n,1)*WTUVTEMP5(42,1)
     WTUV(n,64) = Rpq(n,3)*WTUVTEMP5(39,1)
     WTUV(n,65) = Rpq(n,2)*WTUVTEMP5(41,1)
     WTUV(n,66) = 2*WTUVTEMP4(30,1) + Rpq(n,1)*WTUVTEMP5(45,1)
     WTUV(n,67) = WTUVTEMP4(31,1) + Rpq(n,1)*WTUVTEMP5(46,1)
     WTUV(n,68) = Rpq(n,3)*WTUVTEMP5(42,1)
     WTUV(n,69) = WTUVTEMP4(33,1) + Rpq(n,1)*WTUVTEMP5(48,1)
     WTUV(n,70) = Rpq(n,2)*WTUVTEMP5(45,1)
     WTUV(n,71) = WTUVTEMP4(35,1) + Rpq(n,1)*WTUVTEMP5(50,1)
     WTUV(n,72) = Rpq(n,1)*WTUVTEMP5(51,1)
     WTUV(n,73) = Rpq(n,1)*WTUVTEMP5(52,1)
     WTUV(n,74) = Rpq(n,1)*WTUVTEMP5(53,1)
     WTUV(n,75) = Rpq(n,1)*WTUVTEMP5(54,1)
     WTUV(n,76) = Rpq(n,1)*WTUVTEMP5(55,1)
     WTUV(n,77) = Rpq(n,1)*WTUVTEMP5(56,1)
     WTUV(n,78) = 5*WTUVTEMP4(31,1) + Rpq(n,2)*WTUVTEMP5(51,1)
     WTUV(n,79) = Rpq(n,3)*WTUVTEMP5(51,1)
     WTUV(n,80) = 3*WTUVTEMP4(33,1) + Rpq(n,2)*WTUVTEMP5(53,1)
     WTUV(n,81) = 2*WTUVTEMP4(34,1) + Rpq(n,2)*WTUVTEMP5(54,1)
     WTUV(n,82) = WTUVTEMP4(35,1) + Rpq(n,2)*WTUVTEMP5(55,1)
     WTUV(n,83) = Rpq(n,2)*WTUVTEMP5(56,1)
     WTUV(n,84) = 5*WTUVTEMP4(35,1) + Rpq(n,3)*WTUVTEMP5(56,1)
     
     DO M=1,2
        WTUVTEMP6(57,M) = 5*WTUVTEMP4(21,M+1) + Rpq(n,1)*WTUVTEMP5(36,M+1)
        !     WTUVTEMP6(58,M) = Rpq(n,2)*WTUVTEMP5(36,M+1)
        WTUVTEMP6(59,M) = Rpq(n,3)*WTUVTEMP5(36,M+1)
        WTUVTEMP6(60,M) = 3*WTUVTEMP4(24,M+1) + Rpq(n,1)*WTUVTEMP5(39,M+1)
        !     WTUVTEMP6(61,M) = Rpq(n,2)*WTUVTEMP5(38,M+1)
        WTUVTEMP6(62,M) = 3*WTUVTEMP4(26,M+1) + Rpq(n,1)*WTUVTEMP5(41,M+1)
        WTUVTEMP6(63,M) = 2*WTUVTEMP4(27,M+1) + Rpq(n,1)*WTUVTEMP5(42,M+1)
        !     WTUVTEMP6(64,M) = Rpq(n,3)*WTUVTEMP5(39,M+1)
        !     WTUVTEMP6(65,M) = Rpq(n,2)*WTUVTEMP5(41,M+1)
        WTUVTEMP6(66,M) = 2*WTUVTEMP4(30,M+1) + Rpq(n,1)*WTUVTEMP5(45,M+1)
        WTUVTEMP6(67,M) = WTUVTEMP4(31,M+1) + Rpq(n,1)*WTUVTEMP5(46,M+1)
        !     WTUVTEMP6(68,M) = Rpq(n,3)*WTUVTEMP5(42,M+1)
        WTUVTEMP6(69,M) = WTUVTEMP4(33,M+1) + Rpq(n,1)*WTUVTEMP5(48,M+1)
        !     WTUVTEMP6(70,M) = Rpq(n,M+1)*WTUVTEMP5(45,M+1)
        WTUVTEMP6(71,M) = WTUVTEMP4(35,M+1) + Rpq(n,1)*WTUVTEMP5(50,M+1)
        WTUVTEMP6(72,M) = Rpq(n,1)*WTUVTEMP5(51,M+1)
        !     WTUVTEMP6(73,M) = Rpq(n,1)*WTUVTEMP5(52,M+1)
        WTUVTEMP6(74,M) = Rpq(n,1)*WTUVTEMP5(53,M+1)
        WTUVTEMP6(75,M) = Rpq(n,1)*WTUVTEMP5(54,M+1)
        !     WTUVTEMP6(76,M) = Rpq(n,1)*WTUVTEMP5(55,M+1)
        WTUVTEMP6(77,M) = Rpq(n,1)*WTUVTEMP5(56,M+1)
        WTUVTEMP6(78,M) = 5*WTUVTEMP4(31,M+1) + Rpq(n,2)*WTUVTEMP5(51,M+1)
        WTUVTEMP6(79,M) = Rpq(n,3)*WTUVTEMP5(51,M+1)
        WTUVTEMP6(80,M) = 3*WTUVTEMP4(33,M+1) + Rpq(n,2)*WTUVTEMP5(53,M+1)
        WTUVTEMP6(81,M) = 2*WTUVTEMP4(34,M+1) + Rpq(n,2)*WTUVTEMP5(54,M+1)
        WTUVTEMP6(82,M) = WTUVTEMP4(35,M+1) + Rpq(n,2)*WTUVTEMP5(55,M+1)
        WTUVTEMP6(83,M) = Rpq(n,2)*WTUVTEMP5(56,M+1)
        WTUVTEMP6(84,M) = 5*WTUVTEMP4(35,M+1) + Rpq(n,3)*WTUVTEMP5(56,M+1)
     ENDDO

     WTUV(n, 85) = 6*WTUVTEMP5(36,1) + Rpq(n,1)*WTUVTEMP6(57,1)
     WTUV(n, 86) = Rpq(n,2)*WTUVTEMP6(57,1)
     WTUV(n, 87) = Rpq(n,3)*WTUVTEMP6(57,1)
     WTUV(n, 88) = 4*WTUVTEMP5(39,1) + Rpq(n,1)*WTUVTEMP6(60,1)
     WTUV(n, 89) = Rpq(n,2)*WTUVTEMP6(59,1)
     WTUV(n, 90) = 4*WTUVTEMP5(41,1) + Rpq(n,1)*WTUVTEMP6(62,1)
     WTUV(n, 91) = 3*WTUVTEMP5(42,1) + Rpq(n,1)*WTUVTEMP6(63,1)
     WTUV(n, 92) = Rpq(n,3)*WTUVTEMP6(60,1)
     WTUV(n, 93) = Rpq(n,2)*WTUVTEMP6(62,1)
     WTUV(n, 94) = 3*WTUVTEMP5(45,1) + Rpq(n,1)*WTUVTEMP6(66,1)
     WTUV(n, 95) = 2*WTUVTEMP5(46,1) + Rpq(n,1)*WTUVTEMP6(67,1)
     WTUV(n, 96) = Rpq(n,3)*WTUVTEMP6(63,1)
     WTUV(n, 97) = 2*WTUVTEMP5(48,1) + Rpq(n,1)*WTUVTEMP6(69,1)
     WTUV(n, 98) = Rpq(n,2)*WTUVTEMP6(66,1)
     WTUV(n, 99) = 2*WTUVTEMP5(50,1) + Rpq(n,1)*WTUVTEMP6(71,1)
     WTUV(n,100) = WTUVTEMP5(51,1) + Rpq(n,1)*WTUVTEMP6(72,1)
     WTUV(n,101) = Rpq(n,3)*WTUVTEMP6(67,1)
     WTUV(n,102) = WTUVTEMP5(53,1) + Rpq(n,1)*WTUVTEMP6(74,1)
     WTUV(n,103) = WTUVTEMP5(54,1) + Rpq(n,1)*WTUVTEMP6(75,1)
     WTUV(n,104) = Rpq(n,2)*WTUVTEMP6(71,1)
     WTUV(n,105) = WTUVTEMP5(56,1) + Rpq(n,1)*WTUVTEMP6(77,1)
     WTUV(n,106) = Rpq(n,1)*WTUVTEMP6(78,1)
     WTUV(n,107) = Rpq(n,1)*WTUVTEMP6(79,1)
     WTUV(n,108) = Rpq(n,1)*WTUVTEMP6(80,1)
     WTUV(n,109) = Rpq(n,1)*WTUVTEMP6(81,1)
     WTUV(n,110) = Rpq(n,1)*WTUVTEMP6(82,1)
     WTUV(n,111) = Rpq(n,1)*WTUVTEMP6(83,1)
     WTUV(n,112) = Rpq(n,1)*WTUVTEMP6(84,1)
     WTUV(n,113) = 6*WTUVTEMP5(51,1) + Rpq(n,2)*WTUVTEMP6(78,1)
     WTUV(n,114) = Rpq(n,3)*WTUVTEMP6(78,1)
     WTUV(n,115) = 4*WTUVTEMP5(53,1) + Rpq(n,2)*WTUVTEMP6(80,1)
     WTUV(n,116) = 3*WTUVTEMP5(54,1) + Rpq(n,2)*WTUVTEMP6(81,1)
     WTUV(n,117) = 2*WTUVTEMP5(55,1) + Rpq(n,2)*WTUVTEMP6(82,1)
     WTUV(n,118) = WTUVTEMP5(56,1) + Rpq(n,2)*WTUVTEMP6(83,1)
     WTUV(n,119) = Rpq(n,2)*WTUVTEMP6(84,1)
     WTUV(n,120) = 6*WTUVTEMP5(56,1) + Rpq(n,3)*WTUVTEMP6(84,1)

     WTUVTEMP7( 85,1) = 6*WTUVTEMP5(36,2) + Rpq(n,1)*WTUVTEMP6(57,2)
!     WTUVTEMP7( 86,1) = Rpq(n,2)*WTUVTEMP6(57,2)
     WTUVTEMP7( 87,1) = Rpq(n,3)*WTUVTEMP6(57,2)
     WTUVTEMP7( 88,1) = 4*WTUVTEMP5(39,2) + Rpq(n,1)*WTUVTEMP6(60,2)
!     WTUVTEMP7( 89,1) = Rpq(n,2)*WTUVTEMP6(59,2)
     WTUVTEMP7( 90,1) = 4*WTUVTEMP5(41,2) + Rpq(n,1)*WTUVTEMP6(62,2)
     WTUVTEMP7( 91,1) = 3*WTUVTEMP5(42,2) + Rpq(n,1)*WTUVTEMP6(63,2)
!     WTUVTEMP7( 92,1) = Rpq(n,3)*WTUVTEMP6(60,2)
!     WTUVTEMP7( 93,1) = Rpq(n,2)*WTUVTEMP6(62,2)
     WTUVTEMP7( 94,1) = 3*WTUVTEMP5(45,2) + Rpq(n,1)*WTUVTEMP6(66,2)
     WTUVTEMP7( 95,1) = 2*WTUVTEMP5(46,2) + Rpq(n,1)*WTUVTEMP6(67,2)
!     WTUVTEMP7( 96,1) = Rpq(n,3)*WTUVTEMP6(63,2)
     WTUVTEMP7( 97,1) = 2*WTUVTEMP5(48,2) + Rpq(n,1)*WTUVTEMP6(69,2)
!     WTUVTEMP7( 98,1) = Rpq(n,2)*WTUVTEMP6(66,2)
     WTUVTEMP7( 99,1) = 2*WTUVTEMP5(50,2) + Rpq(n,1)*WTUVTEMP6(71,2)
     WTUVTEMP7(100,1) = WTUVTEMP5(51,2) + Rpq(n,1)*WTUVTEMP6(72,2)
!     WTUVTEMP7(101,1) = Rpq(n,3)*WTUVTEMP6(67,2)
     WTUVTEMP7(102,1) = WTUVTEMP5(53,2) + Rpq(n,1)*WTUVTEMP6(74,2)
     WTUVTEMP7(103,1) = WTUVTEMP5(54,2) + Rpq(n,1)*WTUVTEMP6(75,2)
!     WTUVTEMP7(104,1) = Rpq(n,2)*WTUVTEMP6(71,2)
     WTUVTEMP7(105,1) = WTUVTEMP5(56,2) + Rpq(n,1)*WTUVTEMP6(77,2)
     WTUVTEMP7(106,1) = Rpq(n,1)*WTUVTEMP6(78,2)
!     WTUVTEMP7(107,1) = Rpq(n,1)*WTUVTEMP6(79,2)
     WTUVTEMP7(108,1) = Rpq(n,1)*WTUVTEMP6(80,2)
     WTUVTEMP7(109,1) = Rpq(n,1)*WTUVTEMP6(81,2)
     WTUVTEMP7(110,1) = Rpq(n,1)*WTUVTEMP6(82,2)
!     WTUVTEMP7(111,1) = Rpq(n,1)*WTUVTEMP6(83,2)
     WTUVTEMP7(112,1) = Rpq(n,1)*WTUVTEMP6(84,2)
     WTUVTEMP7(113,1) = 6*WTUVTEMP5(51,2) + Rpq(n,2)*WTUVTEMP6(78,2)
     WTUVTEMP7(114,1) = Rpq(n,3)*WTUVTEMP6(78,2)
     WTUVTEMP7(115,1) = 4*WTUVTEMP5(53,2) + Rpq(n,2)*WTUVTEMP6(80,2)
     WTUVTEMP7(116,1) = 3*WTUVTEMP5(54,2) + Rpq(n,2)*WTUVTEMP6(81,2)
     WTUVTEMP7(117,1) = 2*WTUVTEMP5(55,2) + Rpq(n,2)*WTUVTEMP6(82,2)
     WTUVTEMP7(118,1) = WTUVTEMP5(56,2) + Rpq(n,2)*WTUVTEMP6(83,2)
     WTUVTEMP7(119,1) = Rpq(n,2)*WTUVTEMP6(84,2)
     WTUVTEMP7(120,1) = 6*WTUVTEMP5(56,2) + Rpq(n,3)*WTUVTEMP6(84,2)

     WTUV(n, 121) = 7*WTUVTEMP6(  57,1) + Rpq(n,1)*WTUVTEMP7(  85,1)
     WTUV(n, 122) = Rpq(n,2)*WTUVTEMP7(  85,1)
     WTUV(n, 123) = Rpq(n,3)*WTUVTEMP7(  85,1)
     WTUV(n, 124) = 5*WTUVTEMP6(  60,1) + Rpq(n,1)*WTUVTEMP7(  88,1)
     WTUV(n, 125) = Rpq(n,2)*WTUVTEMP7(  87,1)
     WTUV(n, 126) = 5*WTUVTEMP6(  62,1) + Rpq(n,1)*WTUVTEMP7(  90,1)
     WTUV(n, 127) = 4*WTUVTEMP6(  63,1) + Rpq(n,1)*WTUVTEMP7(  91,1)
     WTUV(n, 128) = Rpq(n,3)*WTUVTEMP7(  88,1)
     WTUV(n, 129) = Rpq(n,2)*WTUVTEMP7(  90,1)
     WTUV(n, 130) = 4*WTUVTEMP6(  66,1) + Rpq(n,1)*WTUVTEMP7(  94,1)
     WTUV(n, 131) = 3*WTUVTEMP6(  67,1) + Rpq(n,1)*WTUVTEMP7(  95,1)
     WTUV(n, 132) = Rpq(n,3)*WTUVTEMP7(  91,1)
     WTUV(n, 133) = 3*WTUVTEMP6(  69,1) + Rpq(n,1)*WTUVTEMP7(  97,1)
     WTUV(n, 134) = Rpq(n,2)*WTUVTEMP7(  94,1)
     WTUV(n, 135) = 3*WTUVTEMP6(  71,1) + Rpq(n,1)*WTUVTEMP7(  99,1)
     WTUV(n, 136) = 2*WTUVTEMP6(  72,1) + Rpq(n,1)*WTUVTEMP7( 100,1)
     WTUV(n, 137) = Rpq(n,3)*WTUVTEMP7(  95,1)
     WTUV(n, 138) = 2*WTUVTEMP6(  74,1) + Rpq(n,1)*WTUVTEMP7( 102,1)
     WTUV(n, 139) = 2*WTUVTEMP6(  75,1) + Rpq(n,1)*WTUVTEMP7( 103,1)
     WTUV(n, 140) = Rpq(n,2)*WTUVTEMP7(  99,1)
     WTUV(n, 141) = 2*WTUVTEMP6(  77,1) + Rpq(n,1)*WTUVTEMP7( 105,1)
     WTUV(n, 142) = WTUVTEMP6(  78,1) + Rpq(n,1)*WTUVTEMP7( 106,1)
     WTUV(n, 143) = Rpq(n,3)*WTUVTEMP7( 100,1)
     WTUV(n, 144) = WTUVTEMP6(  80,1) + Rpq(n,1)*WTUVTEMP7( 108,1)
     WTUV(n, 145) = WTUVTEMP6(  81,1) + Rpq(n,1)*WTUVTEMP7( 109,1)
     WTUV(n, 146) = WTUVTEMP6(  82,1) + Rpq(n,1)*WTUVTEMP7( 110,1)
     WTUV(n, 147) = Rpq(n,2)*WTUVTEMP7( 105,1)
     WTUV(n, 148) = WTUVTEMP6(  84,1) + Rpq(n,1)*WTUVTEMP7( 112,1)
     WTUV(n, 149) = Rpq(n,1)*WTUVTEMP7( 113,1)
     WTUV(n, 150) = Rpq(n,1)*WTUVTEMP7( 114,1)
     WTUV(n, 151) = Rpq(n,1)*WTUVTEMP7( 115,1)
     WTUV(n, 152) = Rpq(n,1)*WTUVTEMP7( 116,1)
     WTUV(n, 153) = Rpq(n,1)*WTUVTEMP7( 117,1)
     WTUV(n, 154) = Rpq(n,1)*WTUVTEMP7( 118,1)
     WTUV(n, 155) = Rpq(n,1)*WTUVTEMP7( 119,1)
     WTUV(n, 156) = Rpq(n,1)*WTUVTEMP7( 120,1)
     WTUV(n, 157) = 7*WTUVTEMP6(  78,1) + Rpq(n,2)*WTUVTEMP7( 113,1)
     WTUV(n, 158) = Rpq(n,3)*WTUVTEMP7( 113,1)
     WTUV(n, 159) = 5*WTUVTEMP6(  80,1) + Rpq(n,2)*WTUVTEMP7( 115,1)
     WTUV(n, 160) = 4*WTUVTEMP6(  81,1) + Rpq(n,2)*WTUVTEMP7( 116,1)
     WTUV(n, 161) = 3*WTUVTEMP6(  82,1) + Rpq(n,2)*WTUVTEMP7( 117,1)
     WTUV(n, 162) = 2*WTUVTEMP6(  83,1) + Rpq(n,2)*WTUVTEMP7( 118,1)
     WTUV(n, 163) = WTUVTEMP6(  84,1) + Rpq(n,2)*WTUVTEMP7( 119,1)
     WTUV(n, 164) = Rpq(n,2)*WTUVTEMP7( 120,1)
     WTUV(n, 165) = 7*WTUVTEMP6(  84,1) + Rpq(n,3)*WTUVTEMP7( 120,1)
  ENDDO
END SUBROUTINE wtuvRecurrenceJMIN0JMAX8




!> \brief wrapper to Overlap hermite integral
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param SJ000 \f$ S^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE buildSJ000(SJ000,PQ,nPrim,Jmax,LUPRI,IPRINT)
 IMPLICIT NONE
 INTEGER,intent(in)          :: nPrim,LUPRI,IPRINT,Jmax
 TYPE(Integrand),intent(in)  :: PQ
 real(realk),intent(inout)   :: SJ000(0:Jmax,nPrim)
! INTEGER  :: J,I
 CALL buildSJ000_OP(SJ000,Nprim,Jmax,PQ%reducedExponents&
      &,PQ%exponents,PQ%squaredDistance)

! IF(IPRINT .GT. 25) THEN
!    DO I=1,NPrim
!       DO J=0,Jmax
!          WRITE(LUPRI,'(2X,A6,I4,A1,I2,A2,ES16.9)')'SJ000(',I,',',J,')=',SJ000(J,I)
!       ENDDO
!    ENDDO
! ENDIF
END SUBROUTINE buildSJ000

!> \brief Overlap hermite integral \f[ S^{j}_{000} = \left(-2 \alpha \right)^{j} \left( \frac{\gamma}{2\pi} \right)^{\frac{3}{2}} \exp(-\gamma R_{pq}^{2}) \f]
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param SJ000 \f$ S^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE buildSJ000_OP(SJ000,nPrim,Jmax,Alpha,P,R2)
use memory_handling
 IMPLICIT NONE
 INTEGER,intent(in)        :: nPrim,Jmax
 real(realk),intent(inout) :: SJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)    :: ALPHA(nPrim),P(nPrim),R2(nprim)
 !
 Real(realk), parameter :: PI=3.14159265358979323846E0_realk
 Real(realk), parameter :: OneHalf=1.5E0_realk, Two = 2.0E0_realk
 REAL(REALK),pointer    :: TEMP(:),TEMP2(:)
 INTEGER                :: J,I
 call mem_alloc(TEMP,nprim)
 call mem_alloc(TEMP2,nprim)
 DO I=1,NPrim
    TEMP(I)=Alpha(I)*R2(I)
 ENDDO
 DO I=1,NPrim
    TEMP2(I)=EXP(-TEMP(I))
 ENDDO
 DO I=1,NPrim
    DO J=0,Jmax
       SJ000(J,I)=((-Two*Alpha(I))**J)&
            &*((PI/P(I))**OneHalf)*TEMP2(I)
    ENDDO
 ENDDO
 call mem_dealloc(TEMP)
 call mem_dealloc(TEMP2)
END SUBROUTINE buildSJ000_OP

!> \brief Gaussian geminal overlap 
!> \author W. Klopper and S. Reine
!> \date 2010-10-13
!> \param GJ000 \f$ G^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) !> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param Coeff Gaussian geminal expansion coefficients
!> \param GGemexponent Gaussian geminal exponents
!> \param n Number of Gaussian geminals in expansion of Slater Geminal
SUBROUTINE buildGJ000(GJ000,nPrim,Jmax,Prefactor,Alpha,R2,Coeff,GGemexponent,n)
use memory_handling
 IMPLICIT NONE
 INTEGER,intent(in)        :: nPrim,Jmax,n
 real(realk),intent(inout) :: GJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)    :: Prefactor(nPrim),ALPHA(nPrim),R2(nPrim)
 REAL(REALK),intent(in)    :: Coeff(n),GGemexponent(n)
 !
 Real(realk), parameter :: PI=3.14159265358979323846E0_realk, OneOverTwoPi=0.159154943091895335768883763372514E0_realk
 Real(realk), parameter :: One = 1.0E0_realk, MinusTwo = -2.0E0_realk
 REAL(REALK),pointer    :: TEMP1(:),TEMP2(:),TEMP3(:),RhoTilde(:),Ooor(:)
 REAL(REALK)            :: DumFact
 INTEGER                :: I,J,K

 call mem_alloc(TEMP1,nPrim)
 call mem_alloc(TEMP2,nPrim)
 call mem_alloc(TEMP3,nPrim)
 call mem_alloc(RhoTilde,nPrim)
 call mem_alloc(Ooor,nPrim)
 DumFact = Coeff(1)*OneOverTwoPi
 DO I=1,NPrim
    Ooor(I) = One/(Alpha(I)+GGemexponent(1))
    RhoTilde(I) = GGemexponent(1)*Ooor(I) 
    TEMP3(I) = Alpha(I)*RhoTilde(I)
    TEMP1(I) = -R2(I)*TEMP3(I)
    TEMP3(I) = TEMP3(I)*MinusTwo
    Ooor(I) = PI*Ooor(I)
    TEMP2(I) = Ooor(I) * SQRT(Ooor(I)) * EXP(TEMP1(I))
    GJ000(0,I) = DumFact*TEMP2(I)
 ENDDO
 DO J=1,Jmax
    DO I=1,NPrim
       TEMP2(I) = TEMP2(I)*TEMP3(I)
       GJ000(J,I) =DumFact*TEMP2(I)
    ENDDO
 ENDDO
 DO K=2,n
    DumFact = Coeff(K)*OneOverTwoPi
    DO I=1,NPrim
       Ooor(I) = One/(Alpha(I)+GGemexponent(K))
       RhoTilde(I) = GGemexponent(K)*Ooor(I) 
       TEMP3(I) = Alpha(I)*RhoTilde(I)
       TEMP1(I) = -R2(I)*TEMP3(I)
       TEMP3(I) = TEMP3(I)*MinusTwo
       Ooor(I) = PI*Ooor(I)
       TEMP2(I) = Ooor(I) * SQRT(Ooor(I)) * EXP(TEMP1(I))
       GJ000(0,I) = GJ000(0,I) + DumFact*TEMP2(I)
    END DO
    DO J=1,Jmax
       DO I=1,NPrim
          TEMP2(I) = TEMP2(I)*TEMP3(I)
          GJ000(J,I) = GJ000(J,I) + DumFact*TEMP2(I)
       ENDDO
    ENDDO
 ENDDO
 DO I=1,NPrim
    GJ000(:,I) = GJ000(:,I)*Alpha(I)*Prefactor(I)
 ENDDO
 call mem_dealloc(TEMP1)
 call mem_dealloc(TEMP2)
 call mem_dealloc(TEMP3)
 call mem_dealloc(RhoTilde)
 call mem_dealloc(Ooor)
END SUBROUTINE buildGJ000

!> \brief Gaussian geminal integral over double commutator with kinetic energy
!> \author W. Klopper and S. Reine
!> \date 2010-10-13
!> \param GJ000 \f$ G^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param Coeff Gaussian geminal expansion coefficients
!> \param GGemexponent Gaussian geminal exponents
!> \param GGemexpoproduct product of two Gaussian geminal exponents
!> \param n Number of Gaussian geminals in expansion of Slater Geminal
SUBROUTINE buildGJ000Grad(GJ000,nPrim,Jmax,Prefactor,Alpha,R2,Coeff,GGemexponent,GGemexpoproduct,n)
use memory_handling
 IMPLICIT NONE
 INTEGER,intent(in)        :: nPrim,Jmax,n
 real(realk),intent(inout) :: GJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)    :: Prefactor(nPrim),ALPHA(nPrim),R2(nPrim)
 REAL(REALK),intent(in)    :: Coeff(n),GGemexponent(n),GGemexpoproduct(n)
 !
 Real(realk), parameter :: TwoSqrtPi=3.54490770181103205459633496668229E0_realk
 Real(realk), parameter :: Anderthalb = 1.5E0_realk, One = 1.0E0_realk, MinusTwo = -2.0E0_realk
 REAL(REALK),pointer    :: TEMP1(:),TEMP2(:),TEMP3(:),RhoHat(:),RhoTilde(:),Ooor(:)
 REAL(REALK)            :: DumFact,InvGGemexponent
 INTEGER                :: I,J,K
 call mem_alloc(TEMP1,nPrim)
 call mem_alloc(TEMP2,nPrim)
 call mem_alloc(TEMP3,nPrim)
 call mem_alloc(RhoHat,nPrim)
 call mem_alloc(RhoTilde,nPrim)
 call mem_alloc(Ooor,nPrim)
 DumFact = TwoSqrtPi*Coeff(1)*GGemexpoproduct(1)
 InvGGemexponent = One/GGemexponent(1)
 DO I=1,NPrim
    Ooor(I) = One/(Alpha(I)+GGemexponent(1))
    RhoHat(I) = Alpha(I)*Ooor(I)
    RhoTilde(I) = GGemexponent(1)*Ooor(I)
    TEMP3(I) = Alpha(I)*RhoTilde(I)
    TEMP1(I) = -R2(I)*TEMP3(I) 
    TEMP3(I) = TEMP3(I)*MinusTwo
    TEMP2(I) = Ooor(I)**2 * SQRT(Ooor(I)) * EXP(TEMP1(I))
    Ooor(I) = Alpha(I)*R2(I)
    TEMP1(I) = Anderthalb + RhoHat(I)*Ooor(I)
    GJ000(0,I) = DumFact*TEMP1(I)*TEMP2(I)
 ENDDO
 DO J=1,Jmax
    DO I=1,NPrim
       TEMP1(I) = TEMP1(I) - Alpha(I)*InvGGemexponent
       TEMP2(I) = TEMP2(I)*TEMP3(I)
       GJ000(J,I) = DumFact*TEMP1(I)*TEMP2(I)
    ENDDO
 ENDDO
 DO K=2,n
    DumFact = TwoSqrtPi*Coeff(K)*GGemexpoproduct(K)
    InvGGemexponent = One/GGemexponent(K)
    DO I=1,NPrim
       Ooor(I) = One/(Alpha(I)+GGemexponent(K))
       RhoHat(I) = Alpha(I)*Ooor(I)
       RhoTilde(I) = GGemexponent(K)*Ooor(I)
       TEMP3(I) = Alpha(I)*RhoTilde(I)
       TEMP1(I) = -R2(I)*TEMP3(I) 
       TEMP3(I) = TEMP3(I)*MinusTwo
       TEMP2(I) = Ooor(I)**2 * SQRT(Ooor(I)) * EXP(TEMP1(I))
       Ooor(I) = Alpha(I)*R2(I)
       TEMP1(I) = Anderthalb + RhoHat(I)*Ooor(I)
       GJ000(0,I) = GJ000(0,I) + DumFact*TEMP1(I)*TEMP2(I)
    END DO
    DO J=1,Jmax
       DO I=1,NPrim
          TEMP1(I) = TEMP1(I) - Alpha(I)*InvGGemexponent
          TEMP2(I) = TEMP2(I)*TEMP3(I)
          GJ000(J,I) = GJ000(J,I) + DumFact*TEMP1(I)*TEMP2(I)
       ENDDO
    ENDDO
 ENDDO
 DO I=1,NPrim
    GJ000(:,I) = GJ000(:,I)*Alpha(I)*Prefactor(I)
 ENDDO
 call mem_dealloc(TEMP1)
 call mem_dealloc(TEMP2)
 call mem_dealloc(TEMP3)
 call mem_dealloc(RhoHat)
 call mem_dealloc(RhoTilde)
 call mem_dealloc(Ooor)
END SUBROUTINE buildGJ000Grad

!> \brief Gaussian geminal Coulomb integrals
!> \author W. Klopper and S. Reine
!> \date 2010-10-13
!> \param GJ000 \f$ G^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param Coeff Gaussian geminal expansion coefficients
!> \param GGemexponent Gaussian geminal exponents
!> \param n Number of Gaussian geminals in expansion of Slater Geminal
!> \param integral containing tabulated boys function \f$ F_{n} \f$ 
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildGJ000Coulomb(GJ000,nPrim,Jmax,Prefactor,Alpha,R2,Coeff,GGemexponent,n,integral,HIGH,LUPRI,IPRINT)
use memory_handling
 IMPLICIT NONE
 INTEGER,intent(in)        :: nPrim,Jmax,n,LUPRI,IPRINT
 real(realk),intent(inout) :: GJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)    :: Prefactor(nPrim),ALPHA(nPrim),R2(nPrim)
 REAL(REALK),intent(in)    :: Coeff(n),GGemexponent(n)
 TYPE(integralitem),intent(in) :: integral
 LOGICAL,intent(in)  :: HIGH
 !
 Real(realk), parameter :: PI=3.14159265358979323846E0_realk
 Real(realk), parameter :: One = 1.0E0_realk, MinusTwo = -2E0_realk, MinusHalf = -0.5E0_realk
 REAL(REALK),pointer    :: TEMP1(:),TEMP2(:),TEMP3(:),RhoTilde(:),BinomArray(:,:),Ooor(:),RhoHat(:),RJ000(:,:)
 INTEGER                :: I,J,K,M
 call mem_alloc(TEMP1,nPrim)
 call mem_alloc(TEMP2,nPrim)
 call mem_alloc(TEMP3,nPrim)
 call mem_alloc(RhoTilde,nPrim)
 call mem_alloc(RhoHat,nPrim)
 call mem_alloc(Ooor,nPrim)
 call mem_alloc(BinomArray,Jmax,Jmax,.TRUE.,.TRUE.)
 call mem_alloc(RJ000,Jmax,nPrim,.TRUE.,.FALSE.)
 !
 call binominit(BinomArray,Jmax)

 !
 DO I=1,NPrim
    Ooor(I) = One/(Alpha(I)+GGemexponent(1))
    RhoHat(I) = Alpha(I)*Ooor(I)
    RhoTilde(I) = GGemexponent(1)*Ooor(I)
    TEMP3(I) = Alpha(I)*RhoTilde(I)
    TEMP1(I) = -R2(I)*TEMP3(I)
    TEMP3(I) = MinusHalf/TEMP3(I)
    TEMP2(I) = Ooor(I)*EXP(TEMP1(I))
    TEMP1(I) = Alpha(I)*RhoHat(I)
 ENDDO
 IF (HIGH) THEN
    CALL buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &TEMP1,R2,Prefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ELSE
    CALL buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &TEMP1,R2,Prefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ENDIF
 DO I=1,NPrim
    GJ000(0,I) = Coeff(1)*RJ000(0,I)*TEMP2(I)
 ENDDO
 DO J=1,Jmax
    DO I=1,NPrim
       Ooor(I) = RhoTilde(I)**J
       TEMP1(I) = RJ000(0,I)*Ooor(I) 
       DO M=1,J
          Ooor(I) = Ooor(I)*TEMP3(I)
          TEMP1(I) = TEMP1(I) + BinomArray(M,J)*RJ000(M,I)*Ooor(I)
       ENDDO
       TEMP2(I) = TEMP2(I)*Alpha(I)*MinusTwo
       GJ000(J,I) = Coeff(1)*TEMP1(I)*TEMP2(I)
    ENDDO
 ENDDO
 DO K=2,n
    DO I=1,NPrim
       Ooor(I) = One/(Alpha(I)+GGemexponent(K))
       RhoHat(I) = Alpha(I)*Ooor(I)
       RhoTilde(I) = GGemexponent(K)*Ooor(I)
       TEMP3(I) = Alpha(I)*RhoTilde(I)
       TEMP1(I) = -R2(I)*TEMP3(I)
       TEMP3(I) = MinusHalf/TEMP3(I)
       TEMP2(I) = Ooor(I)*EXP(TEMP1(I))
       TEMP1(I) = Alpha(I)*RhoHat(I)
    ENDDO
    IF (HIGH) THEN
       CALL buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
            &TEMP1,R2,Prefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
    ELSE
       CALL buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
            &TEMP1,R2,Prefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
    ENDIF
    DO I=1,NPrim
       GJ000(0,I) = GJ000(0,I) + Coeff(K)*RJ000(0,I)*TEMP2(I)
    ENDDO
    DO J=1,Jmax
       DO I=1,NPrim
          Ooor(I) = RhoTilde(I)**J
          TEMP1(I) = RJ000(0,I)*Ooor(I)
          DO M=1,J
             Ooor(I) = Ooor(I)*TEMP3(I)
             TEMP1(I) = TEMP1(I) + BinomArray(M,J)*RJ000(M,I)*Ooor(I)
          ENDDO
          TEMP2(I) = TEMP2(I)*Alpha(I)*MinusTwo
          GJ000(J,I) = GJ000(J,I) + Coeff(K)*TEMP1(I)*TEMP2(I)
       ENDDO
    ENDDO
 ENDDO
 DO I=1,NPrim
    GJ000(:,I) = GJ000(:,I)*Alpha(I)
 ENDDO
 call mem_dealloc(TEMP1)
 call mem_dealloc(TEMP2)
 call mem_dealloc(TEMP3)
 call mem_dealloc(RhoTilde)
 call mem_dealloc(RhoHat)
 call mem_dealloc(Ooor)
 call mem_dealloc(BinomArray)
 call mem_dealloc(RJ000)
END SUBROUTINE buildGJ000Coulomb

!> \brief wrapper to 2 electron Coulomb hermite integral \f$ R^{n}_{000}(\alpha,\textbf{R}_{PQ}) =  \left(-2 \alpha \right)^{n}  F_{n}(\alpha R_{PQ}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ R^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param integral containing tabulated boys function \f$ F_{n} \f$ 
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildRJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,HIGH)
 IMPLICIT NONE
 INTEGER,intent(in)         :: nPrim,LUPRI,IPRINT,Jmax
 REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrim)
 TYPE(Integrand),intent(in) :: PQ
 TYPE(integralitem),intent(in) :: integral
 LOGICAL,intent(in)  :: HIGH

 IF (HIGH) THEN
    CALL buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &PQ%reducedExponents,PQ%squaredDistance,&
         &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ELSE
    CALL buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &PQ%reducedExponents,PQ%squaredDistance,&
         &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ENDIF

END SUBROUTINE buildRJ000

!> \brief wrapper to 1 electron hermite Coulomb integral \f$ R^{n}_{000}(p,\textbf{R}_{PC}) =  \left(-2 p \right)^{n}  F_{n}(p R_{PC}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ R^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param integral containing tabulated boys function \f$ F_{n} \f$ 
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildNuclearRJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,HIGH)
 IMPLICIT NONE
 INTEGER,intent(in)        :: nPrim,LUPRI,IPRINT,Jmax
 REAL(REALK),intent(inout) :: RJ000(0:Jmax,nPrim)
 TYPE(Integrand),intent(in):: PQ
 TYPE(integralitem),intent(in) :: integral
 LOGICAL,intent(in) :: HIGH

 IF (HIGH) THEN
    CALL buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &PQ%Exponents,PQ%squaredDistance,&
         &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ELSE
    CALL buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &PQ%Exponents,PQ%squaredDistance,&
         &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ENDIF

END SUBROUTINE buildNuclearRJ000

!> \brief Low accuracy version of electron hermite Coulomb integral \f$ R^{n}_{000}(\alpha,\textbf{R}_{PQ}) =  \left(-2 \alpha \right)^{n}  F_{n}(\alpha R_{PQ}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine, B. Jansik  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ R^{0:jmax}_{000} \f$
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param alpha reduced exponents
!> \param squared distance
!> \param Prefactor (Coulomb: \f$ \frac{2 \pi^{\frac{5}{2}}}{pq\sqrt{p+q}} \f$, Nuclear att: \f$ \frac{2 \pi}{p}\f$)
!> \param TABFJW tabulated boys function \f$ F_{n} \f$ 
!> \param nTABFJW1 dimension 1 of TABFJW
!> \param nTABFJW2 dimension 2 of TABFJW
SUBROUTINE buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
 IMPLICIT NONE
 INTEGER,intent(in)        :: nPrim,Jmax,Lupri,Iprint,nTABFJW1,nTABFJW2
 REAL(REALK),intent(inout) :: RJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)    :: alpha(nPrim),R2(nprim),Prefactor(nprim)
 REAL(REALK),intent(in)    :: TABFJW(0:nTABFJW1,0:nTABFJW2)
!
 INTEGER         :: I
 REAL(REALK)     :: D2JP36,WVAL!,WVALS(Nprim,3),WVALU(Nprim)
 REAL(REALK),PARAMETER :: HALF =0.5E0_realk,D1=1E0_realk,D2 = 2E0_realk, D4 = 4E0_realk, D10=10E0_realk,D100=100E0_realk
 Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk
 REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
 REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk
 Integer :: IPNT,J
 Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
 REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
 REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
 REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
 REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
 Real(realk), parameter :: PI=3.14159265358979323846E0_realk
 REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
 REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
 REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
 Real(realk) :: R!W2,W3,W4,W5,W6
 REAL(REALK), PARAMETER :: SMALL = 1E-15_realk!0.000001E0_realk

 D2JP36 = 2*JMAX + 36

 DO I = 1, Nprim
    WVAL = alpha(I)*R2(I)
    !  0 < WVAL < 0.000001
    IF (WVAL .LT. SMALL) THEN         
       RJ000(0,I) = D1
       DO J=1,JMAX
          RJ000(J,I)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
       ENDDO
       !  0 < WVAL < 12 
    ELSE IF (WVAL .LT. D12) THEN
       IPNT = NINT(D100*WVAL)
       WDIFF = WVAL - TENTH*IPNT
       !    W2    = WDIFF*WDIFF
       !    W3    = W2*WDIFF
       !    W4    = W3*WDIFF
       !    W5    = W4*WDIFF
       !    W6    = W5*WDIFF

       !    W2    = W2*COEF2
       !    W3    = W3*COEF3
       !    W4    = W4*COEF4
       !    W5    = W5*COEF5
       !    W6    = W6*COEF6

       DO J=0,JMAX
          R = TABFJW(J,IPNT)
          R = R -TABFJW(J+1,IPNT)*WDIFF
          !    R = R + TABFJW(J+2,IPNT)*W2
          !    R = R + TABFJW(J+3,IPNT)*W3
          !    R = R + TABFJW(J+4,IPNT)*W4
          !    R = R + TABFJW(J+5,IPNT)*W5
          !    R = R + TABFJW(J+6,IPNT)*W6
          RJ000(J,I) = R
       ENDDO

       !  12 < WVAL <= (2J+36) 
    ELSE IF (WVAL.LE.D2JP36) THEN
       REXPW = HALF*EXP(-WVAL)
       RWVAL = D1/WVAL
       GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
       RJ000(0,I) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
       DO J=1,JMAX
          RJ000(J,I) = RWVAL*((J - HALF)*RJ000(J-1,I)-REXPW)
       ENDDO
       !  (2J+36) < WVAL 
    ELSE
       RWVAL = PID4/WVAL
       RJ000(0,I) = SQRT(RWVAL)
       RWVAL = RWVAL*PID4I
       DO J = 1, JMAX
          RJ000(J,I) = RWVAL*(J - HALF)*RJ000(J-1,I)
       ENDDO
    END IF
 ENDDO

 ! Scaling
 DO I=1,nPrim
    PREF = Prefactor(I)
    RJ000(0,I) = PREF*RJ000(0,I)
    IF (jmax.GT. 0) THEN
       D2MALPHA = -2*alpha(I)
       DO j=1,jmax
          PREF = PREF*D2MALPHA
          RJ000(J,I) = PREF*RJ000(J,I)
       ENDDO
    ENDIF
 ENDDO

END SUBROUTINE buildRJ000_OP_LA

!> \brief High accuracy version of electron hermite Coulomb integral \f$ R^{n}_{000}(\alpha,\textbf{R}_{PQ}) =  \left(-2 \alpha \right)^{n}  F_{n}(\alpha R_{PQ}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine, B. Jansik  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ R^{0:jmax}_{000} \f$
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param alpha reduced exponents
!> \param squared distance
!> \param Prefactor (Coulomb: \f$ \frac{2 \pi^{\frac{5}{2}}}{pq\sqrt{p+q}} \f$, Nuclear att: \f$ \frac{2 \pi}{p}\f$)
!> \param TABFJW tabulated boys function \f$ F_{n} \f$ 
!> \param nTABFJW1 dimension 1 of TABFJW
!> \param nTABFJW2 dimension 2 of TABFJW
SUBROUTINE buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
 IMPLICIT NONE
 INTEGER,intent(in)         :: nPrim,Jmax,Lupri,Iprint,nTABFJW1,nTABFJW2
 REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)     :: alpha(nPrim),R2(nprim),Prefactor(nprim)
 REAL(REALK),intent(in)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)
!
 INTEGER         :: I,I2
 REAL(REALK)     :: D2JP36,WVAL!,WVALS(Nprim,3),WVALU(Nprim)
 REAL(REALK),PARAMETER :: HALF =0.5E0_realk,D1=1E0_realk,D2 = 2E0_realk, D4 = 4E0_realk, D10=10E0_realk,D100=100E0_realk
 Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk
 REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
 REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk
 Integer :: IPNT,J
 Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
 REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
 REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
 REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
 REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
 Real(realk), parameter :: PI=3.14159265358979323846E0_realk
 REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
 REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
 REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
 Real(realk) :: W2,W3,R!W4,W5,W6,
! LOGICAL :: maxJgt0
 REAL(REALK), PARAMETER :: SMALL = 1E-15_realk!0.000001E0_realk
 integer :: NVAL1(Nprim),NVAL2(Nprim),NVAL3(Nprim),NVAL4(Nprim),N1,N2,N3,N4
 real(realk) :: WVAL2(Nprim),WVAL3(Nprim),WVAL4(Nprim)
 N1 = 0
 N2 = 0
 N3 = 0
 N4 = 0
 D2JP36 = 2*JMAX + 36
 DO I = 1, Nprim
    WVAL = alpha(I)*R2(I)
    !  0 < WVAL < 0.000001
    IF (ABS(WVAL) .LT. SMALL) THEN         
       N1 = N1 + 1
       NVAL1(N1) = I
    !  0 < WVAL < 12 
    ELSE IF (WVAL .LT. D12) THEN
       N2 = N2 + 1
       NVAL2(N2) = I
       WVAL2(N2) = WVAL
    !  12 < WVAL <= (2J+36) 
    ELSE IF (WVAL.LE.D2JP36) THEN
       N3 = N3 + 1
       NVAL3(N3) = I
       WVAL3(N3) = WVAL
       !  (2J+36) < WVAL 
    ELSE
       N4 = N4 + 1
       NVAL4(N4) = I
       WVAL4(N4) = WVAL
    ENDIF
 ENDDO

 IF(N1.GT.0)THEN
    DO I = 1, N1
       I2 = NVAL1(I)
       RJ000(0,I2) = D1
       DO J=1,JMAX
          RJ000(J,I2)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
       ENDDO
    ENDDO
 ENDIF
 IF(N2.GT.0)THEN
    DO I = 1, N2
       I2 = NVAL2(I)
       IPNT = NINT(D100*WVAL2(I))
       WDIFF = WVAL2(I) - TENTH*IPNT
       W2    = WDIFF*WDIFF
       W3    = W2*WDIFF
       !    W4    = W3*WDIFF
       !    W5    = W4*WDIFF
       !    W6    = W5*WDIFF
       W2    = W2*COEF2
       W3    = W3*COEF3
       !    W4    = W4*COEF4
       !    W5    = W5*COEF5
       !    W6    = W6*COEF6
       DO J=0,JMAX
          R = TABFJW(J,IPNT)
          R = R -TABFJW(J+1,IPNT)*WDIFF
          R = R + TABFJW(J+2,IPNT)*W2
          R = R + TABFJW(J+3,IPNT)*W3
          !    R = R + TABFJW(J+4,IPNT)*W4
          !    R = R + TABFJW(J+5,IPNT)*W5
          !    R = R + TABFJW(J+6,IPNT)*W6
          RJ000(J,I2) = R
       ENDDO
    ENDDO
 ENDIF
 IF (N3.GT.0) THEN
    CALL ls_vdinv2(N3,WVAL3,WVAL2) !use WVAL2 as temp array to store D1/WVAL3
    DO I = 1, N3
       I2 = NVAL3(I)
       REXPW = HALF*EXP(-WVAL3(I))
       RWVAL = WVAL2(I)!D1/WVAL3(I)
       GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
       RJ000(0,I2) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
       DO J=1,JMAX
          RJ000(J,I2) = RWVAL*((J - HALF)*RJ000(J-1,I2)-REXPW)
       ENDDO
    ENDDO
 ENDIF
 IF(N4.GT.0)THEN
    CALL ls_vdinv2(N4,WVAL4,WVAL2) !use WVAL2 as temp array to store D1/WVAL3
    DO I = 1, N4
       I2 = NVAL4(I)
       RWVAL = PID4*WVAL2(I) !PID4/WVAL4(I)
       RJ000(0,I2) = SQRT(RWVAL)
!       RWVAL = RWVAL*PID4I = WVAL2(I)
       RWVAL = WVAL2(I)
       DO J = 1, JMAX
          RJ000(J,I2) = RWVAL*(J - HALF)*RJ000(J-1,I2)
       ENDDO
    ENDDO
 END IF

 ! Scaling
 DO I=1,nPrim
    PREF = Prefactor(I)
    RJ000(0,I) = PREF*RJ000(0,I)
    IF (jmax.GT. 0) THEN
       D2MALPHA = -2*alpha(I)
       DO j=1,jmax
          PREF = PREF*D2MALPHA
          RJ000(J,I) = PREF*RJ000(J,I)
       ENDDO
    ENDIF
 ENDDO
 CONTAINS
   subroutine ls_vdinv2(N,A,invA)
     implicit none
     integer,intent(in) :: N
     real(realk),intent(in) :: A(N)
     real(realk),intent(inout) :: invA(N)
     !
!#ifdef VAR_MKL
!     call vdinv(N,A,invA)
!#else
     integer :: i
     DO i=1, N
        invA(i) = 1.0E0_realk/A(i)
     ENDDO
!#endif
   end subroutine ls_vdinv2
   
END SUBROUTINE buildRJ000_OP_HA

!> \brief tabulation of the boys function or incomplete gamma function \f$ F_{n} \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param LUPRI logical unit number for output printing
!> \param JMX tabulate to this order of angular momentum
!> \param sharedTUV stores the tabulated values
!>
!> ***** Tabulation of incomplete gamma function *****
!>   For J = JMX a power series expansion is used, see for
!>   example Eq.(39) given by V. Saunders in "Computational
!>   Techniques in Quantum Chemistry and Molecular Physics",
!>   Reidel 1975.  For J < JMX the values are calculated
!>   using downward recursion in J.
!>
SUBROUTINE gammaTabulation(LUPRI,JMX,sharedTUV)
 use TYPEDEF
 IMPLICIT NONE
 TYPE(TUVitem),intent(inout) :: sharedTUV
 INTEGER,intent(in)       :: JMX,LUPRI
!
 INTEGER           :: MAXJ0,IADR,IPOINT,IORDER,JADR,JMAX,J
 REAL(REALK)       :: DENOM,D2MAX1,R2MAX1,TERM,SUM,REXPW,WVAL,D2WAL
 REAL(REALK), PARAMETER :: HALF = 0.5E0_realk,  TEN6 = 1.0E6_realk
 REAL(REALK), PARAMETER :: D1 = 1E0_realk, D10 = 10E0_realk
 REAL(REALK), PARAMETER :: D2 = 2E0_realk, D4 = 4E0_realk, D12 = 12E0_realk, TENTH = 0.01E0_realk
 REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
 REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk

 REAL(REALK), PARAMETER :: PI    = 3.14159265358979323846E00_realk
 REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
 REAL(REALK), PARAMETER :: R2PI52 = 5.91496717279561287782E00_realk
 REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
 REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI

 REAL(REALK), PARAMETER :: GFAC30 =  .4999489092E0_realk 
 REAL(REALK), PARAMETER :: GFAC31 = -.2473631686E0_realk
 REAL(REALK), PARAMETER :: GFAC32 =  .321180909E0_realk
 REAL(REALK), PARAMETER :: GFAC33 = -.3811559346E0_realk
 REAL(REALK), PARAMETER :: GFAC20 = .4998436875E0_realk
 REAL(REALK), PARAMETER :: GFAC21 = -.24249438E0_realk
 REAL(REALK), PARAMETER :: GFAC22 =  .24642845E0_realk
 REAL(REALK), PARAMETER :: GFAC10 =  .499093162E0_realk
 REAL(REALK), PARAMETER :: GFAC11 = -.2152832E0_realk
 REAL(REALK), PARAMETER :: GFAC00 =  .490E0_realk

 JMAX = JMX + 6
 MAXJ0 = JMAX
 !
 !     WVAL = 0.0
 !
 IADR = 1
 DENOM = D1
 DO J = 0,JMAX
    SHAREDTUV%TABFJW(J,0) = D1/DENOM
    IADR = IADR + 1201
    DENOM = DENOM + D2
 ENDDO
 !
 !     WVAL = 0.1, 0.2, 0.3,... 12.0
 !
 IADR = IADR - 1201
 !D2MAX1 = DFLOAT(2*JMAX + 1)
 D2MAX1 = 2*JMAX + 1
 R2MAX1 = D1/D2MAX1
 DO IPOINT = 1,1200
    !  WVAL = TENTH*DFLOAT(IPOINT)
    WVAL = TENTH*IPOINT
    D2WAL = WVAL + WVAL
    IADR = IADR + 1
    TERM = R2MAX1
    SUM = TERM
    DENOM = D2MAX1
    DO IORDER = 2,200
       DENOM = DENOM + D2
       TERM = TERM*D2WAL/DENOM
       SUM = SUM + TERM
       IF (TERM .LE. 1.0E-15_realk) EXIT
    ENDDO
    REXPW = EXP(-WVAL)
    SHAREDTUV%TABFJW(JMAX,IPOINT) = REXPW*SUM
    DENOM = D2MAX1
    JADR = IADR
    DO J = 1,JMAX
       DENOM = DENOM - D2
       SHAREDTUV%TABFJW(JMAX-J,IPOINT) = (SHAREDTUV%TABFJW(JMAX-J+1,IPOINT)*D2WAL + REXPW)/DENOM
       JADR = JADR - 1201
    ENDDO
 ENDDO
END SUBROUTINE gammaTabulation

!> \brief initialization routine to init the TUVitem and 
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param sharedTUV contain tabulated boys function and tuvindexing
!> \param Integral contain temporary arrays and sharedTUV
!> \param Input contain info about the requested integral 
!> \param OD_LHS contain the left hand side overlapdistribution 
!> \param OD_RHS contain the right hand side overlapdistribution 
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE initTUVitem(sharedTUV,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
 implicit none
 TYPE(TUVitem),intent(inout),target :: sharedTUV
 TYPE(IntegralInput),intent(in)     :: Input
 TYPE(ODITEM),intent(in)            :: OD_LHS,OD_RHS
 Integer,intent(in)                 :: LUPRI,IPRINT
 !
 LOGICAL              :: TABULATE_BOYS
 Integer,parameter :: MAXDER=4
! Integer,parameter :: MAXJ=4*MAXAOJ
 !Integer,parameter :: MAXAOJ=12,MAXDER=2
 !Integer,parameter :: MAXJ=12
 Integer,parameter :: NODES=1201,ORDER=7
! Integer,parameter :: NTABFJ000=NODES*(MAXJ+ORDER)
! Integer,parameter :: NTUVMAX = (MAXJ+1)*(MAXJ+2)*(MAXJ+3)/6
 Integer :: NTABFJ000
 Integer :: NTUVMAX,maxj
 !
 maxJ = OD_LHS%maxJ+OD_RHS%maxJ + Input%geoderivorder &
      & + input%MMorder + input%CMorder + Input%magderorderP+Input%magderorderQ
 IF(INPUT%Operator.EQ.KineticOperator)THEN
    maxJ = maxJ+2
 ENDIF
 IF(INPUT%PropRequireBoys.GT.-1)maxJ = maxJ+INPUT%PropRequireBoys
 NTABFJ000=NODES*(MAXJ+ORDER)
 NTUVMAX = (MAXJ+1)*(MAXJ+2)*(MAXJ+3)/6
 TABULATE_BOYS = .FALSE.
 IF(INPUT%PropRequireBoys.GT.-1) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. CoulombOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator.EQ.NucpotOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. CAMOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. ErfcOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. ErfOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. GGemCouOperator) TABULATE_BOYS = .TRUE.
 IF(TABULATE_BOYS)THEN
    !Pretabulate Gammafunction/Boysfunction 
    NULLIFY(sharedTUV%TABFJW)  !MXQN=13     MAXJ = 4*(MXQN - 1) + 2 
    call mem_alloc(SharedTUV%TABFJW,MAXJ+6,1200,.TRUE.,.TRUE.)     !121*(MAXJ + 7) = 6897
    SharedTUV%nTABFJW1=MAXJ+6
    SharedTUV%nTABFJW2=1200
    CALL GAMMATABULATION(lupri,MAXJ,sharedTUV)    !Jmax is set to 10 is that ok?
 ENDIF
 call mem_alloc(SharedTUV%TUVindex,MAXJ,MAXJ,MAXJ,.TRUE.,.TRUE.,.TRUE.)
 call mem_alloc(SHAREDTUV%Tindex,NTUVMAX)
 call mem_alloc(SHAREDTUV%Uindex,NTUVMAX)
 call mem_alloc(SHAREDTUV%Vindex,NTUVMAX)

 CALL integralTUVindex(SharedTUV,MAXJ,LUPRI,IPRINT)

 call mem_alloc(SHAREDTUV%SYMindex,NTUVMAX)
 CALL integralTUVSYMindex(SharedTUV,MAXJ,Input%iPQxyz,LUPRI,IPRINT)
 SHAREDTUV%iPQxyz = Input%iPQxyz
 CALL BuildPrecalculatedSphmat(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)
! SHAREDTUV%maxJ = maxJ
! SHAREDTUV%NTUVMAX = NTUVMAX
END SUBROUTINE initTUVitem

!> \brief build spherical transformation matrices and attach to sharedTUV 
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param sharedTUV to store the matrices
!> \param OD_LHS contain the left hand side overlapdistribution 
!> \param OD_RHS contain the right hand side overlapdistribution 
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE buildPrecalculatedSphmat(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)
 IMPLICIT NONE
 TYPE(ODITEM),intent(in)     :: OD_LHS,OD_RHS
 TYPE(TUVitem),intent(inout) :: SharedTUV
 INTEGER,intent(in)          :: LUPRI,IPRINT
 !
 INTEGER  :: I,MAXANGMOM

 MAXANGMOM = 0
 DO I=1,OD_LHS%nbatches
    MAXANGMOM = MAX(MAXANGMOM,OD_LHS%BATCH(I)%AO(1)%p%maxangmom)
    MAXANGMOM = MAX(MAXANGMOM,OD_LHS%BATCH(I)%AO(2)%p%maxangmom)
 ENDDO
 DO I=1,OD_RHS%nbatches
    MAXANGMOM = MAX(MAXANGMOM,OD_RHS%BATCH(I)%AO(1)%p%maxangmom)
    MAXANGMOM = MAX(MAXANGMOM,OD_RHS%BATCH(I)%AO(2)%p%maxangmom)
 ENDDO

 SharedTUV%nSPHMAT = MAXANGMOM + 3
 CALL setPrecalculatedSphmat(SharedTUV%SPH_MAT,SharedTUV%nSPHMAT,LUPRI,IPRINT)

END SUBROUTINE buildPrecalculatedSphmat

!> \brief set the TUVindex,Tindex,Uindex and Vindex  
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param sharedTUV to store the info
!> \param MAXJ the maximum angular momentum for this calculation
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE integralTUVindex(sharedTUV,MAXJ,LUPRI,IPRINT)
 implicit none
 TYPE(TUVitem),intent(inout)  :: sharedTUV
 Integer,intent(in)        :: LUPRI,IPRINT,MAXJ
!
 Integer :: TUV,T,U,V,J
!
 TUV=0
 DO J = 0, MAXJ
    DO T = J,0,-1       
       DO U = J-T,0,-1
          TUV=TUV+1
          V=J-T-U
          sharedTUV%TUVindex(T,U,V)=TUV
          sharedTUV%Tindex(TUV)=T
          sharedTUV%Uindex(TUV)=U
          sharedTUV%Vindex(TUV)=V
       ENDDO
    ENDDO
 ENDDO
END SUBROUTINE integralTUVindex

!> \brief set the TUVSYMindex
!> \author \latexonly T. Kj{\ae}rgaard \endlatexonly
!> \date 2011
!> \param sharedTUV to store the info
!> \param MAXJ the maximum angular momentum for this calculation
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE IntegralTUVSYMindex(SharedTUV,MAXJ,IPQxyz,LUPRI,IPRINT)
 implicit none
 TYPE(TUVitem),intent(inout)  :: sharedTUV
 Integer,intent(in)        :: LUPRI,IPRINT,MAXJ,IPQxyz
!
 Integer :: TUV,T,U,V,J,IBTUV
!
 TUV=0
 DO J = 0, MAXJ
    DO T = J,0,-1       
       DO U = J-T,0,-1
          TUV=TUV+1
          V=J-T-U
          IBTUV = MOD(T,2) + 2*MOD(U,2) + 4*MOD(V,2)
          sharedTUV%SYMindex(TUV)=IAND(IBTUV,IPQxyz)
!          print*,'SYM(TUV=',TUV,',ISYM=',IPQxyz,')=',IAND(IBTUV,IPQxyz)
       ENDDO
    ENDDO
 ENDDO

END SUBROUTINE INTEGRALTUVSYMINDEX

!> \brief dealloc memory for the TUVitem
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param sharedTUV to store the info
!> \param Input contain info about the requested integral 
SUBROUTINE freeTUVitem(sharedTUV,Input)
 implicit none
 TYPE(TUVItem),intent(inout)  :: sharedTUV
 TYPE(IntegralInput),intent(in) :: Input
!
 LOGICAL              :: TABULATE_BOYS
!
 TABULATE_BOYS = .FALSE.
 IF(INPUT%PropRequireBoys.GT.-1) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. CoulombOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator.EQ.NucpotOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. CAMOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. ErfcOperator) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. ErfOperator ) TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator .EQ. GGemCouOperator) TABULATE_BOYS = .TRUE.

 IF(TABULATE_BOYS) CALL MEM_DEALLOC(sharedTUV%TABFJW)
 CALL MEM_DEALLOC(SharedTUV%TUVindex)
 CALL MEM_DEALLOC(SHAREDTUV%Tindex)
 CALL MEM_DEALLOC(SHAREDTUV%Uindex)
 CALL MEM_DEALLOC(SHAREDTUV%Vindex)
 CALL freePrecalculatedSphmat(SharedTUV%SPH_MAT,SharedTUV%nSPHMAT)
 CALL MEM_DEALLOC(SHAREDTUV%SYMindex)

END SUBROUTINE freeTUVitem

!> \brief multipole moment hermite integrals 
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param INTEGRAL store temporary arrays
!> \param Input contain info about the requested integral 
!> \param PQ Integrand containing info about overlapdistribution P and Q
!> \param nPrim number of primitive functions
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param CARTORDER cartesian order of the multipole moments
!>
!> The final MOM(T,U,V) integrals are arranged as follows:
!> S(000)
!> S(100) S(010) S(001)
!> S(200) S(110) S(101) S(020) S(011) S(002)
!> S(300) S(210) S(201) S(120) S(111) S(102) S(030) S(021) S(012) S(003)
!>
SUBROUTINE getMomtuv(INTEGRAL,INPUT,PQ,nprim,LUPRI,IPRINT,CARTORDER)
 implicit none
 TYPE(Integrand),intent(in)         :: PQ
 TYPE(Integralitem),intent(inout)   :: integral
 TYPE(IntegralInput),intent(in)     :: Input
 INTEGER,intent(in)                 :: LUPRI,nprim,IPRINT,CARTORDER
!
 INTEGER                 :: J,T,U,V
 INTEGER                 :: ntuv,I,idir,ituv,tuvstart
 INTEGER                 :: Jmax,Jstart,e,f,g,efg,nefg
 INTEGER                 :: Xorder,Yorder,Zorder,C
 real(realk)             :: DISTANCE(nprim,3)
 Real(realk), parameter  :: OneHalf=1.5E0_realk
 Real(realk), parameter  :: PI=3.14159265358979323846E0_realk
 Real(realk), parameter  :: thresh=1E-14_realk

 Xorder=CARTORDER
 Yorder=CARTORDER
 Zorder=CARTORDER

! NPrim=PQ%nPrimitives
 JMAX=PQ%endAngmom!+CARTORDER
 Jstart=PQ%startAngmom !always zero
 Jstart=PQ%startAngmom  !always zero;  ! not for DF !! \Andreas 
!ANDREAS: added JSTART terms
 ntuv=(JMAX+1)*(JMAX+2)*(JMAX+3)/6-Jstart*(jstart+1)*(Jstart+2)/6
!ANDREAS
 nEFG=(CARTORDER+1)*(CARTORDER+2)*(CARTORDER+3)/6
 INTEGRAL%nTUV=ntuv
 INTEGRAL%nEFG=nefg

 IF(INPUT%DO_MULMOM .AND. INPUT%nMultipoleMomentComp .NE. nEFG) THEN
    WRITE(*,*)'ATTENTION in GETMOMTUV: NMMCOMP .NE. nEFG', INPUT%nMultipoleMomentComp,nEFG
    CALL LSQUIT('in GETMOMTUV: NMMCOMP .NE. nEFG',-1)
 END IF

 IF (JMAX+CARTORDER .EQ. 0) THEN
    !Special case JMAX = 0 A SIMPEL OVERLAP MATRIX WITH ONLY S ORBITALS
    !WRITE(LUPRI,*)'SPECIAL CASE OF CARTISIAN MULTIPOLE MOMENTS - OVERLAP'
    DO I=1,NPrim
       Integral%Rtuv(I)=(PI/PQ%P%p%exponents(I))**OneHalf
    ENDDO
 ELSE
    IF(CARTORDER .EQ. 0) THEN
!!$      !SIMPEL OVERLAP
!!$      NULLIFY(INTEGRAL%Wtuv)
!!$      call mem_alloc(INTEGRAL%Wtuv,Nprim,nTUV)
!!$      DO J = Jstart, PQ%endAngmom
!!$         DO T = J,0,-1
!!$            DO U = J-T,0,-1
!!$               V=J-T-U
!!$               DO I=1,nPrim
!!$                  INTEGRAL%Wtuv(I,INTEGRAL%TUVINDEX(T,U,V))=INTTEMP2(I,INTEGRAL%TUVindex(T,U,V))
!!$               ENDDO
!!$            ENDDO
!!$         ENDDO
!!$      ENDDO
    ELSE
       IF(INPUT%OD_MOM)THEN
          DO I = 1,Nprim
             DISTANCE(I,1) = PQ%P%p%center(1+(I-1)*3)-PQ%P%p%ODcenter(1)
             DISTANCE(I,2) = PQ%P%p%center(2+(I-1)*3)-PQ%P%p%ODcenter(2)
             DISTANCE(I,3) = PQ%P%p%center(3+(I-1)*3)-PQ%P%p%ODcenter(3)
             !        WRITE(LUPRI,'(2X,A,3ES16.9)')'DISTANCE =',DISTANCE(I,1),DISTANCE(I,2),DISTANCE(I,3)
          ENDDO
       ELSE
          DO I = 1,Nprim
             DISTANCE(I,1) = PQ%P%p%center(1+(I-1)*3)-INPUT%MOM_CENTER(1)
             DISTANCE(I,2) = PQ%P%p%center(2+(I-1)*3)-INPUT%MOM_CENTER(2)
             DISTANCE(I,3) = PQ%P%p%center(3+(I-1)*3)-INPUT%MOM_CENTER(3)
             !        WRITE(LUPRI,'(2X,A,3ES16.9)')'DISTANCE =',DISTANCE(I,1),DISTANCE(I,2),DISTANCE(I,3)
          ENDDO
       ENDIF
       DO idir=1,3
         DO I=1,Nprim
           IF (abs(DISTANCE(I,idir)).LE.thresh) DISTANCE(I,idir) = 0E0_realk
         ENDDO
       ENDDO
       !      CALL LS_DZERO(INTEGRAL%Rtuv,nTUV*nPrim*nEFG)
       CALL MultipoleRecurrence(INTEGRAL%Rtuv,JMAX,JSTART,CARTORDER,&
            & distance(:,1),distance(:,2),distance(:,3),&
            & Xorder,Yorder,Zorder,NPrim,nTUV,nEFG,&
            & PQ%P%p%Exponents,lupri)
       !   ENDIF
    ENDIF
 ENDIF
 !    Print section
 !    =============
 IF (IPRINT .GE. 10) THEN
    CALL LSHEADER(LUPRI,'Output from ERITUV')
    WRITE (LUPRI,'(2X,A13,I10)') 'MAX angmom ', JMAX
    WRITE (LUPRI,'(2X,A13,I10)') 'NPrim  ', NPrim
    WRITE (LUPRI,'(2X,A13,I10)') 'NTUV  ', nTUV
    WRITE (LUPRI,'(2X,A13,I10)') 'NEFG  ', nEFG
    IF (IPRINT .GE. 20) THEN
       CALL LSHEADER(LUPRI,'Hermite integrals M(t,u,v)')
       efg=0
       DO C=0,CARTORDER
          DO e = C,0,-1
             DO f = C-e,0,-1
                g=C-e-f
                efg=efg+1
                WRITE(LUPRI,*)'MULTIPOLE MOMENT LEVEL = (',e,',',f,',',g,')   EFG=',efg
                DO J = JSTART, JMAX
                   DO T = J,0,-1
                      DO U = J-T,0,-1
                         V=J-T-U
                         ITUV     = INTEGRAL%TUV%TUVINDEX(T,U,V)-Jstart*(Jstart+1)*(Jstart+2)/6
                         TUVSTART = (ITUV-1)*nPrim + (efg-1)*nPrim*nTUV
                         WRITE(LUPRI,*)'MULTIPOLE MOMENT LEVEL = (',T,',',U,',',V,')   TUV=',&
                              &INTEGRAL%TUV%TUVINDEX(T,U,V)-JSTART*(JSTART+1)*(JSTART+2)/6
                         WRITE (LUPRI,'(2X,A7,I1,A1,I1,A1,I1,A1,2X,5F28.17/,(12X,5F28.17))')&
                              & 'CARMOM(',T,',',U,',',V,')', &
                              &(INTEGRAL%Rtuv(I + TUVSTART),I=1,NPrim)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    END IF
 END IF
 Integral%nPrim=nPrim

END SUBROUTINE getMomtuv

!> \brief multipole moment hermite integral recurrence relations \f$ M^{e+1}_{t} = t M^{e}_{t-1}+X_{PC}M^{e}_{t}+\frac{1}{2p}M^{e}_{t+1} \f$ 
!> \author \latexonly T. Kj{\ae}rgaard \endlatexonly
!> \date 2009-02-05
!> \param OUTPUT of the recurrence relations
!> \param Jstart start angular momentum
!> \param JMAX maximum angular momentum
!> \param CARTORDER cartesian order of the multipole moments
!> \param Xpq distance in X axis
!> \param Ypq distance in Y axis
!> \param Zpq distance in Z axis
!> \param Xorder multipole moment expansion order in X
!> \param Yorder multipole moment expansion order in Y
!> \param Zorder multipole moment expansion order in Z
!> \param nPrim number of primitive functions
!> \param nTUV number of TUV components
!> \param nEFG number of cartesian moment components
!> \param TUVindex index array 
!> \param EXPONENTS
!> \param LUPRI logical unit number for output printing
!>
!> Se section 9.5.3 in the Book (page 356)
!> \f[ M^{efg}_{tuv} = M^{e}_{t}  M^{f}_{u}  M^{g}_{v} \f]
!> \f[ M^{e+1}_{t} = t M^{e}_{t-1}+X_{PC}M^{e}_{t}+\frac{1}{2p}M^{e}_{t+1} \f] 
!> starting with \f$ M^{0}_{000} = \sqrt{\frac{\pi}{p}} \f$
!>
SUBROUTINE multipoleRecurrence(OUTPUT,JMAX,JSTART,CARTORDER,&
    & Xpq,Ypq,Zpq,Xorder,Yorder,Zorder,NPrim,nTUV,nEFG,&
    & EXPONENTS,lupri)
 IMPLICIT NONE
 INTEGER,intent(in)       :: Xorder,Yorder,Zorder,nTUV,nEFG,Jmax,nprim,JSTART
INTEGER,intent(in)       :: CARTORDER
 real(realk),intent(inout):: OUTPUT(NTUV*nPrim*nEFG)
 real(realk),intent(in)   :: Xpq(nPrim),Ypq(nPrim),Zpq(nPrim),EXPONENTS(nPrim)
!
 INTEGER                 :: J,T,U,V,TUV
 real(realk)               :: P(nPrim),M000(nPrim)
 INTEGER                 :: I,e,f,g,efg
 INTEGER                 :: lupri,C,offset
 real(realk)             :: MX(nprim,0:Xorder,0:Xorder)
 real(realk)             :: MY(nprim,0:Yorder,0:Yorder)
 real(realk)             :: MZ(nprim,0:Zorder,0:Zorder)
 Real(realk), parameter  :: PI=3.14159265358979323846E0_realk
!#ifdef VAR_MKL
! real(realk)             :: TEMP(nPrim)
!#endif

 !CALL LSHEADER(LUPRI,'Multipole recurrence')
 DO I = 1, Nprim
    P(I)=1/(2E0_realk*EXPONENTS(I))   
 ENDDO
!#ifdef VAR_MKL
! DO I=1,NPrim
!    TEMP(I)=PI/EXPONENTS(I)
! ENDDO
! call vdsqrt(nprim,TEMP,M000)
!#else
 DO I=1,NPrim
    M000(I)=sqrt(PI/EXPONENTS(I))
 ENDDO
!#endif

 CALL LS_DZERO(MX,nPrim*(Xorder+1)*(Xorder+1))
 IF(Xorder .GE. 1)THEN
    IF(Xorder .GT. 1)THEN
       CALL STANDARDLOOPX(M000,MX,Xpq,Xorder,nPrim,P,JMAX,lupri)
    ELSEIF(Xorder .EQ. 1)THEN
       CALL ORDER_EQ_ONE_LOOP(M000,MX,Xpq,Xorder,nPrim,JMAX,lupri)
    ENDIF
 ELSE
    CALL DCOPY(nPrim,M000,1,MX(:,0,0),1)
 ENDIF

 CALL LS_DZERO(MY,nPrim*(Yorder+1)*(Yorder+1))
 IF(Yorder .GE. 1)THEN
    IF(Yorder .GT. 1)THEN
       CALL STANDARDLOOPX(M000,MY,Ypq,Yorder,nPrim,P,JMAX,lupri)
    ELSEIF(Yorder .EQ. 1)THEN
       CALL ORDER_EQ_ONE_LOOP(M000,MY,Ypq,Yorder,nPrim,JMAX,lupri)
    ENDIF
 ELSE
    CALL DCOPY(nPrim,M000,1,MY(:,0,0),1)
 ENDIF

 CALL LS_DZERO(MZ,nPrim*(Zorder+1)*(Zorder+1))
 IF(Zorder .GE. 1)THEN
    IF(Zorder .GT. 1)THEN
       CALL STANDARDLOOPX(M000,MZ,Zpq,Zorder,nPrim,P,JMAX,lupri)
    ELSEIF(Zorder .EQ. 1)THEN
       CALL ORDER_EQ_ONE_LOOP(M000,MZ,Zpq,Zorder,nPrim,JMAX,lupri)
    ENDIF
 ELSE
    CALL DCOPY(nPrim,M000,1,MZ(:,0,0),1)
 ENDIF
 efg=0
 !Nefg=(CARTORDER+1)*(CARTORDER+2)*(CARTORDER+3)/6
 DO C=0,CARTORDER
    DO e = C,0,-1
       DO f = C-e,0,-1
          g=C-e-f
          efg=efg+1
          TUV=0
          DO J = JSTART, JMAX
             DO T = J,0,-1
                DO U = J-T,0,-1
                   V=J-T-U
                   TUV=TUV+1
                   IF(T .LE. e)THEN
                      IF(U .LE. f)THEN
                         IF(V .LE. g)THEN
                            offset = (TUV-1)*nPrim+(efg-1)*nTUV*nPrim
                            DO I=1,nPrim
                               OUTPUT(I + offset)=MX(I,e,T)*MY(I,f,U)*MZ(I,g,V)
                            ENDDO
                         ELSE
                            offset = (TUV-1)*nPrim+(efg-1)*nTUV*nPrim
                            DO I=1,nPrim
                               OUTPUT(I + offset)=0E0_realk
                            ENDDO
                         ENDIF
                      ELSE
                         offset = (TUV-1)*nPrim+(efg-1)*nTUV*nPrim
                         DO I=1,nPrim
                            OUTPUT(I + offset)=0E0_realk
                         ENDDO
                      ENDIF
                   ELSE
                      offset = (TUV-1)*nPrim+(efg-1)*nTUV*nPrim
                      DO I=1,nPrim
                         OUTPUT(I + offset ) = 0.0E0_realk
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO

END SUBROUTINE multipoleRecurrence

!> \brief single direction multipole moment hermite integral recurrence relations \f$ M^{e+1}_{t} = t M^{e}_{t-1}+X_{PC}M^{e}_{t}+\frac{1}{2p}M^{e}_{t+1} \f$ 
!> \author \latexonly T. Kj{\ae}rgaard \endlatexonly
!> \date 2009-02-05
!> \param INPUT \f$ M^{0}_{0}\f$
!> \param M output \f$ M^{e}_{t}\f$ 
!> \param Xpq distance in given direction
!> \param Xorder expansion order in given direction
!> \param nPrim number of primitive functions
!> \param P exponents
!> \param LIMIT the maximum angular moment for the recurrence
!> \param LUPRI logical unit number for output printing
!> 
!> written for X loop (therefore Xpq and Xorder) but general for any direction
!>
SUBROUTINE StandardLoopX(INPUT,M,Xpq,Xorder,nPrim,P,LIMIT,lupri)
 IMPLICIT NONE
 INTEGER,intent(in)       :: nPrim,LIMIT,Xorder,lupri
 real(realk),intent(in)   :: Xpq(nPrim),P(nPrim),INPUT(nPrim)
 real(realk),intent(inout):: M(nPrim,0:Xorder,0:Xorder)
!
 INTEGER   :: C,T,I

 CALL DCOPY(nPrim,INPUT,1,M(1,0,0),1)

 !M(prim,e,t)
 DO C=1,Xorder-1
    IF(C .EQ. 1 )THEN
       DO I = 1, Nprim
          M(I,1,0)=Xpq(I)*INPUT(I) 
       ENDDO
       DO I = 1, Nprim
          M(I,1,1)=INPUT(I)              
       ENDDO
    ELSEIF(C .EQ. 2)THEN
       DO I = 1, Nprim
          M(I,2,0)=Xpq(I)*M(I,1,0)+P(I)*M(I,1,1)     
       ENDDO
       DO I = 1, Nprim
          M(I,2,1)=M(I,1,0)+Xpq(I)*M(I,1,1)          
       ENDDO
       DO I = 1, Nprim
          M(I,2,2)=2*M(I,1,1)                      
       ENDDO
    ELSE
       DO I = 1, Nprim
          M(I,C,0)= Xpq(I)*M(I,C-1,0)+P(I)*M(I,C-1,1)     
       ENDDO
       DO T = 1,C-2
          DO I = 1, Nprim
             M(I,C,T)=T*M(I,C-1,T-1)+Xpq(I)*M(I,C-1,T)+P(I)*M(I,C-1,T+1)
          ENDDO
       ENDDO
       !        T=C-1
       DO I = 1, Nprim
          M(I,C,C-1)=(C-1)*M(I,C-1,C-2)+Xpq(I)*M(I,C-1,C-1)
       ENDDO
       !        T=C
       DO I = 1, Nprim
          M(I,C,C)=C*M(I,C-1,C-1)
       ENDDO
    ENDIF
 ENDDO
 C=Xorder
 IF(C .EQ. 2)THEN
    DO I = 1, Nprim
       M(I,2,0)=Xpq(I)*M(I,1,0)+P(I)*M(I,1,1)     
    ENDDO
    IF(LIMIT .GT. 0)THEN
       DO I = 1, Nprim
          M(I,2,1)=M(I,1,0)+Xpq(I)*M(I,1,1)          
       ENDDO
       IF(LIMIT .GT. 1)THEN
          DO I = 1, Nprim
             M(I,2,2)=2*M(I,1,1)                      
          ENDDO
       ENDIF
    ENDIF
 ELSE
    DO I = 1, Nprim
       M(I,C,0)= Xpq(I)*M(I,C-1,0)+P(I)*M(I,C-1,1)     
    ENDDO
    DO T = 1,MIN(LIMIT,C-2)
       DO I = 1, Nprim
          M(I,C,T)=T*M(I,C-1,T-1)+Xpq(I)*M(I,C-1,T)+P(I)*M(I,C-1,T+1)
       ENDDO
    ENDDO
    IF(LIMIT .GT. C-2)THEN
       !        T=C-1
       DO I = 1, Nprim
          M(I,C,C-1)=(C-1)*M(I,C-1,C-2)+Xpq(I)*M(I,C-1,C-1)
       ENDDO
       IF(LIMIT .GT. C-1)THEN
          !        T=C
          DO I = 1, Nprim
             M(I,C,C)=C*M(I,C-1,C-1)
          ENDDO
       ENDIF
    ENDIF
 ENDIF

END SUBROUTINE StandardLoopX

!> \brief special multipole moment recurrence relations for expansion order equal 1
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param INPUT \f$ M^{0}_{0}\f$
!> \param M output \f$ M^{e}_{t}\f$ 
!> \param Xpq distance in given direction
!> \param Xorder expansion order in given direction
!> \param nPrim number of primitive functions
!> \param LUPRI logical unit number for output printing
!> 
!> written for X loop (therefore Xpq and Xorder) but general for any direction
!>
SUBROUTINE order_eq_one_loop(INPUT,M,Xpq,Xorder,nPrim,LIMIT,lupri)
 IMPLICIT NONE
 INTEGER,intent(in)     :: nPrim,Xorder,lupri,LIMIT
 real(realk),intent(inout) :: M(nPrim,0:Xorder,0:Xorder)
 real(realk),intent(in) :: Xpq(nPrim),INPUT(nPrim)
!
 INTEGER                :: I
 CALL DCOPY(nPrim,INPUT,1,M(1,0,0),1)
 DO I = 1, Nprim
    M(I,1,0)=Xpq(I)*INPUT(I)  
 ENDDO
 IF(LIMIT .GT. 0)THEN
    DO I = 1, Nprim
       M(I,1,1)=INPUT(I)                        
    ENDDO
 ENDIF
END SUBROUTINE order_eq_one_loop

!> \brief wrapper which copy integralprefactors and reducedexponents and then call secondary wrapper to 2 electron complementary error function modified Coulomb hermite integral
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ \tilde{R}^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param integral containing tabulated boys function \f$ F_{n} \f$ 
!> \param Omega argument in the complementary error function erfc\f$(\omega r_{12})\f$
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildErfRJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,OMEGA,HIGH)
 IMPLICIT NONE
 INTEGER,intent(in)           :: nPrim,LUPRI,IPRINT,Jmax
 REAL(REALK),intent(inout)    :: RJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)       :: OMEGA
 TYPE(Integrand),intent(in)   :: PQ
 TYPE(integralitem),intent(in):: integral
 LOGICAL,intent(in)           :: HIGH
!
 REAL(REALK) :: Prefactor(nPrim),alpha(nPrim)

 CALL DCOPY(nPrim,PQ%integralPrefactor,1,Prefactor,1)
 CALL DCOPY(nPrim,PQ%reducedExponents,1,alpha,1)
 !COPY BECAUSE BUILD_ERF2_RJ000 changes the values of prefactor and ALPHA
 CALL buildErf2RJ000(OMEGA,RJ000,nPrim,Jmax,LUPRI,IPRINT,&
      &alpha,PQ%squaredDistance,&
      &Prefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2,HIGH)

END SUBROUTINE buildErfRJ000

!> \brief secondary wrapper to 2 electron complementary error function modified Coulomb hermite integral 
!> \f$ \tilde{R}^{n}_{000}(\alpha,\textbf{R}_{PQ}) =  \left(-2 \alpha \frac{\omega^{2}}{\alpha + \omega} \right)^{n}  \sqrt{\frac{\omega^{2}}{\alpha + \omega}} F_{n}(\alpha \frac{\omega^{2}}{\alpha + \omega} R_{PQ}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param OMEGA argument in the complementary error function erfc\f$(\omega r_{12})\f$
!> \param RJ000 \f$ \tilde{R}^{0:jmax}_{000} \f$
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param alpha reduced exponent (to be modified to \f$ \frac{\omega^{2}}{\alpha + \omega} \f$)
!> \param R2 squared distance
!> \param Prefactor \f$ \frac{2 \pi^{\frac{5}{2}}}{pq\sqrt{p+q}} \f$ to be modified to \f$ \sqrt{\frac{\omega^{2}}{\alpha + \omega}} \frac{2 \pi^{\frac{5}{2}}}{pq\sqrt{p+q}} \f$
!> \param TABFJW tabulated boys function \f$ F_{n} \f$ 
!> \param nTABFJW1 dimension 1 of TABFJW
!> \param nTABFJW2 dimension 2 of TABFJW
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildErf2RJ000(OMEGA,RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2,HIGH)
 IMPLICIT NONE
 INTEGER,intent(in)         :: nPrim,Jmax,Lupri,Iprint,nTABFJW1,nTABFJW2
 REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrim),alpha(nPrim),Prefactor(nprim)
 REAL(REALK),intent(in)     :: R2(nprim),TABFJW(0:nTABFJW1,0:nTABFJW2),OMEGA
 LOGICAL,intent(in)         :: HIGH
!
 Real(realk) :: OMEGA2,BETAAT,SQBETA
 INTEGER     :: I

 OMEGA2=OMEGA*OMEGA
 DO I = 1, NPrim
    BETAAT = OMEGA2 / (ALPHA(I) + OMEGA2)  ! beta = omega^2/(alpha+omega^2) 
    SQBETA = SQRT(BETAAT)
    Prefactor(I) = Prefactor(I)*SQBETA
    ALPHA(I) = ALPHA(I)*BETAAT         
 ENDDO

 IF (HIGH) THEN
    !HIGH ACCURACY
    CALL buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
 ELSE
    !LOW ACCURACY
    CALL buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
 ENDIF

END SUBROUTINE buildErf2RJ000

!> \brief Distribute the primitive hermite integrals stored using tuvPQ into OUTPUT(tuvQ,tuvP,nPrimP,nPrimQ) using tuvP and tuvQ
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param tuvTUV output (ntuvQ,ntuvP,nPrimPQ) = (ntuvQ,ntuvP,nPrimP,nPrimQ)
!> \param WTUV input (ntuvPQ,nPrimPQ) 
!> \param TUVindex indexing for the tuv components
!> \param nPrim number of primitive orbitals = nPrimP*nPrimQ
!> \param startP start angular momentum for the P overlap distribution
!> \param endP end angular momentum for the P overlap distribution
!> \param startQ start angular momentum for the Q overlap distribution
!> \param endQ end angular momentum for the Q overlap distribution
!> \param ioffPQ ofset for the tuvPQ index
!> \param nTUVEFGPQ number of tuv components and derivative efg components
!> \param ntuvPQ  number of tuv components for PQ 
!> \param ntuvP number of tuv components for P overlap distribution
!> \param ntuvQ  number of tuv components for Q overlap distribution
!> \param ideriv derivative index
!> \param lupri logical unit number for output printing
SUBROUTINE distributeHermiteQ1(tuvTUV,WTUV,TUVindex,nPrim,startP,endP,startQ,endQ,&
    &                         ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,lupri)
 implicit none
 Integer,intent(in) :: nPrim,startP,endP,startQ,endQ,ioffPQ!,ioffP,ioffQ
 Integer,intent(in) :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv
 Integer,pointer    :: TUVindex(:,:,:)
 Real(realk),intent(inout) :: tuvTUV(ntuvQ,ntuvP,nPrim)
 Real(realk),intent(in)    :: WTUV(ntuvEFGPQ,nPrim)
 !
 Integer     :: TUVPQindex(nTUVP*nTUVQ)
 Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ
 Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ

 iPrimPQ=1
 iOFF=(ideriv-1)*ntuvPQ-ioffPQ
 ituvP = 0
 inTUVPQ = 0
 DO jP = startP,endP
    DO tP=jP,0,-1
       DO uP=jP-tP,0,-1
          vP=jP-tP-uP
          ituvP = ituvP+1
          ituvQ = 0
          DO jQ = startQ,endQ
             DO tQ=jQ,0,-1
                DO uQ=jQ-tQ,0,-1
                   vQ=jQ-tQ-uQ
                   ituvQ=ituvQ+1
                   ituvPQ=TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
                   inTUVPQ = inTUVPQ+1
                   TUVPQindex(inTUVPQ)=ituvPQ 
                   tuvTUV(ituvQ,ituvP,iPrimPQ) = WTUV(ituvPQ,iPrimPQ)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO
 DO iPrimPQ=2,nPrim
    inTUVPQ = 0
    DO iituvP = 1,ituvP
       DO iituvQ = 1,ituvQ
          inTUVPQ = inTUVPQ+1
          ituvPQ = TUVPQindex(inTUVPQ)
          tuvTUV(iituvQ,iituvP,iPrimPQ) = WTUV(ituvPQ,iPrimPQ)
       ENDDO
    ENDDO
 ENDDO

END SUBROUTINE distributeHermiteQ1

!> \brief Print the Distribute primitive hermite integrals stored using (tuvQ,tuvP,nPrimPQ)
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param tuvTUV to be printed
!> \param nPrim number of primitive orbitals = nPrimP*nPrimQ
!> \param startP start angular momentum for the P overlap distribution
!> \param endP end angular momentum for the P overlap distribution
!> \param startQ start angular momentum for the Q overlap distribution
!> \param endQ end angular momentum for the Q overlap distribution
!> \param ntuvP number of tuv components for P overlap distribution
!> \param ntuvQ  number of tuv components for Q overlap distribution
!> \param lupri logical unit number for output printing
SUBROUTINE printHermitePQ(tuvTUV,nPrim,startP,endP,startQ,endQ,&
    &                   ntuvP,ntuvQ,lupri)
 implicit none
 Integer,intent(in)     :: nPrim,startP,endP,startQ,endQ
 Integer,intent(in)     :: ntuvP,ntuvQ,lupri
 Real(realk),intent(in) :: tuvTUV(ntuvQ,ntuvP,nPrim)
 !
 Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ
 Integer     :: iPrimPQ

 ituvP = 0
 WRITE(LUPRI,'(3X,A)') '***************************************************************'
 WRITE(LUPRI,'(3X,A)') '***                         HermitePQ'
 WRITE(LUPRI,'(3X,A)') '***************************************************************'
 DO jP = startP,endP
    DO tP=jP,0,-1
       DO uP=jP-tP,0,-1
          vP=jP-tP-uP
          ituvP = ituvP+1
          ituvQ = 0
          DO jQ = startQ,endQ
             DO tQ=jQ,0,-1
                DO uQ=jQ-tQ,0,-1
                   vQ=jQ-tQ-uQ
                   ituvQ=ituvQ+1
                   WRITE(LUPRI,'(5X,A,I1,A,I1,A,I1,A,I1,A,I1,A,I1,A)') &
                        &            'W(',tP,',',uP,',',vP,'|',tQ,',',uQ,',',vQ,') ='

                   WRITE(LUPRI,'(5X,6ES10.4/,(5X,6ES10.4))') &
                        &           (tuvTUV(ituvQ,ituvP,iPrimPQ),iPrimPQ=1,nPrim)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO

END SUBROUTINE printHermitePQ

!> \brief wrapper to add calculated contracted integrals to temporary array order PA
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param TUVQ store the integrals until the Q overlap is fully contracted
!> \param OldTUVQ completed contracted integrals
!> \param totOrbQ total number of Orbitals on overlap distribution Q
!> \param startOrb start orbital to place in TUVQ
!> \param endOrb end orbital to place in TUVQ
!> \param nTUVP size of 2. dimension number of tuv components
!> \param nPrimP number of primitive orbitals on overlap distribution P
!> \param nCompQ number of orbital components on overlap distribution Q
!> \param nContQ number of contracted orbitals on overlap distribution Q
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE addToTUVQ_PA(TUVQ,OldTUVQ,totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
 implicit none
 Integer,intent(in)     :: totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
 Real(realk),intent(inout) :: TUVQ(nPrimP,nTUVP,totOrbQ)
 Real(realk),intent(in) :: OldTUVQ(nContQ,nPrimP,nTUVP,nCompQ)

 CALL addToTUVQ1_PA(TUVQ(:,:,startOrb:endOrb),OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)

END SUBROUTINE addToTUVQ_PA

!> \brief routine to add calculated contracted integrals to temporary array order PA
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param TUVQ store the integrals until the Q overlap is fully contracted
!> \param OldTUVQ completed contracted integrals
!> \param nTUVP size of 2. dimension number of tuv components
!> \param nPrimP number of primitive orbitals on overlap distribution P
!> \param nCompQ number of orbital components on overlap distribution Q
!> \param nContQ number of contracted orbitals on overlap distribution Q
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE addToTUVQ1_PA(TUVQ,OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
 implicit none
 Integer,intent(in)     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
 Real(realk),intent(inout) :: TUVQ(nPrimP,nTUVP,nCompQ,nContQ)
 Real(realk),intent(in) :: OldTUVQ(nContQ,nPrimP,nTUVP,nCompQ)
 !
 Integer :: iTUVP,iPrimP,iCompQ,iContQ
 !
 DO iContQ=1,nContQ
    DO iCompQ=1,nCompQ
       DO iTUVP=1,nTUVP
          DO iPrimP=1,nPrimP
             TUVQ(iPrimP,iTUVP,iCompQ,iContQ) = OldTUVQ(iContQ,iPrimP,iTUVP,iCompQ)
          ENDDO
       ENDDO
    ENDDO
 ENDDO

CALL PrintTuvQ_old(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)

END SUBROUTINE addToTUVQ1_PA

!> \brief print routine for TUVQ the temporary array to store integrals fully Q contracted
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param TUVQ the multi dimension array to be printet 
!> \param nTUVP size of 2. dimension
!> \param nPrimP size of 1. dimension
!> \param nCompQ size of 3. dimension
!> \param nContQ size of 4. dimension
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE PrintTuvQ_old(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
 implicit none
 Integer,intent(in)     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
 Real(realk),intent(in) :: TUVQ(nPrimP,nTUVP,nCompQ,nContQ)
 !
 Integer :: iTUVP,iPrimP,iCompQ,iContQ
 !
 IF (IPRINT.GT. 10) THEN
    CALL LSHEADER(LUPRI,'TUVQ')
    DO iContQ=1,nContQ
       DO iCompQ=1,nCompQ
          DO iTUVP=1,nTUVP
             WRITE(LUPRI,'(5X,A,I3,A,I3,A,I3)') 'iTUVP =',iTUVP,' iCompQ =',iCompQ,' iContQ =',iContQ
             WRITE(LUPRI,'(5X,5ES10.4)')  (TUVQ(iPrimP,iTUVP,iCompQ,iContQ), iPrimP=1,nPrimP)
          ENDDO
       ENDDO
    ENDDO
 ENDIF
END SUBROUTINE PrintTuvQ_old

END MODULE Thermite_integrals
