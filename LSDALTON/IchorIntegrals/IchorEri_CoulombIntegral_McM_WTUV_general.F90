!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriCoulombintegralCPUMcMGeneralWTUVMod
use IchorPrecisionModule
private 
public :: IchorwtuvRecurrenceJMIN0JMAX0,IchorwtuvRecurrenceJMIN0JMAX1,&
     & IchorwtuvRecurrenceJMIN0JMAX2,IchorwtuvRecurrenceJMIN0JMAX3,&
     & IchorwtuvRecurrenceCurrent, PrintWTUV, IchorTUVindexFuncFull,&
     & IchorTUVindexFunc,IchorTUVindexFuncJ,IchorwtuvRecurrenceJMIN0JMAX3J

CONTAINS
  subroutine IchorwtuvRecurrenceCurrent(OLD,CUR,J,JMAX,nPrimQ,nPrimP,nPasses,&
       & ntuv,Rpq,RJ000)
    implicit none 
    integer,intent(in)     :: J,nPrimQ,nPrimP,nPasses,ntuv,JMAX
    real(realk),intent(in) :: Rpq(nPrimQ*nPrimP*nPasses,3)
    real(realk),intent(in) :: RJ000(0:JMAX,nPrimQ*nPrimP*nPasses)
    Real(realk),intent(in) :: OLD(nPrimQ*nPrimP*nPasses,ntuv) 
    Real(realk),intent(inout) :: CUR(nPrimQ*nPrimP*nPasses,ntuv) 
    !
    !local variables
    integer :: n,iPrimQ,iPrimP,iPassP
    integer :: ituv,k,TUVauxJm1,TUVauxJm2,t,u,v,tm1,tm2,m1,ituvm1,ituvm2
    integer :: um1,vm1,um2,vm2,dir
    real(realk) :: Xtemp,Ytemp,Ztemp,WJ,WJ1,WJ2,WJ3
    real(realk) :: W100,W010,W001,W200,W110,W101,W020,W011,W002
    !$OMP DO PRIVATE(n,Xtemp,Ytemp,Ztemp,WJ,WJ1,WJ2,WJ3,W100,W010,W001,&
    !$OMP            W200,W110,W101,W020,W011,W002)
    DO n=1,nPrimQ*nPrimP*nPasses
       Xtemp = Rpq(n,1)
       Ytemp = Rpq(n,2)
       Ztemp = Rpq(n,3)
       WJ   = RJ000(j,n)
       WJ1  = RJ000(j+1,n)
       WJ2  = RJ000(j+2,n)
       WJ3  = RJ000(j+3,n)
       W100 = Xtemp*WJ2
       W010 = Ytemp*WJ2
       W001 = Ztemp*WJ2
       W200 = WJ2 + Xtemp*Xtemp*WJ3
       W110 =       Xtemp*Ytemp*WJ3
       W101 =       Xtemp*Ztemp*WJ3
       W020 = WJ2 + Ytemp*Ytemp*WJ3
       W011 =       Ytemp*Ztemp*WJ3
       W002 = WJ2 + Ztemp*Ztemp*WJ3
       CUR(n,1)  = WJ                            !000
       CUR(n,2)  = Xtemp*WJ1                     !100
       CUR(n,3)  = Ytemp*WJ1                     !010
       CUR(n,4)  = Ztemp*WJ1                     !001
       CUR(n,5)  = WJ1 + Xtemp*W100              !200
       CUR(n,6)  =       Ytemp*W100              !110
       CUR(n,7)  =       Ztemp*W100              !101
       CUR(n,8)  = WJ1 + Ytemp*W010              !020
       CUR(n,9)  =       Ztemp*W010              !011
       CUR(n,10) = WJ1 + Ztemp*W001              !002
       CUR(n,11) = 2.0E0_realk*W100 + Xtemp*W200 !300
       CUR(n,12) =       W010 + Xtemp*W110       !210
       CUR(n,13) =       W001 + Xtemp*W101       !201
       CUR(n,14) =       W100 + Ytemp*W110       !120
       CUR(n,15) =              Xtemp*W011       !111
       CUR(n,16) =       W100 + Ztemp*W101       !102
       CUR(n,17) = 2.0E0_realk*W010 + Ytemp*W020 !030
       CUR(n,18) =       W001 + Ytemp*W011       !021
       CUR(n,19) =       W010 + Ztemp*W011       !012
       CUR(n,20) = 2.0E0_realk*W001 + Ztemp*W002 !003
    ENDDO
    !$OMP END DO
    ituv = 20
    DO k=4,JMAX-j
       TUVauxJm1 = IchorTUVindexFuncJ(k-1)
       TUVauxJm2 = IchorTUVindexFuncJ(k-2)
       DO t=k,k-k/2,-1 !special case where we know (t.ge.u).AND.(t.ge.v)
          tm1 = t-1
          tm2 = t-2
          m1  = t-1
          DO u=k-t,0,-1
             v=k-t-u
             ituv   = ituv + 1
             !Thanks to Ulf Ekstroem for calculation of TUV index 
             !ituvm1=TUVauxJm1+ichorfactorial(1+U+V)/(ichorfactorial(U+V-1)*ichorfactorial(2))+V+1
             !ituvm2=TUVauxJm2+ichorfactorial(1+U+V)/(ichorfactorial(U+V-1)*ichorfactorial(2))+V+1
             ituvm1 = TUVauxJm1 + 1 + V + ((U+V+1)*(U+V))/2
             ituvm2 = TUVauxJm2 + 1 + V + ((U+V+1)*(U+V))/2
!$OMP DO PRIVATE(n)
             DO n=1,nPrimQ*nPrimP*nPasses
                CUR(n,ituv) = m1*OLD(n,ituvm2) + Rpq(n,1)*OLD(n,ituvm1)
             ENDDO
!$OMP END DO
          ENDDO
       ENDDO
       DO t=k-k/2-1,0,-1
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
             !Thanks to Ulf Ekstroem for IchorTUVindexFunc
             ituvm1 = IchorTUVindexFunc(um1,vm1,TUVauxJm1)
             ituvm2 = IchorTUVindexFunc(um2,vm2,TUVauxJm2)
!$OMP DO PRIVATE(n)
             DO n=1,nPrimQ*nPrimP*nPasses
                CUR(n,ituv) = m1*OLD(n,ituvm2) + Rpq(n,dir)*OLD(n,ituvm1)
             ENDDO
!$OMP END DO
          ENDDO
       ENDDO
    ENDDO
  end subroutine IchorwtuvRecurrenceCurrent

  !Calculates the TUV index using binomial coefficient function
  !The componentes are orderd using
  !graded reverse lexicographical ordering:
  !(000)
  !(100)
  !(010)
  !(001)
  !(200)
  !(110)
  !(101)
  !(020)
  !(011)
  !(002)
  !(300)
  !(210)
  !(201)
  !etc.
  ! where the first denote the T value, second the U values, etc. 
  ! Huge thanks to Ulf Ekstroem for this calculation!
  integer function IchorTUVindexFuncFull(T,U,V)
    implicit none
    integer,intent(in) :: T,U,V
!    Old Version Using Ichornchoosek
!    IchorTUVindexFuncFull = Ichornchoosek(2+T+U+V,T+U+V-1) + Ichornchoosek(1+U+V,U+V-1) + V + 1
!    Old Version that manual inline Ichornchoosek
!    IchorTUVindexFuncFull = ichorfactorial(2+T+U+V)/(ichorfactorial(T+U+V-1)*ichorfactorial(3)) + &
!         & ichorfactorial(1+U+V)/(ichorfactorial(U+V-1)*ichorfactorial(2)) + V + 1
!    Version that use the cheaper calculation
    IchorTUVindexFuncFull = 1 + V + ((U+V+1)*(U+V))/2+((T+U+V+2)*(T+U+V+1)*(T+U+V))/6
  end function IchorTUVindexFuncFull

  !This function assumes that the first term have been precalculated 
  integer function IchorTUVindexFunc(U,V,TUVauxJ)
    implicit none
    integer,intent(in) :: U,V,TUVauxJ
!    Old Version Using Ichornchoosek
!    IchorTUVindexFunc = TUVauxJ + Ichornchoosek(1+U+V,U+V-1) + V + 1
!    Old Version that manual inline Ichornchoosek
!    IchorTUVindexFunc = TUVauxJ + ichorfactorial(1+U+V)/(ichorfactorial(U+V-1)*ichorfactorial(2)) + V + 1    
!    Version that use the cheaper calculation
    IchorTUVindexFunc = TUVauxJ + 1 + V + ((U+V+1)*(U+V))/2
  end function IchorTUVindexFunc

  !Function to precalculated the first term
  integer function IchorTUVindexFuncJ(J)
    implicit none
    integer,intent(in) :: J
!    Old Version Using Ichornchoosek
!    IchorTUVindexFunc = Ichornchoosek(2+J,J-1)
!    Old Version that manual inline Ichornchoosek
!    IchorTUVindexFuncJ = ichorfactorial(2+J)/(ichorfactorial(J-1)*ichorfactorial(3))
!    Version that use the cheaper calculation
    IchorTUVindexFuncJ = ((J+2)*(J+1)*J)/6 ! = ((T+U+V+2)*(T+U+V+1)*(T+U+V))/6
  end function IchorTUVindexFuncJ

  integer function Ichornchoosek(n,k)
    implicit none
    integer,intent(in) :: n,k
    Ichornchoosek = ichorfactorial(n)/(ichorfactorial(k)*ichorfactorial(n-k))
  end function Ichornchoosek

  integer function Ichorfactorial(n)
    implicit none
    integer,intent(in) :: n
    integer :: n2,i
    n2 = 1
    do i=1,n
       n2 = n2*i
    enddo
    Ichorfactorial = n2
  end function Ichorfactorial

  SUBROUTINE IchorwtuvRecurrenceJMIN0JMAX0(WJ000,WTUV,nPrim)
    implicit none
    Integer,intent(in)           :: nPrim
    Real(realk),intent(in)       :: WJ000(0:0,nPrim)
    Real(realk),intent(inout)    :: WTUV(nPrim,1)
    !
    integer :: n
!$OMP DO PRIVATE(n)
    DO n=1,nPrim
       WTUV(n,1) = WJ000(0,n)
    ENDDO
!$OMP END DO
  END SUBROUTINE ICHORWTUVRECURRENCEJMIN0JMAX0

  SUBROUTINE ichorwtuvrecurrenceJMIN0JMAX1(WJ000,WTUV,Rpq,nPrim)
    implicit none
    Integer,intent(in)           :: nPrim
    Real(realk),intent(in)       :: WJ000(0:1,nPrim),Rpq(nPrim,3)
    Real(realk),intent(inout)    :: WTUV(nPrim,4)
    !
    integer :: n
!$OMP DO PRIVATE(n)
    DO n=1,nPrim
       WTUV(n,1) = WJ000(0,n)
       WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
       WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
       WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
    ENDDO
!$OMP END DO
  END SUBROUTINE ICHORWTUVRECURRENCEJMIN0JMAX1
  
SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX2(WJ000,WTUV,Rpq,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:2,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,10)
  !
  integer :: n
  !000
!$OMP DO PRIVATE(n)
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
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
!$OMP END DO
end SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX2

SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX3(WJ000,WTUV,Rpq,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:3,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,20)
  !
  Real(realk) :: WTUVTEMP_21,WTUVTEMP_22
  Real(realk) :: WTUVTEMP_31,WTUVTEMP_32
  Real(realk) :: WTUVTEMP_41,WTUVTEMP_42
  Real(realk) :: WTUVTEMP25,WTUVTEMP28
  Real(realk) :: WTUVTEMP210,WTUVTEMP29
  integer :: n
!$OMP DO PRIVATE(n)
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0,n)
     WTUV(n,2)   = Rpq(n,1)*WJ000(1,n)
     WTUVTEMP_21 = Rpq(n,1)*WJ000(2,n)
     WTUVTEMP_22 = Rpq(n,1)*WJ000(3,n)
     WTUV(n,3)   = Rpq(n,2)*WJ000(1,n)
     WTUVTEMP_31 = Rpq(n,2)*WJ000(2,n)
     WTUVTEMP_32 = Rpq(n,2)*WJ000(3,n)
     WTUV(n,4)   = Rpq(n,3)*WJ000(1,n)
     WTUVTEMP_41 = Rpq(n,3)*WJ000(2,n)
     WTUVTEMP_42 = Rpq(n,3)*WJ000(3,n)
     WTUV(n,6)   = Rpq(n,1)*WTUVTEMP_31
     WTUV(n,7)   = Rpq(n,1)*WTUVTEMP_41
     WTUV(n,9)   = Rpq(n,2)*WTUVTEMP_41
     WTUVTEMP29  = Rpq(n,2)*WTUVTEMP_42
     WTUV(n,5)   = WJ000(1,n) + Rpq(n,1)*WTUVTEMP_21
     WTUV(n,8)   = WJ000(1,n) + Rpq(n,2)*WTUVTEMP_31
     WTUV(n,10)  = WJ000(1,n) + Rpq(n,3)*WTUVTEMP_41
     WTUVTEMP25  = WJ000(2,n) + Rpq(n,1)*WTUVTEMP_22
     WTUVTEMP28  = WJ000(2,n) + Rpq(n,2)*WTUVTEMP_32
     WTUVTEMP210 = WJ000(2,n) + Rpq(n,3)*WTUVTEMP_42
     WTUV(n,14) = Rpq(n,1)*WTUVTEMP28
     WTUV(n,15) = Rpq(n,1)*WTUVTEMP29
     WTUV(n,16) = Rpq(n,1)*WTUVTEMP210
     WTUV(n,12) = Rpq(n,2)*WTUVTEMP25
     WTUV(n,19) = Rpq(n,2)*WTUVTEMP210
     WTUV(n,13) = Rpq(n,3)*WTUVTEMP25
     WTUV(n,18) = Rpq(n,3)*WTUVTEMP28
     WTUV(n,11) = 2*WTUVTEMP_21 + Rpq(n,1)*WTUVTEMP25
     WTUV(n,17) = 2*WTUVTEMP_31 + Rpq(n,2)*WTUVTEMP28
     WTUV(n,20) = 2*WTUVTEMP_41 + Rpq(n,3)*WTUVTEMP210
  ENDDO
!$OMP END DO
END SUBROUTINE ichorWTUVRECURRENCEJMIN0JMAX3

SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX3J(WJ000,WTUV,Rpq,nPrim,maxJ,J)
  implicit none
  Integer,intent(in)           :: nPrim,maxJ,J
  Real(realk),intent(in)       :: WJ000(0:maxJ,nPrim),Rpq(nPrim,3)
  Real(realk),intent(inout)    :: WTUV(nPrim,20)
  !
  Real(realk) :: WTUVTEMP_21,WTUVTEMP_22
  Real(realk) :: WTUVTEMP_31,WTUVTEMP_32
  Real(realk) :: WTUVTEMP_41,WTUVTEMP_42
  Real(realk) :: WTUVTEMP25,WTUVTEMP28
  Real(realk) :: WTUVTEMP210,WTUVTEMP29
  integer :: n
!$OMP DO PRIVATE(n)
  DO n=1,nPrim
     WTUV(n,1) = WJ000(0+J,n)
     WTUV(n,2)   = Rpq(n,1)*WJ000(1+J,n)
     WTUVTEMP_21 = Rpq(n,1)*WJ000(2+J,n)
     WTUVTEMP_22 = Rpq(n,1)*WJ000(3+J,n)
     WTUV(n,3)   = Rpq(n,2)*WJ000(1+J,n)
     WTUVTEMP_31 = Rpq(n,2)*WJ000(2+J,n)
     WTUVTEMP_32 = Rpq(n,2)*WJ000(3+J,n)
     WTUV(n,4)   = Rpq(n,3)*WJ000(1+J,n)
     WTUVTEMP_41 = Rpq(n,3)*WJ000(2+J,n)
     WTUVTEMP_42 = Rpq(n,3)*WJ000(3+J,n)
     WTUV(n,6)   = Rpq(n,1)*WTUVTEMP_31
     WTUV(n,7)   = Rpq(n,1)*WTUVTEMP_41
     WTUV(n,9)   = Rpq(n,2)*WTUVTEMP_41
     WTUVTEMP29  = Rpq(n,2)*WTUVTEMP_42
     WTUV(n,5)   = WJ000(1+J,n) + Rpq(n,1)*WTUVTEMP_21
     WTUV(n,8)   = WJ000(1+J,n) + Rpq(n,2)*WTUVTEMP_31
     WTUV(n,10)  = WJ000(1+J,n) + Rpq(n,3)*WTUVTEMP_41
     WTUVTEMP25  = WJ000(2+J,n) + Rpq(n,1)*WTUVTEMP_22
     WTUVTEMP28  = WJ000(2+J,n) + Rpq(n,2)*WTUVTEMP_32
     WTUVTEMP210 = WJ000(2+J,n) + Rpq(n,3)*WTUVTEMP_42
     WTUV(n,14) = Rpq(n,1)*WTUVTEMP28
     WTUV(n,15) = Rpq(n,1)*WTUVTEMP29
     WTUV(n,16) = Rpq(n,1)*WTUVTEMP210
     WTUV(n,12) = Rpq(n,2)*WTUVTEMP25
     WTUV(n,19) = Rpq(n,2)*WTUVTEMP210
     WTUV(n,13) = Rpq(n,3)*WTUVTEMP25
     WTUV(n,18) = Rpq(n,3)*WTUVTEMP28
     WTUV(n,11) = 2*WTUVTEMP_21 + Rpq(n,1)*WTUVTEMP25
     WTUV(n,17) = 2*WTUVTEMP_31 + Rpq(n,2)*WTUVTEMP28
     WTUV(n,20) = 2*WTUVTEMP_41 + Rpq(n,3)*WTUVTEMP210
  ENDDO
!$OMP END DO
END SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX3J

! WARNING CHANGE ORDER OF WTUV 
!!$
!!$SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX4(WJ000,WTUV,Rpq,nPrim)
!!$  implicit none
!!$  Integer,intent(in)           :: nPrim
!!$  Real(realk),intent(in)       :: WJ000(0:4,nPrim),Rpq(nPrim,3)
!!$  Real(realk),intent(inout)    :: WTUV(35,nPrim)
!!$  !
!!$  Integer     :: n,M
!!$  Real(realk) :: WTUVTEMP(2:4,3)
!!$  Real(realk) :: WTUVTEMP2(5:10,2)
!!$  Real(realk) :: WTUVTEMP3(11:20)
!!$  !000
!!$  DO n=1,nPrim
!!$     WTUV(1,n) = WJ000(0,n)
!!$  ENDDO
!!$  !it can reuse registers with this collection of this in 1 loop
!!$  DO n=1,nPrim
!!$     WTUV(2,n) = Rpq(n,1)*WJ000(1,n)
!!$     WTUV(3,n) = Rpq(n,2)*WJ000(1,n)
!!$     WTUV(4,n) = Rpq(n,3)*WJ000(1,n)
!!$     WTUVTEMP(2,1) = Rpq(n,1)*WJ000(2,n)
!!$     WTUVTEMP(3,1) = Rpq(n,2)*WJ000(2,n)
!!$     WTUVTEMP(4,1) = Rpq(n,3)*WJ000(2,n)     
!!$     WTUVTEMP(2,2) = Rpq(n,1)*WJ000(3,n)
!!$     WTUVTEMP(3,2) = Rpq(n,2)*WJ000(3,n)
!!$     WTUVTEMP(4,2) = Rpq(n,3)*WJ000(3,n)     
!!$     WTUVTEMP(2,3) = Rpq(n,1)*WJ000(4,n)
!!$     WTUVTEMP(3,3) = Rpq(n,2)*WJ000(4,n)
!!$     WTUVTEMP(4,3) = Rpq(n,3)*WJ000(4,n)
!!$     WTUV(5,n) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
!!$     WTUV(6,n) =              Rpq(n,1)*WTUVTEMP(3,1)
!!$     WTUV(7,n) =              Rpq(n,1)*WTUVTEMP(4,1)
!!$     WTUV(8,n) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
!!$     WTUV(9,n) =              Rpq(n,2)*WTUVTEMP(4,1)
!!$     WTUV(10,n)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
!!$     WTUVTEMP2(5,1) = WJ000(2,n) + Rpq(n,1)*WTUVTEMP(2,2)
!!$     WTUVTEMP2(6,1) =              Rpq(n,1)*WTUVTEMP(3,2)
!!$     WTUVTEMP2(7,1) =              Rpq(n,1)*WTUVTEMP(4,2)
!!$     WTUVTEMP2(8,1) = WJ000(2,n) + Rpq(n,2)*WTUVTEMP(3,2)
!!$     WTUVTEMP2(9,1) =              Rpq(n,2)*WTUVTEMP(4,2)
!!$     WTUVTEMP2(10,1)= WJ000(2,n) + Rpq(n,3)*WTUVTEMP(4,2)     
!!$     WTUVTEMP2(5,2) = WJ000(3,n) + Rpq(n,1)*WTUVTEMP(2,3)
!!$!     WTUVTEMP2(6,2) =              Rpq(n,1)*WTUVTEMP(3,3)
!!$!     WTUVTEMP2(7,2) =              Rpq(n,1)*WTUVTEMP(4,3)
!!$     WTUVTEMP2(8,2) = WJ000(3,n) + Rpq(n,2)*WTUVTEMP(3,3)
!!$     WTUVTEMP2(9,2) =              Rpq(n,2)*WTUVTEMP(4,3)
!!$     WTUVTEMP2(10,2)= WJ000(3,n) + Rpq(n,3)*WTUVTEMP(4,3)
!!$     WTUV(11,n)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
!!$     WTUV(12,n)=  Rpq(n,2)*WTUVTEMP2(5,1)
!!$     WTUV(13,n)= Rpq(n,3)*WTUVTEMP2(5,1)
!!$     WTUV(14,n)= Rpq(n,1)*WTUVTEMP2(8,1)
!!$     WTUV(15,n)= Rpq(n,1)*WTUVTEMP2(9,1)
!!$     WTUV(16,n)= Rpq(n,1)*WTUVTEMP2(10,1)
!!$     WTUV(17,n)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
!!$     WTUV(18,n)= Rpq(n,3)*WTUVTEMP2(8,1)
!!$     WTUV(19,n)= Rpq(n,2)*WTUVTEMP2(10,1)
!!$     WTUV(20,n)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
!!$     WTUVTEMP3(11)= 2*WTUVTEMP(2,2) + Rpq(n,1)*WTUVTEMP2(5,2)
!!$    ! WTUVTEMP3(12)=  Rpq(n,2)*WTUVTEMP2(5,2)
!!$     WTUVTEMP3(13)= Rpq(n,3)*WTUVTEMP2(5,2)
!!$     WTUVTEMP3(14)= Rpq(n,1)*WTUVTEMP2(8,2)
!!$    ! WTUVTEMP3(15)= Rpq(n,1)*WTUVTEMP2(9,2)
!!$     WTUVTEMP3(16)= Rpq(n,1)*WTUVTEMP2(10,2)
!!$     WTUVTEMP3(17)= 2*WTUVTEMP(3,2) + Rpq(n,2)*WTUVTEMP2(8,2)
!!$     WTUVTEMP3(18)= Rpq(n,3)*WTUVTEMP2(8,2)
!!$     WTUVTEMP3(19)= Rpq(n,2)*WTUVTEMP2(10,2)
!!$     WTUVTEMP3(20)= 2*WTUVTEMP(4,2) + Rpq(n,3)*WTUVTEMP2(10,2)
!!$     WTUV(21,n)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11)
!!$     WTUV(22,n)= Rpq(n,2)*WTUVTEMP3(11)
!!$     WTUV(23,n)= Rpq(n,3)*WTUVTEMP3(11)
!!$     WTUV(24,n)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14)
!!$     WTUV(25,n)= Rpq(n,2)*WTUVTEMP3(13)
!!$     WTUV(26,n)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16)
!!$     WTUV(27,n)= Rpq(n,1)*WTUVTEMP3(17)
!!$     WTUV(28,n)= Rpq(n,1)*WTUVTEMP3(18)
!!$     WTUV(29,n)= Rpq(n,1)*WTUVTEMP3(19)
!!$     WTUV(30,n)= Rpq(n,1)*WTUVTEMP3(20)
!!$     WTUV(31,n)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17) 
!!$     WTUV(32,n)= Rpq(n,3)*WTUVTEMP3(17)
!!$     WTUV(33,n)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19)
!!$     WTUV(34,n)= Rpq(n,2)*WTUVTEMP3(20)
!!$     WTUV(35,n)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20)
!!$  ENDDO
!!$END SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX4
!!$
!!$SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX5(WJ000,WTUV,Rpq,nPrim)
!!$  implicit none
!!$  Integer,intent(in)           :: nPrim
!!$  Real(realk),intent(in)       :: WJ000(0:5,nPrim),Rpq(nPrim,3)
!!$  Real(realk),intent(inout)    :: WTUV(nPrim,56)
!!$  !
!!$  Integer     :: n,M
!!$  Real(realk) :: WTUVTEMP(2:4,4)
!!$  Real(realk) :: WTUVTEMP2(5:10,3)
!!$  Real(realk) :: WTUVTEMP3(11:20,2)
!!$  Real(realk) :: WTUVTEMP4(21:35)
!!$  !000
!!$  DO n=1,nPrim
!!$     WTUV(n,1) = WJ000(0,n)
!!$  ENDDO
!!$  !it can reuse registers with this collection of this in 1 loop
!!$  DO n=1,nPrim
!!$     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
!!$     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
!!$     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
!!$     DO M=1,4
!!$        WTUVTEMP(2,M) = Rpq(n,1)*WJ000(M+1,n)
!!$        WTUVTEMP(3,M) = Rpq(n,2)*WJ000(M+1,n)
!!$        WTUVTEMP(4,M) = Rpq(n,3)*WJ000(M+1,n)     
!!$     ENDDO
!!$     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
!!$     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
!!$     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
!!$     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
!!$     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
!!$     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
!!$     DO M=1,3
!!$        WTUVTEMP2(5,M) = WJ000(M+1,n) + Rpq(n,1)*WTUVTEMP(2,M+1)
!!$!        WTUVTEMP2(6,M) =                Rpq(n,1)*WTUVTEMP(3,M+1)
!!$!        WTUVTEMP2(7,M) =                Rpq(n,1)*WTUVTEMP(4,M+1)
!!$        WTUVTEMP2(8,M) = WJ000(M+1,n) + Rpq(n,2)*WTUVTEMP(3,M+1)
!!$        WTUVTEMP2(9,M) =                Rpq(n,2)*WTUVTEMP(4,M+1)
!!$        WTUVTEMP2(10,M)= WJ000(M+1,n) + Rpq(n,3)*WTUVTEMP(4,M+1)     
!!$     ENDDO
!!$     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
!!$     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
!!$     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
!!$     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
!!$     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
!!$     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
!!$     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
!!$     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
!!$     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
!!$     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
!!$     DO M=1,2
!!$        WTUVTEMP3(11,M)= 2*WTUVTEMP(2,M+1) + Rpq(n,1)*WTUVTEMP2(5,M+1)
!!$!        WTUVTEMP3(12,M)=  Rpq(n,2)*WTUVTEMP2(5,M+1)
!!$        WTUVTEMP3(13,M)= Rpq(n,3)*WTUVTEMP2(5,M+1)
!!$        WTUVTEMP3(14,M)= Rpq(n,1)*WTUVTEMP2(8,M+1)
!!$!        WTUVTEMP3(15,M)= Rpq(n,1)*WTUVTEMP2(9,M+1)
!!$        WTUVTEMP3(16,M)= Rpq(n,1)*WTUVTEMP2(10,M+1)
!!$        WTUVTEMP3(17,M)= 2*WTUVTEMP(3,M+1) + Rpq(n,2)*WTUVTEMP2(8,M+1)
!!$        WTUVTEMP3(18,M)= Rpq(n,3)*WTUVTEMP2(8,M+1)
!!$        WTUVTEMP3(19,M)= Rpq(n,2)*WTUVTEMP2(10,M+1)
!!$        WTUVTEMP3(20,M)= 2*WTUVTEMP(4,M+1) + Rpq(n,3)*WTUVTEMP2(10,M+1)
!!$     ENDDO
!!$
!!$     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11,1)
!!$     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11,1)
!!$     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11,1)
!!$     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14,1)
!!$     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13,1)
!!$     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16,1)
!!$     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17,1)
!!$     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18,1)
!!$     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19,1)
!!$     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20,1)
!!$     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17,1) 
!!$     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17,1)
!!$     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19,1)
!!$     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20,1)
!!$     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20,1)
!!$
!!$     WTUVTEMP4(21)= 3*WTUVTEMP2(5,2) + Rpq(n,1)*WTUVTEMP3(11,2)
!!$!     WTUVTEMP4(22)= Rpq(n,2)*WTUVTEMP3(11,2)
!!$     WTUVTEMP4(23)= Rpq(n,3)*WTUVTEMP3(11,2)
!!$     WTUVTEMP4(24)= WTUVTEMP2(8,2) + Rpq(n,1)*WTUVTEMP3(14,2)
!!$!     WTUVTEMP4(25)= Rpq(n,2)*WTUVTEMP3(13,2)
!!$     WTUVTEMP4(26)= WTUVTEMP2(10,2) + Rpq(n,1)*WTUVTEMP3(16,2)
!!$     WTUVTEMP4(27)= Rpq(n,1)*WTUVTEMP3(17,2)
!!$!     WTUVTEMP4(28)= Rpq(n,1)*WTUVTEMP3(18,2)
!!$!     WTUVTEMP4(29)= Rpq(n,1)*WTUVTEMP3(19,2)
!!$     WTUVTEMP4(30)= Rpq(n,1)*WTUVTEMP3(20,2)
!!$     WTUVTEMP4(31)= 3*WTUVTEMP2(8,2) + Rpq(n,2)*WTUVTEMP3(17,2) 
!!$     WTUVTEMP4(32)= Rpq(n,3)*WTUVTEMP3(17,2)
!!$     WTUVTEMP4(33)= WTUVTEMP2(10,2) + Rpq(n,2)*WTUVTEMP3(19,2)
!!$     WTUVTEMP4(34)= Rpq(n,2)*WTUVTEMP3(20,2)
!!$     WTUVTEMP4(35)= 3*WTUVTEMP2(10,2) + Rpq(n,3)*WTUVTEMP3(20,2)
!!$
!!$     WTUV(n,  36) = 4*WTUVTEMP3(11,1) + Rpq(n,1)*WTUVTEMP4(21)
!!$     WTUV(n,  37) = Rpq(n,2)*WTUVTEMP4(21)
!!$     WTUV(n,  38) = Rpq(n,3)*WTUVTEMP4(21)
!!$     WTUV(n,  39) = 2*WTUVTEMP3(14,1) + Rpq(n,1)*WTUVTEMP4(24)
!!$     WTUV(n,  40) = Rpq(n,2)*WTUVTEMP4(23)
!!$     WTUV(n,  41) = 2*WTUVTEMP3(16,1) + Rpq(n,1)*WTUVTEMP4(26)
!!$     WTUV(n,  42) = WTUVTEMP3(17,1) + Rpq(n,1)*WTUVTEMP4(27)
!!$     WTUV(n,  43) = Rpq(n,3)*WTUVTEMP4(24)
!!$     WTUV(n,  44) = Rpq(n,2)*WTUVTEMP4(26)
!!$     WTUV(n,  45) = WTUVTEMP3(20,1) + Rpq(n,1)*WTUVTEMP4(30)
!!$     WTUV(n,  46) = Rpq(n,1)*WTUVTEMP4(31)
!!$     WTUV(n,  47) = Rpq(n,1)*WTUVTEMP4(32)
!!$     WTUV(n,  48) = Rpq(n,1)*WTUVTEMP4(33)
!!$     WTUV(n,  49) = Rpq(n,1)*WTUVTEMP4(34)
!!$     WTUV(n,  50) = Rpq(n,1)*WTUVTEMP4(35)
!!$     WTUV(n,  51) = 4*WTUVTEMP3(17,1) + Rpq(n,2)*WTUVTEMP4(31)
!!$     WTUV(n,  52) = Rpq(n,3)*WTUVTEMP4(31)
!!$     WTUV(n,  53) = 2*WTUVTEMP3(19,1) + Rpq(n,2)*WTUVTEMP4(33)
!!$     WTUV(n,  54) = WTUVTEMP3(20,1) + Rpq(n,2)*WTUVTEMP4(34)
!!$     WTUV(n,  55) = Rpq(n,2)*WTUVTEMP4(35)
!!$     WTUV(n,  56) = 4*WTUVTEMP3(20,1) + Rpq(n,3)*WTUVTEMP4(35)
!!$  ENDDO
!!$END SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX5
!!$
!!$SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX6(WJ000,WTUV,Rpq,nPrim)
!!$  implicit none
!!$  Integer,intent(in)           :: nPrim
!!$  Real(realk),intent(in)       :: WJ000(0:6,nPrim),Rpq(nPrim,3)
!!$  Real(realk),intent(inout)    :: WTUV(nPrim,84)
!!$  !
!!$  Integer     :: n,M
!!$  Real(realk) :: WTUVTEMP(2:4,5)
!!$  Real(realk) :: WTUVTEMP2(5:10,4)
!!$  Real(realk) :: WTUVTEMP3(11:20,3)
!!$  Real(realk) :: WTUVTEMP4(21:35,2)
!!$  Real(realk) :: WTUVTEMP5(36:56,1)
!!$  !000
!!$  !ERROR IN THIS SUBROUTINE SOMEWHERE
!!$  DO n=1,nPrim
!!$     WTUV(n,1) = WJ000(0,n)
!!$  ENDDO
!!$  !it can reuse registers with this collection of this in 1 loop
!!$  DO n=1,nPrim
!!$     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
!!$     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
!!$     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
!!$     DO M=1,5
!!$        WTUVTEMP(2,M) = Rpq(n,1)*WJ000(M+1,n)
!!$        WTUVTEMP(3,M) = Rpq(n,2)*WJ000(M+1,n)
!!$        WTUVTEMP(4,M) = Rpq(n,3)*WJ000(M+1,n)     
!!$     ENDDO
!!$     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
!!$     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
!!$     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
!!$     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
!!$     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
!!$     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
!!$     DO M=1,4
!!$        WTUVTEMP2(5,M) = WJ000(M+1,n) + Rpq(n,1)*WTUVTEMP(2,M+1)
!!$!        WTUVTEMP2(6,M) =                Rpq(n,1)*WTUVTEMP(3,M+1)
!!$!        WTUVTEMP2(7,M) =                Rpq(n,1)*WTUVTEMP(4,M+1)
!!$        WTUVTEMP2(8,M) = WJ000(M+1,n) + Rpq(n,2)*WTUVTEMP(3,M+1)
!!$        WTUVTEMP2(9,M) =                Rpq(n,2)*WTUVTEMP(4,M+1)
!!$        WTUVTEMP2(10,M)= WJ000(M+1,n) + Rpq(n,3)*WTUVTEMP(4,M+1)     
!!$     ENDDO
!!$     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
!!$     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
!!$     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
!!$     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
!!$     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
!!$     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
!!$     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
!!$     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
!!$     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
!!$     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
!!$     DO M=1,3
!!$        WTUVTEMP3(11,M)= 2*WTUVTEMP(2,M+1) + Rpq(n,1)*WTUVTEMP2(5,M+1)
!!$!        WTUVTEMP3(12,M)=  Rpq(n,2)*WTUVTEMP2(5,M+1)
!!$        WTUVTEMP3(13,M)= Rpq(n,3)*WTUVTEMP2(5,M+1)
!!$        WTUVTEMP3(14,M)= Rpq(n,1)*WTUVTEMP2(8,M+1)
!!$!        WTUVTEMP3(15,M)= Rpq(n,1)*WTUVTEMP2(9,M+1)
!!$        WTUVTEMP3(16,M)= Rpq(n,1)*WTUVTEMP2(10,M+1)
!!$        WTUVTEMP3(17,M)= 2*WTUVTEMP(3,M+1) + Rpq(n,2)*WTUVTEMP2(8,M+1)
!!$        WTUVTEMP3(18,M)= Rpq(n,3)*WTUVTEMP2(8,M+1)
!!$        WTUVTEMP3(19,M)= Rpq(n,2)*WTUVTEMP2(10,M+1)
!!$        WTUVTEMP3(20,M)= 2*WTUVTEMP(4,M+1) + Rpq(n,3)*WTUVTEMP2(10,M+1)
!!$     ENDDO
!!$
!!$     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11,1)
!!$     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11,1)
!!$     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11,1)
!!$     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14,1)
!!$     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13,1)
!!$     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16,1)
!!$     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17,1)
!!$     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18,1)
!!$     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19,1)
!!$     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20,1)
!!$     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17,1) 
!!$     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17,1)
!!$     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19,1)
!!$     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20,1)
!!$     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20,1)
!!$
!!$     DO M=1,2
!!$        WTUVTEMP4(21,M)= 3*WTUVTEMP2(5,M+1) + Rpq(n,1)*WTUVTEMP3(11,M+1)
!!$!        WTUVTEMP4(22,M)= Rpq(n,2)*WTUVTEMP3(11,M+1) !not used
!!$        WTUVTEMP4(23,M)= Rpq(n,3)*WTUVTEMP3(11,M+1)
!!$        WTUVTEMP4(24,M)= WTUVTEMP2(8,M+1) + Rpq(n,1)*WTUVTEMP3(14,M+1)
!!$!        WTUVTEMP4(25,M)= Rpq(n,2)*WTUVTEMP3(13,M+1) !not used
!!$        WTUVTEMP4(26,M)= WTUVTEMP2(10,M+1) + Rpq(n,1)*WTUVTEMP3(16,M+1)
!!$        WTUVTEMP4(27,M)= Rpq(n,1)*WTUVTEMP3(17,M+1)
!!$!        WTUVTEMP4(28,M)= Rpq(n,1)*WTUVTEMP3(18,M+1) !not used
!!$!        WTUVTEMP4(29,M)= Rpq(n,1)*WTUVTEMP3(19,M+1) !not used
!!$        WTUVTEMP4(30,M)= Rpq(n,1)*WTUVTEMP3(20,M+1) 
!!$        WTUVTEMP4(31,M)= 3*WTUVTEMP2(8,M+1) + Rpq(n,2)*WTUVTEMP3(17,M+1) 
!!$        WTUVTEMP4(32,M)= Rpq(n,3)*WTUVTEMP3(17,M+1)
!!$        WTUVTEMP4(33,M)= WTUVTEMP2(10,M+1) + Rpq(n,2)*WTUVTEMP3(19,M+1)
!!$        WTUVTEMP4(34,M)= Rpq(n,2)*WTUVTEMP3(20,M+1)
!!$        WTUVTEMP4(35,M)= 3*WTUVTEMP2(10,M+1) + Rpq(n,3)*WTUVTEMP3(20,M+1)
!!$     ENDDO
!!$
!!$     WTUV(n,  36) = 4*WTUVTEMP3(11,1) + Rpq(n,1)*WTUVTEMP4(21,1)
!!$     WTUV(n,  37) = Rpq(n,2)*WTUVTEMP4(21,1)
!!$     WTUV(n,  38) = Rpq(n,3)*WTUVTEMP4(21,1)
!!$     WTUV(n,  39) = 2*WTUVTEMP3(14,1) + Rpq(n,1)*WTUVTEMP4(24,1)
!!$     WTUV(n,  40) = Rpq(n,2)*WTUVTEMP4(23,1)
!!$     WTUV(n,  41) = 2*WTUVTEMP3(16,1) + Rpq(n,1)*WTUVTEMP4(26,1)
!!$     WTUV(n,  42) = WTUVTEMP3(17,1) + Rpq(n,1)*WTUVTEMP4(27,1)
!!$     WTUV(n,  43) = Rpq(n,3)*WTUVTEMP4(24,1)
!!$     WTUV(n,  44) = Rpq(n,2)*WTUVTEMP4(26,1)
!!$     WTUV(n,  45) = WTUVTEMP3(20,1) + Rpq(n,1)*WTUVTEMP4(30,1)
!!$     WTUV(n,  46) = Rpq(n,1)*WTUVTEMP4(31,1)
!!$     WTUV(n,  47) = Rpq(n,1)*WTUVTEMP4(32,1)
!!$     WTUV(n,  48) = Rpq(n,1)*WTUVTEMP4(33,1)
!!$     WTUV(n,  49) = Rpq(n,1)*WTUVTEMP4(34,1)
!!$     WTUV(n,  50) = Rpq(n,1)*WTUVTEMP4(35,1)
!!$     WTUV(n,  51) = 4*WTUVTEMP3(17,1) + Rpq(n,2)*WTUVTEMP4(31,1)
!!$     WTUV(n,  52) = Rpq(n,3)*WTUVTEMP4(31,1)
!!$     WTUV(n,  53) = 2*WTUVTEMP3(19,1) + Rpq(n,2)*WTUVTEMP4(33,1)
!!$     WTUV(n,  54) = WTUVTEMP3(20,1) + Rpq(n,2)*WTUVTEMP4(34,1)
!!$     WTUV(n,  55) = Rpq(n,2)*WTUVTEMP4(35,1)
!!$     WTUV(n,  56) = 4*WTUVTEMP3(20,1) + Rpq(n,3)*WTUVTEMP4(35,1)
!!$
!!$     WTUVTEMP5(36,1) = 4*WTUVTEMP3(11,2) + Rpq(n,1)*WTUVTEMP4(21,2)
!!$!     WTUVTEMP5(37,1) = Rpq(n,2)*WTUVTEMP4(21,2) !not used
!!$     WTUVTEMP5(38,1) = Rpq(n,3)*WTUVTEMP4(21,2)
!!$     WTUVTEMP5(39,1) = 2*WTUVTEMP3(14,2) + Rpq(n,1)*WTUVTEMP4(24,2)
!!$!     WTUVTEMP5(40,1) = Rpq(n,2)*WTUVTEMP4(23,2)!not used
!!$     WTUVTEMP5(41,1) = 2*WTUVTEMP3(16,2) + Rpq(n,1)*WTUVTEMP4(26,2)
!!$     WTUVTEMP5(42,1) = WTUVTEMP3(17,2) + Rpq(n,1)*WTUVTEMP4(27,2)
!!$!     WTUVTEMP5(43,1) = Rpq(n,3)*WTUVTEMP4(24,2)!not used
!!$!     WTUVTEMP5(44,1) = Rpq(n,2)*WTUVTEMP4(26,2)!not used
!!$     WTUVTEMP5(45,1) = WTUVTEMP3(20,2) + Rpq(n,1)*WTUVTEMP4(30,2)
!!$     WTUVTEMP5(46,1) = Rpq(n,1)*WTUVTEMP4(31,2)
!!$!     WTUVTEMP5(47,1) = Rpq(n,1)*WTUVTEMP4(32,2)!not used
!!$     WTUVTEMP5(48,1) = Rpq(n,1)*WTUVTEMP4(33,2)
!!$!     WTUVTEMP5(49,1) = Rpq(n,1)*WTUVTEMP4(34,2)!not used
!!$     WTUVTEMP5(50,1) = Rpq(n,1)*WTUVTEMP4(35,2)
!!$     WTUVTEMP5(51,1) = 4*WTUVTEMP3(17,2) + Rpq(n,2)*WTUVTEMP4(31,2)
!!$     WTUVTEMP5(52,1) = Rpq(n,3)*WTUVTEMP4(31,2)
!!$     WTUVTEMP5(53,1) = 2*WTUVTEMP3(19,2) + Rpq(n,2)*WTUVTEMP4(33,2)
!!$     WTUVTEMP5(54,1) = WTUVTEMP3(20,2) + Rpq(n,2)*WTUVTEMP4(34,2)
!!$     WTUVTEMP5(55,1) = Rpq(n,2)*WTUVTEMP4(35,2)
!!$     WTUVTEMP5(56,1) = 4*WTUVTEMP3(20,2) + Rpq(n,3)*WTUVTEMP4(35,2)
!!$
!!$     WTUV(n,57) = 5*WTUVTEMP4(21,1) + Rpq(n,1)*WTUVTEMP5(36,1)
!!$     WTUV(n,58) = Rpq(n,2)*WTUVTEMP5(36,1)
!!$     WTUV(n,59) = Rpq(n,3)*WTUVTEMP5(36,1)
!!$     WTUV(n,60) = 3*WTUVTEMP4(24,1) + Rpq(n,1)*WTUVTEMP5(39,1)
!!$     WTUV(n,61) = Rpq(n,2)*WTUVTEMP5(38,1)
!!$     WTUV(n,62) = 3*WTUVTEMP4(26,1) + Rpq(n,1)*WTUVTEMP5(41,1)
!!$     WTUV(n,63) = 2*WTUVTEMP4(27,1) + Rpq(n,1)*WTUVTEMP5(42,1)
!!$     WTUV(n,64) = Rpq(n,3)*WTUVTEMP5(39,1)
!!$     WTUV(n,65) = Rpq(n,2)*WTUVTEMP5(41,1)
!!$     WTUV(n,66) = 2*WTUVTEMP4(30,1) + Rpq(n,1)*WTUVTEMP5(45,1)
!!$     WTUV(n,67) = WTUVTEMP4(31,1) + Rpq(n,1)*WTUVTEMP5(46,1)
!!$     WTUV(n,68) = Rpq(n,3)*WTUVTEMP5(42,1)
!!$     WTUV(n,69) = WTUVTEMP4(33,1) + Rpq(n,1)*WTUVTEMP5(48,1)
!!$     WTUV(n,70) = Rpq(n,2)*WTUVTEMP5(45,1)
!!$     WTUV(n,71) = WTUVTEMP4(35,1) + Rpq(n,1)*WTUVTEMP5(50,1)
!!$     WTUV(n,72) = Rpq(n,1)*WTUVTEMP5(51,1)
!!$     WTUV(n,73) = Rpq(n,1)*WTUVTEMP5(52,1)
!!$     WTUV(n,74) = Rpq(n,1)*WTUVTEMP5(53,1)
!!$     WTUV(n,75) = Rpq(n,1)*WTUVTEMP5(54,1)
!!$     WTUV(n,76) = Rpq(n,1)*WTUVTEMP5(55,1)
!!$     WTUV(n,77) = Rpq(n,1)*WTUVTEMP5(56,1)
!!$     WTUV(n,78) = 5*WTUVTEMP4(31,1) + Rpq(n,2)*WTUVTEMP5(51,1)
!!$     WTUV(n,79) = Rpq(n,3)*WTUVTEMP5(51,1)
!!$     WTUV(n,80) = 3*WTUVTEMP4(33,1) + Rpq(n,2)*WTUVTEMP5(53,1)
!!$     WTUV(n,81) = 2*WTUVTEMP4(34,1) + Rpq(n,2)*WTUVTEMP5(54,1)
!!$     WTUV(n,82) = WTUVTEMP4(35,1) + Rpq(n,2)*WTUVTEMP5(55,1)
!!$     WTUV(n,83) = Rpq(n,2)*WTUVTEMP5(56,1)
!!$     WTUV(n,84) = 5*WTUVTEMP4(35,1) + Rpq(n,3)*WTUVTEMP5(56,1)
!!$  ENDDO
!!$END SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX6
!!$
!!$SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX7(WJ000,WTUV,Rpq,nPrim)
!!$  implicit none
!!$  Integer,intent(in)           :: nPrim
!!$  Real(realk),intent(in)       :: WJ000(0:7,nPrim),Rpq(nPrim,3)
!!$  Real(realk),intent(inout)    :: WTUV(nPrim,120)
!!$  !
!!$  Integer     :: n,M
!!$  Real(realk) :: WTUVTEMP(2:4,6)
!!$  Real(realk) :: WTUVTEMP2(5:10,5)
!!$  Real(realk) :: WTUVTEMP3(11:20,4)
!!$  Real(realk) :: WTUVTEMP4(21:35,3)
!!$  Real(realk) :: WTUVTEMP5(36:56,2)
!!$  Real(realk) :: WTUVTEMP6(57:84,1)
!!$  !000
!!$  DO n=1,nPrim
!!$     WTUV(n,1) = WJ000(0,n)
!!$  ENDDO
!!$  !it can reuse registers with this collection of this in 1 loop
!!$  DO n=1,nPrim
!!$     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
!!$     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
!!$     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
!!$     DO M=1,6
!!$        WTUVTEMP(2,M) = Rpq(n,1)*WJ000(M+1,n)
!!$        WTUVTEMP(3,M) = Rpq(n,2)*WJ000(M+1,n)
!!$        WTUVTEMP(4,M) = Rpq(n,3)*WJ000(M+1,n)     
!!$     ENDDO
!!$     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
!!$     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
!!$     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
!!$     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
!!$     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
!!$     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
!!$     DO M=1,5
!!$        WTUVTEMP2(5,M) = WJ000(M+1,n) + Rpq(n,1)*WTUVTEMP(2,M+1)
!!$!        WTUVTEMP2(6,M) =                Rpq(n,1)*WTUVTEMP(3,M+1)
!!$!        WTUVTEMP2(7,M) =                Rpq(n,1)*WTUVTEMP(4,M+1)
!!$        WTUVTEMP2(8,M) = WJ000(M+1,n) + Rpq(n,2)*WTUVTEMP(3,M+1)
!!$        WTUVTEMP2(9,M) =                Rpq(n,2)*WTUVTEMP(4,M+1)
!!$        WTUVTEMP2(10,M)= WJ000(M+1,n) + Rpq(n,3)*WTUVTEMP(4,M+1)     
!!$     ENDDO
!!$     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
!!$     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
!!$     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
!!$     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
!!$     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
!!$     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
!!$     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
!!$     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
!!$     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
!!$     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
!!$     DO M=1,4
!!$        WTUVTEMP3(11,M)= 2*WTUVTEMP(2,M+1) + Rpq(n,1)*WTUVTEMP2(5,M+1)
!!$!        WTUVTEMP3(12,M)=  Rpq(n,2)*WTUVTEMP2(5,M+1)
!!$        WTUVTEMP3(13,M)= Rpq(n,3)*WTUVTEMP2(5,M+1)
!!$        WTUVTEMP3(14,M)= Rpq(n,1)*WTUVTEMP2(8,M+1)
!!$!        WTUVTEMP3(15,M)= Rpq(n,1)*WTUVTEMP2(9,M+1)
!!$        WTUVTEMP3(16,M)= Rpq(n,1)*WTUVTEMP2(10,M+1)
!!$        WTUVTEMP3(17,M)= 2*WTUVTEMP(3,M+1) + Rpq(n,2)*WTUVTEMP2(8,M+1)
!!$        WTUVTEMP3(18,M)= Rpq(n,3)*WTUVTEMP2(8,M+1)
!!$        WTUVTEMP3(19,M)= Rpq(n,2)*WTUVTEMP2(10,M+1)
!!$        WTUVTEMP3(20,M)= 2*WTUVTEMP(4,M+1) + Rpq(n,3)*WTUVTEMP2(10,M+1)
!!$     ENDDO
!!$
!!$     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11,1)
!!$     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11,1)
!!$     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11,1)
!!$     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14,1)
!!$     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13,1)
!!$     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16,1)
!!$     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17,1)
!!$     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18,1)
!!$     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19,1)
!!$     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20,1)
!!$     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17,1) 
!!$     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17,1)
!!$     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19,1)
!!$     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20,1)
!!$     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20,1)
!!$
!!$     DO M=1,3
!!$        WTUVTEMP4(21,M)= 3*WTUVTEMP2(5,M+1) + Rpq(n,1)*WTUVTEMP3(11,M+1)
!!$        !     WTUVTEMP4(22,M)= Rpq(n,2)*WTUVTEMP3(11,M+1)
!!$        WTUVTEMP4(23,M)= Rpq(n,3)*WTUVTEMP3(11,M+1)
!!$        WTUVTEMP4(24,M)= WTUVTEMP2(8,M+1) + Rpq(n,1)*WTUVTEMP3(14,M+1)
!!$        !     WTUVTEMP4(25,M)= Rpq(n,2)*WTUVTEMP3(13,M+1)
!!$        WTUVTEMP4(26,M)= WTUVTEMP2(10,M+1) + Rpq(n,1)*WTUVTEMP3(16,M+1)
!!$        WTUVTEMP4(27,M)= Rpq(n,1)*WTUVTEMP3(17,M+1)
!!$        !     WTUVTEMP4(28,M)= Rpq(n,1)*WTUVTEMP3(18,M+1)
!!$        !     WTUVTEMP4(29,M)= Rpq(n,1)*WTUVTEMP3(19,M+1)
!!$        WTUVTEMP4(30,M)= Rpq(n,1)*WTUVTEMP3(20,M+1)
!!$        WTUVTEMP4(31,M)= 3*WTUVTEMP2(8,M+1) + Rpq(n,2)*WTUVTEMP3(17,M+1) 
!!$        WTUVTEMP4(32,M)= Rpq(n,3)*WTUVTEMP3(17,M+1)
!!$        WTUVTEMP4(33,M)= WTUVTEMP2(10,M+1) + Rpq(n,2)*WTUVTEMP3(19,M+1)
!!$        WTUVTEMP4(34,M)= Rpq(n,2)*WTUVTEMP3(20,M+1)
!!$        WTUVTEMP4(35,M)= 3*WTUVTEMP2(10,M+1) + Rpq(n,3)*WTUVTEMP3(20,M+1)
!!$     ENDDO
!!$
!!$     WTUV(n,  36) = 4*WTUVTEMP3(11,1) + Rpq(n,1)*WTUVTEMP4(21,1)
!!$     WTUV(n,  37) = Rpq(n,2)*WTUVTEMP4(21,1)
!!$     WTUV(n,  38) = Rpq(n,3)*WTUVTEMP4(21,1)
!!$     WTUV(n,  39) = 2*WTUVTEMP3(14,1) + Rpq(n,1)*WTUVTEMP4(24,1)
!!$     WTUV(n,  40) = Rpq(n,2)*WTUVTEMP4(23,1)
!!$     WTUV(n,  41) = 2*WTUVTEMP3(16,1) + Rpq(n,1)*WTUVTEMP4(26,1)
!!$     WTUV(n,  42) = WTUVTEMP3(17,1) + Rpq(n,1)*WTUVTEMP4(27,1)
!!$     WTUV(n,  43) = Rpq(n,3)*WTUVTEMP4(24,1)
!!$     WTUV(n,  44) = Rpq(n,2)*WTUVTEMP4(26,1)
!!$     WTUV(n,  45) = WTUVTEMP3(20,1) + Rpq(n,1)*WTUVTEMP4(30,1)
!!$     WTUV(n,  46) = Rpq(n,1)*WTUVTEMP4(31,1)
!!$     WTUV(n,  47) = Rpq(n,1)*WTUVTEMP4(32,1)
!!$     WTUV(n,  48) = Rpq(n,1)*WTUVTEMP4(33,1)
!!$     WTUV(n,  49) = Rpq(n,1)*WTUVTEMP4(34,1)
!!$     WTUV(n,  50) = Rpq(n,1)*WTUVTEMP4(35,1)
!!$     WTUV(n,  51) = 4*WTUVTEMP3(17,1) + Rpq(n,2)*WTUVTEMP4(31,1)
!!$     WTUV(n,  52) = Rpq(n,3)*WTUVTEMP4(31,1)
!!$     WTUV(n,  53) = 2*WTUVTEMP3(19,1) + Rpq(n,2)*WTUVTEMP4(33,1)
!!$     WTUV(n,  54) = WTUVTEMP3(20,1) + Rpq(n,2)*WTUVTEMP4(34,1)
!!$     WTUV(n,  55) = Rpq(n,2)*WTUVTEMP4(35,1)
!!$     WTUV(n,  56) = 4*WTUVTEMP3(20,1) + Rpq(n,3)*WTUVTEMP4(35,1)
!!$
!!$     DO M=1,2
!!$        WTUVTEMP5(36,M) = 4*WTUVTEMP3(11,M+1) + Rpq(n,1)*WTUVTEMP4(21,M+1)
!!$        !     WTUVTEMP5(37,M) = Rpq(n,2)*WTUVTEMP4(21,M+1)
!!$        WTUVTEMP5(38,M) = Rpq(n,3)*WTUVTEMP4(21,M+1)
!!$        WTUVTEMP5(39,M) = 2*WTUVTEMP3(14,M+1) + Rpq(n,1)*WTUVTEMP4(24,M+1)
!!$        !     WTUVTEMP5(40,M) = Rpq(n,2)*WTUVTEMP4(23,M+1)
!!$        WTUVTEMP5(41,M) = 2*WTUVTEMP3(16,M+1) + Rpq(n,1)*WTUVTEMP4(26,M+1)
!!$        WTUVTEMP5(42,M) = WTUVTEMP3(17,M+1) + Rpq(n,1)*WTUVTEMP4(27,M+1)
!!$        !     WTUVTEMP5(43,M) = Rpq(n,3)*WTUVTEMP4(24,M+1)
!!$        !     WTUVTEMP5(44,M) = Rpq(n,2)*WTUVTEMP4(26,M+1)
!!$        WTUVTEMP5(45,M) = WTUVTEMP3(20,M+1) + Rpq(n,1)*WTUVTEMP4(30,M+1)
!!$        WTUVTEMP5(46,M) = Rpq(n,1)*WTUVTEMP4(31,M+1)
!!$        !     WTUVTEMP5(47,M) = Rpq(n,1)*WTUVTEMP4(32,M+1)
!!$        WTUVTEMP5(48,M) = Rpq(n,1)*WTUVTEMP4(33,M+1)
!!$        !     WTUVTEMP5(49,M) = Rpq(n,1)*WTUVTEMP4(34,M+1)
!!$        WTUVTEMP5(50,M) = Rpq(n,1)*WTUVTEMP4(35,M+1)
!!$        WTUVTEMP5(51,M) = 4*WTUVTEMP3(17,M+1) + Rpq(n,2)*WTUVTEMP4(31,M+1)
!!$        WTUVTEMP5(52,M) = Rpq(n,3)*WTUVTEMP4(31,M+1)
!!$        WTUVTEMP5(53,M) = 2*WTUVTEMP3(19,M+1) + Rpq(n,2)*WTUVTEMP4(33,M+1)
!!$        WTUVTEMP5(54,M) = WTUVTEMP3(20,M+1) + Rpq(n,2)*WTUVTEMP4(34,M+1)
!!$        WTUVTEMP5(55,M) = Rpq(n,2)*WTUVTEMP4(35,M+1)
!!$        WTUVTEMP5(56,M) = 4*WTUVTEMP3(20,M+1) + Rpq(n,3)*WTUVTEMP4(35,M+1)
!!$     ENDDO
!!$     WTUV(n,57) = 5*WTUVTEMP4(21,1) + Rpq(n,1)*WTUVTEMP5(36,1)
!!$     WTUV(n,58) = Rpq(n,2)*WTUVTEMP5(36,1)
!!$     WTUV(n,59) = Rpq(n,3)*WTUVTEMP5(36,1)
!!$     WTUV(n,60) = 3*WTUVTEMP4(24,1) + Rpq(n,1)*WTUVTEMP5(39,1)
!!$     WTUV(n,61) = Rpq(n,2)*WTUVTEMP5(38,1)
!!$     WTUV(n,62) = 3*WTUVTEMP4(26,1) + Rpq(n,1)*WTUVTEMP5(41,1)
!!$     WTUV(n,63) = 2*WTUVTEMP4(27,1) + Rpq(n,1)*WTUVTEMP5(42,1)
!!$     WTUV(n,64) = Rpq(n,3)*WTUVTEMP5(39,1)
!!$     WTUV(n,65) = Rpq(n,2)*WTUVTEMP5(41,1)
!!$     WTUV(n,66) = 2*WTUVTEMP4(30,1) + Rpq(n,1)*WTUVTEMP5(45,1)
!!$     WTUV(n,67) = WTUVTEMP4(31,1) + Rpq(n,1)*WTUVTEMP5(46,1)
!!$     WTUV(n,68) = Rpq(n,3)*WTUVTEMP5(42,1)
!!$     WTUV(n,69) = WTUVTEMP4(33,1) + Rpq(n,1)*WTUVTEMP5(48,1)
!!$     WTUV(n,70) = Rpq(n,2)*WTUVTEMP5(45,1)
!!$     WTUV(n,71) = WTUVTEMP4(35,1) + Rpq(n,1)*WTUVTEMP5(50,1)
!!$     WTUV(n,72) = Rpq(n,1)*WTUVTEMP5(51,1)
!!$     WTUV(n,73) = Rpq(n,1)*WTUVTEMP5(52,1)
!!$     WTUV(n,74) = Rpq(n,1)*WTUVTEMP5(53,1)
!!$     WTUV(n,75) = Rpq(n,1)*WTUVTEMP5(54,1)
!!$     WTUV(n,76) = Rpq(n,1)*WTUVTEMP5(55,1)
!!$     WTUV(n,77) = Rpq(n,1)*WTUVTEMP5(56,1)
!!$     WTUV(n,78) = 5*WTUVTEMP4(31,1) + Rpq(n,2)*WTUVTEMP5(51,1)
!!$     WTUV(n,79) = Rpq(n,3)*WTUVTEMP5(51,1)
!!$     WTUV(n,80) = 3*WTUVTEMP4(33,1) + Rpq(n,2)*WTUVTEMP5(53,1)
!!$     WTUV(n,81) = 2*WTUVTEMP4(34,1) + Rpq(n,2)*WTUVTEMP5(54,1)
!!$     WTUV(n,82) = WTUVTEMP4(35,1) + Rpq(n,2)*WTUVTEMP5(55,1)
!!$     WTUV(n,83) = Rpq(n,2)*WTUVTEMP5(56,1)
!!$     WTUV(n,84) = 5*WTUVTEMP4(35,1) + Rpq(n,3)*WTUVTEMP5(56,1)
!!$     
!!$     WTUVTEMP6(57,1) = 5*WTUVTEMP4(21,2) + Rpq(n,1)*WTUVTEMP5(36,2)
!!$!     WTUVTEMP6(58,1) = Rpq(n,2)*WTUVTEMP5(36,2)
!!$     WTUVTEMP6(59,1) = Rpq(n,3)*WTUVTEMP5(36,2)
!!$     WTUVTEMP6(60,1) = 3*WTUVTEMP4(24,2) + Rpq(n,1)*WTUVTEMP5(39,2)
!!$!     WTUVTEMP6(61,1) = Rpq(n,2)*WTUVTEMP5(38,2)
!!$     WTUVTEMP6(62,1) = 3*WTUVTEMP4(26,2) + Rpq(n,1)*WTUVTEMP5(41,2)
!!$     WTUVTEMP6(63,1) = 2*WTUVTEMP4(27,2) + Rpq(n,1)*WTUVTEMP5(42,2)
!!$!     WTUVTEMP6(64,1) = Rpq(n,3)*WTUVTEMP5(39,2)
!!$!     WTUVTEMP6(65,1) = Rpq(n,2)*WTUVTEMP5(41,2)
!!$     WTUVTEMP6(66,1) = 2*WTUVTEMP4(30,2) + Rpq(n,1)*WTUVTEMP5(45,2)
!!$     WTUVTEMP6(67,1) = WTUVTEMP4(31,2) + Rpq(n,1)*WTUVTEMP5(46,2)
!!$!     WTUVTEMP6(68,1) = Rpq(n,3)*WTUVTEMP5(42,2)
!!$     WTUVTEMP6(69,1) = WTUVTEMP4(33,2) + Rpq(n,1)*WTUVTEMP5(48,2)
!!$!     WTUVTEMP6(70,1) = Rpq(n,2)*WTUVTEMP5(45,2)
!!$     WTUVTEMP6(71,1) = WTUVTEMP4(35,2) + Rpq(n,1)*WTUVTEMP5(50,2)
!!$     WTUVTEMP6(72,1) = Rpq(n,1)*WTUVTEMP5(51,2)
!!$!     WTUVTEMP6(73,1) = Rpq(n,1)*WTUVTEMP5(52,2)
!!$     WTUVTEMP6(74,1) = Rpq(n,1)*WTUVTEMP5(53,2)
!!$     WTUVTEMP6(75,1) = Rpq(n,1)*WTUVTEMP5(54,2)
!!$!     WTUVTEMP6(76,1) = Rpq(n,1)*WTUVTEMP5(55,2)
!!$     WTUVTEMP6(77,1) = Rpq(n,1)*WTUVTEMP5(56,2)
!!$     WTUVTEMP6(78,1) = 5*WTUVTEMP4(31,2) + Rpq(n,2)*WTUVTEMP5(51,2)
!!$     WTUVTEMP6(79,1) = Rpq(n,3)*WTUVTEMP5(51,2)
!!$     WTUVTEMP6(80,1) = 3*WTUVTEMP4(33,2) + Rpq(n,2)*WTUVTEMP5(53,2)
!!$     WTUVTEMP6(81,1) = 2*WTUVTEMP4(34,2) + Rpq(n,2)*WTUVTEMP5(54,2)
!!$     WTUVTEMP6(82,1) = WTUVTEMP4(35,2) + Rpq(n,2)*WTUVTEMP5(55,2)
!!$     WTUVTEMP6(83,1) = Rpq(n,2)*WTUVTEMP5(56,2)
!!$     WTUVTEMP6(84,1) = 5*WTUVTEMP4(35,2) + Rpq(n,3)*WTUVTEMP5(56,2)
!!$
!!$     WTUV(n, 85) = 6*WTUVTEMP5(36,1) + Rpq(n,1)*WTUVTEMP6(57,1)
!!$     WTUV(n, 86) = Rpq(n,2)*WTUVTEMP6(57,1)
!!$     WTUV(n, 87) = Rpq(n,3)*WTUVTEMP6(57,1)
!!$     WTUV(n, 88) = 4*WTUVTEMP5(39,1) + Rpq(n,1)*WTUVTEMP6(60,1)
!!$     WTUV(n, 89) = Rpq(n,2)*WTUVTEMP6(59,1)
!!$     WTUV(n, 90) = 4*WTUVTEMP5(41,1) + Rpq(n,1)*WTUVTEMP6(62,1)
!!$     WTUV(n, 91) = 3*WTUVTEMP5(42,1) + Rpq(n,1)*WTUVTEMP6(63,1)
!!$     WTUV(n, 92) = Rpq(n,3)*WTUVTEMP6(60,1)
!!$     WTUV(n, 93) = Rpq(n,2)*WTUVTEMP6(62,1)
!!$     WTUV(n, 94) = 3*WTUVTEMP5(45,1) + Rpq(n,1)*WTUVTEMP6(66,1)
!!$     WTUV(n, 95) = 2*WTUVTEMP5(46,1) + Rpq(n,1)*WTUVTEMP6(67,1)
!!$     WTUV(n, 96) = Rpq(n,3)*WTUVTEMP6(63,1)
!!$     WTUV(n, 97) = 2*WTUVTEMP5(48,1) + Rpq(n,1)*WTUVTEMP6(69,1)
!!$     WTUV(n, 98) = Rpq(n,2)*WTUVTEMP6(66,1)
!!$     WTUV(n, 99) = 2*WTUVTEMP5(50,1) + Rpq(n,1)*WTUVTEMP6(71,1)
!!$     WTUV(n,100) = WTUVTEMP5(51,1) + Rpq(n,1)*WTUVTEMP6(72,1)
!!$     WTUV(n,101) = Rpq(n,3)*WTUVTEMP6(67,1)
!!$     WTUV(n,102) = WTUVTEMP5(53,1) + Rpq(n,1)*WTUVTEMP6(74,1)
!!$     WTUV(n,103) = WTUVTEMP5(54,1) + Rpq(n,1)*WTUVTEMP6(75,1)
!!$     WTUV(n,104) = Rpq(n,2)*WTUVTEMP6(71,1)
!!$     WTUV(n,105) = WTUVTEMP5(56,1) + Rpq(n,1)*WTUVTEMP6(77,1)
!!$     WTUV(n,106) = Rpq(n,1)*WTUVTEMP6(78,1)
!!$     WTUV(n,107) = Rpq(n,1)*WTUVTEMP6(79,1)
!!$     WTUV(n,108) = Rpq(n,1)*WTUVTEMP6(80,1)
!!$     WTUV(n,109) = Rpq(n,1)*WTUVTEMP6(81,1)
!!$     WTUV(n,110) = Rpq(n,1)*WTUVTEMP6(82,1)
!!$     WTUV(n,111) = Rpq(n,1)*WTUVTEMP6(83,1)
!!$     WTUV(n,112) = Rpq(n,1)*WTUVTEMP6(84,1)
!!$     WTUV(n,113) = 6*WTUVTEMP5(51,1) + Rpq(n,2)*WTUVTEMP6(78,1)
!!$     WTUV(n,114) = Rpq(n,3)*WTUVTEMP6(78,1)
!!$     WTUV(n,115) = 4*WTUVTEMP5(53,1) + Rpq(n,2)*WTUVTEMP6(80,1)
!!$     WTUV(n,116) = 3*WTUVTEMP5(54,1) + Rpq(n,2)*WTUVTEMP6(81,1)
!!$     WTUV(n,117) = 2*WTUVTEMP5(55,1) + Rpq(n,2)*WTUVTEMP6(82,1)
!!$     WTUV(n,118) = WTUVTEMP5(56,1) + Rpq(n,2)*WTUVTEMP6(83,1)
!!$     WTUV(n,119) = Rpq(n,2)*WTUVTEMP6(84,1)
!!$     WTUV(n,120) = 6*WTUVTEMP5(56,1) + Rpq(n,3)*WTUVTEMP6(84,1)
!!$  ENDDO
!!$END SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX7
!!$
!!$SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX8(WJ000,WTUV,Rpq,nPrim)
!!$  implicit none
!!$  Integer,intent(in)           :: nPrim
!!$  Real(realk),intent(in)       :: WJ000(0:8,nPrim),Rpq(nPrim,3)
!!$  Real(realk),intent(inout)    :: WTUV(nPrim,165)
!!$  !
!!$  Integer     :: n,M
!!$  Real(realk) :: WTUVTEMP(2:4,7)
!!$  Real(realk) :: WTUVTEMP2(5:10,6)
!!$  Real(realk) :: WTUVTEMP3(11:20,5)
!!$  Real(realk) :: WTUVTEMP4(21:35,4)
!!$  Real(realk) :: WTUVTEMP5(36:56,3)
!!$  Real(realk) :: WTUVTEMP6(57:84,2)
!!$  Real(realk) :: WTUVTEMP7(85:120,1)
!!$  !000
!!$  DO n=1,nPrim
!!$     WTUV(n,1) = WJ000(0,n)
!!$  ENDDO
!!$  !it can reuse registers with this collection of this in 1 loop
!!$  DO n=1,nPrim
!!$     WTUV(n,2) = Rpq(n,1)*WJ000(1,n)
!!$     WTUV(n,3) = Rpq(n,2)*WJ000(1,n)
!!$     WTUV(n,4) = Rpq(n,3)*WJ000(1,n)
!!$     DO M=1,7
!!$        WTUVTEMP(2,M) = Rpq(n,1)*WJ000(M+1,n)
!!$        WTUVTEMP(3,M) = Rpq(n,2)*WJ000(M+1,n)
!!$        WTUVTEMP(4,M) = Rpq(n,3)*WJ000(M+1,n)     
!!$     ENDDO
!!$     WTUV(n,5) = WJ000(1,n) + Rpq(n,1)*WTUVTEMP(2,1)
!!$     WTUV(n,6) =              Rpq(n,1)*WTUVTEMP(3,1)
!!$     WTUV(n,7) =              Rpq(n,1)*WTUVTEMP(4,1)
!!$     WTUV(n,8) = WJ000(1,n) + Rpq(n,2)*WTUVTEMP(3,1)
!!$     WTUV(n,9) =              Rpq(n,2)*WTUVTEMP(4,1)
!!$     WTUV(n,10)= WJ000(1,n) + Rpq(n,3)*WTUVTEMP(4,1)
!!$     DO M=1,6
!!$        WTUVTEMP2(5,M) = WJ000(M+1,n) + Rpq(n,1)*WTUVTEMP(2,M+1)
!!$!        WTUVTEMP2(6,M) =                Rpq(n,1)*WTUVTEMP(3,M+1)
!!$!        WTUVTEMP2(7,M) =                Rpq(n,1)*WTUVTEMP(4,M+1)
!!$        WTUVTEMP2(8,M) = WJ000(M+1,n) + Rpq(n,2)*WTUVTEMP(3,M+1)
!!$        WTUVTEMP2(9,M) =                Rpq(n,2)*WTUVTEMP(4,M+1)
!!$        WTUVTEMP2(10,M)= WJ000(M+1,n) + Rpq(n,3)*WTUVTEMP(4,M+1)     
!!$     ENDDO
!!$     WTUV(n,11)= 2*WTUVTEMP(2,1) + Rpq(n,1)*WTUVTEMP2(5,1)
!!$     WTUV(n,12)=  Rpq(n,2)*WTUVTEMP2(5,1)
!!$     WTUV(n,13)= Rpq(n,3)*WTUVTEMP2(5,1)
!!$     WTUV(n,14)= Rpq(n,1)*WTUVTEMP2(8,1)
!!$     WTUV(n,15)= Rpq(n,1)*WTUVTEMP2(9,1)
!!$     WTUV(n,16)= Rpq(n,1)*WTUVTEMP2(10,1)
!!$     WTUV(n,17)= 2*WTUVTEMP(3,1) + Rpq(n,2)*WTUVTEMP2(8,1)
!!$     WTUV(n,18)= Rpq(n,3)*WTUVTEMP2(8,1)
!!$     WTUV(n,19)= Rpq(n,2)*WTUVTEMP2(10,1)
!!$     WTUV(n,20)= 2*WTUVTEMP(4,1) + Rpq(n,3)*WTUVTEMP2(10,1)
!!$     DO M=1,5
!!$        WTUVTEMP3(11,M)= 2*WTUVTEMP(2,M+1) + Rpq(n,1)*WTUVTEMP2(5,M+1)
!!$!        WTUVTEMP3(12,M)=  Rpq(n,2)*WTUVTEMP2(5,M+1)
!!$        WTUVTEMP3(13,M)= Rpq(n,3)*WTUVTEMP2(5,M+1)
!!$        WTUVTEMP3(14,M)= Rpq(n,1)*WTUVTEMP2(8,M+1)
!!$!        WTUVTEMP3(15,M)= Rpq(n,1)*WTUVTEMP2(9,M+1)
!!$        WTUVTEMP3(16,M)= Rpq(n,1)*WTUVTEMP2(10,M+1)
!!$        WTUVTEMP3(17,M)= 2*WTUVTEMP(3,M+1) + Rpq(n,2)*WTUVTEMP2(8,M+1)
!!$        WTUVTEMP3(18,M)= Rpq(n,3)*WTUVTEMP2(8,M+1)
!!$        WTUVTEMP3(19,M)= Rpq(n,2)*WTUVTEMP2(10,M+1)
!!$        WTUVTEMP3(20,M)= 2*WTUVTEMP(4,M+1) + Rpq(n,3)*WTUVTEMP2(10,M+1)
!!$     ENDDO
!!$
!!$     WTUV(n,21)= 3*WTUVTEMP2(5,1) + Rpq(n,1)*WTUVTEMP3(11,1)
!!$     WTUV(n,22)= Rpq(n,2)*WTUVTEMP3(11,1)
!!$     WTUV(n,23)= Rpq(n,3)*WTUVTEMP3(11,1)
!!$     WTUV(n,24)= WTUVTEMP2(8,1) + Rpq(n,1)*WTUVTEMP3(14,1)
!!$     WTUV(n,25)= Rpq(n,2)*WTUVTEMP3(13,1)
!!$     WTUV(n,26)= WTUVTEMP2(10,1) + Rpq(n,1)*WTUVTEMP3(16,1)
!!$     WTUV(n,27)= Rpq(n,1)*WTUVTEMP3(17,1)
!!$     WTUV(n,28)= Rpq(n,1)*WTUVTEMP3(18,1)
!!$     WTUV(n,29)= Rpq(n,1)*WTUVTEMP3(19,1)
!!$     WTUV(n,30)= Rpq(n,1)*WTUVTEMP3(20,1)
!!$     WTUV(n,31)= 3*WTUVTEMP2(8,1) + Rpq(n,2)*WTUVTEMP3(17,1) 
!!$     WTUV(n,32)= Rpq(n,3)*WTUVTEMP3(17,1)
!!$     WTUV(n,33)= WTUVTEMP2(10,1) + Rpq(n,2)*WTUVTEMP3(19,1)
!!$     WTUV(n,34)= Rpq(n,2)*WTUVTEMP3(20,1)
!!$     WTUV(n,35)= 3*WTUVTEMP2(10,1) + Rpq(n,3)*WTUVTEMP3(20,1)
!!$
!!$     DO M=1,4
!!$        WTUVTEMP4(21,M)= 3*WTUVTEMP2(5,M+1) + Rpq(n,1)*WTUVTEMP3(11,M+1)
!!$        !     WTUVTEMP4(22,M)= Rpq(n,2)*WTUVTEMP3(11,M+1)
!!$        WTUVTEMP4(23,M)= Rpq(n,3)*WTUVTEMP3(11,M+1)
!!$        WTUVTEMP4(24,M)= WTUVTEMP2(8,M+1) + Rpq(n,1)*WTUVTEMP3(14,M+1)
!!$        !     WTUVTEMP4(25,M)= Rpq(n,2)*WTUVTEMP3(13,M+1)
!!$        WTUVTEMP4(26,M)= WTUVTEMP2(10,M+1) + Rpq(n,1)*WTUVTEMP3(16,M+1)
!!$        WTUVTEMP4(27,M)= Rpq(n,1)*WTUVTEMP3(17,M+1)
!!$        !     WTUVTEMP4(28,M)= Rpq(n,1)*WTUVTEMP3(18,M+1)
!!$        !     WTUVTEMP4(29,M)= Rpq(n,1)*WTUVTEMP3(19,M+1)
!!$        WTUVTEMP4(30,M)= Rpq(n,1)*WTUVTEMP3(20,M+1)
!!$        WTUVTEMP4(31,M)= 3*WTUVTEMP2(8,M+1) + Rpq(n,2)*WTUVTEMP3(17,M+1) 
!!$        WTUVTEMP4(32,M)= Rpq(n,3)*WTUVTEMP3(17,M+1)
!!$        WTUVTEMP4(33,M)= WTUVTEMP2(10,M+1) + Rpq(n,2)*WTUVTEMP3(19,M+1)
!!$        WTUVTEMP4(34,M)= Rpq(n,2)*WTUVTEMP3(20,M+1)
!!$        WTUVTEMP4(35,M)= 3*WTUVTEMP2(10,M+1) + Rpq(n,3)*WTUVTEMP3(20,M+1)
!!$     ENDDO
!!$
!!$     WTUV(n,  36) = 4*WTUVTEMP3(11,1) + Rpq(n,1)*WTUVTEMP4(21,1)
!!$     WTUV(n,  37) = Rpq(n,2)*WTUVTEMP4(21,1)
!!$     WTUV(n,  38) = Rpq(n,3)*WTUVTEMP4(21,1)
!!$     WTUV(n,  39) = 2*WTUVTEMP3(14,1) + Rpq(n,1)*WTUVTEMP4(24,1)
!!$     WTUV(n,  40) = Rpq(n,2)*WTUVTEMP4(23,1)
!!$     WTUV(n,  41) = 2*WTUVTEMP3(16,1) + Rpq(n,1)*WTUVTEMP4(26,1)
!!$     WTUV(n,  42) = WTUVTEMP3(17,1) + Rpq(n,1)*WTUVTEMP4(27,1)
!!$     WTUV(n,  43) = Rpq(n,3)*WTUVTEMP4(24,1)
!!$     WTUV(n,  44) = Rpq(n,2)*WTUVTEMP4(26,1)
!!$     WTUV(n,  45) = WTUVTEMP3(20,1) + Rpq(n,1)*WTUVTEMP4(30,1)
!!$     WTUV(n,  46) = Rpq(n,1)*WTUVTEMP4(31,1)
!!$     WTUV(n,  47) = Rpq(n,1)*WTUVTEMP4(32,1)
!!$     WTUV(n,  48) = Rpq(n,1)*WTUVTEMP4(33,1)
!!$     WTUV(n,  49) = Rpq(n,1)*WTUVTEMP4(34,1)
!!$     WTUV(n,  50) = Rpq(n,1)*WTUVTEMP4(35,1)
!!$     WTUV(n,  51) = 4*WTUVTEMP3(17,1) + Rpq(n,2)*WTUVTEMP4(31,1)
!!$     WTUV(n,  52) = Rpq(n,3)*WTUVTEMP4(31,1)
!!$     WTUV(n,  53) = 2*WTUVTEMP3(19,1) + Rpq(n,2)*WTUVTEMP4(33,1)
!!$     WTUV(n,  54) = WTUVTEMP3(20,1) + Rpq(n,2)*WTUVTEMP4(34,1)
!!$     WTUV(n,  55) = Rpq(n,2)*WTUVTEMP4(35,1)
!!$     WTUV(n,  56) = 4*WTUVTEMP3(20,1) + Rpq(n,3)*WTUVTEMP4(35,1)
!!$
!!$     DO M=1,3
!!$        WTUVTEMP5(36,M) = 4*WTUVTEMP3(11,M+1) + Rpq(n,1)*WTUVTEMP4(21,M+1)
!!$        !     WTUVTEMP5(37,M) = Rpq(n,2)*WTUVTEMP4(21,M+1)
!!$        WTUVTEMP5(38,M) = Rpq(n,3)*WTUVTEMP4(21,M+1)
!!$        WTUVTEMP5(39,M) = 2*WTUVTEMP3(14,M+1) + Rpq(n,1)*WTUVTEMP4(24,M+1)
!!$        !     WTUVTEMP5(40,M) = Rpq(n,2)*WTUVTEMP4(23,M+1)
!!$        WTUVTEMP5(41,M) = 2*WTUVTEMP3(16,M+1) + Rpq(n,1)*WTUVTEMP4(26,M+1)
!!$        WTUVTEMP5(42,M) = WTUVTEMP3(17,M+1) + Rpq(n,1)*WTUVTEMP4(27,M+1)
!!$        !     WTUVTEMP5(43,M) = Rpq(n,3)*WTUVTEMP4(24,M+1)
!!$        !     WTUVTEMP5(44,M) = Rpq(n,2)*WTUVTEMP4(26,M+1)
!!$        WTUVTEMP5(45,M) = WTUVTEMP3(20,M+1) + Rpq(n,1)*WTUVTEMP4(30,M+1)
!!$        WTUVTEMP5(46,M) = Rpq(n,1)*WTUVTEMP4(31,M+1)
!!$        !     WTUVTEMP5(47,M) = Rpq(n,1)*WTUVTEMP4(32,M+1)
!!$        WTUVTEMP5(48,M) = Rpq(n,1)*WTUVTEMP4(33,M+1)
!!$        !     WTUVTEMP5(49,M) = Rpq(n,1)*WTUVTEMP4(34,M+1)
!!$        WTUVTEMP5(50,M) = Rpq(n,1)*WTUVTEMP4(35,M+1)
!!$        WTUVTEMP5(51,M) = 4*WTUVTEMP3(17,M+1) + Rpq(n,2)*WTUVTEMP4(31,M+1)
!!$        WTUVTEMP5(52,M) = Rpq(n,3)*WTUVTEMP4(31,M+1)
!!$        WTUVTEMP5(53,M) = 2*WTUVTEMP3(19,M+1) + Rpq(n,2)*WTUVTEMP4(33,M+1)
!!$        WTUVTEMP5(54,M) = WTUVTEMP3(20,M+1) + Rpq(n,2)*WTUVTEMP4(34,M+1)
!!$        WTUVTEMP5(55,M) = Rpq(n,2)*WTUVTEMP4(35,M+1)
!!$        WTUVTEMP5(56,M) = 4*WTUVTEMP3(20,M+1) + Rpq(n,3)*WTUVTEMP4(35,M+1)
!!$     ENDDO
!!$     WTUV(n,57) = 5*WTUVTEMP4(21,1) + Rpq(n,1)*WTUVTEMP5(36,1)
!!$     WTUV(n,58) = Rpq(n,2)*WTUVTEMP5(36,1)
!!$     WTUV(n,59) = Rpq(n,3)*WTUVTEMP5(36,1)
!!$     WTUV(n,60) = 3*WTUVTEMP4(24,1) + Rpq(n,1)*WTUVTEMP5(39,1)
!!$     WTUV(n,61) = Rpq(n,2)*WTUVTEMP5(38,1)
!!$     WTUV(n,62) = 3*WTUVTEMP4(26,1) + Rpq(n,1)*WTUVTEMP5(41,1)
!!$     WTUV(n,63) = 2*WTUVTEMP4(27,1) + Rpq(n,1)*WTUVTEMP5(42,1)
!!$     WTUV(n,64) = Rpq(n,3)*WTUVTEMP5(39,1)
!!$     WTUV(n,65) = Rpq(n,2)*WTUVTEMP5(41,1)
!!$     WTUV(n,66) = 2*WTUVTEMP4(30,1) + Rpq(n,1)*WTUVTEMP5(45,1)
!!$     WTUV(n,67) = WTUVTEMP4(31,1) + Rpq(n,1)*WTUVTEMP5(46,1)
!!$     WTUV(n,68) = Rpq(n,3)*WTUVTEMP5(42,1)
!!$     WTUV(n,69) = WTUVTEMP4(33,1) + Rpq(n,1)*WTUVTEMP5(48,1)
!!$     WTUV(n,70) = Rpq(n,2)*WTUVTEMP5(45,1)
!!$     WTUV(n,71) = WTUVTEMP4(35,1) + Rpq(n,1)*WTUVTEMP5(50,1)
!!$     WTUV(n,72) = Rpq(n,1)*WTUVTEMP5(51,1)
!!$     WTUV(n,73) = Rpq(n,1)*WTUVTEMP5(52,1)
!!$     WTUV(n,74) = Rpq(n,1)*WTUVTEMP5(53,1)
!!$     WTUV(n,75) = Rpq(n,1)*WTUVTEMP5(54,1)
!!$     WTUV(n,76) = Rpq(n,1)*WTUVTEMP5(55,1)
!!$     WTUV(n,77) = Rpq(n,1)*WTUVTEMP5(56,1)
!!$     WTUV(n,78) = 5*WTUVTEMP4(31,1) + Rpq(n,2)*WTUVTEMP5(51,1)
!!$     WTUV(n,79) = Rpq(n,3)*WTUVTEMP5(51,1)
!!$     WTUV(n,80) = 3*WTUVTEMP4(33,1) + Rpq(n,2)*WTUVTEMP5(53,1)
!!$     WTUV(n,81) = 2*WTUVTEMP4(34,1) + Rpq(n,2)*WTUVTEMP5(54,1)
!!$     WTUV(n,82) = WTUVTEMP4(35,1) + Rpq(n,2)*WTUVTEMP5(55,1)
!!$     WTUV(n,83) = Rpq(n,2)*WTUVTEMP5(56,1)
!!$     WTUV(n,84) = 5*WTUVTEMP4(35,1) + Rpq(n,3)*WTUVTEMP5(56,1)
!!$     
!!$     DO M=1,2
!!$        WTUVTEMP6(57,M) = 5*WTUVTEMP4(21,M+1) + Rpq(n,1)*WTUVTEMP5(36,M+1)
!!$        !     WTUVTEMP6(58,M) = Rpq(n,2)*WTUVTEMP5(36,M+1)
!!$        WTUVTEMP6(59,M) = Rpq(n,3)*WTUVTEMP5(36,M+1)
!!$        WTUVTEMP6(60,M) = 3*WTUVTEMP4(24,M+1) + Rpq(n,1)*WTUVTEMP5(39,M+1)
!!$        !     WTUVTEMP6(61,M) = Rpq(n,2)*WTUVTEMP5(38,M+1)
!!$        WTUVTEMP6(62,M) = 3*WTUVTEMP4(26,M+1) + Rpq(n,1)*WTUVTEMP5(41,M+1)
!!$        WTUVTEMP6(63,M) = 2*WTUVTEMP4(27,M+1) + Rpq(n,1)*WTUVTEMP5(42,M+1)
!!$        !     WTUVTEMP6(64,M) = Rpq(n,3)*WTUVTEMP5(39,M+1)
!!$        !     WTUVTEMP6(65,M) = Rpq(n,2)*WTUVTEMP5(41,M+1)
!!$        WTUVTEMP6(66,M) = 2*WTUVTEMP4(30,M+1) + Rpq(n,1)*WTUVTEMP5(45,M+1)
!!$        WTUVTEMP6(67,M) = WTUVTEMP4(31,M+1) + Rpq(n,1)*WTUVTEMP5(46,M+1)
!!$        !     WTUVTEMP6(68,M) = Rpq(n,3)*WTUVTEMP5(42,M+1)
!!$        WTUVTEMP6(69,M) = WTUVTEMP4(33,M+1) + Rpq(n,1)*WTUVTEMP5(48,M+1)
!!$        !     WTUVTEMP6(70,M) = Rpq(n,M+1)*WTUVTEMP5(45,M+1)
!!$        WTUVTEMP6(71,M) = WTUVTEMP4(35,M+1) + Rpq(n,1)*WTUVTEMP5(50,M+1)
!!$        WTUVTEMP6(72,M) = Rpq(n,1)*WTUVTEMP5(51,M+1)
!!$        !     WTUVTEMP6(73,M) = Rpq(n,1)*WTUVTEMP5(52,M+1)
!!$        WTUVTEMP6(74,M) = Rpq(n,1)*WTUVTEMP5(53,M+1)
!!$        WTUVTEMP6(75,M) = Rpq(n,1)*WTUVTEMP5(54,M+1)
!!$        !     WTUVTEMP6(76,M) = Rpq(n,1)*WTUVTEMP5(55,M+1)
!!$        WTUVTEMP6(77,M) = Rpq(n,1)*WTUVTEMP5(56,M+1)
!!$        WTUVTEMP6(78,M) = 5*WTUVTEMP4(31,M+1) + Rpq(n,2)*WTUVTEMP5(51,M+1)
!!$        WTUVTEMP6(79,M) = Rpq(n,3)*WTUVTEMP5(51,M+1)
!!$        WTUVTEMP6(80,M) = 3*WTUVTEMP4(33,M+1) + Rpq(n,2)*WTUVTEMP5(53,M+1)
!!$        WTUVTEMP6(81,M) = 2*WTUVTEMP4(34,M+1) + Rpq(n,2)*WTUVTEMP5(54,M+1)
!!$        WTUVTEMP6(82,M) = WTUVTEMP4(35,M+1) + Rpq(n,2)*WTUVTEMP5(55,M+1)
!!$        WTUVTEMP6(83,M) = Rpq(n,2)*WTUVTEMP5(56,M+1)
!!$        WTUVTEMP6(84,M) = 5*WTUVTEMP4(35,M+1) + Rpq(n,3)*WTUVTEMP5(56,M+1)
!!$     ENDDO
!!$
!!$     WTUV(n, 85) = 6*WTUVTEMP5(36,1) + Rpq(n,1)*WTUVTEMP6(57,1)
!!$     WTUV(n, 86) = Rpq(n,2)*WTUVTEMP6(57,1)
!!$     WTUV(n, 87) = Rpq(n,3)*WTUVTEMP6(57,1)
!!$     WTUV(n, 88) = 4*WTUVTEMP5(39,1) + Rpq(n,1)*WTUVTEMP6(60,1)
!!$     WTUV(n, 89) = Rpq(n,2)*WTUVTEMP6(59,1)
!!$     WTUV(n, 90) = 4*WTUVTEMP5(41,1) + Rpq(n,1)*WTUVTEMP6(62,1)
!!$     WTUV(n, 91) = 3*WTUVTEMP5(42,1) + Rpq(n,1)*WTUVTEMP6(63,1)
!!$     WTUV(n, 92) = Rpq(n,3)*WTUVTEMP6(60,1)
!!$     WTUV(n, 93) = Rpq(n,2)*WTUVTEMP6(62,1)
!!$     WTUV(n, 94) = 3*WTUVTEMP5(45,1) + Rpq(n,1)*WTUVTEMP6(66,1)
!!$     WTUV(n, 95) = 2*WTUVTEMP5(46,1) + Rpq(n,1)*WTUVTEMP6(67,1)
!!$     WTUV(n, 96) = Rpq(n,3)*WTUVTEMP6(63,1)
!!$     WTUV(n, 97) = 2*WTUVTEMP5(48,1) + Rpq(n,1)*WTUVTEMP6(69,1)
!!$     WTUV(n, 98) = Rpq(n,2)*WTUVTEMP6(66,1)
!!$     WTUV(n, 99) = 2*WTUVTEMP5(50,1) + Rpq(n,1)*WTUVTEMP6(71,1)
!!$     WTUV(n,100) = WTUVTEMP5(51,1) + Rpq(n,1)*WTUVTEMP6(72,1)
!!$     WTUV(n,101) = Rpq(n,3)*WTUVTEMP6(67,1)
!!$     WTUV(n,102) = WTUVTEMP5(53,1) + Rpq(n,1)*WTUVTEMP6(74,1)
!!$     WTUV(n,103) = WTUVTEMP5(54,1) + Rpq(n,1)*WTUVTEMP6(75,1)
!!$     WTUV(n,104) = Rpq(n,2)*WTUVTEMP6(71,1)
!!$     WTUV(n,105) = WTUVTEMP5(56,1) + Rpq(n,1)*WTUVTEMP6(77,1)
!!$     WTUV(n,106) = Rpq(n,1)*WTUVTEMP6(78,1)
!!$     WTUV(n,107) = Rpq(n,1)*WTUVTEMP6(79,1)
!!$     WTUV(n,108) = Rpq(n,1)*WTUVTEMP6(80,1)
!!$     WTUV(n,109) = Rpq(n,1)*WTUVTEMP6(81,1)
!!$     WTUV(n,110) = Rpq(n,1)*WTUVTEMP6(82,1)
!!$     WTUV(n,111) = Rpq(n,1)*WTUVTEMP6(83,1)
!!$     WTUV(n,112) = Rpq(n,1)*WTUVTEMP6(84,1)
!!$     WTUV(n,113) = 6*WTUVTEMP5(51,1) + Rpq(n,2)*WTUVTEMP6(78,1)
!!$     WTUV(n,114) = Rpq(n,3)*WTUVTEMP6(78,1)
!!$     WTUV(n,115) = 4*WTUVTEMP5(53,1) + Rpq(n,2)*WTUVTEMP6(80,1)
!!$     WTUV(n,116) = 3*WTUVTEMP5(54,1) + Rpq(n,2)*WTUVTEMP6(81,1)
!!$     WTUV(n,117) = 2*WTUVTEMP5(55,1) + Rpq(n,2)*WTUVTEMP6(82,1)
!!$     WTUV(n,118) = WTUVTEMP5(56,1) + Rpq(n,2)*WTUVTEMP6(83,1)
!!$     WTUV(n,119) = Rpq(n,2)*WTUVTEMP6(84,1)
!!$     WTUV(n,120) = 6*WTUVTEMP5(56,1) + Rpq(n,3)*WTUVTEMP6(84,1)
!!$
!!$     WTUVTEMP7( 85,1) = 6*WTUVTEMP5(36,2) + Rpq(n,1)*WTUVTEMP6(57,2)
!!$!     WTUVTEMP7( 86,1) = Rpq(n,2)*WTUVTEMP6(57,2)
!!$     WTUVTEMP7( 87,1) = Rpq(n,3)*WTUVTEMP6(57,2)
!!$     WTUVTEMP7( 88,1) = 4*WTUVTEMP5(39,2) + Rpq(n,1)*WTUVTEMP6(60,2)
!!$!     WTUVTEMP7( 89,1) = Rpq(n,2)*WTUVTEMP6(59,2)
!!$     WTUVTEMP7( 90,1) = 4*WTUVTEMP5(41,2) + Rpq(n,1)*WTUVTEMP6(62,2)
!!$     WTUVTEMP7( 91,1) = 3*WTUVTEMP5(42,2) + Rpq(n,1)*WTUVTEMP6(63,2)
!!$!     WTUVTEMP7( 92,1) = Rpq(n,3)*WTUVTEMP6(60,2)
!!$!     WTUVTEMP7( 93,1) = Rpq(n,2)*WTUVTEMP6(62,2)
!!$     WTUVTEMP7( 94,1) = 3*WTUVTEMP5(45,2) + Rpq(n,1)*WTUVTEMP6(66,2)
!!$     WTUVTEMP7( 95,1) = 2*WTUVTEMP5(46,2) + Rpq(n,1)*WTUVTEMP6(67,2)
!!$!     WTUVTEMP7( 96,1) = Rpq(n,3)*WTUVTEMP6(63,2)
!!$     WTUVTEMP7( 97,1) = 2*WTUVTEMP5(48,2) + Rpq(n,1)*WTUVTEMP6(69,2)
!!$!     WTUVTEMP7( 98,1) = Rpq(n,2)*WTUVTEMP6(66,2)
!!$     WTUVTEMP7( 99,1) = 2*WTUVTEMP5(50,2) + Rpq(n,1)*WTUVTEMP6(71,2)
!!$     WTUVTEMP7(100,1) = WTUVTEMP5(51,2) + Rpq(n,1)*WTUVTEMP6(72,2)
!!$!     WTUVTEMP7(101,1) = Rpq(n,3)*WTUVTEMP6(67,2)
!!$     WTUVTEMP7(102,1) = WTUVTEMP5(53,2) + Rpq(n,1)*WTUVTEMP6(74,2)
!!$     WTUVTEMP7(103,1) = WTUVTEMP5(54,2) + Rpq(n,1)*WTUVTEMP6(75,2)
!!$!     WTUVTEMP7(104,1) = Rpq(n,2)*WTUVTEMP6(71,2)
!!$     WTUVTEMP7(105,1) = WTUVTEMP5(56,2) + Rpq(n,1)*WTUVTEMP6(77,2)
!!$     WTUVTEMP7(106,1) = Rpq(n,1)*WTUVTEMP6(78,2)
!!$!     WTUVTEMP7(107,1) = Rpq(n,1)*WTUVTEMP6(79,2)
!!$     WTUVTEMP7(108,1) = Rpq(n,1)*WTUVTEMP6(80,2)
!!$     WTUVTEMP7(109,1) = Rpq(n,1)*WTUVTEMP6(81,2)
!!$     WTUVTEMP7(110,1) = Rpq(n,1)*WTUVTEMP6(82,2)
!!$!     WTUVTEMP7(111,1) = Rpq(n,1)*WTUVTEMP6(83,2)
!!$     WTUVTEMP7(112,1) = Rpq(n,1)*WTUVTEMP6(84,2)
!!$     WTUVTEMP7(113,1) = 6*WTUVTEMP5(51,2) + Rpq(n,2)*WTUVTEMP6(78,2)
!!$     WTUVTEMP7(114,1) = Rpq(n,3)*WTUVTEMP6(78,2)
!!$     WTUVTEMP7(115,1) = 4*WTUVTEMP5(53,2) + Rpq(n,2)*WTUVTEMP6(80,2)
!!$     WTUVTEMP7(116,1) = 3*WTUVTEMP5(54,2) + Rpq(n,2)*WTUVTEMP6(81,2)
!!$     WTUVTEMP7(117,1) = 2*WTUVTEMP5(55,2) + Rpq(n,2)*WTUVTEMP6(82,2)
!!$     WTUVTEMP7(118,1) = WTUVTEMP5(56,2) + Rpq(n,2)*WTUVTEMP6(83,2)
!!$     WTUVTEMP7(119,1) = Rpq(n,2)*WTUVTEMP6(84,2)
!!$     WTUVTEMP7(120,1) = 6*WTUVTEMP5(56,2) + Rpq(n,3)*WTUVTEMP6(84,2)
!!$
!!$     WTUV(n, 121) = 7*WTUVTEMP6(  57,1) + Rpq(n,1)*WTUVTEMP7(  85,1)
!!$     WTUV(n, 122) = Rpq(n,2)*WTUVTEMP7(  85,1)
!!$     WTUV(n, 123) = Rpq(n,3)*WTUVTEMP7(  85,1)
!!$     WTUV(n, 124) = 5*WTUVTEMP6(  60,1) + Rpq(n,1)*WTUVTEMP7(  88,1)
!!$     WTUV(n, 125) = Rpq(n,2)*WTUVTEMP7(  87,1)
!!$     WTUV(n, 126) = 5*WTUVTEMP6(  62,1) + Rpq(n,1)*WTUVTEMP7(  90,1)
!!$     WTUV(n, 127) = 4*WTUVTEMP6(  63,1) + Rpq(n,1)*WTUVTEMP7(  91,1)
!!$     WTUV(n, 128) = Rpq(n,3)*WTUVTEMP7(  88,1)
!!$     WTUV(n, 129) = Rpq(n,2)*WTUVTEMP7(  90,1)
!!$     WTUV(n, 130) = 4*WTUVTEMP6(  66,1) + Rpq(n,1)*WTUVTEMP7(  94,1)
!!$     WTUV(n, 131) = 3*WTUVTEMP6(  67,1) + Rpq(n,1)*WTUVTEMP7(  95,1)
!!$     WTUV(n, 132) = Rpq(n,3)*WTUVTEMP7(  91,1)
!!$     WTUV(n, 133) = 3*WTUVTEMP6(  69,1) + Rpq(n,1)*WTUVTEMP7(  97,1)
!!$     WTUV(n, 134) = Rpq(n,2)*WTUVTEMP7(  94,1)
!!$     WTUV(n, 135) = 3*WTUVTEMP6(  71,1) + Rpq(n,1)*WTUVTEMP7(  99,1)
!!$     WTUV(n, 136) = 2*WTUVTEMP6(  72,1) + Rpq(n,1)*WTUVTEMP7( 100,1)
!!$     WTUV(n, 137) = Rpq(n,3)*WTUVTEMP7(  95,1)
!!$     WTUV(n, 138) = 2*WTUVTEMP6(  74,1) + Rpq(n,1)*WTUVTEMP7( 102,1)
!!$     WTUV(n, 139) = 2*WTUVTEMP6(  75,1) + Rpq(n,1)*WTUVTEMP7( 103,1)
!!$     WTUV(n, 140) = Rpq(n,2)*WTUVTEMP7(  99,1)
!!$     WTUV(n, 141) = 2*WTUVTEMP6(  77,1) + Rpq(n,1)*WTUVTEMP7( 105,1)
!!$     WTUV(n, 142) = WTUVTEMP6(  78,1) + Rpq(n,1)*WTUVTEMP7( 106,1)
!!$     WTUV(n, 143) = Rpq(n,3)*WTUVTEMP7( 100,1)
!!$     WTUV(n, 144) = WTUVTEMP6(  80,1) + Rpq(n,1)*WTUVTEMP7( 108,1)
!!$     WTUV(n, 145) = WTUVTEMP6(  81,1) + Rpq(n,1)*WTUVTEMP7( 109,1)
!!$     WTUV(n, 146) = WTUVTEMP6(  82,1) + Rpq(n,1)*WTUVTEMP7( 110,1)
!!$     WTUV(n, 147) = Rpq(n,2)*WTUVTEMP7( 105,1)
!!$     WTUV(n, 148) = WTUVTEMP6(  84,1) + Rpq(n,1)*WTUVTEMP7( 112,1)
!!$     WTUV(n, 149) = Rpq(n,1)*WTUVTEMP7( 113,1)
!!$     WTUV(n, 150) = Rpq(n,1)*WTUVTEMP7( 114,1)
!!$     WTUV(n, 151) = Rpq(n,1)*WTUVTEMP7( 115,1)
!!$     WTUV(n, 152) = Rpq(n,1)*WTUVTEMP7( 116,1)
!!$     WTUV(n, 153) = Rpq(n,1)*WTUVTEMP7( 117,1)
!!$     WTUV(n, 154) = Rpq(n,1)*WTUVTEMP7( 118,1)
!!$     WTUV(n, 155) = Rpq(n,1)*WTUVTEMP7( 119,1)
!!$     WTUV(n, 156) = Rpq(n,1)*WTUVTEMP7( 120,1)
!!$     WTUV(n, 157) = 7*WTUVTEMP6(  78,1) + Rpq(n,2)*WTUVTEMP7( 113,1)
!!$     WTUV(n, 158) = Rpq(n,3)*WTUVTEMP7( 113,1)
!!$     WTUV(n, 159) = 5*WTUVTEMP6(  80,1) + Rpq(n,2)*WTUVTEMP7( 115,1)
!!$     WTUV(n, 160) = 4*WTUVTEMP6(  81,1) + Rpq(n,2)*WTUVTEMP7( 116,1)
!!$     WTUV(n, 161) = 3*WTUVTEMP6(  82,1) + Rpq(n,2)*WTUVTEMP7( 117,1)
!!$     WTUV(n, 162) = 2*WTUVTEMP6(  83,1) + Rpq(n,2)*WTUVTEMP7( 118,1)
!!$     WTUV(n, 163) = WTUVTEMP6(  84,1) + Rpq(n,2)*WTUVTEMP7( 119,1)
!!$     WTUV(n, 164) = Rpq(n,2)*WTUVTEMP7( 120,1)
!!$     WTUV(n, 165) = 7*WTUVTEMP6(  84,1) + Rpq(n,3)*WTUVTEMP7( 120,1)
!!$  ENDDO
!!$END SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX8

subroutine PrintWTUV(WTUV,AngmomPQ,nPrimQP,nPasses,nTUV,lupri)
  implicit none
  integer,intent(in) :: AngmomPQ,nPrimQP,nPasses,nTUV,lupri
  real(realk),intent(in) :: WTUV(nPrimQP,nPasses,nTUV) 
  !
  integer :: TUV,J,T,U,V,I,iPass
  WRITE(LUPRI,*)'Output from WTUVrecurrence'
  WRITE (LUPRI,'(2X,A,I10)') 'JMAX   :', AngmomPQ
  WRITE (LUPRI,'(2X,A,I10)') 'NPrimQP:', nPrimQP
  WRITE (LUPRI,'(2X,A,I10)') 'nPasses:', nPasses
  WRITE (LUPRI,'(2X,A,I10)') 'NTUV   :', nTUV
  WRITE(LUPRI,*)'Hermite integrals S(t,u,v)'
  DO iPass = 1,nPasses
     WRITE (LUPRI,*) ' '
     WRITE (LUPRI,'(2X,A,I5)') 'iPass = ',iPass
     TUV=0
     DO J = 0, AngmomPQ
        DO T = J,0,-1
           DO U = J-T,0,-1
              V=J-T-U
              TUV=TUV+1
              WRITE (LUPRI,'(2X,A2,I3,A1,I3,A1,I3,A1,2X,5ES16.8/,(18X,5ES16.8))')&
                   & 'W(',T,',',U,',',V,')', (WTUV(I,iPass,TUV),I=1,nPrimQP)
              WRITE (LUPRI,*) ' '
           ENDDO
        ENDDO
     ENDDO
  ENDDO
end subroutine PrintWTUV

end MODULE IchorEriCoulombintegralCPUMcMGeneralWTUVMod
