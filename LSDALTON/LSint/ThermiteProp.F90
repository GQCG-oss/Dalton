!> @file
!> Contains the orbital and overlap types, and subroutines to set these up
MODULE Thermite_prop
use typedef
use thermite_OD
use thermite_integrals
use OD_type
use ls_util
use OverlapType
CONTAINS
SUBROUTINE OVERLAPINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               INTEGRAL(iPrimP,ijk) = SX0*SY0*SZ0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE OVERLAPINTGEN

SUBROUTINE OVERLAPINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1
real(realk),parameter :: D0=0.0E0_realk
DO ijk=1,ijkP
   INTEGRAL(ijk)=D0
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1 = AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               INTEGRAL(ijk)=INTEGRAL(ijk) + SX0*SY0*SZ0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE OVERLAPINTSEG

SUBROUTINE DIPLENINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,POdist,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),POdist(nPrimP,3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DX0 = POdist(iPrimP,1)*SX0
               DY0 = POdist(iPrimP,2)*SY0
               DZ0 = POdist(iPrimP,3)*SZ0
               IF(iP1+iP2.GT.0) DX0 = DX0 + SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,1,0,0,1)
               IF(jP1+jP2.GT.0) DY0 = DY0 + SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,1,0,0,2)
               IF(kP1+kP2.GT.0) DZ0 = DZ0 + SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,1,0,0,3)
               INTEGRAL(iPrimP,1,ijk) = DX0*SY0*SZ0
               INTEGRAL(iPrimP,2,ijk) = SX0*DY0*SZ0
               INTEGRAL(iPrimP,3,ijk) = SX0*SY0*DZ0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE DIPLENINTGEN

SUBROUTINE DIPLENINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,POdist,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),POdist(nPrimP,3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1
real(realk),parameter :: D0=0.0E0_realk
DO ijk=1,ijkP
   INTEGRAL(1,ijk)=D0
   INTEGRAL(2,ijk)=D0
   INTEGRAL(3,ijk)=D0
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1 = AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DX0 = POdist(iPrimP,1)*SX0
               DY0 = POdist(iPrimP,2)*SY0
               DZ0 = POdist(iPrimP,3)*SZ0
               IF(iP1+iP2.GT.0) DX0 = DX0 + SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,1,0,0,1)
               IF(jP1+jP2.GT.0) DY0 = DY0 + SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,1,0,0,2)
               IF(kP1+kP2.GT.0) DZ0 = DZ0 + SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,1,0,0,3)
               INTEGRAL(1,ijk)=INTEGRAL(1,ijk)+ DX0*SY0*SZ0
               INTEGRAL(2,ijk)=INTEGRAL(2,ijk)+ SX0*DY0*SZ0
               INTEGRAL(3,ijk)=INTEGRAL(3,ijk)+ SX0*SY0*DZ0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE DIPLENINTSEG

SUBROUTINE DIPVELINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,1,0,1)
               DY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,1,0,2)
               DZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,1,0,3)
               INTEGRAL(iPrimP,1,ijk) = DX0*SY0*SZ0
               INTEGRAL(iPrimP,2,ijk) = SX0*DY0*SZ0
               INTEGRAL(iPrimP,3,ijk) = SX0*SY0*DZ0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE DIPVELINTGEN

SUBROUTINE DIPVELINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1
real(realk),parameter :: D0=0.0E0_realk
DO ijk=1,ijkP
   INTEGRAL(1,ijk)=D0
   INTEGRAL(2,ijk)=D0
   INTEGRAL(3,ijk)=D0
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1 = AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,1,0,1)
               DY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,1,0,2)
               DZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,1,0,3)
               INTEGRAL(1,ijk)=INTEGRAL(1,ijk)+ DX0*SY0*SZ0
               INTEGRAL(2,ijk)=INTEGRAL(2,ijk)+ SX0*DY0*SZ0
               INTEGRAL(3,ijk)=INTEGRAL(3,ijk)+ SX0*SY0*DZ0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE DIPVELINTSEG

SUBROUTINE THETAINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX1,DY1,DZ1,DX2,DY2,DZ2
real(realk) :: X2,Y2,Z2,R2
real(realk),parameter :: D0=0.0E0_realk,D2=2.0E0_realk,D1P5=1.5E0_realk,D4=4.0E0_realk
ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)
               DY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
               DZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
               DX2 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,2,1)
               DY2 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,2,2)
               DZ2 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,2,3)
               X2 = DX2*SY0*SZ0
               Y2 = SX0*DY2*SZ0
               Z2 = SX0*SY0*DZ2
               R2 = (X2 + Y2 + Z2)/D2               
               INTEGRAL(iPrimP,1,ijk) = D1P5*X2 - R2
               INTEGRAL(iPrimP,2,ijk) = D1P5*DX1*DY1*SZ0
               INTEGRAL(iPrimP,3,ijk) = D1P5*DX1*SY0*DZ1
               INTEGRAL(iPrimP,4,ijk) = D1P5*Y2 - R2
               INTEGRAL(iPrimP,5,ijk) = D1P5*SX0*DY1*DZ1
               INTEGRAL(iPrimP,6,ijk) = D1P5*Z2 - R2
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE THETAINTGEN

SUBROUTINE THETAINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX1,DY1,DZ1,DX2,DY2,DZ2
real(realk) :: X2,Y2,Z2,R2
real(realk),parameter :: D0=0.0E0_realk,D2=2.0E0_realk,D1P5=1.5E0_realk,D4=4.0E0_realk

DO ijk=1,ijkP
   INTEGRAL(1,ijk)=D0
   INTEGRAL(2,ijk)=D0
   INTEGRAL(3,ijk)=D0
   INTEGRAL(4,ijk)=D0
   INTEGRAL(5,ijk)=D0
   INTEGRAL(6,ijk)=D0
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1 = AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)           
               DY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
               DZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
               DX2 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,2,1)           
               DY2 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,2,2)
               DZ2 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,2,3)
               X2 = DX2*SY0*SZ0
               Y2 = SX0*DY2*SZ0
               Z2 = SX0*SY0*DZ2
               R2 = (X2 + Y2 + Z2)/D2               
               INTEGRAL(1,ijk)=INTEGRAL(1,ijk) + D1P5*X2 - R2
               INTEGRAL(2,ijk)=INTEGRAL(2,ijk) + D1P5*DX1*DY1*SZ0
               INTEGRAL(3,ijk)=INTEGRAL(3,ijk) + D1P5*DX1*SY0*DZ1
               INTEGRAL(4,ijk)=INTEGRAL(4,ijk) + D1P5*Y2 - R2
               INTEGRAL(5,ijk)=INTEGRAL(5,ijk) + D1P5*SX0*DY1*DZ1
               INTEGRAL(6,ijk)=INTEGRAL(6,ijk) + D1P5*Z2 - R2
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE THETAINTSEG

SUBROUTINE ROTSTRINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1,OVERLAP,DX1,DY1,DZ1
real(realk),parameter :: D0=0.0E0_realk,D2=2.0E0_realk

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1 = AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,1,0,1)
               DY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,1,0,2)
               DZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,1,0,3)
               SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)
               SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
               SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
               DX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,1,1,1)
               DY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,1,1,2)
               DZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,1,1,3)
               OVERLAP = SX0*SY0*SZ0
               INTEGRAL(iPrimP,1,ijk)=D2*DX1*SY0*SZ0 - OVERLAP
               INTEGRAL(iPrimP,2,ijk)=(SX1*DY0 + DX0*SY1)*SZ0
               INTEGRAL(iPrimP,3,ijk)=(SX1*DZ0 + DX0*SZ1)*SY0
               INTEGRAL(iPrimP,4,ijk)=D2*SX0*DY1*SZ0 - OVERLAP
               INTEGRAL(iPrimP,5,ijk)=SX0*(DY0*SZ1 + SY1*DZ0)
               INTEGRAL(iPrimP,6,ijk)=D2*SX0*SY0*DZ1 - OVERLAP
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE ROTSTRINTGEN

SUBROUTINE ROTSTRINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1,OVERLAP,DX1,DY1,DZ1
real(realk),parameter :: D0=0.0E0_realk,D2=2.0E0_realk
DO ijk=1,ijkP
   INTEGRAL(1,ijk)=D0
   INTEGRAL(2,ijk)=D0
   INTEGRAL(3,ijk)=D0
   INTEGRAL(4,ijk)=D0
   INTEGRAL(5,ijk)=D0
   INTEGRAL(6,ijk)=D0
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1 = AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,1,0,1)
               DY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,1,0,2)
               DZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,1,0,3)
               SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)
               SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
               SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
               DX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,1,1,1)
               DY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,1,1,2)
               DZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,1,1,3)
               OVERLAP = SX0*SY0*SZ0
               INTEGRAL(1,ijk)=INTEGRAL(1,ijk) + D2*DX1*SY0*SZ0 - OVERLAP
               INTEGRAL(2,ijk)=INTEGRAL(2,ijk) + (SX1*DY0 + DX0*SY1)*SZ0
               INTEGRAL(3,ijk)=INTEGRAL(3,ijk) + (SX1*DZ0 + DX0*SZ1)*SY0
               INTEGRAL(4,ijk)=INTEGRAL(4,ijk) + D2*SX0*DY1*SZ0 - OVERLAP
               INTEGRAL(5,ijk)=INTEGRAL(5,ijk) + SX0*(DY0*SZ1 + SY1*DZ0)
               INTEGRAL(6,ijk)=INTEGRAL(6,ijk) + D2*SX0*SY0*DZ1 - OVERLAP
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE ROTSTRINTSEG

SUBROUTINE ANGMOMINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1
ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)           
               DX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,1,0,1)
               SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               DY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,1,0,2)
               SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,1,0,3)
               SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
               INTEGRAL(iPrimP,1,ijk) = -(SY1*DZ0 - DY0*SZ1)*SX0
               INTEGRAL(iPrimP,2,ijk) = -(SZ1*DX0 - DZ0*SX1)*SY0
               INTEGRAL(iPrimP,3,ijk) = -(SX1*DY0 - DX0*SY1)*SZ0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE ANGMOMINTGEN

SUBROUTINE ANGMOMINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,DX0,DY0,DZ0,SX1,SY1,SZ1
real(realk),parameter :: D0=0.0E0_realk
DO ijk=1,ijkP
   INTEGRAL(1,ijk)=D0
   INTEGRAL(2,ijk)=D0
   INTEGRAL(3,ijk)=D0
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1 = AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1) 
               DX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,1,0,1)
               SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)               
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               DY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,1,0,2)
               SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               DZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,1,0,3)
               SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
               INTEGRAL(1,ijk)=INTEGRAL(1,ijk)-(SY1*DZ0-DY0*SZ1)*SX0
               INTEGRAL(2,ijk)=INTEGRAL(2,ijk)-(SZ1*DX0-DZ0*SX1)*SY0
               INTEGRAL(3,ijk)=INTEGRAL(3,ijk)-(SX1*DY0-DX0*SY1)*SZ0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE ANGMOMINTSEG

SUBROUTINE LONMOM1INTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,CDdist,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),CDdist(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,SX1,SY1,SZ1,DX2,DY2,DZ2,SX21,SY21,SZ21
real(realk),parameter :: D0=0.0E0_realk
real(realk),parameter :: D4INV = 0.250E0_realk, D2INV = 0.50E0_realk

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1 = AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    DO iPrimP=1,nPrimP
     SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1) 
     SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
     SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
     
     SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)               
     SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
     SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
     
     SX21 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,2,1,1)               
     SY21 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,2,1,2)
     SZ21 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,2,1,3)
     
     DX2 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,2,0,1)
     DY2 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,2,0,2)
     DZ2 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,2,0,3)

     INTEGRAL(iPrimP,1,ijk)=&
          & - D4INV *(CDdist(2)*(DX2*SY0*SZ1 + SX0*DY2*SZ1 + SX0*SY0*SZ21) &
          & - CDdist(3)*(DX2*SY1*SZ0 + SX0*SY21*SZ0 + SX0*SY1*DZ2))
     INTEGRAL(iPRimP,2,ijk)=&
          & - D4INV *(CDdist(3)*(SX21*SY0*SZ0 + SX1*DY2*SZ0 + SX1*SY0*DZ2) &
          & - CDdist(1)*(DX2*SY0*SZ1 + SX0*DY2*SZ1 + SX0*SY0*SZ21))
     INTEGRAL(iPrimP,3,ijk)=&
          & - D4INV *(CDdist(1)*(DX2*SY1*SZ0 + SX0*SY21*SZ0 + SX0*SY1*DZ2) &
          & - CDdist(2)*(SX21*SY0*SZ0 + SX1*DY2*SZ0 + SX1*SY0*DZ2))
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE LONMOM1INTGEN

SUBROUTINE LONMOM1INTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,CDdist,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,JMAX
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),CDdist(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,SX1,SY1,SZ1,DX2,DY2,DZ2,SX21,SY21,SZ21
real(realk),parameter :: D0=0.0E0_realk
real(realk),parameter :: D4INV = 0.250E0_realk, D2INV = 0.50E0_realk
!
integer :: MAXT,MAXU,MAXV,IU,IV,IT,IUMAX,ITMAX,ituv
real(realk) :: EV,FV,EU,FU,ET,FT,FX,FY,FZ,ATUV

DO ijk=1,ijkP
   INTEGRAL(1,ijk)=D0
   INTEGRAL(2,ijk)=D0
   INTEGRAL(3,ijk)=D0
ENDDO
!       Kinetic energy contribution
ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1 = AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    DO iPrimP=1,nPrimP
     SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1) 
     SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
     SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
     
     SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)               
     SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
     SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
     
     SX21 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,2,1,1)               
     SY21 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,2,1,2)
     SZ21 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,2,1,3)
     
     DX2 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,2,0,1)
     DY2 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,2,0,2)
     DZ2 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,2,0,3)

     INTEGRAL(1,ijk)=INTEGRAL(1,ijk) &
          & - D4INV *(CDdist(2)*(DX2*SY0*SZ1 + SX0*DY2*SZ1 + SX0*SY0*SZ21) &
          & - CDdist(3)*(DX2*SY1*SZ0 + SX0*SY21*SZ0 + SX0*SY1*DZ2))
     INTEGRAL(2,ijk)=INTEGRAL(2,ijk) &
          & - D4INV *(CDdist(3)*(SX21*SY0*SZ0 + SX1*DY2*SZ0 + SX1*SY0*DZ2) &
          & - CDdist(1)*(DX2*SY0*SZ1 + SX0*DY2*SZ1 + SX0*SY0*SZ21))
     INTEGRAL(3,ijk)=INTEGRAL(3,ijk) &
          & - D4INV *(CDdist(1)*(DX2*SY1*SZ0 + SX0*SY21*SZ0 + SX0*SY1*DZ2) &
          & - CDdist(2)*(SX21*SY0*SZ0 + SX1*DY2*SZ0 + SX1*SY0*DZ2))

    ENDDO
   ENDDO 
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE LONMOM1INTSEG

SUBROUTINE LONMOM2INTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,CDdist,RTUV,nTUV,SharedTUV,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,JMAX
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),CDdist(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,SX1,SY1,SZ1,DX2,DY2,DZ2,SX21,SY21,SZ21
real(realk),parameter :: D0=0.0E0_realk
real(realk),parameter :: D4INV = 0.250E0_realk, D2INV = 0.50E0_realk
!
integer :: MAXT,MAXU,MAXV,IU,IV,IT,IUMAX,ITMAX,ituv,X
real(realk) :: EV,FV,EU,FU,ET,FT,FX,FY,FZ,ATUV,signijk
real(realk),parameter :: D1=1.0E0_realk

DO ijk=1,ijkP
  DO X=1,3
     DO iPrimP=1,nPrimP
        INTEGRAL(iPrimP,X,ijk)=D0
     ENDDO
  ENDDO
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1 = AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    MAXT=iP1+iP2
    MAXU=jP1+jP2
    MAXV=kP1+kP2
    DO iPrimP=1,nPrimP
     DO IV = 0, MAXV + 1
      EV = ODC(iPrimP,kP1,kP2,IV,0,0,3) 
      FV = ODC(iPrimP,kP1,kP2,IV,0,1,3) 
      IUMAX = MAXU + 1
      IF(IV .GT. MAXV) IUMAX = MAXU
      DO IU = 0, IUMAX
       EU = ODC(iPrimP,jP1,jP2,IU,0,0,2) 
       FU = ODC(iPrimP,jP1,jP2,IU,0,1,2) 
       ITMAX = MAXT + 1
       IF((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
       DO IT = 0, ITMAX
        ET = ODC(iPrimP,iP1,iP2,IT,0,0,1) 
        FT = ODC(iPrimP,iP1,iP2,IT,0,1,1) 
        FX = FT*EU*EV
        FY = ET*FU*EV
        FZ = ET*EU*FV
        iTUV = sharedTUV%TUVindex(it,iu,iv)
        ATUV = -D2INV*RTUV(iPrimP,iTUV)
        INTEGRAL(iPrimP,1,ijk)=INTEGRAL(iPrimP,1,ijk) + ATUV*(CDdist(2)*FZ-CDdist(3)*FY) 
        INTEGRAL(iPrimP,2,ijk)=INTEGRAL(iPrimP,2,ijk) + ATUV*(CDdist(3)*FX-CDdist(1)*FZ)
        INTEGRAL(iPrimP,3,ijk)=INTEGRAL(iPrimP,3,ijk) + ATUV*(CDdist(1)*FY-CDdist(2)*FX)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO 
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE LONMOM2INTGEN

SUBROUTINE LONMOM2INTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,CDdist,RTUV,nTUV,SharedTUV,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,JMAX
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),CDdist(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nTUV)
!real(realk) :: RTUV(0:JMAX,nPrimP)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,SX1,SY1,SZ1,DX2,DY2,DZ2,SX21,SY21,SZ21
real(realk),parameter :: D0=0.0E0_realk
real(realk),parameter :: D4INV = 0.250E0_realk, D2INV = 0.50E0_realk
!
integer :: MAXT,MAXU,MAXV,IU,IV,IT,IUMAX,ITMAX,ituv
real(realk) :: EV,FV,EU,FU,ET,FT,FX,FY,FZ,ATUV,signijk
real(realk),parameter :: D1=1.0E0_realk

DO ijk=1,ijkP
   INTEGRAL(1,ijk)=D0
   INTEGRAL(2,ijk)=D0
   INTEGRAL(3,ijk)=D0
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1 = AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    MAXT=iP1+iP2
    MAXU=jP1+jP2
    MAXV=kP1+kP2
    DO iPrimP=1,nPrimP
     DO IV = 0, MAXV + 1
      EV = ODC(iPrimP,kP1,kP2,IV,0,0,3) 
      FV = ODC(iPrimP,kP1,kP2,IV,0,1,3) 
      IUMAX = MAXU + 1
      IF(IV .GT. MAXV) IUMAX = MAXU
      DO IU = 0, IUMAX
       EU = ODC(iPrimP,jP1,jP2,IU,0,0,2) 
       FU = ODC(iPrimP,jP1,jP2,IU,0,1,2) 
       ITMAX = MAXT + 1
       IF((IU .GT. MAXU) .OR. (IV .GT. MAXV)) ITMAX = MAXT
       DO IT = 0, ITMAX
        ET = ODC(iPrimP,iP1,iP2,IT,0,0,1) 
        FT = ODC(iPrimP,iP1,iP2,IT,0,1,1) 
        FX = FT*EU*EV
        FY = ET*FU*EV
        FZ = ET*EU*FV
        iTUV = sharedTUV%TUVindex(it,iu,iv)
        ATUV = -D2INV*RTUV(iPrimP,iTUV)
        INTEGRAL(1,ijk)=INTEGRAL(1,ijk) + ATUV*(CDdist(2)*FZ-CDdist(3)*FY) 
        INTEGRAL(2,ijk)=INTEGRAL(2,ijk) + ATUV*(CDdist(3)*FX-CDdist(1)*FZ)
        INTEGRAL(3,ijk)=INTEGRAL(3,ijk) + ATUV*(CDdist(1)*FY-CDdist(2)*FX)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO 
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE LONMOM2INTSEG

SUBROUTINE POTINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,RTUV,nTUV,SharedTUV,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: EV,EU,ET,EEE,signijk
integer :: MAXT,MAXU,MAXV,IV,IU,IT,ituv
real(realk),parameter :: D1=1.0E0_realk
real(realk),parameter :: D0=0.0E0_realk

DO ijk=1,ijkP
   DO iPrimP=1,nPrimP
      INTEGRAL(iPrimP,ijk) = D0
   ENDDO
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            MAXT = iP1 + iP2
            MAXU = jP1 + jP2
            MAXV = kP1 + kP2
            DO iPrimP=1,nPrimP
               DO IV = 0, MAXV
                  EV = ODC(iPrimP,kP1,kP2,IV,0,0,3)
                  DO IU = 0, MAXU
                     EU = ODC(iPrimP,jP1,jP2,IU,0,0,2)
                     DO IT = 0, MAXT
                        ET = ODC(iPrimP,iP1,iP2,IT,0,0,1)
                        EEE = ET*EU*EV
                        ituv = sharedTUV%TUVindex(it,iu,iv)
                        INTEGRAL(iPrimP,ijk) = INTEGRAL(iPrimP,ijk) + EEE*RTUV(iPrimP,ituv)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE POTINTGEN

SUBROUTINE POTINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,RTUV,nTUV,SharedTUV,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: EV,EU,ET,EEE,signijk
integer :: MAXT,MAXU,MAXV,IV,IU,IT,ituv
real(realk),parameter :: D0=0.0E0_realk
real(realk),parameter :: D1=1.0E0_realk

DO ijk=1,ijkP
   INTEGRAL(ijk)=D0
   INTEGRAL(ijk)=D0
   INTEGRAL(ijk)=D0
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            MAXT = iP1 + iP2
            MAXU = jP1 + jP2
            MAXV = kP1 + kP2
            DO iPrimP=1,nPrimP
               DO IV = 0, MAXV
                  EV = ODC(iPrimP,kP1,kP2,IV,0,0,3)
                  DO IU = 0, MAXU
                     EU = ODC(iPrimP,jP1,jP2,IU,0,0,2)
                     DO IT = 0, MAXT
                        ET = ODC(iPrimP,iP1,iP2,IT,0,0,1)
                        EEE = ET*EU*EV
                        ituv = sharedTUV%TUVindex(it,iu,iv)
                        INTEGRAL(ijk) = INTEGRAL(ijk) + EEE*RTUV(iPrimP,ituv)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE POTINTseg

SUBROUTINE PSOINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,RTUV,nTUV,nPrimPQ,nPassQ,SharedTUV,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA
integer     :: AngmomB,ijkP,nOperatorComp,lupri,nPrimPQ,nPassQ
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,nPassQ,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nPassQ,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: EV,EU,ET,FV,FU,EE,FE,EF,FT,FEE,EFE,EEF,AH0T,AH0U,AH0V
integer :: MAXT,MAXU,MAXV,IV,IU,IT,ituv,ituv2,ituv3,ituv4,iatomc,X
real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk

DO ijk=1,ijkP
   DO IATOMC = 1,nPassQ
      DO X=1,3
         DO iPrimP=1,nPrimP
            INTEGRAL(iPrimP,X,iAtomC,ijk)=D0
         ENDDO
      ENDDO
   ENDDO
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            MAXT = iP1 + iP2 + 1 
            MAXU = jP1 + jP2 + 1 
            MAXV = kP1 + kP2 + 1 
            DO iPrimP=1,nPrimP
               DO IV = 0, MAXV
                  EV = ODC(iPrimP,kP1,kP2,IV,0,0,3)
                  FV = ODC(iPrimP,kP1,kP2,IV,1,0,3)
                  DO IU = 0, MAXU
                     EU = ODC(iPrimP,jP1,jP2,IU,0,0,2)
                     FU = ODC(iPrimP,jP1,jP2,IU,1,0,2)
                     EE = EU*EV
                     FE = FU*EV
                     EF = EU*FV
                     DO IT = 0, MAXT
                        ET = ODC(iPrimP,iP1,iP2,IT,0,0,1)
                        FT = ODC(iPrimP,iP1,iP2,IT,1,0,1)
                        FEE = FT*EE
                        EFE = ET*FE
                        EEF = ET*EF
                        DO IATOMC = 1,nPassQ
                           ituv2 = sharedTUV%TUVindex(it+1,iu,iv)
                           ituv3 = sharedTUV%TUVindex(it,iu+1,iv)
                           ituv4 = sharedTUV%TUVindex(it,iu,iv+1)
                           AH0T = RTUV(iPrimP,iAtomc,ituv2)
                           AH0U = RTUV(iPrimP,iAtomc,ituv3)
                           AH0V = RTUV(iPrimP,iAtomc,ituv4)
                           INTEGRAL(iPrimP,1,iAtomC,ijk) = INTEGRAL(iPrimP,1,iAtomC,ijk) + EFE*AH0V-EEF*AH0U
                           INTEGRAL(iPrimP,2,iAtomC,ijk) = INTEGRAL(iPrimP,2,iAtomC,ijk) + EEF*AH0T-FEE*AH0V
                           INTEGRAL(iPrimP,3,iAtomC,ijk) = INTEGRAL(iPrimP,3,iAtomC,ijk) + FEE*AH0U-EFE*AH0T
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE PSOINTGEN

SUBROUTINE PSOINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,RTUV,nTUV,nPrimPQ,nPassQ,SharedTUV,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP
integer     :: nOperatorComp,lupri,nPrimPQ,nPassQ
real(realk) :: SHGTF(nPrimP)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,nPassQ,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nPassQ,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: EV,EU,ET,FV,FU,EE,FE,EF,FT,FEE,EFE,EEF,signijk
real(realk) :: AH0T,AH0U,AH0V
integer :: MAXT,MAXU,MAXV,IV,IU,IT,ituv,ituv2,ituv3,ituv4,iatomC
real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk

DO ijk=1,ijkP
   DO IATOMC = 1,nPassQ
      INTEGRAL(1,iAtomC,ijk)=D0
      INTEGRAL(2,iAtomC,ijk)=D0
      INTEGRAL(3,iAtomC,ijk)=D0
   ENDDO
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            MAXT = iP1 + iP2 + 1 
            MAXU = jP1 + jP2 + 1 
            MAXV = kP1 + kP2 + 1 
            DO iPrimP=1,nPrimP
               DO IV = 0, MAXV
                  EV = ODC(iPrimP,kP1,kP2,IV,0,0,3)
                  FV = ODC(iPrimP,kP1,kP2,IV,1,0,3)
                  DO IU = 0, MAXU
                     EU = ODC(iPrimP,jP1,jP2,IU,0,0,2)
                     FU = ODC(iPrimP,jP1,jP2,IU,1,0,2)
                     EE = EU*EV
                     FE = FU*EV
                     EF = EU*FV
                     DO IT = 0, MAXT
                        ET = ODC(iPrimP,iP1,iP2,IT,0,0,1)
                        FT = ODC(iPrimP,iP1,iP2,IT,1,0,1)
                        FEE = FT*EE
                        EFE = ET*FE
                        EEF = ET*EF
                        DO IATOMC = 1,nPassQ
                           ituv2 = sharedTUV%TUVindex(it+1,iu,iv)
                           ituv3 = sharedTUV%TUVindex(it,iu+1,iv)
                           ituv4 = sharedTUV%TUVindex(it,iu,iv+1)
                           AH0T = RTUV(iPrimP,iAtomc,ituv2)
                           AH0U = RTUV(iPrimP,iAtomc,ituv3)
                           AH0V = RTUV(iPrimP,iAtomc,ituv4)
                           INTEGRAL(1,iAtomC,ijk) = INTEGRAL(1,iAtomC,ijk) + EFE*AH0V-EEF*AH0U
                           INTEGRAL(2,iAtomC,ijk) = INTEGRAL(2,iAtomC,ijk) + EEF*AH0T-FEE*AH0V
                           INTEGRAL(3,iAtomC,ijk) = INTEGRAL(3,iAtomC,ijk) + FEE*AH0U-EFE*AH0T
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE PSOINTseg

SUBROUTINE NSTLONINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,RTUV,nTUV,nPrimPQ,nPassQ,SharedTUV,CDdist,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP
integer     :: nOperatorComp,lupri,nPrimPQ,nPassQ
real(realk) :: SHGTF(nPrimP),CDdist(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,nPassQ,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nPassQ,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: EV,EU,ET,FV,FU,EE,FE,EF,FT,FEE,EFE,EEF,signijk
real(realk) :: DV,GV,DU,GU,DT,GT,AH0T,AH0U,AH0V
integer :: MAXT,MAXU,MAXV,IV,IU,IT,ituv,ituv2,ituv3,ituv4,iatomC,X
real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
real(realk),parameter :: D2INV = 0.50E0_realk

DO ijk=1,ijkP
   DO IATOMC = 1,nPassQ
      DO X=1,nOperatorComp
         DO iPrimP=1,nPrimP
            INTEGRAL(iPrimP,X,iAtomC,ijk)=D0
         ENDDO
      ENDDO
   ENDDO
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1=AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    MAXT = iP1 + iP2 + 2
    MAXU = jP1 + jP2 + 2 
    MAXV = kP1 + kP2 + 2 
    DO iPrimP=1,nPrimP
     DO IV = 0, MAXV
      DV = ODC(iPrimP,kP1,kP2,IV,1,1,3)
      EV = ODC(iPrimP,kP1,kP2,IV,0,0,3)
      FV = ODC(iPrimP,kP1,kP2,IV,1,0,3)
      GV = ODC(iPrimP,kP1,kP2,IV,0,1,3)
      DO IU = 0, MAXU
       DU = ODC(iPrimP,jP1,jP2,IU,1,1,2)
       EU = ODC(iPrimP,jP1,jP2,IU,0,0,2)
       FU = ODC(iPrimP,jP1,jP2,IU,1,0,2)
       GU = ODC(iPrimP,jP1,jP2,IU,0,1,2)
       DO IT = 0, MAXT
         DT = ODC(iPrimP,iP1,iP2,IT,1,1,1)
         ET = ODC(iPrimP,iP1,iP2,IT,0,0,1)
         FT = ODC(iPrimP,iP1,iP2,IT,1,0,1)
         GT = ODC(iPrimP,iP1,iP2,IT,0,1,1)
         DO IATOMC = 1,nPassQ
!            IF(IT+IU+IV+1.LT.AngmomA+AngmomB+2)THEN
               ituv2 = sharedTUV%TUVindex(it+1,iu,iv)
               ituv3 = sharedTUV%TUVindex(it,iu+1,iv)
               ituv4 = sharedTUV%TUVindex(it,iu,iv+1)
               AH0T = RTUV(iPrimP,iAtomc,ituv2)
               AH0U = RTUV(iPrimP,iAtomc,ituv3)
               AH0V = RTUV(iPrimP,iAtomc,ituv4)
!            ELSE
 !              AH0T = D0
  !             AH0U = D0
   !            AH0V = D0
    !        ENDIF
            INTEGRAL(iPrimP,1,iAtomC,ijk) = INTEGRAL(iPrimP,1,iAtomC,ijk) - (CDdist(2)*ET*(FU*GV*AH0V - EU*DV*AH0U)&
                 & + (GU*FV*AH0U - DU*EV*AH0V)*CDdist(3)*ET)*D2INV
            INTEGRAL(iPrimP,2,iAtomC,ijk) = INTEGRAL(iPrimP,2,iAtomC,ijk) - (CDdist(3)*GT*(FU*EV*AH0V - EU*FV*AH0U)&
                 & + (EU*DV*AH0U - FU*GV*AH0V)*CDdist(1)*ET)*D2INV
            INTEGRAL(iPrimP,3,iAtomC,ijk) = INTEGRAL(iPrimP,3,iAtomC,ijk) - (CDdist(1)*ET*(DU*EV*AH0V - GU*FV*AH0U)&
                 & + (EU*FV*AH0U - FU*EV*AH0V)*CDdist(2)*GT)*D2INV
            INTEGRAL(iPrimP,4,iAtomC,ijk) = INTEGRAL(iPrimP,4,iAtomC,ijk) + (CDdist(2)*EU*(FT*GV*AH0V - ET*DV*AH0T)&
                 & + (ET*FV*AH0T - FT*EV*AH0V)*CDdist(3)*GU)*D2INV
            INTEGRAL(iPrimP,5,iAtomC,ijk) = INTEGRAL(iPrimP,5,iAtomC,ijk) + (CDdist(3)*EU*(DT*EV*AH0V - GT*FV*AH0T)&
                 & + (ET*DV*AH0T - FT*GV*AH0V)*CDdist(1)*EU)*D2INV
            INTEGRAL(iPrimP,6,iAtomC,ijk) = INTEGRAL(iPrimP,6,iAtomC,ijk) + (CDdist(1)*GU*(FT*EV*AH0V - ET*FV*AH0T)&
                 & + (GT*FV*AH0T - DT*EV*AH0V)*CDdist(2)*EU)*D2INV
            INTEGRAL(iPrimP,7,iAtomC,ijk) = INTEGRAL(iPrimP,7,iAtomC,ijk) - (CDdist(2)*GV*(FT*EU*AH0U - ET*FU*AH0T)&
                 & + (ET*DU*AH0T - FT*GU*AH0U)*CDdist(3)*EV)*D2INV
            INTEGRAL(iPrimP,8,iAtomC,ijk) = INTEGRAL(iPrimP,8,iAtomC,ijk) - (CDdist(3)*EV*(DT*EU*AH0U - GT*FU*AH0T)&
                 & + (ET*FU*AH0T - FT*EU*AH0U)*CDdist(1)*GV)*D2INV
            INTEGRAL(iPrimP,9,iAtomC,ijk) = INTEGRAL(iPrimP,9,iAtomC,ijk) - (CDdist(1)*EV*(FT*GU*AH0U - ET*DU*AH0T)&
                 & + (GT*FU*AH0T - DT*EU*AH0U)*CDdist(2)*EV)*D2INV
         ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE NSTLONINTgen

SUBROUTINE NSTLONINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,RTUV,nTUV,nPrimPQ,nPassQ,SharedTUV,CDdist,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP
integer     :: nOperatorComp,lupri,nPrimPQ,nPassQ
real(realk) :: SHGTF(nPrimP),CDdist(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,nPassQ,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nPassQ,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: EV,EU,ET,FV,FU,EE,FE,EF,FT,FEE,EFE,EEF,signijk
real(realk) :: DV,GV,DU,GU,DT,GT,AH0T,AH0U,AH0V
integer :: MAXT,MAXU,MAXV,IV,IU,IT,ituv,ituv2,ituv3,ituv4,iatomC
real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
real(realk),parameter :: D2INV = 0.50E0_realk

DO ijk=1,ijkP
   DO IATOMC = 1,nPassQ
      INTEGRAL(1,iAtomC,ijk)=D0
      INTEGRAL(2,iAtomC,ijk)=D0
      INTEGRAL(3,iAtomC,ijk)=D0
      INTEGRAL(4,iAtomC,ijk)=D0
      INTEGRAL(5,iAtomC,ijk)=D0
      INTEGRAL(6,iAtomC,ijk)=D0
      INTEGRAL(7,iAtomC,ijk)=D0
      INTEGRAL(8,iAtomC,ijk)=D0
      INTEGRAL(9,iAtomC,ijk)=D0
   ENDDO
ENDDO


ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1=AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    MAXT = iP1 + iP2 + 2
    MAXU = jP1 + jP2 + 2 
    MAXV = kP1 + kP2 + 2 
    DO iPrimP=1,nPrimP
     DO IV = 0, MAXV
      DV = ODC(iPrimP,kP1,kP2,IV,1,1,3)
      EV = ODC(iPrimP,kP1,kP2,IV,0,0,3)
      FV = ODC(iPrimP,kP1,kP2,IV,1,0,3)
      GV = ODC(iPrimP,kP1,kP2,IV,0,1,3)
      DO IU = 0, MAXU
       DU = ODC(iPrimP,jP1,jP2,IU,1,1,2)
       EU = ODC(iPrimP,jP1,jP2,IU,0,0,2)
       FU = ODC(iPrimP,jP1,jP2,IU,1,0,2)
       GU = ODC(iPrimP,jP1,jP2,IU,0,1,2)
       DO IT = 0, MAXT
         DT = ODC(iPrimP,iP1,iP2,IT,1,1,1)
         ET = ODC(iPrimP,iP1,iP2,IT,0,0,1)
         FT = ODC(iPrimP,iP1,iP2,IT,1,0,1)
         GT = ODC(iPrimP,iP1,iP2,IT,0,1,1)
         DO IATOMC = 1,nPassQ
!            IF(IT+IU+IV+1.LT.AngmomA+AngmomB+2)THEN
               ituv2 = sharedTUV%TUVindex(it+1,iu,iv)
               ituv3 = sharedTUV%TUVindex(it,iu+1,iv)
               ituv4 = sharedTUV%TUVindex(it,iu,iv+1)
               AH0T = RTUV(iPrimP,iAtomc,ituv2)
               AH0U = RTUV(iPrimP,iAtomc,ituv3)
               AH0V = RTUV(iPrimP,iAtomc,ituv4)
!            ELSE
 !              AH0T = D0
  !             AH0U = D0
   !            AH0V = D0
    !        ENDIF
            INTEGRAL(1,iAtomC,ijk) = INTEGRAL(1,iAtomC,ijk) - (CDdist(2)*ET*(FU*GV*AH0V - EU*DV*AH0U) &
                 &+ (GU*FV*AH0U - DU*EV*AH0V)*CDdist(3)*ET)*D2INV
            INTEGRAL(2,iAtomC,ijk) = INTEGRAL(2,iAtomC,ijk) - (CDdist(3)*GT*(FU*EV*AH0V - EU*FV*AH0U) &
                 &+ (EU*DV*AH0U - FU*GV*AH0V)*CDdist(1)*ET)*D2INV
            INTEGRAL(3,iAtomC,ijk) = INTEGRAL(3,iAtomC,ijk) - (CDdist(1)*ET*(DU*EV*AH0V - GU*FV*AH0U) &
                 &+ (EU*FV*AH0U - FU*EV*AH0V)*CDdist(2)*GT)*D2INV
            INTEGRAL(4,iAtomC,ijk) = INTEGRAL(4,iAtomC,ijk) + (CDdist(2)*EU*(FT*GV*AH0V - ET*DV*AH0T) &
                 &+ (ET*FV*AH0T - FT*EV*AH0V)*CDdist(3)*GU)*D2INV
            INTEGRAL(5,iAtomC,ijk) = INTEGRAL(5,iAtomC,ijk) + (CDdist(3)*EU*(DT*EV*AH0V - GT*FV*AH0T) &
                 &+ (ET*DV*AH0T - FT*GV*AH0V)*CDdist(1)*EU)*D2INV
            INTEGRAL(6,iAtomC,ijk) = INTEGRAL(6,iAtomC,ijk) + (CDdist(1)*GU*(FT*EV*AH0V - ET*FV*AH0T) &
                 &+ (GT*FV*AH0T - DT*EV*AH0V)*CDdist(2)*EU)*D2INV
            INTEGRAL(7,iAtomC,ijk) = INTEGRAL(7,iAtomC,ijk) - (CDdist(2)*GV*(FT*EU*AH0U - ET*FU*AH0T) &
                 &+ (ET*DU*AH0T - FT*GU*AH0U)*CDdist(3)*EV)*D2INV
            INTEGRAL(8,iAtomC,ijk) = INTEGRAL(8,iAtomC,ijk) - (CDdist(3)*EV*(DT*EU*AH0U - GT*FU*AH0T) &
                 &+ (ET*FU*AH0T - FT*EU*AH0U)*CDdist(1)*GV)*D2INV
            INTEGRAL(9,iAtomC,ijk) = INTEGRAL(9,iAtomC,ijk) - (CDdist(1)*EV*(FT*GU*AH0U - ET*DU*AH0T) &
                 &+ (GT*FU*AH0T - DT*EU*AH0U)*CDdist(2)*EV)*D2INV
         ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE NSTLONINTseg

SUBROUTINE NSTNOLINTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,RTUV,nTUV,nPrimPQ,nPassQ,SharedTUV,CDdist,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP
integer     :: nOperatorComp,lupri,nPrimPQ,nPassQ
real(realk) :: SHGTF(nPrimP),CDdist(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,nPassQ,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nPassQ,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: EV,EU,ET,FV,FU,EE,FE,EF,FT,FEE,EFE,EEF,signijk
real(realk) :: DV,GV,DU,GU,DT,GT,AH0T,AH0U,AH0V
integer :: MAXT,MAXU,MAXV,IV,IU,IT,ituv,ituv2,ituv3,ituv4,iatomC,X
real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
real(realk),parameter :: D2INV = 0.50E0_realk

DO ijk=1,ijkP
   DO IATOMC = 1,nPassQ
      DO X=1,nOperatorComp
         DO iPrimP=1,nPrimP
            INTEGRAL(iPrimP,X,iAtomC,ijk)=D0
         ENDDO
      ENDDO
   ENDDO
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1=AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    MAXT = iP1 + iP2 + 2
    MAXU = jP1 + jP2 + 2 
    MAXV = kP1 + kP2 + 2 
    DO iPrimP=1,nPrimP
     DO IV = 0, MAXV
      EV = ODC(iPrimP,kP1,kP2,IV,0,0,3)
      FV = ODC(iPrimP,kP1,kP2,IV,0,1,3)
      DO IU = 0, MAXU
       EU = ODC(iPrimP,jP1,jP2,IU,0,0,2)
       FU = ODC(iPrimP,jP1,jP2,IU,0,1,2)
       DO IT = 0, MAXT
         ET = ODC(iPrimP,iP1,iP2,IT,0,0,1)
         FT = ODC(iPrimP,iP1,iP2,IT,0,1,1)
         FEE = FT*EU*EV
         EFE = ET*FU*EV
         EEF = ET*EU*FV
         DO IATOMC = 1,nPassQ
            ituv2 = sharedTUV%TUVindex(it+1,iu,iv)
            ituv3 = sharedTUV%TUVindex(it,iu+1,iv)
            ituv4 = sharedTUV%TUVindex(it,iu,iv+1)
            AH0T = RTUV(iPrimP,iAtomc,ituv2)
            AH0U = RTUV(iPrimP,iAtomc,ituv3)
            AH0V = RTUV(iPrimP,iAtomc,ituv4)
            INTEGRAL(iPrimP,1,iAtomC,ijk) = INTEGRAL(iPrimP,1,iAtomC,ijk)- D2INV*(EFE*AH0U + EEF*AH0V)
            INTEGRAL(iPrimP,2,iAtomC,ijk) = INTEGRAL(iPrimP,2,iAtomC,ijk)+ D2INV*FEE*AH0U
            INTEGRAL(iPrimP,3,iAtomC,ijk) = INTEGRAL(iPrimP,3,iAtomC,ijk)+ D2INV*FEE*AH0V
            INTEGRAL(iPrimP,4,iAtomC,ijk) = INTEGRAL(iPrimP,4,iAtomC,ijk)+ D2INV*EFE*AH0T
            INTEGRAL(iPrimP,5,iAtomC,ijk) = INTEGRAL(iPrimP,5,iAtomC,ijk)- D2INV*(FEE*AH0T + EEF*AH0V)
            INTEGRAL(iPrimP,6,iAtomC,ijk) = INTEGRAL(iPrimP,6,iAtomC,ijk)+ D2INV*EFE*AH0V
            INTEGRAL(iPrimP,7,iAtomC,ijk) = INTEGRAL(iPrimP,7,iAtomC,ijk)+ D2INV*EEF*AH0T
            INTEGRAL(iPrimP,8,iAtomC,ijk) = INTEGRAL(iPrimP,8,iAtomC,ijk)+ D2INV*EEF*AH0U
            INTEGRAL(iPrimP,9,iAtomC,ijk) = INTEGRAL(iPrimP,9,iAtomC,ijk)- D2INV*(FEE*AH0T + EFE*AH0U)
         ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE NSTNOLINTgen

SUBROUTINE NSTNOLINTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,RTUV,nTUV,nPrimPQ,nPassQ,SharedTUV,CDdist,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP
integer     :: nOperatorComp,lupri,nPrimPQ,nPassQ
real(realk) :: SHGTF(nPrimP),CDdist(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,nPassQ,ijkP)
!
integer     :: nTUV
real(realk) :: RTUV(nPrimP,nPassQ,nTUV)
TYPE(TUVitem),intent(in)     :: SharedTUV
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: EV,EU,ET,FV,FU,EE,FE,EF,FT,FEE,EFE,EEF,signijk
real(realk) :: DV,GV,DU,GU,DT,GT,AH0T,AH0U,AH0V
integer :: MAXT,MAXU,MAXV,IV,IU,IT,ituv,ituv2,ituv3,ituv4,iatomC
real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
real(realk),parameter :: D2INV = 0.50E0_realk

DO ijk=1,ijkP
   DO IATOMC = 1,nPassQ
      INTEGRAL(1,iAtomC,ijk)=D0
      INTEGRAL(2,iAtomC,ijk)=D0
      INTEGRAL(3,iAtomC,ijk)=D0
      INTEGRAL(4,iAtomC,ijk)=D0
      INTEGRAL(5,iAtomC,ijk)=D0
      INTEGRAL(6,iAtomC,ijk)=D0
      INTEGRAL(7,iAtomC,ijk)=D0
      INTEGRAL(8,iAtomC,ijk)=D0
      INTEGRAL(9,iAtomC,ijk)=D0
   ENDDO
ENDDO

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1=AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    MAXT = iP1 + iP2 + 2
    MAXU = jP1 + jP2 + 2 
    MAXV = kP1 + kP2 + 2 
    DO iPrimP=1,nPrimP
     DO IV = 0, MAXV
      EV = ODC(iPrimP,kP1,kP2,IV,0,0,3)
      FV = ODC(iPrimP,kP1,kP2,IV,0,1,3)
      DO IU = 0, MAXU
       EU = ODC(iPrimP,jP1,jP2,IU,0,0,2)
       FU = ODC(iPrimP,jP1,jP2,IU,0,1,2)
       DO IT = 0, MAXT
         ET = ODC(iPrimP,iP1,iP2,IT,0,0,1)
         FT = ODC(iPrimP,iP1,iP2,IT,0,1,1)
         FEE = FT*EU*EV
         EFE = ET*FU*EV
         EEF = ET*EU*FV
!         WRITE(lupri,*)'FEE',FEE
!         WRITE(lupri,*)'EFE',EFE
!         WRITE(lupri,*)'EEF',EEF
!         WRITE(lupri,*)'IT,IU,IV',IT,IV,IU
         DO IATOMC = 1,nPassQ
!            WRITE(lupri,*)'IATOMC',IATOMC
            ituv2 = sharedTUV%TUVindex(it+1,iu,iv)
            ituv3 = sharedTUV%TUVindex(it,iu+1,iv)
            ituv4 = sharedTUV%TUVindex(it,iu,iv+1)
            AH0T = RTUV(iPrimP,iAtomc,ituv2)
            AH0U = RTUV(iPrimP,iAtomc,ituv3)
            AH0V = RTUV(iPrimP,iAtomc,ituv4)
!            WRITE(lupri,*)'AH0T',AH0T
!            WRITE(lupri,*)'AH0T',AH0U
!            WRITE(lupri,*)'AH0T',AH0V

            INTEGRAL(1,iAtomC,ijk) = INTEGRAL(1,iAtomC,ijk) - D2INV*(EFE*AH0U + EEF*AH0V)
            INTEGRAL(2,iAtomC,ijk) = INTEGRAL(2,iAtomC,ijk) + D2INV*FEE*AH0U
            INTEGRAL(3,iAtomC,ijk) = INTEGRAL(3,iAtomC,ijk) + D2INV*FEE*AH0V
            INTEGRAL(4,iAtomC,ijk) = INTEGRAL(4,iAtomC,ijk) + D2INV*EFE*AH0T
            INTEGRAL(5,iAtomC,ijk) = INTEGRAL(5,iAtomC,ijk) - D2INV*(FEE*AH0T + EEF*AH0V)
            INTEGRAL(6,iAtomC,ijk) = INTEGRAL(6,iAtomC,ijk) + D2INV*EFE*AH0V
            INTEGRAL(7,iAtomC,ijk) = INTEGRAL(7,iAtomC,ijk) + D2INV*EEF*AH0T
            INTEGRAL(8,iAtomC,ijk) = INTEGRAL(8,iAtomC,ijk) + D2INV*EEF*AH0U
            INTEGRAL(9,iAtomC,ijk) = INTEGRAL(9,iAtomC,ijk) - D2INV*(FEE*AH0T + EFE*AH0U)
!            WRITE(lupri,*)'INT1',- D2INV*(EFE*AH0U + EEF*AH0V)
!            WRITE(lupri,*)'INT2',D2INV*FEE*AH0U
!            WRITE(lupri,*)'INT3',D2INV*FEE*AH0V
!            WRITE(lupri,*)'INT4',D2INV*EFE*AH0T
!            WRITE(lupri,*)'INT5',D2INV*(FEE*AH0T + EEF*AH0V)
!            WRITE(lupri,*)'INT6',D2INV*EFE*AH0V
!            WRITE(lupri,*)'INT7',D2INV*EEF*AH0T
!            WRITE(lupri,*)'INT8',D2INV*EEF*AH0U
!            WRITE(lupri,*)'INT9',D2INV*(FEE*AH0T + EFE*AH0U)
         ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE NSTNOLINTseg

SUBROUTINE DCM1INTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,distanceAB,AODIST,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),distanceAB(3),AODIST(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,SX1,SY1,SZ1,DX1,DY1,DZ1,DX2,DY2,DZ2
real(realk) :: SX2,SY2,SZ2,DIFABX,DIFABY,DIFABZ,ADX,ADY,ADZ,DPLX,DPLY,DPLZ
real(realk),PARAMETER :: D0 = 0.0E0_realk
real(realk),PARAMETER :: D2 = 2.0E0_realk
real(realk),PARAMETER :: D2I = 1.0E0_realk/D2
DIFABX = distanceAB(1)
DIFABY = distanceAB(2)
DIFABZ = distanceAB(3)
ADX = AODIST(1)
ADY = AODIST(2)
ADZ = AODIST(3)
ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)
               SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
               SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
               SX2 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,2,1)
               SY2 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,2,2)
               SZ2 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,2,3)
               DX1 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,0,1) + &
                    & ADX * ODC(iPrimP,iP1,iP2,0,0,0,1))
               DY1 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,0,2) + &
                    & ADY * ODC(iPrimP,jP1,jP2,0,0,0,2))
               DZ1 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,0,3) + &
                    & ADZ * ODC(iPrimP,kP1,kP2,0,0,0,3))
               DX2 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,1,1) + &
                    & ADX * ODC(iPrimP,iP1,iP2,0,0,1,1))
               DY2 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,1,2) + &
                    & ADY * ODC(iPrimP,jP1,jP2,0,0,1,2))
               DZ2 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,1,3) + &
                    & ADZ * ODC(iPrimP,kP1,kP2,0,0,1,3))
               DPLX = DX2*SY0*SZ0
               DPLY = DX1*SY1*SZ0
               DPLZ = DX1*SY0*SZ1
               !X-CM1 X,Y,Z
               INTEGRAL(iPrimP,1,ijk) =  D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
               INTEGRAL(iPrimP,2,ijk) =  D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
               INTEGRAL(iPrimP,3,ijk) =  D2I*(DIFABX*DPLY - DIFABY*DPLX)
               DPLX = SX1*DY1*SZ0
               DPLY = SX0*DY2*SZ0
               DPLZ = SX0*DY1*SZ1
               !Y-CM1 X,Y,Z
               INTEGRAL(iPrimP,4,ijk) =  D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
               INTEGRAL(iPrimP,5,ijk) =  D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
               INTEGRAL(iPrimP,6,ijk) =  D2I*(DIFABX*DPLY - DIFABY*DPLX)
               DPLX = SX1*SY0*DZ1
               DPLY = SX0*SY1*DZ1
               DPLZ = SX0*SY0*DZ2
               !Z-CM1 X,Y,Z
               INTEGRAL(iPrimP,7,ijk) =  D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
               INTEGRAL(iPrimP,8,ijk) =  D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
               INTEGRAL(iPrimP,9,ijk) =  D2I*(DIFABX*DPLY - DIFABY*DPLX)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE DCM1INTGEN

SUBROUTINE DCM1INTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,distanceAB,AODIST,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),distanceAB(3),AODIST(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,SX1,SY1,SZ1,DX1,DY1,DZ1,DX2,DY2,DZ2
real(realk) :: SX2,SY2,SZ2,DIFABX,DIFABY,DIFABZ,ADX,ADY,ADZ,DPLX,DPLY,DPLZ
real(realk),PARAMETER :: D0 = 0.0E0_realk
real(realk),PARAMETER :: D2 = 2.0E0_realk
real(realk),PARAMETER :: D2I = 1.0E0_realk/D2
DO ijk=1,ijkP
   do ip1=1,nOperatorComp
      INTEGRAL(ip1,ijk)=D0
   enddo
ENDDO
DIFABX = distanceAB(1)
DIFABY = distanceAB(2)
DIFABZ = distanceAB(3)
ADX = AODIST(1)
ADY = AODIST(2)
ADZ = AODIST(3)

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
   DO jP2=P2-iP2,0,-1
      kP2=P2-iP2-jP2
      P1=AngmomA
      DO iP1=P1,0,-1
         DO jP1=P1-iP1,0,-1
            kP1=P1-iP1-jP1
            ijk=ijk+1
            DO iPrimP=1,nPrimP
               SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)
               SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
               SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
               SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)
               SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
               SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
               SX2 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,2,1)
               SY2 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,2,2)
               SZ2 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,2,3)
               DX1 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,0,1) + &
                    & ADX * ODC(iPrimP,iP1,iP2,0,0,0,1))
               DY1 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,0,2) + &
                    & ADY * ODC(iPrimP,jP1,jP2,0,0,0,2))
               DZ1 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,0,3) + &
                    & ADZ * ODC(iPrimP,kP1,kP2,0,0,0,3))
               DX2 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,1,1) + &
                    & ADX * ODC(iPrimP,iP1,iP2,0,0,1,1))
               DY2 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,1,2) + &
                    & ADY * ODC(iPrimP,jP1,jP2,0,0,1,2))
               DZ2 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,1,3) + &
                    & ADZ * ODC(iPrimP,kP1,kP2,0,0,1,3))
               DPLX = DX2*SY0*SZ0
               DPLY = DX1*SY1*SZ0
               DPLZ = DX1*SY0*SZ1
               !X-CM1 X,Y,Z
               INTEGRAL(1,ijk)=INTEGRAL(1,ijk) + D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
               INTEGRAL(2,ijk)=INTEGRAL(2,ijk) + D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
               INTEGRAL(3,ijk)=INTEGRAL(3,ijk) + D2I*(DIFABX*DPLY - DIFABY*DPLX)
               DPLX = SX1*DY1*SZ0
               DPLY = SX0*DY2*SZ0
               DPLZ = SX0*DY1*SZ1
               !Y-CM1 X,Y,Z
               INTEGRAL(4,ijk)=INTEGRAL(4,ijk) + D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
               INTEGRAL(5,ijk)=INTEGRAL(5,ijk) + D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
               INTEGRAL(6,ijk)=INTEGRAL(6,ijk) + D2I*(DIFABX*DPLY - DIFABY*DPLX)
               DPLX = SX1*SY0*DZ1
               DPLY = SX0*SY1*DZ1
               DPLZ = SX0*SY0*DZ2
               !Z-CM1 X,Y,Z
               INTEGRAL(7,ijk)=INTEGRAL(7,ijk) + D2I*(DIFABY*DPLZ - DIFABZ*DPLY)
               INTEGRAL(8,ijk)=INTEGRAL(8,ijk) + D2I*(DIFABZ*DPLX - DIFABX*DPLZ)
               INTEGRAL(9,ijk)=INTEGRAL(9,ijk) + D2I*(DIFABX*DPLY - DIFABY*DPLX)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE DCM1INTSEG

SUBROUTINE DCM2INTgen(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,distanceAB,AODIST,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),distanceAB(3),AODIST(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nPrimP,nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,SX1,SY1,SZ1,DX1,DY1,DZ1,DX2,DY2,DZ2
real(realk) :: SX2,SY2,SZ2,DIFABX,DIFABY,DIFABZ,ADX,ADY,ADZ,DPLX,DPLY,DPLZ
real(realk) :: DPLXX,DPLXY,DPLXZ,DPLYY,DPLYZ,DPLZZ,DX3,DY3,DZ3
real(realk),PARAMETER :: D0 = 0.0E0_realk, D4INV = 0.25E0_realk
real(realk),PARAMETER :: D2 = 2.0E0_realk
real(realk),PARAMETER :: D2I = 1.0E0_realk/D2
DIFABX = distanceAB(1)
DIFABY = distanceAB(2)
DIFABZ = distanceAB(3)
ADX = AODIST(1)
ADY = AODIST(2)
ADZ = AODIST(3)
ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1=AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    DO iPrimP=1,nPrimP
     SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)
     SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
     SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
     SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)
     SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
     SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
     SX2 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,2,1)
     SY2 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,2,2)
     SZ2 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,2,3)
     DX1 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,0,1) + &
          & ADX * ODC(iPrimP,iP1,iP2,0,0,0,1))
     DY1 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,0,2) + &
          & ADY * ODC(iPrimP,jP1,jP2,0,0,0,2))
     DZ1 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,0,3) + &
          & ADZ * ODC(iPrimP,kP1,kP2,0,0,0,3))
     DX2 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,1,1) + &
          & ADX * ODC(iPrimP,iP1,iP2,0,0,1,1))
     DY2 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,1,2) + &
          & ADY * ODC(iPrimP,jP1,jP2,0,0,1,2))
     DZ2 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,1,3) + &
          & ADZ * ODC(iPrimP,kP1,kP2,0,0,1,3))
     DX3 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,2,1) + &
          & ADX * ODC(iPrimP,iP1,iP2,0,0,2,1))
     DY3 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,2,2) + &
          & ADY * ODC(iPrimP,jP1,jP2,0,0,2,2))
     DZ3 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,2,3) + &
          & ADZ * ODC(iPrimP,kP1,kP2,0,0,2,3))
     !X-CM2 XX,XY,XZ,YY,YZ,ZZ
     DPLXX = DX3*SY0*SZ0
     DPLXY = DX2*SY1*SZ0
     DPLXZ = DX2*SY0*SZ1
     DPLYY = DX1*SY2*SZ0
     DPLYZ = DX1*SY1*SZ1
     DPLZZ = DX1*SY0*SZ2
     INTEGRAL(iPrimP,1,ijk)=(D2*DIFABZ*DIFABY*DPLYZ-DIFABZ*DIFABZ*DPLYY-DIFABY*DIFABY*DPLZZ)*D4INV
     INTEGRAL(iPrimP,2,ijk)=(DIFABZ*DIFABZ*DPLXY-DIFABY*DIFABZ*DPLXZ-DIFABZ*DIFABX*DPLYZ+DIFABY*DIFABX*DPLZZ)*D4INV
     INTEGRAL(iPrimP,3,ijk)=(DIFABY*DIFABY*DPLXZ-DIFABX*DIFABY*DPLYZ+DIFABZ*DIFABX*DPLYY-DIFABY*DIFABZ*DPLXY)*D4INV
     INTEGRAL(iPrimP,4,ijk)=(D2*DIFABZ*DIFABX*DPLXZ-DIFABZ*DIFABZ*DPLXX-DIFABX*DIFABX*DPLZZ)*D4INV
     INTEGRAL(iPrimP,5,ijk)=(DIFABZ*DIFABY*DPLXX-DIFABZ*DIFABX*DPLXY-DIFABY*DIFABX*DPLXZ+DIFABX*DIFABX*DPLYZ)*D4INV
     INTEGRAL(iPrimP,6,ijk)=(D2*DIFABX*DIFABY*DPLXY-DIFABY*DIFABY*DPLXX-DIFABX*DIFABX*DPLYY)*D4INV
     !Y-CM2 XX,XY,XZ,YY,YZ,ZZ
     DPLXX = SX2*DY1*SZ0
     DPLXY = SX1*DY2*SZ0
     DPLXZ = SX1*DY1*SZ1
     DPLYY = SX0*DY3*SZ0
     DPLYZ = SX0*DY2*SZ1
     DPLZZ = SX0*DY1*SZ2
     INTEGRAL(iPrimP,7,ijk)=(D2*DIFABZ*DIFABY*DPLYZ-DIFABZ*DIFABZ*DPLYY-DIFABY*DIFABY*DPLZZ)*D4INV
     INTEGRAL(iPrimP,8,ijk)=(DIFABZ*DIFABZ*DPLXY-DIFABY*DIFABZ*DPLXZ-DIFABZ*DIFABX*DPLYZ+DIFABY*DIFABX*DPLZZ)*D4INV
     INTEGRAL(iPrimP,9,ijk)=(DIFABY*DIFABY*DPLXZ-DIFABX*DIFABY*DPLYZ+DIFABZ*DIFABX*DPLYY-DIFABY*DIFABZ*DPLXY)*D4INV
     INTEGRAL(iPrimP,10,ijk)=(D2*DIFABZ*DIFABX*DPLXZ-DIFABZ*DIFABZ*DPLXX-DIFABX*DIFABX*DPLZZ)*D4INV
     INTEGRAL(iPrimP,11,ijk)=(DIFABZ*DIFABY*DPLXX-DIFABZ*DIFABX*DPLXY-DIFABY*DIFABX*DPLXZ+DIFABX*DIFABX*DPLYZ)*D4INV
     INTEGRAL(iPrimP,12,ijk)=(D2*DIFABX*DIFABY*DPLXY-DIFABY*DIFABY*DPLXX-DIFABX*DIFABX*DPLYY)*D4INV
     !Z-CM2 XX,XY,XZ,YY,YZ,ZZ
     DPLXX = SX2*SY0*DZ1
     DPLXY = SX1*SY1*DZ1
     DPLXZ = SX1*SY0*DZ2
     DPLYY = SX0*SY2*DZ1
     DPLYZ = SX0*SY1*DZ2
     DPLZZ = SX0*SY0*DZ3
     INTEGRAL(iPrimP,13,ijk)=(D2*DIFABZ*DIFABY*DPLYZ-DIFABZ*DIFABZ*DPLYY-DIFABY*DIFABY*DPLZZ)*D4INV
     INTEGRAL(iPrimP,14,ijk)=(DIFABZ*DIFABZ*DPLXY-DIFABY*DIFABZ*DPLXZ-DIFABZ*DIFABX*DPLYZ+DIFABY*DIFABX*DPLZZ)*D4INV
     INTEGRAL(iPrimP,15,ijk)=(DIFABY*DIFABY*DPLXZ-DIFABX*DIFABY*DPLYZ+DIFABZ*DIFABX*DPLYY-DIFABY*DIFABZ*DPLXY)*D4INV
     INTEGRAL(iPrimP,16,ijk)=(D2*DIFABZ*DIFABX*DPLXZ-DIFABZ*DIFABZ*DPLXX-DIFABX*DIFABX*DPLZZ)*D4INV
     INTEGRAL(iPrimP,17,ijk)=(DIFABZ*DIFABY*DPLXX-DIFABZ*DIFABX*DPLXY-DIFABY*DIFABX*DPLXZ+DIFABX*DIFABX*DPLYZ)*D4INV
     INTEGRAL(iPrimP,18,ijk)=(D2*DIFABX*DIFABY*DPLXY-DIFABY*DIFABY*DPLXX-DIFABX*DIFABX*DPLYY)*D4INV
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE DCM2INTGEN

SUBROUTINE DCM2INTseg(INTEGRAL,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,SHGTF,&
     & AngmomA,AngmomB,ijkP,nOperatorComp,distanceAB,AODIST,lupri)
implicit none
integer     :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
integer     :: AngmomA,AngmomB,ijkP,nOperatorComp,lupri
real(realk) :: SHGTF(nPrimP),distanceAB(3),AODIST(3)
real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
real(realk) :: INTEGRAL(nOperatorComp,ijkP)
!
integer :: ijk,P2,P1,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
real(realk) :: SX0,SY0,SZ0,SX1,SY1,SZ1,DX1,DY1,DZ1,DX2,DY2,DZ2
real(realk) :: SX2,SY2,SZ2,DIFABX,DIFABY,DIFABZ,ADX,ADY,ADZ,DPLX,DPLY,DPLZ
real(realk) :: DPLXX,DPLXY,DPLXZ,DPLYY,DPLYZ,DPLZZ,DX3,DY3,DZ3
real(realk),PARAMETER :: D0 = 0.0E0_realk, D4INV = 0.25E0_realk
real(realk),PARAMETER :: D2 = 2.0E0_realk
real(realk),PARAMETER :: D2I = 1.0E0_realk/D2
DO ijk=1,ijkP
   do ip1=1,nOperatorComp
      INTEGRAL(ip1,ijk)=D0
   enddo
ENDDO
DIFABX = distanceAB(1)
DIFABY = distanceAB(2)
DIFABZ = distanceAB(3)
ADX = AODIST(1)
ADY = AODIST(2)
ADZ = AODIST(3)

ijk=0
P2 = AngmomB
DO iP2=P2,0,-1
 DO jP2=P2-iP2,0,-1
  kP2=P2-iP2-jP2
  P1=AngmomA
  DO iP1=P1,0,-1
   DO jP1=P1-iP1,0,-1
    kP1=P1-iP1-jP1
    ijk=ijk+1
    DO iPrimP=1,nPrimP
     SX0 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,0,1)
     SY0 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,0,2)
     SZ0 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,0,3)
     SX1 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,1,1)
     SY1 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,1,2)
     SZ1 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,1,3)
     SX2 = SHGTF(iPrimP)*ODC(iPrimP,iP1,iP2,0,0,2,1)
     SY2 = SHGTF(iPrimP)*ODC(iPrimP,jP1,jP2,0,0,2,2)
     SZ2 = SHGTF(iPrimP)*ODC(iPrimP,kP1,kP2,0,0,2,3)
     DX1 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,0,1) + &
          & ADX * ODC(iPrimP,iP1,iP2,0,0,0,1))
     DY1 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,0,2) + &
          & ADY * ODC(iPrimP,jP1,jP2,0,0,0,2))
     DZ1 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,0,3) + &
          & ADZ * ODC(iPrimP,kP1,kP2,0,0,0,3))
     DX2 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,1,1) + &
          & ADX * ODC(iPrimP,iP1,iP2,0,0,1,1))
     DY2 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,1,2) + &
          & ADY * ODC(iPrimP,jP1,jP2,0,0,1,2))
     DZ2 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,1,3) + &
          & ADZ * ODC(iPrimP,kP1,kP2,0,0,1,3))
     DX3 = SHGTF(iPrimP)*(ODC(iPrimP,iP1+1,iP2,0,0,2,1) + &
          & ADX * ODC(iPrimP,iP1,iP2,0,0,2,1))
     DY3 = SHGTF(iPrimP)*(ODC(iPrimP,jP1+1,jP2,0,0,2,2) + &
          & ADY * ODC(iPrimP,jP1,jP2,0,0,2,2))
     DZ3 = SHGTF(iPrimP)*(ODC(iPrimP,kP1+1,kP2,0,0,2,3) + &
          & ADZ * ODC(iPrimP,kP1,kP2,0,0,2,3))
     !X-CM2 XX,XY,XZ,YY,YZ,ZZ
     DPLXX = DX3*SY0*SZ0
     DPLXY = DX2*SY1*SZ0
     DPLXZ = DX2*SY0*SZ1
     DPLYY = DX1*SY2*SZ0
     DPLYZ = DX1*SY1*SZ1
     DPLZZ = DX1*SY0*SZ2
     INTEGRAL(1,ijk)=INTEGRAL(1,ijk)+(D2*DIFABZ*DIFABY*DPLYZ-DIFABZ*DIFABZ*DPLYY-DIFABY*DIFABY*DPLZZ)*D4INV
     INTEGRAL(2,ijk)=INTEGRAL(2,ijk)+(DIFABZ*DIFABZ*DPLXY-DIFABY*DIFABZ*DPLXZ-DIFABZ*DIFABX*DPLYZ+DIFABY*DIFABX*DPLZZ)*D4INV
     INTEGRAL(3,ijk)=INTEGRAL(3,ijk)+(DIFABY*DIFABY*DPLXZ-DIFABX*DIFABY*DPLYZ+DIFABZ*DIFABX*DPLYY-DIFABY*DIFABZ*DPLXY)*D4INV
     INTEGRAL(4,ijk)=INTEGRAL(4,ijk)+(D2*DIFABZ*DIFABX*DPLXZ-DIFABZ*DIFABZ*DPLXX-DIFABX*DIFABX*DPLZZ)*D4INV
     INTEGRAL(5,ijk)=INTEGRAL(5,ijk)+(DIFABZ*DIFABY*DPLXX-DIFABZ*DIFABX*DPLXY-DIFABY*DIFABX*DPLXZ+DIFABX*DIFABX*DPLYZ)*D4INV
     INTEGRAL(6,ijk)=INTEGRAL(6,ijk)+(D2*DIFABX*DIFABY*DPLXY-DIFABY*DIFABY*DPLXX-DIFABX*DIFABX*DPLYY)*D4INV
     !Y-CM2 XX,XY,XZ,YY,YZ,ZZ
     DPLXX = SX2*DY1*SZ0
     DPLXY = SX1*DY2*SZ0
     DPLXZ = SX1*DY1*SZ1
     DPLYY = SX0*DY3*SZ0
     DPLYZ = SX0*DY2*SZ1
     DPLZZ = SX0*DY1*SZ2
     INTEGRAL(7,ijk)=INTEGRAL(7,ijk)+(D2*DIFABZ*DIFABY*DPLYZ-DIFABZ*DIFABZ*DPLYY-DIFABY*DIFABY*DPLZZ)*D4INV
     INTEGRAL(8,ijk)=INTEGRAL(8,ijk)+(DIFABZ*DIFABZ*DPLXY-DIFABY*DIFABZ*DPLXZ-DIFABZ*DIFABX*DPLYZ+DIFABY*DIFABX*DPLZZ)*D4INV
     INTEGRAL(9,ijk)=INTEGRAL(9,ijk)+(DIFABY*DIFABY*DPLXZ-DIFABX*DIFABY*DPLYZ+DIFABZ*DIFABX*DPLYY-DIFABY*DIFABZ*DPLXY)*D4INV
     INTEGRAL(10,ijk)=INTEGRAL(10,ijk)+(D2*DIFABZ*DIFABX*DPLXZ-DIFABZ*DIFABZ*DPLXX-DIFABX*DIFABX*DPLZZ)*D4INV
     INTEGRAL(11,ijk)=INTEGRAL(11,ijk)+(DIFABZ*DIFABY*DPLXX-DIFABZ*DIFABX*DPLXY-DIFABY*DIFABX*DPLXZ+DIFABX*DIFABX*DPLYZ)*D4INV
     INTEGRAL(12,ijk)=INTEGRAL(12,ijk)+(D2*DIFABX*DIFABY*DPLXY-DIFABY*DIFABY*DPLXX-DIFABX*DIFABX*DPLYY)*D4INV
     !Z-CM2 XX,XY,XZ,YY,YZ,ZZ
     DPLXX = SX2*SY0*DZ1
     DPLXY = SX1*SY1*DZ1
     DPLXZ = SX1*SY0*DZ2
     DPLYY = SX0*SY2*DZ1
     DPLYZ = SX0*SY1*DZ2
     DPLZZ = SX0*SY0*DZ3
     INTEGRAL(13,ijk)=INTEGRAL(13,ijk)+(D2*DIFABZ*DIFABY*DPLYZ-DIFABZ*DIFABZ*DPLYY-DIFABY*DIFABY*DPLZZ)*D4INV
     INTEGRAL(14,ijk)=INTEGRAL(14,ijk)+(DIFABZ*DIFABZ*DPLXY-DIFABY*DIFABZ*DPLXZ-DIFABZ*DIFABX*DPLYZ+DIFABY*DIFABX*DPLZZ)*D4INV
     INTEGRAL(15,ijk)=INTEGRAL(15,ijk)+(DIFABY*DIFABY*DPLXZ-DIFABX*DIFABY*DPLYZ+DIFABZ*DIFABX*DPLYY-DIFABY*DIFABZ*DPLXY)*D4INV
     INTEGRAL(16,ijk)=INTEGRAL(16,ijk)+(D2*DIFABZ*DIFABX*DPLXZ-DIFABZ*DIFABZ*DPLXX-DIFABX*DIFABX*DPLZZ)*D4INV
     INTEGRAL(17,ijk)=INTEGRAL(17,ijk)+(DIFABZ*DIFABY*DPLXX-DIFABZ*DIFABX*DPLXY-DIFABY*DIFABX*DPLXZ+DIFABX*DIFABX*DPLYZ)*D4INV
     INTEGRAL(18,ijk)=INTEGRAL(18,ijk)+(D2*DIFABX*DIFABY*DPLXY-DIFABY*DIFABY*DPLXX-DIFABX*DIFABX*DPLYY)*D4INV
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE DCM2INTSEG

SUBROUTINE buildpropEcoeff(ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,&
     & JMAXAP,JMAXBP,JMAXTP,DIFODC,KINODC,EXPA,EXPB,preexpfacinv3,EXPPI,&
     & PADIST,PBDIST,DOMOM1,ORIGIN,PROPTYPE,AODIST,BODIST,IPRINT,LUPRI)
implicit none
Integer             :: JmaxA,JmaxB,JmaxP,JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,JMAXT
Integer             :: nPrimP,IPRINT,LUPRI,PROPTYPE
logical             :: DIFODC,KINODC,DOMOM1
!maybe EXPA(nprimA) should be expa(nprimP) ??
REAL(REALK)         :: ORIGIN(3),AODIST(3),BODIST(3)!passes modifies this
REAL(REALK)         :: EXPPI(nPrimP),preexpfacinv3(nPrimP),FAC,TEXPB1,EXPA(nPrimP),EXPB(nPrimP),FACB,CORAX,CORAY,CORAZ,SHGTF
REAL(REALK)         :: INT,OAX,OAY,OAZ,PAX,PAY,PAZ,EXPP,PBX,PBY,PBZ,EXPPIH,PADIST(nPrimP,3),PBDIST(nPrimP,3)
Integer             :: A,B,T,IC,AB,IT,D,X,IM,M,IA,IB,ID,ITEX,n1,n2,iPrimP
Real(realk),pointer :: ODCUND(:,:,:,:,:,:,:)
Real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
integer (kind=long) :: nsize

n1=nPrimP*3*(JMAXA+1)*(JMAXB+1)*(JMAXT+1)*(JMAXD+1)*(JMAXM+1)
CALL LS_DZERO(ODC,n1)
!=========================================
allocate(ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3))
nsize = size(ODCUND)*mem_realsize
call mem_allocated_mem_real(nsize)
!=========================================

n2=nPrimP*3*(JMAXAP+4)*(JMAXBP+3)*(JMAXTP+3)*(JMAXD+1)*(JMAXM+1)
CALL LS_DZERO(ODCUND,n2)

CALL ONEODC(ODCUND,JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,preexpfacinv3,nPrimP,&
     & EXPPI,PADIST,PBDIST,IPRINT,LUPRI)

CALL COPYODC1(ODC,ODCUND,nPrimP,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,&
     & JMAXTP,JMAXD,JMAXM,lupri)
!
!  Expansion coefficients for derivatives
!
IF(DIFODC) THEN
   IF((PROPTYPE.EQ.19.OR.PROPTYPE.EQ.190).OR.(PROPTYPE.EQ.27.OR.PROPTYPE.EQ.26))THEN
      CALL DODCB(ODC,ODCUND,nPrimP,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,&
           &                    JMAXTP,JMAXD,JMAXM,EXPB,IPRINT,LUPRI)
      CALL COPYODC2(ODC,ODCUND,nPrimP,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,&
          & JMAXTP,JMAXD,JMAXM)
   ELSE
      CALL DODCA(ODC,ODCUND,nPrimP,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,&
           &                    JMAXTP,JMAXD,JMAXM,EXPA,IPRINT,LUPRI)
      CALL COPYODC2(ODC,ODCUND,nPrimP,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,&
           & JMAXTP,JMAXD,JMAXM)
   ENDIF
ENDIF

!
!  Expansion coeffiecients for moments or electric derivatives
!
IF (JMAXM .GT. 0) THEN
   IF((PROPTYPE .EQ. 49) .OR. (PROPTYPE .EQ. 51)) THEN
!      CALL EFODC(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,&
!      &                 JMAXTP,JMAXD,JMAXM,EXPA,EXPB,IPRINT)
   ELSEIF(PROPTYPE .EQ. 63 .OR. PROPTYPE .EQ. 83) THEN
!      CALL DODCB2(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,&
!      &                  JMAXTP,JMAXD,JMAXM,EXPB,IPRINT)
   ELSEIF(DOMOM1) THEN
      IF(PROPTYPE.EQ.55)THEN
         CALL MODCB(ODC,ODCUND,nPrimP,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,&
              &           JMAXTP,JMAXD,JMAXM,BODIST,IPRINT,LUPRI)
      ELSE
         CALL MODCA(ODC,ODCUND,nPrimP,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,&
              &           JMAXTP,JMAXD,JMAXM,AODIST,IPRINT,LUPRI)
      END IF
      CALL COPYODC3(ODC,ODCUND,nPrimP,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,&
           & JMAXTP,JMAXD,JMAXM)
   ENDIF
END IF

IF(IPRINT.GT.100)THEN
   WRITE(lupri,*)'buildpropEcoeff'
   CALL WRITE_ODC(ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,LUPRI)
ENDIF
!================================================
nsize = size(ODCUND)*mem_realsize
call mem_deallocated_mem_real(nsize)
deallocate(ODCUND)
!================================================

end SUBROUTINE buildpropEcoeff

SUBROUTINE ONEODC(ODCUND,JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,FAC,&
     &  nPrimP,EXPPI,PADIST,PBDIST,IPRINT,LUPRI)
implicit none
Integer     :: JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,nPrimP,IPRINT,LUPRI
REAL(REALK) :: EXPPI(nPrimP),FAC(nPrimP),PADIST(nPrimP,3),PBDIST(nPrimP,3)
Real(realk) :: ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
!
Integer     :: A,B,T,AB,X,iprimP
Real(realk) :: Tplus1,EXPPIH(nPrimP)
Real(realk),parameter :: D05=0.5E0_realk

DO IPRIMP = 1,nPrimP
   EXPPIH(IPRIMP) = EXPPI(IPRIMP)*D05
ENDDO

DO X=1,3
   A=0
   DO IPRIMP = 1,nPrimP   
      ODCUND(IprimP,0,0,0,0,0,X) = FAC(iPrimP)
   ENDDO
   DO B = 1, JMAXBP
      !   AB = B
      DO T = 0, B!AB
         Tplus1 = (T+1)
         DO IPRIMP = 1,nPrimP
            ODCUND(IPRIMP,0,B,T,0,0,X) = EXPPIH(IPRIMP)*ODCUND(IPRIMP,0,B-1,T-1,0,0,X)&
                 & + PBDIST(IPRIMP,X)*ODCUND(IPRIMP,0,B-1,T,0,0,X) &
                 & + Tplus1*ODCUND(IPRIMP,0,B-1,T+1,0,0,X)
         ENDDO
      ENDDO
   ENDDO
   DO A=1,JMAXAP
      DO T = 0, A
         Tplus1 = (T+1)
         DO IPRIMP = 1,nPrimP
            ODCUND(IPRIMP,A,0,T,0,0,X) = EXPPIH(IPRIMP)*ODCUND(IPRIMP,A-1,0,T-1,0,0,X)&
           & + PADIST(IPRIMP,X)*ODCUND(IPRIMP,A-1,0,T ,0,0,X)+  Tplus1*ODCUND(IPRIMP,A-1,0,T+1,0,0,X)
         ENDDO
      ENDDO
      DO B = 1, JMAXBP
         AB = A + B
         DO T = 0, AB
            Tplus1 = (T+1)
            DO IPRIMP = 1,nPrimP
               ODCUND(IPRIMP,A,B,T,0,0,X) = EXPPIH(IPRIMP)*ODCUND(IPRIMP,A,B-1,T-1,0,0,X)&
                    & + PBDIST(IPRIMP,X)*ODCUND(IPRIMP,A,B-1,T  ,0,0,X) + Tplus1*ODCUND(IPRIMP,A,B-1,T+1,0,0,X)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE ONEODC

SUBROUTINE COPYODC1(ODC,ODCUND,nPrimP,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,&
     & JMAXTP,JMAXD,JMAXM,lupri)
implicit none
Integer             :: JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM,nPrimP,lupri
Real(realk) :: ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
Real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
!
INTEGER :: IC,IT,IB,IA,iPrimP

DO IC = 1, 3
   DO IT = 0, JMAXT
      DO IB = 0, JMAXB
         DO IA = 0, JMAXA
            DO IPRIMP = 1,nPrimP
               ODC(IPRIMP,IA,IB,IT,0,0,IC) = ODCUND(IPRIMP,IA,IB,IT,0,0,IC)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE COPYODC1

SUBROUTINE WRITE_ODC(ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,LUPRI)
implicit none
integer :: nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,LUPRI
Real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
!
Integer :: IX,M,D,T,B,A,iprimp
DO M=0,JMAXM
 DO D=0,JMAXD
  IF(M.EQ.0.AND.D.EQ.0)THEN
     WRITE(lupri,'(3X,A)')'Undifferentiated Ecoefficients'
  ELSE 
     WRITE(lupri,'(I3,A)') D, '. order differentiated Ecoefficients'
     WRITE(lupri,'(I3,A)') M, '. order moment Ecoefficients  '
  ENDIF
  DO IX=1,3
   WRITE(lupri,'(A,I3)') 'Cartesian Component=',IX
   WRITE(lupri,'(A3,A3,A3)')'A','B','T'
   DO T=0,JMAXT
    DO B=0,JMAXB
     DO A=0,JMAXA
        WRITE(lupri,'(I3,I3,I3,5ES15.5/,(9X,5ES15.5))') A,B,T,(ODC(iPrimP,A,B,T,D,M,IX),iPrimP=1,nPrimP)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO               
END SUBROUTINE WRITE_ODC

!Eq 9.3.28
SUBROUTINE DODCA(ODC,ODCUND,nPrimP,JMAXA,JMAXAP,JMAXBP,&
     & JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM,EXPA,IPRINT,LUPRI)
  implicit none
  integer :: JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM
  integer :: IPRINT,LUPRI,nPrimP
  real(realk) :: EXPA(nPrimP)
  real(realk) :: ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
  real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
  !
  INTEGER A, B, D, T, X,iPrimP
!WARNING STATIC ALLOCATION
  REAL(REALK) :: EXPA2(nPrimP),FACA
  DO iPrimP=1,nPrimP
     EXPA2(iPrimP) = EXPA(iPrimP) + EXPA(iPrimP)
  ENDDO
  DO X = 1, 3
   DO D = 1, JMAXD
    DO A = 0, JMAXAP
     FACA = A
     DO B = 0, JMAXBP
      DO T = 0, A + B + D
       DO iPrimP=1,nPrimP
        ODCUND(iPrimP,A,B,T,D,0,X)=EXPA2(iPrimP)*ODCUND(iPrimP,A+1,B,T,D-1,0,X)&
             & - FACA*ODCUND(iPrimP,A-1,B,T,D-1,0,X)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
END SUBROUTINE DODCA

SUBROUTINE DODCB(ODC,ODCUND,nPrimP,JMAXA,JMAXAP,JMAXBP,&
     & JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM,EXPB,IPRINT,LUPRI)
  implicit none
  integer :: JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM
  integer :: IPRINT,LUPRI,nPrimP
  real(realk) :: EXPB(nPrimP)
  real(realk) :: ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
  real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
  !
  INTEGER A, B, D, T, X,iPrimP
!WARNING STATIC ALLOCATION
  REAL(REALK) :: EXPB2(nPrimP),FACB
  DO iPrimP=1,nPrimP
     EXPB2(iPrimP) = EXPB(iPrimP) + EXPB(iPrimP)
  ENDDO
  DO X = 1, 3
   DO D = 1, JMAXD
    DO A = 0, JMAXAP
     DO B = 0, JMAXBP
      FACB = B
      DO T = 0, A + B + D
       DO iPrimP=1,nPrimP
          ODCUND(iPrimP,A,B,T,D,0,X) = EXPB2(iPrimP)*ODCUND(iPrimP,A,B+1,T,D-1,0,X)&
               & - FACB*ODCUND(iPrimP,A,B - 1,T,D-1,0,X)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
END SUBROUTINE DODCB

SUBROUTINE COPYODC2(ODC,ODCUND,nPrimP,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,&
     & JMAXTP,JMAXD,JMAXM)
implicit none
Integer             :: JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM,nPrimP
Real(realk) :: ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
Real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
!
INTEGER :: IC,IT,IB,IA,ID,iPrimP
DO IC = 1, 3
   DO ID = 0, JMAXD
      DO IT = 0, JMAXT
         DO IB = 0, JMAXB
            DO IA = 0, JMAXA
               DO iPrimP=1,nPrimP
                  ODC(iPrimP,IA,IB,IT,ID,0,IC) = ODCUND(iPrimP,IA,IB,IT,ID,0,IC)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE COPYODC2

!Eq 9.3.16 
subroutine MODCA(ODC,ODCUND,nPrimP,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,&
  & JMAXTP,JMAXD,JMAXM,AODIST,IPRINT,LUPRI)
  implicit none
  Integer     :: JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM,nPrimP,IPRINT,LUPRI
  real(realk) :: AODIST(3)
  Real(realk) :: ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
  Real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
  !
  INTEGER A, B, T, D, X, M,iPrimP
  REAL(REALK) :: FAC
  DO X=1,3
     FAC = AODIST(X)
     DO M = 1, JMAXM
        DO A = 0, JMAXAP
           DO B = 0, JMAXBP
              DO D = 0, JMAXD
                 DO T = 0, A + B + D + M
                    DO iPrimP=1,nPrimP
                       ODCUND(iPrimP,A,B,T,D,M,X) = ODCUND(iPrimP,A + 1,B,T,D,M-1,X) &
                     & + FAC*ODCUND(iPrimP,A,B,T,D,M-1,X)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
END subroutine MODCA

subroutine MODCB(ODC,ODCUND,nPrimP,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,&
  & JMAXTP,JMAXD,JMAXM,BODIST,IPRINT,LUPRI)
  implicit none
  Integer     :: JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM,nPrimP,IPRINT,LUPRI
  real(realk) :: BODIST(3)
  Real(realk) :: ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
  Real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
  !
  INTEGER A, B, T, D, X, M,iPrimP
  REAL(REALK) :: FAC
  DO X=1,3
     FAC = BODIST(X)
     DO M = 1, JMAXM
        DO A = 0, JMAXAP
           DO B = 0, JMAXBP
              DO D = 0, JMAXD
                 DO T = 0, A + B + D + M
                    DO iPrimP=1,nPrimP
                       ODCUND(iPrimP,A,B,T,D,M,X) = ODCUND(iPrimP,A,B+1,T,D,M-1,X) &
                     & + FAC*ODCUND(iPrimP,A,B,T,D,M-1,X)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
END subroutine MODCB

SUBROUTINE COPYODC3(ODC,ODCUND,nPrimP,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,&
     & JMAXTP,JMAXD,JMAXM)
implicit none
Integer             :: JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,JMAXD,JMAXM,nPrimP
Real(realk) :: ODCUND(nPrimP,-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
Real(realk) :: ODC(nPrimP,0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
!
INTEGER :: IC,IT,IB,IA,ID,IM,iPrimP
DO IC = 1, 3
   DO IM = 0, JMAXM
      DO ID = 0, JMAXD
         DO IT = 0, JMAXT
            DO IB = 0, JMAXB
               DO IA = 0, JMAXA
                  DO iPrimP=1,nPrimP
                     ODC(iPrimP,IA,IB,IT,ID,IM,IC) = ODCUND(iPrimP,IA,IB,IT,ID,IM,IC)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE COPYODC3

!> \brief new distributePQ to lstensor
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param RES contains the result lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE DistributePropIntegrals(RES,QPmat2,nA,nB,nOperatorComp,PQ,Input,Lsoutput,LUPRI,IPRINT)
implicit none 
Type(integrand)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput !not used
Integer,intent(in)              :: LUPRI,IPRINT,nA,nB,nOperatorComp
REAL(REALK),target,intent(in)   :: QPMAT2(nA,nB,nOperatorComp)  !nA,nB
TYPE(lstensor),intent(inout)    :: RES
!
TYPE(overlap),pointer :: P,Q
logical :: SameLHSaos,antiAB,add,permuteAB,AntipermuteAB
integer :: iPassP,atomA,batchA,atomB,batchB,indAB,indBA,atomC,n1,n2,sA,sB
integer :: nFullOperatorComp,startOp,nPass,nOperatorComp2,iPassQ,maxbat,maxang
real(realk),pointer :: AB(:),BA(:)
P => PQ%P%p
Q => PQ%Q%p

SameLHSaos     = INPUT%SameLHSaos
antiAB=Input%PropAnti
startOP = 0
add = INPUT%AddToIntegral
nFullOperatorComp = Lsoutput%ndim(5)
!nOperatorComp is actually nPassQ*nOperatorComp
IF(Q%orbital1%TYPE_Nucleus.AND.(.NOT.add))then
   nPass = Q%nPasses
   nOperatorComp2 = nOperatorComp/nPass
ELSE
   nPass = 1
   nOperatorComp2 = nOperatorComp
ENDIF

atomA  = P%orb1atom(1)
batchA = P%orb1batch(1)
atomB  = P%orb2atom(1)
batchB = P%orb2batch(1)
permuteAB = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
AntipermuteAB  = permuteAB.AND.antiAB
   
indAB = res%index(atomA,atomB,1,1)
!AB => res%LSAO(indAB)%batch(batchA,batchB,1,1)%elms
AB => res%LSAO(indAB)%elms
n1 = res%LSAO(indAB)%nLocal(1)
n2 = res%LSAO(indAB)%nLocal(2)
maxBat = res%LSAO(indAB)%maxbat
maxAng = res%LSAO(indAB)%maxAng
sA = res%LSAO(indAB)%startLocalOrb(1+(batchA-1)*maxAng) - 1
sB = res%LSAO(indAB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat) - 1 
IF (permuteAB) THEN
   indBA = res%index(atomB,atomA,1,1)
   BA => res%LSAO(indBA)%elms
ENDIF
DO iPassQ = 1,nPass
   IF(Q%orbital1%TYPE_Nucleus.AND.(.NOT.add))then
      atomC  = Q%orb1atom(1)
      startOP = (atomC-1)*nOperatorComp2
      IF(nFullOperatorComp.EQ.nOperatorComp2.AND.(atomC.GT.1))THEN
         !for setting%scheme%AONuclearSpecID .NE. 0 we only calculate the 
         !property integral for a specific atom (therefor the nFullOperatorComp=nOperatorComp)
         startOP=0
      ENDIF
   ENDIF
   IF (antipermuteAB) THEN
      CALL intAB(AB,n1,n2,sA,sB,nA,nB,QPmat2,add,nOperatorComp2,nPass,iPassQ,nFullOperatorComp,startOP,lupri)
      CALL AntiIntBA(BA,n1,n2,sA,sB,nA,nB,QPmat2,add,nOperatorComp2,nPass,iPassQ,nFullOperatorComp,startOP,lupri)
   ELSE IF (permuteAB) THEN
      CALL intAB(AB,n1,n2,sA,sB,nA,nB,QPmat2,add,nOperatorComp2,nPass,iPassQ,nFullOperatorComp,startOP,lupri)
      CALL intBA(BA,n1,n2,sA,sB,nA,nB,QPmat2,add,nOperatorComp2,nPass,iPassQ,nFullOperatorComp,startOP,lupri)
   ELSE
      CALL intAB(AB,n1,n2,sA,sB,nA,nB,QPmat2,add,nOperatorComp2,nPass,iPassQ,nFullOperatorComp,startOP,lupri)
   ENDIF
ENDDO
CONTAINS
SUBROUTINE intAB(AB,n1,n2,sA,sB,nA,nB,integralsAB,add,nOperatorComp2,nPass,iPass,nFullOperatorComp,startOP,lupri)
implicit none
Integer,intent(IN)       :: nA,nB,nOperatorComp2,lupri,nPass,iPass,nFullOperatorComp,startOP
Integer,intent(IN)       :: n1,n2,sA,sB
Real(realk),intent(IN)   :: integralsAB(nA,nB,nOperatorComp2,nPass)
Real(realk),intent(INOUT):: AB(n1,n2,nFullOperatorComp)!(nA*nB*nFullOperatorComp)
Logical,intent(IN)       :: add
!
Integer :: X,iB,iA
IF (.not.add) THEN
   do X=1,nOperatorComp2
      do iB=1,nB
         do iA=1,nA
            AB(sA+iA,sB+iB,startOP+X) = integralsAB(iA,iB,X,iPass)
         enddo
      enddo
   enddo
ELSE
!$OMP CRITICAL (PropdistributeOMP)
   do X=1,nOperatorComp2
      do iB=1,nB
         do iA=1,nA
            AB(sA+iA,sB+iB,startOP+X) = AB(sA+iA,sB+iB,startOP+X) + integralsAB(iA,iB,X,iPass)
         enddo
      enddo
   enddo
!$OMP END CRITICAL (PropdistributeOMP)
ENDIF
END SUBROUTINE intAB

SUBROUTINE intBA(BA,n1,n2,sA,sB,nA,nB,integralsAB,add,nOperatorComp2,nPass,iPass,nFullOperatorComp,startOP,lupri)
implicit none
Integer,intent(IN)       :: nA,nB,nOperatorComp2,lupri,nPass,iPass,nFullOperatorComp,startOP
Integer,intent(IN)       :: n1,n2,sA,sB
Real(realk),intent(IN)   :: integralsAB(nA,nB,nOperatorComp2,nPass)
Real(realk),intent(INOUT):: BA(n2,n1,nFullOperatorComp)!(nB,nA,nFullOperatorComp)
Logical,intent(IN)       :: add
!
Integer :: X,iB,iA
IF (.not.add) THEN
   do X=1,nOperatorComp2
      do iB=1,nB
         do iA=1,nA
            BA(sB+iB,sA+iA,startOP+X) = integralsAB(iA,iB,X,iPass)
         ENDDO
      ENDDO
   ENDDO
ELSE
!$OMP CRITICAL (PropdistributeOMP)
   do X=1,nOperatorComp2
      do iB=1,nB
         do iA=1,nA
            BA(sB+iB,sA+iA,startOP+X) = BA(sB+iB,sA+iA,startOP+X) + integralsAB(iA,iB,X,iPass)
         ENDDO
      ENDDO
   ENDDO
!$OMP END CRITICAL (PropdistributeOMP)
ENDIF
END SUBROUTINE intBA

SUBROUTINE AntiintBA(BA,n1,n2,sA,sB,nA,nB,integralsAB,add,nOperatorComp2,nPass,iPass,nFullOperatorComp,startOP,lupri)
implicit none
Integer,intent(IN)       :: nA,nB,nOperatorComp2,lupri,nPass,iPass,nFullOperatorComp,startOP
Integer,intent(IN)       :: n1,n2,sA,sB
Real(realk),intent(IN)   :: integralsAB(nA,nB,nOperatorComp2,nPass)
Real(realk),intent(INOUT):: BA(n2,n1,nFullOperatorComp)!(nB,nA,nFullOperatorComp)
Logical,intent(IN)       :: add
!
Integer :: X,iB,iA
IF (.not.add) THEN
   do X=1,nOperatorComp2
      do iB=1,nB
         do iA=1,nA
            BA(sB+iB,sA+iA,startOP+X) = -integralsAB(iA,iB,X,iPass)
         ENDDO
      ENDDO
   ENDDO
ELSE
!$OMP CRITICAL (PropdistributeOMP)
   do X=1,nOperatorComp2
      do iB=1,nB
         do iA=1,nA
            BA(sB+iB,sA+iA,startOP+X) = BA(sB+iB,sA+iA,startOP+X) - integralsAB(iA,iB,X,iPass)
         ENDDO
      ENDDO
   ENDDO
!$OMP END CRITICAL (PropdistributeOMP)
ENDIF
END SUBROUTINE AntiintBA

END SUBROUTINE DistributePropIntegrals

!> \brief new distributePQ to lstensor
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param RES contains the result lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE DistributePropExpVal(expval,QPmat2,Dlhs,nA,nB,nOperatorComp,PQ,Input,Lsoutput,LUPRI,IPRINT)
implicit none 
Type(integrand)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput !not used
Integer,intent(in)              :: LUPRI,IPRINT,nA,nB,nOperatorComp
REAL(REALK),target,intent(in)   :: QPMAT2(nA,nB,nOperatorComp)  !nA,nB
TYPE(lstensor),intent(in)       :: Dlhs
Real(Realk),intent(inout)       :: expval(:)
!
TYPE(overlap),pointer :: P,Q
logical :: SameLHSaos,antiAB,permuteAB,AntipermuteAB
integer :: iPassP,atomA,batchA,atomB,batchB,DmatindAB,DmatindBA,atomC,n1,n2,sA,sB
integer :: nFullOperatorComp,startOp,nPass,nOperatorComp2,iPassQ,I,ndmat,idmat
integer :: maxbat,maxAng
real(realk),pointer :: DmatAB(:),DmatBA(:)
real(realk),pointer :: TmpExpVal(:,:)
P => PQ%P%p
Q => PQ%Q%p
ndmat = Input%NDMAT_LHS
SameLHSaos     = INPUT%SameLHSaos
antiAB=Input%PropAnti
!startOP = 0
nFullOperatorComp = Lsoutput%ndim(5)/ndmat
!nOperatorComp is actually nPassQ*nOperatorComp
nPass = Q%nPasses
nOperatorComp2 = nOperatorComp/nPass

!WRITE(lupri,*)'DistributePropIntegrals (1,1) nOperatorComp',nOperatorComp
!DO ipassQ=1,nOperatorComp
!   WRITE(lupri,*)'DistributePropIntegrals (1,1) iOperatorComp',ipassQ
!   call ls_output(QPMAT2(:,:,ipassQ),1,nA,1,nB,nA,nB,1,lupri)
!ENDDO

atomA  = P%orb1atom(1)
batchA = P%orb1batch(1)
atomB  = P%orb2atom(1)
batchB = P%orb2batch(1)
permuteAB = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
AntipermuteAB  = permuteAB.AND.antiAB
call mem_alloc(TmpExpVal,nOperatorComp2,ndmat)
DmatindAB = Dlhs%index(atomA,atomB,1,1)
DmatAB => Dlhs%LSAO(DmatindAB)%elms
n1 = Dlhs%LSAO(DmatindAB)%nLocal(1)
n2 = Dlhs%LSAO(DmatindAB)%nLocal(2)
maxBat = Dlhs%LSAO(DmatindAB)%maxbat
maxAng = Dlhs%LSAO(DmatindAB)%maxAng
sA = Dlhs%LSAO(DmatindAB)%startLocalOrb(1+(batchA-1)*maxAng) - 1
sB = Dlhs%LSAO(DmatindAB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat) - 1 
IF (permuteAB) THEN
   DmatindBA = Dlhs%index(atomB,atomA,1,1)
   DmatBA => Dlhs%LSAO(DmatindBA)%elms
ENDIF
DO iPassQ = 1,nPass
   DO IDMAT=1,NDMAT
      DO I = 1,nOperatorComp2
         TmpExpVal(I,IDMAT) = 0E0_realk
      ENDDO
   ENDDO
   IF (antipermuteAB) THEN
      CALL intAB(TmpExpval,DmatAB,n1,n2,sA,sB,ndmat,nA,nB,QPmat2,nOperatorComp2,nPass,iPassQ,lupri)
      CALL AntiIntBA(TmpExpval,DmatBA,n1,n2,sA,sB,ndmat,nA,nB,QPmat2,nOperatorComp2,nPass,iPassQ,lupri)
   ELSE IF (permuteAB) THEN
      CALL intAB(TmpExpval,DmatAB,n1,n2,sA,sB,ndmat,nA,nB,QPmat2,nOperatorComp2,nPass,iPassQ,lupri)
      CALL intBA(TmpExpval,DmatBA,n1,n2,sA,sB,ndmat,nA,nB,QPmat2,nOperatorComp2,nPass,iPassQ,lupri)
   ELSE
      CALL intAB(TmpExpval,DmatAB,n1,n2,sA,sB,ndmat,nA,nB,QPmat2,nOperatorComp2,nPass,iPassQ,lupri)
   ENDIF   
   atomC  = Q%orb1atom(iPassQ)
!$OMP CRITICAL (ProPExpVal)
   DO IDMAT=1,NDMAT
      startOP = (atomC-1)*nOperatorComp2 + (IDMAT-1)*nFullOperatorComp
      IF(nFullOperatorComp.EQ.nOperatorComp2.AND.(atomC.GT.1))THEN
         !for setting%scheme%AONuclearSpecID .NE. 0 we only calculate the 
         !property integral for a specific atom (therefor the nFullOperatorComp=nOperatorComp)
         startOP=(IDMAT-1)*nFullOperatorComp
      ENDIF
      DO I = 1,nOperatorComp2
         ExpVal(startOP+I) = ExpVal(startOP+I) + TmpExpVal(I,idmat)
      ENDDO
   ENDDO
!$OMP END CRITICAL (ProPExpVal)
ENDDO
call mem_dealloc(TmpExpVal)
CONTAINS
SUBROUTINE intAB(ExpVal,DmatAB,n1,n2,sA,sB,ndmat,nA,nB,integralsAB,nOperatorComp2,nPass,iPass,lupri)
implicit none
Integer,intent(IN)       :: ndmat,nA,nB,nOperatorComp2,lupri,nPass,iPass,n1,n2,sA,sB
Real(realk),intent(IN)   :: integralsAB(nA,nB,nOperatorComp2,nPass)
Real(realk),intent(IN)   :: DmatAB(n1,n2,ndmat) !(nA*nB,ndmat)
Real(realk),intent(INOUT):: Expval(nOperatorComp2,ndmat)
!
Integer :: X,IA,iB,idmat
real(realk) :: res
do Idmat=1,ndmat
   do X=1,nOperatorComp2
      res = 0.0E0_realk
      do IB=1,nB
       do IA=1,nA
         res = res + DmatAB(sA+IA,sB+IB,idmat)*integralsAB(IA,IB,X,iPass)
       ENDDO
      ENDDO
      Expval(X,idmat) = res
   enddo
enddo
END SUBROUTINE intAB

SUBROUTINE intBA(ExpVal,DmatBA,n1,n2,sA,sB,ndmat,nA,nB,integralsAB,nOperatorComp2,nPass,iPass,lupri)
implicit none
Integer,intent(IN)       :: ndmat,nA,nB,nOperatorComp2,lupri,nPass,iPass,n1,n2,sA,sB
Real(realk),intent(IN)   :: integralsAB(nA,nB,nOperatorComp2,nPass)
Real(realk),intent(IN)   :: DmatBA(n2,n1,ndmat) !(nB,nA,ndmat)
Real(realk),intent(INOUT):: Expval(nOperatorComp2,ndmat)
!
Integer :: X,iB,iA,idmat
real(realk) :: res
do Idmat=1,ndmat
   do X=1,nOperatorComp2
      res = 0.0E0_realk
      do iB=1,nB
         do iA=1,nA
            res = res + DmatBA(sB+iB,sA+iA,idmat)*integralsAB(iA,iB,X,iPass)
         ENDDO
      ENDDO
      Expval(X,idmat) = res
   enddo
ENDDO
END SUBROUTINE intBA

SUBROUTINE AntiintBA(ExpVal,DmatBA,n1,n2,sA,sB,ndmat,nA,nB,integralsAB,nOperatorComp2,nPass,iPass,lupri)
implicit none
Integer,intent(IN)       :: ndmat,nA,nB,nOperatorComp2,lupri,nPass,iPass,n1,n2,sA,sB
Real(realk),intent(IN)   :: integralsAB(nA,nB,nOperatorComp2,nPass)
Real(realk),intent(IN)   :: DmatBA(n2,n1,ndmat) !(nB,nA,ndmat)
Real(realk),intent(INOUT):: Expval(nOperatorComp2,ndmat)
!
Integer :: X,iB,iA,idmat
real(realk) :: res
do Idmat=1,ndmat
   do X=1,nOperatorComp2
      res = 0.0E0_realk
      do iB=1,nB
         do iA=1,nA
            res = res - DmatBA(sB+iB,sA+iA,idmat)*integralsAB(iA,iB,X,iPass)
         enddo
      enddo
      Expval(X,idmat) = res
   enddo
enddo
END SUBROUTINE AntiintBA

END SUBROUTINE DistributePropExpVal

END MODULE Thermite_prop
