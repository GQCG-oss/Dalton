MODULE AGC_GPU_OBS_HorizontalRecurrenceRHSModCtoD
 use IchorPrecisionMod
  
 CONTAINS

!Transfer angmom from C to D
subroutine HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri,iASync)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,    4)
  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,    2:    4,    1:    1)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  integer :: iP,ilmP,iTUVC
!  real(realk) :: Tmp(nTUVA,nTUVB) ordering
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVC,ilmP) &
!$ACC PRESENT(nPasses,ThetaP,ThetaP2) ASYNC(iASync)
  DO iP = 1,nContPQ*nPasses
    DO iTUVC=  2,  4
     DO ilmP = 1,nlmP
        ThetaP(IP,ilmP,iTUVC,1) = ThetaP2(IP,ilmP,iTUVC)
     ENDDO
    ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_RHS_Q1C1D0CtoD

!Transfer angmom from C to D
subroutine HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri,iASync)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,   10)
  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,    2:    4,    2:    4)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVC
  real(realk) :: Xcd,Ycd,Zcd
!  real(realk) :: Tmp(nTUVA,nTUVB) ordering
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,&
!$ACC         iTUVC,ilmP,Xcd,Ycd,Zcd) &
!$ACC PRESENT(nPasses,Qdistance12,ThetaP,ThetaP2) ASYNC(iASync)
  DO iP = 1,nContPQ*nPasses
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
    DO ilmP = 1,nlmP
     ThetaP(IP,ilmP,2,2) = ThetaP2(IP,ilmP,5) + Xcd*ThetaP2(iP,ilmP,2) 
     ThetaP(IP,ilmP,3,2) = ThetaP2(IP,ilmP,6) + Xcd*ThetaP2(iP,ilmP,3) 
     ThetaP(IP,ilmP,4,2) = ThetaP2(IP,ilmP,7) + Xcd*ThetaP2(iP,ilmP,4) 
     ThetaP(IP,ilmP,2,3) = ThetaP2(IP,ilmP,6) + Ycd*ThetaP2(iP,ilmP,2) 
     ThetaP(IP,ilmP,3,3) = ThetaP2(IP,ilmP,8) + Ycd*ThetaP2(iP,ilmP,3) 
     ThetaP(IP,ilmP,4,3) = ThetaP2(IP,ilmP,9) + Ycd*ThetaP2(iP,ilmP,4) 
     ThetaP(IP,ilmP,2,4) = ThetaP2(IP,ilmP,7) + Zcd*ThetaP2(iP,ilmP,2) 
     ThetaP(IP,ilmP,3,4) = ThetaP2(IP,ilmP,9) + Zcd*ThetaP2(iP,ilmP,3) 
     ThetaP(IP,ilmP,4,4) = ThetaP2(IP,ilmP,10) + Zcd*ThetaP2(iP,ilmP,4) 
    ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_RHS_Q2C1D1CtoD

!Transfer angmom from C to D
subroutine HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri,iASync)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,   10)
  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,    5:   10,    1:    1)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  integer :: iP,ilmP,iTUVC
!  real(realk) :: Tmp(nTUVA,nTUVB) ordering
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVC,ilmP) &
!$ACC PRESENT(nPasses,ThetaP,ThetaP2) ASYNC(iASync)
  DO iP = 1,nContPQ*nPasses
    DO iTUVC=  5, 10
     DO ilmP = 1,nlmP
        ThetaP(IP,ilmP,iTUVC,1) = ThetaP2(IP,ilmP,iTUVC)
     ENDDO
    ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_RHS_Q2C2D0CtoD

!Transfer angmom from C to D
subroutine HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri,iASync)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,   20)
  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,    5:   10,    2:    4)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVC
  real(realk) :: Xcd,Ycd,Zcd
!  real(realk) :: Tmp(nTUVA,nTUVB) ordering
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,&
!$ACC         iTUVC,ilmP,Xcd,Ycd,Zcd) &
!$ACC PRESENT(nPasses,Qdistance12,ThetaP,ThetaP2) ASYNC(iASync)
  DO iP = 1,nContPQ*nPasses
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
    DO ilmP = 1,nlmP
     ThetaP(IP,ilmP,5,2) = ThetaP2(IP,ilmP,11) + Xcd*ThetaP2(iP,ilmP,5) 
     ThetaP(IP,ilmP,6,2) = ThetaP2(IP,ilmP,12) + Xcd*ThetaP2(iP,ilmP,6) 
     ThetaP(IP,ilmP,7,2) = ThetaP2(IP,ilmP,13) + Xcd*ThetaP2(iP,ilmP,7) 
     ThetaP(IP,ilmP,8,2) = ThetaP2(IP,ilmP,14) + Xcd*ThetaP2(iP,ilmP,8) 
     ThetaP(IP,ilmP,9,2) = ThetaP2(IP,ilmP,15) + Xcd*ThetaP2(iP,ilmP,9) 
     ThetaP(IP,ilmP,10,2) = ThetaP2(IP,ilmP,16) + Xcd*ThetaP2(iP,ilmP,10) 
     ThetaP(IP,ilmP,5,3) = ThetaP2(IP,ilmP,12) + Ycd*ThetaP2(iP,ilmP,5) 
     ThetaP(IP,ilmP,6,3) = ThetaP2(IP,ilmP,14) + Ycd*ThetaP2(iP,ilmP,6) 
     ThetaP(IP,ilmP,7,3) = ThetaP2(IP,ilmP,15) + Ycd*ThetaP2(iP,ilmP,7) 
     ThetaP(IP,ilmP,8,3) = ThetaP2(IP,ilmP,17) + Ycd*ThetaP2(iP,ilmP,8) 
     ThetaP(IP,ilmP,9,3) = ThetaP2(IP,ilmP,18) + Ycd*ThetaP2(iP,ilmP,9) 
     ThetaP(IP,ilmP,10,3) = ThetaP2(IP,ilmP,19) + Ycd*ThetaP2(iP,ilmP,10) 
     ThetaP(IP,ilmP,5,4) = ThetaP2(IP,ilmP,13) + Zcd*ThetaP2(iP,ilmP,5) 
     ThetaP(IP,ilmP,6,4) = ThetaP2(IP,ilmP,15) + Zcd*ThetaP2(iP,ilmP,6) 
     ThetaP(IP,ilmP,7,4) = ThetaP2(IP,ilmP,16) + Zcd*ThetaP2(iP,ilmP,7) 
     ThetaP(IP,ilmP,8,4) = ThetaP2(IP,ilmP,18) + Zcd*ThetaP2(iP,ilmP,8) 
     ThetaP(IP,ilmP,9,4) = ThetaP2(IP,ilmP,19) + Zcd*ThetaP2(iP,ilmP,9) 
     ThetaP(IP,ilmP,10,4) = ThetaP2(IP,ilmP,20) + Zcd*ThetaP2(iP,ilmP,10) 
    ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_RHS_Q3C2D1CtoD

!Transfer angmom from C to D
subroutine HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri,iASync)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,   35)
  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,    5:   10,    5:   10)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVC
  real(realk) :: Xcd,Ycd,Zcd
  real(realk) :: Tmp1(  5: 20,  2:  4)
!  real(realk) :: Tmp(nTUVA,nTUVB) ordering
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,&
!$ACC         Tmp1,&
!$ACC         iTUVC,ilmP,Xcd,Ycd,Zcd) &
!$ACC PRESENT(nPasses,Qdistance12,ThetaP,ThetaP2) ASYNC(iASync)
  DO iP = 1,nContPQ*nPasses
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
    DO ilmP = 1,nlmP
     Tmp1(5,2) = ThetaP2(IP,ilmP,11) + Xcd*ThetaP2(iP,ilmP,5) 
     Tmp1(6,2) = ThetaP2(IP,ilmP,12) + Xcd*ThetaP2(iP,ilmP,6) 
     Tmp1(7,2) = ThetaP2(IP,ilmP,13) + Xcd*ThetaP2(iP,ilmP,7) 
     Tmp1(8,2) = ThetaP2(IP,ilmP,14) + Xcd*ThetaP2(iP,ilmP,8) 
     Tmp1(9,2) = ThetaP2(IP,ilmP,15) + Xcd*ThetaP2(iP,ilmP,9) 
     Tmp1(10,2) = ThetaP2(IP,ilmP,16) + Xcd*ThetaP2(iP,ilmP,10) 
     Tmp1(11,2) = ThetaP2(IP,ilmP,21) + Xcd*ThetaP2(iP,ilmP,11) 
     Tmp1(12,2) = ThetaP2(IP,ilmP,22) + Xcd*ThetaP2(iP,ilmP,12) 
     Tmp1(13,2) = ThetaP2(IP,ilmP,23) + Xcd*ThetaP2(iP,ilmP,13) 
     Tmp1(14,2) = ThetaP2(IP,ilmP,24) + Xcd*ThetaP2(iP,ilmP,14) 
     Tmp1(15,2) = ThetaP2(IP,ilmP,25) + Xcd*ThetaP2(iP,ilmP,15) 
     Tmp1(16,2) = ThetaP2(IP,ilmP,26) + Xcd*ThetaP2(iP,ilmP,16) 
     Tmp1(17,2) = ThetaP2(IP,ilmP,27) + Xcd*ThetaP2(iP,ilmP,17) 
     Tmp1(18,2) = ThetaP2(IP,ilmP,28) + Xcd*ThetaP2(iP,ilmP,18) 
     Tmp1(19,2) = ThetaP2(IP,ilmP,29) + Xcd*ThetaP2(iP,ilmP,19) 
     Tmp1(20,2) = ThetaP2(IP,ilmP,30) + Xcd*ThetaP2(iP,ilmP,20) 
     Tmp1(5,3) = ThetaP2(IP,ilmP,12) + Ycd*ThetaP2(iP,ilmP,5) 
     Tmp1(6,3) = ThetaP2(IP,ilmP,14) + Ycd*ThetaP2(iP,ilmP,6) 
     Tmp1(7,3) = ThetaP2(IP,ilmP,15) + Ycd*ThetaP2(iP,ilmP,7) 
     Tmp1(8,3) = ThetaP2(IP,ilmP,17) + Ycd*ThetaP2(iP,ilmP,8) 
     Tmp1(9,3) = ThetaP2(IP,ilmP,18) + Ycd*ThetaP2(iP,ilmP,9) 
     Tmp1(10,3) = ThetaP2(IP,ilmP,19) + Ycd*ThetaP2(iP,ilmP,10) 
     Tmp1(11,3) = ThetaP2(IP,ilmP,22) + Ycd*ThetaP2(iP,ilmP,11) 
     Tmp1(12,3) = ThetaP2(IP,ilmP,24) + Ycd*ThetaP2(iP,ilmP,12) 
     Tmp1(13,3) = ThetaP2(IP,ilmP,25) + Ycd*ThetaP2(iP,ilmP,13) 
     Tmp1(14,3) = ThetaP2(IP,ilmP,27) + Ycd*ThetaP2(iP,ilmP,14) 
     Tmp1(15,3) = ThetaP2(IP,ilmP,28) + Ycd*ThetaP2(iP,ilmP,15) 
     Tmp1(16,3) = ThetaP2(IP,ilmP,29) + Ycd*ThetaP2(iP,ilmP,16) 
     Tmp1(17,3) = ThetaP2(IP,ilmP,31) + Ycd*ThetaP2(iP,ilmP,17) 
     Tmp1(18,3) = ThetaP2(IP,ilmP,32) + Ycd*ThetaP2(iP,ilmP,18) 
     Tmp1(19,3) = ThetaP2(IP,ilmP,33) + Ycd*ThetaP2(iP,ilmP,19) 
     Tmp1(20,3) = ThetaP2(IP,ilmP,34) + Ycd*ThetaP2(iP,ilmP,20) 
     Tmp1(5,4) = ThetaP2(IP,ilmP,13) + Zcd*ThetaP2(iP,ilmP,5) 
     Tmp1(6,4) = ThetaP2(IP,ilmP,15) + Zcd*ThetaP2(iP,ilmP,6) 
     Tmp1(7,4) = ThetaP2(IP,ilmP,16) + Zcd*ThetaP2(iP,ilmP,7) 
     Tmp1(8,4) = ThetaP2(IP,ilmP,18) + Zcd*ThetaP2(iP,ilmP,8) 
     Tmp1(9,4) = ThetaP2(IP,ilmP,19) + Zcd*ThetaP2(iP,ilmP,9) 
     Tmp1(10,4) = ThetaP2(IP,ilmP,20) + Zcd*ThetaP2(iP,ilmP,10) 
     Tmp1(11,4) = ThetaP2(IP,ilmP,23) + Zcd*ThetaP2(iP,ilmP,11) 
     Tmp1(12,4) = ThetaP2(IP,ilmP,25) + Zcd*ThetaP2(iP,ilmP,12) 
     Tmp1(13,4) = ThetaP2(IP,ilmP,26) + Zcd*ThetaP2(iP,ilmP,13) 
     Tmp1(14,4) = ThetaP2(IP,ilmP,28) + Zcd*ThetaP2(iP,ilmP,14) 
     Tmp1(15,4) = ThetaP2(IP,ilmP,29) + Zcd*ThetaP2(iP,ilmP,15) 
     Tmp1(16,4) = ThetaP2(IP,ilmP,30) + Zcd*ThetaP2(iP,ilmP,16) 
     Tmp1(17,4) = ThetaP2(IP,ilmP,32) + Zcd*ThetaP2(iP,ilmP,17) 
     Tmp1(18,4) = ThetaP2(IP,ilmP,33) + Zcd*ThetaP2(iP,ilmP,18) 
     Tmp1(19,4) = ThetaP2(IP,ilmP,34) + Zcd*ThetaP2(iP,ilmP,19) 
     Tmp1(20,4) = ThetaP2(IP,ilmP,35) + Zcd*ThetaP2(iP,ilmP,20) 
     ThetaP(IP,ilmP,5,5) = Tmp1(11,2) + Xcd*Tmp1(5,2) 
     ThetaP(IP,ilmP,6,5) = Tmp1(12,2) + Xcd*Tmp1(6,2) 
     ThetaP(IP,ilmP,7,5) = Tmp1(13,2) + Xcd*Tmp1(7,2) 
     ThetaP(IP,ilmP,8,5) = Tmp1(14,2) + Xcd*Tmp1(8,2) 
     ThetaP(IP,ilmP,9,5) = Tmp1(15,2) + Xcd*Tmp1(9,2) 
     ThetaP(IP,ilmP,10,5) = Tmp1(16,2) + Xcd*Tmp1(10,2) 
     ThetaP(IP,ilmP,5,6) = Tmp1(11,3) + Xcd*Tmp1(5,3) 
     ThetaP(IP,ilmP,6,6) = Tmp1(12,3) + Xcd*Tmp1(6,3) 
     ThetaP(IP,ilmP,7,6) = Tmp1(13,3) + Xcd*Tmp1(7,3) 
     ThetaP(IP,ilmP,8,6) = Tmp1(14,3) + Xcd*Tmp1(8,3) 
     ThetaP(IP,ilmP,9,6) = Tmp1(15,3) + Xcd*Tmp1(9,3) 
     ThetaP(IP,ilmP,10,6) = Tmp1(16,3) + Xcd*Tmp1(10,3) 
     ThetaP(IP,ilmP,5,7) = Tmp1(11,4) + Xcd*Tmp1(5,4) 
     ThetaP(IP,ilmP,6,7) = Tmp1(12,4) + Xcd*Tmp1(6,4) 
     ThetaP(IP,ilmP,7,7) = Tmp1(13,4) + Xcd*Tmp1(7,4) 
     ThetaP(IP,ilmP,8,7) = Tmp1(14,4) + Xcd*Tmp1(8,4) 
     ThetaP(IP,ilmP,9,7) = Tmp1(15,4) + Xcd*Tmp1(9,4) 
     ThetaP(IP,ilmP,10,7) = Tmp1(16,4) + Xcd*Tmp1(10,4) 
     ThetaP(IP,ilmP,5,8) = Tmp1(12,3) + Ycd*Tmp1(5,3) 
     ThetaP(IP,ilmP,6,8) = Tmp1(14,3) + Ycd*Tmp1(6,3) 
     ThetaP(IP,ilmP,7,8) = Tmp1(15,3) + Ycd*Tmp1(7,3) 
     ThetaP(IP,ilmP,8,8) = Tmp1(17,3) + Ycd*Tmp1(8,3) 
     ThetaP(IP,ilmP,9,8) = Tmp1(18,3) + Ycd*Tmp1(9,3) 
     ThetaP(IP,ilmP,10,8) = Tmp1(19,3) + Ycd*Tmp1(10,3) 
     ThetaP(IP,ilmP,5,9) = Tmp1(12,4) + Ycd*Tmp1(5,4) 
     ThetaP(IP,ilmP,6,9) = Tmp1(14,4) + Ycd*Tmp1(6,4) 
     ThetaP(IP,ilmP,7,9) = Tmp1(15,4) + Ycd*Tmp1(7,4) 
     ThetaP(IP,ilmP,8,9) = Tmp1(17,4) + Ycd*Tmp1(8,4) 
     ThetaP(IP,ilmP,9,9) = Tmp1(18,4) + Ycd*Tmp1(9,4) 
     ThetaP(IP,ilmP,10,9) = Tmp1(19,4) + Ycd*Tmp1(10,4) 
     ThetaP(IP,ilmP,5,10) = Tmp1(13,4) + Zcd*Tmp1(5,4) 
     ThetaP(IP,ilmP,6,10) = Tmp1(15,4) + Zcd*Tmp1(6,4) 
     ThetaP(IP,ilmP,7,10) = Tmp1(16,4) + Zcd*Tmp1(7,4) 
     ThetaP(IP,ilmP,8,10) = Tmp1(18,4) + Zcd*Tmp1(8,4) 
     ThetaP(IP,ilmP,9,10) = Tmp1(19,4) + Zcd*Tmp1(9,4) 
     ThetaP(IP,ilmP,10,10) = Tmp1(20,4) + Zcd*Tmp1(10,4) 
    ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_RHS_Q4C2D2CtoD
end module
