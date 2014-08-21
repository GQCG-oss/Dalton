MODULE AGC_GPU_OBS_HorizontalRecurrenceRHSModDtoC
 use IchorPrecisionModule
  
 CONTAINS

!Transfer angmom from D to C
subroutine HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,    4)
  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,1,    2:    4)
  !Local variables
  integer :: iP,ilmP,iTUVD
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVD,ilmP) &
!$ACC PRESENT(nContPQ,nPasses,ThetaP,ThetaP2)
  DO iP = 1,nContPQ*nPasses
    DO iTUVD=  2,  4
     DO ilmP = 1,nlmP
        ThetaP(IP,ilmP,1,iTUVD) = ThetaP2(IP,ilmP,iTUVD)
     ENDDO
    ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_RHS_Q1C0D1DtoC

!Transfer angmom from D to C
subroutine HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,   10)
  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,1,    5:   10)
  !Local variables
  integer :: iP,ilmP,iTUVD
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVD,ilmP) &
!$ACC PRESENT(nContPQ,nPasses,ThetaP,ThetaP2)
  DO iP = 1,nContPQ*nPasses
    DO iTUVD=  5, 10
     DO ilmP = 1,nlmP
        ThetaP(IP,ilmP,1,iTUVD) = ThetaP2(IP,ilmP,iTUVD)
     ENDDO
    ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_RHS_Q2C0D2DtoC

!Transfer angmom from D to C
subroutine HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,   20)
  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,    2:    4,    5:   10)
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVD
  real(realk) :: Xcd,Ycd,Zcd
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,&
!$ACC         iTUVD,ilmP,Xcd,Ycd,Zcd) &
!$ACC PRESENT(nContPQ,nPasses,Qdistance12,ThetaP,ThetaP2)
  DO iP = 1,nContPQ*nPasses
   Xcd = -Qdistance12(1)
   Ycd = -Qdistance12(2)
   Zcd = -Qdistance12(3)
    DO ilmP = 1,nlmP
     ThetaP(iP,ilmP,2,5) = ThetaP2(iP,ilmP,11) + Xcd*ThetaP2(iP,ilmP,5) 
     ThetaP(iP,ilmP,2,6) = ThetaP2(iP,ilmP,12) + Xcd*ThetaP2(iP,ilmP,6) 
     ThetaP(iP,ilmP,2,7) = ThetaP2(iP,ilmP,13) + Xcd*ThetaP2(iP,ilmP,7) 
     ThetaP(iP,ilmP,2,8) = ThetaP2(iP,ilmP,14) + Xcd*ThetaP2(iP,ilmP,8) 
     ThetaP(iP,ilmP,2,9) = ThetaP2(iP,ilmP,15) + Xcd*ThetaP2(iP,ilmP,9) 
     ThetaP(iP,ilmP,2,10) = ThetaP2(iP,ilmP,16) + Xcd*ThetaP2(iP,ilmP,10) 
     ThetaP(iP,ilmP,3,5) = ThetaP2(iP,ilmP,12) + Ycd*ThetaP2(iP,ilmP,5) 
     ThetaP(iP,ilmP,3,6) = ThetaP2(iP,ilmP,14) + Ycd*ThetaP2(iP,ilmP,6) 
     ThetaP(iP,ilmP,3,7) = ThetaP2(iP,ilmP,15) + Ycd*ThetaP2(iP,ilmP,7) 
     ThetaP(iP,ilmP,3,8) = ThetaP2(iP,ilmP,17) + Ycd*ThetaP2(iP,ilmP,8) 
     ThetaP(iP,ilmP,3,9) = ThetaP2(iP,ilmP,18) + Ycd*ThetaP2(iP,ilmP,9) 
     ThetaP(iP,ilmP,3,10) = ThetaP2(iP,ilmP,19) + Ycd*ThetaP2(iP,ilmP,10) 
     ThetaP(iP,ilmP,4,5) = ThetaP2(iP,ilmP,13) + Zcd*ThetaP2(iP,ilmP,5) 
     ThetaP(iP,ilmP,4,6) = ThetaP2(iP,ilmP,15) + Zcd*ThetaP2(iP,ilmP,6) 
     ThetaP(iP,ilmP,4,7) = ThetaP2(iP,ilmP,16) + Zcd*ThetaP2(iP,ilmP,7) 
     ThetaP(iP,ilmP,4,8) = ThetaP2(iP,ilmP,18) + Zcd*ThetaP2(iP,ilmP,8) 
     ThetaP(iP,ilmP,4,9) = ThetaP2(iP,ilmP,19) + Zcd*ThetaP2(iP,ilmP,9) 
     ThetaP(iP,ilmP,4,10) = ThetaP2(iP,ilmP,20) + Zcd*ThetaP2(iP,ilmP,10) 
    ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_RHS_Q3C1D2DtoC
end module
