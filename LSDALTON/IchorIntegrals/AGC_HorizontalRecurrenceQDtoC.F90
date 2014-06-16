MODULE AGC_OBS_HorizontalRecurrenceRHSModDtoC
 use IchorPrecisionModule
  
 CONTAINS

!Transfer angmom from D to C
subroutine HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nlmP,    4,nContPQ*nPasses)
  real(realk),intent(inout) :: ThetaP(nlmP,1,    2:    4,nContPQ*nPasses)
  !Local variables
  integer :: iP,ilmP,iTUVD
!$OMP DO&
!$OMP PRIVATE(iP,iTUVD,ilmP) 
  DO iP = 1,nContPQ*nPasses
    DO iTUVD=  2,  4
     DO ilmP = 1,nlmP
        ThetaP(ilmP,1,iTUVD,IP) = ThetaP2(ilmP,iTUVD,IP)
     ENDDO
    ENDDO
   ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_RHS_Q1C0D1DtoC

!Transfer angmom from D to C
subroutine HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nlmP,   10,nContPQ*nPasses)
  real(realk),intent(inout) :: ThetaP(nlmP,1,    5:   10,nContPQ*nPasses)
  !Local variables
  integer :: iP,ilmP,iTUVD
!$OMP DO&
!$OMP PRIVATE(iP,iTUVD,ilmP) 
  DO iP = 1,nContPQ*nPasses
    DO iTUVD=  5, 10
     DO ilmP = 1,nlmP
        ThetaP(ilmP,1,iTUVD,IP) = ThetaP2(ilmP,iTUVD,IP)
     ENDDO
    ENDDO
   ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_RHS_Q2C0D2DtoC

!Transfer angmom from D to C
subroutine HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nlmP,   20,nContPQ*nPasses)
  real(realk),intent(inout) :: ThetaP(nlmP,    2:    4,    5:   10,nContPQ*nPasses)
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVD
  real(realk) :: Xcd,Ycd,Zcd
!$OMP DO&
!$OMP PRIVATE(iP,&
!$OMP         iTUVD,ilmP,Xcd,Ycd,Zcd) 
  DO iP = 1,nContPQ*nPasses
   Xcd = -Qdistance12(1)
   Ycd = -Qdistance12(2)
   Zcd = -Qdistance12(3)
    DO ilmP = 1,nlmP
     ThetaP(ilmP, 2, 5,IP) = ThetaP2(ilmP,11,IP) + Xcd*ThetaP2(ilmP, 5, IP) 
     ThetaP(ilmP, 2, 6,IP) = ThetaP2(ilmP,12,IP) + Xcd*ThetaP2(ilmP, 6, IP) 
     ThetaP(ilmP, 2, 7,IP) = ThetaP2(ilmP,13,IP) + Xcd*ThetaP2(ilmP, 7, IP) 
     ThetaP(ilmP, 2, 8,IP) = ThetaP2(ilmP,14,IP) + Xcd*ThetaP2(ilmP, 8, IP) 
     ThetaP(ilmP, 2, 9,IP) = ThetaP2(ilmP,15,IP) + Xcd*ThetaP2(ilmP, 9, IP) 
     ThetaP(ilmP, 2,10,IP) = ThetaP2(ilmP,16,IP) + Xcd*ThetaP2(ilmP,10, IP) 
     ThetaP(ilmP, 3, 5,IP) = ThetaP2(ilmP,12,IP) + Ycd*ThetaP2(ilmP, 5, IP) 
     ThetaP(ilmP, 3, 6,IP) = ThetaP2(ilmP,14,IP) + Ycd*ThetaP2(ilmP, 6, IP) 
     ThetaP(ilmP, 3, 7,IP) = ThetaP2(ilmP,15,IP) + Ycd*ThetaP2(ilmP, 7, IP) 
     ThetaP(ilmP, 3, 8,IP) = ThetaP2(ilmP,17,IP) + Ycd*ThetaP2(ilmP, 8, IP) 
     ThetaP(ilmP, 3, 9,IP) = ThetaP2(ilmP,18,IP) + Ycd*ThetaP2(ilmP, 9, IP) 
     ThetaP(ilmP, 3,10,IP) = ThetaP2(ilmP,19,IP) + Ycd*ThetaP2(ilmP,10, IP) 
     ThetaP(ilmP, 4, 5,IP) = ThetaP2(ilmP,13,IP) + Zcd*ThetaP2(ilmP, 5, IP) 
     ThetaP(ilmP, 4, 6,IP) = ThetaP2(ilmP,15,IP) + Zcd*ThetaP2(ilmP, 6, IP) 
     ThetaP(ilmP, 4, 7,IP) = ThetaP2(ilmP,16,IP) + Zcd*ThetaP2(ilmP, 7, IP) 
     ThetaP(ilmP, 4, 8,IP) = ThetaP2(ilmP,18,IP) + Zcd*ThetaP2(ilmP, 8, IP) 
     ThetaP(ilmP, 4, 9,IP) = ThetaP2(ilmP,19,IP) + Zcd*ThetaP2(ilmP, 9, IP) 
     ThetaP(ilmP, 4,10,IP) = ThetaP2(ilmP,20,IP) + Zcd*ThetaP2(ilmP,10, IP) 
    ENDDO
   ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_RHS_Q3C1D2DtoC
#ifdef VAR_OPENACC

!Transfer angmom from D to C
subroutine HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(realk),intent(in) :: Qdistance12(3)
  real(realk),intent(in) :: ThetaP2(nlmP,    4,nContPQ*nPasses)
  real(realk),intent(inout) :: ThetaP(nlmP,1,    2:    4,nContPQ*nPasses)
  !Local variables
  integer :: iP,ilmP,iTUVD
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVD,ilmP) &
!$ACC PRESENT(nlmP,nContPQ,nPasses,ThetaP,ThetaP2)
  DO iP = 1,nContPQ*nPasses
    DO iTUVD=  2,  4
     DO ilmP = 1,nlmP
        ThetaP(ilmP,1,iTUVD,IP) = ThetaP2(ilmP,iTUVD,IP)
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
  real(realk),intent(in) :: ThetaP2(nlmP,   10,nContPQ*nPasses)
  real(realk),intent(inout) :: ThetaP(nlmP,1,    5:   10,nContPQ*nPasses)
  !Local variables
  integer :: iP,ilmP,iTUVD
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVD,ilmP) &
!$ACC PRESENT(nlmP,nContPQ,nPasses,ThetaP,ThetaP2)
  DO iP = 1,nContPQ*nPasses
    DO iTUVD=  5, 10
     DO ilmP = 1,nlmP
        ThetaP(ilmP,1,iTUVD,IP) = ThetaP2(ilmP,iTUVD,IP)
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
  real(realk),intent(in) :: ThetaP2(nlmP,   20,nContPQ*nPasses)
  real(realk),intent(inout) :: ThetaP(nlmP,    2:    4,    5:   10,nContPQ*nPasses)
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVD
  real(realk) :: Xcd,Ycd,Zcd
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,&
!$ACC         iTUVD,ilmP,Xcd,Ycd,Zcd) &
!$ACC PRESENT(nlmP,nContPQ,nPasses,Qdistance12,ThetaP,ThetaP2)
  DO iP = 1,nContPQ*nPasses
   Xcd = -Qdistance12(1)
   Ycd = -Qdistance12(2)
   Zcd = -Qdistance12(3)
    DO ilmP = 1,nlmP
     ThetaP(ilmP, 2, 5,IP) = ThetaP2(ilmP,11,IP) + Xcd*ThetaP2(ilmP, 5, IP) 
     ThetaP(ilmP, 2, 6,IP) = ThetaP2(ilmP,12,IP) + Xcd*ThetaP2(ilmP, 6, IP) 
     ThetaP(ilmP, 2, 7,IP) = ThetaP2(ilmP,13,IP) + Xcd*ThetaP2(ilmP, 7, IP) 
     ThetaP(ilmP, 2, 8,IP) = ThetaP2(ilmP,14,IP) + Xcd*ThetaP2(ilmP, 8, IP) 
     ThetaP(ilmP, 2, 9,IP) = ThetaP2(ilmP,15,IP) + Xcd*ThetaP2(ilmP, 9, IP) 
     ThetaP(ilmP, 2,10,IP) = ThetaP2(ilmP,16,IP) + Xcd*ThetaP2(ilmP,10, IP) 
     ThetaP(ilmP, 3, 5,IP) = ThetaP2(ilmP,12,IP) + Ycd*ThetaP2(ilmP, 5, IP) 
     ThetaP(ilmP, 3, 6,IP) = ThetaP2(ilmP,14,IP) + Ycd*ThetaP2(ilmP, 6, IP) 
     ThetaP(ilmP, 3, 7,IP) = ThetaP2(ilmP,15,IP) + Ycd*ThetaP2(ilmP, 7, IP) 
     ThetaP(ilmP, 3, 8,IP) = ThetaP2(ilmP,17,IP) + Ycd*ThetaP2(ilmP, 8, IP) 
     ThetaP(ilmP, 3, 9,IP) = ThetaP2(ilmP,18,IP) + Ycd*ThetaP2(ilmP, 9, IP) 
     ThetaP(ilmP, 3,10,IP) = ThetaP2(ilmP,19,IP) + Ycd*ThetaP2(ilmP,10, IP) 
     ThetaP(ilmP, 4, 5,IP) = ThetaP2(ilmP,13,IP) + Zcd*ThetaP2(ilmP, 5, IP) 
     ThetaP(ilmP, 4, 6,IP) = ThetaP2(ilmP,15,IP) + Zcd*ThetaP2(ilmP, 6, IP) 
     ThetaP(ilmP, 4, 7,IP) = ThetaP2(ilmP,16,IP) + Zcd*ThetaP2(ilmP, 7, IP) 
     ThetaP(ilmP, 4, 8,IP) = ThetaP2(ilmP,18,IP) + Zcd*ThetaP2(ilmP, 8, IP) 
     ThetaP(ilmP, 4, 9,IP) = ThetaP2(ilmP,19,IP) + Zcd*ThetaP2(ilmP, 9, IP) 
     ThetaP(ilmP, 4,10,IP) = ThetaP2(ilmP,20,IP) + Zcd*ThetaP2(ilmP,10, IP) 
    ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_RHS_Q3C1D2DtoC
#endif
end module
