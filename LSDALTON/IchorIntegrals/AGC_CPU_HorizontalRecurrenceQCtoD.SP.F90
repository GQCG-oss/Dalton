module SPAGC_CPU_OBS_HorizontalRecurrenceRHSModCtoD
 use IchorPrecisionMod
  
 CONTAINS

!Transfer angmom from C to D
subroutine SPHorizontalRR_CPU_RHS_Q1C1D0CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(reals),intent(in) :: Qdistance12(3)
  real(reals),intent(in) :: ThetaP2(nlmP,    4,nContPQ*nPasses)
  real(reals),intent(inout) :: ThetaP(nlmP,    2:    4,    1:    1,nContPQ*nPasses)
  !Local variables
  integer :: iP,ilmP,iTUVC
!  real(reals) :: Tmp(nTUVA,nTUVB) ordering
!$OMP DO &
!$OMP PRIVATE(iP,iTUVC,ilmP) 
  DO iP = 1,nContPQ*nPasses
    DO iTUVC=  2,  4
     DO ilmP = 1,nlmP
        ThetaP(ilmP,iTUVC,1,IP) = ThetaP2(ilmP,iTUVC,IP)
     ENDDO
    ENDDO
  ENDDO
!$OMP END DO
end subroutine SPHorizontalRR_CPU_RHS_Q1C1D0CtoD

!Transfer angmom from C to D
subroutine SPHorizontalRR_CPU_RHS_Q2C1D1CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(reals),intent(in) :: Qdistance12(3)
  real(reals),intent(in) :: ThetaP2(nlmP,   10,nContPQ*nPasses)
  real(reals),intent(inout) :: ThetaP(nlmP,    2:    4,    2:    4,nContPQ*nPasses)
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVC
  real(reals) :: Xcd,Ycd,Zcd
!  real(reals) :: Tmp(nTUVA,nTUVB) ordering
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iTUVC,ilmP,Xcd,Ycd,Zcd) 
  DO iP = 1,nContPQ*nPasses
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
    DO ilmP = 1,nlmP
     ThetaP(ilmP,2,2,IP) = ThetaP2(ilmP,5,IP) + Xcd*ThetaP2(ilmP,2,IP) 
     ThetaP(ilmP,3,2,IP) = ThetaP2(ilmP,6,IP) + Xcd*ThetaP2(ilmP,3,IP) 
     ThetaP(ilmP,4,2,IP) = ThetaP2(ilmP,7,IP) + Xcd*ThetaP2(ilmP,4,IP) 
     ThetaP(ilmP,2,3,IP) = ThetaP2(ilmP,6,IP) + Ycd*ThetaP2(ilmP,2,IP) 
     ThetaP(ilmP,3,3,IP) = ThetaP2(ilmP,8,IP) + Ycd*ThetaP2(ilmP,3,IP) 
     ThetaP(ilmP,4,3,IP) = ThetaP2(ilmP,9,IP) + Ycd*ThetaP2(ilmP,4,IP) 
     ThetaP(ilmP,2,4,IP) = ThetaP2(ilmP,7,IP) + Zcd*ThetaP2(ilmP,2,IP) 
     ThetaP(ilmP,3,4,IP) = ThetaP2(ilmP,9,IP) + Zcd*ThetaP2(ilmP,3,IP) 
     ThetaP(ilmP,4,4,IP) = ThetaP2(ilmP,10,IP) + Zcd*ThetaP2(ilmP,4,IP) 
    ENDDO
  ENDDO
!$OMP END DO
end subroutine SPHorizontalRR_CPU_RHS_Q2C1D1CtoD

!Transfer angmom from C to D
subroutine SPHorizontalRR_CPU_RHS_Q2C2D0CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(reals),intent(in) :: Qdistance12(3)
  real(reals),intent(in) :: ThetaP2(nlmP,   10,nContPQ*nPasses)
  real(reals),intent(inout) :: ThetaP(nlmP,    5:   10,    1:    1,nContPQ*nPasses)
  !Local variables
  integer :: iP,ilmP,iTUVC
!  real(reals) :: Tmp(nTUVA,nTUVB) ordering
!$OMP DO &
!$OMP PRIVATE(iP,iTUVC,ilmP) 
  DO iP = 1,nContPQ*nPasses
    DO iTUVC=  5, 10
     DO ilmP = 1,nlmP
        ThetaP(ilmP,iTUVC,1,IP) = ThetaP2(ilmP,iTUVC,IP)
     ENDDO
    ENDDO
  ENDDO
!$OMP END DO
end subroutine SPHorizontalRR_CPU_RHS_Q2C2D0CtoD

!Transfer angmom from C to D
subroutine SPHorizontalRR_CPU_RHS_Q3C2D1CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(reals),intent(in) :: Qdistance12(3)
  real(reals),intent(in) :: ThetaP2(nlmP,   20,nContPQ*nPasses)
  real(reals),intent(inout) :: ThetaP(nlmP,    5:   10,    2:    4,nContPQ*nPasses)
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVC
  real(reals) :: Xcd,Ycd,Zcd
!  real(reals) :: Tmp(nTUVA,nTUVB) ordering
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iTUVC,ilmP,Xcd,Ycd,Zcd) 
  DO iP = 1,nContPQ*nPasses
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
    DO ilmP = 1,nlmP
     ThetaP(ilmP,5,2,IP) = ThetaP2(ilmP,11,IP) + Xcd*ThetaP2(ilmP,5,IP) 
     ThetaP(ilmP,6,2,IP) = ThetaP2(ilmP,12,IP) + Xcd*ThetaP2(ilmP,6,IP) 
     ThetaP(ilmP,7,2,IP) = ThetaP2(ilmP,13,IP) + Xcd*ThetaP2(ilmP,7,IP) 
     ThetaP(ilmP,8,2,IP) = ThetaP2(ilmP,14,IP) + Xcd*ThetaP2(ilmP,8,IP) 
     ThetaP(ilmP,9,2,IP) = ThetaP2(ilmP,15,IP) + Xcd*ThetaP2(ilmP,9,IP) 
     ThetaP(ilmP,10,2,IP) = ThetaP2(ilmP,16,IP) + Xcd*ThetaP2(ilmP,10,IP) 
     ThetaP(ilmP,5,3,IP) = ThetaP2(ilmP,12,IP) + Ycd*ThetaP2(ilmP,5,IP) 
     ThetaP(ilmP,6,3,IP) = ThetaP2(ilmP,14,IP) + Ycd*ThetaP2(ilmP,6,IP) 
     ThetaP(ilmP,7,3,IP) = ThetaP2(ilmP,15,IP) + Ycd*ThetaP2(ilmP,7,IP) 
     ThetaP(ilmP,8,3,IP) = ThetaP2(ilmP,17,IP) + Ycd*ThetaP2(ilmP,8,IP) 
     ThetaP(ilmP,9,3,IP) = ThetaP2(ilmP,18,IP) + Ycd*ThetaP2(ilmP,9,IP) 
     ThetaP(ilmP,10,3,IP) = ThetaP2(ilmP,19,IP) + Ycd*ThetaP2(ilmP,10,IP) 
     ThetaP(ilmP,5,4,IP) = ThetaP2(ilmP,13,IP) + Zcd*ThetaP2(ilmP,5,IP) 
     ThetaP(ilmP,6,4,IP) = ThetaP2(ilmP,15,IP) + Zcd*ThetaP2(ilmP,6,IP) 
     ThetaP(ilmP,7,4,IP) = ThetaP2(ilmP,16,IP) + Zcd*ThetaP2(ilmP,7,IP) 
     ThetaP(ilmP,8,4,IP) = ThetaP2(ilmP,18,IP) + Zcd*ThetaP2(ilmP,8,IP) 
     ThetaP(ilmP,9,4,IP) = ThetaP2(ilmP,19,IP) + Zcd*ThetaP2(ilmP,9,IP) 
     ThetaP(ilmP,10,4,IP) = ThetaP2(ilmP,20,IP) + Zcd*ThetaP2(ilmP,10,IP) 
    ENDDO
  ENDDO
!$OMP END DO
end subroutine SPHorizontalRR_CPU_RHS_Q3C2D1CtoD

!Transfer angmom from C to D
subroutine SPHorizontalRR_CPU_RHS_Q4C2D2CtoD(nContPQ,nPasses,nlmP,&
         & Qdistance12,ThetaP2,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri
  real(reals),intent(in) :: Qdistance12(3)
  real(reals),intent(in) :: ThetaP2(nlmP,   35,nContPQ*nPasses)
  real(reals),intent(inout) :: ThetaP(nlmP,    5:   10,    5:   10,nContPQ*nPasses)
  !Local variables
  integer :: iP,iC,iPassQ,ilmP,iTUVC
  real(reals) :: Xcd,Ycd,Zcd
  real(reals) :: Tmp1(  5: 20,  2:  4)
!  real(reals) :: Tmp(nTUVA,nTUVB) ordering
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         Tmp1,&
!$OMP         iTUVC,ilmP,Xcd,Ycd,Zcd) 
  DO iP = 1,nContPQ*nPasses
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
    DO ilmP = 1,nlmP
     Tmp1(5,2) = ThetaP2(ilmP,11,IP) + Xcd*ThetaP2(ilmP,5,IP) 
     Tmp1(6,2) = ThetaP2(ilmP,12,IP) + Xcd*ThetaP2(ilmP,6,IP) 
     Tmp1(7,2) = ThetaP2(ilmP,13,IP) + Xcd*ThetaP2(ilmP,7,IP) 
     Tmp1(8,2) = ThetaP2(ilmP,14,IP) + Xcd*ThetaP2(ilmP,8,IP) 
     Tmp1(9,2) = ThetaP2(ilmP,15,IP) + Xcd*ThetaP2(ilmP,9,IP) 
     Tmp1(10,2) = ThetaP2(ilmP,16,IP) + Xcd*ThetaP2(ilmP,10,IP) 
     Tmp1(11,2) = ThetaP2(ilmP,21,IP) + Xcd*ThetaP2(ilmP,11,IP) 
     Tmp1(12,2) = ThetaP2(ilmP,22,IP) + Xcd*ThetaP2(ilmP,12,IP) 
     Tmp1(13,2) = ThetaP2(ilmP,23,IP) + Xcd*ThetaP2(ilmP,13,IP) 
     Tmp1(14,2) = ThetaP2(ilmP,24,IP) + Xcd*ThetaP2(ilmP,14,IP) 
     Tmp1(15,2) = ThetaP2(ilmP,25,IP) + Xcd*ThetaP2(ilmP,15,IP) 
     Tmp1(16,2) = ThetaP2(ilmP,26,IP) + Xcd*ThetaP2(ilmP,16,IP) 
     Tmp1(17,2) = ThetaP2(ilmP,27,IP) + Xcd*ThetaP2(ilmP,17,IP) 
     Tmp1(18,2) = ThetaP2(ilmP,28,IP) + Xcd*ThetaP2(ilmP,18,IP) 
     Tmp1(19,2) = ThetaP2(ilmP,29,IP) + Xcd*ThetaP2(ilmP,19,IP) 
     Tmp1(20,2) = ThetaP2(ilmP,30,IP) + Xcd*ThetaP2(ilmP,20,IP) 
     Tmp1(5,3) = ThetaP2(ilmP,12,IP) + Ycd*ThetaP2(ilmP,5,IP) 
     Tmp1(6,3) = ThetaP2(ilmP,14,IP) + Ycd*ThetaP2(ilmP,6,IP) 
     Tmp1(7,3) = ThetaP2(ilmP,15,IP) + Ycd*ThetaP2(ilmP,7,IP) 
     Tmp1(8,3) = ThetaP2(ilmP,17,IP) + Ycd*ThetaP2(ilmP,8,IP) 
     Tmp1(9,3) = ThetaP2(ilmP,18,IP) + Ycd*ThetaP2(ilmP,9,IP) 
     Tmp1(10,3) = ThetaP2(ilmP,19,IP) + Ycd*ThetaP2(ilmP,10,IP) 
     Tmp1(11,3) = ThetaP2(ilmP,22,IP) + Ycd*ThetaP2(ilmP,11,IP) 
     Tmp1(12,3) = ThetaP2(ilmP,24,IP) + Ycd*ThetaP2(ilmP,12,IP) 
     Tmp1(13,3) = ThetaP2(ilmP,25,IP) + Ycd*ThetaP2(ilmP,13,IP) 
     Tmp1(14,3) = ThetaP2(ilmP,27,IP) + Ycd*ThetaP2(ilmP,14,IP) 
     Tmp1(15,3) = ThetaP2(ilmP,28,IP) + Ycd*ThetaP2(ilmP,15,IP) 
     Tmp1(16,3) = ThetaP2(ilmP,29,IP) + Ycd*ThetaP2(ilmP,16,IP) 
     Tmp1(17,3) = ThetaP2(ilmP,31,IP) + Ycd*ThetaP2(ilmP,17,IP) 
     Tmp1(18,3) = ThetaP2(ilmP,32,IP) + Ycd*ThetaP2(ilmP,18,IP) 
     Tmp1(19,3) = ThetaP2(ilmP,33,IP) + Ycd*ThetaP2(ilmP,19,IP) 
     Tmp1(20,3) = ThetaP2(ilmP,34,IP) + Ycd*ThetaP2(ilmP,20,IP) 
     Tmp1(5,4) = ThetaP2(ilmP,13,IP) + Zcd*ThetaP2(ilmP,5,IP) 
     Tmp1(6,4) = ThetaP2(ilmP,15,IP) + Zcd*ThetaP2(ilmP,6,IP) 
     Tmp1(7,4) = ThetaP2(ilmP,16,IP) + Zcd*ThetaP2(ilmP,7,IP) 
     Tmp1(8,4) = ThetaP2(ilmP,18,IP) + Zcd*ThetaP2(ilmP,8,IP) 
     Tmp1(9,4) = ThetaP2(ilmP,19,IP) + Zcd*ThetaP2(ilmP,9,IP) 
     Tmp1(10,4) = ThetaP2(ilmP,20,IP) + Zcd*ThetaP2(ilmP,10,IP) 
     Tmp1(11,4) = ThetaP2(ilmP,23,IP) + Zcd*ThetaP2(ilmP,11,IP) 
     Tmp1(12,4) = ThetaP2(ilmP,25,IP) + Zcd*ThetaP2(ilmP,12,IP) 
     Tmp1(13,4) = ThetaP2(ilmP,26,IP) + Zcd*ThetaP2(ilmP,13,IP) 
     Tmp1(14,4) = ThetaP2(ilmP,28,IP) + Zcd*ThetaP2(ilmP,14,IP) 
     Tmp1(15,4) = ThetaP2(ilmP,29,IP) + Zcd*ThetaP2(ilmP,15,IP) 
     Tmp1(16,4) = ThetaP2(ilmP,30,IP) + Zcd*ThetaP2(ilmP,16,IP) 
     Tmp1(17,4) = ThetaP2(ilmP,32,IP) + Zcd*ThetaP2(ilmP,17,IP) 
     Tmp1(18,4) = ThetaP2(ilmP,33,IP) + Zcd*ThetaP2(ilmP,18,IP) 
     Tmp1(19,4) = ThetaP2(ilmP,34,IP) + Zcd*ThetaP2(ilmP,19,IP) 
     Tmp1(20,4) = ThetaP2(ilmP,35,IP) + Zcd*ThetaP2(ilmP,20,IP) 
     ThetaP(ilmP,5,5,IP) = Tmp1(11,2) + Xcd*Tmp1(5,2) 
     ThetaP(ilmP,6,5,IP) = Tmp1(12,2) + Xcd*Tmp1(6,2) 
     ThetaP(ilmP,7,5,IP) = Tmp1(13,2) + Xcd*Tmp1(7,2) 
     ThetaP(ilmP,8,5,IP) = Tmp1(14,2) + Xcd*Tmp1(8,2) 
     ThetaP(ilmP,9,5,IP) = Tmp1(15,2) + Xcd*Tmp1(9,2) 
     ThetaP(ilmP,10,5,IP) = Tmp1(16,2) + Xcd*Tmp1(10,2) 
     ThetaP(ilmP,5,6,IP) = Tmp1(11,3) + Xcd*Tmp1(5,3) 
     ThetaP(ilmP,6,6,IP) = Tmp1(12,3) + Xcd*Tmp1(6,3) 
     ThetaP(ilmP,7,6,IP) = Tmp1(13,3) + Xcd*Tmp1(7,3) 
     ThetaP(ilmP,8,6,IP) = Tmp1(14,3) + Xcd*Tmp1(8,3) 
     ThetaP(ilmP,9,6,IP) = Tmp1(15,3) + Xcd*Tmp1(9,3) 
     ThetaP(ilmP,10,6,IP) = Tmp1(16,3) + Xcd*Tmp1(10,3) 
     ThetaP(ilmP,5,7,IP) = Tmp1(11,4) + Xcd*Tmp1(5,4) 
     ThetaP(ilmP,6,7,IP) = Tmp1(12,4) + Xcd*Tmp1(6,4) 
     ThetaP(ilmP,7,7,IP) = Tmp1(13,4) + Xcd*Tmp1(7,4) 
     ThetaP(ilmP,8,7,IP) = Tmp1(14,4) + Xcd*Tmp1(8,4) 
     ThetaP(ilmP,9,7,IP) = Tmp1(15,4) + Xcd*Tmp1(9,4) 
     ThetaP(ilmP,10,7,IP) = Tmp1(16,4) + Xcd*Tmp1(10,4) 
     ThetaP(ilmP,5,8,IP) = Tmp1(12,3) + Ycd*Tmp1(5,3) 
     ThetaP(ilmP,6,8,IP) = Tmp1(14,3) + Ycd*Tmp1(6,3) 
     ThetaP(ilmP,7,8,IP) = Tmp1(15,3) + Ycd*Tmp1(7,3) 
     ThetaP(ilmP,8,8,IP) = Tmp1(17,3) + Ycd*Tmp1(8,3) 
     ThetaP(ilmP,9,8,IP) = Tmp1(18,3) + Ycd*Tmp1(9,3) 
     ThetaP(ilmP,10,8,IP) = Tmp1(19,3) + Ycd*Tmp1(10,3) 
     ThetaP(ilmP,5,9,IP) = Tmp1(12,4) + Ycd*Tmp1(5,4) 
     ThetaP(ilmP,6,9,IP) = Tmp1(14,4) + Ycd*Tmp1(6,4) 
     ThetaP(ilmP,7,9,IP) = Tmp1(15,4) + Ycd*Tmp1(7,4) 
     ThetaP(ilmP,8,9,IP) = Tmp1(17,4) + Ycd*Tmp1(8,4) 
     ThetaP(ilmP,9,9,IP) = Tmp1(18,4) + Ycd*Tmp1(9,4) 
     ThetaP(ilmP,10,9,IP) = Tmp1(19,4) + Ycd*Tmp1(10,4) 
     ThetaP(ilmP,5,10,IP) = Tmp1(13,4) + Zcd*Tmp1(5,4) 
     ThetaP(ilmP,6,10,IP) = Tmp1(15,4) + Zcd*Tmp1(6,4) 
     ThetaP(ilmP,7,10,IP) = Tmp1(16,4) + Zcd*Tmp1(7,4) 
     ThetaP(ilmP,8,10,IP) = Tmp1(18,4) + Zcd*Tmp1(8,4) 
     ThetaP(ilmP,9,10,IP) = Tmp1(19,4) + Zcd*Tmp1(9,4) 
     ThetaP(ilmP,10,10,IP) = Tmp1(20,4) + Zcd*Tmp1(10,4) 
    ENDDO
  ENDDO
!$OMP END DO
end subroutine SPHorizontalRR_CPU_RHS_Q4C2D2CtoD
end module
