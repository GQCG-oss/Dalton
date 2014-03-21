MODULE AGC_OBS_HorizontalRecurrenceLHSModBtoA
 use IchorPrecisionModule
  
 CONTAINS

subroutine HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(    4,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    1:    1,    2:    4,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iTUVB) 
  DO iP = 1,nTUVQ*nContQP*nPasses
     DO iTUVB=  2,  4
        ThetaP(1,iTUVB,iP) = AuxCont(iTUVB,iP)
     ENDDO
  ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_LHS_P1A0B1BtoA

subroutine HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(   10,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    1:    1,    5:   10,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iTUVB) 
  DO iP = 1,nTUVQ*nContQP*nPasses
     DO iTUVB=  5, 10
        ThetaP(1,iTUVB,iP) = AuxCont(iTUVB,iP)
     ENDDO
  ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_LHS_P2A0B2BtoA

subroutine HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(   20,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    2:    4,    5:   10,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iPassP,iTUVB,iAtomA,iAtomB,Xab,Yab,Zab) 
  DO iP = 1,nTUVQ*nContQP*nPasses
   iPassP = (iP-1)/(nTUVQ*nContQP)+1
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = -Pdistance12(1,iAtomA,iAtomB)
   Yab = -Pdistance12(2,iAtomA,iAtomB)
   Zab = -Pdistance12(3,iAtomA,iAtomB)
     ThetaP( 2, 5      ,iP) = AuxCont(11      ,iP)+ Xab*AuxCont( 5,      ip) 
     ThetaP( 2, 6      ,iP) = AuxCont(12      ,iP)+ Xab*AuxCont( 6,      ip) 
     ThetaP( 2, 7      ,iP) = AuxCont(13      ,iP)+ Xab*AuxCont( 7,      ip) 
     ThetaP( 2, 8      ,iP) = AuxCont(14      ,iP)+ Xab*AuxCont( 8,      ip) 
     ThetaP( 2, 9      ,iP) = AuxCont(15      ,iP)+ Xab*AuxCont( 9,      ip) 
     ThetaP( 2,10      ,iP) = AuxCont(16      ,iP)+ Xab*AuxCont(10,      ip) 
     ThetaP( 3, 5      ,iP) = AuxCont(12      ,iP)+ Yab*AuxCont( 5,      ip) 
     ThetaP( 3, 6      ,iP) = AuxCont(14      ,iP)+ Yab*AuxCont( 6,      ip) 
     ThetaP( 3, 7      ,iP) = AuxCont(15      ,iP)+ Yab*AuxCont( 7,      ip) 
     ThetaP( 3, 8      ,iP) = AuxCont(17      ,iP)+ Yab*AuxCont( 8,      ip) 
     ThetaP( 3, 9      ,iP) = AuxCont(18      ,iP)+ Yab*AuxCont( 9,      ip) 
     ThetaP( 3,10      ,iP) = AuxCont(19      ,iP)+ Yab*AuxCont(10,      ip) 
     ThetaP( 4, 5      ,iP) = AuxCont(13      ,iP)+ Zab*AuxCont( 5,      ip) 
     ThetaP( 4, 6      ,iP) = AuxCont(15      ,iP)+ Zab*AuxCont( 6,      ip) 
     ThetaP( 4, 7      ,iP) = AuxCont(16      ,iP)+ Zab*AuxCont( 7,      ip) 
     ThetaP( 4, 8      ,iP) = AuxCont(18      ,iP)+ Zab*AuxCont( 8,      ip) 
     ThetaP( 4, 9      ,iP) = AuxCont(19      ,iP)+ Zab*AuxCont( 9,      ip) 
     ThetaP( 4,10      ,iP) = AuxCont(20      ,iP)+ Zab*AuxCont(10,      ip) 
  ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_LHS_P3A1B2BtoA
#ifdef VAR_OPENACC

subroutine HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(    4,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    1:    1,    2:    4,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,&
!$ACC         iTUVB) &
!$ACC PRESENT(nTUVQ,nContQP,nPasses,&
!$ACC         AuxCont,ThetaP)
  DO iP = 1,nTUVQ*nContQP*nPasses
     DO iTUVB=  2,  4
        ThetaP(1,iTUVB,iP) = AuxCont(iTUVB,iP)
     ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_LHS_P1A0B1BtoA

subroutine HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(   10,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    1:    1,    5:   10,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,&
!$ACC         iTUVB) &
!$ACC PRESENT(nTUVQ,nContQP,nPasses,&
!$ACC         AuxCont,ThetaP)
  DO iP = 1,nTUVQ*nContQP*nPasses
     DO iTUVB=  5, 10
        ThetaP(1,iTUVB,iP) = AuxCont(iTUVB,iP)
     ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_LHS_P2A0B2BtoA

subroutine HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(   20,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    2:    4,    5:   10,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,&
!$ACC         iPassP,iTUVB,iAtomA,iAtomB,Xab,Yab,Zab) &
!$ACC PRESENT(nTUVQ,nContQP,nPasses,&
!$ACC         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP)
  DO iP = 1,nTUVQ*nContQP*nPasses
   iPassP = (iP-1)/(nTUVQ*nContQP)+1
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = -Pdistance12(1,iAtomA,iAtomB)
   Yab = -Pdistance12(2,iAtomA,iAtomB)
   Zab = -Pdistance12(3,iAtomA,iAtomB)
     ThetaP( 2, 5      ,iP) = AuxCont(11      ,iP)+ Xab*AuxCont( 5,      ip) 
     ThetaP( 2, 6      ,iP) = AuxCont(12      ,iP)+ Xab*AuxCont( 6,      ip) 
     ThetaP( 2, 7      ,iP) = AuxCont(13      ,iP)+ Xab*AuxCont( 7,      ip) 
     ThetaP( 2, 8      ,iP) = AuxCont(14      ,iP)+ Xab*AuxCont( 8,      ip) 
     ThetaP( 2, 9      ,iP) = AuxCont(15      ,iP)+ Xab*AuxCont( 9,      ip) 
     ThetaP( 2,10      ,iP) = AuxCont(16      ,iP)+ Xab*AuxCont(10,      ip) 
     ThetaP( 3, 5      ,iP) = AuxCont(12      ,iP)+ Yab*AuxCont( 5,      ip) 
     ThetaP( 3, 6      ,iP) = AuxCont(14      ,iP)+ Yab*AuxCont( 6,      ip) 
     ThetaP( 3, 7      ,iP) = AuxCont(15      ,iP)+ Yab*AuxCont( 7,      ip) 
     ThetaP( 3, 8      ,iP) = AuxCont(17      ,iP)+ Yab*AuxCont( 8,      ip) 
     ThetaP( 3, 9      ,iP) = AuxCont(18      ,iP)+ Yab*AuxCont( 9,      ip) 
     ThetaP( 3,10      ,iP) = AuxCont(19      ,iP)+ Yab*AuxCont(10,      ip) 
     ThetaP( 4, 5      ,iP) = AuxCont(13      ,iP)+ Zab*AuxCont( 5,      ip) 
     ThetaP( 4, 6      ,iP) = AuxCont(15      ,iP)+ Zab*AuxCont( 6,      ip) 
     ThetaP( 4, 7      ,iP) = AuxCont(16      ,iP)+ Zab*AuxCont( 7,      ip) 
     ThetaP( 4, 8      ,iP) = AuxCont(18      ,iP)+ Zab*AuxCont( 8,      ip) 
     ThetaP( 4, 9      ,iP) = AuxCont(19      ,iP)+ Zab*AuxCont( 9,      ip) 
     ThetaP( 4,10      ,iP) = AuxCont(20      ,iP)+ Zab*AuxCont(10,      ip) 
  ENDDO
end subroutine HorizontalRR_GPU_LHS_P3A1B2BtoA
#endif
end module
