MODULE AGC_GPU_PrimContractSegQMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
  
CONTAINS
  

   subroutine PrimitiveContractionAGPUSegQ(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nTUVP,nTUVQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB*nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nPrimB*nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iContA,iPrimA
    real(realk) :: TMP
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContA,iPrimA,TMP) &
!$ACC PRESENT(nPasses,nPrimA,nContA,nPrimB,&
!$ACC         ACC,AUXarrayCont,AUXarray2) ASYNC(iASync)
    do iP=1,nPasses*nTUVP*nTUVQ*nPrimB
      do iContA=1,nContA
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         TMP = TMP + ACC(iPrimA,iContA)*AUXarray2(iPrimA,iP)
        enddo
        AUXarrayCont(iContA,iP) = TMP
      enddo
    enddo
   end subroutine PrimitiveContractionAGPUSegQ


   subroutine PrimitiveContractionBGPUSegQ(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nTUVP,nTUVQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nContA,nPrimB,nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iContB,iContA,iPrimB
    real(realk) :: TMP
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContA,iContB,iPrimB,TMP) &
!$ACC PRESENT(nPasses,nPrimB,nContB,nContA,&
!$ACC         BCC,AUXarrayCont,AUXarray2) ASYNC(iASync)
    do iP=1,nPasses*nTUVP*nTUVQ
     do iContB=1,nContB
      do iContA=1,nContA
       AUXarrayCont(iContA,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
       do iContA=1,nContA
        AUXarrayCont(iContA,iContB,iP) = AUXarrayCont(iContA,iContB,iP) + TMP*AUXarray2(iContA,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
   end subroutine PrimitiveContractionBGPUSegQ

END MODULE AGC_GPU_PrimContractSegQMod
