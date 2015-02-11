MODULE AGC_GPU_PrimContractSegPMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
  
CONTAINS
  

   subroutine PrimitiveContractionCGPUSegP(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ,nTUVP,nTUVQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iContC,iPrimC
    real(realk) :: TMP
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContC,iPrimC,TMP) &
!$ACC PRESENT(nContC,nPrimC,nPasses,nPrimD,&
!$ACC        CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)
    do iP = 1,nPasses*nTUVP*nTUVQ*nPrimD
     do iContC=1,nContC
      tmp = 0.0E0_realk
      do iPrimC=1,nPrimC
       tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
      enddo
      AUXarrayCont(iContC,iP) = tmp
     enddo
    enddo
   end subroutine PrimitiveContractionCGPUSegP

   subroutine PrimitiveContractionDGPUSegP(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ,nTUVP,nTUVQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iContC,iContD,iPrimD
    real(realk) :: TMP
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContD,iPrimD,iContC,TMP) &
!$ACC PRESENT(nContC,nPrimD,nPasses,nContD,&
!$ACC         DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)
    do iP = 1,nPasses*nTUVP*nTUVQ
     do iContD=1,nContD
      do iContC=1,nContC
       AUXarrayCont(iContC,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       tmp = DCC(nPrimD,nContD)
       do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iP) = AUXarrayCont(iContC,iContD,iP) + TMP*AUXarray2(iContC,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
   end subroutine PrimitiveContractionDGPUSegP
END MODULE AGC_GPU_PrimContractSegPMod
