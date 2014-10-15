MODULE AGC_GPU_PrimContractGenMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
  
CONTAINS
  
   subroutine PrimitiveContractionCGPUGen(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ,nTUVP,nTUVQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPrimA*nPrimB*nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPrimA*nPrimB*nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iContC,iPrimC
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContC,iPrimC,TMP) &
!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!$ACC        CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)
   do iP=1,nTUVP*nTUVQ*nPrimA*nPrimB*nPrimD*nPasses
    do iContC=1,nContC
     TMP = 0.0E0_realk
     do iPrimC=1,nPrimC
      TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
     enddo
     AUXarrayCont(iContC,iP)=TMP
    enddo
   enddo
   end subroutine PrimitiveContractionCGPUGen

   subroutine PrimitiveContractionDGPUGen(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ,nTUVP,nTUVQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPrimA*nPrimB*nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPrimA*nPrimB*nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iContC,iContD,iPrimD
    real(realk) :: TMP
!Scaling p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContC,iContD,iPrimD,TMP) &
!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!$ACC         DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)
   do iP=1,nTUVP*nTUVQ*nPrimA*nPrimB*nPasses
    do iContD=1,nContD
     do iContC=1,nContC
      AUXarrayCont(iContC,iContD,iP)=0.0E0_realk
     enddo
     do iPrimD=1,nPrimD
      TMP = DCC(iPrimD,iContD)
      do iContC=1,nContC
       AUXarrayCont(iContC,iContD,iP)=AUXarrayCont(iContC,iContD,iP)+TMP*AUXarray2(iContC,iPrimD,iP)
      enddo
     enddo
    enddo
   enddo
   end subroutine PrimitiveContractionDGPUGen

   subroutine PrimitiveContractionAGPUGen(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ,nTUVP,nTUVQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(nContC*nContD,nPrimA,nPrimB*nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD,nContA,nPrimB*nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iC,iContA,iPrimA
    real(realk) :: TMP
!Scaling p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContC*nContD*nContA*nTUV*nPassQ
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContA,iPrimA,TMP,iC) &
!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!$ACC         ACC,AUXarrayCont,AUXarray2) ASYNC(iASync)
   do iP=1,nTUVP*nTUVQ*nPrimB*nPasses
    do iContA=1,nContA
     do iC=1,nContC*nContD
      AUXarrayCont(iC,iContA,iP)=0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iC=1,nContC*nContD
       AUXarrayCont(iC,iContA,iP)=AUXarrayCont(iC,iContA,iP)+TMP*AUXarray2(iC,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
   end subroutine PrimitiveContractionAGPUGen

   subroutine PrimitiveContractionBGPUGen(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ,nTUVP,nTUVQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nContC*nContD*nContA,nPrimB,nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD*nContA,nContB,nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iC,iContB,iPrimB
    real(realk) :: TMP
    !Scaling p*c**4*nTUV*nPassQ: nPrimB*nContC*nContD*nContA*nContB*nTUV*nPassQ
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContB,iPrimB,TMP,iC) &
!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!$ACC        BCC,AUXarrayCont,AUXarray2) ASYNC(iASync)
   do iP=1,nTUVP*nTUVQ*nPasses
    do iContB=1,nContB
     do iC=1,nContC*nContD*nContA
      AUXarrayCont(iC,iContB,iP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iC=1,nContC*nContD*nContA
       AUXarrayCont(iC,iContB,iP)=AUXarrayCont(iC,iContB,iP)+TMP*AUXarray2(iC,iPrimB,iP)
      enddo
     enddo
    enddo
   enddo
   end subroutine PrimitiveContractionBGPUGen
END MODULE AGC_GPU_PrimContractGenMod
