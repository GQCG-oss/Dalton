MODULE AGC_CPU_PrimContractGenMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
  
CONTAINS
  

   subroutine PrimitiveContractionCCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iContC,iPrimC
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      TMP = 0.0E0_realk
      do iPrimC=1,nPrimC
       TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
      enddo
      AUXarrayCont(iContC,iP)=TMP
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen1

   subroutine PrimitiveContractionDCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iContD,iContC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContD,iContC,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
      do iContC=1,nContC
       AUXarrayCont(iContC,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
       do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iP)=AUXarrayCont(iContC,iContD,iP)+TMP*AUXarray2(iContC,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen1

   subroutine PrimitiveContractionACPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iC,iP,iContA,iPrimA
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iC,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
      do iC=1,nContC*nContD
       AUXarrayCont(iC,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
       do iC=1,nContC*nContD
        AUXarrayCont(iC,iContA,iP)=AUXarrayCont(iC,iContA,iP)+TMP*AUXarray2(iC,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen1

   subroutine PrimitiveContractionBCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iC,iP,iContB,iPrimB
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iC,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
      do iC=1,nContC*nContD*nContA
       AUXarrayCont(iC,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
       do iC=1,nContC*nContD*nContA
        AUXarrayCont(iC,iContB,iP)=AUXarrayCont(iC,iContB,iP)+TMP*AUXarray2(iC,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen1

   subroutine PrimitiveContractionCCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(    4,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,    4
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,    4
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen4

   subroutine PrimitiveContractionDCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(    4*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,4*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,4*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen4

   subroutine PrimitiveContractionACPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(    4*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,4*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,4*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen4

   subroutine PrimitiveContractionBCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(    4*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,4*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,4*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen4

   subroutine PrimitiveContractionCCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   10,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   10
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   10
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen10

   subroutine PrimitiveContractionDCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   10*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,10*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,10*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen10

   subroutine PrimitiveContractionACPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   10*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,10*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,10*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen10

   subroutine PrimitiveContractionBCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   10*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,10*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,10*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen10

   subroutine PrimitiveContractionCCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   20,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   20
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   20
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen20

   subroutine PrimitiveContractionDCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   20*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,20*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,20*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen20

   subroutine PrimitiveContractionACPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   20*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,20*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,20*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen20

   subroutine PrimitiveContractionBCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   20*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,20*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,20*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen20

   subroutine PrimitiveContractionCCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   35,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   35
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   35
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen35

   subroutine PrimitiveContractionDCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   35*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,35*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,35*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen35

   subroutine PrimitiveContractionACPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   35*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,35*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,35*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen35

   subroutine PrimitiveContractionBCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   35*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,35*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,35*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen35

   subroutine PrimitiveContractionCCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   16,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   16
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   16
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen16

   subroutine PrimitiveContractionDCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   16*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,16*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,16*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen16

   subroutine PrimitiveContractionACPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   16*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,16*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,16*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen16

   subroutine PrimitiveContractionBCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,16*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,16*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen16

   subroutine PrimitiveContractionCCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   40,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   40
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   40
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen40

   subroutine PrimitiveContractionDCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   40*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,40*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,40*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen40

   subroutine PrimitiveContractionACPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   40*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,40*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,40*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen40

   subroutine PrimitiveContractionBCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   40*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,40*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,40*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen40

   subroutine PrimitiveContractionCCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   80,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   80
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   80
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen80

   subroutine PrimitiveContractionDCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   80*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,80*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,80*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen80

   subroutine PrimitiveContractionACPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   80*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,80*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,80*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen80

   subroutine PrimitiveContractionBCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   80*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,80*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,80*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen80

   subroutine PrimitiveContractionCCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  140,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  140
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  140
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen140

   subroutine PrimitiveContractionDCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  140*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,140*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,140*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen140

   subroutine PrimitiveContractionACPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  140*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,140*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,140*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen140

   subroutine PrimitiveContractionBCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  140*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,140*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,140*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen140

   subroutine PrimitiveContractionCCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  100,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  100
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  100
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen100

   subroutine PrimitiveContractionDCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  100*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,100*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,100*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen100

   subroutine PrimitiveContractionACPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  100*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,100*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,100*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen100

   subroutine PrimitiveContractionBCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,100*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,100*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen100

   subroutine PrimitiveContractionCCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  200,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  200
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  200
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen200

   subroutine PrimitiveContractionDCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  200*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,200*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,200*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen200

   subroutine PrimitiveContractionACPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  200*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,200*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,200*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen200

   subroutine PrimitiveContractionBCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  200*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,200*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,200*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen200

   subroutine PrimitiveContractionCCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  350,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  350
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  350
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen350

   subroutine PrimitiveContractionDCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  350*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,350*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,350*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen350

   subroutine PrimitiveContractionACPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  350*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,350*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,350*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen350

   subroutine PrimitiveContractionBCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  350*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,350*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,350*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen350

   subroutine PrimitiveContractionCCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  400,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  400
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  400
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen400

   subroutine PrimitiveContractionDCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  400*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,400*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,400*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen400

   subroutine PrimitiveContractionACPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  400*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,400*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,400*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen400

   subroutine PrimitiveContractionBCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,400*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,400*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen400

   subroutine PrimitiveContractionCCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  700,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  700
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  700
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen700

   subroutine PrimitiveContractionDCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  700*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,700*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,700*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen700

   subroutine PrimitiveContractionACPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  700*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,700*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,700*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen700

   subroutine PrimitiveContractionBCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  700*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,700*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,700*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen700

   subroutine PrimitiveContractionCCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1, 1225
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1, 1225
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen1225

   subroutine PrimitiveContractionDCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2( 1225*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,1225*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,1225*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen1225

   subroutine PrimitiveContractionACPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2( 1225*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,1225*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,1225*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen1225

   subroutine PrimitiveContractionBCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,1225*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,1225*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen1225
END MODULE AGC_CPU_PrimContractGenMod
