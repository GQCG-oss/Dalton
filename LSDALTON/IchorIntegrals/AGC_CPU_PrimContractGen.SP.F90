module SPAGC_CPU_PrimContractGenMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
  
CONTAINS
  

   subroutine SPPrimitiveContractionCCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iContC,iPrimC
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      TMP = 0.0E0_reals
      do iPrimC=1,nPrimC
       TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
      enddo
      AUXarrayCont(iContC,iP)=TMP
     enddo
    enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionCCPUGen1

   subroutine SPPrimitiveContractionDCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iContD,iContC,iPrimD
    integer :: iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContD,iContC,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
      do iContC=1,nContC
       AUXarrayCont(iContC,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen1

   subroutine SPPrimitiveContractionACPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iC,iP,iContA,iPrimA
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iC,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
      do iC=1,nContC*nContD
       AUXarrayCont(iC,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen1

   subroutine SPPrimitiveContractionBCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iC,iP,iContB,iPrimB
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iC,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
      do iC=1,nContC*nContD*nContA
       AUXarrayCont(iC,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen1

   subroutine SPPrimitiveContractionCCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(    4,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(    4,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,    4
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen4

   subroutine SPPrimitiveContractionDCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(    4*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(    4*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,4*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen4

   subroutine SPPrimitiveContractionACPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(    4*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(    4*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,4*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen4

   subroutine SPPrimitiveContractionBCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(    4*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(    4*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,4*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen4

   subroutine SPPrimitiveContractionCCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   10,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   10,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   10
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen10

   subroutine SPPrimitiveContractionDCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   10*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   10*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,10*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen10

   subroutine SPPrimitiveContractionACPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   10*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   10*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,10*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen10

   subroutine SPPrimitiveContractionBCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   10*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   10*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,10*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen10

   subroutine SPPrimitiveContractionCCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   20,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   20,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   20
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen20

   subroutine SPPrimitiveContractionDCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   20*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   20*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,20*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen20

   subroutine SPPrimitiveContractionACPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   20*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   20*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,20*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen20

   subroutine SPPrimitiveContractionBCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   20*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   20*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,20*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen20

   subroutine SPPrimitiveContractionCCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   35,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   35,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   35
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen35

   subroutine SPPrimitiveContractionDCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   35*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   35*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,35*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen35

   subroutine SPPrimitiveContractionACPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   35*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   35*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,35*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen35

   subroutine SPPrimitiveContractionBCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   35*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   35*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,35*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen35

   subroutine SPPrimitiveContractionCCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   16,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   16,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   16
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen16

   subroutine SPPrimitiveContractionDCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   16*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   16*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,16*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen16

   subroutine SPPrimitiveContractionACPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   16*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   16*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,16*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen16

   subroutine SPPrimitiveContractionBCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   16*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   16*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,16*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen16

   subroutine SPPrimitiveContractionCCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   40,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   40,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   40
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen40

   subroutine SPPrimitiveContractionDCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   40*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   40*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,40*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen40

   subroutine SPPrimitiveContractionACPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   40*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   40*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,40*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen40

   subroutine SPPrimitiveContractionBCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   40*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   40*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,40*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen40

   subroutine SPPrimitiveContractionCCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   80,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   80,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   80
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen80

   subroutine SPPrimitiveContractionDCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   80*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   80*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,80*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen80

   subroutine SPPrimitiveContractionACPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   80*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   80*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,80*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen80

   subroutine SPPrimitiveContractionBCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   80*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   80*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,80*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen80

   subroutine SPPrimitiveContractionCCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  140,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  140,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  140
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen140

   subroutine SPPrimitiveContractionDCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  140*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  140*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,140*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen140

   subroutine SPPrimitiveContractionACPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  140*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  140*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,140*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen140

   subroutine SPPrimitiveContractionBCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  140*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  140*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,140*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen140

   subroutine SPPrimitiveContractionCCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  100,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  100,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  100
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen100

   subroutine SPPrimitiveContractionDCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  100*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  100*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,100*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen100

   subroutine SPPrimitiveContractionACPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  100*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  100*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,100*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen100

   subroutine SPPrimitiveContractionBCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  100*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  100*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,100*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen100

   subroutine SPPrimitiveContractionCCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  200,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  200,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  200
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen200

   subroutine SPPrimitiveContractionDCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  200*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  200*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,200*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen200

   subroutine SPPrimitiveContractionACPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  200*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  200*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,200*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen200

   subroutine SPPrimitiveContractionBCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  200*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  200*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,200*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen200

   subroutine SPPrimitiveContractionCCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  350,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  350,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  350
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen350

   subroutine SPPrimitiveContractionDCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  350*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  350*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,350*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen350

   subroutine SPPrimitiveContractionACPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  350*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  350*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,350*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen350

   subroutine SPPrimitiveContractionBCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  350*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  350*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,350*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen350

   subroutine SPPrimitiveContractionCCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  400,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  400,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  400
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen400

   subroutine SPPrimitiveContractionDCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  400*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  400*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,400*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen400

   subroutine SPPrimitiveContractionACPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  400*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  400*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,400*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen400

   subroutine SPPrimitiveContractionBCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  400*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  400*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,400*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen400

   subroutine SPPrimitiveContractionCCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  700,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  700,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  700
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen700

   subroutine SPPrimitiveContractionDCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  700*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  700*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,700*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen700

   subroutine SPPrimitiveContractionACPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  700*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  700*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,700*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen700

   subroutine SPPrimitiveContractionBCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  700*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  700*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,700*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen700

   subroutine SPPrimitiveContractionCCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont( 1225,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(reals) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1, 1225
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUGen1225

   subroutine SPPrimitiveContractionDCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2( 1225*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont( 1225*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(reals) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,1225*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUGen1225

   subroutine SPPrimitiveContractionACPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2( 1225*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont( 1225*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(reals) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,1225*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionACPUGen1225

   subroutine SPPrimitiveContractionBCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2( 1225*nContC*nContD*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont( 1225*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(reals) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,1225*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionBCPUGen1225
END module SPAGC_CPU_PrimContractGenMod
