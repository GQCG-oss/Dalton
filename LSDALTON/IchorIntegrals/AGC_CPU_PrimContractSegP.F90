MODULE AGC_CPU_PrimContractSegPMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
  
CONTAINS
  

   subroutine PrimitiveContractionCCPUSegP1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      tmp = 0.0E0_realk
      do iPrimC=1,nPrimC
       tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
      enddo
      AUXarrayCont(iContC,iP) = tmp
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP1

   subroutine PrimitiveContractionDCPUSegP1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iContC
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iContC,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
      do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iPassP) = AUXarrayCont(iContC,iContD,iPassP) + tmp*AUXarray2(iContC,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP1

   subroutine PrimitiveContractionCCPUSegP4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(    4,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,    4
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,    4
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP4

   subroutine PrimitiveContractionDCPUSegP4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(    4*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,    4*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,    4*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP4

   subroutine PrimitiveContractionCCPUSegP10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   10,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   10
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   10
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP10

   subroutine PrimitiveContractionDCPUSegP10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   10*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   10*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   10*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP10

   subroutine PrimitiveContractionCCPUSegP20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   20,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   20
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   20
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP20

   subroutine PrimitiveContractionDCPUSegP20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   20*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   20*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   20*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP20

   subroutine PrimitiveContractionCCPUSegP35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   35,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   35
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   35
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP35

   subroutine PrimitiveContractionDCPUSegP35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   35*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   35*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   35*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP35

   subroutine PrimitiveContractionCCPUSegP16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   16,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   16
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   16
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP16

   subroutine PrimitiveContractionDCPUSegP16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   16*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   16*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   16*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP16

   subroutine PrimitiveContractionCCPUSegP40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   40,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   40
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   40
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP40

   subroutine PrimitiveContractionDCPUSegP40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   40*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   40*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   40*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP40

   subroutine PrimitiveContractionCCPUSegP80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   80,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   80
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   80
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP80

   subroutine PrimitiveContractionDCPUSegP80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   80*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   80*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   80*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP80

   subroutine PrimitiveContractionCCPUSegP140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  140,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  140
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  140
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP140

   subroutine PrimitiveContractionDCPUSegP140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  140*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  140*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  140*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP140

   subroutine PrimitiveContractionCCPUSegP100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  100,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  100
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  100
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP100

   subroutine PrimitiveContractionDCPUSegP100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  100*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  100*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  100*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP100

   subroutine PrimitiveContractionCCPUSegP200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  200,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  200
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  200
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP200

   subroutine PrimitiveContractionDCPUSegP200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  200*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  200*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  200*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP200

   subroutine PrimitiveContractionCCPUSegP350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  350,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  350
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  350
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP350

   subroutine PrimitiveContractionDCPUSegP350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  350*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  350*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  350*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP350

   subroutine PrimitiveContractionCCPUSegP400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  400,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  400
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  400
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP400

   subroutine PrimitiveContractionDCPUSegP400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  400*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  400*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  400*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP400

   subroutine PrimitiveContractionCCPUSegP700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  700,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  700
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  700
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP700

   subroutine PrimitiveContractionDCPUSegP700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  700*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  700*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  700*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP700

   subroutine PrimitiveContractionCCPUSegP1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1, 1225
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1, 1225
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP1225

   subroutine PrimitiveContractionDCPUSegP1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2( 1225*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1, 1225*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1, 1225*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP1225
END MODULE AGC_CPU_PrimContractSegPMod
