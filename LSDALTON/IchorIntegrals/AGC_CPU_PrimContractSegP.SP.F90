module SPAGC_CPU_PrimContractSegPMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
  
CONTAINS
  

   subroutine SPPrimitiveContractionCCPUSegP1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      tmp = 0.0E0_reals
      do iPrimC=1,nPrimC
       tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
      enddo
      AUXarrayCont(iContC,iP) = tmp
     enddo
    enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionCCPUSegP1

   subroutine SPPrimitiveContractionDCPUSegP1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iContC
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iContC,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
      do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP1

   subroutine SPPrimitiveContractionCCPUSegP4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(    4,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(    4,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,    4
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP4

   subroutine SPPrimitiveContractionDCPUSegP4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(    4*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(    4*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,    4*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP4

   subroutine SPPrimitiveContractionCCPUSegP10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   10,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   10,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   10
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP10

   subroutine SPPrimitiveContractionDCPUSegP10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   10*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   10*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   10*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP10

   subroutine SPPrimitiveContractionCCPUSegP20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   20,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   20,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   20
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP20

   subroutine SPPrimitiveContractionDCPUSegP20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   20*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   20*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   20*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP20

   subroutine SPPrimitiveContractionCCPUSegP35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   35,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   35,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   35
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP35

   subroutine SPPrimitiveContractionDCPUSegP35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   35*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   35*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   35*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP35

   subroutine SPPrimitiveContractionCCPUSegP16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   16,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   16,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   16
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP16

   subroutine SPPrimitiveContractionDCPUSegP16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   16*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   16*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   16*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP16

   subroutine SPPrimitiveContractionCCPUSegP40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   40,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   40,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   40
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP40

   subroutine SPPrimitiveContractionDCPUSegP40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   40*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   40*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   40*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP40

   subroutine SPPrimitiveContractionCCPUSegP80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(   80,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   80,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   80
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP80

   subroutine SPPrimitiveContractionDCPUSegP80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(   80*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   80*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   80*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP80

   subroutine SPPrimitiveContractionCCPUSegP140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  140,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  140,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  140
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP140

   subroutine SPPrimitiveContractionDCPUSegP140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  140*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  140*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  140*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP140

   subroutine SPPrimitiveContractionCCPUSegP100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  100,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  100,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  100
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP100

   subroutine SPPrimitiveContractionDCPUSegP100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  100*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  100*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  100*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP100

   subroutine SPPrimitiveContractionCCPUSegP200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  200,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  200,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  200
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP200

   subroutine SPPrimitiveContractionDCPUSegP200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  200*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  200*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  200*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP200

   subroutine SPPrimitiveContractionCCPUSegP350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  350,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  350,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  350
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP350

   subroutine SPPrimitiveContractionDCPUSegP350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  350*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  350*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  350*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP350

   subroutine SPPrimitiveContractionCCPUSegP400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  400,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  400,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  400
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP400

   subroutine SPPrimitiveContractionDCPUSegP400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  400*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  400*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  400*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP400

   subroutine SPPrimitiveContractionCCPUSegP700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2(  700,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  700,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  700
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP700

   subroutine SPPrimitiveContractionDCPUSegP700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2(  700*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  700*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  700*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP700

   subroutine SPPrimitiveContractionCCPUSegP1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: CCC(nPrimC,nContC)
    real(reals),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD*nPasses)
    real(reals),intent(inout) :: AUXarrayCont( 1225,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1, 1225
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionCCPUSegP1225

   subroutine SPPrimitiveContractionDCPUSegP1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(reals),intent(in) :: DCC(nPrimD,nContD)
    real(reals),intent(in) :: AUXarray2( 1225*nContC,nPrimD,nPasses)
    real(reals),intent(inout) :: AUXarrayCont( 1225*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1, 1225*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_reals
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
   end subroutine SPPrimitiveContractionDCPUSegP1225
END module SPAGC_CPU_PrimContractSegPMod
