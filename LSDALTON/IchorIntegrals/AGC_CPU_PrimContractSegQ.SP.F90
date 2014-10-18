module SPAGC_CPU_PrimContractSegQMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
  
CONTAINS
  

   subroutine SPPrimitiveContractionACPUSegQ1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     TMP = 0.0E0_reals
     do iPrimA=1,nPrimA
      TMP = TMP + ACC(iPrimA,iContA)*AUXarray2(iPrimA,iP)
     enddo
     AUXarrayCont(iContA,iP) = tmp
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ1


   subroutine SPPrimitiveContractionBCPUSegQ1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iContA
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iContA,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iContA=1,nContA
       AUXarrayCont(iContA,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iContA=1,nContA
       AUXarrayCont(iContA,iContB,iPassP)=AUXarrayCont(iContA,iContB,iPassP)+tmp*AUXarray2(iContA,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ1


   subroutine SPPrimitiveContractionACPUSegQ4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(    4,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(    4,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,    4
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,    4
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ4


   subroutine SPPrimitiveContractionBCPUSegQ4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(    4*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(    4*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,    4*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,    4*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ4


   subroutine SPPrimitiveContractionACPUSegQ10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   10,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   10,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   10
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   10
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ10


   subroutine SPPrimitiveContractionBCPUSegQ10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   10*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   10*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   10*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   10*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ10


   subroutine SPPrimitiveContractionACPUSegQ20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   20,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   20,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   20
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   20
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ20


   subroutine SPPrimitiveContractionBCPUSegQ20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   20*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   20*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   20*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   20*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ20


   subroutine SPPrimitiveContractionACPUSegQ35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   35,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   35,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   35
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   35
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ35


   subroutine SPPrimitiveContractionBCPUSegQ35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   35*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   35*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   35*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   35*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ35


   subroutine SPPrimitiveContractionACPUSegQ16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   16,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   16,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   16
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   16
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ16


   subroutine SPPrimitiveContractionBCPUSegQ16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   16*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   16*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   16*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   16*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ16


   subroutine SPPrimitiveContractionACPUSegQ40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   40,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   40,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   40
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   40
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ40


   subroutine SPPrimitiveContractionBCPUSegQ40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   40*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   40*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   40*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   40*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ40


   subroutine SPPrimitiveContractionACPUSegQ80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(   80,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   80,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   80
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   80
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ80


   subroutine SPPrimitiveContractionBCPUSegQ80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(   80*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(   80*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   80*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   80*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ80


   subroutine SPPrimitiveContractionACPUSegQ140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  140,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  140,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  140
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  140
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ140


   subroutine SPPrimitiveContractionBCPUSegQ140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  140*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  140*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  140*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  140*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ140


   subroutine SPPrimitiveContractionACPUSegQ100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  100,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  100,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  100
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  100
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ100


   subroutine SPPrimitiveContractionBCPUSegQ100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  100*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  100*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  100*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  100*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ100


   subroutine SPPrimitiveContractionACPUSegQ200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  200,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  200,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  200
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  200
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ200


   subroutine SPPrimitiveContractionBCPUSegQ200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  200*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  200*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  200*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  200*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ200


   subroutine SPPrimitiveContractionACPUSegQ350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  350,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  350,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  350
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  350
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ350


   subroutine SPPrimitiveContractionBCPUSegQ350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  350*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  350*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  350*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  350*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ350


   subroutine SPPrimitiveContractionACPUSegQ400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  400,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  400,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  400
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  400
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ400


   subroutine SPPrimitiveContractionBCPUSegQ400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  400*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  400*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  400*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  400*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ400


   subroutine SPPrimitiveContractionACPUSegQ700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2(  700,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  700,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  700
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  700
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ700


   subroutine SPPrimitiveContractionBCPUSegQ700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2(  700*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont(  700*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  700*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  700*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ700


   subroutine SPPrimitiveContractionACPUSegQ1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: ACC(nPrimA,nContA)
    real(reals),intent(in) :: AUXarray2( 1225,nPrimA,nPrimB*nPasses)
    real(reals),intent(inout) :: AUXarrayCont( 1225,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1, 1225
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_reals
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1, 1225
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionACPUSegQ1225


   subroutine SPPrimitiveContractionBCPUSegQ1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(reals),intent(in) :: BCC(nPrimB,nContB)
    real(reals),intent(in) :: AUXarray2( 1225*nContA,nPrimB,nPasses)
    real(reals),intent(inout) :: AUXarrayCont( 1225*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(reals) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1, 1225*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_reals
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1, 1225*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine SPPrimitiveContractionBCPUSegQ1225

END module SPAGC_CPU_PrimContractSegQMod
