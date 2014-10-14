MODULE AGC_CPU_PrimContractSegQMod
!Automatic Generated Code (AGC) by runPrimContraction.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
  
CONTAINS
  

   subroutine PrimitiveContractionACPUSegQ1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     TMP = 0.0E0_realk
     do iPrimA=1,nPrimA
      TMP = TMP + ACC(iPrimA,iContA)*AUXarray2(iPrimA,iP)
     enddo
     AUXarrayCont(iContA,iP) = tmp
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ1


   subroutine PrimitiveContractionBCPUSegQ1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iContA
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iContA,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iContA=1,nContA
       AUXarrayCont(iContA,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ1


   subroutine PrimitiveContractionACPUSegQ4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(    4,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,    4
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ4


   subroutine PrimitiveContractionBCPUSegQ4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(    4*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,    4*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ4


   subroutine PrimitiveContractionACPUSegQ10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   10,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   10
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ10


   subroutine PrimitiveContractionBCPUSegQ10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   10*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   10*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ10


   subroutine PrimitiveContractionACPUSegQ20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   20,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   20
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ20


   subroutine PrimitiveContractionBCPUSegQ20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   20*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   20*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ20


   subroutine PrimitiveContractionACPUSegQ35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   35,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   35
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ35


   subroutine PrimitiveContractionBCPUSegQ35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   35*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   35*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ35


   subroutine PrimitiveContractionACPUSegQ16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   16,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   16
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ16


   subroutine PrimitiveContractionBCPUSegQ16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   16*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ16


   subroutine PrimitiveContractionACPUSegQ40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   40,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   40
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ40


   subroutine PrimitiveContractionBCPUSegQ40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   40*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   40*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ40


   subroutine PrimitiveContractionACPUSegQ80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   80,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   80
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ80


   subroutine PrimitiveContractionBCPUSegQ80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   80*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   80*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ80


   subroutine PrimitiveContractionACPUSegQ140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  140,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  140
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ140


   subroutine PrimitiveContractionBCPUSegQ140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  140*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  140*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ140


   subroutine PrimitiveContractionACPUSegQ100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  100,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  100
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ100


   subroutine PrimitiveContractionBCPUSegQ100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  100*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ100


   subroutine PrimitiveContractionACPUSegQ200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  200,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  200
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ200


   subroutine PrimitiveContractionBCPUSegQ200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  200*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  200*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ200


   subroutine PrimitiveContractionACPUSegQ350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  350,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  350
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ350


   subroutine PrimitiveContractionBCPUSegQ350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  350*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  350*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ350


   subroutine PrimitiveContractionACPUSegQ400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  400,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  400
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ400


   subroutine PrimitiveContractionBCPUSegQ400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  400*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ400


   subroutine PrimitiveContractionACPUSegQ700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  700,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  700
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ700


   subroutine PrimitiveContractionBCPUSegQ700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  700*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  700*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ700


   subroutine PrimitiveContractionACPUSegQ1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1, 1225
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
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
   end subroutine PrimitiveContractionACPUSegQ1225


   subroutine PrimitiveContractionBCPUSegQ1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1, 1225*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
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
   end subroutine PrimitiveContractionBCPUSegQ1225

END MODULE AGC_CPU_PrimContractSegQMod
