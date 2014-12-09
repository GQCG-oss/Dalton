MODULE IchorEriGabPrimMod
!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
 
CONTAINS
  subroutine GabPrimitiveContractionGen1A(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont(nPrimB,nPrimB,nContA)
    !
    integer :: iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: TMP,TMPACC
    !$OMP DO COLLAPSE(3) &
    !$OMP PRIVATE(iPrimC,iPrimD,iPrimA,iPrimB,iContD,TMP,TMPACC) 
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         TMPACC = ACC(iPrimA,iContC)
         do iPrimC=1,nPrimA
          TMP = TMP + TMPACC*ACC(iPrimC,iContC)*AUXarray2(iPrimC,iPrimD,iPrimA,iPrimB)
         enddo
        enddo
        AUXarrayCont(iPrimD,iPrimB,iContC) = TMP
       enddo
      enddo
     enddo
     !$OMP ENDDO 
  end subroutine GabPrimitiveContractionGen1A

  subroutine GabPrimitiveContractionGen1B(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nPrimB,nPrimB,nContA)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB)
    !
    integer :: iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: TMP,TMPBCC
    !$OMP DO COLLAPSE(2) PRIVATE(iContC,iPrimD,iPrimB,iContD,TMP,TMPBCC) 
     do iContC=1,nContA
      do iContD=1,nContB
       TMP = 0.0E0_realk
       do iPrimB=1,nPrimB
        TMPBCC = BCC(iPrimB,iContD)
        do iPrimD=1,nPrimB
         TMP = TMP + TMPBCC*BCC(iPrimD,iContD)*AUXarray2(iPrimD,iPrimB,iContC)
        enddo
       enddo
       AUXarrayCont(iContC,iContD) = TMP
      enddo
     enddo
     !$OMP END DO
  end subroutine GabPrimitiveContractionGen1B

  subroutine GabPrimitiveContractionGen16A(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   16,nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont(   16,nPrimB,nPrimB,nContA)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD,iTUV
    real(realk) :: TMP,ACCTMP
!$OMP DO COLLAPSE(4) &     
!$OMP PRIVATE(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iContD,iContC,TMP,ACCTMP)
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,   16
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         AUXarrayCont(iTUV,iPrimD,iPrimB,iContC) = TMP
        enddo
       enddo
      enddo
     enddo
!$OMP END DO
  end subroutine GabPrimitiveContractionGen16A

  subroutine GabPrimitiveContractionGen16B(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16,nPrimB,nPrimB,nContA)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContA,nContB)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD,iTUV
    real(realk) :: TMP,BCCTMP
!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iTUV,iPrimD,iPrimB,iContD,iContC,TMP,BCCTMP) 
     do iContC=1,nContA
      do iContD=1,nContB
       do iTUV=1,   16
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         BCCTMP = BCC(iPrimB,iContD)
         do iPrimD=1,nPrimB
          TMP = TMP + BCC(iPrimD,iContD)*BCCTMP*AUXarray2(iTUV,iPrimD,iPrimB,iContC)
         enddo
        enddo
        AUXarrayCont(iTUV,iContC,iContD) = TMP
       enddo
      enddo
     enddo
!$OMP END DO
  end subroutine GabPrimitiveContractionGen16B

  subroutine GabPrimitiveContractionGen100A(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  100,nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont(  100,nPrimB,nPrimB,nContA)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD,iTUV
    real(realk) :: TMP,ACCTMP
!$OMP DO COLLAPSE(4) &     
!$OMP PRIVATE(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iContD,iContC,TMP,ACCTMP)
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,  100
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         AUXarrayCont(iTUV,iPrimD,iPrimB,iContC) = TMP
        enddo
       enddo
      enddo
     enddo
!$OMP END DO
  end subroutine GabPrimitiveContractionGen100A

  subroutine GabPrimitiveContractionGen100B(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100,nPrimB,nPrimB,nContA)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContA,nContB)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD,iTUV
    real(realk) :: TMP,BCCTMP
!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iTUV,iPrimD,iPrimB,iContD,iContC,TMP,BCCTMP) 
     do iContC=1,nContA
      do iContD=1,nContB
       do iTUV=1,  100
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         BCCTMP = BCC(iPrimB,iContD)
         do iPrimD=1,nPrimB
          TMP = TMP + BCC(iPrimD,iContD)*BCCTMP*AUXarray2(iTUV,iPrimD,iPrimB,iContC)
         enddo
        enddo
        AUXarrayCont(iTUV,iContC,iContD) = TMP
       enddo
      enddo
     enddo
!$OMP END DO
  end subroutine GabPrimitiveContractionGen100B

  subroutine GabPrimitiveContractionGen400A(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  400,nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont(  400,nPrimB,nPrimB,nContA)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD,iTUV
    real(realk) :: TMP,ACCTMP
!$OMP DO COLLAPSE(4) &     
!$OMP PRIVATE(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iContD,iContC,TMP,ACCTMP)
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,  400
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         AUXarrayCont(iTUV,iPrimD,iPrimB,iContC) = TMP
        enddo
       enddo
      enddo
     enddo
!$OMP END DO
  end subroutine GabPrimitiveContractionGen400A

  subroutine GabPrimitiveContractionGen400B(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400,nPrimB,nPrimB,nContA)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContA,nContB)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD,iTUV
    real(realk) :: TMP,BCCTMP
!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iTUV,iPrimD,iPrimB,iContD,iContC,TMP,BCCTMP) 
     do iContC=1,nContA
      do iContD=1,nContB
       do iTUV=1,  400
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         BCCTMP = BCC(iPrimB,iContD)
         do iPrimD=1,nPrimB
          TMP = TMP + BCC(iPrimD,iContD)*BCCTMP*AUXarray2(iTUV,iPrimD,iPrimB,iContC)
         enddo
        enddo
        AUXarrayCont(iTUV,iContC,iContD) = TMP
       enddo
      enddo
     enddo
!$OMP END DO
  end subroutine GabPrimitiveContractionGen400B

  subroutine GabPrimitiveContractionGen1225A(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nPrimB,nPrimB,nContA)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD,iTUV
    real(realk) :: TMP,ACCTMP
!$OMP DO COLLAPSE(4) &     
!$OMP PRIVATE(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iContD,iContC,TMP,ACCTMP)
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1, 1225
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         AUXarrayCont(iTUV,iPrimD,iPrimB,iContC) = TMP
        enddo
       enddo
      enddo
     enddo
!$OMP END DO
  end subroutine GabPrimitiveContractionGen1225A

  subroutine GabPrimitiveContractionGen1225B(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimB,nPrimB,nContA)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContA,nContB)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD,iTUV
    real(realk) :: TMP,BCCTMP
!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iTUV,iPrimD,iPrimB,iContD,iContC,TMP,BCCTMP) 
     do iContC=1,nContA
      do iContD=1,nContB
       do iTUV=1, 1225
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         BCCTMP = BCC(iPrimB,iContD)
         do iPrimD=1,nPrimB
          TMP = TMP + BCC(iPrimD,iContD)*BCCTMP*AUXarray2(iTUV,iPrimD,iPrimB,iContC)
         enddo
        enddo
        AUXarrayCont(iTUV,iContC,iContD) = TMP
       enddo
      enddo
     enddo
!$OMP END DO
  end subroutine GabPrimitiveContractionGen1225B

  subroutine ExtractGabElmP1Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP
    real(realk) :: TMP
    !$OMP SINGLE
     TMP = ABS(AUXarray(1))
     do iContP=2,nContP
      IF(ABS(AUXarray(iContP)).GT.TMP)THEN
       TMP = ABS(AUXarray(iContP))
      ENDIF
     enddo
     Output(1) = SQRT(TMP)
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine ExtractGabElmP1Gen


  subroutine ExtractGabElmP3Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(    3,    3,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(    3)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    3
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    3
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine ExtractGabElmP3Gen

  subroutine ExtractGabElmP5Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(    5,    5,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(    5)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    5
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    5
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine ExtractGabElmP5Gen

  subroutine ExtractGabElmP9Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(    9,    9,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(    9)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    9
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    9
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine ExtractGabElmP9Gen

  subroutine ExtractGabElmP15Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(   15,   15,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(   15)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,   15
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,   15
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine ExtractGabElmP15Gen

  subroutine ExtractGabElmP25Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(   25,   25,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(   25)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,   25
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,   25
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine ExtractGabElmP25Gen
END MODULE IchorEriGabPrimMod
