MODULE IchorEriGabintegralOBSGeneralModSeg
!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_OBS_VERTICALRECURRENCEMODA
use AGC_OBS_VERTICALRECURRENCEMODB
use AGC_OBS_VERTICALRECURRENCEMODD
use AGC_OBS_VERTICALRECURRENCEMODC
use AGC_OBS_VERTICALRECURRENCEMODASeg
use AGC_OBS_VERTICALRECURRENCEMODBSeg
use AGC_OBS_VERTICALRECURRENCEMODDSeg
use AGC_OBS_VERTICALRECURRENCEMODCSeg
use AGC_OBS_TRANSFERRECURRENCEMODAtoCSeg
use AGC_OBS_TRANSFERRECURRENCEMODAtoDSeg
use AGC_OBS_TRANSFERRECURRENCEMODBtoCSeg
use AGC_OBS_TRANSFERRECURRENCEMODBtoDSeg
use AGC_OBS_TRANSFERRECURRENCEMODCtoASeg
use AGC_OBS_TRANSFERRECURRENCEMODDtoASeg
use AGC_OBS_TRANSFERRECURRENCEMODCtoBSeg
use AGC_OBS_TRANSFERRECURRENCEMODDtoBSeg
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorGabIntegral_OBS_Seg,IchorGabIntegral_OBS_general_sizeSeg  
  
CONTAINS
  
  
  subroutine IchorGabIntegral_OBS_Seg(nPrimA,nPrimB,&
       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContP,pexp,ACC,BCC,&
       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,Pdistance12,PQorder,CDAB,Acenter,Bcenter,&
       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    implicit none
    integer,intent(in) :: nPrimP,nPasses,nPrimA,nPrimB
    integer,intent(in) :: MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContP,nTABFJW1,nTABFJW2
    integer,intent(in) :: AngmomA,AngmomB
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)
    logical,intent(in)     :: Psegmented
    real(realk),intent(in) :: pexp(nPrimP)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(inout) :: CDAB(nPasses)
    real(realk),intent(in) :: integralPrefactor(nPrimP*nPrimP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimP)
    real(realk),intent(in) :: reducedExponents(nPrimP*nPrimP)
    !reducedExponents(nPrimP,nPrimP)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
!   Local variables 
    integer :: AngmomP,I,J,nContQP,la,lb,lc,ld,nsize,angmomid
    
    !Setup combined Angmom info
    AngmomP = AngmomA+AngmomB
!    nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
!    nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
!    nlmA = 2*AngmomA+1
!    nlmB = 2*AngmomB+1
    AngmomID = 10*AngmomA+AngmomB
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSeg0(nPasses,nPrimP,nPrimP,&
               & reducedExponents,TABFJW,Pcent,Pcent,integralPrefactor,&
               & PpreExpFac,PpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        call ExtractGabElmP1Seg(TMParray2,CDAB,nPasses)
    CASE(  10)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoCSeg(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,4,Pdistance12,TMParray1,TMParray2,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,3,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP3Seg(TMParray1,CDAB,nPasses)
    CASE(  11)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSeg(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,10,Pdistance12,TMParray1,TMParray2,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,9,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP9Seg(TMParray1,CDAB,nPasses)
    CASE(  20)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSeg(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,10,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nPasses,TMParray2,TMParray1)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,5,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nPasses,TMParray2,TMParray1)
        call ExtractGabElmP5Seg(TMParray1,CDAB,nPasses)
    CASE(  21)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoCSeg(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,20,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nPasses,TMParray2,TMParray1)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,15,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nPasses,TMParray2,TMParray1)
        call ExtractGabElmP15Seg(TMParray1,CDAB,nPasses)
    CASE(  22)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoCSeg(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,35,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,nPasses,TMParray2,TMParray1)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,25,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nPasses,TMParray2,TMParray1)
        call ExtractGabElmP25Seg(TMParray1,CDAB,nPasses)
    CASE(   1)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoDSeg(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,4,Pdistance12,TMParray1,TMParray2,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,3,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP3Seg(TMParray1,CDAB,nPasses)
    CASE(   2)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoDSeg(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,10,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nPasses,TMParray2,TMParray1)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,5,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nPasses,TMParray2,TMParray1)
        call ExtractGabElmP5Seg(TMParray1,CDAB,nPasses)
    CASE(  12)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoDSeg(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,20,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nPasses,TMParray2,TMParray1)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,15,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nPasses,TMParray2,TMParray1)
        call ExtractGabElmP15Seg(TMParray1,CDAB,nPasses)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorGabIntegral_OBS_Seg',-1)
    END SELECT
  end subroutine IchorGabIntegral_OBS_Seg
  
  subroutine IchorGabIntegral_OBS_general_sizeSeg(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,nPrimP,nContP)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: AngmomA,AngmomB
    integer,intent(in) :: nPrimP,nContP
    ! local variables
    integer :: AngmomID
    
    AngmomID = 10*AngmomA+AngmomB
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,1)
    CASE(   1)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,16)
       TMParray2maxSize = MAX(TMParray2maxSize,12)
       TMParray1maxSize = MAX(TMParray1maxSize,9)
    CASE(   2)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
       TMParray1maxSize = MAX(TMParray1maxSize,50)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
       TMParray1maxSize = MAX(TMParray1maxSize,25)
    CASE(  10)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,16)
       TMParray2maxSize = MAX(TMParray2maxSize,12)
       TMParray1maxSize = MAX(TMParray1maxSize,9)
    CASE(  11)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,90)
       TMParray1maxSize = MAX(TMParray1maxSize,81)
    CASE(  12)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,400)
       TMParray2maxSize = MAX(TMParray2maxSize,360)
       TMParray1maxSize = MAX(TMParray1maxSize,300)
       TMParray2maxSize = MAX(TMParray2maxSize,270)
       TMParray1maxSize = MAX(TMParray1maxSize,225)
    CASE(  20)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
       TMParray1maxSize = MAX(TMParray1maxSize,50)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
       TMParray1maxSize = MAX(TMParray1maxSize,25)
    CASE(  21)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,400)
       TMParray2maxSize = MAX(TMParray2maxSize,360)
       TMParray1maxSize = MAX(TMParray1maxSize,300)
       TMParray2maxSize = MAX(TMParray2maxSize,270)
       TMParray1maxSize = MAX(TMParray1maxSize,225)
    CASE(  22)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225)
       TMParray2maxSize = MAX(TMParray2maxSize,1260)
       TMParray1maxSize = MAX(TMParray1maxSize,875)
       TMParray2maxSize = MAX(TMParray2maxSize,900)
       TMParray1maxSize = MAX(TMParray1maxSize,625)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorGabIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorGabIntegral_OBS_general_sizeSeg
  
  subroutine ExtractGabElmP1Seg(AUXarray,Output,nPasses)
    implicit none
    integer,intent(in) :: nPasses
    real(realk),intent(in) :: AUXarray(nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ
    do iPassQ = 1,nPasses
     Output(iPassQ) = SQRT(ABS(AUXarray(iPassQ)))
    enddo
  end subroutine ExtractGabElmP1Seg

  subroutine ExtractGabElmP3Seg(AUXarray,Output,nPasses)
    implicit none
    integer,intent(in) :: nPasses
    real(realk),intent(in) :: AUXarray(    3,    3,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,i
    real(realk) :: TMP(    3)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,    3
      TMP(i) = ABS(AUXarray(i,i,iPassQ))
     enddo
     Output(iPassQ) = SQRT(MAXVAL(TMP))
    enddo
  end subroutine ExtractGabElmP3Seg

  subroutine ExtractGabElmP5Seg(AUXarray,Output,nPasses)
    implicit none
    integer,intent(in) :: nPasses
    real(realk),intent(in) :: AUXarray(    5,    5,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,i
    real(realk) :: TMP(    5)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,    5
      TMP(i) = ABS(AUXarray(i,i,iPassQ))
     enddo
     Output(iPassQ) = SQRT(MAXVAL(TMP))
    enddo
  end subroutine ExtractGabElmP5Seg

  subroutine ExtractGabElmP9Seg(AUXarray,Output,nPasses)
    implicit none
    integer,intent(in) :: nPasses
    real(realk),intent(in) :: AUXarray(    9,    9,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,i
    real(realk) :: TMP(    9)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,    9
      TMP(i) = ABS(AUXarray(i,i,iPassQ))
     enddo
     Output(iPassQ) = SQRT(MAXVAL(TMP))
    enddo
  end subroutine ExtractGabElmP9Seg

  subroutine ExtractGabElmP15Seg(AUXarray,Output,nPasses)
    implicit none
    integer,intent(in) :: nPasses
    real(realk),intent(in) :: AUXarray(   15,   15,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,i
    real(realk) :: TMP(   15)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,   15
      TMP(i) = ABS(AUXarray(i,i,iPassQ))
     enddo
     Output(iPassQ) = SQRT(MAXVAL(TMP))
    enddo
  end subroutine ExtractGabElmP15Seg

  subroutine ExtractGabElmP25Seg(AUXarray,Output,nPasses)
    implicit none
    integer,intent(in) :: nPasses
    real(realk),intent(in) :: AUXarray(   25,   25,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,i
    real(realk) :: TMP(   25)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,   25
      TMP(i) = ABS(AUXarray(i,i,iPassQ))
     enddo
     Output(iPassQ) = SQRT(MAXVAL(TMP))
    enddo
  end subroutine ExtractGabElmP25Seg
END MODULE IchorEriGabintegralOBSGeneralModSeg
