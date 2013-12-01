MODULE IchorEriCoulombintegralOBSGeneralModSegP
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted LHS and a General Contracted RHS Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_OBS_VERTICALRECURRENCEMODASegP
use AGC_OBS_VERTICALRECURRENCEMODBSegP
use AGC_OBS_VERTICALRECURRENCEMODDSegP
use AGC_OBS_VERTICALRECURRENCEMODCSegP
use AGC_OBS_VERTICALRECURRENCEMODA
use AGC_OBS_VERTICALRECURRENCEMODB
use AGC_OBS_VERTICALRECURRENCEMODC
use AGC_OBS_VERTICALRECURRENCEMODD
use AGC_OBS_TRANSFERRECURRENCEMODAtoCSegP
use AGC_OBS_TRANSFERRECURRENCEMODAtoDSegP
use AGC_OBS_TRANSFERRECURRENCEMODBtoCSegP
use AGC_OBS_TRANSFERRECURRENCEMODBtoDSegP
use AGC_OBS_TRANSFERRECURRENCEMODCtoASegP
use AGC_OBS_TRANSFERRECURRENCEMODDtoASegP
use AGC_OBS_TRANSFERRECURRENCEMODCtoBSegP
use AGC_OBS_TRANSFERRECURRENCEMODDtoBSegP
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_OBS_SegP,IchorCoulombIntegral_OBS_general_sizeSegP  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsC,nAtomsD
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3),Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
!   Local variables 
    integer :: AngmomPQ,AngmomP,AngmomQ,I,J,nContQP,la,lb,lc,ld,nsize,angmomid
    
    !Setup combined Angmom info
    AngmomP = AngmomA+AngmomB
    AngmomQ = AngmomC+AngmomD
    AngmomPQ  = AngmomP + AngmomQ
!    nTUV = (AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6
!    nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
!    nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
!    nTUVC = (AngmomC+1)*(AngmomC+2)*(AngmomC+3)/6
!    nTUVD = (AngmomD+1)*(AngmomD+2)*(AngmomD+3)/6
!    nlmA = 2*AngmomA+1
!    nlmB = 2*AngmomB+1
!    nlmC = 2*AngmomC+1
!    nlmD = 2*AngmomD+1
    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegP0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP1(TMParray2,CDAB     ,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegP1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(1,nContQ*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(35,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegP3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(1,nContQ*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(35,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegP4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(1,nContQ*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(25,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(25,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nContQ*nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceSegP1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrenceSegP2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(1,nContQ*nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceSegP1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrenceSegP2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrenceSegP3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(1,nContQ*nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrenceSegP2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(1,nContQ*nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrenceSegP3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(1,nContQ*nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrenceSegP4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(1,nContQ*nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegP1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegP2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(1,nContQ*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(35,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(9,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(9,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(9,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(9,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(9,nContQ*nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegP3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegP20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(1,nContQ*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(35,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(25,nContQ*nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegP700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(25,nContQ*nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_SegP',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_SegP
  
  
  subroutine IchorCoulombIntegral_OBS_general_sizeSegP(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         &nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP
    ! local variables
    integer :: AngmomID
    
    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimQ)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQ)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQ)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQ)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQ)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQ)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQ)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQ)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQ)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQ)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQ)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQ)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,315*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,324*nContQ)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQ)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQ)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQ)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQ)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQ)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQ)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQ)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,165*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContQ)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_general_sizeSegP

  subroutine PrimitiveContractionSegP1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    !Due to P being segmented the P contraction have already been done and we need to 
    !go from nPrimQ to nContQ
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: tmp,tmpC(nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       tmp = 0.0E0_realk
       do iPrimC=1,nPrimC
        tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iPrimD,iPassQ)
       enddo
       tmpC(iPrimD) = tmp
      enddo
      do iContD=1,nContD
       tmp = 0.0E0_realk
       do iPrimD=1,nPrimD
        tmp = tmp + DCC(iPrimD,iContD)*tmpC(iPrimD)
       enddo
       AUXarrayCont(iContC,iContD,iPassQ) = tmp
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP1


  subroutine PrimitiveContractionSegP4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(    4,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(    4)
    real(realk) :: TMPC(    4,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,    4
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,    4
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,    4
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,    4
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,    4
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,    4
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP4

  subroutine PrimitiveContractionSegP10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   10,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(   10)
    real(realk) :: TMPC(   10,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   10
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,   10
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,   10
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   10
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,   10
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,   10
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP10

  subroutine PrimitiveContractionSegP20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   20,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(   20)
    real(realk) :: TMPC(   20,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   20
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,   20
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,   20
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   20
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,   20
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,   20
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP20

  subroutine PrimitiveContractionSegP35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   35,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(   35)
    real(realk) :: TMPC(   35,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   35
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,   35
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,   35
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   35
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,   35
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,   35
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP35

  subroutine PrimitiveContractionSegP16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   16,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(   16)
    real(realk) :: TMPC(   16,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   16
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,   16
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,   16
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   16
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,   16
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,   16
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP16

  subroutine PrimitiveContractionSegP40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   40,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(   40)
    real(realk) :: TMPC(   40,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   40
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,   40
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,   40
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   40
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,   40
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,   40
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP40

  subroutine PrimitiveContractionSegP80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   80,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(   80)
    real(realk) :: TMPC(   80,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   80
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,   80
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,   80
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   80
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,   80
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,   80
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP80

  subroutine PrimitiveContractionSegP140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  140,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(  140)
    real(realk) :: TMPC(  140,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  140
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,  140
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,  140
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  140
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,  140
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,  140
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP140

  subroutine PrimitiveContractionSegP100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  100,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(  100)
    real(realk) :: TMPC(  100,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  100
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,  100
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,  100
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  100
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,  100
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,  100
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP100

  subroutine PrimitiveContractionSegP200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  200,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(  200)
    real(realk) :: TMPC(  200,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  200
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,  200
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,  200
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  200
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,  200
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,  200
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP200

  subroutine PrimitiveContractionSegP350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  350,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(  350)
    real(realk) :: TMPC(  350,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  350
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,  350
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,  350
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  350
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,  350
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,  350
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP350

  subroutine PrimitiveContractionSegP400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  400,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(  400)
    real(realk) :: TMPC(  400,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  400
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,  400
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,  400
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  400
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,  400
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,  400
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP400

  subroutine PrimitiveContractionSegP700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  700,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray(  700)
    real(realk) :: TMPC(  700,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  700
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1,  700
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1,  700
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  700
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1,  700
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1,  700
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP700

  subroutine PrimitiveContractionSegP1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMPArray( 1225)
    real(realk) :: TMPC( 1225,nPrimD),CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1, 1225
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimC=1,nPrimC
        CCCTMP = CCC(iPrimC,iContC)
        do iTUV=1, 1225
         tmparray(iTUV) = tmparray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
       enddo
       do iTUV=1, 1225
        tmpC(iTUV,iPrimD) = tmparray(iTUV)
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1, 1225
        tmparray(iTUV) = 0.0E0_realk
       enddo
       do iPrimD=1,nPrimD
        DCCTMP = DCC(iPrimD,iContD)
        do iTUV=1, 1225
         tmparray(iTUV) = tmparray(iTUV) + DCCTMP*tmpC(iTUV,iPrimD)
        enddo
       enddo
       do iTUV=1, 1225
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmparray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP1225
END MODULE IchorEriCoulombintegralOBSGeneralModSegP
