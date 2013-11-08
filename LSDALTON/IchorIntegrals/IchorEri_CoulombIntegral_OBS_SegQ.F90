MODULE IchorEriCoulombintegralOBSGeneralModSegQ
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for a General Contracted LHS Segmented contracted RHS and Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_OBS_VERTICALRECURRENCEMODASegQ
use AGC_OBS_VERTICALRECURRENCEMODBSegQ
use AGC_OBS_VERTICALRECURRENCEMODDSegQ
use AGC_OBS_VERTICALRECURRENCEMODCSegQ
use AGC_OBS_TRANSFERRECURRENCEMODAtoCSegQ
use AGC_OBS_TRANSFERRECURRENCEMODAtoDSegQ
use AGC_OBS_TRANSFERRECURRENCEMODBtoCSegQ
use AGC_OBS_TRANSFERRECURRENCEMODBtoDSegQ
use AGC_OBS_TRANSFERRECURRENCEMODCtoASegQ
use AGC_OBS_TRANSFERRECURRENCEMODDtoASegQ
use AGC_OBS_TRANSFERRECURRENCEMODCtoBSegQ
use AGC_OBS_TRANSFERRECURRENCEMODDtoBSegQ
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_OBS_SegQ,IchorCoulombIntegral_OBS_general_sizeSegQ  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&
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
        call VerticalRecurrenceSegQ0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ1(TMParray2,CDAB     ,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(25,nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(25,nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceSegQ1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrenceSegQ2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(1,nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceSegQ1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrenceSegQ2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrenceSegQ3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(1,nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrenceSegQ2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(1,nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrenceSegQ3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(1,nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrenceSegQ4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(1,nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(9,nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(9,nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(9,nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(9,nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(9,nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionSegQ20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(25,nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
         call PrimitiveContractionSegQ700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(25,nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_SegQ',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_SegQ
  
  
  subroutine IchorCoulombIntegral_OBS_general_sizeSegQ(TMParray1maxsize,&
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
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,1*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,1*nContP)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,4*nContP)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,10*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,6*nContP)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,4*nContP)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,10*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,6*nContP)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,35*nContP)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,3*nContP)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContP)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,5*nContP)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContP)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,3*nContP)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContP)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,54*nContP)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,54*nContP)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,15*nContP)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContP)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContP)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContP)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContP)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,525*nContP)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,5*nContP)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContP)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,15*nContP)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContP)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContP)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContP)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContP)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,525*nContP)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,25*nContP)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContP)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,250*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContP)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,500*nContP)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,250*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,500*nContP)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,165*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,875*nContP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_general_sizeSegQ
  subroutine PrimitiveContractionSegQ1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    !Due to Q being segmented the Q contraction have already been done and we need to 
    !go from nPrimP to nContP
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: tmp,tmpA(nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       tmp = 0.0E0_realk
       do iPrimA=1,nPrimA
        tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iPrimA,iPrimB,iPassQ)
       enddo
       tmpA(iPrimB) = tmp
      enddo
      do iContB=1,nContB
       do iPrimB=1,nPrimB
        tmp = tmp + BCC(iPrimB,iContB)*tmpA(iPrimB)
       enddo
       AUXarrayCont(iContA,iContB,iPassQ) = tmp
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ1


  subroutine PrimitiveContractionSegQ4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(    4,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(    4)
    real(realk) :: TMPA(    4,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,    4
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,    4
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,    4
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,    4
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,    4
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,    4
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ4


  subroutine PrimitiveContractionSegQ10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   10,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(   10)
    real(realk) :: TMPA(   10,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   10
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,   10
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,   10
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   10
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,   10
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,   10
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ10


  subroutine PrimitiveContractionSegQ20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   20,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(   20)
    real(realk) :: TMPA(   20,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   20
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,   20
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,   20
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   20
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,   20
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,   20
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ20


  subroutine PrimitiveContractionSegQ35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   35,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(   35)
    real(realk) :: TMPA(   35,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   35
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,   35
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,   35
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   35
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,   35
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,   35
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ35


  subroutine PrimitiveContractionSegQ16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(   16)
    real(realk) :: TMPA(   16,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   16
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,   16
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,   16
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   16
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,   16
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,   16
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ16


  subroutine PrimitiveContractionSegQ40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   40,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(   40)
    real(realk) :: TMPA(   40,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   40
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,   40
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,   40
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   40
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,   40
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,   40
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ40


  subroutine PrimitiveContractionSegQ80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   80,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(   80)
    real(realk) :: TMPA(   80,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   80
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,   80
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,   80
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   80
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,   80
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,   80
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ80


  subroutine PrimitiveContractionSegQ140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  140,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(  140)
    real(realk) :: TMPA(  140,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  140
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,  140
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,  140
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  140
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,  140
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,  140
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ140


  subroutine PrimitiveContractionSegQ100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(  100)
    real(realk) :: TMPA(  100,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  100
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,  100
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,  100
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  100
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,  100
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,  100
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ100


  subroutine PrimitiveContractionSegQ200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  200,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(  200)
    real(realk) :: TMPA(  200,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  200
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,  200
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,  200
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  200
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,  200
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,  200
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ200


  subroutine PrimitiveContractionSegQ350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  350,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(  350)
    real(realk) :: TMPA(  350,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  350
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,  350
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,  350
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  350
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,  350
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,  350
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ350


  subroutine PrimitiveContractionSegQ400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(  400)
    real(realk) :: TMPA(  400,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  400
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,  400
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,  400
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  400
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,  400
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,  400
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ400


  subroutine PrimitiveContractionSegQ700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  700,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray(  700)
    real(realk) :: TMPA(  700,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  700
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1,  700
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1,  700
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  700
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1,  700
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1,  700
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ700


  subroutine PrimitiveContractionSegQ1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: ACCTMP,BCCTMP,TMPArray( 1225)
    real(realk) :: TMPA( 1225,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1, 1225
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimA=1,nPrimA
        ACCTMP = ACC(iPrimA,iContA)
        do iTUV=1, 1225
         tmpArray(iTUV) = tmpArray(iTUV) + ACCTMP*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
       enddo
       do iTUV=1, 1225
        tmpA(iTUV,iPrimB) = tmpArray(iTUV)
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1, 1225
        TMPArray(iTUV) = 0.0E0_realk
       enddo
       do iPrimB=1,nPrimB
        BCCTMP = BCC(iPrimB,iContB)
        do iTUV=1, 1225
         tmpArray(iTUV) = tmpArray(iTUV) + BCC(iPrimB,iContB)*tmpA(iTUV,iPrimB)
        enddo
       enddo
       do iTUV=1, 1225
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmpArray(iTUV)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ1225

END MODULE IchorEriCoulombintegralOBSGeneralModSegQ
