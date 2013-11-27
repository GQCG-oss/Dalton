MODULE IchorEriCoulombintegralOBSGeneralModGen
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_OBS_VERTICALRECURRENCEMODA
use AGC_OBS_VERTICALRECURRENCEMODB
use AGC_OBS_VERTICALRECURRENCEMODD
use AGC_OBS_VERTICALRECURRENCEMODC
use AGC_OBS_TRANSFERRECURRENCEMODAtoCGen
use AGC_OBS_TRANSFERRECURRENCEMODAtoDGen
use AGC_OBS_TRANSFERRECURRENCEMODBtoCGen
use AGC_OBS_TRANSFERRECURRENCEMODBtoDGen
use AGC_OBS_TRANSFERRECURRENCEMODCtoAGen
use AGC_OBS_TRANSFERRECURRENCEMODDtoAGen
use AGC_OBS_TRANSFERRECURRENCEMODCtoBGen
use AGC_OBS_TRANSFERRECURRENCEMODDtoBGen
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_OBS_Gen,IchorCoulombIntegral_OBS_general_sizeGen  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
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
        call VerticalRecurrence0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen1(TMParray2,CDAB     ,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrence1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(1,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrence1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call PrimitiveContractionGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
         call PrimitiveContractionGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(25,nContQP*nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_Gen',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_Gen
  
  
  subroutine IchorCoulombIntegral_OBS_general_sizeGen(TMParray1maxsize,&
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
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimQP)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQP)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQP)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQP)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQP)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQP)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQP)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQP)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQP)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQP)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQP)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQP)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,315*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,324*nContQP)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQP)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQP)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQP)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQP)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQP)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,165*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContQP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_general_sizeGen

  subroutine PrimitiveContractionGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: TMP
    real(realk) :: tmpArray1(nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(nPrimA,nPrimB)
    real(realk) :: tmpArray3(nPrimB)
    !Scaling p**4*c*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nPassQ
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         TMP = 0.0E0_realk
         do iPrimC=1,nPrimC
          TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
         enddo
         tmpArray1(iPrimD,iPrimA,iPrimB) = TMP
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         TMP = 0.0E0_realk
         do iPrimD=1,nPrimD
          TMP = TMP + DCC(iPrimD,iContD)*tmpArray1(iPrimD,iPrimA,iPrimB)
         enddo
         tmpArray2(iPrimA,iPrimB) = TMP
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          TMP = TMP + ACC(iPrimA,iContA)*tmpArray2(iPrimA,iPrimB)
         enddo
         tmpArray3(iPrimB) = TMP
        enddo
        do iContB=1,nContB
         TMP = 0.0E0_realk
         do iPrimB=1,nPrimB
          TMP = TMP + BCC(iPrimB,iContB)*tmpArray3(iPrimB)
         enddo
         AUXarrayCont(iContC,iContD,iContA,iContB,iPassQ) = TMP
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen1

  subroutine PrimitiveContractionGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(    4,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(    4)
    real(realk) :: tmpArray1(    4,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(    4,nPrimA,nPrimB)
    real(realk) :: tmpArray3(    4,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,    4
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,    4
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,    4
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,    4
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,    4
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,    4
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,    4
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,    4
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,    4
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,    4
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,    4
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,    4
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen4

  subroutine PrimitiveContractionGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   10,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(   10)
    real(realk) :: tmpArray1(   10,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(   10,nPrimA,nPrimB)
    real(realk) :: tmpArray3(   10,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   10
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,   10
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,   10
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   10
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   10
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   10
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,   10
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   10
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   10
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,   10
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   10
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,   10
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen10

  subroutine PrimitiveContractionGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   20,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(   20)
    real(realk) :: tmpArray1(   20,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(   20,nPrimA,nPrimB)
    real(realk) :: tmpArray3(   20,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   20
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,   20
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,   20
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   20
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   20
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   20
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,   20
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   20
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   20
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,   20
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   20
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,   20
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen20

  subroutine PrimitiveContractionGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   35,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(   35)
    real(realk) :: tmpArray1(   35,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(   35,nPrimA,nPrimB)
    real(realk) :: tmpArray3(   35,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   35
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,   35
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,   35
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   35
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   35
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   35
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,   35
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   35
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   35
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,   35
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   35
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,   35
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen35

  subroutine PrimitiveContractionGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   16,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(   16)
    real(realk) :: tmpArray1(   16,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(   16,nPrimA,nPrimB)
    real(realk) :: tmpArray3(   16,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   16
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,   16
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,   16
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   16
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   16
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   16
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,   16
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   16
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   16
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,   16
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   16
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,   16
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen16

  subroutine PrimitiveContractionGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   40,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(   40)
    real(realk) :: tmpArray1(   40,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(   40,nPrimA,nPrimB)
    real(realk) :: tmpArray3(   40,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   40
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,   40
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,   40
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   40
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   40
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   40
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,   40
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   40
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   40
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,   40
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   40
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,   40
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen40

  subroutine PrimitiveContractionGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   80,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(   80)
    real(realk) :: tmpArray1(   80,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(   80,nPrimA,nPrimB)
    real(realk) :: tmpArray3(   80,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   80
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,   80
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,   80
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   80
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   80
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   80
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,   80
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   80
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,   80
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,   80
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,   80
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,   80
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen80

  subroutine PrimitiveContractionGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  140,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(  140)
    real(realk) :: tmpArray1(  140,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(  140,nPrimA,nPrimB)
    real(realk) :: tmpArray3(  140,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  140
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,  140
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,  140
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  140
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  140
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  140
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,  140
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  140
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  140
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,  140
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  140
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,  140
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen140

  subroutine PrimitiveContractionGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  100,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(  100)
    real(realk) :: tmpArray1(  100,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(  100,nPrimA,nPrimB)
    real(realk) :: tmpArray3(  100,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  100
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,  100
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,  100
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  100
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  100
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  100
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,  100
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  100
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  100
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,  100
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  100
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,  100
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen100

  subroutine PrimitiveContractionGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  200,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(  200)
    real(realk) :: tmpArray1(  200,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(  200,nPrimA,nPrimB)
    real(realk) :: tmpArray3(  200,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  200
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,  200
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,  200
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  200
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  200
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  200
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,  200
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  200
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  200
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,  200
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  200
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,  200
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen200

  subroutine PrimitiveContractionGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  350,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(  350)
    real(realk) :: tmpArray1(  350,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(  350,nPrimA,nPrimB)
    real(realk) :: tmpArray3(  350,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  350
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,  350
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,  350
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  350
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  350
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  350
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,  350
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  350
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  350
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,  350
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  350
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,  350
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen350

  subroutine PrimitiveContractionGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  400,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(  400)
    real(realk) :: tmpArray1(  400,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(  400,nPrimA,nPrimB)
    real(realk) :: tmpArray3(  400,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  400
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,  400
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,  400
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  400
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  400
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  400
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,  400
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  400
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  400
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,  400
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  400
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,  400
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen400

  subroutine PrimitiveContractionGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  700,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray(  700)
    real(realk) :: tmpArray1(  700,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2(  700,nPrimA,nPrimB)
    real(realk) :: tmpArray3(  700,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  700
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1,  700
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1,  700
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  700
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  700
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  700
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1,  700
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  700
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1,  700
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1,  700
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1,  700
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1,  700
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen700

  subroutine PrimitiveContractionGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMPArray( 1225)
    real(realk) :: tmpArray1( 1225,nPrimD,nPrimA,nPrimB)
    real(realk) :: tmpArray2( 1225,nPrimA,nPrimB)
    real(realk) :: tmpArray3( 1225,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
    do iPassQ = 1,nPasses
     iContQP = 0
     do iContC=1,nContC
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1, 1225
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimC=1,nPrimC
          CCCTMP = CCC(iPrimC,iContC)
          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iTUV=1, 1225
           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
          enddo
         enddo
         do iTUV=1, 1225
          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
      do iContD=1,nContD
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1, 1225
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimD=1,nPrimD
          DCCTMP = DCC(iPrimD,iContD)
          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iTUV=1, 1225
           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1, 1225
          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
       enddo
       do iContA=1,nContA
        do iPrimB=1,nPrimB
         do iTUV=1, 1225
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContA)
          !Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iTUV=1, 1225
           TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)
          enddo
         enddo
         do iTUV=1, 1225
          tmpArray3(iTUV,iPrimB) = TMPArray(iTUV)
         enddo
        enddo
        do iContB=1,nContB
         do iTUV=1, 1225
          TMPArray(iTUV) = 0.0E0_realk
         enddo
         do iPrimB=1,nPrimB
          BCCTMP = BCC(iPrimB,iContB)
          !Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iTUV=1, 1225
           TMPArray(iTUV) = TMPArray(iTUV) + BCCTMP*tmpArray3(iTUV,iPrimB)
          enddo
         enddo
         do iTUV=1, 1225
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassQ) = TMPArray(iTUV)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionGen1225
END MODULE IchorEriCoulombintegralOBSGeneralModGen
