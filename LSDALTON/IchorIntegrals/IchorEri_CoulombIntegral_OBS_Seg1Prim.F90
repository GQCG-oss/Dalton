MODULE IchorEriCoulombintegralOBSGeneralModSeg1Prim
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset containing a single primitive
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_OBS_VERTICALRECURRENCEMODASeg1Prim
use AGC_OBS_VERTICALRECURRENCEMODBSeg1Prim
use AGC_OBS_VERTICALRECURRENCEMODDSeg1Prim
use AGC_OBS_VERTICALRECURRENCEMODCSeg1Prim
use AGC_OBS_VERTICALRECURRENCEMODA
use AGC_OBS_VERTICALRECURRENCEMODB
use AGC_OBS_VERTICALRECURRENCEMODC
use AGC_OBS_VERTICALRECURRENCEMODD
use AGC_OBS_TRANSFERRECURRENCEMODAtoCSeg1Prim
use AGC_OBS_TRANSFERRECURRENCEMODAtoDSeg1Prim
use AGC_OBS_TRANSFERRECURRENCEMODBtoCSeg1Prim
use AGC_OBS_TRANSFERRECURRENCEMODBtoDSeg1Prim
use AGC_OBS_TRANSFERRECURRENCEMODCtoASeg1Prim
use AGC_OBS_TRANSFERRECURRENCEMODDtoASeg1Prim
use AGC_OBS_TRANSFERRECURRENCEMODCtoBSeg1Prim
use AGC_OBS_TRANSFERRECURRENCEMODDtoBSeg1Prim
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_OBS_Seg1Prim,IchorCoulombIntegral_OBS_general_sizeSeg1Prim  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&
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
        call VerticalRecurrenceSeg1Prim0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(25,nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(25,nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceSeg1Prim1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrenceSeg1Prim2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(1,nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceSeg1Prim1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrenceSeg1Prim2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrenceSeg1Prim3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(1,nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrenceSeg1Prim2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(1,nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrenceSeg1Prim3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(1,nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrenceSeg1Prim4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(1,nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(9,nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(9,nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(9,nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(9,nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(9,nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(1,nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nPasses,35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(35,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(1,nPasses,25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(25,nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nPasses,20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(1,nPasses,25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(25,nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_Seg1Prim',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_Seg1Prim
  
  subroutine IchorCoulombIntegral_OBS_general_sizeSeg1Prim(TMParray1maxsize,&
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
       TMParray2maxSize = MAX(TMParray2maxSize,1)
       TMParray1maxSize = MAX(TMParray1maxSize,1)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
       TMParray1maxSize = MAX(TMParray1maxSize,4)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
       TMParray1maxSize = MAX(TMParray1maxSize,4)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35)
       TMParray1maxSize = MAX(TMParray1maxSize,35)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
       TMParray1maxSize = MAX(TMParray1maxSize,3)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16)
       TMParray2maxSize = MAX(TMParray2maxSize,12)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
       TMParray1maxSize = MAX(TMParray1maxSize,18)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16)
       TMParray2maxSize = MAX(TMParray2maxSize,12)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
       TMParray1maxSize = MAX(TMParray1maxSize,18)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140)
       TMParray2maxSize = MAX(TMParray2maxSize,105)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10)
       TMParray1maxSize = MAX(TMParray1maxSize,6)
       TMParray2maxSize = MAX(TMParray2maxSize,5)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
       TMParray2maxSize = MAX(TMParray2maxSize,24)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
       TMParray1maxSize = MAX(TMParray1maxSize,50)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
       TMParray2maxSize = MAX(TMParray2maxSize,24)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
       TMParray2maxSize = MAX(TMParray2maxSize,120)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
       TMParray1maxSize = MAX(TMParray1maxSize,50)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
       TMParray2maxSize = MAX(TMParray2maxSize,120)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350)
       TMParray2maxSize = MAX(TMParray2maxSize,210)
       TMParray1maxSize = MAX(TMParray1maxSize,175)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
       TMParray1maxSize = MAX(TMParray1maxSize,3)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16)
       TMParray2maxSize = MAX(TMParray2maxSize,12)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
       TMParray1maxSize = MAX(TMParray1maxSize,18)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16)
       TMParray2maxSize = MAX(TMParray2maxSize,12)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
       TMParray1maxSize = MAX(TMParray1maxSize,18)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140)
       TMParray2maxSize = MAX(TMParray2maxSize,105)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20)
       TMParray1maxSize = MAX(TMParray1maxSize,15)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700)
       TMParray2maxSize = MAX(TMParray2maxSize,525)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10)
       TMParray1maxSize = MAX(TMParray1maxSize,6)
       TMParray2maxSize = MAX(TMParray2maxSize,5)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
       TMParray2maxSize = MAX(TMParray2maxSize,24)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
       TMParray1maxSize = MAX(TMParray1maxSize,50)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40)
       TMParray2maxSize = MAX(TMParray2maxSize,24)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
       TMParray2maxSize = MAX(TMParray2maxSize,120)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
       TMParray1maxSize = MAX(TMParray1maxSize,50)
       TMParray2maxSize = MAX(TMParray2maxSize,30)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
       TMParray2maxSize = MAX(TMParray2maxSize,120)
       TMParray1maxSize = MAX(TMParray1maxSize,100)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350)
       TMParray2maxSize = MAX(TMParray2maxSize,210)
       TMParray1maxSize = MAX(TMParray1maxSize,175)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20)
       TMParray1maxSize = MAX(TMParray1maxSize,15)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700)
       TMParray2maxSize = MAX(TMParray2maxSize,525)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35)
       TMParray1maxSize = MAX(TMParray1maxSize,25)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350)
       TMParray2maxSize = MAX(TMParray2maxSize,250)
       TMParray1maxSize = MAX(TMParray1maxSize,150)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700)
       TMParray2maxSize = MAX(TMParray2maxSize,500)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350)
       TMParray2maxSize = MAX(TMParray2maxSize,250)
       TMParray1maxSize = MAX(TMParray1maxSize,150)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700)
       TMParray2maxSize = MAX(TMParray2maxSize,500)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,165*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225)
       TMParray2maxSize = MAX(TMParray2maxSize,875)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_general_sizeSeg1Prim
  
END MODULE IchorEriCoulombintegralOBSGeneralModSeg1Prim
