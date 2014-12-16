MODULE IchorEriCoulombintegralGPUOBSGeneralModSeg1Prim
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset containing a single primitive
use IchorprecisionMod
use IchorCommonMod
use IchorMemory
use AGC_GPU_OBS_BUILDRJ000ModGen
use AGC_GPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_GPU_OBS_VERTICALRECURRENCEMODASeg1Prim
use AGC_GPU_OBS_VERTICALRECURRENCEMODBSeg1Prim
use AGC_GPU_OBS_VERTICALRECURRENCEMODDSeg1Prim
use AGC_GPU_OBS_VERTICALRECURRENCEMODCSeg1Prim
use AGC_GPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_GPU_OBS_TRMODAtoCSeg1Prim1
use AGC_GPU_OBS_TRMODAtoCSeg1Prim2
use AGC_GPU_OBS_TRMODAtoDSeg1Prim1
use AGC_GPU_OBS_TRMODAtoDSeg1Prim2
use AGC_GPU_OBS_TRMODBtoCSeg1Prim1
use AGC_GPU_OBS_TRMODBtoDSeg1Prim1
use AGC_GPU_OBS_TRMODCtoASeg1Prim
use AGC_GPU_OBS_TRMODDtoASeg1Prim
use AGC_GPU_OBS_TRMODCtoBSeg1Prim
use AGC_GPU_OBS_TRMODDtoBSeg1Prim
use AGC_GPU_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_GPU_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_GPU_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_GPU_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_GPU_OBS_Sphcontract1Mod
use AGC_GPU_OBS_Sphcontract2Mod
  
private   
public :: ICI_GPU_OBS_Seg1Prim
  
CONTAINS
  
  
  subroutine ICI_GPU_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & IatomAPass,iatomBPass,iASync)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsA,nAtomsB
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD
    integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
    integer,intent(in) :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP*nAtomsA*nAtomsB)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ)           !qcent(3,nPrimQ)
    real(realk),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP*nAtomsA*nAtomsB)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    integer,intent(in) :: localintsmaxsize
    real(realk),intent(inout) :: LOCALINTS(localintsmaxsize)
    real(realk),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3) !Ccenter-Dcenter
    !Qdistance12(3)
    real(realk),intent(in) :: Pdistance12(3*nAtomsA*nAtomsB)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3,nAtomsA),Bcenter(3,nAtomsB),Ccenter(3),Dcenter(3)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
    integer(kind=acckind),intent(in) :: iASync
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
    IF(UseGeneralCode) AngmomID = AngmomID + 10000 !force to use general code
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceGPUSeg1Prim0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,LOCALINTS(1),iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceGPUSeg1Prim1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2(1),iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call BuildRJ000GPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP1Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,3,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call BuildRJ000GPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP1Q2CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,3,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call BuildRJ000GPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call BuildRJ000GPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP2Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,9,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call BuildRJ000GPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,9,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call BuildRJ000GPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & TMParray2(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(1,nPasses,TMParray2(1),&
            & LOCALINTS(1),iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call BuildRJ000GPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP2Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(4,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,5,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call BuildRJ000GPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,5,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call BuildRJ000GPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,5,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(5,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call BuildRJ000GPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP2Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(20,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,5,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(5,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call BuildRJ000GPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP2Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(35,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,5,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(5,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call BuildRJ000GPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & TMParray2(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(1,nPasses,TMParray2(1),&
            & LOCALINTS(1),iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call BuildRJ000GPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP3Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(4,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,15,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call BuildRJ000GPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP3Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,15,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call BuildRJ000GPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP3Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,15,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(15,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call BuildRJ000GPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP3Q3AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(20,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,15,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(15,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call BuildRJ000GPUSeg1Prim7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP3Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(35,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,15,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(15,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call BuildRJ000GPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & TMParray2(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(1,nPasses,TMParray2(1),&
            & LOCALINTS(1),iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call BuildRJ000GPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP4Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(4,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,25,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call BuildRJ000GPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP4Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,25,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call BuildRJ000GPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP4Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,25,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(25,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call BuildRJ000GPUSeg1Prim7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP4Q3AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(20,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,25,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(25,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call BuildRJ000GPUSeg1Prim8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1),iASync)
        call VerticalRecurrenceGPUSeg1Prim8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1),iASync)
        call TransferRecurrenceGPUP4Q4AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1),iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(35,nPasses,TMParray1(1),&
            & TMParray2(1),iASync)
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,25,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(25,nPasses,TMParray1(1),&
            & LOCALINTS(1),iASync)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in ICI_GPU_OBS_Seg1Prim',-1)
    END SELECT
  end subroutine ICI_GPU_OBS_Seg1Prim
  
END MODULE IchorEriCoulombintegralGPUOBSGeneralModSeg1Prim
