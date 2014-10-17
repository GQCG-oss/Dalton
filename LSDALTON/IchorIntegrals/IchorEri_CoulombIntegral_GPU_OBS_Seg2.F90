MODULE IchorEriCoulombintegralGPUOBSGeneralModSeg2
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_GPU_OBS_BUILDRJ000ModGen
use AGC_GPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_GPU_OBS_VERTICALRECURRENCEMODASeg
use AGC_GPU_OBS_VERTICALRECURRENCEMODBSeg
use AGC_GPU_OBS_VERTICALRECURRENCEMODDSeg
use AGC_GPU_OBS_VERTICALRECURRENCEMODCSeg
use AGC_GPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_GPU_OBS_TRMODAtoCSeg1
use AGC_GPU_OBS_TRMODAtoCSeg2
use AGC_GPU_OBS_TRMODAtoCSeg3
use AGC_GPU_OBS_TRMODAtoDSeg1
use AGC_GPU_OBS_TRMODAtoDSeg2
use AGC_GPU_OBS_TRMODBtoCSeg1
use AGC_GPU_OBS_TRMODBtoDSeg1
use AGC_GPU_OBS_TRMODCtoASeg
use AGC_GPU_OBS_TRMODDtoASeg
use AGC_GPU_OBS_TRMODCtoBSeg
use AGC_GPU_OBS_TRMODDtoBSeg
use AGC_GPU_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_GPU_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_GPU_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_GPU_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_GPU_OBS_Sphcontract1Mod
use AGC_GPU_OBS_Sphcontract2Mod
  
private   
public :: ICI_GPU_OBS_Seg2
  
CONTAINS
  
  
  subroutine ICI_GPU_OBS_Seg2(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceGPUSeg1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUSeg2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(1,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceGPUSeg1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUSeg2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,1,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUSeg3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(1,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUSeg2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(1,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUSeg3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(1,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUSeg4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(1,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceGPUSeg1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q1BtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q2DtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q1BtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q2CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q3DtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q2CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q3CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q4CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUSeg2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(1,nPasses,TMParray2,&
            & LOCALINTS,iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q1BtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(4,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q2BtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(5,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q1BtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(4,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q2BtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q3DtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(20,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(5,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q2BtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q3CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(20,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q4CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(35,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q1AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q2DtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q3DtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q2CtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q3CtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP1Q4CtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q1AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q2AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(9,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q3DtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(9,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q2AtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(9,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q3CtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(9,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q4CtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(9,nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUSeg3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(1,nPasses,TMParray2,&
            & LOCALINTS,iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q1BtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(4,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q2BtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(15,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q1BtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(4,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q2BtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q3BtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(20,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(15,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q2BtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q3BtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(20,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q4CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(35,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q1AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q2AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(5,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP2Q3DtoASeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(5,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q1AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q2AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(15,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP3Q3AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(15,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP4Q1AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(1,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP4Q2AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(1,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(25,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call VerticalRecurrenceGPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call TransferRecurrenceGPUP4Q3AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        !Primitive Contraction have already been done
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2,iASync)
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(1,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(25,nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in ICI_GPU_OBS_Seg',-1)
    END SELECT
  end subroutine ICI_GPU_OBS_Seg2
  
END MODULE IchorEriCoulombintegralGPUOBSGeneralModSeg2
