module SPIchorEriCoulombintegralGPUOBSGeneralModSegP2
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted LHS and a General Contracted RHS Basisset 
use IchorprecisionMod
use IchorCommonMod
use IchorMemory
use SPAGC_GPU_OBS_BUILDRJ000ModGen
use SPAGC_GPU_OBS_BUILDRJ000ModSeg1Prim
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODASegP
use SPAGC_GPU_PrimContractSegPMod
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODBSegP
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODDSegP
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODCSegP
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODAGen
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODBGen
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODCGen
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODDGen
use SPAGC_GPU_OBS_TRMODAtoCSegP1
use SPAGC_GPU_OBS_TRMODAtoCSegP2
use SPAGC_GPU_OBS_TRMODAtoCSegP3
use SPAGC_GPU_OBS_TRMODAtoDSegP1
use SPAGC_GPU_OBS_TRMODAtoDSegP2
use SPAGC_GPU_OBS_TRMODBtoCSegP1
use SPAGC_GPU_OBS_TRMODBtoDSegP1
use SPAGC_GPU_OBS_TRMODCtoASegP
use SPAGC_GPU_OBS_TRMODDtoASegP
use SPAGC_GPU_OBS_TRMODCtoBSegP
use SPAGC_GPU_OBS_TRMODDtoBSegP
use SPAGC_GPU_OBS_HorizontalRecurrenceLHSModAtoB
use SPAGC_GPU_OBS_HorizontalRecurrenceLHSModBtoA
use SPAGC_GPU_OBS_HorizontalRecurrenceRHSModCtoD
use SPAGC_GPU_OBS_HorizontalRecurrenceRHSModDtoC
use SPAGC_GPU_OBS_Sphcontract1Mod
use SPAGC_GPU_OBS_Sphcontract2Mod
  
private   
public :: SPICI_GPU_OBS_SegP2
  
CONTAINS
  
  
  subroutine SPICI_GPU_OBS_SegP2(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    real(reals),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(reals),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(reals),intent(in) :: pcent(3*nPrimP*nAtomsA*nAtomsB)           !qcent(3,nPrimP)
    real(reals),intent(in) :: qcent(3*nPrimQ)           !qcent(3,nPrimQ)
    real(reals),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP*nAtomsA*nAtomsB)
    real(reals),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(reals),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(reals),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(reals) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(reals) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    integer,intent(in) :: localintsmaxsize
    real(reals),intent(inout) :: LOCALINTS(localintsmaxsize)
    real(reals),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(reals),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(reals),intent(in) :: Qdistance12(3) !Ccenter-Dcenter
    !Qdistance12(3)
    real(reals),intent(in) :: Pdistance12(3*nAtomsA*nAtomsB)           !Acenter-Bcenter 
    real(reals),intent(in) :: Acenter(3,nAtomsA),Bcenter(3,nAtomsB),Ccenter(3),Dcenter(3)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(reals),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
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
        call SPVerticalRecurrenceGPUSegP1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,4,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUSegP2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call SPVerticalRecurrenceGPUSegP1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,4,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUSegP2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUSegP3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,20,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUSegP2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUSegP3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,20,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUSegP4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,35,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,35,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call SPVerticalRecurrenceGPUSegP1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,1,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,1,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q2DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,35,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,35,iASync)
        call SPHorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUSegP2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,1,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,1,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
        call SPHorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA0(35,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q2DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,35,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,35,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUSegP3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
         call SPPrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,1,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,1,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q3BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q3BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,35,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,35,iASync)
        call SPHorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA1(35,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP4Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,4,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,4,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP4Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,10,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,10,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC0(25,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call SPBuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP4Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
         call SPPrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,20,iASync)
         call SPPrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,20,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC1(25,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE DEFAULT
        call ichorquit('Unknown Case in ICI_GPU_OBS_SegP',-1)
    END SELECT
  end subroutine SPICI_GPU_OBS_SegP2
  
END module SPIchorEriCoulombintegralGPUOBSGeneralModSegP2
