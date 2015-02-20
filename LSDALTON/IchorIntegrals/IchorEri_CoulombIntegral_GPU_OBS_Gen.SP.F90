module SPIchorEriCoulombintegralGPUOBSGeneralModGen
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
use IchorCommonMod
use IchorMemory
use SPAGC_GPU_OBS_BUILDRJ000ModGen
use SPAGC_GPU_OBS_BUILDRJ000ModSeg1Prim
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODAGen
use SPAGC_GPU_PrimContractGenMod
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODBGen
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODDGen
use SPAGC_GPU_OBS_VERTICALRECURRENCEMODCGen
use SPAGC_GPU_OBS_TRMODAtoCGen1
use SPAGC_GPU_OBS_TRMODAtoCGen2
use SPAGC_GPU_OBS_TRMODAtoDGen1
use SPAGC_GPU_OBS_TRMODAtoDGen2
use SPAGC_GPU_OBS_TRMODBtoCGen1
use SPAGC_GPU_OBS_TRMODBtoDGen1
use SPAGC_GPU_OBS_TRMODCtoAGen
use SPAGC_GPU_OBS_TRMODDtoAGen
use SPAGC_GPU_OBS_TRMODCtoBGen
use SPAGC_GPU_OBS_TRMODDtoBGen
use SPAGC_GPU_OBS_HorizontalRecurrenceLHSModAtoB
use SPAGC_GPU_OBS_HorizontalRecurrenceLHSModBtoA
use SPAGC_GPU_OBS_HorizontalRecurrenceRHSModCtoD
use SPAGC_GPU_OBS_HorizontalRecurrenceRHSModDtoC
use SPAGC_GPU_OBS_Sphcontract1Mod
use SPAGC_GPU_OBS_Sphcontract2Mod
  
private   
public :: SPICI_GPU_OBS_Gen
  
CONTAINS
  
  
  subroutine SPICI_GPU_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
        call SPVerticalRecurrenceGPUGen0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,1,1,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,1,1,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,1,1,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,LOCALINTS,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,1,1,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call SPVerticalRecurrenceGPUGen1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,1,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,1,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,1,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,1,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,4,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,4,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,4,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,4,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,10,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,10,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,10,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,4,10,iASync)
        call SPHorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,1,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,1,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,1,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,1,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,4,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call SPBuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,1,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,1,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,1,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,1,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(1,nContQP*nPasses,TMParray2,&
            & LOCALINTS,iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,4,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,4,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(4,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,10,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(5,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,20,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,20,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,20,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,20,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(20,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(5,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,35,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,35,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,35,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,10,35,iASync)
        call SPHorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP2_maxAngA2(35,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(5,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call SPBuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,1,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,1,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,1,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,1,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(1,nContQP*nPasses,TMParray2,&
            & LOCALINTS,iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,4,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,4,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,4,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,4,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(4,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,10,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,10,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,10,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(15,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,20,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,20,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,20,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,20,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(20,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(15,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP3Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,35,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,35,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,35,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,20,35,iASync)
        call SPHorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP3_maxAngA2(35,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(15,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call SPBuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,1,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,1,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,1,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,1,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(1,nContQP*nPasses,TMParray2,&
            & LOCALINTS,iASync)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call SPBuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP4Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,4,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,4,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,4,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,4,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(4,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,10,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,10,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,10,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,10,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call SPBuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,10,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,10,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,10,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,10,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ2_maxAngC2(25,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call SPBuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP4Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,20,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,20,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,20,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,20,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(20,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ3_maxAngC2(25,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call SPBuildRJ000GPUGen8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
        call SPVerticalRecurrenceGPUGen8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        call SPTransferRecurrenceGPUP4Q4AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,35,iASync)
         call SPPrimitiveContractionDGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,35,iASync)
         call SPPrimitiveContractionAGPUGen(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,35,iASync)
         call SPPrimitiveContractionBGPUGen(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
               & nContC,nPrimD,nContD,35,35,iASync)
        call SPHorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS1_GPU_maxAngP4_maxAngA2(35,nContQP*nPasses,TMParray1,&
            & TMParray2,iASync)
        call SPHorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
        call SPSphericalContractOBS2_GPU_maxAngQ4_maxAngC2(25,nContQP*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE DEFAULT
        call ichorquit('Unknown Case in ICI_GPU_OBS_Gen',-1)
    END SELECT
  end subroutine SPICI_GPU_OBS_Gen
  
END module SPIchorEriCoulombintegralGPUOBSGeneralModGen
