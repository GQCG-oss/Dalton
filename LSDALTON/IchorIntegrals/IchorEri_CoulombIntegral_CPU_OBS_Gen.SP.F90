module SPIchorEriCoulombintegralCPUOBSGeneralModGen
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
use IchorCommonMod
use IchorMemory
use SPAGC_CPU_OBS_BUILDRJ000ModGen
use SPAGC_CPU_OBS_BUILDRJ000ModSeg1Prim
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use SPAGC_CPU_PrimContractGenMod
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use SPAGC_CPU_OBS_TRMODAtoCGen1
use SPAGC_CPU_OBS_TRMODAtoCGen2
use SPAGC_CPU_OBS_TRMODAtoDGen1
use SPAGC_CPU_OBS_TRMODAtoDGen2
use SPAGC_CPU_OBS_TRMODBtoCGen1
use SPAGC_CPU_OBS_TRMODBtoDGen1
use SPAGC_CPU_OBS_TRMODCtoAGen
use SPAGC_CPU_OBS_TRMODDtoAGen
use SPAGC_CPU_OBS_TRMODCtoBGen
use SPAGC_CPU_OBS_TRMODDtoBGen
use SPAGC_CPU_OBS_HorizontalRecurrenceLHSModAtoB
use SPAGC_CPU_OBS_HorizontalRecurrenceLHSModBtoA
use SPAGC_CPU_OBS_HorizontalRecurrenceRHSModCtoD
use SPAGC_CPU_OBS_HorizontalRecurrenceRHSModDtoC
use SPAGC_CPU_OBS_Sphcontract1Mod
use SPAGC_CPU_OBS_Sphcontract2Mod
  
private   
public :: SPICI_CPU_OBS_Gen
  
CONTAINS
  
  
  subroutine SPICI_CPU_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
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
       & IatomAPass,iatomBPass)
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
        call SPVerticalRecurrenceCPUGen0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen1(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen1(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen1(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen1(TMParray1,LOCALINTS,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call SPVerticalRecurrenceCPUGen1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call SPBuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP1Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call SPBuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call SPBuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call SPBuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call SPBuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call SPHorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call SPBuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2(1,nContQP*nPasses,TMParray2,&
            & LOCALINTS)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call SPBuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call SPBuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call SPBuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call SPBuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call SPBuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2(35,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call SPBuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2(1,nContQP*nPasses,TMParray2,&
            & LOCALINTS)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call SPBuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP3Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call SPBuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call SPBuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call SPBuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP3Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call SPBuildRJ000CPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP3Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2(35,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call SPBuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2(1,nContQP*nPasses,TMParray2,&
            & LOCALINTS)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call SPBuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP4Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call SPBuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call SPBuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ2_maxAngC2(25,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call SPBuildRJ000CPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP4Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ3_maxAngC2(25,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call SPBuildRJ000CPUGen8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
        call SPVerticalRecurrenceCPUGen8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        call SPTransferRecurrenceCPUP4Q4AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
         call SPPrimitiveContractionCCPUGen1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionDCPUGen1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionACPUGen1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call SPPrimitiveContractionBCPUGen1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call SPHorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2(35,nContQP*nPasses,TMParray1,&
            & TMParray2)
        call SPHorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
        call SPSphericalContractOBS2_CPU_maxAngQ4_maxAngC2(25,nContQP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE DEFAULT
        call ichorquit('Unknown Case in ICI_CPU_OBS_Gen',-1)
    END SELECT
  end subroutine SPICI_CPU_OBS_Gen
  
END module SPIchorEriCoulombintegralCPUOBSGeneralModGen
