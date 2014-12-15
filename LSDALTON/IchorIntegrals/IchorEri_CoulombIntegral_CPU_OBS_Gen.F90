MODULE IchorEriCoulombintegralCPUOBSGeneralModGen
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionMod
use IchorCommonMod
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000ModGen
use AGC_CPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_CPU_PrimContractGenMod
use AGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_CPU_OBS_TRMODAtoCGen1
use AGC_CPU_OBS_TRMODAtoCGen2
use AGC_CPU_OBS_TRMODAtoDGen1
use AGC_CPU_OBS_TRMODAtoDGen2
use AGC_CPU_OBS_TRMODBtoCGen1
use AGC_CPU_OBS_TRMODBtoDGen1
use AGC_CPU_OBS_TRMODCtoAGen
use AGC_CPU_OBS_TRMODDtoAGen
use AGC_CPU_OBS_TRMODCtoBGen
use AGC_CPU_OBS_TRMODDtoBGen
use AGC_CPU_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_CPU_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_CPU_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_CPU_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_CPU_OBS_Sphcontract1Mod
use AGC_CPU_OBS_Sphcontract2Mod
  
private   
public :: ICI_CPU_OBS_Gen
  
CONTAINS
  
  
  subroutine ICI_CPU_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPasses,MaxPasses,IntPrint,lupri,&
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
    integer,intent(in) :: MaxPasses,IntPrint,lupri
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
    real(realk),intent(in) :: integralPrefactor(nPrimQ*nPrimP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQ*nPrimP)
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
!   Local variables 
    integer :: AngmomPQ,AngmomP,AngmomQ,I,J,la,lb,lc,ld,nsize,angmomid
    
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
        call VerticalRecurrenceCPUGen0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,TMParray2(1))
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionCCPUGen1(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen1(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen1(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen1(TMParray1,LOCALINTS,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceCPUGen1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2(1))
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionCCPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,3,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,3,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,9,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,9,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(1,nContQ*nContP*nPasses,TMParray2(1),&
            & LOCALINTS(1))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,5,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,5,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,5,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,5,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(35,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,5,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionCCPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(1,nContQ*nContP*nPasses,TMParray2(1),&
            & LOCALINTS(1))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,15,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,15,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,15,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,15,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call BuildRJ000CPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(35,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,15,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
         call PrimitiveContractionCCPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(1,nContQ*nContP*nPasses,TMParray2(1),&
            & LOCALINTS(1))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP4Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,25,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,25,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,25,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(25,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call BuildRJ000CPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP4Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,25,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(25,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call BuildRJ000CPUGen8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUGen8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP4Q4AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
         call PrimitiveContractionCCPUGen1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUGen1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionACPUGen1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
         call PrimitiveContractionBCPUGen1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(35,nContQ*nContP*nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,25,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(25,nContQ*nContP*nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in ICI_CPU_OBS_Gen',-1)
    END SELECT
  end subroutine ICI_CPU_OBS_Gen
  
END MODULE IchorEriCoulombintegralCPUOBSGeneralModGen
