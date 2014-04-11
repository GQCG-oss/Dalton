MODULE IchorEriCoulombintegralGPUOBSGeneralModGen
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_GPU_OBS_BUILDRJ000ModGen
use AGC_GPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_GPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_GPU_OBS_TRMODAtoCGen
use AGC_GPU_OBS_TRMODAtoDGen
use AGC_GPU_OBS_TRMODBtoCGen
use AGC_GPU_OBS_TRMODBtoDGen
use AGC_GPU_OBS_TRMODCtoAGen
use AGC_GPU_OBS_TRMODDtoAGen
use AGC_GPU_OBS_TRMODCtoBGen
use AGC_GPU_OBS_TRMODDtoBGen
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_GPU_OBS_Gen,IchorCoulombIntegral_GPU_OBS_general_sizeGen  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_GPU_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
       & BasisCont1,BasisCont2,BasisCont3,IatomAPass,iatomBPass)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsA,nAtomsB
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
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
    integer,intent(in) :: BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize
    real(realk) :: BasisCont1(BasisCont1maxsize) 
    real(realk) :: BasisCont2(BasisCont2maxsize) 
    real(realk) :: BasisCont3(BasisCont3maxsize) 
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
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*1.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*1.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen1(TMParray2,LOCALINTS,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*12.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*12),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,9,Qdistance12,TMParray2(1:nContQP*nPasses*36),&
            & LOCALINTS(1:nContQP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*81.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,9,Qdistance12,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*81),lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*6.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*6),&
            & LOCALINTS(1:nContQP*nPasses*5))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*24.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*24),&
            & TMParray1(1:nContQP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*20),&
            & LOCALINTS(1:nContQP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*50),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*50),&
            & TMParray2(1:nContQP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*120),&
            & TMParray1(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*100),&
            & TMParray2(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*210.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*175.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(35,nContQP*nPasses,TMParray2(1:nContQP*nPasses*210),&
            & TMParray1(1:nContQP*nPasses*175))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*175),&
            & TMParray2(1:nContQP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*180),&
            & LOCALINTS(1:nContQP*nPasses*125))
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*72.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*72),&
            & TMParray1(1:nContQP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*60),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*180),&
            & TMParray1(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*150),&
            & LOCALINTS(1:nContQP*nPasses*135),lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*180),&
            & TMParray1(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*150),&
            & TMParray2(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*300.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*360),&
            & TMParray1(1:nContQP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*270.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*300),&
            & TMParray2(1:nContQP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*270),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*630.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*525.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(35,nContQP*nPasses,TMParray2(1:nContQP*nPasses*630),&
            & TMParray1(1:nContQP*nPasses*525))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*540.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*525),&
            & TMParray2(1:nContQP*nPasses*540),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*540),&
            & LOCALINTS(1:nContQP*nPasses*375))
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*36),&
            & LOCALINTS(1:nContQP*nPasses*25))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*144.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*144),&
            & TMParray1(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray1(1:nContQP*nPasses*100),&
            & LOCALINTS(1:nContQP*nPasses*75),lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*250.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*360),&
            & TMParray1(1:nContQP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray1(1:nContQP*nPasses*250),&
            & LOCALINTS(1:nContQP*nPasses*225),lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*250.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*360),&
            & TMParray1(1:nContQP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray1(1:nContQP*nPasses*250),&
            & TMParray2(1:nContQP*nPasses*150),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(25,nContQP*nPasses,TMParray2(1:nContQP*nPasses*150),&
            & LOCALINTS(1:nContQP*nPasses*125))
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*720.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*500.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*720),&
            & TMParray1(1:nContQP*nPasses*500))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*450.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray1(1:nContQP*nPasses*500),&
            & TMParray2(1:nContQP*nPasses*450),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(25,nContQP*nPasses,TMParray2(1:nContQP*nPasses*450),&
            & LOCALINTS(1:nContQP*nPasses*375))
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*9.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*165.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q4AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*1260.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*875.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(35,nContQP*nPasses,TMParray2(1:nContQP*nPasses*1260),&
            & TMParray1(1:nContQP*nPasses*875))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*900.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,25,Qdistance12,TMParray1(1:nContQP*nPasses*875),&
            & TMParray2(1:nContQP*nPasses*900),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*625.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(25,nContQP*nPasses,TMParray2(1:nContQP*nPasses*900),&
            & LOCALINTS(1:nContQP*nPasses*625))
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,1,Qdistance12,TMParray1(1:nContQP*nPasses*4),&
            & LOCALINTS(1:nContQP*nPasses*3),lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*6.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,1,Qdistance12,TMParray2(1:nContQP*nPasses*10),&
            & TMParray1(1:nContQP*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*6),&
            & LOCALINTS(1:nContQP*nPasses*5))
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,1,Qdistance12,TMParray1(1:nContQP*nPasses*4),&
            & LOCALINTS(1:nContQP*nPasses*3),lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,1,Qdistance12,TMParray2(1:nContQP*nPasses*10),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,1,Qdistance12,TMParray2(1:nContQP*nPasses*20),&
            & TMParray1(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*6.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,1,Qdistance12,TMParray2(1:nContQP*nPasses*10),&
            & TMParray1(1:nContQP*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*6),&
            & LOCALINTS(1:nContQP*nPasses*5))
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,1,Qdistance12,TMParray2(1:nContQP*nPasses*20),&
            & TMParray1(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,1,Qdistance12,TMParray2(1:nContQP*nPasses*35),&
            & TMParray1(1:nContQP*nPasses*36),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*36),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*12.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*12),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*30),&
            & TMParray1(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*12.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*12),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q3DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*30),&
            & TMParray1(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q3CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*105.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*108.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*105),&
            & TMParray1(1:nContQP*nPasses*108),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*108),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*6.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*6),&
            & LOCALINTS(1:nContQP*nPasses*5))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*24.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*24),&
            & TMParray1(1:nContQP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*20),&
            & LOCALINTS(1:nContQP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*50),&
            & TMParray2(1:nContQP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*24.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*24),&
            & TMParray1(1:nContQP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*20),&
            & LOCALINTS(1:nContQP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*50),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*120),&
            & TMParray1(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*100),&
            & TMParray2(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*50),&
            & TMParray2(1:nContQP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*120),&
            & TMParray1(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*100),&
            & TMParray2(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*210.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*175.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(35,nContQP*nPasses,TMParray2(1:nContQP*nPasses*210),&
            & TMParray1(1:nContQP*nPasses*175))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*175),&
            & TMParray2(1:nContQP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*180),&
            & LOCALINTS(1:nContQP*nPasses*125))
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*12.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*12),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*30),&
            & TMParray1(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*30),&
            & TMParray1(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*105.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*108.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,3,Qdistance12,TMParray2(1:nContQP*nPasses*105),&
            & TMParray1(1:nContQP*nPasses*108),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(3,nContQP*nPasses,TMParray1(1:nContQP*nPasses*108),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,9,Qdistance12,TMParray2(1:nContQP*nPasses*36),&
            & LOCALINTS(1:nContQP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,9,Qdistance12,TMParray2(1:nContQP*nPasses*90),&
            & TMParray1(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(9,nContQP*nPasses,TMParray1(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*162.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,9,Qdistance12,TMParray2(1:nContQP*nPasses*180),&
            & TMParray1(1:nContQP*nPasses*162),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(9,nContQP*nPasses,TMParray1(1:nContQP*nPasses*162),&
            & LOCALINTS(1:nContQP*nPasses*135))
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,9,Qdistance12,TMParray2(1:nContQP*nPasses*90),&
            & TMParray1(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(9,nContQP*nPasses,TMParray1(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*162.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,9,Qdistance12,TMParray2(1:nContQP*nPasses*180),&
            & TMParray1(1:nContQP*nPasses*162),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(9,nContQP*nPasses,TMParray1(1:nContQP*nPasses*162),&
            & LOCALINTS(1:nContQP*nPasses*135))
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*315.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*324.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,9,Qdistance12,TMParray2(1:nContQP*nPasses*315),&
            & TMParray1(1:nContQP*nPasses*324),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(9,nContQP*nPasses,TMParray1(1:nContQP*nPasses*324),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(1,nContQP*nPasses,TMParray1(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*72.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*72),&
            & TMParray1(1:nContQP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*60),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*180),&
            & TMParray1(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*150),&
            & TMParray2(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*72.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*72),&
            & TMParray1(1:nContQP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*60),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*180),&
            & TMParray1(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*150),&
            & LOCALINTS(1:nContQP*nPasses*135),lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q3BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*300.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*360),&
            & TMParray1(1:nContQP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*270.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*300),&
            & TMParray2(1:nContQP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*270),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*180),&
            & TMParray1(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*150),&
            & TMParray2(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q3BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*300.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*360),&
            & TMParray1(1:nContQP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*270.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*300),&
            & TMParray2(1:nContQP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*270),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*630.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*525.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(35,nContQP*nPasses,TMParray2(1:nContQP*nPasses*630),&
            & TMParray1(1:nContQP*nPasses*525))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*540.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*525),&
            & TMParray2(1:nContQP*nPasses*540),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*540),&
            & LOCALINTS(1:nContQP*nPasses*375))
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*24.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*24),&
            & TMParray1(1:nContQP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*20),&
            & LOCALINTS(1:nContQP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*60),&
            & TMParray1(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*50),&
            & TMParray2(1:nContQP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*120),&
            & TMParray1(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray1(1:nContQP*nPasses*100),&
            & TMParray2(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(5,nContQP*nPasses,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*72.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*72),&
            & TMParray1(1:nContQP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*60),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*180),&
            & TMParray1(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*150),&
            & TMParray2(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q3AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*300.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*360),&
            & TMParray1(1:nContQP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*270.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray1(1:nContQP*nPasses*300),&
            & TMParray2(1:nContQP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(15,nContQP*nPasses,TMParray2(1:nContQP*nPasses*270),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*144.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(4,nContQP*nPasses,TMParray2(1:nContQP*nPasses*144),&
            & TMParray1(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQP,nPasses,25,Qdistance12,TMParray1(1:nContQP*nPasses*100),&
            & LOCALINTS(1:nContQP*nPasses*75),lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*250.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray2(1:nContQP*nPasses*360),&
            & TMParray1(1:nContQP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQP,nPasses,25,Qdistance12,TMParray1(1:nContQP*nPasses*250),&
            & TMParray2(1:nContQP*nPasses*150),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(25,nContQP*nPasses,TMParray2(1:nContQP*nPasses*150),&
            & LOCALINTS(1:nContQP*nPasses*125))
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q3AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nContP*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nContP*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionGPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*720.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*500.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(20,nContQP*nPasses,TMParray2(1:nContQP*nPasses*720),&
            & TMParray1(1:nContQP*nPasses*500))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*450.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQP,nPasses,25,Qdistance12,TMParray1(1:nContQP*nPasses*500),&
            & TMParray2(1:nContQP*nPasses*450),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(25,nContQP*nPasses,TMParray2(1:nContQP*nPasses*450),&
            & LOCALINTS(1:nContQP*nPasses*375))
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_GPU_OBS_Gen',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_GPU_OBS_Gen
  
  
  subroutine IchorCoulombIntegral_GPU_OBS_general_sizeGen(TMParray1maxsize,&
         & TMParray2maxsize,BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
         & AngmomA,AngmomB,AngmomC,AngmomD,nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(inout) :: BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP
    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD
    ! local variables
    integer :: AngmomID
    
    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    BasisCont1maxsize = 1
    BasisCont2maxsize = 1
    BasisCont3maxsize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
    BasisCont1maxsize =     1*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =     1*nPrimA*nPrimB
    BasisCont3maxsize =     1*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimQP)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
    BasisCont1maxsize =     4*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =     4*nPrimA*nPrimB
    BasisCont3maxsize =     4*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQP)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
    BasisCont1maxsize =    10*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    10*nPrimA*nPrimB
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,10*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,6*nContQP)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
    BasisCont1maxsize =     4*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =     4*nPrimA*nPrimB
    BasisCont3maxsize =     4*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQP)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
    BasisCont1maxsize =    10*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    10*nPrimA*nPrimB
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,10*nContQP)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
    BasisCont1maxsize =    20*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    20*nPrimA*nPrimB
    BasisCont3maxsize =    20*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,18*nContQP)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
    BasisCont1maxsize =    10*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    10*nPrimA*nPrimB
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,10*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,6*nContQP)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
    BasisCont1maxsize =    20*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    20*nPrimA*nPrimB
    BasisCont3maxsize =    20*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,18*nContQP)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
    BasisCont1maxsize =    35*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    35*nPrimA*nPrimB
    BasisCont3maxsize =    35*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,35*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQP)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
    BasisCont1maxsize =     4*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =     4*nPrimA*nPrimB
    BasisCont3maxsize =     4*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQP)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
    BasisCont1maxsize =    16*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    16*nPrimA*nPrimB
    BasisCont3maxsize =    16*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,12*nContQP)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,18*nContQP)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
    BasisCont1maxsize =    16*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    16*nPrimA*nPrimB
    BasisCont3maxsize =    16*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,12*nContQP)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
    BasisCont1maxsize =    80*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    80*nPrimA*nPrimB
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,54*nContQP)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,18*nContQP)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
    BasisCont1maxsize =    80*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    80*nPrimA*nPrimB
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,54*nContQP)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
    BasisCont1maxsize =   140*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   140*nPrimA*nPrimB
    BasisCont3maxsize =   140*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,105*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,108*nContQP)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
    BasisCont1maxsize =    10*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    10*nPrimA*nPrimB
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,10*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,6*nContQP)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,24*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,50*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,24*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,50*nContQP)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,120*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,50*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,120*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
    BasisCont1maxsize =   350*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   350*nPrimA*nPrimB
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,210*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,175*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
    BasisCont1maxsize =     4*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =     4*nPrimA*nPrimB
    BasisCont3maxsize =     4*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQP)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
    BasisCont1maxsize =    16*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    16*nPrimA*nPrimB
    BasisCont3maxsize =    16*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,12*nContQP)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,18*nContQP)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
    BasisCont1maxsize =    16*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    16*nPrimA*nPrimB
    BasisCont3maxsize =    16*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,12*nContQP)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
    BasisCont1maxsize =    80*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    80*nPrimA*nPrimB
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,54*nContQP)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,18*nContQP)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
    BasisCont1maxsize =    80*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    80*nPrimA*nPrimB
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,54*nContQP)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
    BasisCont1maxsize =   140*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   140*nPrimA*nPrimB
    BasisCont3maxsize =   140*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,105*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,108*nContQP)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
    BasisCont1maxsize =    10*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    10*nPrimA*nPrimB
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,10*nContQP)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQP)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,54*nContQP)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQP)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,162*nContQP)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,54*nContQP)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,162*nContQP)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
    BasisCont1maxsize =   350*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   350*nPrimA*nPrimB
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,315*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,324*nContQP)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
    BasisCont1maxsize =    20*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    20*nPrimA*nPrimB
    BasisCont3maxsize =    20*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,18*nContQP)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
    BasisCont1maxsize =    80*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    80*nPrimA*nPrimB
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,72*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
    BasisCont1maxsize =    80*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    80*nPrimA*nPrimB
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,72*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
    BasisCont1maxsize =   400*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   400*nPrimA*nPrimB
    BasisCont3maxsize =   400*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,270*nContQP)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
    BasisCont1maxsize =   400*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   400*nPrimA*nPrimB
    BasisCont3maxsize =   400*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,270*nContQP)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
    BasisCont1maxsize =   700*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   700*nPrimA*nPrimB
    BasisCont3maxsize =   700*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,630*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,525*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,540*nContQP)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
    BasisCont1maxsize =    10*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    10*nPrimA*nPrimB
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,10*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,6*nContQP)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,24*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,50*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
    BasisCont1maxsize =    40*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    40*nPrimA*nPrimB
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,24*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,50*nContQP)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,120*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
    BasisCont1maxsize =   100*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   100*nPrimA*nPrimB
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,50*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContQP)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,120*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
    BasisCont1maxsize =   350*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   350*nPrimA*nPrimB
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,210*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,175*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
    BasisCont1maxsize =    20*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    20*nPrimA*nPrimB
    BasisCont3maxsize =    20*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,18*nContQP)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
    BasisCont1maxsize =    80*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    80*nPrimA*nPrimB
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,72*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
    BasisCont1maxsize =    80*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    80*nPrimA*nPrimB
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,72*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
    BasisCont1maxsize =   400*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   400*nPrimA*nPrimB
    BasisCont3maxsize =   400*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,270*nContQP)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
    BasisCont1maxsize =   200*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   200*nPrimA*nPrimB
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,180*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContQP)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
    BasisCont1maxsize =   400*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   400*nPrimA*nPrimB
    BasisCont3maxsize =   400*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,270*nContQP)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
    BasisCont1maxsize =   700*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   700*nPrimA*nPrimB
    BasisCont3maxsize =   700*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,630*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,525*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,540*nContQP)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
    BasisCont1maxsize =    35*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =    35*nPrimA*nPrimB
    BasisCont3maxsize =    35*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,35*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQP)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
    BasisCont1maxsize =   140*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   140*nPrimA*nPrimB
    BasisCont3maxsize =   140*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,144*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
    BasisCont1maxsize =   350*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   350*nPrimA*nPrimB
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,250*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
    BasisCont1maxsize =   140*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   140*nPrimA*nPrimB
    BasisCont3maxsize =   140*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,144*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContQP)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
    BasisCont1maxsize =   350*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   350*nPrimA*nPrimB
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,250*nContQP)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
    BasisCont1maxsize =   700*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   700*nPrimA*nPrimB
    BasisCont3maxsize =   700*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,720*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,500*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,450*nContQP)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
    BasisCont1maxsize =   350*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   350*nPrimA*nPrimB
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,250*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
    BasisCont1maxsize =   700*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =   700*nPrimA*nPrimB
    BasisCont3maxsize =   700*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,720*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,500*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,450*nContQP)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
    BasisCont1maxsize =  1225*nPrimD*nPrimA*nPrimB
    BasisCont2maxsize =  1225*nPrimA*nPrimB
    BasisCont3maxsize =  1225*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,9*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,165*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1260*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,875*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,900*nContQP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_GPU_OBS_general_sizeGen

  subroutine PrimitiveContractionGPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: TMP
    real(realk) :: BasisCont1(nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(nPrimA,nPrimB)
    real(realk) :: BasisCont3(nPrimB)
    !Scaling p**4*c*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nPassQ
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2,&
!!$OMP        BasisCont1,BasisCont2,BasisCont3)
!!$OMP SINGLE
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iPrimD,iPrimC,TMP,iContC,iPassP)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         TMP = 0.0E0_realk
         do iPrimC=1,nPrimC
          TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
         enddo
         BasisCont1(iPrimD,iPrimA,iPrimB) = TMP
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iPrimD,TMP,iContD,iContC,iPassP)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         TMP = 0.0E0_realk
         do iPrimD=1,nPrimD
          TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iPrimD,iPrimA,iPrimB)
         enddo
         BasisCont2(iPrimA,iPrimB) = TMP
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO PRIVATE(iPrimB,iPrimA,TMP,iContA,iContD,iContC,iPassP)
        do iPrimB=1,nPrimB
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iPrimA,iPrimB)
         enddo
         BasisCont3(iPrimB) = TMP
        enddo
!!$OMP END DO
!!$OMP DO PRIVATE(iPrimB,TMP,iContA,iContB,iContD,iContC,iPassP)
        do iContB=1,nContB
         TMP = 0.0E0_realk
         do iPrimB=1,nPrimB
          TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iPrimB)
         enddo
         AUXarrayCont(iContC,iContD,iContA,iContB,iPassP) = TMP
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END SINGLE
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
  end subroutine PrimitiveContractionGPUGen1

   subroutine PrimitiveContractionGPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(    4,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(    4,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(    4,nPrimA,nPrimB)
    real(realk) :: BasisCont3(    4,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,    4
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,    4
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,    4
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,    4
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen4

   subroutine PrimitiveContractionGPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   10,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(   10,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(   10,nPrimA,nPrimB)
    real(realk) :: BasisCont3(   10,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   10
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   10
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,   10
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,   10
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen10

   subroutine PrimitiveContractionGPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   20,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(   20,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(   20,nPrimA,nPrimB)
    real(realk) :: BasisCont3(   20,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   20
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   20
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,   20
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,   20
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen20

   subroutine PrimitiveContractionGPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   35,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(   35,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(   35,nPrimA,nPrimB)
    real(realk) :: BasisCont3(   35,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   35
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   35
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,   35
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,   35
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen35

   subroutine PrimitiveContractionGPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   16,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(   16,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(   16,nPrimA,nPrimB)
    real(realk) :: BasisCont3(   16,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   16
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   16
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,   16
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,   16
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen16

   subroutine PrimitiveContractionGPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   40,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(   40,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(   40,nPrimA,nPrimB)
    real(realk) :: BasisCont3(   40,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   40
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   40
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,   40
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,   40
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen40

   subroutine PrimitiveContractionGPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   80,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(   80,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(   80,nPrimA,nPrimB)
    real(realk) :: BasisCont3(   80,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,   80
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,   80
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,   80
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,   80
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen80

   subroutine PrimitiveContractionGPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  140,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(  140,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(  140,nPrimA,nPrimB)
    real(realk) :: BasisCont3(  140,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  140
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  140
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,  140
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,  140
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen140

   subroutine PrimitiveContractionGPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  100,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(  100,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(  100,nPrimA,nPrimB)
    real(realk) :: BasisCont3(  100,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  100
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  100
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,  100
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,  100
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen100

   subroutine PrimitiveContractionGPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  200,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(  200,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(  200,nPrimA,nPrimB)
    real(realk) :: BasisCont3(  200,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  200
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  200
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,  200
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,  200
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen200

   subroutine PrimitiveContractionGPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  350,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(  350,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(  350,nPrimA,nPrimB)
    real(realk) :: BasisCont3(  350,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  350
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  350
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,  350
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,  350
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen350

   subroutine PrimitiveContractionGPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  400,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(  400,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(  400,nPrimA,nPrimB)
    real(realk) :: BasisCont3(  400,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  400
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  400
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,  400
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,  400
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen400

   subroutine PrimitiveContractionGPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  700,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1(  700,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2(  700,nPrimA,nPrimB)
    real(realk) :: BasisCont3(  700,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1,  700
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1,  700
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1,  700
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1,  700
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen700

   subroutine PrimitiveContractionGPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD,BasisCont1,BasisCont2,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContC,nContD,nContA,nContB,nPasses)
    !
    integer :: iPassP,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iContQP,iPrimQP
    real(realk) :: TMP
    real(realk) :: BasisCont1( 1225,nPrimD,nPrimA,nPrimB)
    real(realk) :: BasisCont2( 1225,nPrimA,nPrimB)
    real(realk) :: BasisCont3( 1225,nPrimB)
    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP PRIVATE(iTUV,iPassP,iContC,iContD,iContA,iContB,iPrimC,iPrimD,iPrimA,iPrimB,&
!!$OMP         BasisCont1,BasisCont2,BasisCont3,TMP) &
!!$OMP SHARED(nContC,nContD,nContA,nContB,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&
!!$OMP        ACC,BCC,CCC,DCC,AUXarrayCont,AUXarray2)
!$OMP SINGLE
    do iPassP = 1,nPasses
     do iContC=1,nContC
!!$OMP DO COLLAPSE(4) PRIVATE(iPrimB,iPrimA,iPrimD,iTUV,iPrimC,TMP,iPassP,iContC)
      do iPrimB=1,nPrimB
       do iPrimA=1,nPrimA
        do iPrimD=1,nPrimD
         do iTUV=1, 1225
          TMP = 0.0E0_realk
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
          do iPrimC=1,nPrimC
           TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassP)
          enddo
          BasisCont1(iTUV,iPrimD,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
      enddo
!!$OMP END DO
      do iContD=1,nContD
!!$OMP DO COLLAPSE(3) PRIVATE(iPrimB,iPrimA,iTUV,iPrimD,TMP,iPassP,iContC,iContD)
       do iPrimB=1,nPrimB
        do iPrimA=1,nPrimA
         do iTUV=1, 1225
          TMP = 0.0E0_realk
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
          do iPrimD=1,nPrimD
           TMP = TMP + DCC(iPrimD,iContD)*BasisCont1(iTUV,iPrimD,iPrimA,iPrimB)
          enddo
          BasisCont2(iTUV,iPrimA,iPrimB) = TMP
         enddo
        enddo
       enddo
!!$OMP END DO
       do iContA=1,nContA
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iPrimA,iTUV,TMP,iPassP,iContC,iContD,iContA)
        do iPrimB=1,nPrimB
         do iTUV=1, 1225
          TMP = 0.0E0_realk
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContA*nContC*nContD*nTUV*nPassQ
          do iPrimA=1,nPrimA
           TMP = TMP + ACC(iPrimA,iContA)*BasisCont2(iTUV,iPrimA,iPrimB)
          enddo
          BasisCont3(iTUV,iPrimB) = TMP
         enddo
        enddo
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(iPrimB,iTUV,TMP,iPassP,iContC,iContD,iContA,iContB)
        do iContB=1,nContB
         do iTUV=1, 1225
          TMP = 0.0E0_realk
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContB*nContC*nContD*nTUV*nPassQ
          do iPrimB=1,nPrimB
           TMP = TMP + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
          enddo
          AUXarrayCont(iTUV,iContC,iContD,iContA,iContB,iPassP) = TMP
         enddo
        enddo
!!$OMP END DO
       enddo
      enddo
     enddo
    enddo
!!$OMP END PARALLEL DO
!$OMP END SINGLE
!$OMP BARRIER
   end subroutine PrimitiveContractionGPUGen1225
END MODULE IchorEriCoulombintegralGPUOBSGeneralModGen
