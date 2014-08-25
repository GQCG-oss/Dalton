MODULE IchorEriCoulombintegralCPUOBSGeneralModGen
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000ModGen
use AGC_CPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_CPU_OBS_TRMODAtoCGen
use AGC_CPU_OBS_TRMODAtoDGen
use AGC_CPU_OBS_TRMODBtoCGen
use AGC_CPU_OBS_TRMODBtoDGen
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
public :: IchorCoulombIntegral_CPU_OBS_Gen,IchorCoulombIntegral_CPU_OBS_general_sizeGen  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_CPU_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
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
        call VerticalRecurrenceCPUGen0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*1.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen1(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*1.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen1(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*1.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
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
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
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
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*12),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
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
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,9,Qdistance12,TMParray1(1:nContQP*nPasses*36),&
            & LOCALINTS(1:nContQP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*81.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,9,Qdistance12,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*81),lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*6),&
            & LOCALINTS(1:nContQP*nPasses*5))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*24),&
            & TMParray2(1:nContQP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*20),&
            & LOCALINTS(1:nContQP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*50),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*50),&
            & TMParray1(1:nContQP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*120),&
            & TMParray2(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*100),&
            & TMParray1(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(35,nContQP*nPasses,TMParray1(1:nContQP*nPasses*210),&
            & TMParray2(1:nContQP*nPasses*175))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*175),&
            & TMParray1(1:nContQP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*180),&
            & LOCALINTS(1:nContQP*nPasses*125))
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*72),&
            & TMParray2(1:nContQP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*60),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*180),&
            & TMParray2(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*150),&
            & LOCALINTS(1:nContQP*nPasses*135),lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*180),&
            & TMParray2(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*150),&
            & TMParray1(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*360),&
            & TMParray2(1:nContQP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*300),&
            & TMParray1(1:nContQP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*270),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(35,nContQP*nPasses,TMParray1(1:nContQP*nPasses*630),&
            & TMParray2(1:nContQP*nPasses*525))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*525),&
            & TMParray1(1:nContQP*nPasses*540),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*540),&
            & LOCALINTS(1:nContQP*nPasses*375))
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*36),&
            & LOCALINTS(1:nContQP*nPasses*25))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q1AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*144),&
            & TMParray2(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray2(1:nContQP*nPasses*100),&
            & LOCALINTS(1:nContQP*nPasses*75),lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*360),&
            & TMParray2(1:nContQP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray2(1:nContQP*nPasses*250),&
            & LOCALINTS(1:nContQP*nPasses*225),lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*360),&
            & TMParray2(1:nContQP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,25,Qdistance12,TMParray2(1:nContQP*nPasses*250),&
            & TMParray1(1:nContQP*nPasses*150),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(25,nContQP*nPasses,TMParray1(1:nContQP*nPasses*150),&
            & LOCALINTS(1:nContQP*nPasses*125))
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q3AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*720),&
            & TMParray2(1:nContQP*nPasses*500))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,25,Qdistance12,TMParray2(1:nContQP*nPasses*500),&
            & TMParray1(1:nContQP*nPasses*450),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(25,nContQP*nPasses,TMParray1(1:nContQP*nPasses*450),&
            & LOCALINTS(1:nContQP*nPasses*375))
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*9.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*165.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q4AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*1260.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*875.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(35,nContQP*nPasses,TMParray1(1:nContQP*nPasses*1260),&
            & TMParray2(1:nContQP*nPasses*875))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*900.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,25,Qdistance12,TMParray2(1:nContQP*nPasses*875),&
            & TMParray1(1:nContQP*nPasses*900),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*625.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(25,nContQP*nPasses,TMParray1(1:nContQP*nPasses*900),&
            & LOCALINTS(1:nContQP*nPasses*625))
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,1,Qdistance12,TMParray2(1:nContQP*nPasses*4),&
            & LOCALINTS(1:nContQP*nPasses*3),lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,1,Qdistance12,TMParray1(1:nContQP*nPasses*10),&
            & TMParray2(1:nContQP*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*6),&
            & LOCALINTS(1:nContQP*nPasses*5))
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,1,Qdistance12,TMParray2(1:nContQP*nPasses*4),&
            & LOCALINTS(1:nContQP*nPasses*3),lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,1,Qdistance12,TMParray1(1:nContQP*nPasses*10),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,1,Qdistance12,TMParray1(1:nContQP*nPasses*20),&
            & TMParray2(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,1,Qdistance12,TMParray1(1:nContQP*nPasses*10),&
            & TMParray2(1:nContQP*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*6),&
            & LOCALINTS(1:nContQP*nPasses*5))
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,1,Qdistance12,TMParray1(1:nContQP*nPasses*20),&
            & TMParray2(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,1,Qdistance12,TMParray1(1:nContQP*nPasses*35),&
            & TMParray2(1:nContQP*nPasses*36),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*36),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
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
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*12),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*30),&
            & TMParray2(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*12),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*30),&
            & TMParray2(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*105),&
            & TMParray2(1:nContQP*nPasses*108),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*108),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*6),&
            & LOCALINTS(1:nContQP*nPasses*5))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*24),&
            & TMParray2(1:nContQP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*20),&
            & LOCALINTS(1:nContQP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*50),&
            & TMParray1(1:nContQP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*24),&
            & TMParray2(1:nContQP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*20),&
            & LOCALINTS(1:nContQP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*50),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*120),&
            & TMParray2(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*100),&
            & TMParray1(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*50),&
            & TMParray1(1:nContQP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*120),&
            & TMParray2(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*100),&
            & TMParray1(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(35,nContQP*nPasses,TMParray1(1:nContQP*nPasses*210),&
            & TMParray2(1:nContQP*nPasses*175))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*175),&
            & TMParray1(1:nContQP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*180),&
            & LOCALINTS(1:nContQP*nPasses*125))
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*12),&
            & LOCALINTS(1:nContQP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*30),&
            & TMParray2(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*30),&
            & TMParray2(1:nContQP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,3,Qdistance12,TMParray1(1:nContQP*nPasses*105),&
            & TMParray2(1:nContQP*nPasses*108),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(3,nContQP*nPasses,TMParray2(1:nContQP*nPasses*108),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,9,Qdistance12,TMParray1(1:nContQP*nPasses*36),&
            & LOCALINTS(1:nContQP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,9,Qdistance12,TMParray1(1:nContQP*nPasses*90),&
            & TMParray2(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(9,nContQP*nPasses,TMParray2(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,9,Qdistance12,TMParray1(1:nContQP*nPasses*180),&
            & TMParray2(1:nContQP*nPasses*162),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(9,nContQP*nPasses,TMParray2(1:nContQP*nPasses*162),&
            & LOCALINTS(1:nContQP*nPasses*135))
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,9,Qdistance12,TMParray1(1:nContQP*nPasses*90),&
            & TMParray2(1:nContQP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(9,nContQP*nPasses,TMParray2(1:nContQP*nPasses*54),&
            & LOCALINTS(1:nContQP*nPasses*45))
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,9,Qdistance12,TMParray1(1:nContQP*nPasses*180),&
            & TMParray2(1:nContQP*nPasses*162),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(9,nContQP*nPasses,TMParray2(1:nContQP*nPasses*162),&
            & LOCALINTS(1:nContQP*nPasses*135))
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*315.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*324.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,9,Qdistance12,TMParray1(1:nContQP*nPasses*315),&
            & TMParray2(1:nContQP*nPasses*324),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(9,nContQP*nPasses,TMParray2(1:nContQP*nPasses*324),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(1,nContQP*nPasses,TMParray2(1:nContQP*nPasses*18),&
            & LOCALINTS(1:nContQP*nPasses*15))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*72),&
            & TMParray2(1:nContQP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*60),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*180),&
            & TMParray2(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*150),&
            & TMParray1(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*72),&
            & TMParray2(1:nContQP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*60),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*180),&
            & TMParray2(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*150),&
            & LOCALINTS(1:nContQP*nPasses*135),lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*360),&
            & TMParray2(1:nContQP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*300),&
            & TMParray1(1:nContQP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*270),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*180),&
            & TMParray2(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*150),&
            & TMParray1(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoCGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*360),&
            & TMParray2(1:nContQP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*300),&
            & TMParray1(1:nContQP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*270),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q4CtoBGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(35,nContQP*nPasses,TMParray1(1:nContQP*nPasses*630),&
            & TMParray2(1:nContQP*nPasses*525))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*525),&
            & TMParray1(1:nContQP*nPasses*540),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*540),&
            & LOCALINTS(1:nContQP*nPasses*375))
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*24),&
            & TMParray2(1:nContQP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*20),&
            & LOCALINTS(1:nContQP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*60),&
            & TMParray2(1:nContQP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*50),&
            & TMParray1(1:nContQP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*30),&
            & LOCALINTS(1:nContQP*nPasses*25))
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoAGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*120),&
            & TMParray2(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,5,Qdistance12,TMParray2(1:nContQP*nPasses*100),&
            & TMParray1(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(5,nContQP*nPasses,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*72),&
            & TMParray2(1:nContQP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*60),&
            & LOCALINTS(1:nContQP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*180),&
            & TMParray2(1:nContQP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*150),&
            & TMParray1(1:nContQP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*90),&
            & LOCALINTS(1:nContQP*nPasses*75))
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*360),&
            & TMParray2(1:nContQP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,15,Qdistance12,TMParray2(1:nContQP*nPasses*300),&
            & TMParray1(1:nContQP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,nContQP*nPasses,TMParray1(1:nContQP*nPasses*270),&
            & LOCALINTS(1:nContQP*nPasses*225))
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q1AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nContQP*nPasses,TMParray1(1:nContQP*nPasses*144),&
            & TMParray2(1:nContQP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQP,nPasses,25,Qdistance12,TMParray2(1:nContQP*nPasses*100),&
            & LOCALINTS(1:nContQP*nPasses*75),lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQP*nPasses,TMParray1(1:nContQP*nPasses*360),&
            & TMParray2(1:nContQP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQP,nPasses,25,Qdistance12,TMParray2(1:nContQP*nPasses*250),&
            & TMParray1(1:nContQP*nPasses*150),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(25,nContQP*nPasses,TMParray1(1:nContQP*nPasses*150),&
            & LOCALINTS(1:nContQP*nPasses*125))
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q3AtoDGen(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        nContQP = nContQ*nContP
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nPrimA*nPrimB*nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nContC*nContD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUGen700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUGen700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nContQP*nPasses,TMParray1(1:nContQP*nPasses*720),&
            & TMParray2(1:nContQP*nPasses*500))
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQP,nPasses,25,Qdistance12,TMParray2(1:nContQP*nPasses*500),&
            & TMParray1(1:nContQP*nPasses*450),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(25,nContQP*nPasses,TMParray1(1:nContQP*nPasses*450),&
            & LOCALINTS(1:nContQP*nPasses*375))
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_CPU_OBS_Gen',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_CPU_OBS_Gen
  
  
  subroutine IchorCoulombIntegral_CPU_OBS_general_sizeGen(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP
    integer,intent(in) :: nContA,nContB,nContC,nContD
    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD
    ! local variables
    integer :: AngmomID
    
    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,1*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,1*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,1*nContC*nContD*nContA*nPrimB)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD*nContA*nContB)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD*nContA*nContB)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nContA*nContB)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,35*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,35*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,35*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,35*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQP)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD*nContA*nContB)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQP)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD*nContA*nContB)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQP)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nContA*nContB)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQP)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQP)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQP)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQP)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,315*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,324*nContQP)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQP)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQP)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
    TMParray2maxSize = MAX(TMParray2maxSize,35*nContC*nPrimD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,35*nContC*nContD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,35*nContC*nContD*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,35*nContC*nContD*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQP)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQP)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQP)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,9*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,165*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nPrimQP)
    TMParray1maxSize = MAX(TMParray1maxSize,1225*nContC*nPrimD*nPrimA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,1225*nContC*nContD*nPrimA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,1225*nContC*nContD*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,1225*nContC*nContD*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContQP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_CPU_OBS_general_sizeGen

   subroutine PrimitiveContractionCCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iContC,iPrimC
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      TMP = 0.0E0_realk
      do iPrimC=1,nPrimC
       TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
      enddo
      AUXarrayCont(iContC,iP)=TMP
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen1

   subroutine PrimitiveContractionDCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iContD,iContC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContD,iContC,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
      do iContC=1,nContC
       AUXarrayCont(iContC,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
       do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iP)=AUXarrayCont(iContC,iContD,iP)+TMP*AUXarray2(iContC,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen1

   subroutine PrimitiveContractionACPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iC,iP,iContA,iPrimA
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iC,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
      do iC=1,nContC*nContD
       AUXarrayCont(iC,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
       do iC=1,nContC*nContD
        AUXarrayCont(iC,iContA,iP)=AUXarrayCont(iC,iContA,iP)+TMP*AUXarray2(iC,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen1

   subroutine PrimitiveContractionBCPUGen1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iC,iP,iContB,iPrimB
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iC,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
      do iC=1,nContC*nContD*nContA
       AUXarrayCont(iC,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
       do iC=1,nContC*nContD*nContA
        AUXarrayCont(iC,iContB,iP)=AUXarrayCont(iC,iContB,iP)+TMP*AUXarray2(iC,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen1

   subroutine PrimitiveContractionCCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(    4,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,    4
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,    4
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen4

   subroutine PrimitiveContractionDCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(    4*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,4*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,4*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen4

   subroutine PrimitiveContractionACPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(    4*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,4*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,4*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen4

   subroutine PrimitiveContractionBCPUGen4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(    4*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,4*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,4*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen4

   subroutine PrimitiveContractionCCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   10,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   10
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   10
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen10

   subroutine PrimitiveContractionDCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   10*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,10*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,10*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen10

   subroutine PrimitiveContractionACPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   10*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,10*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,10*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen10

   subroutine PrimitiveContractionBCPUGen10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   10*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,10*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,10*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen10

   subroutine PrimitiveContractionCCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   20,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   20
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   20
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen20

   subroutine PrimitiveContractionDCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   20*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,20*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,20*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen20

   subroutine PrimitiveContractionACPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   20*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,20*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,20*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen20

   subroutine PrimitiveContractionBCPUGen20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   20*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,20*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,20*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen20

   subroutine PrimitiveContractionCCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   35,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   35
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   35
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen35

   subroutine PrimitiveContractionDCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   35*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,35*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,35*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen35

   subroutine PrimitiveContractionACPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   35*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,35*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,35*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen35

   subroutine PrimitiveContractionBCPUGen35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   35*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,35*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,35*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen35

   subroutine PrimitiveContractionCCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   16,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   16
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   16
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen16

   subroutine PrimitiveContractionDCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   16*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,16*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,16*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen16

   subroutine PrimitiveContractionACPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   16*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,16*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,16*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen16

   subroutine PrimitiveContractionBCPUGen16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,16*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,16*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen16

   subroutine PrimitiveContractionCCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   40,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   40
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   40
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen40

   subroutine PrimitiveContractionDCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   40*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,40*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,40*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen40

   subroutine PrimitiveContractionACPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   40*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,40*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,40*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen40

   subroutine PrimitiveContractionBCPUGen40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   40*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,40*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,40*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen40

   subroutine PrimitiveContractionCCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   80,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,   80
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,   80
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen80

   subroutine PrimitiveContractionDCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   80*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,80*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,80*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen80

   subroutine PrimitiveContractionACPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   80*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,80*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,80*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen80

   subroutine PrimitiveContractionBCPUGen80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   80*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,80*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,80*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen80

   subroutine PrimitiveContractionCCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  140,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  140
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  140
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen140

   subroutine PrimitiveContractionDCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  140*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,140*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,140*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen140

   subroutine PrimitiveContractionACPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  140*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,140*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,140*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen140

   subroutine PrimitiveContractionBCPUGen140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  140*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,140*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,140*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen140

   subroutine PrimitiveContractionCCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  100,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  100
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  100
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen100

   subroutine PrimitiveContractionDCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  100*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,100*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,100*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen100

   subroutine PrimitiveContractionACPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  100*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,100*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,100*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen100

   subroutine PrimitiveContractionBCPUGen100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,100*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,100*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen100

   subroutine PrimitiveContractionCCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  200,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  200
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  200
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen200

   subroutine PrimitiveContractionDCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  200*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,200*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,200*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen200

   subroutine PrimitiveContractionACPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  200*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,200*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,200*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen200

   subroutine PrimitiveContractionBCPUGen200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  200*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,200*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,200*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen200

   subroutine PrimitiveContractionCCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  350,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  350
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  350
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen350

   subroutine PrimitiveContractionDCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  350*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,350*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,350*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen350

   subroutine PrimitiveContractionACPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  350*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,350*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,350*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen350

   subroutine PrimitiveContractionBCPUGen350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  350*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,350*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,350*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen350

   subroutine PrimitiveContractionCCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  400,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  400
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  400
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen400

   subroutine PrimitiveContractionDCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  400*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,400*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,400*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen400

   subroutine PrimitiveContractionACPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  400*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,400*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,400*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen400

   subroutine PrimitiveContractionBCPUGen400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,400*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,400*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen400

   subroutine PrimitiveContractionCCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  700,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1,  700
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1,  700
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen700

   subroutine PrimitiveContractionDCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  700*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,700*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,700*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen700

   subroutine PrimitiveContractionACPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  700*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,700*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,700*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen700

   subroutine PrimitiveContractionBCPUGen700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  700*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,700*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,700*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen700

   subroutine PrimitiveContractionCCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD*nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContC,nPrimD*nPrimA*nPrimB*nPasses)
    !
    integer :: iP,iPrimC,iContC,iTUV
    real(realk) :: TMP
!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) 
    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses
     do iContC=1,nContC
      do iTUV=1, 1225
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       TMP = CCC(iPrimC,iContC)
       do iTUV=1, 1225
        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUGen1225

   subroutine PrimitiveContractionDCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2( 1225*nContC,nPrimD,nPrimA*nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContC,nContD,nPrimA*nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContD,iPrimD
    real(realk) :: TMP
!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) 
    do iP = 1,nPrimA*nPrimB*nPasses
     do iContD=1,nContD
       do iTUV=1,1225*nContC
       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       TMP = DCC(iPrimD,iContD)
        do iTUV=1,1225*nContC
        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUGen1225

   subroutine PrimitiveContractionACPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2( 1225*nContC*nContD,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContC*nContD,nContA,nPrimB*nPasses)
    !
    integer :: iTUV,iP,iContA,iPrimA
    real(realk) :: TMP
!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
    do iP = 1,nPrimB*nPasses
     do iContA=1,nContA
       do iTUV=1,1225*nContC*nContD
       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
      enddo
      do iPrimA=1,nPrimA
       TMP = ACC(iPrimA,iContA)
        do iTUV=1,1225*nContC*nContD
        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUGen1225

   subroutine PrimitiveContractionBCPUGen1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&
       & nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225*nContC*nContD*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContC*nContD*nContA,nContB,nPasses)
    !
    integer :: iTUV,iP,iContB,iPrimB
    real(realk) :: TMP
!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) 
    do iP = 1,nPasses
     do iContB=1,nContB
       do iTUV=1,1225*nContC*nContD*nContA
       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk
      enddo
      do iPrimB=1,nPrimB
       TMP = BCC(iPrimB,iContB)
        do iTUV=1,1225*nContC*nContD*nContA
        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUGen1225
END MODULE IchorEriCoulombintegralCPUOBSGeneralModGen
