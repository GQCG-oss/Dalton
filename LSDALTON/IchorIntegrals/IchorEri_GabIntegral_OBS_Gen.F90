MODULE IchorEriGabintegralOBSGeneralModGen
!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000MODGen
use AGC_CPU_OBS_BUILDRJ000MODSeg1Prim
use IchorEriCoulombintegralOBSGeneralModGen
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
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorGabIntegral_OBS_Gen,IchorGabIntegral_OBS_general_sizeGen  
  
CONTAINS
  
  
  subroutine IchorGabIntegral_OBS_Gen(nPrimA,nPrimB,&
       & nPrimP,IntPrint,lupri,&
       & nContA,nContB,nContP,pexp,ACC,BCC,&
       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,Pdistance12,PQorder,LOCALINTS,Acenter,Bcenter,&
       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisContmaxsize,BasisCont)
    implicit none
    integer,intent(in) :: nPrimP,nPrimA,nPrimB
    integer,intent(in) :: IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContP,nTABFJW1,nTABFJW2
    integer,intent(in) :: AngmomA,AngmomB
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)
    logical,intent(in)     :: Psegmented
    real(realk),intent(in) :: pexp(nPrimP)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(inout) :: LOCALINTS(1)
    real(realk),intent(in) :: integralPrefactor(nPrimP*nPrimP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimP)
    real(realk),intent(in) :: reducedExponents(nPrimP*nPrimP)
    !reducedExponents(nPrimP,nPrimP)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize,BasisContmaxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    real(realk),intent(inout) :: BasisCont(BasisContmaxsize)
!   Local variables 
    integer :: AngmomP,I,J,nContQP,la,lb,lc,ld,nsize,angmomid,IatomAPass(1),IatomBPass(1)
    
    !Setup combined Angmom info
    AngmomP = AngmomA+AngmomB
    IatomAPass(1) = 1
    IatomBPass(1) = 1
!    nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
!    nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
!    nlmA = 2*AngmomA+1
!    nlmB = 2*AngmomB+1
    AngmomID = 10*AngmomA+AngmomB
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*1.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen0(1,nPrimP,nPrimP,&
               & reducedExponents,TABFJW,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*1.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen1(TMParray2(1:nPrimP*nPrimP*1),&
             & TMParray1(1:nContP*1),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        call ExtractGabElmP1Gen(TMParray1,LOCALINTS,nContP)
    CASE(  10)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoCGen(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen16(TMParray2(1:nPrimP*nPrimP*16),&
             & TMParray1(1:nContP*16),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*12.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP,1,4,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*9.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,1,3,Pdistance12,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP3Gen(TMParray1,LOCALINTS,nContP)
    CASE(  11)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCGen(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen100(TMParray2(1:nPrimP*nPrimP*100),&
             & TMParray1(1:nContP*100),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*90.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*81.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,1,9,Pdistance12,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP9Gen(TMParray1,LOCALINTS,nContP)
    CASE(  20)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCGen(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen100(TMParray2(1:nPrimP*nPrimP*100),&
             & TMParray1(1:nContP*100),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContP,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*50.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContP,TMParray2,&
            & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nContP*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,1,5,Pdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*25.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContP,TMParray2,&
            & TMParray1)
        call ExtractGabElmP5Gen(TMParray1,LOCALINTS,nContP)
    CASE(  21)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoCGen(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen400(TMParray2(1:nPrimP*nPrimP*400),&
             & TMParray1(1:nContP*400),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP,1,20,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*300.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nContP,TMParray2,&
            & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nContP*270.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,1,15,Pdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*225.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContP,TMParray2,&
            & TMParray1)
        call ExtractGabElmP15Gen(TMParray1,LOCALINTS,nContP)
    CASE(  22)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*9.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen8(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*165.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen8A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q4AtoCGen(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen1225(TMParray2(1:nPrimP*nPrimP*1225),&
             & TMParray1(1:nContP*1225),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*1260.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP,1,35,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*875.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,nContP,TMParray2,&
            & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nContP*900.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,1,25,Pdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*625.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nContP,TMParray2,&
            & TMParray1)
        call ExtractGabElmP25Gen(TMParray1,LOCALINTS,nContP)
    CASE(   1)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen2(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen2B(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoDGen(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen16(TMParray2(1:nPrimP*nPrimP*16),&
             & TMParray1(1:nContP*16),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*12.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP,1,4,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*9.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,1,3,Pdistance12,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP3Gen(TMParray1,LOCALINTS,nContP)
    CASE(   2)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen4(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen4B(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoDGen(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen100(TMParray2(1:nPrimP*nPrimP*100),&
             & TMParray1(1:nContP*100),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContP,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*50.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContP,TMParray2,&
            & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nContP*30.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,1,5,Pdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*25.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContP,TMParray2,&
            & TMParray1)
        call ExtractGabElmP5Gen(TMParray1,LOCALINTS,nContP)
    CASE(  12)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call BuildRJ000CPUGen6(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUGen6B(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoDGen(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen400(TMParray2(1:nPrimP*nPrimP*400),&
             & TMParray1(1:nContP*400),nPrimP,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*360.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP,1,20,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*300.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nContP,TMParray2,&
            & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nContP*270.GT.TMParray2maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,1,15,Pdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*225.GT.TMParray1maxsize)THEN
          call ichorquit('nContPtoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContP,TMParray2,&
            & TMParray1)
        call ExtractGabElmP15Gen(TMParray1,LOCALINTS,nContP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorGabIntegral_OBS_Gen',-1)
    END SELECT
  end subroutine IchorGabIntegral_OBS_Gen
  
  
  subroutine IchorGabIntegral_OBS_general_sizeGen(TMParray1maxsize,&
         &TMParray2maxsize,BasisContmaxsize,AngmomA,AngmomB,nPrimP,nContP,nPrimB)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize,BasisContmaxsize
    integer,intent(in) :: AngmomA,AngmomB
    integer,intent(in) :: nPrimP,nContP,nPrimB
    ! local variables
    integer :: AngmomID
    
    AngmomID = 10*AngmomA+AngmomB
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    BasisContmaxsize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
       BasisContmaxsize =     1*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,1*nContP)
    CASE(   1)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       BasisContmaxsize =    16*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,12*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,9*nContP)
    CASE(   2)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       BasisContmaxsize =   100*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,50*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,25*nContP)
    CASE(  10)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       BasisContmaxsize =    16*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,12*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,9*nContP)
    CASE(  11)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       BasisContmaxsize =   100*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,90*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,81*nContP)
    CASE(  12)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       BasisContmaxsize =   400*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,270*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,225*nContP)
    CASE(  20)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       BasisContmaxsize =   100*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,50*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,30*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,25*nContP)
    CASE(  21)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       BasisContmaxsize =   400*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,360*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,300*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,270*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,225*nContP)
    CASE(  22)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       BasisContmaxsize =  1225*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,9*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,165*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,1260*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,875*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,900*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,625*nContP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorGabIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorGabIntegral_OBS_general_sizeGen
  subroutine GabPrimitiveContractionGen1(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB)
    !
    integer :: iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: TMP,TMPACC,TMPBCC
    real(realk) :: BasisCont(nPrimB,nPrimB)
    !Scaling p**4*c: nPrimA*nPrimB*nPrimC*nPrimD*nContC
    !$OMP SINGLE
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         TMPACC = ACC(iPrimA,iContC)
         do iPrimC=1,nPrimA
          TMP = TMP + TMPACC*ACC(iPrimC,iContC)*AUXarray2(iPrimC,iPrimD,iPrimA,iPrimB)
         enddo
        enddo
        BasisCont(iPrimD,iPrimB) = TMP
       enddo
      enddo
      do iContD=1,nContB
       TMP = 0.0E0_realk
       do iPrimB=1,nPrimB
        TMPBCC = BCC(iPrimB,iContD)
        do iPrimD=1,nPrimB
         TMP = TMP + TMPBCC*BCC(iPrimD,iContD)*BasisCont(iPrimD,iPrimB)
        enddo
       enddo
       AUXarrayCont(iContC,iContD) = TMP
      enddo
     enddo
    !$OMP END SINGLE
  end subroutine GabPrimitiveContractionGen1

  subroutine GabPrimitiveContractionGen16(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16,nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContA,nContB)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont(   16,nPrimB,nPrimB)
    real(realk) :: ACCTMP,BCCTMP
    !$OMP SINGLE
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,   16
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         BasisCont(iTUV,iPrimD,iPrimB) = TMP
        enddo
       enddo
      enddo
      do iContD=1,nContB
       do iTUV=1,   16
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         BCCTMP = BCC(iPrimB,iContD)
         do iPrimD=1,nPrimB
          TMP = TMP + BCC(iPrimD,iContD)*BCCTMP*BasisCont(iTUV,iPrimD,iPrimB)
         enddo
        enddo
        AUXarrayCont(iTUV,iContC,iContD) = TMP
       enddo
      enddo
     enddo
    !$OMP END SINGLE
  end subroutine GabPrimitiveContractionGen16

  subroutine GabPrimitiveContractionGen100(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100,nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContA,nContB)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont(  100,nPrimB,nPrimB)
    real(realk) :: ACCTMP,BCCTMP
    !$OMP SINGLE
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,  100
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         BasisCont(iTUV,iPrimD,iPrimB) = TMP
        enddo
       enddo
      enddo
      do iContD=1,nContB
       do iTUV=1,  100
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         BCCTMP = BCC(iPrimB,iContD)
         do iPrimD=1,nPrimB
          TMP = TMP + BCC(iPrimD,iContD)*BCCTMP*BasisCont(iTUV,iPrimD,iPrimB)
         enddo
        enddo
        AUXarrayCont(iTUV,iContC,iContD) = TMP
       enddo
      enddo
     enddo
    !$OMP END SINGLE
  end subroutine GabPrimitiveContractionGen100

  subroutine GabPrimitiveContractionGen400(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400,nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContA,nContB)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont(  400,nPrimB,nPrimB)
    real(realk) :: ACCTMP,BCCTMP
    !$OMP SINGLE
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,  400
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         BasisCont(iTUV,iPrimD,iPrimB) = TMP
        enddo
       enddo
      enddo
      do iContD=1,nContB
       do iTUV=1,  400
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         BCCTMP = BCC(iPrimB,iContD)
         do iPrimD=1,nPrimB
          TMP = TMP + BCC(iPrimD,iContD)*BCCTMP*BasisCont(iTUV,iPrimD,iPrimB)
         enddo
        enddo
        AUXarrayCont(iTUV,iContC,iContD) = TMP
       enddo
      enddo
     enddo
    !$OMP END SINGLE
  end subroutine GabPrimitiveContractionGen400

  subroutine GabPrimitiveContractionGen1225(AUXarray2,AUXarrayCont,nPrimP,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimA,nPrimB,nPrimA,nPrimB)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContA,nContB)
    !
    integer :: iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont( 1225,nPrimB,nPrimB)
    real(realk) :: ACCTMP,BCCTMP
    !$OMP SINGLE
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1, 1225
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB)
          enddo
         enddo
         BasisCont(iTUV,iPrimD,iPrimB) = TMP
        enddo
       enddo
      enddo
      do iContD=1,nContB
       do iTUV=1, 1225
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         BCCTMP = BCC(iPrimB,iContD)
         do iPrimD=1,nPrimB
          TMP = TMP + BCC(iPrimD,iContD)*BCCTMP*BasisCont(iTUV,iPrimD,iPrimB)
         enddo
        enddo
        AUXarrayCont(iTUV,iContC,iContD) = TMP
       enddo
      enddo
     enddo
    !$OMP END SINGLE
  end subroutine GabPrimitiveContractionGen1225

  subroutine ExtractGabElmP1Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP
    real(realk) :: TMP
    !$OMP SINGLE
     TMP = ABS(AUXarray(1))
     do iContP=2,nContP
      IF(ABS(AUXarray(iContP)).GT.TMP)THEN
       TMP = ABS(AUXarray(iContP))
      ENDIF
     enddo
     Output(1) = SQRT(TMP)
    !$OMP END SINGLE
  end subroutine ExtractGabElmP1Gen


  subroutine ExtractGabElmP3Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(    3,    3,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(    3)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    3
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    3
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
  end subroutine ExtractGabElmP3Gen

  subroutine ExtractGabElmP5Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(    5,    5,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(    5)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    5
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    5
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
  end subroutine ExtractGabElmP5Gen

  subroutine ExtractGabElmP9Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(    9,    9,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(    9)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    9
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    9
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
  end subroutine ExtractGabElmP9Gen

  subroutine ExtractGabElmP15Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(   15,   15,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(   15)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,   15
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,   15
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
  end subroutine ExtractGabElmP15Gen

  subroutine ExtractGabElmP25Gen(AUXarray,Output,nContP)
    implicit none
    integer,intent(in) :: nContP
    real(realk),intent(in) :: AUXarray(   25,   25,nContP)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: iContP,i
    real(realk) :: TMP(   25)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,   25
      TMP(i) = ABS(AUXarray(i,i,1))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,   25
       TMP(i) = ABS(AUXarray(i,i,iContP))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(1) = SQRT(TotalMaxValue)
    !$OMP END SINGLE
  end subroutine ExtractGabElmP25Gen
END MODULE IchorEriGabintegralOBSGeneralModGen
