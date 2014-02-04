MODULE IchorEriGabintegralOBSGeneralModGen
!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory
!Contains routines for General Contracted Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use IchorEriCoulombintegralOBSGeneralModGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODA
use AGC_CPU_OBS_VERTICALRECURRENCEMODB
use AGC_CPU_OBS_VERTICALRECURRENCEMODD
use AGC_CPU_OBS_VERTICALRECURRENCEMODC
use AGC_OBS_TRANSFERRECURRENCEMODAtoCGen
use AGC_OBS_TRANSFERRECURRENCEMODAtoDGen
use AGC_OBS_TRANSFERRECURRENCEMODBtoCGen
use AGC_OBS_TRANSFERRECURRENCEMODBtoDGen
use AGC_OBS_TRANSFERRECURRENCEMODCtoAGen
use AGC_OBS_TRANSFERRECURRENCEMODDtoAGen
use AGC_OBS_TRANSFERRECURRENCEMODCtoBGen
use AGC_OBS_TRANSFERRECURRENCEMODDtoBGen
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
       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContP,pexp,ACC,BCC,&
       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,Pdistance12,PQorder,LOCALINTS,Acenter,Bcenter,&
       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisContmaxsize,BasisCont)
    implicit none
    integer,intent(in) :: nPrimP,nPasses,nPrimA,nPrimB
    integer,intent(in) :: MaxPasses,IntPrint,lupri
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
    real(realk),intent(inout) :: LOCALINTS(nPasses)
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
    integer :: AngmomP,I,J,nContQP,la,lb,lc,ld,nsize,angmomid
    
    !Setup combined Angmom info
    AngmomP = AngmomA+AngmomB
!    nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
!    nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
!    nlmA = 2*AngmomA+1
!    nlmB = 2*AngmomB+1
    AngmomID = 10*AngmomA+AngmomB
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*1.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU0(nPasses,nPrimP,nPrimP,&
               & reducedExponents,TABFJW,Pcent,Pcent,integralPrefactor,&
               & PpreExpFac,PpreExpFac,TMParray2(1:nPrimP*nPrimP*nPasses*1))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*1.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen1(TMParray2(1:nPrimP*nPrimP*nPasses*1),&
             & TMParray1(1:nContP*nPasses*1),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        call ExtractGabElmP1Gen(TMParray1(1:nContP*nPasses*1),&
            & LOCALINTS,nContP,nPasses)
    CASE(  10)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU2A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1AtoCGen(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*10),TMParray1(1:nPrimP*nPrimP*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen16(TMParray1(1:nPrimP*nPrimP*nPasses*16),&
             & TMParray2(1:nContP*nPasses*16),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*16),&
            & TMParray1(1:nContP*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,3,Pdistance12,TMParray1(1:nContP*nPasses*12),&
            & TMParray2(1:nContP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP3Gen(TMParray2(1:nContP*nPasses*9),&
            & LOCALINTS,nContP,nPasses)
    CASE(  11)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCGen(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*35),TMParray1(1:nPrimP*nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen100(TMParray1(1:nPrimP*nPrimP*nPasses*100),&
             & TMParray2(1:nContP*nPasses*100),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*90),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*81.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,9,Pdistance12,TMParray1(1:nContP*nPasses*90),&
            & TMParray2(1:nContP*nPasses*81),lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP9Gen(TMParray2(1:nContP*nPasses*81),&
            & LOCALINTS,nContP,nPasses)
    CASE(  20)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCGen(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*35),TMParray1(1:nPrimP*nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen100(TMParray1(1:nPrimP*nPrimP*nPasses*100),&
             & TMParray2(1:nContP*nPasses*100),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*60),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContP*nPasses,TMParray1(1:nContP*nPasses*60),&
            & TMParray2(1:nContP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,5,Pdistance12,TMParray2(1:nContP*nPasses*50),&
            & TMParray1(1:nContP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContP*nPasses,TMParray1(1:nContP*nPasses*30),&
            & TMParray2(1:nContP*nPasses*25))
        call ExtractGabElmP5Gen(TMParray2(1:nContP*nPasses*25),&
            & LOCALINTS,nContP,nPasses)
    CASE(  21)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3AtoCGen(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*84),TMParray1(1:nPrimP*nPrimP*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen400(TMParray1(1:nPrimP*nPrimP*nPasses*400),&
             & TMParray2(1:nContP*nPasses*400),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*400),&
            & TMParray1(1:nContP*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nContP*nPasses,TMParray1(1:nContP*nPasses*360),&
            & TMParray2(1:nContP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,15,Pdistance12,TMParray2(1:nContP*nPasses*300),&
            & TMParray1(1:nContP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContP*nPasses,TMParray1(1:nContP*nPasses*270),&
            & TMParray2(1:nContP*nPasses*225))
        call ExtractGabElmP15Gen(TMParray2(1:nContP*nPasses*225),&
            & LOCALINTS,nContP,nPasses)
    CASE(  22)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*165.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU8A(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*165))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q4AtoCGen(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*165),TMParray1(1:nPrimP*nPrimP*nPasses*1225))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen1225(TMParray1(1:nPrimP*nPrimP*nPasses*1225),&
             & TMParray2(1:nContP*nPasses*1225),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*1260.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,35,Pdistance12,TMParray2(1:nContP*nPasses*1225),&
            & TMParray1(1:nContP*nPasses*1260),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*875.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,nContP*nPasses,TMParray1(1:nContP*nPasses*1260),&
            & TMParray2(1:nContP*nPasses*875))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*900.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,25,Pdistance12,TMParray2(1:nContP*nPasses*875),&
            & TMParray1(1:nContP*nPasses*900),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*625.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nContP*nPasses,TMParray1(1:nContP*nPasses*900),&
            & TMParray2(1:nContP*nPasses*625))
        call ExtractGabElmP25Gen(TMParray2(1:nContP*nPasses*625),&
            & LOCALINTS,nContP,nPasses)
    CASE(   1)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU2B(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1BtoDGen(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*10),TMParray1(1:nPrimP*nPrimP*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen16(TMParray1(1:nPrimP*nPrimP*nPasses*16),&
             & TMParray2(1:nContP*nPasses*16),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*16),&
            & TMParray1(1:nContP*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,3,Pdistance12,TMParray1(1:nContP*nPasses*12),&
            & TMParray2(1:nContP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP3Gen(TMParray2(1:nContP*nPasses*9),&
            & LOCALINTS,nContP,nPasses)
    CASE(   2)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2BtoDGen(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*35),TMParray1(1:nPrimP*nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen100(TMParray1(1:nPrimP*nPrimP*nPasses*100),&
             & TMParray2(1:nContP*nPasses*100),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*60),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContP*nPasses,TMParray1(1:nContP*nPasses*60),&
            & TMParray2(1:nContP*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,5,Pdistance12,TMParray2(1:nContP*nPasses*50),&
            & TMParray1(1:nContP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContP*nPasses,TMParray1(1:nContP*nPasses*30),&
            & TMParray2(1:nContP*nPasses*25))
        call ExtractGabElmP5Gen(TMParray2(1:nContP*nPasses*25),&
            & LOCALINTS,nContP,nPasses)
    CASE(  12)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6B(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3BtoDGen(nPasses,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & TMParray2(1:nPrimP*nPrimP*nPasses*84),TMParray1(1:nPrimP*nPrimP*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call GabPrimitiveContractionGen400(TMParray1(1:nPrimP*nPrimP*nPasses*400),&
             & TMParray2(1:nContP*nPasses*400),nPrimP,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*400),&
            & TMParray1(1:nContP*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nContP*nPasses,TMParray1(1:nContP*nPasses*360),&
            & TMParray2(1:nContP*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,15,Pdistance12,TMParray2(1:nContP*nPasses*300),&
            & TMParray1(1:nContP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContP*nPasses,TMParray1(1:nContP*nPasses*270),&
            & TMParray2(1:nContP*nPasses*225))
        call ExtractGabElmP15Gen(TMParray2(1:nContP*nPasses*225),&
            & LOCALINTS,nContP,nPasses)
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
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContP)
    CASE(   2)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       BasisContmaxsize =   100*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,25*nContP)
    CASE(  10)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       BasisContmaxsize =    16*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContP)
    CASE(  11)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       BasisContmaxsize =   100*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,81*nContP)
    CASE(  12)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       BasisContmaxsize =   400*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,225*nContP)
    CASE(  20)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       BasisContmaxsize =   100*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,25*nContP)
    CASE(  21)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       BasisContmaxsize =   400*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,225*nContP)
    CASE(  22)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       BasisContmaxsize =  1225*nPrimB*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,165*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,625*nContP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorGabIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorGabIntegral_OBS_general_sizeGen
  subroutine GabPrimitiveContractionGen1(AUXarray2,AUXarrayCont,nPrimP,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: TMP,TMPACC,TMPBCC
    real(realk) :: BasisCont(nPrimB,nPrimB)
    !Scaling p**4*c*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nPassQ
    do iPassQ = 1,nPasses
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         TMPACC = ACC(iPrimA,iContC)
         do iPrimC=1,nPrimA
          TMP = TMP + TMPACC*ACC(iPrimC,iContC)*AUXarray2(iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
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
       AUXarrayCont(iContC,iContD,iPassQ) = TMP
      enddo
     enddo
    enddo
  end subroutine GabPrimitiveContractionGen1

  subroutine GabPrimitiveContractionGen16(AUXarray2,AUXarrayCont,nPrimP,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16,nPrimA,nPrimB,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont(   16,nPrimB,nPrimB)
    real(realk) :: ACCTMP,BCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,   16
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
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
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = TMP
       enddo
      enddo
     enddo
    enddo
  end subroutine GabPrimitiveContractionGen16

  subroutine GabPrimitiveContractionGen100(AUXarray2,AUXarrayCont,nPrimP,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100,nPrimA,nPrimB,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont(  100,nPrimB,nPrimB)
    real(realk) :: ACCTMP,BCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,  100
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
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
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = TMP
       enddo
      enddo
     enddo
    enddo
  end subroutine GabPrimitiveContractionGen100

  subroutine GabPrimitiveContractionGen400(AUXarray2,AUXarrayCont,nPrimP,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400,nPrimA,nPrimB,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont(  400,nPrimB,nPrimB)
    real(realk) :: ACCTMP,BCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1,  400
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
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
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = TMP
       enddo
      enddo
     enddo
    enddo
  end subroutine GabPrimitiveContractionGen400

  subroutine GabPrimitiveContractionGen1225(AUXarray2,AUXarrayCont,nPrimP,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimA,nPrimB,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont( 1225,nPrimB,nPrimB)
    real(realk) :: ACCTMP,BCCTMP
    do iPassQ = 1,nPasses
     do iContC=1,nContA
      do iPrimB=1,nPrimB
       do iPrimD=1,nPrimB
        do iTUV=1, 1225
         TMP = 0.0E0_realk
         do iPrimA=1,nPrimA
          ACCTMP = ACC(iPrimA,iContC)
          do iPrimC=1,nPrimA
           TMP = TMP + ACC(iPrimC,iContC)*ACCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)
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
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = TMP
       enddo
      enddo
     enddo
    enddo
  end subroutine GabPrimitiveContractionGen1225

  subroutine ExtractGabElmP1Gen(AUXarray,Output,nContP,nPasses)
    implicit none
    integer,intent(in) :: nContP,nPasses
    real(realk),intent(in) :: AUXarray(nContP,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,iContP
    real(realk) :: TMP
    do iPassQ = 1,nPasses
     TMP = ABS(AUXarray(1,iPassQ))
     do iContP=2,nContP
      IF(ABS(AUXarray(iContP,iPassQ)).GT.TMP)THEN
       TMP = ABS(AUXarray(iContP,iPassQ))
      ENDIF
     enddo
     Output(iPassQ) = SQRT(TMP)
    enddo
  end subroutine ExtractGabElmP1Gen


  subroutine ExtractGabElmP3Gen(AUXarray,Output,nContP,nPasses)
    implicit none
    integer,intent(in) :: nContP,nPasses
    real(realk),intent(in) :: AUXarray(    3,    3,nContP,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,iContP,i
    real(realk) :: TMP(    3)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,    3
      TMP(i) = ABS(AUXarray(i,i,1,iPassQ))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    3
       TMP(i) = ABS(AUXarray(i,i,iContP,iPassQ))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(iPassQ) = SQRT(TotalMaxValue)
    enddo
  end subroutine ExtractGabElmP3Gen

  subroutine ExtractGabElmP5Gen(AUXarray,Output,nContP,nPasses)
    implicit none
    integer,intent(in) :: nContP,nPasses
    real(realk),intent(in) :: AUXarray(    5,    5,nContP,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,iContP,i
    real(realk) :: TMP(    5)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,    5
      TMP(i) = ABS(AUXarray(i,i,1,iPassQ))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    5
       TMP(i) = ABS(AUXarray(i,i,iContP,iPassQ))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(iPassQ) = SQRT(TotalMaxValue)
    enddo
  end subroutine ExtractGabElmP5Gen

  subroutine ExtractGabElmP9Gen(AUXarray,Output,nContP,nPasses)
    implicit none
    integer,intent(in) :: nContP,nPasses
    real(realk),intent(in) :: AUXarray(    9,    9,nContP,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,iContP,i
    real(realk) :: TMP(    9)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,    9
      TMP(i) = ABS(AUXarray(i,i,1,iPassQ))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,    9
       TMP(i) = ABS(AUXarray(i,i,iContP,iPassQ))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(iPassQ) = SQRT(TotalMaxValue)
    enddo
  end subroutine ExtractGabElmP9Gen

  subroutine ExtractGabElmP15Gen(AUXarray,Output,nContP,nPasses)
    implicit none
    integer,intent(in) :: nContP,nPasses
    real(realk),intent(in) :: AUXarray(   15,   15,nContP,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,iContP,i
    real(realk) :: TMP(   15)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,   15
      TMP(i) = ABS(AUXarray(i,i,1,iPassQ))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,   15
       TMP(i) = ABS(AUXarray(i,i,iContP,iPassQ))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(iPassQ) = SQRT(TotalMaxValue)
    enddo
  end subroutine ExtractGabElmP15Gen

  subroutine ExtractGabElmP25Gen(AUXarray,Output,nContP,nPasses)
    implicit none
    integer,intent(in) :: nContP,nPasses
    real(realk),intent(in) :: AUXarray(   25,   25,nContP,nPasses)
    real(realk),intent(inout) :: Output(nPasses)
    !
    integer :: iPassQ,iContP,i
    real(realk) :: TMP(   25)
    real(realk) :: MaxValue,TotalMaxValue
    do iPassQ = 1,nPasses
     do i=1,   25
      TMP(i) = ABS(AUXarray(i,i,1,iPassQ))
     enddo
     TotalMaxValue = MAXVAL(TMP)
     do iContP=2,nContP
      do i=1,   25
       TMP(i) = ABS(AUXarray(i,i,iContP,iPassQ))
      enddo
      maxvalue = MAXVAL(TMP)
      IF(MaxValue.GT.TotalMaxValue)THEN
       TotalMaxValue = MaxValue
      ENDIF
     enddo
     Output(iPassQ) = SQRT(TotalMaxValue)
    enddo
  end subroutine ExtractGabElmP25Gen
END MODULE IchorEriGabintegralOBSGeneralModGen
