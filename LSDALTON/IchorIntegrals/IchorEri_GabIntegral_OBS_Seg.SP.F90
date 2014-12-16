module SPIchorEriGabintegralOBSGeneralModSeg
!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset 
use IchorprecisionMod
use IchorCommonMod
use IchorMemory
use SPAGC_CPU_OBS_BUILDRJ000MODGen
use SPAGC_CPU_OBS_BUILDRJ000MODSeg1Prim
use IchorEriGabintegralCPUMcMGeneralMod
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODASeg
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODBSeg
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODDSeg
use SPAGC_CPU_OBS_VERTICALRECURRENCEMODCSeg
use SPAGC_CPU_OBS_TRMODAtoCSeg1
use SPAGC_CPU_OBS_TRMODAtoCSeg2
use SPAGC_CPU_OBS_TRMODAtoCSeg3
use SPAGC_CPU_OBS_TRMODAtoDSeg1
use SPAGC_CPU_OBS_TRMODAtoDSeg2
use SPAGC_CPU_OBS_TRMODBtoCSeg1
use SPAGC_CPU_OBS_TRMODBtoDSeg1
use SPAGC_CPU_OBS_TRMODCtoASeg
use SPAGC_CPU_OBS_TRMODDtoASeg
use SPAGC_CPU_OBS_TRMODCtoBSeg
use SPAGC_CPU_OBS_TRMODDtoBSeg
use SPAGC_CPU_OBS_HorizontalRecurrenceLHSModAtoB
use SPAGC_CPU_OBS_HorizontalRecurrenceLHSModBtoA
use SPAGC_CPU_OBS_HorizontalRecurrenceRHSModCtoD
use SPAGC_CPU_OBS_HorizontalRecurrenceRHSModDtoC
use SPAGC_CPU_OBS_Sphcontract1Mod
use SPAGC_CPU_OBS_Sphcontract2Mod
  
private   
public :: IGI_OBS_Seg,IGI_OBS_general_sizeSeg  
  
CONTAINS
  
  
  subroutine SPIGI_OBS_Seg(nPrimA,nPrimB,&
       & nPrimP,IntPrint,lupri,&
       & nContA,nContB,nContP,pexp,ACC,BCC,&
       & nOrbCompA,nOrbCompB,nCartOrbCompA,nCartOrbCompB,&
       & nCartOrbCompP,nOrbCompP,nTUVP,nTUV,&
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
    integer,intent(in) :: nOrbCompA,nOrbCompB,nCartOrbCompA,nCartOrbCompB
    integer,intent(in) :: nCartOrbCompP,nOrbCompP,nTUVP,nTUV
    real(reals),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)
    logical,intent(in)     :: Psegmented
    real(reals),intent(in) :: pexp(nPrimP)
    real(reals),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(reals),intent(in) :: PpreExpFac(nPrimP)
    real(reals),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(reals),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(reals) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(reals),intent(inout) :: LOCALINTS(1)
    real(reals),intent(in) :: integralPrefactor(nPrimP*nPrimP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimP)
    real(reals),intent(in) :: reducedExponents(nPrimP*nPrimP)
    !reducedExponents(nPrimP,nPrimP)
    real(reals),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(reals),intent(in) :: Acenter(3),Bcenter(3)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize,BasisContmaxsize
!   TMP variables - allocated outside
    real(reals),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    real(reals),intent(inout) :: BasisCont(BasisContmaxsize)
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
    IF(UseGeneralCode) AngmomID = AngmomID + 10000 !force to use general code
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(1*1.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUSeg0(1,nPrimP,nPrimP,&
               & reducedExponents,TABFJW,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        call SPExtractGabElmP1Seg(TMParray2,LOCALINTS)
    CASE(  10)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPBuildRJ000CPUGen2(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUGen2A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(1*16.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPTransferRecurrenceCPUP1Q1AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*12.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_LHS_P1A1B0AtoB(1,1,4,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(1*9.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_RHS_Q1C1D0CtoD(1,1,3,Pdistance12,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation RHS needed
        call SPExtractGabElmP3Seg(TMParray2,LOCALINTS)
    CASE(  11)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPBuildRJ000CPUGen4(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUGen4A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(1*100.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPTransferRecurrenceCPUP2Q2AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*90.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_LHS_P2A1B1AtoB(1,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(1*81.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_RHS_Q2C1D1CtoD(1,1,9,Pdistance12,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation RHS needed
        call SPExtractGabElmP9Seg(TMParray2,LOCALINTS)
    CASE(  20)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPBuildRJ000CPUGen4(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUGen4A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(1*100.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPTransferRecurrenceCPUP2Q2AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*60.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_LHS_P2A2B0AtoB(1,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*50.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*30.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_RHS_Q2C2D0CtoD(1,1,5,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*25.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,1,TMParray1,&
            & TMParray2)
        call SPExtractGabElmP5Seg(TMParray2,LOCALINTS)
    CASE(  21)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPBuildRJ000CPUGen6(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUGen6A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(1*400.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPTransferRecurrenceCPUP3Q3AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*360.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_LHS_P3A2B1AtoB(1,1,20,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*300.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*270.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_RHS_Q3C2D1CtoD(1,1,15,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*225.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,1,TMParray1,&
            & TMParray2)
        call SPExtractGabElmP15Seg(TMParray2,LOCALINTS)
    CASE(  22)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*9.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPBuildRJ000CPUGen8(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*165.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUGen8A(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(1*1225.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPTransferRecurrenceCPUP4Q4AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*1260.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_LHS_P4A2B2AtoB(1,1,35,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*875.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2(35,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*900.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_RHS_Q4C2D2CtoD(1,1,25,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*625.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS2_CPU_maxAngQ4_maxAngC2(25,1,TMParray1,&
            & TMParray2)
        call SPExtractGabElmP25Seg(TMParray2,LOCALINTS)
    CASE(   1)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPBuildRJ000CPUGen2(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUGen2B(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(1*16.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPTransferRecurrenceCPUP1Q1BtoDSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*12.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_LHS_P1A0B1BtoA(1,1,4,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(1*9.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_RHS_Q1C0D1DtoC(1,1,3,Pdistance12,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation RHS needed
        call SPExtractGabElmP3Seg(TMParray2,LOCALINTS)
    CASE(   2)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPBuildRJ000CPUGen4(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUGen4B(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(1*100.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPTransferRecurrenceCPUP2Q2BtoDSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*60.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_LHS_P2A0B2BtoA(1,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*50.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*30.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_RHS_Q2C0D2DtoC(1,1,5,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*25.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,1,TMParray1,&
            & TMParray2)
        call SPExtractGabElmP5Seg(TMParray2,LOCALINTS)
    CASE(  12)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPBuildRJ000CPUGen6(1,nPrimP,nPrimP,reducedExponents,&
               & TABFJW,Pcent,Pcent,IatomApass,IatomBpass,&
               & 1,1,1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimP*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimPtoo small',-1)
        ENDIF
#endif
        call SPVerticalRecurrenceCPUGen6B(1,nPrimP,nPrimP,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(1*400.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPTransferRecurrenceCPUP3Q3BtoDSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*360.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_LHS_P3A1B2BtoA(1,1,20,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*300.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*270.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPHorizontalRR_CPU_RHS_Q3C1D2DtoC(1,1,15,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*225.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SPSphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,1,TMParray1,&
            & TMParray2)
        call SPExtractGabElmP15Seg(TMParray2,LOCALINTS)
    CASE DEFAULT
        call SPIGI_CPU_McM_general(nPrimA,nPrimB,&
           & nPrimP,IntPrint,lupri,&
           & nContA,nContB,nContP,pexp,ACC,BCC,&
           & nOrbCompA,nOrbCompB,nCartOrbCompA,nCartOrbCompB,&
           & nCartOrbCompP,nOrbCompP,nTUVP,nTUV,&
           & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
           & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&
           & AngmomA,AngmomB,Pdistance12,PQorder,LOCALINTS,&
           & Acenter,Bcenter,spherical,&
           & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    END SELECT
  end subroutine SPIGI_OBS_Seg
  
  subroutine SPIGI_OBS_general_sizeSeg(TMParray1maxsize,&
         &TMParray2maxsize,BasisContmaxsize,AngmomA,AngmomB,nPrimP,nContP,nPrimB)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize,BasisContmaxsize
    integer,intent(in) :: AngmomA,AngmomB
    integer,intent(in) :: nPrimP,nContP,nPrimB
    ! local variables
    integer :: AngmomID
    
    AngmomID = 10*AngmomA+AngmomB
    IF(UseGeneralCode) AngmomID = AngmomID + 10000 !force to use general code
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    BasisContmaxsize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,1)
    CASE(   1)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
       TMParray2maxSize = MAX(TMParray2maxSize,9)
    CASE(   2)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,25)
    CASE(  10)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
       TMParray2maxSize = MAX(TMParray2maxSize,9)
    CASE(  11)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
       TMParray2maxSize = MAX(TMParray2maxSize,81)
    CASE(  12)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
       TMParray2maxSize = MAX(TMParray2maxSize,225)
    CASE(  20)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,25)
    CASE(  21)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
       TMParray2maxSize = MAX(TMParray2maxSize,225)
    CASE(  22)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,9*nPrimP*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,165*nPrimP*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225)
       TMParray1maxSize = MAX(TMParray1maxSize,1260)
       TMParray2maxSize = MAX(TMParray2maxSize,875)
       TMParray1maxSize = MAX(TMParray1maxSize,900)
       TMParray2maxSize = MAX(TMParray2maxSize,625)
    CASE DEFAULT
      call SPIGI_CPU_McM_general_size(TMParray1maxsize,&
          & TMParray2maxsize,AngmomA,AngmomB,&
          & nPrimP,nContP,nPrimB,.TRUE.)
    END SELECT
  end subroutine SPIGI_OBS_general_sizeSeg
  
  subroutine SPExtractGabElmP1Seg(AUXarray,Output)
    implicit none
    real(reals),intent(in) :: AUXarray(1)
    real(reals),intent(inout) :: Output(1)
    !$OMP SINGLE
     Output(1) = SQRT(ABS(AUXarray(1)))
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine SPExtractGabElmP1Seg

  subroutine SPExtractGabElmP3Seg(AUXarray,Output)
    implicit none
    real(reals),intent(in) :: AUXarray(    3,    3)
    real(reals),intent(inout) :: Output(1)
    !
    integer :: i
    real(reals) :: TMP(    3)
    real(reals) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    3
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine SPExtractGabElmP3Seg

  subroutine SPExtractGabElmP5Seg(AUXarray,Output)
    implicit none
    real(reals),intent(in) :: AUXarray(    5,    5)
    real(reals),intent(inout) :: Output(1)
    !
    integer :: i
    real(reals) :: TMP(    5)
    real(reals) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    5
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine SPExtractGabElmP5Seg

  subroutine SPExtractGabElmP9Seg(AUXarray,Output)
    implicit none
    real(reals),intent(in) :: AUXarray(    9,    9)
    real(reals),intent(inout) :: Output(1)
    !
    integer :: i
    real(reals) :: TMP(    9)
    real(reals) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    9
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine SPExtractGabElmP9Seg

  subroutine SPExtractGabElmP15Seg(AUXarray,Output)
    implicit none
    real(reals),intent(in) :: AUXarray(   15,   15)
    real(reals),intent(inout) :: Output(1)
    !
    integer :: i
    real(reals) :: TMP(   15)
    real(reals) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,   15
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine SPExtractGabElmP15Seg

  subroutine SPExtractGabElmP25Seg(AUXarray,Output)
    implicit none
    real(reals),intent(in) :: AUXarray(   25,   25)
    real(reals),intent(inout) :: Output(1)
    !
    integer :: i
    real(reals) :: TMP(   25)
    real(reals) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,   25
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
    !$OMP BARRIER
  end subroutine SPExtractGabElmP25Seg
END module SPIchorEriGabintegralOBSGeneralModSeg
