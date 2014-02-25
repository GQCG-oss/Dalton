MODULE IchorEriGabintegralOBSGeneralModSeg
!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000MODGen
use AGC_CPU_OBS_BUILDRJ000MODSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODASeg
use AGC_CPU_OBS_VERTICALRECURRENCEMODBSeg
use AGC_CPU_OBS_VERTICALRECURRENCEMODDSeg
use AGC_CPU_OBS_VERTICALRECURRENCEMODCSeg
use AGC_CPU_OBS_TRMODAtoCSeg
use AGC_CPU_OBS_TRMODAtoDSeg
use AGC_CPU_OBS_TRMODBtoCSeg
use AGC_CPU_OBS_TRMODBtoDSeg
use AGC_CPU_OBS_TRMODCtoASeg
use AGC_CPU_OBS_TRMODDtoASeg
use AGC_CPU_OBS_TRMODCtoBSeg
use AGC_CPU_OBS_TRMODDtoBSeg
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorGabIntegral_OBS_Seg,IchorGabIntegral_OBS_general_sizeSeg  
  
CONTAINS
  
  
  subroutine IchorGabIntegral_OBS_Seg(nPrimA,nPrimB,&
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
        IF(1*1.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg0(1,nPrimP,nPrimP,&
               & reducedExponents,TABFJW,Pcent,Pcent,integralPrefactor,&
               & IatomApass,IatomBpass,1,1,1,PpreExpFac,PpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        call ExtractGabElmP1Seg(TMParray2,LOCALINTS)
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
        IF(1*16.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*12.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(1,1,4,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(1*9.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(1,1,3,Pdistance12,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP3Seg(TMParray2,LOCALINTS)
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
        IF(1*100.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*90.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(1,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(1*81.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(1,1,9,Pdistance12,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP9Seg(TMParray2,LOCALINTS)
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
        IF(1*100.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*60.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(1,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*50.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*30.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(1,1,5,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*25.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,1,TMParray1,&
            & TMParray2)
        call ExtractGabElmP5Seg(TMParray2,LOCALINTS)
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
        IF(1*400.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*360.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(1,1,20,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*300.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*270.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(1,1,15,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*225.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,1,TMParray1,&
            & TMParray2)
        call ExtractGabElmP15Seg(TMParray2,LOCALINTS)
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
        IF(1*1225.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q4AtoCSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*1260.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(1,1,35,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*875.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*900.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(1,1,25,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*625.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,1,TMParray1,&
            & TMParray2)
        call ExtractGabElmP25Seg(TMParray2,LOCALINTS)
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
        IF(1*16.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoDSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*12.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(1,1,4,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(1*9.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(1,1,3,Pdistance12,TMParray1,&
            & TMParray2,lupri)
        !no Spherical Transformation RHS needed
        call ExtractGabElmP3Seg(TMParray2,LOCALINTS)
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
        IF(1*100.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoDSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*60.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(1,1,10,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*50.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*30.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(1,1,5,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*25.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,1,TMParray1,&
            & TMParray2)
        call ExtractGabElmP5Seg(TMParray2,LOCALINTS)
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
        IF(1*400.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoDSeg(1,nPrimP,nPrimP,reducedExponents,&
               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&
               & 1,1,1,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(1*360.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(1,1,20,Pdistance12,1,1,1,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*300.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,1,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(1*270.GT.TMParray1maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(1,1,15,Pdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(1*225.GT.TMParray2maxsize)THEN
          call ichorquit('1too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,1,TMParray1,&
            & TMParray2)
        call ExtractGabElmP15Seg(TMParray2,LOCALINTS)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorGabIntegral_OBS_Seg',-1)
    END SELECT
  end subroutine IchorGabIntegral_OBS_Seg
  
  subroutine IchorGabIntegral_OBS_general_sizeSeg(TMParray1maxsize,&
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
        CALL ICHORQUIT('Unknown Case in IchorGabIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorGabIntegral_OBS_general_sizeSeg
  
  subroutine ExtractGabElmP1Seg(AUXarray,Output)
    implicit none
    real(realk),intent(in) :: AUXarray(1)
    real(realk),intent(inout) :: Output(1)
    !$OMP SINGLE
     Output(1) = SQRT(ABS(AUXarray(1)))
    !$OMP END SINGLE
  end subroutine ExtractGabElmP1Seg

  subroutine ExtractGabElmP3Seg(AUXarray,Output)
    implicit none
    real(realk),intent(in) :: AUXarray(    3,    3)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: i
    real(realk) :: TMP(    3)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    3
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
  end subroutine ExtractGabElmP3Seg

  subroutine ExtractGabElmP5Seg(AUXarray,Output)
    implicit none
    real(realk),intent(in) :: AUXarray(    5,    5)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: i
    real(realk) :: TMP(    5)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    5
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
  end subroutine ExtractGabElmP5Seg

  subroutine ExtractGabElmP9Seg(AUXarray,Output)
    implicit none
    real(realk),intent(in) :: AUXarray(    9,    9)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: i
    real(realk) :: TMP(    9)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,    9
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
  end subroutine ExtractGabElmP9Seg

  subroutine ExtractGabElmP15Seg(AUXarray,Output)
    implicit none
    real(realk),intent(in) :: AUXarray(   15,   15)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: i
    real(realk) :: TMP(   15)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,   15
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
  end subroutine ExtractGabElmP15Seg

  subroutine ExtractGabElmP25Seg(AUXarray,Output)
    implicit none
    real(realk),intent(in) :: AUXarray(   25,   25)
    real(realk),intent(inout) :: Output(1)
    !
    integer :: i
    real(realk) :: TMP(   25)
    real(realk) :: MaxValue,TotalMaxValue
    !$OMP SINGLE
     do i=1,   25
      TMP(i) = ABS(AUXarray(i,i))
     enddo
     Output(1) = SQRT(MAXVAL(TMP))
    !$OMP END SINGLE
  end subroutine ExtractGabElmP25Seg
END MODULE IchorEriGabintegralOBSGeneralModSeg
