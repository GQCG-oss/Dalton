MODULE IchorEriCoulombintegralCPUOBSGeneralModSegQ
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for a General Contracted LHS Segmented contracted RHS and Basisset 
use IchorEriCoulombintegralCPUMcMGeneralMod
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000ModGen
use AGC_CPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODASegQ
use AGC_CPU_OBS_VERTICALRECURRENCEMODBSegQ
use AGC_CPU_OBS_VERTICALRECURRENCEMODDSegQ
use AGC_CPU_OBS_VERTICALRECURRENCEMODCSegQ
use AGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_CPU_OBS_TRMODAtoCSegQ
use AGC_CPU_OBS_TRMODAtoDSegQ
use AGC_CPU_OBS_TRMODBtoCSegQ
use AGC_CPU_OBS_TRMODBtoDSegQ
use AGC_CPU_OBS_TRMODCtoASegQ
use AGC_CPU_OBS_TRMODDtoASegQ
use AGC_CPU_OBS_TRMODCtoBSegQ
use AGC_CPU_OBS_TRMODDtoBSegQ
use AGC_CPU_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_CPU_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_CPU_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_CPU_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_CPU_OBS_Sphcontract1Mod
use AGC_CPU_OBS_Sphcontract2Mod
  
private   
public :: ICI_CPU_OBS_SegQ,ICI_CPU_OBS_general_sizeSegQ  
  
CONTAINS
  
  
  subroutine ICI_CPU_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    IF(UseGeneralCode) AngmomID = AngmomID + 10000 !force to use general code
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*1.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*1.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ1(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
         call PrimitiveContractionBCPUSegQ1(TMParray1,LOCALINTS,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,1,&
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
        IF(nPrimP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,1,&
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*81.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(35,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(35,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(25,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q3AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(25,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q4AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*1260.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*875.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(35,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*900.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*625.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(25,nContP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,1,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,1,&
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
        IF(nPrimP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(35,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(3,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(9,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(9,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(9,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(9,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*315.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*324.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(9,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(1,nContP*nPasses,TMParray2,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContP,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(35,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(5,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContP,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
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
        IF(nPrimP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(25,nContP*nPasses,TMParray1,&
            & LOCALINTS)
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
        IF(nPrimP*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q3AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nPrimB*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContA*nPrimB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionACPUSegQ700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContA*nContB*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContA*nContB*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionBCPUSegQ700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContP,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nContP*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContP,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(25,nContP*nPasses,TMParray1,&
            & LOCALINTS)
    CASE DEFAULT
        call ICI_CPU_McM_general(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    END SELECT
  end subroutine ICI_CPU_OBS_SegQ
  
  
  subroutine ICI_CPU_OBS_general_sizeSegQ(TMParray1maxsize,&
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
    IF(UseGeneralCode) AngmomID = AngmomID + 10000 !force to use general code
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,1*nContA*nPrimB)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContA*nContB)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContA*nContB)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContA*nContB)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,35*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,35*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContP)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContA*nContB)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContP)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContA*nContB)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContP)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContA*nContB)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContP)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContP)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContP)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContP)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,315*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,324*nContP)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContP)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContP)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimP)
    TMParray2maxSize = MAX(TMParray2maxSize,35*nContA*nPrimB)
    TMParray1maxSize = MAX(TMParray1maxSize,35*nContA*nContB)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContP)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContP)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContP)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContP)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContP)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContP)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,9*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,165*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nPrimP)
    TMParray1maxSize = MAX(TMParray1maxSize,1225*nContA*nPrimB)
    TMParray2maxSize = MAX(TMParray2maxSize,1225*nContA*nContB)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContP)
    CASE DEFAULT
     call ICI_CPU_McM_general_size(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP,&
         & .FALSE.,.TRUE.)
    END SELECT
  end subroutine ICI_CPU_OBS_general_sizeSegQ

   subroutine PrimitiveContractionACPUSegQ1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     TMP = 0.0E0_realk
     do iPrimA=1,nPrimA
      TMP = TMP + ACC(iPrimA,iContA)*AUXarray2(iPrimA,iP)
     enddo
     AUXarrayCont(iContA,iP) = tmp
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ1


   subroutine PrimitiveContractionBCPUSegQ1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iContA
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iContA,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iContA=1,nContA
       AUXarrayCont(iContA,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iContA=1,nContA
       AUXarrayCont(iContA,iContB,iPassP)=AUXarrayCont(iContA,iContB,iPassP)+tmp*AUXarray2(iContA,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ1


   subroutine PrimitiveContractionACPUSegQ4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(    4,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,    4
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,    4
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ4


   subroutine PrimitiveContractionBCPUSegQ4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(    4*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,    4*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,    4*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ4


   subroutine PrimitiveContractionACPUSegQ10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   10,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   10
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   10
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ10


   subroutine PrimitiveContractionBCPUSegQ10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   10*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   10*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   10*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ10


   subroutine PrimitiveContractionACPUSegQ20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   20,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   20
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   20
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ20


   subroutine PrimitiveContractionBCPUSegQ20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   20*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   20*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   20*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ20


   subroutine PrimitiveContractionACPUSegQ35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   35,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   35
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   35
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ35


   subroutine PrimitiveContractionBCPUSegQ35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   35*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   35*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   35*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ35


   subroutine PrimitiveContractionACPUSegQ16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   16,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   16
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   16
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ16


   subroutine PrimitiveContractionBCPUSegQ16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   16*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   16*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ16


   subroutine PrimitiveContractionACPUSegQ40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   40,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   40
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   40
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ40


   subroutine PrimitiveContractionBCPUSegQ40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   40*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   40*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   40*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ40


   subroutine PrimitiveContractionACPUSegQ80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(   80,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,   80
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,   80
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ80


   subroutine PrimitiveContractionBCPUSegQ80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   80*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,   80*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,   80*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ80


   subroutine PrimitiveContractionACPUSegQ140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  140,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  140
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  140
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ140


   subroutine PrimitiveContractionBCPUSegQ140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  140*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  140*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  140*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ140


   subroutine PrimitiveContractionACPUSegQ100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  100,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  100
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  100
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ100


   subroutine PrimitiveContractionBCPUSegQ100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  100*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  100*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ100


   subroutine PrimitiveContractionACPUSegQ200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  200,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  200
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  200
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ200


   subroutine PrimitiveContractionBCPUSegQ200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  200*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  200*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  200*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ200


   subroutine PrimitiveContractionACPUSegQ350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  350,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  350
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  350
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ350


   subroutine PrimitiveContractionBCPUSegQ350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  350*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  350*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  350*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ350


   subroutine PrimitiveContractionACPUSegQ400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  400,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  400
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  400
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ400


   subroutine PrimitiveContractionBCPUSegQ400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  400*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  400*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ400


   subroutine PrimitiveContractionACPUSegQ700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2(  700,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1,  700
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1,  700
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ700


   subroutine PrimitiveContractionBCPUSegQ700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  700*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1,  700*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1,  700*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ700


   subroutine PrimitiveContractionACPUSegQ1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimA,nPrimB*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContA,nPrimB*nPasses)
    !
    integer :: iP,iContA,iPrimA,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) 
   do iP = 1,nPasses*nPrimB
    do iContA=1,nContA
     do iTUV=1, 1225
      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk
     enddo
     do iPrimA=1,nPrimA
      TMP = ACC(iPrimA,iContA)
      do iTUV=1, 1225
       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionACPUSegQ1225


   subroutine PrimitiveContractionBCPUSegQ1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225*nContA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContA,nContB,nPasses)
    !
    integer :: iPassP,iContB,iPrimB,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) 
   do iPassP = 1,nPasses
    do iContB=1,nContB
     do iTUV=1, 1225*nContA
      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk
     enddo
     do iPrimB=1,nPrimB
      TMP = BCC(iPrimB,iContB)
      do iTUV=1, 1225*nContA
       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)
      enddo
     enddo
    enddo
   enddo
!$OMP END DO
   end subroutine PrimitiveContractionBCPUSegQ1225

END MODULE IchorEriCoulombintegralCPUOBSGeneralModSegQ
