MODULE IchorEriCoulombintegralCPUOBSGeneralModSegP
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted LHS and a General Contracted RHS Basisset 
use IchorEriCoulombintegralCPUMcMGeneralMod
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000ModGen
use AGC_CPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODASegP
use AGC_CPU_OBS_VERTICALRECURRENCEMODBSegP
use AGC_CPU_OBS_VERTICALRECURRENCEMODDSegP
use AGC_CPU_OBS_VERTICALRECURRENCEMODCSegP
use AGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_CPU_OBS_TRMODAtoCSegP
use AGC_CPU_OBS_TRMODAtoDSegP
use AGC_CPU_OBS_TRMODBtoCSegP
use AGC_CPU_OBS_TRMODBtoDSegP
use AGC_CPU_OBS_TRMODCtoASegP
use AGC_CPU_OBS_TRMODDtoASegP
use AGC_CPU_OBS_TRMODCtoBSegP
use AGC_CPU_OBS_TRMODDtoBSegP
use AGC_CPU_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_CPU_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_CPU_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_CPU_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_CPU_OBS_Sphcontract1Mod
use AGC_CPU_OBS_Sphcontract2Mod
  
private   
public :: ICI_CPU_OBS_SegP,ICI_CPU_OBS_general_sizeSegP  
  
CONTAINS
  
  
  subroutine ICI_CPU_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*1.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*1.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP1(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
         call PrimitiveContractionDCPUSegP1(TMParray1,LOCALINTS,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,1,&
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
        IF(nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,1,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*81.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(1,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(35,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(1,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(35,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(1,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(25,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q3AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(25,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q4AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP1225(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP1225(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*1260.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*875.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(35,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*900.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*625.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(25,nContQ*nPasses,TMParray1,&
            & LOCALINTS)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,1,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(1,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(1,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(1,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP35(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP35(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP4(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP4(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,1,&
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
        IF(nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP10(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP10(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(1,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(35,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP16(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP16(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(9,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(9,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(9,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(9,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*315.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*324.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(9,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP20(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP20(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(1,nContQ*nPasses,TMParray2,&
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
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(35,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP40(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP40(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP100(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP100(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP80(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP80(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP200(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP200(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP400(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP400(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP140(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP140(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
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
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP350(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP350(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(25,nContQ*nPasses,TMParray1,&
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
        IF(nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCCPUSegP700(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDCPUSegP700(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(25,nContQ*nPasses,TMParray1,&
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
  end subroutine ICI_CPU_OBS_SegP
  
  
  subroutine ICI_CPU_OBS_general_sizeSegP(TMParray1maxsize,&
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
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,1*nContC*nPrimD)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,35*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,35*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQ)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQ)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,4*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,4*nContC*nContD)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,16*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,16*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQ)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQ)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQ)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQ)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQ)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,315*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,324*nContQ)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQ)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,10*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,10*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,40*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,40*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,100*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,100*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,20*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,20*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,80*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,80*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,200*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,200*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,400*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,400*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQ)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQ)
    TMParray2maxSize = MAX(TMParray2maxSize,35*nContC*nPrimD)
    TMParray1maxSize = MAX(TMParray1maxSize,35*nContC*nContD)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQ)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQ)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,140*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,140*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQ)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,350*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,350*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQ)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,700*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,700*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQ)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,9*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,165*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nPrimQ)
    TMParray1maxSize = MAX(TMParray1maxSize,1225*nContC*nPrimD)
    TMParray2maxSize = MAX(TMParray2maxSize,1225*nContC*nContD)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContQ)
    CASE DEFAULT
     call ICI_CPU_McM_general_size(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP,&
         & .TRUE.,.FALSE.)
    END SELECT
  end subroutine ICI_CPU_OBS_general_sizeSegP

   subroutine PrimitiveContractionCCPUSegP1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      tmp = 0.0E0_realk
      do iPrimC=1,nPrimC
       tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
      enddo
      AUXarrayCont(iContC,iP) = tmp
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP1

   subroutine PrimitiveContractionDCPUSegP1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iContC
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iContC,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
      do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iPassP) = AUXarrayCont(iContC,iContD,iPassP) + tmp*AUXarray2(iContC,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP1

   subroutine PrimitiveContractionCCPUSegP4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(    4,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,    4
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,    4
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP4

   subroutine PrimitiveContractionDCPUSegP4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(    4*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,    4*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,    4*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP4

   subroutine PrimitiveContractionCCPUSegP10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   10,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   10
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   10
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP10

   subroutine PrimitiveContractionDCPUSegP10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   10*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   10*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   10*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP10

   subroutine PrimitiveContractionCCPUSegP20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   20,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   20
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   20
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP20

   subroutine PrimitiveContractionDCPUSegP20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   20*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   20*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   20*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP20

   subroutine PrimitiveContractionCCPUSegP35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   35,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   35
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   35
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP35

   subroutine PrimitiveContractionDCPUSegP35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   35*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   35*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   35*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP35

   subroutine PrimitiveContractionCCPUSegP16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   16,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   16
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   16
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP16

   subroutine PrimitiveContractionDCPUSegP16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   16*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   16*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   16*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP16

   subroutine PrimitiveContractionCCPUSegP40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   40,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   40
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   40
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP40

   subroutine PrimitiveContractionDCPUSegP40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   40*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   40*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   40*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP40

   subroutine PrimitiveContractionCCPUSegP80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(   80,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,   80
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,   80
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP80

   subroutine PrimitiveContractionDCPUSegP80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   80*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,   80*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,   80*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP80

   subroutine PrimitiveContractionCCPUSegP140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  140,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  140
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  140
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP140

   subroutine PrimitiveContractionDCPUSegP140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  140*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  140*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  140*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP140

   subroutine PrimitiveContractionCCPUSegP100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  100,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  100
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  100
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP100

   subroutine PrimitiveContractionDCPUSegP100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  100*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  100*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  100*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP100

   subroutine PrimitiveContractionCCPUSegP200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  200,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  200
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  200
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP200

   subroutine PrimitiveContractionDCPUSegP200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  200*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  200*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  200*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP200

   subroutine PrimitiveContractionCCPUSegP350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  350,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  350
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  350
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP350

   subroutine PrimitiveContractionDCPUSegP350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  350*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  350*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  350*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP350

   subroutine PrimitiveContractionCCPUSegP400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  400,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  400
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  400
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP400

   subroutine PrimitiveContractionDCPUSegP400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  400*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  400*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  400*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP400

   subroutine PrimitiveContractionCCPUSegP700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(  700,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1,  700
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1,  700
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP700

   subroutine PrimitiveContractionDCPUSegP700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  700*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1,  700*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1,  700*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP700

   subroutine PrimitiveContractionCCPUSegP1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD*nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContC,nPrimD*nPasses)
    !
    integer :: iP,iContC,iPrimC,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) 
    do iP = 1,nPasses*nPrimD
     do iContC=1,nContC
      do iTUV=1, 1225
       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk
      enddo
      do iPrimC=1,nPrimC
       tmp = CCC(iPrimC,iContC)
       do iTUV=1, 1225
        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionCCPUSegP1225

   subroutine PrimitiveContractionDCPUSegP1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2( 1225*nContC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225*nContC,nContD,nPasses)
    !
    integer :: iPassP,iContD,iPrimD,iTUV
    real(realk) :: TMP
!$OMP DO &
!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) 
    do iPassP = 1,nPasses
     do iContD=1,nContD
       do iTUV=1, 1225*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk
       enddo
      do iPrimD=1,nPrimD
       tmp = DCC(iPrimD,iContD)
       do iTUV=1, 1225*nContC
        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)
       enddo
      enddo
     enddo
    enddo
!$OMP END DO
   end subroutine PrimitiveContractionDCPUSegP1225
END MODULE IchorEriCoulombintegralCPUOBSGeneralModSegP
