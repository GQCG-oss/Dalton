MODULE IchorEriCoulombintegralGPUOBSGeneralModSegP
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted LHS and a General Contracted RHS Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_GPU_OBS_BUILDRJ000ModGen
use AGC_GPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_GPU_OBS_VERTICALRECURRENCEMODASegP
use AGC_GPU_OBS_VERTICALRECURRENCEMODBSegP
use AGC_GPU_OBS_VERTICALRECURRENCEMODDSegP
use AGC_GPU_OBS_VERTICALRECURRENCEMODCSegP
use AGC_GPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_GPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_GPU_OBS_TRMODAtoCSegP
use AGC_GPU_OBS_TRMODAtoDSegP
use AGC_GPU_OBS_TRMODBtoCSegP
use AGC_GPU_OBS_TRMODBtoDSegP
use AGC_GPU_OBS_TRMODCtoASegP
use AGC_GPU_OBS_TRMODDtoASegP
use AGC_GPU_OBS_TRMODCtoBSegP
use AGC_GPU_OBS_TRMODDtoBSegP
use AGC_GPU_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_GPU_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_GPU_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_GPU_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_GPU_OBS_Sphcontract1Mod
use AGC_GPU_OBS_Sphcontract2Mod
  
private   
public :: ICI_GPU_OBS_SegP,ICI_GPU_OBS_general_sizeSegP  
  
CONTAINS
  
  
  subroutine ICI_GPU_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
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
       & IatomAPass,iatomBPass,iASync)
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
    integer(kind=acckind),intent(in) :: iASync
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
        call VerticalRecurrenceGPUSegP0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*1.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,1,iASync)
         call PrimitiveContractionDGPUSegP(TMParray1,LOCALINTS,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,1,iASync)
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
        call VerticalRecurrenceGPUSegP1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & LOCALINTS,lupri,iASync)
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
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & LOCALINTS,lupri,iASync)
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
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*81.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
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
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(35,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
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
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q3AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(35,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
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
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(25,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q3AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(25,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*9.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*165.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q4AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*1260.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*875.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(35,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*900.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*625.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(25,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,4,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,4,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,20,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,10,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,20,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,1,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,1,35,iASync)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & LOCALINTS,lupri,iASync)
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
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
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
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA0(35,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP1Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,4,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,4,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*315.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*324.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(9,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUSegP3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(1,nContQ*nPasses,TMParray2,&
            & LOCALINTS,iASync)
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
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q3BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q3BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,35,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQ,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA1(35,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,10,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP3Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,20,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,4,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q1C0D1DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri,iASync)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,10,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q2C0D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(25,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000GPUGen7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceGPUGen7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceGPUP4Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nPrimD*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nContC*nPrimD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionCGPUSegP(TMParray2,TMParray1,nPrimP,nPrimQ,nPasses,&
               & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,35,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContC*nContD*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContC*nContD*nPasses too small',-1)
        ENDIF
#endif
         call PrimitiveContractionDGPUSegP(TMParray1,TMParray2,nPrimP,nPrimQ,nPasses,&
               & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,35,20,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQ,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_GPU_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1,&
            & TMParray2,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_GPU_RHS_Q3C1D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri,iASync)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(25,nContQ*nPasses,TMParray1,&
            & LOCALINTS,iASync)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in ICI_GPU_OBS_SegP',-1)
    END SELECT
  end subroutine ICI_GPU_OBS_SegP
  
  
  subroutine ICI_GPU_OBS_general_sizeSegP(TMParray1maxsize,&
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
        CALL ICHORQUIT('Unknown Case in ICI_OBS_general_size',-1)
    END SELECT
  end subroutine ICI_GPU_OBS_general_sizeSegP

   subroutine PrimitiveContractionCGPUSegP(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ,nTUVP,nTUVQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iContC,iPrimC
    real(realk) :: TMP
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContC,iPrimC,TMP) &
!$ACC PRESENT(nContC,nPrimC,nPasses,nPrimD,&
!$ACC        CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)
    do iP = 1,nPasses*nTUVP*nTUVQ*nPrimD
     do iContC=1,nContC
      tmp = 0.0E0_realk
      do iPrimC=1,nPrimC
       tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)
      enddo
      AUXarrayCont(iContC,iP) = tmp
     enddo
    enddo
   end subroutine PrimitiveContractionCGPUSegP

   subroutine PrimitiveContractionDGPUSegP(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,nTUVP,nTUVQ,iASync)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ,nTUVP,nTUVQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPasses*nTUVP*nTUVQ)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses*nTUVP*nTUVQ)
    integer(kind=acckind),intent(in) :: iASync
    !
    integer :: iP,iContC,iContD,iPrimD
    real(realk) :: TMP
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iContD,iPrimD,iContC,TMP) &
!$ACC PRESENT(nContC,nPrimD,nPasses,nContD,&
!$ACC         DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)
    do iP = 1,nPasses*nTUVP*nTUVQ
     do iContD=1,nContD
      do iContC=1,nContC
       AUXarrayCont(iContC,iContD,iP) = 0.0E0_realk
      enddo
      do iPrimD=1,nPrimD
       tmp = DCC(nPrimD,nContD)
       do iContC=1,nContC
        AUXarrayCont(iContC,iContD,iP) = AUXarrayCont(iContC,iContD,iP) + TMP*AUXarray2(iContC,iPrimD,iP)
       enddo
      enddo
     enddo
    enddo
   end subroutine PrimitiveContractionDGPUSegP
END MODULE IchorEriCoulombintegralGPUOBSGeneralModSegP
