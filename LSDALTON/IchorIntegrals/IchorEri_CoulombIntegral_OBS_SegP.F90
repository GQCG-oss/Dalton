MODULE IchorEriCoulombintegralOBSGeneralModSegP
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted LHS and a General Contracted RHS Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_VERTICALRECURRENCEMODASegP
use AGC_CPU_OBS_VERTICALRECURRENCEMODBSegP
use AGC_CPU_OBS_VERTICALRECURRENCEMODDSegP
use AGC_CPU_OBS_VERTICALRECURRENCEMODCSegP
use AGC_CPU_OBS_VERTICALRECURRENCEMODA
use AGC_CPU_OBS_VERTICALRECURRENCEMODB
use AGC_CPU_OBS_VERTICALRECURRENCEMODC
use AGC_CPU_OBS_VERTICALRECURRENCEMODD
use AGC_OBS_TRANSFERRECURRENCEMODAtoCSegP
use AGC_OBS_TRANSFERRECURRENCEMODAtoDSegP
use AGC_OBS_TRANSFERRECURRENCEMODBtoCSegP
use AGC_OBS_TRANSFERRECURRENCEMODBtoDSegP
use AGC_OBS_TRANSFERRECURRENCEMODCtoASegP
use AGC_OBS_TRANSFERRECURRENCEMODDtoASegP
use AGC_OBS_TRANSFERRECURRENCEMODCtoBSegP
use AGC_OBS_TRANSFERRECURRENCEMODDtoBSegP
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_OBS_SegP,IchorCoulombIntegral_OBS_general_sizeSegP  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
       & BasisCont1,BasisCont2,BasisCont3)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsA,nAtomsB
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ)           !qcent(3,nPrimQ)
    real(realk),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP)
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
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3,nAtomsA),Bcenter(3,nAtomsB),Ccenter(3),Dcenter(3)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    integer,intent(in) :: BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize
    real(realk) :: BasisCont1(BasisCont1maxsize) 
    real(realk) :: BasisCont2(BasisCont2maxsize) 
    real(realk) :: BasisCont3(BasisCont3maxsize) 
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
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2(1:nPrimQ*nPasses*1))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*1.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP1(TMParray2(1:nPrimQ*nPasses*1),&
            & LOCALINTS(1:nContQ*nPasses*1),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*4))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP4(TMParray2(1:nPrimQ*nPasses*4),&
            & TMParray1(1:nContQ*nPasses*4),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,1,Pdistance12,TMParray1(1:nContQ*nPasses*4),&
            & LOCALINTS(1:nContQ*nPasses*3),lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*10),TMParray1(1:nPrimQ*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP16(TMParray1(1:nPrimQ*nPasses*16),&
            & TMParray2(1:nContQ*nPasses*16),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*16),&
            & TMParray1(1:nContQ*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*12),&
            & LOCALINTS(1:nContQ*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*30),&
            & LOCALINTS(1:nContQ*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP10(TMParray2(1:nPrimQ*nPasses*10),&
            & TMParray1(1:nContQ*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,1,Pdistance12,TMParray1(1:nContQ*nPasses*10),&
            & LOCALINTS(1:nContQ*nPasses*9),lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*36),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1(1:nContQ*nPasses*36),&
            & LOCALINTS(1:nContQ*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*81.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*81),lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP10(TMParray2(1:nPrimQ*nPasses*10),&
            & TMParray1(1:nContQ*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,1,Pdistance12,TMParray1(1:nContQ*nPasses*10),&
            & TMParray2(1:nContQ*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*6),&
            & LOCALINTS(1:nContQ*nPasses*5))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*24),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*24),&
            & TMParray2(1:nContQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*20),&
            & LOCALINTS(1:nContQ*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*50),&
            & LOCALINTS(1:nContQ*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*50),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*30),&
            & LOCALINTS(1:nContQ*nPasses*25))
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*120),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*120),&
            & TMParray2(1:nContQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP350(TMParray1(1:nPrimQ*nPasses*350),&
            & TMParray2(1:nContQ*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,35,Pdistance12,TMParray2(1:nContQ*nPasses*350),&
            & TMParray1(1:nContQ*nPasses*210),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(35,nContQ*nPasses,TMParray1(1:nContQ*nPasses*210),&
            & TMParray2(1:nContQ*nPasses*175))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*175),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*180),&
            & LOCALINTS(1:nContQ*nPasses*125))
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP20(TMParray2(1:nPrimQ*nPasses*20),&
            & TMParray1(1:nContQ*nPasses*20),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,1,Pdistance12,TMParray1(1:nContQ*nPasses*20),&
            & TMParray2(1:nContQ*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*18),&
            & LOCALINTS(1:nContQ*nPasses*15))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP80(TMParray1(1:nPrimQ*nPasses*80),&
            & TMParray2(1:nContQ*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*80),&
            & TMParray1(1:nContQ*nPasses*72),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*72),&
            & TMParray2(1:nContQ*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*60),&
            & LOCALINTS(1:nContQ*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*180),&
            & TMParray2(1:nContQ*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*150),&
            & LOCALINTS(1:nContQ*nPasses*135),lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*180),&
            & TMParray2(1:nContQ*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*150),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP400(TMParray1(1:nPrimQ*nPasses*400),&
            & TMParray2(1:nContQ*nPasses*400),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*400),&
            & TMParray1(1:nContQ*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*360),&
            & TMParray2(1:nContQ*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*300),&
            & TMParray1(1:nContQ*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*270),&
            & LOCALINTS(1:nContQ*nPasses*225))
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*120))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*120),TMParray1(1:nPrimQ*nPasses*700))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP700(TMParray1(1:nPrimQ*nPasses*700),&
            & TMParray2(1:nContQ*nPasses*700),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,35,Pdistance12,TMParray2(1:nContQ*nPasses*700),&
            & TMParray1(1:nContQ*nPasses*630),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(35,nContQ*nPasses,TMParray1(1:nContQ*nPasses*630),&
            & TMParray2(1:nContQ*nPasses*525))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*525),&
            & TMParray1(1:nContQ*nPasses*540),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*540),&
            & LOCALINTS(1:nContQ*nPasses*375))
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP35(TMParray2(1:nPrimQ*nPasses*35),&
            & TMParray1(1:nContQ*nPasses*35),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,1,Pdistance12,TMParray1(1:nContQ*nPasses*35),&
            & TMParray2(1:nContQ*nPasses*36),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*36),&
            & LOCALINTS(1:nContQ*nPasses*25))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q1AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*140))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP140(TMParray1(1:nPrimQ*nPasses*140),&
            & TMParray2(1:nContQ*nPasses*140),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*140),&
            & TMParray1(1:nContQ*nPasses*144),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*144),&
            & TMParray2(1:nContQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,25,Qdistance12,TMParray2(1:nContQ*nPasses*100),&
            & LOCALINTS(1:nContQ*nPasses*75),lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP350(TMParray1(1:nPrimQ*nPasses*350),&
            & TMParray2(1:nContQ*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*350),&
            & TMParray1(1:nContQ*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*360),&
            & TMParray2(1:nContQ*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,25,Qdistance12,TMParray2(1:nContQ*nPasses*250),&
            & LOCALINTS(1:nContQ*nPasses*225),lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP350(TMParray1(1:nPrimQ*nPasses*350),&
            & TMParray2(1:nContQ*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*350),&
            & TMParray1(1:nContQ*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*360),&
            & TMParray2(1:nContQ*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,25,Qdistance12,TMParray2(1:nContQ*nPasses*250),&
            & TMParray1(1:nContQ*nPasses*150),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(25,nContQ*nPasses,TMParray1(1:nContQ*nPasses*150),&
            & LOCALINTS(1:nContQ*nPasses*125))
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*120))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q3AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*120),TMParray1(1:nPrimQ*nPasses*700))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP700(TMParray1(1:nPrimQ*nPasses*700),&
            & TMParray2(1:nContQ*nPasses*700),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*700),&
            & TMParray1(1:nContQ*nPasses*720),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*720),&
            & TMParray2(1:nContQ*nPasses*500))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,25,Qdistance12,TMParray2(1:nContQ*nPasses*500),&
            & TMParray1(1:nContQ*nPasses*450),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(25,nContQ*nPasses,TMParray1(1:nContQ*nPasses*450),&
            & LOCALINTS(1:nContQ*nPasses*375))
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*165.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*165))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q4AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*165),TMParray1(1:nPrimQ*nPasses*1225))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP1225(TMParray1(1:nPrimQ*nPasses*1225),&
            & TMParray2(1:nContQ*nPasses*1225),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*1260.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,35,Pdistance12,TMParray2(1:nContQ*nPasses*1225),&
            & TMParray1(1:nContQ*nPasses*1260),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*875.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(35,nContQ*nPasses,TMParray1(1:nContQ*nPasses*1260),&
            & TMParray2(1:nContQ*nPasses*875))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*900.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,25,Qdistance12,TMParray2(1:nContQ*nPasses*875),&
            & TMParray1(1:nContQ*nPasses*900),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*625.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nContQ*nPasses,TMParray1(1:nContQ*nPasses*900),&
            & LOCALINTS(1:nContQ*nPasses*625))
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*4))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP4(TMParray2(1:nPrimQ*nPasses*4),&
            & TMParray1(1:nContQ*nPasses*4),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,1,Qdistance12,TMParray1(1:nContQ*nPasses*4),&
            & LOCALINTS(1:nContQ*nPasses*3),lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP10(TMParray2(1:nPrimQ*nPasses*10),&
            & TMParray1(1:nContQ*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1(1:nContQ*nPasses*10),&
            & TMParray2(1:nContQ*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*6),&
            & LOCALINTS(1:nContQ*nPasses*5))
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*4))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP4(TMParray2(1:nPrimQ*nPasses*4),&
            & TMParray1(1:nContQ*nPasses*4),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray1(1:nContQ*nPasses*4),&
            & LOCALINTS(1:nContQ*nPasses*3),lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP10(TMParray2(1:nPrimQ*nPasses*10),&
            & TMParray1(1:nContQ*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1(1:nContQ*nPasses*10),&
            & LOCALINTS(1:nContQ*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP20(TMParray2(1:nPrimQ*nPasses*20),&
            & TMParray1(1:nContQ*nPasses*20),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,1,Qdistance12,TMParray1(1:nContQ*nPasses*20),&
            & TMParray2(1:nContQ*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*18),&
            & LOCALINTS(1:nContQ*nPasses*15))
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP10(TMParray2(1:nPrimQ*nPasses*10),&
            & TMParray1(1:nContQ*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,1,Qdistance12,TMParray1(1:nContQ*nPasses*10),&
            & TMParray2(1:nContQ*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*6),&
            & LOCALINTS(1:nContQ*nPasses*5))
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP20(TMParray2(1:nPrimQ*nPasses*20),&
            & TMParray1(1:nContQ*nPasses*20),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,1,Qdistance12,TMParray1(1:nContQ*nPasses*20),&
            & TMParray2(1:nContQ*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*18),&
            & LOCALINTS(1:nContQ*nPasses*15))
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP35(TMParray2(1:nPrimQ*nPasses*35),&
            & TMParray1(1:nContQ*nPasses*35),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,1,Qdistance12,TMParray1(1:nContQ*nPasses*35),&
            & TMParray2(1:nContQ*nPasses*36),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*36),&
            & LOCALINTS(1:nContQ*nPasses*25))
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*4))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP4(TMParray2(1:nPrimQ*nPasses*4),&
            & TMParray1(1:nContQ*nPasses*4),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,1,Pdistance12,TMParray1(1:nContQ*nPasses*4),&
            & LOCALINTS(1:nContQ*nPasses*3),lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*10),TMParray1(1:nPrimQ*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP16(TMParray1(1:nPrimQ*nPasses*16),&
            & TMParray2(1:nContQ*nPasses*16),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*16),&
            & TMParray1(1:nContQ*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*12),&
            & LOCALINTS(1:nContQ*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*30),&
            & TMParray2(1:nContQ*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*18),&
            & LOCALINTS(1:nContQ*nPasses*15))
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*10),TMParray1(1:nPrimQ*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP16(TMParray1(1:nPrimQ*nPasses*16),&
            & TMParray2(1:nContQ*nPasses*16),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*16),&
            & TMParray1(1:nContQ*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*12),&
            & LOCALINTS(1:nContQ*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*30),&
            & LOCALINTS(1:nContQ*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP80(TMParray1(1:nPrimQ*nPasses*80),&
            & TMParray2(1:nContQ*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*80),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*54),&
            & LOCALINTS(1:nContQ*nPasses*45))
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*30),&
            & TMParray2(1:nContQ*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*18),&
            & LOCALINTS(1:nContQ*nPasses*15))
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP80(TMParray1(1:nPrimQ*nPasses*80),&
            & TMParray2(1:nContQ*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*80),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*54),&
            & LOCALINTS(1:nContQ*nPasses*45))
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*140))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP140(TMParray1(1:nPrimQ*nPasses*140),&
            & TMParray2(1:nContQ*nPasses*140),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContQ*nPasses,35,Pdistance12,TMParray2(1:nContQ*nPasses*140),&
            & TMParray1(1:nContQ*nPasses*105),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*105),&
            & TMParray2(1:nContQ*nPasses*108),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*108),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP10(TMParray2(1:nPrimQ*nPasses*10),&
            & TMParray1(1:nContQ*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,1,Pdistance12,TMParray1(1:nContQ*nPasses*10),&
            & TMParray2(1:nContQ*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*6),&
            & LOCALINTS(1:nContQ*nPasses*5))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*24),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*24),&
            & TMParray2(1:nContQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*20),&
            & LOCALINTS(1:nContQ*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*50),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*30),&
            & LOCALINTS(1:nContQ*nPasses*25))
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*24),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*24),&
            & TMParray2(1:nContQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*20),&
            & LOCALINTS(1:nContQ*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*50),&
            & LOCALINTS(1:nContQ*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3DtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*120),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*120),&
            & TMParray2(1:nContQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*50),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*30),&
            & LOCALINTS(1:nContQ*nPasses*25))
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*120),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*120),&
            & TMParray2(1:nContQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP350(TMParray1(1:nPrimQ*nPasses*350),&
            & TMParray2(1:nContQ*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContQ*nPasses,35,Pdistance12,TMParray2(1:nContQ*nPasses*350),&
            & TMParray1(1:nContQ*nPasses*210),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(35,nContQ*nPasses,TMParray1(1:nContQ*nPasses*210),&
            & TMParray2(1:nContQ*nPasses*175))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*175),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*180),&
            & LOCALINTS(1:nContQ*nPasses*125))
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*10),TMParray1(1:nPrimQ*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP16(TMParray1(1:nPrimQ*nPasses*16),&
            & TMParray2(1:nContQ*nPasses*16),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*16),&
            & TMParray1(1:nContQ*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*12),&
            & LOCALINTS(1:nContQ*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*30),&
            & TMParray2(1:nContQ*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*18),&
            & LOCALINTS(1:nContQ*nPasses*15))
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP80(TMParray1(1:nPrimQ*nPasses*80),&
            & TMParray2(1:nContQ*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*80),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*54),&
            & LOCALINTS(1:nContQ*nPasses*45))
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*30),&
            & TMParray2(1:nContQ*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*18),&
            & LOCALINTS(1:nContQ*nPasses*15))
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP80(TMParray1(1:nPrimQ*nPasses*80),&
            & TMParray2(1:nContQ*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*80),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*54),&
            & LOCALINTS(1:nContQ*nPasses*45))
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*140))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP140(TMParray1(1:nPrimQ*nPasses*140),&
            & TMParray2(1:nContQ*nPasses*140),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContQ*nPasses,35,Pdistance12,TMParray2(1:nContQ*nPasses*140),&
            & TMParray1(1:nContQ*nPasses*105),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,3,Qdistance12,TMParray1(1:nContQ*nPasses*105),&
            & TMParray2(1:nContQ*nPasses*108),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nContQ*nPasses,TMParray2(1:nContQ*nPasses*108),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*36),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,9,Qdistance12,TMParray1(1:nContQ*nPasses*36),&
            & LOCALINTS(1:nContQ*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1(1:nContQ*nPasses*90),&
            & TMParray2(1:nContQ*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(9,nContQ*nPasses,TMParray2(1:nContQ*nPasses*54),&
            & LOCALINTS(1:nContQ*nPasses*45))
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,9,Qdistance12,TMParray1(1:nContQ*nPasses*180),&
            & TMParray2(1:nContQ*nPasses*162),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(9,nContQ*nPasses,TMParray2(1:nContQ*nPasses*162),&
            & LOCALINTS(1:nContQ*nPasses*135))
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,9,Qdistance12,TMParray1(1:nContQ*nPasses*90),&
            & TMParray2(1:nContQ*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(9,nContQ*nPasses,TMParray2(1:nContQ*nPasses*54),&
            & LOCALINTS(1:nContQ*nPasses*45))
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,9,Qdistance12,TMParray1(1:nContQ*nPasses*180),&
            & TMParray2(1:nContQ*nPasses*162),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(9,nContQ*nPasses,TMParray2(1:nContQ*nPasses*162),&
            & LOCALINTS(1:nContQ*nPasses*135))
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q4CtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP350(TMParray1(1:nPrimQ*nPasses*350),&
            & TMParray2(1:nContQ*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*315.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContQ*nPasses,35,Pdistance12,TMParray2(1:nContQ*nPasses*350),&
            & TMParray1(1:nContQ*nPasses*315),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*324.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,9,Qdistance12,TMParray1(1:nContQ*nPasses*315),&
            & TMParray2(1:nContQ*nPasses*324),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(9,nContQ*nPasses,TMParray2(1:nContQ*nPasses*324),&
            & LOCALINTS(1:nContQ*nPasses*225))
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegP3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP20(TMParray2(1:nPrimQ*nPasses*20),&
            & TMParray1(1:nContQ*nPasses*20),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,1,Pdistance12,TMParray1(1:nContQ*nPasses*20),&
            & TMParray2(1:nContQ*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(1,nContQ*nPasses,TMParray2(1:nContQ*nPasses*18),&
            & LOCALINTS(1:nContQ*nPasses*15))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q1BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP80(TMParray1(1:nPrimQ*nPasses*80),&
            & TMParray2(1:nContQ*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*80),&
            & TMParray1(1:nContQ*nPasses*72),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*72),&
            & TMParray2(1:nContQ*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*60),&
            & LOCALINTS(1:nContQ*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*180),&
            & TMParray2(1:nContQ*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*150),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q1BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP80(TMParray1(1:nPrimQ*nPasses*80),&
            & TMParray2(1:nContQ*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*80),&
            & TMParray1(1:nContQ*nPasses*72),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*72),&
            & TMParray2(1:nContQ*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*60),&
            & LOCALINTS(1:nContQ*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*180),&
            & TMParray2(1:nContQ*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*150),&
            & LOCALINTS(1:nContQ*nPasses*135),lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3BtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP400(TMParray1(1:nPrimQ*nPasses*400),&
            & TMParray2(1:nContQ*nPasses*400),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*400),&
            & TMParray1(1:nContQ*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*360),&
            & TMParray2(1:nContQ*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*300),&
            & TMParray1(1:nContQ*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*270),&
            & LOCALINTS(1:nContQ*nPasses*225))
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*180),&
            & TMParray2(1:nContQ*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*150),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3BtoCSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP400(TMParray1(1:nPrimQ*nPasses*400),&
            & TMParray2(1:nContQ*nPasses*400),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*400),&
            & TMParray1(1:nContQ*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*360),&
            & TMParray2(1:nContQ*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*300),&
            & TMParray1(1:nContQ*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*270),&
            & LOCALINTS(1:nContQ*nPasses*225))
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*120))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q4CtoBSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*120),TMParray1(1:nPrimQ*nPasses*700))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP700(TMParray1(1:nPrimQ*nPasses*700),&
            & TMParray2(1:nContQ*nPasses*700),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContQ*nPasses,35,Pdistance12,TMParray2(1:nContQ*nPasses*700),&
            & TMParray1(1:nContQ*nPasses*630),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(35,nContQ*nPasses,TMParray1(1:nContQ*nPasses*630),&
            & TMParray2(1:nContQ*nPasses*525))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*525),&
            & TMParray1(1:nContQ*nPasses*540),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*540),&
            & LOCALINTS(1:nContQ*nPasses*375))
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*20),TMParray1(1:nPrimQ*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP40(TMParray1(1:nPrimQ*nPasses*40),&
            & TMParray2(1:nContQ*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*40),&
            & TMParray1(1:nContQ*nPasses*24),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*24),&
            & TMParray2(1:nContQ*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*20),&
            & LOCALINTS(1:nContQ*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP100(TMParray1(1:nPrimQ*nPasses*100),&
            & TMParray2(1:nContQ*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*60),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*60),&
            & TMParray2(1:nContQ*nPasses*50))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*50),&
            & TMParray1(1:nContQ*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*30),&
            & LOCALINTS(1:nContQ*nPasses*25))
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*120),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*120),&
            & TMParray2(1:nContQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,5,Qdistance12,TMParray2(1:nContQ*nPasses*100),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nContQ*nPasses,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*35),TMParray1(1:nPrimQ*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP80(TMParray1(1:nPrimQ*nPasses*80),&
            & TMParray2(1:nContQ*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*80),&
            & TMParray1(1:nContQ*nPasses*72),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*72),&
            & TMParray2(1:nContQ*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*60),&
            & LOCALINTS(1:nContQ*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP200(TMParray1(1:nPrimQ*nPasses*200),&
            & TMParray2(1:nContQ*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*200),&
            & TMParray1(1:nContQ*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*180),&
            & TMParray2(1:nContQ*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*150),&
            & TMParray1(1:nContQ*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*90),&
            & LOCALINTS(1:nContQ*nPasses*75))
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP400(TMParray1(1:nPrimQ*nPasses*400),&
            & TMParray2(1:nContQ*nPasses*400),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*400),&
            & TMParray1(1:nContQ*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*360),&
            & TMParray2(1:nContQ*nPasses*300))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,15,Qdistance12,TMParray2(1:nContQ*nPasses*300),&
            & TMParray1(1:nContQ*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContQ*nPasses,TMParray1(1:nContQ*nPasses*270),&
            & LOCALINTS(1:nContQ*nPasses*225))
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*56),TMParray1(1:nPrimQ*nPasses*140))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP140(TMParray1(1:nPrimQ*nPasses*140),&
            & TMParray2(1:nContQ*nPasses*140),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,4,Pdistance12,TMParray2(1:nContQ*nPasses*140),&
            & TMParray1(1:nContQ*nPasses*144),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nContQ*nPasses,TMParray1(1:nContQ*nPasses*144),&
            & TMParray2(1:nContQ*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQ,nPasses,25,Qdistance12,TMParray2(1:nContQ*nPasses*100),&
            & LOCALINTS(1:nContQ*nPasses*75),lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*84),TMParray1(1:nPrimQ*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP350(TMParray1(1:nPrimQ*nPasses*350),&
            & TMParray2(1:nContQ*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,10,Pdistance12,TMParray2(1:nContQ*nPasses*350),&
            & TMParray1(1:nContQ*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContQ*nPasses,TMParray1(1:nContQ*nPasses*360),&
            & TMParray2(1:nContQ*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2(1:nContQ*nPasses*250),&
            & TMParray1(1:nContQ*nPasses*150),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(25,nContQ*nPasses,TMParray1(1:nContQ*nPasses*150),&
            & LOCALINTS(1:nContQ*nPasses*125))
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimQ*nPasses*120))
#ifdef VAR_DEBUGICHOR
        IF(nPrimQ*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimQ*nPasses*120),TMParray1(1:nPrimQ*nPasses*700))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegP700(TMParray1(1:nPrimQ*nPasses*700),&
            & TMParray2(1:nContQ*nPasses*700),nPrimP,nPrimQ,nPasses,&
              & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContQ*nPasses,20,Pdistance12,TMParray2(1:nContQ*nPasses*700),&
            & TMParray1(1:nContQ*nPasses*720),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nContQ*nPasses,TMParray1(1:nContQ*nPasses*720),&
            & TMParray2(1:nContQ*nPasses*500))
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQ,nPasses,25,Qdistance12,TMParray2(1:nContQ*nPasses*500),&
            & TMParray1(1:nContQ*nPasses*450),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContQ*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContQ*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(25,nContQ*nPasses,TMParray1(1:nContQ*nPasses*450),&
            & LOCALINTS(1:nContQ*nPasses*375))
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_SegP',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_SegP
  
  
  subroutine IchorCoulombIntegral_OBS_general_sizeSegP(TMParray1maxsize,&
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
    BasisCont3maxsize =     1*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimQ)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
    BasisCont3maxsize =     4*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQ)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
    BasisCont3maxsize =    10*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
    BasisCont3maxsize =     4*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQ)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
    BasisCont3maxsize =    10*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
    BasisCont3maxsize =    20*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
    BasisCont3maxsize =    10*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
    BasisCont3maxsize =    20*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
    BasisCont3maxsize =    35*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQ)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
    BasisCont3maxsize =     4*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQ)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
    BasisCont3maxsize =    16*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
    BasisCont3maxsize =    16*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
    BasisCont3maxsize =    80*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
    BasisCont3maxsize =    80*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
    BasisCont3maxsize =   140*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQ)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
    BasisCont3maxsize =    10*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
    BasisCont3maxsize =   350*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
    BasisCont3maxsize =     4*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContQ)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
    BasisCont3maxsize =    16*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
    BasisCont3maxsize =    16*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQ)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
    BasisCont3maxsize =    80*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
    BasisCont3maxsize =    80*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
    BasisCont3maxsize =   140*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQ)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
    BasisCont3maxsize =    10*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQ)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQ)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQ)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQ)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQ)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
    BasisCont3maxsize =   350*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,315*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,324*nContQ)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
    BasisCont3maxsize =    20*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
    BasisCont3maxsize =    80*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
    BasisCont3maxsize =    80*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
    BasisCont3maxsize =   400*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
    BasisCont3maxsize =   400*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
    BasisCont3maxsize =   700*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQ)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
    BasisCont3maxsize =    10*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQ)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
    BasisCont3maxsize =    40*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQ)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
    BasisCont3maxsize =   100*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQ)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
    BasisCont3maxsize =   350*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
    BasisCont3maxsize =    20*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQ)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
    BasisCont3maxsize =    80*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
    BasisCont3maxsize =    80*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQ)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
    BasisCont3maxsize =   400*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
    BasisCont3maxsize =   200*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQ)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
    BasisCont3maxsize =   400*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQ)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
    BasisCont3maxsize =   700*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQ)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
    BasisCont3maxsize =    35*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQ)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQ)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
    BasisCont3maxsize =   140*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
    BasisCont3maxsize =   350*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQ)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
    BasisCont3maxsize =   140*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQ)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
    BasisCont3maxsize =   350*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
    BasisCont3maxsize =   700*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQ)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
    BasisCont3maxsize =   350*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQ)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
    BasisCont3maxsize =   700*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQ)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
    BasisCont3maxsize =  1225*nPrimD
       TMParray2maxSize = MAX(TMParray2maxSize,165*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nPrimQ)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContQ)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContQ)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContQ)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_general_sizeSegP

  subroutine PrimitiveContractionSegP1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    !Due to P being segmented the P contraction have already been done and we need to 
    !go from nPrimQ to nContQ
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: tmp,BasisCont3(nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       tmp = 0.0E0_realk
       do iPrimC=1,nPrimC
        tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iPrimD,iPassQ)
       enddo
       BasisCont3(iPrimD) = tmp
      enddo
      do iContD=1,nContD
       tmp = 0.0E0_realk
       do iPrimD=1,nPrimD
        tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iPrimD)
       enddo
       AUXarrayCont(iContC,iContD,iPassQ) = tmp
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP1


  subroutine PrimitiveContractionSegP4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(    4,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(    4,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,    4
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,    4
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP4

  subroutine PrimitiveContractionSegP10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   10,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(   10,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   10
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   10
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP10

  subroutine PrimitiveContractionSegP20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   20,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(   20,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   20
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   20
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP20

  subroutine PrimitiveContractionSegP35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   35,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(   35,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   35
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   35
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP35

  subroutine PrimitiveContractionSegP16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   16,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(   16,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   16
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   16
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP16

  subroutine PrimitiveContractionSegP40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   40,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(   40,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   40
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   40
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP40

  subroutine PrimitiveContractionSegP80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(   80,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(   80,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,   80
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,   80
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP80

  subroutine PrimitiveContractionSegP140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  140,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(  140,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  140
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  140
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP140

  subroutine PrimitiveContractionSegP100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  100,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(  100,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  100
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  100
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP100

  subroutine PrimitiveContractionSegP200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  200,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(  200,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  200
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  200
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP200

  subroutine PrimitiveContractionSegP350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  350,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(  350,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  350
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  350
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP350

  subroutine PrimitiveContractionSegP400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  400,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(  400,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  400
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  400
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP400

  subroutine PrimitiveContractionSegP700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(  700,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3(  700,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1,  700
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1,  700
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP700

  subroutine PrimitiveContractionSegP1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContQ,CCC,DCC,nPrimC,nContC,nPrimD,nContD,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ
    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimC,nPrimD,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContC,nContD,nPasses)
    !
    integer :: iPassQ,iContC,iContD,iPrimC,iPrimD,iTUV
    real(realk) :: TMP
    real(realk) :: BasisCont3( 1225,nPrimD)
    do iPassQ = 1,nPasses
     do iContC=1,nContC
      do iPrimD=1,nPrimD
       do iTUV=1, 1225
        tmp = 0.0E0_realk
        do iPrimC=1,nPrimC
         tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iTUV,iPrimC,iPrimD,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimD) = tmp
       enddo
      enddo
      do iContD=1,nContD
       do iTUV=1, 1225
        tmp = 0.0E0_realk
        do iPrimD=1,nPrimD
         tmp = tmp + DCC(iPrimD,iContD)*BasisCont3(iTUV,iPrimD)
        enddo
        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegP1225
END MODULE IchorEriCoulombintegralOBSGeneralModSegP
