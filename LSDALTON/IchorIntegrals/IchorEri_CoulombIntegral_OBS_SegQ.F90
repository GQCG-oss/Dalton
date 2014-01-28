MODULE IchorEriCoulombintegralOBSGeneralModSegQ
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for a General Contracted LHS Segmented contracted RHS and Basisset 
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_VERTICALRECURRENCEMODASegQ
use AGC_CPU_OBS_VERTICALRECURRENCEMODBSegQ
use AGC_CPU_OBS_VERTICALRECURRENCEMODDSegQ
use AGC_CPU_OBS_VERTICALRECURRENCEMODCSegQ
use AGC_CPU_OBS_VERTICALRECURRENCEMODA
use AGC_CPU_OBS_VERTICALRECURRENCEMODB
use AGC_CPU_OBS_VERTICALRECURRENCEMODC
use AGC_CPU_OBS_VERTICALRECURRENCEMODD
use AGC_OBS_TRANSFERRECURRENCEMODAtoCSegQ
use AGC_OBS_TRANSFERRECURRENCEMODAtoDSegQ
use AGC_OBS_TRANSFERRECURRENCEMODBtoCSegQ
use AGC_OBS_TRANSFERRECURRENCEMODBtoDSegQ
use AGC_OBS_TRANSFERRECURRENCEMODCtoASegQ
use AGC_OBS_TRANSFERRECURRENCEMODDtoASegQ
use AGC_OBS_TRANSFERRECURRENCEMODCtoBSegQ
use AGC_OBS_TRANSFERRECURRENCEMODDtoBSegQ
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_OBS_SegQ,IchorCoulombIntegral_OBS_general_sizeSegQ  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&
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
        IF(nPrimP*nPasses*1.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2(1:nPrimP*nPasses*1))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*1.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ1(TMParray2(1:nPrimP*nPasses*1),&
            & LOCALINTS(1:nContP*nPasses*1),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*4))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ4(TMParray2(1:nPrimP*nPasses*4),&
            & TMParray1(1:nContP*nPasses*4),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,1,Pdistance12,TMParray1(1:nContP*nPasses*4),&
            & LOCALINTS(1:nContP*nPasses*3),lupri)
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
               & TMParray2(1:nPrimP*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*10),TMParray1(1:nPrimP*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ16(TMParray1(1:nPrimP*nPasses*16),&
            & TMParray2(1:nContP*nPasses*16),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*16),&
            & TMParray1(1:nContP*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*12),&
            & LOCALINTS(1:nContP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*30),&
            & LOCALINTS(1:nContP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ10(TMParray2(1:nPrimP*nPasses*10),&
            & TMParray1(1:nContP*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,1,Pdistance12,TMParray1(1:nContP*nPasses*10),&
            & LOCALINTS(1:nContP*nPasses*9),lupri)
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
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*36),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,9,Qdistance12,TMParray1(1:nContP*nPasses*36),&
            & LOCALINTS(1:nContP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*90),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*81.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,9,Qdistance12,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*81),lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ10(TMParray2(1:nPrimP*nPasses*10),&
            & TMParray1(1:nContP*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContP*nPasses,1,Pdistance12,TMParray1(1:nContP*nPasses*10),&
            & TMParray2(1:nContP*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(1,nContP*nPasses,TMParray2(1:nContP*nPasses*6),&
            & LOCALINTS(1:nContP*nPasses*5))
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
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*24),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nContP*nPasses,TMParray1(1:nContP*nPasses*24),&
            & TMParray2(1:nContP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*20),&
            & LOCALINTS(1:nContP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*50),&
            & LOCALINTS(1:nContP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*50),&
            & TMParray1(1:nContP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContP*nPasses,TMParray1(1:nContP*nPasses*30),&
            & LOCALINTS(1:nContP*nPasses*25))
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*120),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nContP*nPasses,TMParray1(1:nContP*nPasses*120),&
            & TMParray2(1:nContP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nContP*nPasses,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ350(TMParray1(1:nPrimP*nPasses*350),&
            & TMParray2(1:nContP*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContP*nPasses,35,Pdistance12,TMParray2(1:nContP*nPasses*350),&
            & TMParray1(1:nContP*nPasses*210),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(35,nContP*nPasses,TMParray1(1:nContP*nPasses*210),&
            & TMParray2(1:nContP*nPasses*175))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*175),&
            & TMParray1(1:nContP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nContP*nPasses,TMParray1(1:nContP*nPasses*180),&
            & LOCALINTS(1:nContP*nPasses*125))
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ20(TMParray2(1:nPrimP*nPasses*20),&
            & TMParray1(1:nContP*nPasses*20),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP*nPasses,1,Pdistance12,TMParray1(1:nContP*nPasses*20),&
            & TMParray2(1:nContP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(1,nContP*nPasses,TMParray2(1:nContP*nPasses*18),&
            & LOCALINTS(1:nContP*nPasses*15))
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
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ80(TMParray1(1:nPrimP*nPasses*80),&
            & TMParray2(1:nContP*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*80),&
            & TMParray1(1:nContP*nPasses*72),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nContP*nPasses,TMParray1(1:nContP*nPasses*72),&
            & TMParray2(1:nContP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*60),&
            & LOCALINTS(1:nContP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContP*nPasses,TMParray1(1:nContP*nPasses*180),&
            & TMParray2(1:nContP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*150),&
            & LOCALINTS(1:nContP*nPasses*135),lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContP*nPasses,TMParray1(1:nContP*nPasses*180),&
            & TMParray2(1:nContP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*150),&
            & TMParray1(1:nContP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nContP*nPasses,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ400(TMParray1(1:nPrimP*nPasses*400),&
            & TMParray2(1:nContP*nPasses*400),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*300),&
            & TMParray1(1:nContP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContP*nPasses,TMParray1(1:nContP*nPasses*270),&
            & LOCALINTS(1:nContP*nPasses*225))
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*120))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*120),TMParray1(1:nPrimP*nPasses*700))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ700(TMParray1(1:nPrimP*nPasses*700),&
            & TMParray2(1:nContP*nPasses*700),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP*nPasses,35,Pdistance12,TMParray2(1:nContP*nPasses*700),&
            & TMParray1(1:nContP*nPasses*630),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(35,nContP*nPasses,TMParray1(1:nContP*nPasses*630),&
            & TMParray2(1:nContP*nPasses*525))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*525),&
            & TMParray1(1:nContP*nPasses*540),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nContP*nPasses,TMParray1(1:nContP*nPasses*540),&
            & LOCALINTS(1:nContP*nPasses*375))
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ35(TMParray2(1:nPrimP*nPasses*35),&
            & TMParray1(1:nContP*nPasses*35),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,1,Pdistance12,TMParray1(1:nContP*nPasses*35),&
            & TMParray2(1:nContP*nPasses*36),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(1,nContP*nPasses,TMParray2(1:nContP*nPasses*36),&
            & LOCALINTS(1:nContP*nPasses*25))
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
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q1AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*140))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ140(TMParray1(1:nPrimP*nPasses*140),&
            & TMParray2(1:nContP*nPasses*140),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*140),&
            & TMParray1(1:nContP*nPasses*144),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nContP*nPasses,TMParray1(1:nContP*nPasses*144),&
            & TMParray2(1:nContP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,25,Qdistance12,TMParray2(1:nContP*nPasses*100),&
            & LOCALINTS(1:nContP*nPasses*75),lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ350(TMParray1(1:nPrimP*nPasses*350),&
            & TMParray2(1:nContP*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*350),&
            & TMParray1(1:nContP*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContP*nPasses,TMParray1(1:nContP*nPasses*360),&
            & TMParray2(1:nContP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,25,Qdistance12,TMParray2(1:nContP*nPasses*250),&
            & LOCALINTS(1:nContP*nPasses*225),lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ350(TMParray1(1:nPrimP*nPasses*350),&
            & TMParray2(1:nContP*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*350),&
            & TMParray1(1:nContP*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContP*nPasses,TMParray1(1:nContP*nPasses*360),&
            & TMParray2(1:nContP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,25,Qdistance12,TMParray2(1:nContP*nPasses*250),&
            & TMParray1(1:nContP*nPasses*150),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(25,nContP*nPasses,TMParray1(1:nContP*nPasses*150),&
            & LOCALINTS(1:nContP*nPasses*125))
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*120))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q3AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*120),TMParray1(1:nPrimP*nPasses*700))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ700(TMParray1(1:nPrimP*nPasses*700),&
            & TMParray2(1:nContP*nPasses*700),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*700),&
            & TMParray1(1:nContP*nPasses*720),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nContP*nPasses,TMParray1(1:nContP*nPasses*720),&
            & TMParray2(1:nContP*nPasses*500))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,25,Qdistance12,TMParray2(1:nContP*nPasses*500),&
            & TMParray1(1:nContP*nPasses*450),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(25,nContP*nPasses,TMParray1(1:nContP*nPasses*450),&
            & LOCALINTS(1:nContP*nPasses*375))
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*165.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*165))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*1225.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q4AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*165),TMParray1(1:nPrimP*nPasses*1225))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ1225(TMParray1(1:nPrimP*nPasses*1225),&
            & TMParray2(1:nContP*nPasses*1225),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,25,Qdistance12,TMParray2(1:nContP*nPasses*875),&
            & TMParray1(1:nContP*nPasses*900),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*625.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(25,nContP*nPasses,TMParray1(1:nContP*nPasses*900),&
            & LOCALINTS(1:nContP*nPasses*625))
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*4))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ4(TMParray2(1:nPrimP*nPasses*4),&
            & TMParray1(1:nContP*nPasses*4),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,1,Qdistance12,TMParray1(1:nContP*nPasses*4),&
            & LOCALINTS(1:nContP*nPasses*3),lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ10(TMParray2(1:nPrimP*nPasses*10),&
            & TMParray1(1:nContP*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,1,Qdistance12,TMParray1(1:nContP*nPasses*10),&
            & TMParray2(1:nContP*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(1,nContP*nPasses,TMParray2(1:nContP*nPasses*6),&
            & LOCALINTS(1:nContP*nPasses*5))
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*4))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ4(TMParray2(1:nPrimP*nPasses*4),&
            & TMParray1(1:nContP*nPasses*4),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,1,Qdistance12,TMParray1(1:nContP*nPasses*4),&
            & LOCALINTS(1:nContP*nPasses*3),lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ10(TMParray2(1:nPrimP*nPasses*10),&
            & TMParray1(1:nContP*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,1,Qdistance12,TMParray1(1:nContP*nPasses*10),&
            & LOCALINTS(1:nContP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ20(TMParray2(1:nPrimP*nPasses*20),&
            & TMParray1(1:nContP*nPasses*20),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,1,Qdistance12,TMParray1(1:nContP*nPasses*20),&
            & TMParray2(1:nContP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(1,nContP*nPasses,TMParray2(1:nContP*nPasses*18),&
            & LOCALINTS(1:nContP*nPasses*15))
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ10(TMParray2(1:nPrimP*nPasses*10),&
            & TMParray1(1:nContP*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,1,Qdistance12,TMParray1(1:nContP*nPasses*10),&
            & TMParray2(1:nContP*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(1,nContP*nPasses,TMParray2(1:nContP*nPasses*6),&
            & LOCALINTS(1:nContP*nPasses*5))
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ20(TMParray2(1:nPrimP*nPasses*20),&
            & TMParray1(1:nContP*nPasses*20),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,1,Qdistance12,TMParray1(1:nContP*nPasses*20),&
            & TMParray2(1:nContP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(1,nContP*nPasses,TMParray2(1:nContP*nPasses*18),&
            & LOCALINTS(1:nContP*nPasses*15))
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ35(TMParray2(1:nPrimP*nPasses*35),&
            & TMParray1(1:nContP*nPasses*35),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,1,Qdistance12,TMParray1(1:nContP*nPasses*35),&
            & TMParray2(1:nContP*nPasses*36),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(1,nContP*nPasses,TMParray2(1:nContP*nPasses*36),&
            & LOCALINTS(1:nContP*nPasses*25))
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*4))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*4.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ4(TMParray2(1:nPrimP*nPasses*4),&
            & TMParray1(1:nContP*nPasses*4),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,1,Pdistance12,TMParray1(1:nContP*nPasses*4),&
            & LOCALINTS(1:nContP*nPasses*3),lupri)
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
               & TMParray2(1:nPrimP*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*10),TMParray1(1:nPrimP*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ16(TMParray1(1:nPrimP*nPasses*16),&
            & TMParray2(1:nContP*nPasses*16),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*16),&
            & TMParray1(1:nContP*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*12),&
            & LOCALINTS(1:nContP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*30),&
            & TMParray2(1:nContP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nContP*nPasses,TMParray2(1:nContP*nPasses*18),&
            & LOCALINTS(1:nContP*nPasses*15))
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*10),TMParray1(1:nPrimP*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ16(TMParray1(1:nPrimP*nPasses*16),&
            & TMParray2(1:nContP*nPasses*16),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*16),&
            & TMParray1(1:nContP*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*12),&
            & LOCALINTS(1:nContP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*30),&
            & LOCALINTS(1:nContP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q3DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ80(TMParray1(1:nPrimP*nPasses*80),&
            & TMParray2(1:nContP*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*80),&
            & TMParray1(1:nContP*nPasses*60),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*60),&
            & TMParray2(1:nContP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nContP*nPasses,TMParray2(1:nContP*nPasses*54),&
            & LOCALINTS(1:nContP*nPasses*45))
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*30),&
            & TMParray2(1:nContP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nContP*nPasses,TMParray2(1:nContP*nPasses*18),&
            & LOCALINTS(1:nContP*nPasses*15))
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q3CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ80(TMParray1(1:nPrimP*nPasses*80),&
            & TMParray2(1:nContP*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*80),&
            & TMParray1(1:nContP*nPasses*60),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*60),&
            & TMParray2(1:nContP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nContP*nPasses,TMParray2(1:nContP*nPasses*54),&
            & LOCALINTS(1:nContP*nPasses*45))
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*140))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ140(TMParray1(1:nPrimP*nPasses*140),&
            & TMParray2(1:nContP*nPasses*140),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A0B1BtoA(nContP*nPasses,35,Pdistance12,TMParray2(1:nContP*nPasses*140),&
            & TMParray1(1:nContP*nPasses*105),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*105),&
            & TMParray2(1:nContP*nPasses*108),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nContP*nPasses,TMParray2(1:nContP*nPasses*108),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*10))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ10(TMParray2(1:nPrimP*nPasses*10),&
            & TMParray1(1:nContP*nPasses*10),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContP*nPasses,1,Pdistance12,TMParray1(1:nContP*nPasses*10),&
            & TMParray2(1:nContP*nPasses*6),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(1,nContP*nPasses,TMParray2(1:nContP*nPasses*6),&
            & LOCALINTS(1:nContP*nPasses*5))
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
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*24),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nContP*nPasses,TMParray1(1:nContP*nPasses*24),&
            & TMParray2(1:nContP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*20),&
            & LOCALINTS(1:nContP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*50),&
            & TMParray1(1:nContP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContP*nPasses,TMParray1(1:nContP*nPasses*30),&
            & LOCALINTS(1:nContP*nPasses*25))
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*24),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(4,nContP*nPasses,TMParray1(1:nContP*nPasses*24),&
            & TMParray2(1:nContP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*20),&
            & LOCALINTS(1:nContP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*50),&
            & LOCALINTS(1:nContP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3DtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*120),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nContP*nPasses,TMParray1(1:nContP*nPasses*120),&
            & TMParray2(1:nContP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nContP*nPasses,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*50),&
            & TMParray1(1:nContP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(5,nContP*nPasses,TMParray1(1:nContP*nPasses*30),&
            & LOCALINTS(1:nContP*nPasses*25))
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*120),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(20,nContP*nPasses,TMParray1(1:nContP*nPasses*120),&
            & TMParray2(1:nContP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(5,nContP*nPasses,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ350(TMParray1(1:nPrimP*nPasses*350),&
            & TMParray2(1:nContP*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A0B2BtoA(nContP*nPasses,35,Pdistance12,TMParray2(1:nContP*nPasses*350),&
            & TMParray1(1:nContP*nPasses*210),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA0(35,nContP*nPasses,TMParray1(1:nContP*nPasses*210),&
            & TMParray2(1:nContP*nPasses*175))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*175),&
            & TMParray1(1:nContP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(5,nContP*nPasses,TMParray1(1:nContP*nPasses*180),&
            & LOCALINTS(1:nContP*nPasses*125))
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*10))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*16.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*10),TMParray1(1:nPrimP*nPasses*16))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ16(TMParray1(1:nPrimP*nPasses*16),&
            & TMParray2(1:nContP*nPasses*16),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*16),&
            & TMParray1(1:nContP*nPasses*12),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*12),&
            & LOCALINTS(1:nContP*nPasses*9),lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*30),&
            & TMParray2(1:nContP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(3,nContP*nPasses,TMParray2(1:nContP*nPasses*18),&
            & LOCALINTS(1:nContP*nPasses*15))
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ80(TMParray1(1:nPrimP*nPasses*80),&
            & TMParray2(1:nContP*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*80),&
            & TMParray1(1:nContP*nPasses*60),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*60),&
            & TMParray2(1:nContP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(3,nContP*nPasses,TMParray2(1:nContP*nPasses*54),&
            & LOCALINTS(1:nContP*nPasses*45))
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q2CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*30),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*30),&
            & TMParray2(1:nContP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(3,nContP*nPasses,TMParray2(1:nContP*nPasses*18),&
            & LOCALINTS(1:nContP*nPasses*15))
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ80(TMParray1(1:nPrimP*nPasses*80),&
            & TMParray2(1:nContP*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*80),&
            & TMParray1(1:nContP*nPasses*60),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*60),&
            & TMParray2(1:nContP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(3,nContP*nPasses,TMParray2(1:nContP*nPasses*54),&
            & LOCALINTS(1:nContP*nPasses*45))
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP1Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*140))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ140(TMParray1(1:nPrimP*nPasses*140),&
            & TMParray2(1:nContP*nPasses*140),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P1A1B0AtoB(nContP*nPasses,35,Pdistance12,TMParray2(1:nContP*nPasses*140),&
            & TMParray1(1:nContP*nPasses*105),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,3,Qdistance12,TMParray1(1:nContP*nPasses*105),&
            & TMParray2(1:nContP*nPasses*108),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(3,nContP*nPasses,TMParray2(1:nContP*nPasses*108),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*36),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,9,Qdistance12,TMParray1(1:nContP*nPasses*36),&
            & LOCALINTS(1:nContP*nPasses*27),lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*90),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,9,Qdistance12,TMParray1(1:nContP*nPasses*90),&
            & TMParray2(1:nContP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(9,nContP*nPasses,TMParray2(1:nContP*nPasses*54),&
            & LOCALINTS(1:nContP*nPasses*45))
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*180),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,9,Qdistance12,TMParray1(1:nContP*nPasses*180),&
            & TMParray2(1:nContP*nPasses*162),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(9,nContP*nPasses,TMParray2(1:nContP*nPasses*162),&
            & LOCALINTS(1:nContP*nPasses*135))
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*90),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,9,Qdistance12,TMParray1(1:nContP*nPasses*90),&
            & TMParray2(1:nContP*nPasses*54),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(9,nContP*nPasses,TMParray2(1:nContP*nPasses*54),&
            & LOCALINTS(1:nContP*nPasses*45))
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*180),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,9,Qdistance12,TMParray1(1:nContP*nPasses*180),&
            & TMParray2(1:nContP*nPasses*162),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(9,nContP*nPasses,TMParray2(1:nContP*nPasses*162),&
            & LOCALINTS(1:nContP*nPasses*135))
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ350(TMParray1(1:nPrimP*nPasses*350),&
            & TMParray2(1:nContP*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*315.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A1B1AtoB(nContP*nPasses,35,Pdistance12,TMParray2(1:nContP*nPasses*350),&
            & TMParray1(1:nContP*nPasses*315),lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*324.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,9,Qdistance12,TMParray1(1:nContP*nPasses*315),&
            & TMParray2(1:nContP*nPasses*324),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(9,nContP*nPasses,TMParray2(1:nContP*nPasses*324),&
            & LOCALINTS(1:nContP*nPasses*225))
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSegQ3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
        !No reason for the Electron Transfer Recurrence Relation 
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ20(TMParray2(1:nPrimP*nPasses*20),&
            & TMParray1(1:nContP*nPasses*20),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP*nPasses,1,Pdistance12,TMParray1(1:nContP*nPasses*20),&
            & TMParray2(1:nContP*nPasses*18),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(1,nContP*nPasses,TMParray2(1:nContP*nPasses*18),&
            & LOCALINTS(1:nContP*nPasses*15))
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
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q1BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ80(TMParray1(1:nPrimP*nPasses*80),&
            & TMParray2(1:nContP*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*80),&
            & TMParray1(1:nContP*nPasses*72),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nContP*nPasses,TMParray1(1:nContP*nPasses*72),&
            & TMParray2(1:nContP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*60),&
            & LOCALINTS(1:nContP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContP*nPasses,TMParray1(1:nContP*nPasses*180),&
            & TMParray2(1:nContP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*150),&
            & TMParray1(1:nContP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nContP*nPasses,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q1BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ80(TMParray1(1:nPrimP*nPasses*80),&
            & TMParray2(1:nContP*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*80),&
            & TMParray1(1:nContP*nPasses*72),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(4,nContP*nPasses,TMParray1(1:nContP*nPasses*72),&
            & TMParray2(1:nContP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C1D0CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*60),&
            & LOCALINTS(1:nContP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContP*nPasses,TMParray1(1:nContP*nPasses*180),&
            & TMParray2(1:nContP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C1D1CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*150),&
            & LOCALINTS(1:nContP*nPasses*135),lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3BtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ400(TMParray1(1:nPrimP*nPasses*400),&
            & TMParray2(1:nContP*nPasses*400),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*300),&
            & TMParray1(1:nContP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContP*nPasses,TMParray1(1:nContP*nPasses*270),&
            & LOCALINTS(1:nContP*nPasses*225))
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(10,nContP*nPasses,TMParray1(1:nContP*nPasses*180),&
            & TMParray2(1:nContP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C2D0CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*150),&
            & TMParray1(1:nContP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC2(15,nContP*nPasses,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3BtoCSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ400(TMParray1(1:nPrimP*nPasses*400),&
            & TMParray2(1:nContP*nPasses*400),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q3C2D1CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*300),&
            & TMParray1(1:nContP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC2(15,nContP*nPasses,TMParray1(1:nContP*nPasses*270),&
            & LOCALINTS(1:nContP*nPasses*225))
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*120))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q4CtoBSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*120),TMParray1(1:nPrimP*nPasses*700))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ700(TMParray1(1:nPrimP*nPasses*700),&
            & TMParray2(1:nContP*nPasses*700),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A1B2BtoA(nContP*nPasses,35,Pdistance12,TMParray2(1:nContP*nPasses*700),&
            & TMParray1(1:nContP*nPasses*630),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA1(35,nContP*nPasses,TMParray1(1:nContP*nPasses*630),&
            & TMParray2(1:nContP*nPasses*525))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q4C2D2CtoD(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*525),&
            & TMParray1(1:nContP*nPasses*540),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ4_maxAngC2(15,nContP*nPasses,TMParray1(1:nContP*nPasses*540),&
            & LOCALINTS(1:nContP*nPasses*375))
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*40.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*20),TMParray1(1:nPrimP*nPasses*40))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ40(TMParray1(1:nPrimP*nPasses*40),&
            & TMParray2(1:nContP*nPasses*40),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*40),&
            & TMParray1(1:nContP*nPasses*24),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(4,nContP*nPasses,TMParray1(1:nContP*nPasses*24),&
            & TMParray2(1:nContP*nPasses*20))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*20),&
            & LOCALINTS(1:nContP*nPasses*15),lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*100.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ100(TMParray1(1:nPrimP*nPasses*100),&
            & TMParray2(1:nContP*nPasses*100),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*50),&
            & TMParray1(1:nContP*nPasses*30),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(5,nContP*nPasses,TMParray1(1:nContP*nPasses*30),&
            & LOCALINTS(1:nContP*nPasses*25))
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP2Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P2A2B0AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*120),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP2_maxAngA2(20,nContP*nPasses,TMParray1(1:nContP*nPasses*120),&
            & TMParray2(1:nContP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,5,Qdistance12,TMParray2(1:nContP*nPasses*100),&
            & TMParray1(1:nContP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(5,nContP*nPasses,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*35))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*80.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*35),TMParray1(1:nPrimP*nPasses*80))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ80(TMParray1(1:nPrimP*nPasses*80),&
            & TMParray2(1:nContP*nPasses*80),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*80),&
            & TMParray1(1:nContP*nPasses*72),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(4,nContP*nPasses,TMParray1(1:nContP*nPasses*72),&
            & TMParray2(1:nContP*nPasses*60))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*60),&
            & LOCALINTS(1:nContP*nPasses*45),lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*200.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*200))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ200(TMParray1(1:nPrimP*nPasses*200),&
            & TMParray2(1:nContP*nPasses*200),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P3A2B1AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*200),&
            & TMParray1(1:nContP*nPasses*180),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP3_maxAngA2(10,nContP*nPasses,TMParray1(1:nContP*nPasses*180),&
            & TMParray2(1:nContP*nPasses*150))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*150),&
            & TMParray1(1:nContP*nPasses*90),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(15,nContP*nPasses,TMParray1(1:nContP*nPasses*90),&
            & LOCALINTS(1:nContP*nPasses*75))
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*400.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP3Q3AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*400))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ400(TMParray1(1:nPrimP*nPasses*400),&
            & TMParray2(1:nContP*nPasses*400),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
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
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,15,Qdistance12,TMParray2(1:nContP*nPasses*300),&
            & TMParray1(1:nContP*nPasses*270),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(15,nContP*nPasses,TMParray1(1:nContP*nPasses*270),&
            & LOCALINTS(1:nContP*nPasses*225))
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*56))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*140.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q1AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*56),TMParray1(1:nPrimP*nPasses*140))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ140(TMParray1(1:nPrimP*nPasses*140),&
            & TMParray2(1:nContP*nPasses*140),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,4,Pdistance12,TMParray2(1:nContP*nPasses*140),&
            & TMParray1(1:nContP*nPasses*144),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(4,nContP*nPasses,TMParray1(1:nContP*nPasses*144),&
            & TMParray2(1:nContP*nPasses*100))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q1C0D1DtoC(nContP,nPasses,25,Qdistance12,TMParray2(1:nContP*nPasses*100),&
            & LOCALINTS(1:nContP*nPasses*75),lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*84))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*350.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*84),TMParray1(1:nPrimP*nPasses*350))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ350(TMParray1(1:nPrimP*nPasses*350),&
            & TMParray2(1:nContP*nPasses*350),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,10,Pdistance12,TMParray2(1:nContP*nPasses*350),&
            & TMParray1(1:nContP*nPasses*360),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(10,nContP*nPasses,TMParray1(1:nContP*nPasses*360),&
            & TMParray2(1:nContP*nPasses*250))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q2C0D2DtoC(nContP,nPasses,25,Qdistance12,TMParray2(1:nContP*nPasses*250),&
            & TMParray1(1:nContP*nPasses*150),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ2_maxAngC0(25,nContP*nPasses,TMParray1(1:nContP*nPasses*150),&
            & LOCALINTS(1:nContP*nPasses*125))
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPassestoo small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPU7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
               & TMParray2(1:nPrimP*nPasses*120))
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPasses*700.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPassestoo small',-1)
        ENDIF
#endif
        call TransferRecurrenceP4Q3AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2(1:nPrimP*nPasses*120),TMParray1(1:nPrimP*nPasses*700))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
         call PrimitiveContractionSegQ700(TMParray1(1:nPrimP*nPasses*700),&
            & TMParray2(1:nContP*nPasses*700),nPrimP,nPrimQ,nPasses,&
              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_LHS_P4A2B2AtoB(nContP*nPasses,20,Pdistance12,TMParray2(1:nContP*nPasses*700),&
            & TMParray1(1:nContP*nPasses*720),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_maxAngP4_maxAngA2(20,nContP*nPasses,TMParray1(1:nContP*nPasses*720),&
            & TMParray2(1:nContP*nPasses*500))
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call HorizontalRR_RHS_Q3C1D2DtoC(nContP,nPasses,25,Qdistance12,TMParray2(1:nContP*nPasses*500),&
            & TMParray1(1:nContP*nPasses*450),lupri)
#ifdef VAR_DEBUGICHOR
        IF(nContP*nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nContP*nPassestoo small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_maxAngQ3_maxAngC1(25,nContP*nPasses,TMParray1(1:nContP*nPasses*450),&
            & LOCALINTS(1:nContP*nPasses*375))
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_SegQ',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_SegQ
  
  
  subroutine IchorCoulombIntegral_OBS_general_sizeSegQ(TMParray1maxsize,&
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
    BasisCont3maxsize =     1*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimP)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
    BasisCont3maxsize =     4*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContP)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
    BasisCont3maxsize =     4*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContP)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
    BasisCont3maxsize =    20*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
    BasisCont3maxsize =    20*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
    BasisCont3maxsize =    35*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContP)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
    BasisCont3maxsize =     4*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContP)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
    BasisCont3maxsize =    16*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
    BasisCont3maxsize =    16*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
    BasisCont3maxsize =   140*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContP)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
    BasisCont3maxsize =     4*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nContP)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
    BasisCont3maxsize =    16*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
    BasisCont3maxsize =    16*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContP)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
    BasisCont3maxsize =   140*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContP)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContP)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContP)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContP)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContP)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContP)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,315*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,324*nContP)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
    BasisCont3maxsize =    20*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
    BasisCont3maxsize =   400*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
    BasisCont3maxsize =   400*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
    BasisCont3maxsize =   700*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContP)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
    BasisCont3maxsize =    10*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContP)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
    BasisCont3maxsize =    40*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContP)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
    BasisCont3maxsize =   100*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContP)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
    BasisCont3maxsize =    20*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContP)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
    BasisCont3maxsize =    80*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContP)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
    BasisCont3maxsize =   400*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
    BasisCont3maxsize =   200*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContP)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
    BasisCont3maxsize =   400*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContP)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
    BasisCont3maxsize =   700*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContP)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
    BasisCont3maxsize =    35*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContP)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
    BasisCont3maxsize =   140*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContP)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
    BasisCont3maxsize =   140*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContP)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContP)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
    BasisCont3maxsize =   700*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContP)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
    BasisCont3maxsize =   350*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContP)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
    BasisCont3maxsize =   700*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContP)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
    BasisCont3maxsize =  1225*nPrimB
       TMParray2maxSize = MAX(TMParray2maxSize,165*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nPrimP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContP)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContP)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_general_sizeSegQ
  subroutine PrimitiveContractionSegQ1(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    !Due to Q being segmented the Q contraction have already been done and we need to 
    !go from nPrimP to nContP
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    real(realk) :: tmp,BasisCont3(nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       tmp = 0.0E0_realk
       do iPrimA=1,nPrimA
        tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iPrimA,iPrimB,iPassQ)
       enddo
       BasisCont3(iPrimB) = tmp
      enddo
      do iContB=1,nContB
       tmp = 0.0E0_realk
       do iPrimB=1,nPrimB
        tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iPrimB)
       enddo
       AUXarrayCont(iContA,iContB,iPassQ) = tmp
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ1


  subroutine PrimitiveContractionSegQ4(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(    4,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(    4,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(    4,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,    4
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,    4
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ4


  subroutine PrimitiveContractionSegQ10(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   10,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   10,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(   10,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   10
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   10
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ10


  subroutine PrimitiveContractionSegQ20(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   20,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   20,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(   20,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   20
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   20
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ20


  subroutine PrimitiveContractionSegQ35(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   35,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   35,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(   35,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   35
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   35
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ35


  subroutine PrimitiveContractionSegQ16(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   16,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   16,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(   16,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   16
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   16
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ16


  subroutine PrimitiveContractionSegQ40(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   40,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   40,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(   40,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   40
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   40
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ40


  subroutine PrimitiveContractionSegQ80(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(   80,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(   80,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(   80,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,   80
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,   80
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ80


  subroutine PrimitiveContractionSegQ140(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  140,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  140,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(  140,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  140
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  140
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ140


  subroutine PrimitiveContractionSegQ100(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  100,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  100,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(  100,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  100
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  100
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ100


  subroutine PrimitiveContractionSegQ200(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  200,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  200,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(  200,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  200
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  200
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ200


  subroutine PrimitiveContractionSegQ350(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  350,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  350,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(  350,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  350
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  350
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ350


  subroutine PrimitiveContractionSegQ400(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  400,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  400,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(  400,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  400
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  400
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ400


  subroutine PrimitiveContractionSegQ700(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2(  700,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(  700,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3(  700,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1,  700
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1,  700
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ700


  subroutine PrimitiveContractionSegQ1225(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&
       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB,BasisCont3)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: AUXarray2( 1225,nPrimA,nPrimB,nPasses)
    real(realk),intent(inout) :: AUXarrayCont( 1225,nContA,nContB,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iPrimA,iPrimB
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: TMP
    real(realk) :: BasisCont3( 1225,nPrimB)
    do iPassQ = 1,nPasses
     do iContA=1,nContA
      do iPrimB=1,nPrimB
       do iTUV=1, 1225
        TMP = 0.0E0_realk
        do iPrimA=1,nPrimA
         tmp = tmp + ACC(iPrimA,iContA)*AUXarray2(iTUV,iPrimA,iPrimB,iPassQ)
        enddo
        BasisCont3(iTUV,iPrimB) = tmp
       enddo
      enddo
      do iContB=1,nContB
       do iTUV=1, 1225
        TMP = 0.0E0_realk
        do iPrimB=1,nPrimB
         tmp = tmp + BCC(iPrimB,iContB)*BasisCont3(iTUV,iPrimB)
        enddo
        AUXarrayCont(iTUV,iContA,iContB,iPassQ) = tmp
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSegQ1225

END MODULE IchorEriCoulombintegralOBSGeneralModSegQ
