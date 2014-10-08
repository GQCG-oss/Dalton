MODULE IchorEriCoulombintegralCPUOBSGeneralModSeg1Prim
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset containing a single primitive
use IchorEriCoulombintegralCPUMcMGeneralMod
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000ModGen
use AGC_CPU_OBS_BUILDRJ000ModSeg1Prim
use IchorGaussianGeminalMod, only: GGemOperatorCalc
use AGC_CPU_OBS_VERTICALRECURRENCEMODASeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODBSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODDSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODCSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_CPU_OBS_TRMODAtoCSeg1Prim
use AGC_CPU_OBS_TRMODAtoDSeg1Prim
use AGC_CPU_OBS_TRMODBtoCSeg1Prim
use AGC_CPU_OBS_TRMODBtoDSeg1Prim
use AGC_CPU_OBS_TRMODCtoASeg1Prim
use AGC_CPU_OBS_TRMODDtoASeg1Prim
use AGC_CPU_OBS_TRMODCtoBSeg1Prim
use AGC_CPU_OBS_TRMODDtoBSeg1Prim
use AGC_CPU_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_CPU_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_CPU_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_CPU_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_CPU_OBS_Sphcontract1Mod
use AGC_CPU_OBS_Sphcontract2Mod
  
private   
public :: ICI_CPU_OBS_Seg1Prim,ICI_CPU_OBS_general_sizeSeg1Prim  
  
CONTAINS
  
  
  subroutine ICI_CPU_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    IF(GGemOperatorCalc) AngmomID = AngmomID + 10000 !force to use general code
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPasses*1.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
               & PpreExpFac,QpreExpFac,LOCALINTS)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,1,&
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
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,1,&
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
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*81.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(1,nPasses,TMParray2,&
            & LOCALINTS)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(35,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(1,nPasses,TMParray2,&
            & LOCALINTS)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(35,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(1,nPasses,TMParray2,&
            & LOCALINTS)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q1AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(25,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q3AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(25,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*9.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*165.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*1225.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q4AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*1260.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*875.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(35,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*900.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*625.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(25,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(1,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,1,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,1,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(1,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(1,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(1,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*36.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,1,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(1,nPasses,TMParray2,&
            & LOCALINTS)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*3.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,1,&
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
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*5.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(1,nPasses,TMParray2,&
            & LOCALINTS)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*210.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*175.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(35,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*3.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*10.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*16.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*12.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*9.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q2CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP1Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*105.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*108.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,3,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(3,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*36.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*27.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,9,Qdistance12,TMParray1,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(9,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(9,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*54.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(9,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*162.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(9,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*315.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
        !no Spherical Transformation LHS needed
#ifdef VAR_DEBUGICHOR
        IF(nPasses*324.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,9,Qdistance12,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(9,nPasses,TMParray2,&
            & LOCALINTS)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*18.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,1,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1,&
            & TMParray2,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(1,nPasses,TMParray2,&
            & LOCALINTS)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*135.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*630.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(1,nPasses,35,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*525.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(35,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*540.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*4.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*20.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*40.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*24.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*20.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*15.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*50.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*30.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*25.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP2Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,5,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(5,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*5.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*35.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*80.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*72.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*60.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*45.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*200.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*180.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*150.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*90.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*400.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP3Q3AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*300.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*270.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,15,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*225.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(15,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*6.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*56.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*140.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*144.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,4,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*100.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(4,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*75.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(1,nPasses,25,Qdistance12,TMParray2,&
            & LOCALINTS,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*7.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*84.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*350.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*360.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,10,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*250.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(10,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*150.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(1,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*125.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(25,nPasses,TMParray1,&
            & LOCALINTS)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*8.GT.TMParray2maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call BuildRJ000CPUSeg1Prim7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPrimP*nPrimQ*nPasses*120.GT.TMParray1maxsize)THEN
          call ichorquit('nPrimP*nPrimQ*nPasses too small',-1)
        ENDIF
#endif
        call VerticalRecurrenceCPUSeg1Prim7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*700.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call TransferRecurrenceCPUP4Q3AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1,TMParray2)
        !Primitive Contraction have already been done
#ifdef VAR_DEBUGICHOR
        IF(nPasses*720.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(1,nPasses,20,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*500.GT.TMParray2maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(20,nPasses,TMParray1,&
            & TMParray2)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*450.GT.TMParray1maxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(1,nPasses,25,Qdistance12,TMParray2,&
            & TMParray1,lupri)
#ifdef VAR_DEBUGICHOR
        IF(nPasses*375.GT.LOCALINTSmaxsize)THEN
          call ichorquit('nPasses too small',-1)
        ENDIF
#endif
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(25,nPasses,TMParray1,&
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
  end subroutine ICI_CPU_OBS_Seg1Prim
  
  subroutine ICI_CPU_OBS_general_sizeSeg1Prim(TMParray1maxsize,&
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
    IF(GGemOperatorCalc) AngmomID = AngmomID + 10000 !force to use general code
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
    TMParray2maxSize = MAX(TMParray2maxSize,1)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5)
       TMParray1maxSize = MAX(TMParray1maxSize,35)
       TMParray2maxSize = MAX(TMParray2maxSize,36)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140)
       TMParray1maxSize = MAX(TMParray1maxSize,105)
       TMParray2maxSize = MAX(TMParray2maxSize,108)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,24)
       TMParray2maxSize = MAX(TMParray2maxSize,20)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,24)
       TMParray2maxSize = MAX(TMParray2maxSize,20)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,120)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,120)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,210)
       TMParray2maxSize = MAX(TMParray2maxSize,175)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140)
       TMParray1maxSize = MAX(TMParray1maxSize,105)
       TMParray2maxSize = MAX(TMParray2maxSize,108)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,36)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,36)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,162)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,162)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,315)
       TMParray2maxSize = MAX(TMParray2maxSize,324)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,72)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,72)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700)
       TMParray1maxSize = MAX(TMParray1maxSize,630)
       TMParray2maxSize = MAX(TMParray2maxSize,525)
       TMParray1maxSize = MAX(TMParray1maxSize,540)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,24)
       TMParray2maxSize = MAX(TMParray2maxSize,20)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,24)
       TMParray2maxSize = MAX(TMParray2maxSize,20)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,120)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,120)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,210)
       TMParray2maxSize = MAX(TMParray2maxSize,175)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,72)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,72)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700)
       TMParray1maxSize = MAX(TMParray1maxSize,630)
       TMParray2maxSize = MAX(TMParray2maxSize,525)
       TMParray1maxSize = MAX(TMParray1maxSize,540)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5)
       TMParray1maxSize = MAX(TMParray1maxSize,35)
       TMParray2maxSize = MAX(TMParray2maxSize,36)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140)
       TMParray1maxSize = MAX(TMParray1maxSize,144)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,250)
       TMParray1maxSize = MAX(TMParray1maxSize,150)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140)
       TMParray1maxSize = MAX(TMParray1maxSize,144)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,250)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700)
       TMParray1maxSize = MAX(TMParray1maxSize,720)
       TMParray2maxSize = MAX(TMParray2maxSize,500)
       TMParray1maxSize = MAX(TMParray1maxSize,450)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,250)
       TMParray1maxSize = MAX(TMParray1maxSize,150)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700)
       TMParray1maxSize = MAX(TMParray1maxSize,720)
       TMParray2maxSize = MAX(TMParray2maxSize,500)
       TMParray1maxSize = MAX(TMParray1maxSize,450)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,9*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,165*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225)
       TMParray1maxSize = MAX(TMParray1maxSize,1260)
       TMParray2maxSize = MAX(TMParray2maxSize,875)
       TMParray1maxSize = MAX(TMParray1maxSize,900)
    CASE DEFAULT
     call ICI_CPU_McM_general_size(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP,&
         & .TRUE.,.TRUE.)
    END SELECT
  end subroutine ICI_CPU_OBS_general_sizeSeg1Prim
  
END MODULE IchorEriCoulombintegralCPUOBSGeneralModSeg1Prim
