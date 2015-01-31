MODULE IchorEriCoulombintegralCPUOBSGeneralModSeg1Prim2
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset containing a single primitive
use IchorEriCoulombintegralCPUMcMGeneralMod
use IchorprecisionMod
use IchorCommonMod
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000ModGen
use AGC_CPU_OBS_BUILDRJ000ModSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODASeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODBSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODDSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODCSeg1Prim
use AGC_CPU_OBS_VERTICALRECURRENCEMODAGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODBGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODCGen
use AGC_CPU_OBS_VERTICALRECURRENCEMODDGen
use AGC_CPU_OBS_TRMODAtoCSeg1Prim1
use AGC_CPU_OBS_TRMODAtoCSeg1Prim2
use AGC_CPU_OBS_TRMODAtoDSeg1Prim1
use AGC_CPU_OBS_TRMODAtoDSeg1Prim2
use AGC_CPU_OBS_TRMODBtoCSeg1Prim1
use AGC_CPU_OBS_TRMODBtoDSeg1Prim1
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
public :: ICI_CPU_OBS_Seg1Prim2
  
CONTAINS
  
  
  subroutine ICI_CPU_OBS_Seg1Prim2(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPasses,MaxPasses,IntPrint,lupri,&
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
    integer,intent(in) :: MaxPasses,IntPrint,lupri
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
    real(realk),intent(in) :: integralPrefactor(nPrimQ*nPrimP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQ*nPrimP)
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
    integer :: AngmomPQ,AngmomP,AngmomQ,I,J,la,lb,lc,ld,nsize,angmomid
    
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
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceCPUSeg1Prim1D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim2D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceCPUSeg1Prim1C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim2C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceCPUSeg1Prim1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray2(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q2DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q2CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q3DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q2CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q3CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A0B1BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q3DtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q3CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A0B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA0(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call BuildRJ000CPUSeg1Prim2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q2DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q2CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP1Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        !No reason for the Electron Transfer Recurrence Relation 
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray1(1),&
            & TMParray2(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray2(1),&
            & LOCALINTS(1))
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q2BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q1BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C1D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C1D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q3BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q2BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C2D0CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q3BtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C2D1CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call BuildRJ000CPUSeg1Prim7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q4CtoBSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A1B2BtoA(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA1(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q4C2D2CtoD(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call BuildRJ000CPUSeg1Prim3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5D(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP2Q3DtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP2_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call BuildRJ000CPUSeg1Prim4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP3Q3AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP3_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call BuildRJ000CPUSeg1Prim5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP4Q1AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q1C0D1DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & LOCALINTS(1),lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call BuildRJ000CPUSeg1Prim6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP4Q2AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q2C0D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call BuildRJ000CPUSeg1Prim7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&
               & MaxPasses,nAtomsA,nAtomsB,TMParray2(1))
        call VerticalRecurrenceCPUSeg1Prim7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TMParray2(1),Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
               & TMParray1(1))
        call TransferRecurrenceCPUP4Q3AtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&
               & TMParray1(1),TMParray2(1))
        !Primitive Contraction have already been done
        call HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQ,nContP,nPasses,nTUVQ,&
            & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS1_CPU_maxAngP4_maxAngA2(nTUVQ,nContQ,nContP,nPasses,TMParray1(1),&
            & TMParray2(1))
        call HorizontalRR_CPU_RHS_Q3C1D2DtoC(nContQ,nContP,nPasses,nOrbCompP,Qdistance12,TMParray2(1),&
            & TMParray1(1),lupri)
        call SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nOrbCompP,nContQ,nContP,nPasses,TMParray1(1),&
            & LOCALINTS(1))
    CASE DEFAULT
#ifdef VAR_OPENACC
        CALL ICHORQUIT('ICI_CPU_McM_general called with OpenACC',-1)
#endif
        call ICI_CPU_McM_general(nPrimA,nPrimB,nPrimC,nPrimD,&
           & nPrimP,nPrimQ,nPasses,MaxPasses,IntPrint,lupri,&
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
  end subroutine ICI_CPU_OBS_Seg1Prim2
  
END MODULE IchorEriCoulombintegralCPUOBSGeneralModSeg1Prim2
