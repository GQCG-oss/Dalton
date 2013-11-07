MODULE IchorEriCoulombintegralOBSGeneralMod
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_OBS_VERTICALRECURRENCEMODA
use AGC_OBS_VERTICALRECURRENCEMODASegQ
use AGC_OBS_VERTICALRECURRENCEMODASegP
use AGC_OBS_VERTICALRECURRENCEMODASeg
use AGC_OBS_VERTICALRECURRENCEMODASeg1Prim
use AGC_OBS_VERTICALRECURRENCEMODB
use AGC_OBS_VERTICALRECURRENCEMODBSegQ
use AGC_OBS_VERTICALRECURRENCEMODBSegP
use AGC_OBS_VERTICALRECURRENCEMODBSeg
use AGC_OBS_VERTICALRECURRENCEMODBSeg1Prim
use AGC_OBS_VERTICALRECURRENCEMODD
use AGC_OBS_VERTICALRECURRENCEMODDSegQ
use AGC_OBS_VERTICALRECURRENCEMODDSegP
use AGC_OBS_VERTICALRECURRENCEMODDSeg
use AGC_OBS_VERTICALRECURRENCEMODDSeg1Prim
use AGC_OBS_VERTICALRECURRENCEMODC
use AGC_OBS_VERTICALRECURRENCEMODCSegQ
use AGC_OBS_VERTICALRECURRENCEMODCSegP
use AGC_OBS_VERTICALRECURRENCEMODCSeg
use AGC_OBS_VERTICALRECURRENCEMODCSeg1Prim
use AGC_OBS_TRANSFERRECURRENCEMODAtoC
use AGC_OBS_TRANSFERRECURRENCEMODAtoD
use AGC_OBS_TRANSFERRECURRENCEMODBtoC
use AGC_OBS_TRANSFERRECURRENCEMODBtoD
use AGC_OBS_TRANSFERRECURRENCEMODCtoA
use AGC_OBS_TRANSFERRECURRENCEMODDtoA
use AGC_OBS_TRANSFERRECURRENCEMODCtoB
use AGC_OBS_TRANSFERRECURRENCEMODDtoB
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
use AGC_OBS_HorizontalRecurrenceLHSModBtoA
use AGC_OBS_HorizontalRecurrenceRHSModCtoD
use AGC_OBS_HorizontalRecurrenceRHSModDtoC
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_OBS_general,IchorCoulombIntegral_OBS_general_size  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsC,nAtomsD
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3),Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
  
    IF(nAtomsC*nAtomsD.NE.nPasses)Call ichorquit('nPass error',-1)
  
    IF(PQorder)THEN
       call IchorQuit('PQorder OBS general expect to get QP ordering',-1)
    ENDIF
    IF(.NOT.spherical)THEN
       call IchorQuit('cartesian not testet',-1)
    ENDIF
    
   IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN
    call IchorCoulombIntegral_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
   ELSEIF(Psegmented.AND.Qsegmented)THEN
    call IchorCoulombIntegral_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
   ELSEIF(Psegmented)THEN
    call IchorCoulombIntegral_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
   ELSEIF(Qsegmented)THEN
    call IchorCoulombIntegral_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
   ELSE
    call IchorCoulombIntegral_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
   ENDIF
  end subroutine IchorCoulombIntegral_OBS_general
  
  subroutine IchorCoulombIntegral_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsC,nAtomsD
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3),Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
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
        call VerticalRecurrenceSeg1Prim0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceSeg1Prim1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrenceSeg1Prim2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceSeg1Prim1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrenceSeg1Prim2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrenceSeg1Prim3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrenceSeg1Prim2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrenceSeg1Prim3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrenceSeg1Prim4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1Prim3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_Seg1Prim',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_Seg1Prim
  
  subroutine IchorCoulombIntegral_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsC,nAtomsD
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3),Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
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
        call VerticalRecurrenceSeg0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSeg2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceSeg1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrenceSeg2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceSeg1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrenceSeg2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrenceSeg3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrenceSeg2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrenceSeg3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrenceSeg4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSeg1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSeg3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        !Primitive Contraction have already been done
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_Seg',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_Seg
  
  subroutine IchorCoulombIntegral_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsC,nAtomsD
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3),Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
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
        call VerticalRecurrenceSegP0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegP1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegP2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegP3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegP4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceSegP1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrenceSegP2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceSegP1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrenceSegP2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrenceSegP3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrenceSegP2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrenceSegP3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrenceSegP4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegP1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegP2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegP3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegP so need to contract on Q side only - not implemented',-1)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_SegP',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_SegP
  
  subroutine IchorCoulombIntegral_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsC,nAtomsD
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3),Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
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
        call VerticalRecurrenceSegQ0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrenceSegQ1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrenceSegQ2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrenceSegQ1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrenceSegQ2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrenceSegQ3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrenceSegQ2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrenceSegQ3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrenceSegQ4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrenceSegQ3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        CALL ICHORQUIT('SegQ so need to contract on P side only - not implemented',-1)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_SegQ',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_SegQ
  
  subroutine IchorCoulombIntegral_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsC,nAtomsD
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimQP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimQP)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3),Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
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
        call VerticalRecurrence0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,CDAB     ,   1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrence1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,   4,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  35,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q4AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence1D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,   4,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence2D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence1C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,   4,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence2C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  35,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        !no need for LHS Horizontal recurrence relations, it would be a simply copy
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   1,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   1,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
        call VerticalRecurrence1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,   4,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A0B1BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
        call VerticalRecurrence2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A0B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA0(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence3D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence4D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence4C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   3,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   9,nContQP*nPasses,TMParray2,CDAB     )
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
        call VerticalRecurrence3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   1,nContQP*nPasses,TMParray2,CDAB     )
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
        call VerticalRecurrence4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
        call VerticalRecurrence5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
        call VerticalRecurrence6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3BtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A1B2BtoA(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA1(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
        call VerticalRecurrence5D(nPasses,nPrimP,nPrimQ,nAtomsC,nAtomsD,reducedExponents,&
               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
        call VerticalRecurrence5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q1AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q1C0D1DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
        call VerticalRecurrence6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q2AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C0D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC0(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
        call VerticalRecurrence7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP4Q3AtoD(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C1D2DtoC(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC1(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_Gen',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_Gen
  
  
  subroutine IchorCoulombIntegral_OBS_general_size(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         &nPrimQP,nContQP)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nPrimQP,nContQP
    ! local variables
    integer :: AngmomID
    
    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimQP)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,3*nContQP)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,5*nContQP)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,3*nContQP)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContQP)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,5*nContQP)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,25*nContQP)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,3*nContQP)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContQP)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContQP)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,27*nContQP)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,75*nContQP)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,5*nContQP)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,25*nContQP)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,75*nContQP)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,25*nContQP)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,75*nContQP)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,125*nContQP)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,3*nContQP)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContQP)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,16*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,12*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContQP)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,27*nContQP)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,105*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,108*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,75*nContQP)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContQP)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,27*nContQP)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,36*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,27*nContQP)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,81*nContQP)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,135*nContQP)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,54*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,162*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,135*nContQP)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,315*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,324*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,225*nContQP)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,75*nContQP)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,135*nContQP)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,225*nContQP)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,75*nContQP)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,225*nContQP)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,375*nContQP)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,5*nContQP)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,25*nContQP)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,40*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,24*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,20*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,75*nContQP)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,100*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,60*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,50*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,30*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,25*nContQP)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,75*nContQP)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,210*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,175*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,125*nContQP)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,20*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,18*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,15*nContQP)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,75*nContQP)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,80*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,72*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,60*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,45*nContQP)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,135*nContQP)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,225*nContQP)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,200*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,180*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,150*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,90*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,75*nContQP)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,400*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,300*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,270*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,225*nContQP)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,630*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,525*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,540*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,375*nContQP)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,35*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,36*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,25*nContQP)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,75*nContQP)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,125*nContQP)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,56*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,140*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,144*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,75*nContQP)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,225*nContQP)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,375*nContQP)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,84*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,350*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,360*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,250*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,150*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,125*nContQP)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,120*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,700*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,720*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,500*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,450*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,375*nContQP)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,165*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1225*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1260*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,875*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,900*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,625*nContQP)
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general_size',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_general_size

  subroutine PrimitiveContraction(AUXarray2,AUXarrayCont,nTUVfull,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nTUVfull,nPrimP,nPrimQ,nPasses,nContP,nContQ
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nTUVfull,nPrimQ,nPrimP,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nTUVfull,nContQ,nContP,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: B,ABCDTMP,ABDTMP,ABTMP
    !all passes have same ACCs,BCCs,...
    !maybe construct big CC(nPrimQP,nContQP) matrix and call dgemm nPass times
    !the construction of CC should scale as c**4*p**4 and the 
    !dgemm should scale as c**4*p**4*L**6 but hopefully with efficient FLOP count, although not quadratic matrices....
    !special for nContPQ = 1 
    !special for nContP = 1
    !special for nContQ = 1
    !special for nContA = 1 ...
    !memory should be c**4*p**4 + p**4*L**6 which is fine
    !this would be a simple sum for segmentet! or maybe the sum can be moved into the previous electron transfer reccurence
    do iPassQ = 1,nPasses
       do iContB=1,nContB
          do iContA=1,nContA
             iContP = iContA+(iContB-1)*nContA
             do iContD=1,nContD
                do iContC=1,nContC
                   iContQ = iContC+(iContD-1)*nContC
                   do iTUV=1,nTUVfull
                      AUXarrayCont(iTUV,iContQ,iContP,iPassQ) = 0.0E0_realk
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    do iPassQ = 1,nPasses
       do iContB=1,nContB
          do iPrimB=1,nPrimB
             B = BCC(iPrimB,iContB)
             do iContA=1,nContA
                iContP = iContA+(iContB-1)*nContA
                do iPrimA=1,nPrimA
                   iPrimP = iPrimA + (iPrimB-1)*nPrimA
                   ABTMP = ACC(iPrimA,iContA)*B
                   do iContD=1,nContD
                      do iPrimD=1,nPrimD
                         ABDTMP = DCC(iPrimD,iContD)*ABTMP
                         do iContC=1,nContC
                            iContQ = iContC+(iContD-1)*nContC
                            iPrimQ = (iPrimD-1)*nPrimC
                            do iPrimC=1,nPrimC
                               ABCDTMP = CCC(iPrimC,iContC)*ABDTMP
                               iPrimQ = iPrimQ + 1
                               do iTUV=1,nTUVfull
                                  AUXarrayCont(iTUV,iContQ,iContP,iPassQ) = AUXarrayCont(iTUV,iContQ,iContP,iPassQ) + &
                                       & ABCDTMP*AUXarray2(iTUV,iPrimQ,iPrimP,iPassQ)
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine PrimitiveContraction

  subroutine PrimitiveContractionSeg(AUXarray2,AUXarrayCont,nTUVfull,nPrimP,nPrimQ,nPasses)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nTUVfull,nPrimP,nPrimQ,nPasses
    real(realk),intent(in) :: AUXarray2(nTUVfull,nPrimQ,nPrimP,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nTUVfull,nPasses)
    !
    integer :: iPassQ,iTUV,iPrimQ,iPrimP
    !Maybe the sum can be moved into the previous electron transfer reccurence or Vertical 
    do iPassQ = 1,nPasses
     do iTUV=1,nTUVfull
      AUXarrayCont(iTUV,iPassQ) = 0.0E0_realk
     enddo
    enddo
    do iPassQ = 1,nPasses
     do iPrimP = 1,nPrimP
      do iPrimQ = 1,nPrimQ
       do iTUV=1,nTUVfull
        AUXarrayCont(iTUV,iPassQ) = AUXarrayCont(iTUV,iPassQ) &
                                  & + AUXarray2(iTUV,iPrimQ,iPrimP,iPassQ)
       enddo
      enddo
     enddo
    enddo
  end subroutine PrimitiveContractionSeg

  
END MODULE IchorEriCoulombintegralOBSGeneralMod
