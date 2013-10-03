MODULE IchorEriCoulombintegralOBSGeneralMod
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_OBS_VERTICALRECURRENCEMODA
use AGC_OBS_VERTICALRECURRENCEMODD
use AGC_OBS_VERTICALRECURRENCEMODC
use AGC_OBS_TRANSFERRECURRENCEMODAtoC
use AGC_OBS_TRANSFERRECURRENCEMODDtoA
use AGC_OBS_TRANSFERRECURRENCEMODCtoA
use AGC_OBS_HorizontalRecurrenceLHSModAtoB
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
!   Local variables 
    real(realk),pointer :: squaredDistance(:)
    integer :: AngmomPQ,AngmomP,AngmomQ,I,J,nContQP,la,lb,lc,ld,nsize,angmomid
    real(realk),pointer :: RJ000(:),OUTPUTinterest(:)
  
 IF(nAtomsC*nAtomsD.NE.nPasses)Call ichorquit('nPass error',-1)
  
!IF(.TRUE.)THEN
!    call interest_initialize()
!
!    nsize = size(CDAB)*10
!    call mem_ichor_alloc(OUTPUTinterest,nsize)
!    nsize = nPrimQP
!    
!    
!    la = AngmomA+1
!    lb = AngmomB+1
!    lc = AngmomC+1
!    ld = AngmomD+1
!    call interest_eri(OUTPUTinterest,nsize,&
!       & la,Aexp,Acenter(1),Acenter(2),Acenter(3),ACC,&
!         & lb,Bexp,Bcenter(1),Bcenter(2),Bcenter(3),BCC,&
!         & lc,Cexp,Ccenter(1,1),Ccenter(2,1),Ccenter(3,1),CCC,&
!         & ld,Dexp,Dcenter(1,1),Dcenter(2,1),Dcenter(3,1),DCC,&
!         & lupri)!,&
!         !         & .false.)
!    write(lupri,*)'OUTPUTinterest',OUTPUTinterest
!ENDIF
 
    
    IF(PQorder)THEN
       call IchorQuit('PQorder OBS general expect to get QP ordering',-1)
    ENDIF
    IF(.NOT.spherical)THEN
       call IchorQuit('cartesian not testet',-1)
    ENDIF
    
    
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,   1,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,   1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        !no need for LHS Horizontal recurrence relations a simply copy
        !no Spherical Transformation LHS needed
        CDAB = TMParray1
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call VerticalRecurrence1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,   4,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,   4,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        !no Spherical Transformation LHS needed
        CDAB = TMParray2
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q1C1D0CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
        call VerticalRecurrence3C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP1Q2CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P1A1B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   3,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        !no Spherical Transformation LHS needed
        CDAB = TMParray2
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q1AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P2A1B1AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        !no Spherical Transformation LHS needed
        call HorizontalRR_RHS_Q2C1D1CtoD(nContQP,nPasses,   9,Qdistance12,TMParray1,CDAB     ,lupri)
        !no Spherical Transformation RHS needed
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call VerticalRecurrence2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q2C2D0CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
        call VerticalRecurrence5C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q3CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
        call VerticalRecurrence6C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP2Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P2A2B0AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(   5,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call VerticalRecurrence3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q3C2D1CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
        call VerticalRecurrence7C(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        call TransferRecurrenceP3Q4CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P3A2B1AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  15,nContQP*nPasses,TMParray1,CDAB     )
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call VerticalRecurrence4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  35,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  35,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
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
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call HorizontalRR_LHS_P4A2B2AtoB(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call HorizontalRR_RHS_Q4C2D2CtoD(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  25,nContQP*nPasses,TMParray1,CDAB     )
    CASE DEFAULT
        CALL ICHORQUIT('Unknown Case in IchorCoulombIntegral_OBS_general',-1)
    END SELECT
  end subroutine IchorCoulombIntegral_OBS_general
  
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
    TMParray2maxSize = 0
    TMParray1maxSize = 0
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,1*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,1*nPrimQP)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,4*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,3*nContQP)
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
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,9*nContQP)
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
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,10*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nContQP)
       TMParray2maxSize = MAX(TMParray2maxSize,6*nContQP)
       TMParray1maxSize = MAX(TMParray1maxSize,5*nContQP)
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

  SUBROUTINE buildRJ000_general(nPasses,nPrimQ,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&
       & TABFJW,RJ000,JMAX,Pcent,Qcent)
    IMPLICIT NONE
    INTEGER,intent(in)         :: nPrimP,nPrimQ,Jmax,nTABFJW1,nTABFJW2,nPasses
    REAL(REALK),intent(in)     :: reducedExponents(nPrimQ,nPrimP)
    REAL(REALK),intent(in)     :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
    REAL(REALK),intent(in)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrimQ*nPrimP,nPasses)
    !
    REAL(REALK)     :: D2JP36,WVAL
    REAL(REALK),PARAMETER :: HALF =0.5E0_realk,D1=1E0_realk,D2 = 2E0_realk, D4 = 4E0_realk, D100=100E0_realk
    Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk
    REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
    REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk
    Integer :: IPNT,J
    Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
    REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
    REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
    REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
    REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
    Real(realk), parameter :: PI=3.14159265358979323846E0_realk
    REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
    REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
    REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
    Real(realk) :: W2,W3,R
    REAL(REALK), PARAMETER :: SMALL = 1E-15_realk
    REAL(REALK) :: PX,PY,PZ,PQX,PQY,PQZ,squaredDistance
    integer :: iPassQ,iPQ,iPrimP,iPrimQ
    !make for different values of JMAX => loop unroll  
    !sorting? 
    D2JP36 = 2*JMAX + 36
    DO iPassQ=1, nPasses
      ipq = 0
      DO iPrimP=1, nPrimP
        px = Pcent(1,iPrimP)
        py = Pcent(2,iPrimP)
        pz = Pcent(3,iPrimP)
        DO iPrimQ=1, nPrimQ
          ipq = ipq + 1
          pqx = px - Qcent(1,iPrimQ,iPassQ)
          pqy = py - Qcent(2,iPrimQ,iPassQ)
          pqz = pz - Qcent(3,iPrimQ,iPassQ)
          squaredDistance = pqx*pqx+pqy*pqy+pqz*pqz
          WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
          !  0 < WVAL < 0.000001
          IF (ABS(WVAL) .LT. SMALL) THEN
             RJ000(0,ipq,ipassq) = D1
             DO J=1,JMAX
                RJ000(J,ipq,ipassq)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
             ENDDO
             !  0 < WVAL < 12 
          ELSE IF (WVAL .LT. D12) THEN
             IPNT = NINT(D100*WVAL)
             WDIFF = WVAL - TENTH*IPNT
             W2    = WDIFF*WDIFF
             W3    = W2*WDIFF
             W2    = W2*COEF2
             W3    = W3*COEF3
             DO J=0,JMAX
                R = TABFJW(J,IPNT)
                R = R -TABFJW(J+1,IPNT)*WDIFF
                R = R + TABFJW(J+2,IPNT)*W2
                R = R + TABFJW(J+3,IPNT)*W3
                RJ000(J,ipq,ipassq) = R
             ENDDO
             !  12 < WVAL <= (2J+36) 
          ELSE IF (WVAL.LE.D2JP36) THEN
             REXPW = HALF*EXP(-WVAL)
             RWVAL = D1/WVAL
             GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
             RJ000(0,ipq,ipassq) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
             DO J=1,JMAX
                RJ000(J,ipq,ipassq) = RWVAL*((J - HALF)*RJ000(J-1,ipq,ipassq)-REXPW)
             ENDDO
             !  (2J+36) < WVAL 
          ELSE
             RWVAL = PID4/WVAL
             RJ000(0,ipq,ipassq) = SQRT(RWVAL)
             RWVAL = RWVAL*PID4I
             DO J = 1, JMAX
                RJ000(J,ipq,ipassq) = RWVAL*(J - HALF)*RJ000(J-1,ipq,ipassq)
             ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE buildRJ000_general
  
END MODULE IchorEriCoulombintegralOBSGeneralMod
