MODULE IchorEriCoulombintegralOBSGeneralMod
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_OBS_VERTICALRECURRENCEMOD
use AGC_OBS_TRANSFERRECURRENCEMOD
use AGC_OBS_HorizontalRecurrenceLHSMod
use AGC_OBS_HorizontalRecurrenceRHSMod
use AGC_OBS_Sphcontract1Mod
use AGC_OBS_Sphcontract2Mod
  
private   
public :: IchorCoulombIntegral_OBS_general  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical)
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
!   Local variables 
    real(realk),pointer :: squaredDistance(:)
    integer :: AngmomPQ,AngmomP,AngmomQ,I,J,nContQP,la,lb,lc,ld,nsize
    real(realk),pointer :: RJ000(:),TMParray1(:),TMParray2(:),OUTPUTinterest(:)
  
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
    IF(AngmomA.EQ.  0)THEN
     IF(AngmomB.EQ.  0)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        !This is the Angmom(A= 0,B= 0,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,1*nPrimQP*nPasses)
        call VerticalRecurrence0(nPasses,nPrimP,nPrimQ,&
               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&
               & PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,1*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,   1,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,   1,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !no need for LHS Horizontal recurrence relations a simply copy
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 0,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 0,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 0,C= 1,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 0,C= 1,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 0,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 0,C= 2,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 0,C= 2,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 0,C= 2,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ENDIF ! C if statement
     ELSEIF(AngmomB.EQ.  1)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 0,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 1,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 1,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 2,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 2,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 1,C= 2,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ENDIF ! C if statement
     ELSEIF(AngmomB.EQ.  2)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 0,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 1,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 1,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 2,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 2,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 0,B= 2,C= 2,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ENDIF ! C if statement
     ENDIF ! B if statement
    ELSEIF(AngmomA.EQ.  1)THEN
     IF(AngmomB.EQ.  0)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        !This is the Angmom(A= 1,B= 0,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,4*nPrimQP*nPasses)
        call VerticalRecurrence1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,4*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,   4,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,   4,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,3*nContQP*nPasses)
        call HorizontalRR_LHS_P1A1B0(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 0,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 0,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        !This is the Angmom(A= 1,B= 0,C= 1,D= 0) combi
        call mem_ichor_alloc(TMParray2,10*nPrimQP*nPasses)
        call VerticalRecurrence2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,16*nPrimQP*nPasses)
        call TransferRecurrenceP1Q1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,16*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  16,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,12*nContQP*nPasses)
        call HorizontalRR_LHS_P1A1B0(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,9*nContQP*nPasses)
        call HorizontalRR_RHS_Q1C1D0(nContQP,nPasses,   3,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 0,C= 1,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 0,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 0,C= 2,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 0,C= 2,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 0,C= 2,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ENDIF ! C if statement
     ELSEIF(AngmomB.EQ.  1)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        !This is the Angmom(A= 1,B= 1,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,10*nPrimQP*nPasses)
        call VerticalRecurrence2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,10*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,9*nContQP*nPasses)
        call HorizontalRR_LHS_P2A1B1(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 1,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 1,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        !This is the Angmom(A= 1,B= 1,C= 1,D= 0) combi
        call mem_ichor_alloc(TMParray2,20*nPrimQP*nPasses)
        call VerticalRecurrence3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,40*nPrimQP*nPasses)
        call TransferRecurrenceP2Q1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,40*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,36*nContQP*nPasses)
        call HorizontalRR_LHS_P2A1B1(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,27*nContQP*nPasses)
        call HorizontalRR_RHS_Q1C1D0(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
       ELSEIF(AngmomD.EQ.  1)THEN
        !This is the Angmom(A= 1,B= 1,C= 1,D= 1) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,100*nPrimQP*nPasses)
        call TransferRecurrenceP2Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,100*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,90*nContQP*nPasses)
        call HorizontalRR_LHS_P2A1B1(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,81*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C1D1(nContQP,nPasses,   9,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 1,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 1,C= 2,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 1,C= 2,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 1,C= 2,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ENDIF ! C if statement
     ELSEIF(AngmomB.EQ.  2)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 0,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 1,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 1,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 2,D= 0 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 2,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 1,B= 2,C= 2,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ENDIF ! C if statement
     ENDIF ! B if statement
    ELSEIF(AngmomA.EQ.  2)THEN
     IF(AngmomB.EQ.  0)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,10*nPrimQP*nPasses)
        call VerticalRecurrence2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,10*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,6*nContQP*nPasses)
        call HorizontalRR_LHS_P2A2B0(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray1,5*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   1,nContQP*nPasses,TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 0,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,10*nPrimQP*nPasses)
        call VerticalRecurrence2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,10*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  10,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,6*nContQP*nPasses)
        call HorizontalRR_LHS_P2A2B0(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 0,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 0,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call mem_ichor_alloc(TMParray2,20*nPrimQP*nPasses)
        call VerticalRecurrence3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,40*nPrimQP*nPasses)
        call TransferRecurrenceP2Q1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,40*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,24*nContQP*nPasses)
        call HorizontalRR_LHS_P2A2B0(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,20*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP2_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,15*nContQP*nPasses)
        call HorizontalRR_RHS_Q1C1D0(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 0,C= 1,D= 0) combi
        call mem_ichor_alloc(TMParray2,20*nPrimQP*nPasses)
        call VerticalRecurrence3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,40*nPrimQP*nPasses)
        call TransferRecurrenceP2Q1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,40*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  40,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,24*nContQP*nPasses)
        call HorizontalRR_LHS_P2A2B0(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,18*nContQP*nPasses)
        call HorizontalRR_RHS_Q1C1D0(nContQP,nPasses,   6,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,100*nPrimQP*nPasses)
        call TransferRecurrenceP2Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,100*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,60*nContQP*nPasses)
        call HorizontalRR_LHS_P2A2B0(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,50*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,45*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C1D1(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 0,C= 1,D= 1) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,100*nPrimQP*nPasses)
        call TransferRecurrenceP2Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,100*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,60*nContQP*nPasses)
        call HorizontalRR_LHS_P2A2B0(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,54*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C1D1(nContQP,nPasses,   6,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 0,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,100*nPrimQP*nPasses)
        call TransferRecurrenceP2Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,100*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,60*nContQP*nPasses)
        call HorizontalRR_LHS_P2A2B0(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,50*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP2_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,30*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C2D0(nContQP,nPasses,   5,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation RHS
        call mem_ichor_alloc(TMParray2,25*nContQP*nPasses)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(   5,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ELSE
        !This is the Angmom(A= 2,B= 0,C= 2,D= 0) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,100*nPrimQP*nPasses)
        call TransferRecurrenceP2Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,100*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 100,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,60*nContQP*nPasses)
        call HorizontalRR_LHS_P2A2B0(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,36*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C2D0(nContQP,nPasses,   6,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 0,C= 2,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 0,C= 2,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ENDIF ! C if statement
     ELSEIF(AngmomB.EQ.  1)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,20*nPrimQP*nPasses)
        call VerticalRecurrence3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,20*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,18*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray1,15*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   1,nContQP*nPasses,TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 1,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,20*nPrimQP*nPasses)
        call VerticalRecurrence3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,20*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  20,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,18*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 1,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 1,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,80*nPrimQP*nPasses)
        call TransferRecurrenceP3Q1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,80*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,72*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,60*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP3_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,45*nContQP*nPasses)
        call HorizontalRR_RHS_Q1C1D0(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 1,C= 1,D= 0) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,80*nPrimQP*nPasses)
        call TransferRecurrenceP3Q1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,80*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,  80,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,72*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,54*nContQP*nPasses)
        call HorizontalRR_RHS_Q1C1D0(nContQP,nPasses,  18,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call mem_ichor_alloc(TMParray2,56*nPrimQP*nPasses)
        call VerticalRecurrence5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,200*nPrimQP*nPasses)
        call TransferRecurrenceP3Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,200*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,180*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,150*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,135*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C1D1(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 1,C= 1,D= 1) combi
        call mem_ichor_alloc(TMParray2,56*nPrimQP*nPasses)
        call VerticalRecurrence5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,200*nPrimQP*nPasses)
        call TransferRecurrenceP3Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,200*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,180*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,162*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C1D1(nContQP,nPasses,  18,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 1,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call mem_ichor_alloc(TMParray2,56*nPrimQP*nPasses)
        call VerticalRecurrence5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,200*nPrimQP*nPasses)
        call TransferRecurrenceP3Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,200*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,180*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,150*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,90*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C2D0(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation RHS
        call mem_ichor_alloc(TMParray2,75*nContQP*nPasses)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  15,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ELSE
        !This is the Angmom(A= 2,B= 1,C= 2,D= 0) combi
        call mem_ichor_alloc(TMParray2,56*nPrimQP*nPasses)
        call VerticalRecurrence5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,200*nPrimQP*nPasses)
        call TransferRecurrenceP3Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,200*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 200,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,180*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,108*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C2D0(nContQP,nPasses,  18,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call mem_ichor_alloc(TMParray2,84*nPrimQP*nPasses)
        call VerticalRecurrence6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,400*nPrimQP*nPasses)
        call TransferRecurrenceP3Q3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,400*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,360*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,300*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP3_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,270*nContQP*nPasses)
        call HorizontalRR_RHS_Q3C2D1(nContQP,nPasses,  15,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation RHS
        call mem_ichor_alloc(TMParray2,225*nContQP*nPasses)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  15,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ELSE
        !This is the Angmom(A= 2,B= 1,C= 2,D= 1) combi
        call mem_ichor_alloc(TMParray2,84*nPrimQP*nPasses)
        call VerticalRecurrence6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,400*nPrimQP*nPasses)
        call TransferRecurrenceP3Q3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,400*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 400,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,360*nContQP*nPasses)
        call HorizontalRR_LHS_P3A2B1(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,324*nContQP*nPasses)
        call HorizontalRR_RHS_Q3C2D1(nContQP,nPasses,  18,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 1,C= 2,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ENDIF ! C if statement
     ELSEIF(AngmomB.EQ.  2)THEN
      IF(AngmomC.EQ.  0)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,35*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  35,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  35,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,36*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray1,25*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   1,nContQP*nPasses,TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 2,C= 0,D= 0) combi
        call mem_ichor_alloc(TMParray2,35*nPrimQP*nPasses)
        call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !No reason for the Electron Transfer Recurrence Relation 
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray1,35*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray2,TMParray1,  35,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray2,TMParray1,  35,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray2)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,36*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,   1,Pdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation LHS needed
        !no need for RHS Horizontal recurrence relations 
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 2,C= 0,D= 1 combi not implemented",-1)
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 2,C= 0,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  1)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call mem_ichor_alloc(TMParray2,56*nPrimQP*nPasses)
        call VerticalRecurrence5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,140*nPrimQP*nPasses)
        call TransferRecurrenceP4Q1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,140*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,144*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,100*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP4_maxAngA2(   4,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,75*nContQP*nPasses)
        call HorizontalRR_RHS_Q1C1D0(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 2,C= 1,D= 0) combi
        call mem_ichor_alloc(TMParray2,56*nPrimQP*nPasses)
        call VerticalRecurrence5(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,140*nPrimQP*nPasses)
        call TransferRecurrenceP4Q1(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,140*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 140,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,144*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,   4,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,108*nContQP*nPasses)
        call HorizontalRR_RHS_Q1C1D0(nContQP,nPasses,  36,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call mem_ichor_alloc(TMParray2,84*nPrimQP*nPasses)
        call VerticalRecurrence6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,350*nPrimQP*nPasses)
        call TransferRecurrenceP4Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,350*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,360*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,250*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,225*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C1D1(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation RHS needed
        CDAB = TMParray1
        call mem_ichor_dealloc(TMParray1)
        ELSE
        !This is the Angmom(A= 2,B= 2,C= 1,D= 1) combi
        call mem_ichor_alloc(TMParray2,84*nPrimQP*nPasses)
        call VerticalRecurrence6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,350*nPrimQP*nPasses)
        call TransferRecurrenceP4Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,350*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,360*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,324*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C1D1(nContQP,nPasses,  36,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  2)THEN
        CALL ICHORQUIT("Angmom A= 2,B= 2,C= 1,D= 2 combi not implemented",-1)
       ENDIF ! D if statement
      ELSEIF(AngmomC.EQ.  2)THEN
       IF(AngmomD.EQ.  0)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call mem_ichor_alloc(TMParray2,84*nPrimQP*nPasses)
        call VerticalRecurrence6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,350*nPrimQP*nPasses)
        call TransferRecurrenceP4Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,350*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,360*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,250*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  10,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,150*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C2D0(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation RHS
        call mem_ichor_alloc(TMParray2,125*nContQP*nPasses)
        call SphericalContractOBS2_maxAngQ2_maxAngC2(  25,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ELSE
        !This is the Angmom(A= 2,B= 2,C= 2,D= 0) combi
        call mem_ichor_alloc(TMParray2,84*nPrimQP*nPasses)
        call VerticalRecurrence6(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,350*nPrimQP*nPasses)
        call TransferRecurrenceP4Q2(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,350*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 350,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,360*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,  10,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,216*nContQP*nPasses)
        call HorizontalRR_RHS_Q2C2D0(nContQP,nPasses,  36,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  1)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call mem_ichor_alloc(TMParray2,120*nPrimQP*nPasses)
        call VerticalRecurrence7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,700*nPrimQP*nPasses)
        call TransferRecurrenceP4Q3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,700*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,720*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,500*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  20,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,450*nContQP*nPasses)
        call HorizontalRR_RHS_Q3C2D1(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation RHS
        call mem_ichor_alloc(TMParray2,375*nContQP*nPasses)
        call SphericalContractOBS2_maxAngQ3_maxAngC2(  25,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ELSE
        !This is the Angmom(A= 2,B= 2,C= 2,D= 1) combi
        call mem_ichor_alloc(TMParray2,120*nPrimQP*nPasses)
        call VerticalRecurrence7(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,700*nPrimQP*nPasses)
        call TransferRecurrenceP4Q3(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,700*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2, 700,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,720*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,  20,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,648*nContQP*nPasses)
        call HorizontalRR_RHS_Q3C2D1(nContQP,nPasses,  36,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ELSEIF(AngmomD.EQ.  2)THEN
        IF(spherical)THEN
        !This is the Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call mem_ichor_alloc(TMParray2,165*nPrimQP*nPasses)
        call VerticalRecurrence8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,1225*nPrimQP*nPasses)
        call TransferRecurrenceP4Q4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,1225*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,1260*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation LHS
        call mem_ichor_alloc(TMParray2,875*nContQP*nPasses)
        call SphericalContractOBS1_maxAngP4_maxAngA2(  35,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,900*nContQP*nPasses)
        call HorizontalRR_RHS_Q4C2D2(nContQP,nPasses,  25,Qdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !Spherical Transformation RHS
        call mem_ichor_alloc(TMParray2,625*nContQP*nPasses)
        call SphericalContractOBS2_maxAngQ4_maxAngC2(  25,nContQP*nPasses,TMParray1,TMParray2)
        call mem_ichor_dealloc(TMParray1)
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ELSE
        !This is the Angmom(A= 2,B= 2,C= 2,D= 2) combi
        call mem_ichor_alloc(TMParray2,165*nPrimQP*nPasses)
        call VerticalRecurrence8(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,TMParray2)
        !Electron Transfer Recurrence Relation 
        call mem_ichor_alloc(TMParray1,1225*nPrimQP*nPasses)
        call TransferRecurrenceP4Q4(nPasses,nPrimP,nPrimQ,reducedExponents,&
               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
               & TMParray2,TMParray1)
        call mem_ichor_dealloc(TMParray2)
        nContQP = nContQ*nContP
        call mem_ichor_alloc(TMParray2,1225*nContQP*nPasses)
        IF(Qsegmented.AND.Psegmented)THEN
         call PrimitiveContractionSeg(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses)
        ELSE
         call PrimitiveContraction(TMParray1,TMParray2,1225,nPrimP,nPrimQ,nPasses,&
              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&
              & nContC,nPrimD,nContD)
        ENDIF
        call mem_ichor_dealloc(TMParray1)
        !LHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray1,1260*nContQP*nPasses)
        call HorizontalRR_LHS_P4A2B2(nContQP*nPasses,  35,Pdistance12,TMParray2,TMParray1,lupri)
        call mem_ichor_dealloc(TMParray2)
        !no Spherical Transformation LHS needed
        !RHS Horizontal recurrence relations 
        call mem_ichor_alloc(TMParray2,1296*nContQP*nPasses)
        call HorizontalRR_RHS_Q4C2D2(nContQP,nPasses,  36,Qdistance12,TMParray1,TMParray2,lupri)
        call mem_ichor_dealloc(TMParray1)
        !no Spherical Transformation RHS needed
        CDAB = TMParray2
        call mem_ichor_dealloc(TMParray2)
        ENDIF
       ENDIF ! D if statement
      ENDIF ! C if statement
     ENDIF ! B if statement
    ENDIF ! A if statement
  end subroutine IchorCoulombIntegral_OBS_general
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
                               print*,'ABCDTMP',ABCDTMP
                               iPrimQ = iPrimQ + 1
                               do iTUV=1,nTUVfull
                                  AUXarrayCont(iTUV,iContQ,iContP,iPassQ) = AUXarrayCont(iTUV,iContQ,iContP,iPassQ) + &
                                       & ABCDTMP*AUXarray2(iTUV,iPrimQ,iPrimP,iPassQ)
                                  !WRITE(lupri,'(A,I6,A,ES16.8)')'cont AUXarrayCont(',iTUV,',iContQ,iContP,iPassQ)',&
                                  !     & ABCDTMP*AUXarray2(iTUV,iPrimQ,iPrimP,iPassQ)
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
