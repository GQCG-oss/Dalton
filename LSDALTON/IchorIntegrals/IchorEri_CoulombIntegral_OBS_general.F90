MODULE IchorEriCoulombintegralOBSGeneralMod
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
use IchorEriCoulombintegralOBSGeneralModGen
use IchorEriCoulombintegralOBSGeneralModSegQ
use IchorEriCoulombintegralOBSGeneralModSegP
use IchorEriCoulombintegralOBSGeneralModSeg
use IchorEriCoulombintegralOBSGeneralModSeg1Prim
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000ModGen
use AGC_CPU_OBS_BUILDRJ000ModSeg1Prim
public :: IchorCoulombIntegral_OBS_general,IchorCoulombIntegral_OBS_general_size  
  
CONTAINS
  
  
  subroutine IchorCoulombIntegral_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
       & BasisCont1,BasisCont2,BasisCont3,IatomAPass,iatomBPass)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: nAtomsA,nAtomsB
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP*nAtomsA*nAtomsB)   !pcent(3,nPrimP)
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
    real(realk),intent(in) :: Pdistance12(3*nAtomsA*nAtomsB)  !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3,nAtomsA),Bcenter(3,nAtomsB),Ccenter(3),Dcenter(3)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    real(realk),intent(inout) :: BasisCont1(BasisCont1maxsize) !nTUVP*nTUVQ*nPrimD*nPrimA*nPrimB
    real(realk),intent(inout) :: BasisCont2(BasisCont1maxsize) !nTUVP*nTUVQ*nPrimA*nPrimB
    real(realk),intent(inout) :: BasisCont3(BasisCont1maxsize) !nTUVP*nTUVQ*nPrimB or nTUVP*nTUVQ*nPrimD
    integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  
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
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
       & BasisCont1,BasisCont2,BasisCont3,IatomAPass,iatomBPass)
   ELSEIF(Psegmented.AND.Qsegmented)THEN
    call IchorCoulombIntegral_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
       & BasisCont1,BasisCont2,BasisCont3,IatomAPass,iatomBPass)
   ELSEIF(Psegmented)THEN
    call IchorCoulombIntegral_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
       & BasisCont1,BasisCont2,BasisCont3,IatomAPass,iatomBPass)
   ELSEIF(Qsegmented)THEN
    call IchorCoulombIntegral_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
       & BasisCont1,BasisCont2,BasisCont3,IatomAPass,iatomBPass)
   ELSE
    call IchorCoulombIntegral_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&
       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&
       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize,&
       & BasisCont1,BasisCont2,BasisCont3,IatomAPass,iatomBPass)
   ENDIF
  end subroutine IchorCoulombIntegral_OBS_general
  
  
  subroutine IchorCoulombIntegral_OBS_general_size(TMParray1maxsize,&
         & TMParray2maxsize,BasisCont1maxsize,BasisCont2maxsize,&
         & BasisCont3maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,&
         & nContP,nContQ,nPrimQP,nContQP,Psegmented,Qsegmented)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(inout) :: BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP
    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD
    logical,intent(in) :: Psegmented,Qsegmented
    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN
     call IchorCoulombIntegral_OBS_general_sizeSeg1Prim(TMParray1maxsize,&
         & TMParray2maxsize,BasisCont1maxsize,BasisCont2maxsize,&
         & BasisCont3maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ELSEIF(Psegmented.AND.Qsegmented)THEN
     call IchorCoulombIntegral_OBS_general_sizeSeg(TMParray1maxsize,&
         & TMParray2maxsize,BasisCont1maxsize,BasisCont2maxsize,&
         & BasisCont3maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ELSEIF(Psegmented)THEN
     call IchorCoulombIntegral_OBS_general_sizeSegP(TMParray1maxsize,&
         & TMParray2maxsize,BasisCont1maxsize,BasisCont2maxsize,&
         & BasisCont3maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ELSEIF(Qsegmented)THEN
     call IchorCoulombIntegral_OBS_general_sizeSegQ(TMParray1maxsize,&
         & TMParray2maxsize,BasisCont1maxsize,BasisCont2maxsize,&
         & BasisCont3maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ELSE
     call IchorCoulombIntegral_OBS_general_sizeGen(TMParray1maxsize,&
         & TMParray2maxsize,BasisCont1maxsize,BasisCont2maxsize,&
         & BasisCont3maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ENDIF
  end subroutine IchorCoulombIntegral_OBS_general_size
  
END MODULE IchorEriCoulombintegralOBSGeneralMod
