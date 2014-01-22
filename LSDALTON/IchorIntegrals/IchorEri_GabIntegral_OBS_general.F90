MODULE IchorEriGabintegralOBSGeneralMod
!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory
use IchorEriGabintegralOBSGeneralModGen
use IchorEriGabintegralOBSGeneralModSeg
use IchorprecisionModule
use IchorCommonModule
use IchorMemory
public :: IchorGabIntegral_OBS_general,IchorGabIntegral_OBS_general_size  
  
CONTAINS
  
  
  subroutine IchorGabIntegral_OBS_general(nPrimA,nPrimB,&
       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContP,pexp,ACC,BCC,&
       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,Pdistance12,PQorder,LOCALINTS,Acenter,Bcenter,&
       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisContmaxsize,BasisCont)
    implicit none
    integer,intent(in) :: nPrimP,nPasses,nPrimA,nPrimB
    integer,intent(in) :: MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContP,nTABFJW1,nTABFJW2
    integer,intent(in) :: AngmomA,AngmomB
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)
    logical,intent(in)     :: Psegmented
    real(realk),intent(in) :: pexp(nPrimP)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(inout) :: LOCALINTS(nPasses)
    real(realk),intent(in) :: integralPrefactor(nPrimP*nPrimP)
    logical,intent(in) :: PQorder
    !integralPrefactor(nPrimP,nPrimP)
    real(realk),intent(in) :: reducedExponents(nPrimP*nPrimP)
    !reducedExponents(nPrimP,nPrimP)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),intent(in) :: Acenter(3),Bcenter(3)
    logical,intent(in) :: spherical
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize,BasisContmaxsize
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    real(realk),intent(inout) :: BasisCont(BasisContmaxsize)
    IF(PQorder)THEN
       call IchorQuit('PQorder OBS general expect to get QP ordering',-1)
    ENDIF
    IF(.NOT.spherical)THEN
       call IchorQuit('cartesian not testet',-1)
    ENDIF
    
   IF(Psegmented)THEN
    call IchorGabIntegral_OBS_Seg(nPrimA,nPrimB,&
       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContP,pexp,ACC,BCC,&
       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,Pdistance12,PQorder,LOCALINTS,Acenter,Bcenter,&
       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisContmaxsize,BasisCont)
   ELSE
    call IchorGabIntegral_OBS_Gen(nPrimA,nPrimB,&
       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContP,pexp,ACC,BCC,&
       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,Pdistance12,PQorder,LOCALINTS,Acenter,Bcenter,&
       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
       & BasisContmaxsize,BasisCont)
   ENDIF
  end subroutine IchorGabIntegral_OBS_general
  
  
  subroutine IchorGabIntegral_OBS_general_size(TMParray1maxsize,&
         & TMParray2maxsize,BasisContmaxsize,AngmomA,AngmomB,nPrimP,&
         & nContP,nPrimB,Psegmented)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(inout) :: BasisContmaxsize
    integer,intent(in) :: AngmomA,AngmomB
    integer,intent(in) :: nPrimP,nContP,nPrimB
    logical,intent(in) :: Psegmented
    IF(Psegmented)THEN
     call IchorGabIntegral_OBS_general_sizeSeg(TMParray1maxsize,&
         & TMParray2maxsize,BasisContmaxsize,AngmomA,AngmomB,&
         & nPrimP,nContP,nPrimB)
    ELSE
     call IchorGabIntegral_OBS_general_sizeGen(TMParray1maxsize,&
         &TMParray2maxsize,BasisContmaxsize,AngmomA,AngmomB,&
         & nPrimP,nContP,nPrimB)
    ENDIF
  end subroutine IchorGabIntegral_OBS_general_size
  
END MODULE IchorEriGabintegralOBSGeneralMod
