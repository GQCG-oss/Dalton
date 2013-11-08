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
public :: IchorCoulombIntegral_OBS_general,IchorCoulombIntegral_OBS_general_size  
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
  
  
  subroutine IchorCoulombIntegral_OBS_general_size(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         &nPrimQP,nContQP)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nPrimQP,nContQP
    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN
     call IchorCoulombIntegral_OBS_general_sizeSeg1Prim(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         &nPrimQP,nContQP)
    ELSEIF(Psegmented.AND.Qsegmented)THEN
     call IchorCoulombIntegral_OBS_general_sizeSeg(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         &nPrimQP,nContQP)
    ELSEIF(Psegmented)THEN
     call IchorCoulombIntegral_OBS_general_sizeSegP(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         &nPrimQP,nContQP)
    ELSEIF(Qsegmented)THEN
     call IchorCoulombIntegral_OBS_general_sizeSegQ(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         &nPrimQP,nContQP)
    ELSE
     call IchorCoulombIntegral_OBS_general_sizeGen(TMParray1maxsize,&
         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         &nPrimQP,nContQP)
    ENDIF
  end subroutine IchorCoulombIntegral_OBS_general_size
  
