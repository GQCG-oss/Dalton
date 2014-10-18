MODULE IchorEriCoulombintegralCPUOBSGeneralMod
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
use IchorEriCoulombintegralCPUOBSGeneralModGen
use IchorEriCoulombintegralCPUOBSGeneralModSegQ
use IchorEriCoulombintegralCPUOBSGeneralModSegP
use IchorEriCoulombintegralCPUOBSGeneralModSeg
use IchorEriCoulombintegralCPUOBSGeneralModSeg1Prim
use IchorEriCoulombintegralCPUOBSGeneralModGen2
use IchorEriCoulombintegralCPUOBSGeneralModSegQ2
use IchorEriCoulombintegralCPUOBSGeneralModSegP2
use IchorEriCoulombintegralCPUOBSGeneralModSeg2
use IchorEriCoulombintegralCPUOBSGeneralModSeg1Prim2
use IchorEriCoulombintegralCPUOBSGeneralModGenSize
use IchorEriCoulombintegralCPUOBSGeneralModSegQSize
use IchorEriCoulombintegralCPUOBSGeneralModSegPSize
use IchorEriCoulombintegralCPUOBSGeneralModSegSize
use IchorEriCoulombintegralCPUOBSGeneralModSeg1PrimSize
use SPIchorEriCoulombintegralCPUOBSGeneralModGen
use SPIchorEriCoulombintegralCPUOBSGeneralModSegQ
use SPIchorEriCoulombintegralCPUOBSGeneralModSegP
use SPIchorEriCoulombintegralCPUOBSGeneralModSeg
use SPIchorEriCoulombintegralCPUOBSGeneralModSeg1Prim
use SPIchorEriCoulombintegralCPUOBSGeneralModGen2
use SPIchorEriCoulombintegralCPUOBSGeneralModSegQ2
use SPIchorEriCoulombintegralCPUOBSGeneralModSegP2
use SPIchorEriCoulombintegralCPUOBSGeneralModSeg2
use SPIchorEriCoulombintegralCPUOBSGeneralModSeg1Prim2
use IchorprecisionMod
use IchorCommonMod
use IchorMemory
use AGC_CPU_OBS_BUILDRJ000ModGen
use AGC_CPU_OBS_BUILDRJ000ModSeg1Prim
public :: ICI_CPU_OBS_general
  
CONTAINS
  
  
  subroutine ICI_CPU_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,&
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
       & IatomAPass,iatomBPass,UseSP)
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
    logical,intent(in)     :: Qsegmented,Psegmented,UseSP
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
!   TMP variables - allocated outside
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
    real(reals),allocatable :: Aexp_SP(:),Bexp_SP(:),Cexp_SP(:),Dexp_SP(:)
    real(reals),allocatable :: pexp_SP(:),qexp_SP(:),pcent_SP(:),qcent_SP(:)
    real(reals),allocatable :: QpreExpFac_SP(:),PpreExpFac_SP(:),TABFJW_SP(:,:)
    real(reals),allocatable :: ACC_SP(:,:),BCC_SP(:,:),CCC_SP(:,:),DCC_SP(:,:)
    real(reals),allocatable :: LOCALINTS_SP(:),integralPrefactor_SP(:)
    real(reals),allocatable :: reducedExponents_SP(:)
    real(reals) :: Qdistance12_SP(3),Ccenter_SP(3),Dcenter_SP(3)
    real(reals),allocatable :: Pdistance12_SP(:)
    real(reals),allocatable :: Acenter_SP(:,:),Bcenter_SP(:,:)
    real(reals),allocatable :: TmpArray1_SP(:),TmpArray2_SP(:)
  
    IF(PQorder)THEN
       call IchorQuit('PQorder OBS general expect to get QP ordering',-1)
    ENDIF
    IF(.NOT.spherical)THEN
       call IchorQuit('cartesian not testet',-1)
    ENDIF
    
  IF(.NOT.UseSP)THEN
   IF(.NOT.UseGeneralCode.AND.(((AngmomA.LE.2).AND.(AngmomA.GE.AngmomB)).AND.&
       ((AngmomA.GE.AngmomC).AND.(AngmomC.GE.AngmomD))))THEN
    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN
     call ICI_CPU_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ELSEIF(Psegmented.AND.Qsegmented)THEN
     call ICI_CPU_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ELSEIF(Psegmented)THEN
     call ICI_CPU_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ELSEIF(Qsegmented)THEN
     call ICI_CPU_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ELSE
     call ICI_CPU_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ENDIF
   ELSE
    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN
     call ICI_CPU_OBS_Seg1Prim2(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ELSEIF(Psegmented.AND.Qsegmented)THEN
     call ICI_CPU_OBS_Seg2(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ELSEIF(Psegmented)THEN
     call ICI_CPU_OBS_SegP2(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ELSEIF(Qsegmented)THEN
     call ICI_CPU_OBS_SegQ2(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ELSE
     call ICI_CPU_OBS_Gen2(nPrimA,nPrimB,nPrimC,nPrimD,&
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
    ENDIF
   ENDIF
   ELSE !Single Precision Routines!
  allocate(Aexp_SP(nPrimA))
  allocate(Bexp_SP(nPrimB))
  allocate(Cexp_SP(nPrimC))
  allocate(Dexp_SP(nPrimD))
  allocate(pexp_SP(nPrimP))
  allocate(qexp_SP(nPrimQ))
  allocate(pcent_SP(3*nPrimP*nAtomsA*nAtomsB))
  allocate(qcent_SP(3*nPrimQ))
  allocate(QpreExpFac_SP(nPrimQ))
  allocate(PpreExpFac_SP(nPrimP*nAtomsA*nAtomsB))
  allocate(TABFJW_SP(0:nTABFJW1,0:nTABFJW2))
  allocate(ACC_SP(nPrimA,nContA))
  allocate(BCC_SP(nPrimB,nContB))
  allocate(CCC_SP(nPrimC,nContC))
  allocate(DCC_SP(nPrimD,nContD))
  allocate(LOCALINTS_SP(localintsmaxsize))
  allocate(integralPrefactor_SP(nPrimQP))
  allocate(reducedExponents_SP(nPrimQP))
  allocate(Pdistance12_SP(3*nAtomsA*nAtomsB))
  allocate(Acenter_SP(3,nAtomsA))
  allocate(Bcenter_SP(3,nAtomsB))
  allocate(TmpArray1_SP(TMParray1maxsize))
  allocate(TmpArray2_SP(TMParray2maxsize))
  Aexp_SP = Real(Aexp,kind=reals)
  Bexp_SP = Real(Bexp,kind=reals)
  Cexp_SP = Real(Cexp,kind=reals)
  Dexp_SP = Real(Dexp,kind=reals)
  pexp_SP = Real(pexp,kind=reals)
  qexp_SP = Real(qexp,kind=reals)
  pcent_SP = Real(pcent,kind=reals)
  qcent_SP = Real(qcent,kind=reals)
  QpreExpFac_SP = Real(QpreExpFac,kind=reals)
  PpreExpFac_SP = Real(PpreExpFac,kind=reals)
  TABFJW_SP = Real(TABFJW,kind=reals)
  ACC_SP = Real(ACC,kind=reals)
  BCC_SP = Real(BCC,kind=reals)
  CCC_SP = Real(CCC,kind=reals)
  DCC_SP = Real(DCC,kind=reals)
  integralPrefactor_SP = Real(integralPrefactor,kind=reals)
  reducedExponents_SP = Real(reducedExponents,kind=reals)
  Pdistance12_SP = Real(Pdistance12,kind=reals)
  Acenter_SP = Real(Acenter,kind=reals)
  Bcenter_SP = Real(Bcenter,kind=reals)
  Qdistance12_SP = Real(Qdistance12,kind=reals)
  Ccenter_SP = Real(Ccenter,kind=reals)
  Dcenter_SP = Real(Dcenter,kind=reals)
   IF(.NOT.UseGeneralCode.AND.(((AngmomA.LE.2).AND.(AngmomA.GE.AngmomB)).AND.&
       ((AngmomA.GE.AngmomC).AND.(AngmomC.GE.AngmomD))))THEN
    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN
     call SPICI_CPU_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ELSEIF(Psegmented.AND.Qsegmented)THEN
     call SPICI_CPU_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ELSEIF(Psegmented)THEN
     call SPICI_CPU_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ELSEIF(Qsegmented)THEN
     call SPICI_CPU_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ELSE
     call SPICI_CPU_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ENDIF
   ELSE
    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN
     call SPICI_CPU_OBS_Seg1Prim2(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ELSEIF(Psegmented.AND.Qsegmented)THEN
     call SPICI_CPU_OBS_Seg2(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ELSEIF(Psegmented)THEN
     call SPICI_CPU_OBS_SegP2(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ELSEIF(Qsegmented)THEN
     call SPICI_CPU_OBS_SegQ2(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ELSE
     call SPICI_CPU_OBS_Gen2(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&
       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&
       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&
       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&
       & PQorder,LOCALINTS_SP,localintsmaxsize,&
       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&
       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&
       & IatomAPass,iatomBPass)
    ENDIF
   ENDIF
   !COPY TO DP
   LOCALINTS  = REAL(LOCALINTS_SP,KIND=realk)
  deallocate(Aexp_SP)
  deallocate(Bexp_SP)
  deallocate(Cexp_SP)
  deallocate(Dexp_SP)
  deallocate(pexp_SP)
  deallocate(qexp_SP)
  deallocate(pcent_SP)
  deallocate(qcent_SP)
  deallocate(QpreExpFac_SP)
  deallocate(PpreExpFac_SP)
  deallocate(TABFJW_SP)
  deallocate(ACC_SP)
  deallocate(BCC_SP)
  deallocate(CCC_SP)
  deallocate(DCC_SP)
  deallocate(LOCALINTS_SP)
  deallocate(integralPrefactor_SP)
  deallocate(reducedExponents_SP)
  deallocate(Pdistance12_SP)
  deallocate(Acenter_SP)
  deallocate(Bcenter_SP)
  deallocate(TmpArray1_SP)
  deallocate(TmpArray2_SP)
  ENDIF
  end subroutine ICI_CPU_OBS_general
  
  
  subroutine ICI_CPU_OBS_general_size(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,&
         & nContP,nContQ,nPrimQP,nContQP,Psegmented,Qsegmented)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP
    integer,intent(in) :: nContA,nContB,nContC,nContD
    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD
    logical,intent(in) :: Psegmented,Qsegmented
    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN
     call ICI_CPU_OBS_general_sizeSeg1Prim(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ELSEIF(Psegmented.AND.Qsegmented)THEN
     call ICI_CPU_OBS_general_sizeSeg(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ELSEIF(Psegmented)THEN
     call ICI_CPU_OBS_general_sizeSegP(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ELSEIF(Qsegmented)THEN
     call ICI_CPU_OBS_general_sizeSegQ(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ELSE
     call ICI_CPU_OBS_general_sizeGen(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    ENDIF
  end subroutine ICI_CPU_OBS_general_size
  
END MODULE IchorEriCoulombintegralCPUOBSGeneralMod
