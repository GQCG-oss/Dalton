!> @file
!> Contains the main Ichor screening integral drivers for calculation electron repulsion screening integrals
!> Ichor is the "Integral Code Hand Optimized for Rapid evaluation" 
!> Ichor is the ethereal golden fluid that is the blood of the greek gods

!> \brief Main Ichor drivers for the calculation of screening integrals 
!> based on the Obara Saika(OS)/Head-Gordon-Pople(HGP)
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorGabLoopmod
  use IchorprecisionMod
  use IchorCommonMod
  use IchorEriGabintegralOBSGeneralMod, only: IGI_OBS_general, &
       & IGI_OBS_general_size
  use IchorEriCoulombintegralCPUMcMGeneralMod, only: TmpArray3,TmpArray4,&
       & DetermineSizeTmpArray34,precalcichorsphmat,freeichorsphmat,&
       & nTmpArray3,nTmpArray4
  use IchorMemory
  use IchorGaussianGeminalMod, only: GGemOperatorCalc
  use IchorParametersMod
public :: GabIntLoop
private
CONTAINS

subroutine GabIntLoop(nPrimA,nPrimB,nPrimP,intprint,lupri,nContA,&
     & nAtomsA,nAtomsB,nContB,nContP,nTABFJW1,nTABFJW2,AngmomA,AngmomB,&
     & TMParray1maxsize,TMParray2maxsize,iBatchIndexOfTypeA,&
     & iBatchIndexOfTypeB,OutputDim1,OutputDim2,expB,ContractCoeffB,&
     & nCartOrbCompA,nCartOrbCompB,nOrbCompA,nOrbCompB,nCartOrbCompP,&
     & nOrbCompP,nTUVP,nTUV,&
     & Bcenter,expA,ContractCoeffA,Acenter,expP,TABFJW,&
     & reducedExponents,integralPrefactor,pcent,PpreExpFac,&
     & OutputStorage,Psegmented,PQorder,Spherical,TriangularLHSAtomLoop)
  implicit none
  integer,intent(in) :: nPrimA,nPrimB,nPrimP,intprint,lupri,nContA,nAtomsA
  integer,intent(in) :: nAtomsB,nContB,nContP,nTABFJW1,nTABFJW2,AngmomA,AngmomB,TMParray1maxsize,TMParray2maxsize
  integer,intent(in) :: iBatchIndexOfTypeA,iBatchIndexOfTypeB
  integer,intent(in) :: OutputDim1,OutputDim2
  integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nOrbCompA,nOrbCompB,nCartOrbCompP
  integer,intent(in) :: nOrbCompP,nTUVP,nTUV
  real(realk),intent(in) :: expB(nPrimB),ContractCoeffB(nPrimB,nContB),Bcenter(3,nAtomsB)
  real(realk),intent(in) :: expA(nPrimA),ContractCoeffA(nPrimA,nContA),Acenter(3,nAtomsA)
  real(realk),intent(in) :: expP(nPrimP),TABFJW(0:nTABFJW1,0:nTABFJW2)
  real(realk),intent(in) :: reducedExponents(nPrimP*nPrimP)
  real(realk),intent(in) :: integralPrefactor(nPrimP*nPrimP)
  real(realk),intent(inout) :: pcent(3*nPrimP),PpreExpFac(nPrimP)
!  real(realk),intent(inout) :: TmpArray1(TMParray1maxsize)
!  real(realk),intent(inout) :: TmpArray2(TMParray2maxsize)
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2)
  logical,intent(in) :: Psegmented,PQorder,Spherical,TriangularLHSAtomLoop
  !local variables
  real(realk) :: CDAB(1),BcenterSpec(3),AcenterSpec(3),Pdistance12(3)
  integer :: IatomB,iBatchA,IatomA,iBatchB,nPasses,MaxPasses,iPass,nPass
  real(realk),allocatable :: TmpArray1(:),TmpArray2(:)
  integer,allocatable :: iPassA(:),iPassB(:)
  IF(TriangularLHSAtomLoop)THEN
     allocate(iPassA(nAtomsA*nAtomsB+1/2))
     call mem_ichor_alloc(iPassA)
     allocate(iPassB(nAtomsA*nAtomsB+1/2))
     call mem_ichor_alloc(iPassB)
     iPass = 0
     DO IatomB = 1,nAtomsB
        DO IatomA = IatomB,nAtomsA
           iPass = iPass+1
           iPassA(iPass) = IatomA
           iPassB(iPass) = IatomB
        ENDDO
     ENDDO
     nPass = iPass
  ELSE
     nPass = nAtomsA*nAtomsB
  ENDIF
  allocate(TmpArray1(TMParray1maxsize))
  allocate(TmpArray2(TMParray2maxsize))
  call mem_ichor_alloc(TmpArray1) 
  call mem_ichor_alloc(TmpArray2) 
  MaxPasses = 1
  nPasses = 1
  IF(MAX(AngmomA,AngmomB).GT.MaxSpecialAngmom.OR.GGemOperatorCalc)THEN
     call DetermineSizeTmpArray34(nTUVP,nCartOrbCompP,nPrimP,nTUVP,&
          & nCartOrbCompP,nPrimP,1,&
          & AngmomA,AngmomB,AngmomA,AngmomB,AngmomA+AngmomB,AngmomA+AngmomB,&
          & AngmomA+AngmomB+AngmomA+AngmomB)
     allocate(TmpArray3(nTmpArray3))
     call mem_ichor_alloc(TmpArray3)
     allocate(TmpArray4(nTmpArray4))
     call mem_ichor_alloc(TmpArray4)
     CALL PreCalciChorSPHMAT(MAX(AngmomA,AngmomB))
  ENDIF
!$OMP PARALLEL DEFAULT(none) &
!$OMP PRIVATE(iPass,iAtomA,iAtomB,iBatchB,BcenterSpec,iBatchA,AcenterSpec) &
!$OMP SHARED(nPrimA,nPrimB,nPrimP,intprint,lupri,nContA,nAtomsA,&
!$OMP        nAtomsB,nContB,nContP,nTABFJW1,nTABFJW2,AngmomA,AngmomB,&
!$OMP        TMParray1maxsize,TMParray2maxsize,iBatchIndexOfTypeA,&
!$OMP        iBatchIndexOfTypeB,OutputDim1,OutputDim2,&
!$OMP        expB,ContractCoeffB,Bcenter,expA,ContractCoeffA,Acenter,&
!$OMP        expP,TABFJW,reducedExponents,integralPrefactor,nPass,&
!$OMP        pcent,PpreExpFac,OutputStorage,iPassA,iPassB,CDAB,&
!$OMP        Psegmented,PQorder,Spherical,TriangularLHSAtomLoop,&
!$OMP        nCartOrbCompA,nCartOrbCompB,nOrbCompA,nOrbCompB,&
!$OMP        nCartOrbCompP,nOrbCompP,nTUVP,nTUV,&
!$OMP        Pdistance12,nPasses,MaxPasses,TmpArray1,TmpArray2)
  DO IPass = 1,nPass
     IF(TriangularLHSAtomLoop)THEN
        IatomA = iPassA(iPass)
        IatomB = iPassB(iPass)
     ELSE
        IatomA = iPass - ((iPass-1)/nAtomsA)*nAtomsA
        IatomB = (iPass-1)/nAtomsA+1
     ENDIF
     iBatchB = iBatchIndexOfTypeB + IatomB
     BcenterSpec(1) = Bcenter(1,IatomB)
     BcenterSpec(2) = Bcenter(2,IatomB)
     BcenterSpec(3) = Bcenter(3,IatomB)
     iBatchA = iBatchIndexOfTypeA+IatomA
     AcenterSpec(1) = Acenter(1,IatomA)
     AcenterSpec(2) = Acenter(2,IatomA)
     AcenterSpec(3) = Acenter(3,IatomA)
!$OMP SINGLE
     CALL Build_qcent_Qdistance12_QpreExpFac(nPrimA,nPrimB,&
          & nContA,nContB,expA,expB,AcenterSpec,BcenterSpec,ContractCoeffA,&
          & ContractCoeffB,PSegmented,&
          & pcent,Pdistance12,PpreExpFac,INTPRINT)
!$OMP END SINGLE
!$OMP BARRIER
     call IGI_OBS_general(nPrimA,nPrimB,nPrimP,&
          & intprint,lupri,nContA,nContB,nContP,expP,&
          & ContractCoeffA,ContractCoeffB,&
          & nOrbCompA,nOrbCompB,nCartOrbCompA,nCartOrbCompB,&
          & nCartOrbCompP,nOrbCompP,nTUVP,nTUV,&
          & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
          & expA,expB,Psegmented,reducedExponents,integralPrefactor,&
          & AngmomA,AngmomB,Pdistance12,PQorder,CDAB,&
          & AcenterSpec,BcenterSpec,Spherical,&
          & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)
!$OMP SINGLE
     OutputStorage(iBatchA,iBatchB) = CDAB(1)
!$OMP END SINGLE
  ENDDO
!$OMP END PARALLEL
  IF(MAX(AngmomA,AngmomB).GT.MaxSpecialAngmom.OR.GGemOperatorCalc)THEN
    call mem_ichor_dealloc(TmpArray3)
    deallocate(TmpArray3)
    call mem_ichor_dealloc(TmpArray4)
    deallocate(TmpArray4)
    call FreeIchorSPHMAT()
  ENDIF
  call mem_ichor_dealloc(TmpArray1)
  deallocate(TmpArray1)
  call mem_ichor_dealloc(TmpArray2)
  deallocate(TmpArray2)
  IF(TriangularLHSAtomLoop)THEN
     call mem_ichor_dealloc(iPassA)
     deallocate(iPassA)
     call mem_ichor_dealloc(iPassB)
     deallocate(iPassB)
  ENDIF
end subroutine GabIntLoop

END MODULE IchorGabLoopmod
