!> @file
!> Contains the main Ichor integral drivers for calculation electron repulsion integrals
!> Ichor is the "Integral Code Hand Optimized for Rapid evaluation" 
!> Ichor is the ethereal golden fluid that is the blood of the greek gods
!> The code was written by Thomas Kjaergaard with inspiration from the 
!> McMurchie-Davidson algorithm based Thermite code written by Simen Reine and Thomas Kjaergaard and 
!> as well as the Obara-Saika algorithm based Interest code written by Michal Repisky.
!> This code is based on Obara-Saika (and Head-Gordon-Pople) 
!> The Interest code was also used for debugging purposes
!> \brief Main Ichor drivers for the calculation of integrals 
!> based on the Obara Saika(OS)/Head-Gordon-Pople(HGP)
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriLinkmod
  use IchorprecisionMod
  use IchorCommonMod
  use IchorEriCoulombintegralCPUMcMGeneralMod, only: TmpArray3,TmpArray4,DetermineSizeTmpArray34,&
       & precalcichorsphmat, freeichorsphmat,nTmpArray3,nTmpArray4
  use IchorEriCoulombintegralCPUOBSGeneralMod, only: ICI_CPU_OBS_general, &
       & ICI_CPU_OBS_general_size
  use IchorEriCoulombintegralGPUOBSGeneralMod, only: ICI_GPU_OBS_general, &
       & ICI_GPU_OBS_general_size
  use ICI_seg_seg_SSSS_mod, only: ICI_seg_seg_SSSS
  use IchorMemory
  use IchorEriToolsmod
  use IchorEriDistMod

CONTAINS

subroutine IchorTypeLinKLoop(nAtomsA,nPrimA,nContA,nOrbCompA,&
     & expA,ContractCoeffA,AngmomA,Acenter,nOrbA,&
     & nAtomsB,nPrimB,nContB,nOrbCompB,&
     & expB,ContractCoeffB,AngmomB,Bcenter,nOrbB,&
     & nAtomsC,nPrimC,nContC,nOrbCompC,&
     & expC,ContractCoeffC,AngmomC,Ccenter,nOrbC,&
     & nAtomsD,nPrimD,nContD,nOrbCompD,&
     & expD,ContractCoeffD,AngmomD,Dcenter,nOrbD,&
     & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
     & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
     & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
     & qcent,Ppreexpfac,Qpreexpfac,&
     & Qiprim1,Qiprim2,&
     & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
     & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
     & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
     & reducedExponents,integralPrefactor,&
     & PcentPass,Pdistance12Pass,PpreExpFacPass,&
     & Qdistance12,PQorder,&
     & Spherical,TMParray1maxsize,nLocalInt,&
     & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri,&
     & DmatBD,DmatAD,DmatBC,&
     & DmatAC,KmatBD,KmatAD,KmatBC,KmatAC,nDimA,nDimB,nDimC,nDimD,nDmat,&
     & ReducedDmatBD,ReducedDmatBC,AtomGAB,AtomGCD,&
     & nKetList,KetList,nBraList,BraList,nBraketList,BraketList,&
     & SameRHSaos,SameLHSaos,SameODs)
  implicit none
  logical,intent(in) :: TriangularRHSAtomLoop,TriangularLHSAtomLoop,PermuteLHSTypes
  logical,intent(in) :: Qsegmented,Psegmented,PQorder,Spherical,CSscreen
  logical,intent(in) :: TriangularODAtomLoop
  logical,intent(in) :: SameRHSaos,SameLHSaos,SameODs
  integer,intent(in) :: nTABFJW1,nTABFJW2,lupri,nDimA,nDimB,nDimC,nDimD,nDmat
  integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
  integer,intent(in) :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV
  real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2),THRESHOLD_CS
  integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)
  integer,intent(in) :: nKetList(nAtomsD),KetList(nAtomsC,nAtomsD)
  integer,intent(in) :: nBraList(nAtomsB),BraList(nAtomsA,nAtomsB)
  integer,intent(in) :: nBraketList(nAtomsD),BraketList(nAtomsB,nAtomsD)
  real(realk),intent(in) :: DmatBD(ndimB,ndimD,nDmat),DmatAD(ndimA,ndimD,nDmat)
  real(realk),intent(in) :: DmatBC(ndimB,ndimC,nDmat),DmatAC(ndimA,ndimC,nDmat)
  real(realk),intent(inout) :: KmatBD(ndimB,ndimD,nDmat),KmatAD(ndimA,ndimD,nDmat)
  real(realk),intent(inout) :: KmatBC(ndimB,ndimC,nDmat),KmatAC(ndimA,ndimC,nDmat)
  real(realk),intent(in) :: ReducedDmatBD(nAtomsB,nAtomsD),ReducedDmatBC(nAtomsB,nAtomsC)
  real(realk),intent(in) :: AtomGAB(nAtomsA,nAtomsB),AtomGCD(nAtomsC,nAtomsD)
  !D
  integer,intent(in) :: nAtomsD,nPrimD,nContD,nOrbCompD,AngmomD,nOrbD
  real(realk),intent(in) :: Dcenter(3,nAtomsD),expD(nPrimD),ContractCoeffD(nPrimD,nContD)
  !C
  integer,intent(in) :: nAtomsC,nPrimC,nContC,nOrbCompC,AngmomC,nOrbC
  real(realk),intent(in) :: Ccenter(3,nAtomsC),expC(nPrimC),ContractCoeffC(nPrimC,nContC)
  logical,intent(in) :: noScreenCD2(nAtomsC,nAtomsD)
  logical,intent(in) :: noScreenAB(nAtomsA,nAtomsB)
  !Q
  integer,intent(in) :: nContQ,nPrimQ
  real(realk),intent(inout) :: Qcent(3,nPrimQ),Qdistance12(3),QpreExpFac(nPrimQ),expQ(nPrimQ)
  !P
  integer,intent(in) :: nContP,nPrimP
  real(realk),intent(inout) :: PpreExpFac(nPrimP)
  real(realk),intent(in) :: PcentPass(3,nPrimP,natomsA*natomsB)
  real(realk),intent(in) :: Pdistance12Pass(3,natomsA*natomsB)
  real(realk),intent(in) :: PpreExpFacPass(nPrimP,natomsA*natomsB),expP(nPrimP)
  !A & B
  integer,intent(in) :: AngmomB,AngmomA,nOrbA,nOrbB
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbCompA,nOrbCompB,nContA,nContB
  integer,intent(in) :: nPrimA,nPrimB
  real(realk),intent(in) :: Acenter(3,nAtomsA),expA(nPrimA),ContractCoeffA(nPrimA,nContA)
  real(realk),intent(in) :: Bcenter(3,nAtomsB),expB(nPrimB),ContractCoeffB(nPrimB,nContB)
  !collected
  integer,intent(in) :: TotalAngmom
!  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
!  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  integer,intent(in) :: TMParray1maxsize,TMParray2maxsize,nLocalInt
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: integralPrefactor(nPrimQ,nPrimP)
  !local variables
  integer :: A,B,C,D,iC,iB,iA,IatomA,IatomB,IatomC,IatomD
  integer :: startC,nPasses,startD,intprint,iPass,IatomBend
  real(realk) :: DcenterSpec(3),CcenterSpec(3),AcenterSpec(3),BcenterSpec(3),GABELM
  logical :: PermuteRHS
  real(realk),allocatable :: TmpArray1(:),TmpArray2(:)
  real(realk),allocatable :: LocalIntPass1(:),LocalIntPass2(:)
  integer,allocatable :: IatomAPass(:),IatomBPass(:)
  logical,allocatable :: DoInt(:,:)
  logical :: NOELEMENTSADDED
  integer :: iOrbQ,iOrbB,iOrbA,iOrbD,iOrbC,I4,I3,I2
  integer :: startA,startB,ndim,nOrbQ,MaxPasses
  integer :: TMParray1maxsizePass,TMParray2maxsizePass,nLocalIntPass
#ifdef VAR_OMP
  integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif

  !FIXME: Determine MaxPasses 
  MaxPasses = nAtomsA*nAtomsB

  TMParray1maxsizePass = TMParray1maxsize*MaxPasses
  TMParray2maxsizePass = TMParray2maxsize*MaxPasses
  
  ndim = nOrbA*nOrbB*nAtomsA*nAtomsB
  nOrbQ = nOrbC*nOrbD
  allocate(TmpArray1(TMParray1maxsize*MaxPasses))
  call mem_ichor_alloc(TmpArray1)
  allocate(TmpArray2(TMParray2maxsize*MaxPasses))     
  call mem_ichor_alloc(TmpArray2)
  allocate(IatomAPass(MaxPasses))
  call mem_ichor_alloc(IatomAPass)  
  allocate(IatomBPass(MaxPasses))
  call mem_ichor_alloc(IatomBPass)
    
  nLocalIntPass = nLocalint*MaxPasses
  allocate(LocalIntPass1(nLocalIntPass))
  CALL Mem_ichor_alloc(LocalIntPass1)
  allocate(LocalIntPass2(nLocalint*nAtomsA*nAtomsB))
  CALL Mem_ichor_alloc(LocalIntPass2)

  allocate(DoINT(nAtomsA,nAtomsB))
  CALL Mem_ichor_alloc(DoINT)
  !Link Procedure
  IF(UseGeneralCode)THEN
    !$OMP MASTER
     call DetermineSizeTmpArray34(nTUVQ,nCartOrbCompQ,nPrimQ,nTUVP,nCartOrbCompP,nPrimP,MaxPasses,&
          & AngmomA,AngmomB,AngmomC,AngmomD,AngmomA+AngmomB,AngmomC+AngmomD,TotalAngmom)
     allocate(TmpArray3(nTmpArray3))
     call mem_ichor_alloc(TmpArray3)
     allocate(TmpArray4(nTmpArray4))
     call mem_ichor_alloc(TmpArray4)
     CALL PreCalciChorSPHMAT(MAX(AngmomA,AngmomB,AngmomC,AngmomD))
     !$OMP END MASTER
     !$OMP BARRIER 
  ENDIF
  DO iatomD = 1,natomsD
   DcenterSpec(1) = Dcenter(1,IAtomD)
   DcenterSpec(2) = Dcenter(2,IAtomD)
   DcenterSpec(3) = Dcenter(3,IAtomD)
   startD = (IatomD-1)*nOrbD 
   DO iC = 1,nKetList(IatomD) !include IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
    IatomC = KetList(iC,IatomD)

    PermuteRHS = TriangularRHSAtomLoop.AND.IatomD.LT.IatomC
    DoINT = .FALSE.
    !making List of {AB} which contribute for CD: K(AC) = sum_{BD} (A,B|C,D) Dmat(B,D)
    LOOPB: DO iB = 1,nBraketList(IatomD)
     IatomB = BraKetList(iB,IatomD)
     NOELEMENTSADDED = .TRUE.
     LOOPA: DO iA = 1,nBraList(IatomB)
      IatomA = BraList(iA,IatomB)
      IF(ReducedDmatBD(IatomB,IatomD)*atomGAB(IatomA,IatomB)*atomGCD(IatomC,IatomD) .LE. THRESHOLD_CS) EXIT LOOPA
      DoINT(IatomA,IatomB) = .TRUE. 
      NOELEMENTSADDED = .FALSE.
     ENDDO LOOPA
     IF(NOELEMENTSADDED)EXIT LOOPB
    ENDDO LOOPB
    IF(SameRHSaos)THEN !permutational symmetry (A,B|D,C) = (A,B|C,D)
     !making List of {AB} which contribute for CD: K(AD) = sum_{BC} (A,B|D,C) Dmat(B,C)
     LOOPB2: DO iB = 1,nBraketList(IatomC) 
      IatomB = BraKetList(iB,IatomC)
      NOELEMENTSADDED = .TRUE.
      LOOPA2: DO iA = 1,nBraList(IatomB)
       IatomA = BraList(iA,IatomB)
       IF( ReducedDmatBC(IatomB,IatomC)*AtomGAB(IatomA,IatomB)*AtomGCD(IatomC,IatomD) .LE. THRESHOLD_CS) EXIT LOOPA2
       DoINT(IatomA,IatomB) = .TRUE. 
       NOELEMENTSADDED = .FALSE.
      ENDDO LOOPA2
      IF(NOELEMENTSADDED)EXIT LOOPB2
     ENDDO LOOPB2
    ENDIF

    !Make IatomAPass,IatomBPass,nPasses from  DoInt(IatomA,IatomB)
    iPass=0
    IatomBend = nAtomsB
    DO IatomA = 1,nAtomsA
       IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
       DO IatomB = 1,IatomBend
          IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
             IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
                IF(DoINT(IatomA,IatomB))THEN
                   iPass = iPass + 1
                   IatomAPass(iPass) = IatomA
                   IatomBPass(iPass) = IatomB
                ENDIF
             ENDIF
          ELSE
             IF(DoINT(IatomA,IatomB))THEN
                iPass = iPass + 1
                IatomAPass(iPass) = IatomA
                IatomBPass(iPass) = IatomB
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    nPasses = iPass

    CcenterSpec(1) = Ccenter(1,IAtomC)
    CcenterSpec(2) = Ccenter(2,IAtomC)
    CcenterSpec(3) = Ccenter(3,IAtomC)

    !FIXME SHOULD NOT BE NEEDED
    IF(nPasses.NE.nAtomsA*nAtomsB)THEN
       do I4 = 1,nLocalIntPass
          LocalIntPass1(I4) = 0.0E0_realk
       enddo
       do I4 = 1,nLocalint*nAtomsA*nAtomsB
          LocalIntPass2(I4) = 0.0E0_realk
       enddo
    ENDIF

    !output: Qcent,Qdistance12,QpreExpFac
    CALL Build_qcent_Qdistance12_QpreExpFac(nPrimC,nPrimD,nContC,nContD,&
         & expC,expD,CcenterSpec,DcenterSpec,ContractCoeffC,ContractCoeffD,&
         & Qsegmented,Qcent,Qdistance12,QpreExpFac,INTPRINT)
    !Unique for each iPassQ (iAtomC,iAtomD) iteration: qcent,qdistance12,qpreexpfac, Qiprim1(nPrimQ), output:
    !LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,MaxPasses)
    !IatomAPass,iatomBPass changes and 
    !     IF(iAtomC.EQ.1.AND.iAtomD.EQ.1)INTPRINT=1000
    call ICI_CPU_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
         & nPrimQ,nPrimP*nPrimQ,nPasses,MaxPasses,intprint,lupri,&
         & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
         & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
         & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
         & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
         & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
         & pcentPass,qcent,PpreexpfacPass,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
         & Qiprim1,Qiprim2,expA,expB,expC,expD,&
         & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
         & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,Qdistance12,PQorder,&
         & LocalIntPass1,nLocalIntPass,Acenter,Bcenter,CcenterSpec,DcenterSpec,&
         & nAtomsA,nAtomsB,Spherical,TmpArray1,TMParray1maxsizePass,TmpArray2,&
         & TMParray2maxsizePass,IatomAPass,iatomBPass,useSP)
    !output private LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
    !reorder (including LHS permute) to LocalIntPass(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
    !this can be done on the accelerator.
     IF(.TRUE.)THEN !nPrimLast
        call MainTriDistributetoLocalIntPass2CPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
             & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
             & MaxPasses,nPasses,TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1,&
             & LocalIntPass2,IatomAPass,iatomBPass,nContQ,nContP)
     ELSE
        call MainTriDistributetoLocalIntPass2GPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
             & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
             & MaxPasses,nPasses,TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1,&
             & LocalIntPass2,IatomAPass,iatomBPass,nContQ,nContP)
     ENDIF
    !TriangularLHSAtomLoop have ensured that LocalIntPass have full set of AtomA and AtomB
    !Contract LocalIntPass(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD) With D(B,D) and D()
    !K(AC) = sum_{BD} (A,B|C,D) Dmat(B,D) = sum_{BD} (A,B|C,D) Dmat(B,D) 
    !K(AD) = sum_{BC} (A,B|D,C) Dmat(B,C) = sum_{BC} (A,B|C,D) Dmat(B,C) if PermuteRHS
    !K(BC) = sum_{AD} (B,A|C,D) Dmat(A,D) = sum_{AD} (A,B|C,D) Dmat(A,D) if PermuteLHSTypes
    !K(BD) = sum_{AC} (B,A|D,C) Dmat(A,C) = sum_{AC} (A,B|C,D) Dmat(A,C) if PermuteLHSTypes and PermuteRHS
    startC = (IatomC-1)*nOrbC 
    call DistributeLink_KmatAC(nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,nDmat,DmatBD,&
         & LocalIntPass2,KmatAC,startC,startD,lupri)

!     IF(PermuteLHSTypes)THEN
!      IF(PermuteRHS)THEN

!    CALL LinkDistribution()


   ENDDO !IatomC
  ENDDO !iAtomD
  IF(UseGeneralCode)THEN
    !$OMP MASTER
    call mem_ichor_dealloc(TmpArray3)
    deallocate(TmpArray3)
    call mem_ichor_dealloc(TmpArray4)
    deallocate(TmpArray4)
    !FIXME MUCH LATER
    call FreeIchorSPHMAT()
    !$OMP END MASTER
    !$OMP BARRIER 
  ENDIF

  call mem_ichor_dealloc(TmpArray1)
  deallocate(TmpArray1)
  call mem_ichor_dealloc(TmpArray2)
  deallocate(TmpArray2)    
  call mem_ichor_dealloc(IatomBPass) 
  deallocate(IatomBPass) 
  call mem_ichor_dealloc(IatomAPass) 
  deallocate(IatomAPass) 
  CALL Mem_ichor_dealloc(LocalIntPass1)
  deallocate(LocalIntPass1)
  CALL Mem_ichor_dealloc(LocalIntPass2)
  deallocate(LocalIntPass2)
  CALL Mem_ichor_dealloc(DoINT)
  deallocate(DoINT)

end subroutine IchorTypeLinKLoop

subroutine DistributeLink_KmatAC(nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,nDmat,&
     & DmatBD,LocalIntPass2,KmatAC,startC,startD,lupri)
  implicit none
  integer,intent(in) :: nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,nDmat
  integer,intent(in) :: startC,startD,lupri
  real(realk),intent(in) :: DmatBD(ndimB,ndimD,nDmat)
  real(realk),intent(in) :: LocalIntPass2(ndimA,ndimB,nOrbC,nOrbD)
  real(realk),intent(inout) :: KmatAC(ndimA,ndimC,nDmat)
  !
  integer :: iOrbD,iOrbC,B,A,idmat
  real(realk) :: TMPD
  DO idmat = 1,nDmat
     DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
           DO B = 1,ndimB !nAtomsB*nOrbB
              TMPD = DmatBD(B,startD + iOrbD,idmat)
              DO A = 1,ndimA !nAtomsA*nOrbA
                 KmatAC(A,startC+iOrbC,idmat) = KmatAC(A,startC+iOrbC,idmat) + &
                      & LocalIntPass2(A,B,iOrbC,iOrbD)*TMPD
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
end subroutine DistributeLink_KmatAC

subroutine IchorTypeMOtransLoop(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
     & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
     & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
     & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
     & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
     & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
     & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
     & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
     & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
     & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
     & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
     & qcent,Ppreexpfac,Qpreexpfac,&
     & Qiprim1,Qiprim2,&
     & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
     & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
     & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
     & reducedExponents,integralPrefactor,&
     & PcentPass,Pdistance12Pass,PpreExpFacPass,&
     & Qdistance12,PQorder,&
     & BATCHGCD,BATCHGAB,Spherical,TMParray1maxsize,nLocalInt,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,OutputStorage,&
     & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri,&
     & nCMO1,nCMO2,nCMO3,nCMO4,CMO1A,CMO2B,CMO3C,CMO4D,&
     & CMO1B,CMO2A,CMO3D,CMO4C,PermuteRHSTypes)
  implicit none
  logical,intent(in) :: TriangularRHSAtomLoop,TriangularLHSAtomLoop,PermuteLHSTypes
  logical,intent(in) :: Qsegmented,Psegmented,PQorder,Spherical,CSscreen
  logical,intent(in) :: TriangularODAtomLoop,PermuteRHSTypes
  integer,intent(in) :: nTABFJW1,nTABFJW2,lupri,nCMO1,nCMO2,nCMO3,nCMO4
  integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
  integer,intent(in) :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV
  real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2),THRESHOLD_CS
  real(realk),intent(in) :: CMO1A(nOrbA*nAtomsA,nCMO1),CMO2B(nOrbB*nAtomsB,nCMO2)
  real(realk),intent(in) :: CMO3C(nOrbC*nAtomsC,nCMO3),CMO4D(nOrbD*nAtomsD,nCMO4)
  real(realk),intent(in) :: CMO2A(nOrbA*nAtomsA,nCMO2),CMO1B(nOrbB*nAtomsB,nCMO1)
  real(realk),intent(in) :: CMO4C(nOrbC*nAtomsC,nCMO4),CMO3D(nOrbD*nAtomsD,nCMO3)

  integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)
  !D
  integer,intent(in) :: nAtomsD,nPrimD,nContD,nOrbCompD,AngmomD,nBatchD,nOrbD
  integer,intent(in) :: iBatchIndexOfTypeD
  integer,intent(in) :: startOrbitalD(nAtomsD)
  real(realk),intent(in) :: Dcenter(3,nAtomsD),expD(nPrimD),ContractCoeffD(nPrimD,nContD)
  !C
  integer,intent(in) :: nAtomsC,nPrimC,nContC,nOrbCompC,AngmomC,nBatchC,nOrbC
  integer,intent(in) :: iBatchIndexOfTypeC
  integer,intent(in) :: startOrbitalC(nAtomsC)
  real(realk),intent(in) :: Ccenter(3,nAtomsC),expC(nPrimC),ContractCoeffC(nPrimC,nContC)
  logical,intent(in) :: noScreenCD2(nAtomsC,nAtomsD)
  logical,intent(in) :: noScreenAB(nAtomsA,nAtomsB)
  !Q
  integer,intent(in) :: nContQ,nPrimQ
  real(realk),intent(inout) :: Qcent(3,nPrimQ),Qdistance12(3),QpreExpFac(nPrimQ),expQ(nPrimQ)
  !P
  integer,intent(in) :: nContP,nPrimP
  real(realk),intent(inout) :: PpreExpFac(nPrimP)
  real(realk),intent(in) :: PcentPass(3,nPrimP,natomsA*natomsB)
  real(realk),intent(in) :: Pdistance12Pass(3,natomsA*natomsB)
  real(realk),intent(in) :: PpreExpFacPass(nPrimP,natomsA*natomsB),expP(nPrimP)
  !
  real(realk),intent(in) :: BATCHGCD(nBatchC,nBatchD)
  !A & B
  integer,intent(in) :: iBatchIndexOfTypeA,iBatchIndexOfTypeB,AngmomB,AngmomA,nOrbA,nOrbB
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbCompA,nOrbCompB,nBatchA,nBatchB,nContA,nContB
  integer,intent(in) :: nPrimA,nPrimB
  integer,intent(in) :: startOrbitalA(nAtomsA)
  integer,intent(in) :: startOrbitalB(nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA),expA(nPrimA),ContractCoeffA(nPrimA,nContA)
  real(realk),intent(in) :: Bcenter(3,nAtomsB),expB(nPrimB),ContractCoeffB(nPrimB,nContB)
  real(realk),intent(in) ::  BATCHGAB(nBatchA*nBatchB)
  !collected
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,TotalAngmom
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  integer,intent(in) :: TMParray1maxsize,TMParray2maxsize,nLocalInt
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: integralPrefactor(nPrimQ,nPrimP)
  !local variables
  integer :: iBatchD,IatomD,IatomC,IatomB,nPasses,intprint
  real(realk) :: DcenterSpec(3),CcenterSpec(3),AcenterSpec(3),BcenterSpec(3),GABELM
  logical :: PermuteRHS
  real(realk),allocatable :: TmpArray1(:),TmpArray2(:)
  real(realk),allocatable :: LocalIntPass1(:),LocalIntPass2(:)
  real(realk),allocatable :: OutputA(:,:,:)
  real(realk),allocatable :: OutputC(:,:,:),OutputCD(:,:,:)
  integer,allocatable :: IatomAPass(:),IatomBPass(:)
  integer :: iOrbQ,iOrbB,iOrbA,iOrbD,iOrbC,I4,I3,I2
  integer :: ndim,nOrbQ,MaxPasses
  integer :: TMParray1maxsizePass,TMParray2maxsizePass,nLocalIntPass
  integer :: nDimA,nDimB,nDimC,nDimD
#ifdef VAR_OMP
  integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
  nLocalIntPass = nOrbA*nAtomsA*nOrbB*nAtomsB*nOrbC*nOrbD
  nDimA = nOrbA*nAtomsA
  nDimB = nOrbB*nAtomsB
  nDimC = nOrbC*nAtomsC
  nDimD = nOrbD*nAtomsD
  call DetermineMaxPasses(nAtomsD,iBatchIndexOfTypeD,nAtomsC,nAtomsA,nAtomsB,&
       & iBatchIndexOfTypeC,iBatchIndexOfTypeA,nBatchB,nBatchA,iBatchIndexOfTypeB,&
       & TriangularRHSAtomLoop,CSscreen,TriangularLHSAtomLoop,TriangularODAtomLoop,&
       & noScreenCD2,BATCHGAB,THRESHOLD_CS,noScreenAB,BATCHGCD,nBatchC,nBatchD,MaxPasses)
  TMParray1maxsizePass = TMParray1maxsize*MaxPasses
  TMParray2maxsizePass = TMParray2maxsize*MaxPasses

  ndim = nOrbA*nOrbB*nAtomsA*nAtomsB
  nOrbQ = nOrbC*nOrbD
  allocate(TmpArray1(TMParray1maxsize*MaxPasses))
  call mem_ichor_alloc(TmpArray1)
  allocate(TmpArray2(TMParray2maxsize*MaxPasses))     
  call mem_ichor_alloc(TmpArray2)
  allocate(IatomAPass(MaxPasses))
  call mem_ichor_alloc(IatomAPass)  
  allocate(IatomBPass(MaxPasses))
  call mem_ichor_alloc(IatomBPass)

  nLocalIntPass = nLocalint*MaxPasses
  allocate(LocalIntPass1(nLocalIntPass))
  CALL Mem_ichor_alloc(LocalIntPass1)
  allocate(LocalIntPass2(nLocalint*nAtomsA*nAtomsB))
  CALL Mem_ichor_alloc(LocalIntPass2)

  allocate(OutputA(nCMO1,MAX(ndimA,ndimB),nOrbC*nOrbD))
  CALL Mem_ichor_alloc(OutputA)
  allocate(OutputCD(nCMO1,nCMO2,ndimC*nDimD))
  CALL Mem_ichor_alloc(OutputCD)
  !is this necessary
  call ichorzero2(OutputCD,nCMO1*nCMO2,ndimC*nDimD)

  IF(UseGeneralCode)THEN
     call DetermineSizeTmpArray34(nTUVQ,nCartOrbCompQ,nPrimQ,nTUVP,nCartOrbCompP,nPrimP,MaxPasses,&
          & AngmomA,AngmomB,AngmomC,AngmomD,AngmomA+AngmomB,AngmomC+AngmomD,TotalAngmom)
     allocate(TmpArray3(nTmpArray3))
     call mem_ichor_alloc(TmpArray3)
     allocate(TmpArray4(nTmpArray4))
     call mem_ichor_alloc(TmpArray4)
     CALL PreCalciChorSPHMAT(MAX(AngmomA,AngmomB,AngmomC,AngmomD))
  ENDIF
!$OMP PARALLEL DEFAULT(none) &
!$OMP PRIVATE(iAtomD,iAtomC,GABELM,iBatchD,DcenterSpec,PermuteRHS,&
!$OMP         CcenterSpec,iOrbQ,I3,I4,iOrbD,iOrbC,iAtomB) &
!$OMP SHARED(nAtomsD,iBatchIndexOfTypeD,Dcenter,nAtomsC,&
!$OMP        TriangularRHSAtomLoop,noScreenCD2,Ccenter,&
!$OMP        nAtomsA,nAtomsB,BATCHGCD,iBatchIndexOfTypeC,CSscreen,&
!$OMP        nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
!$OMP        BATCHGAB,THRESHOLD_CS,IatomAPass,IatomBPass,&
!$OMP        MaxPasses,TriangularLHSAtomLoop,TriangularODAtomLoop,nOrbB,&
!$OMP        Qsegmented,nPasses,noScreenAB,nLocalInt,TotalAngmom,nOrbA,&
!$OMP        nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,intprint,lupri,&
!$OMP        nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,TABFJW,&
!$OMP        ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
!$OMP        pcentPass,qcent,PpreexpfacPass,Qpreexpfac,nTABFJW1,nTABFJW2,&
!$OMP        Qiprim1,Qiprim2,expA,expB,expC,expD,Psegmented,reducedExponents,&
!$OMP        integralPrefactor,AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,&
!$OMP        Qdistance12,PQorder,LocalIntPass1,LocalIntPass2,nLocalIntPass,&
!$OMP        Spherical,TmpArray1,TMParray1maxsizePass,TmpArray2,Bcenter,nOrbQ,&
!$OMP        TMParray2maxsizePass,Acenter,TmpArray3,TmpArray4,&
!$OMP        nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,PermuteLHSTypes,nOrbD,nOrbC,&
!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputStorage,&
!$OMP        nDimA,nDimB,nDimC,nDimD,CMO1A,CMO1B,CMO2A,CMO2B,nCMO1,nCMO2,&
!$OMP        OutputA,OutputCD,PermuteRHSTypes,nTmpArray3,nTmpArray4,&
!$OMP        nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
!$OMP        nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV)
  DO IatomD = 1,nAtomsD
   GABELM = 0.0E0_realk 
   iBatchD = iBatchIndexOfTypeD + IatomD
   DcenterSpec(1) = Dcenter(1,IAtomD)
   DcenterSpec(2) = Dcenter(2,IAtomD)
   DcenterSpec(3) = Dcenter(3,IAtomD)
   DO IatomC = 1,nAtomsC
    IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
    PermuteRHS = TriangularRHSAtomLoop.AND.IatomD.LT.IatomC
    IF(noScreenCD2(IatomC,IatomD))THEN
     CcenterSpec(1) = Ccenter(1,IAtomC)
     CcenterSpec(2) = Ccenter(2,IAtomC)
     CcenterSpec(3) = Ccenter(3,IAtomC)
     IF(CSscreen)GABELM = BATCHGCD(iBatchIndexOfTypeC+IatomC,iBatchD)
     !output: IatomAPass,IatomBPass,nPasses
!$OMP SINGLE
     nPasses = nAtomsA*nAtomsB
     CALL BUILD_noScreen2(CSscreen,nAtomsA,nAtomsB,&
          & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
          & BATCHGAB,THRESHOLD_CS,GABELM,nPasses,IatomAPass,IatomBPass,MaxPasses,&
          & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB) 
!$OMP END SINGLE
!$OMP BARRIER
     IF(nPasses.EQ.0)CYCLE
     IF(nPasses.NE.nAtomsA*nAtomsB)THEN
!!$OMP DO PRIVATE(I4)
!        do I4 = 1,nLocalIntPass
!           LocalIntPass1(I4) = 0.0E0_realk
!        enddo
!!$OMP END DO NOWAIT
!$OMP DO PRIVATE(I4)
        do I4 = 1,nLocalint*nAtomsA*nAtomsB
           LocalIntPass2(I4) = 0.0E0_realk
        enddo
!$OMP END DO
     ENDIF
!$OMP SINGLE
     !output: Qcent,Qdistance12,QpreExpFac
     CALL Build_qcent_Qdistance12_QpreExpFac(nPrimC,nPrimD,nContC,nContD,&
          & expC,expD,CcenterSpec,DcenterSpec,ContractCoeffC,ContractCoeffD,&
          & Qsegmented,Qcent,Qdistance12,QpreExpFac,INTPRINT)
!$OMP END SINGLE
!$OMP BARRIER
     !Unique for each iPassQ (iAtomC,iAtomD) iteration: qcent,qdistance12,qpreexpfac, Qiprim1(nPrimQ), output:
     !LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,MaxPasses)
     !IatomAPass,iatomBPass changes and 
!     IF(iAtomC.EQ.1.AND.iAtomD.EQ.1)INTPRINT=1000

     call ICI_CPU_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
          & nPrimQ,nPrimP*nPrimQ,nPasses,MaxPasses,intprint,lupri,&
          & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
          & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
          & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
          & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
          & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
          & pcentPass,qcent,PpreexpfacPass,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
          & Qiprim1,Qiprim2,expA,expB,expC,expD,&
          & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
          & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,Qdistance12,PQorder,&
          & LocalIntPass1,nLocalIntPass,Acenter,Bcenter,CcenterSpec,DcenterSpec,&
          & nAtomsA,nAtomsB,Spherical,TmpArray1,TMParray1maxsizePass,TmpArray2,&
          & TMParray2maxsizePass,IatomAPass,iatomBPass,useSP)
     !output private LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
     !reorder (including LHS permute) to LocalIntPass(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
     !this can be done on the accelerator
     call MainTriDistributetoLocalIntPass2CPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
          & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
          & MaxPasses,nPasses,TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1,&
          & LocalIntPass2,IatomAPass,iatomBPass,nContQ,nContP)

     !OutputA(nCMO1,ndimB,nOrbC,nOrbD)=LocalIntPass2(ndimA,ndimB,nOrb,nOrbD)*CMO1A(ndimA,nCMO1)
     IF(nDimA.NE.1)THEN
        CALL DistributeLink_MOtransformA(nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,&
             & CMO1A,LocalIntPass2,OutputA,lupri,nCMO1)        
     ELSE
        CALL DistributeLink_MOtransformA1(nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,&
             & CMO1A,LocalIntPass2,OutputA,lupri,nCMO1)
     ENDIF
     !OutputCD1(nCMO1,nCMO2,ndimC,ndimD) = OutputCD1(nCMO1,nCMO2,ndimC,ndimD) 
     !                                   + OutputA(nCMO1,ndimB,nOrbC,nOrbD)*CMO2B(ndimB,nCMO2)
     CALL DistributeLink_MOtransformB(nOrbD,nOrbC,ndimB,ndimC,ndimD,&
          & CMO2B,OutputA,OutputCD,lupri,nCMO1,nCMO2,iAtomC,iAtomD,PermuteRHS)

     IF(PermuteLHSTypes)THEN
        !OutputA(nCMO1,ndimA,nOrbC,nOrbD) = LocalIntPass2(ndimA,ndimB,nOrb,nOrbD)*CMO1B(ndimB,nCMO1)
        CALL DistributeLink_MOtransformA2(nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,&
             & CMO1B,LocalIntPass2,OutputA,lupri,nCMO1)        
        !OutputCD(nCMO1,nCMO2,ndimC,ndimD) = OutputCD(nCMO1,nCMO2,ndimC,ndimD) 
        !                                  + OutputA(nCMO1,ndimA,nOrbC,nOrbD)*CMO2A(ndimA,nCMO2)
        CALL DistributeLink_MOtransformB2(nOrbD,nOrbC,ndimA,ndimC,ndimD,&
             & CMO2A,OutputA,OutputCD,lupri,nCMO1,nCMO2,iAtomC,iAtomD,PermuteRHS)
     ENDIF
    ENDIF !noscreenCD2
   ENDDO !IatomC
  ENDDO !iAtomD
!$OMP END PARALLEL
  IF(UseGeneralCode)THEN
    call mem_ichor_dealloc(TmpArray3)
    deallocate(TmpArray3)
    call mem_ichor_dealloc(TmpArray4)
    deallocate(TmpArray4)
    call FreeIchorSPHMAT()
  ENDIF

  !SYMMETRY 
  !(ndimA,ndimB,ndimC,ndimD) = (ndimC,ndimD,ndimA,ndimB) =
  !(ndimB,ndimA,ndimC,ndimD) = (ndimC,ndimD,ndimB,ndimA) = 
  !(ndimA,ndimB,ndimD,ndimC) = (ndimD,ndimC,ndimA,ndimB) = 
  !(ndimB,ndimA,ndimD,ndimC) = (ndimD,ndimC,ndimB,ndimA) =

  !LHS                        1  LocalIntPass2           2 different OutputA          1 different OutputCD        1 different OutputC                            
  !LHS,RHS                    1  LocalIntPass2           2 different OutputA          1 different OutputCD        2 different OutputC                            
  !LHS and RHS                1  LocalIntPass2           2 different OutputA          1 different OutputCD        2 different OutputC                            
  !                           1  LocalIntPass2           4 different OutputA          4 different OutputCD        4 different OutputC                            
  !CMO1A*CMO2B*CMO3C*CMO4D*(ndimA,ndimB,nOrbC,nOrbD) - (nCMO1,ndimB,nOrbC,nOrbD) - (nCMO1,nCMO2,ndimC,ndimD) - (nCMO1,nCMO2,nCMO3,ndimD)    FULL
  !CMO1B*CMO2A*CMO3C*CMO4D*(ndimB,ndimA,nOrbC,nOrbD) - (nCMO1,ndimA,nOrbC,nOrbD) - (nCMO1,nCMO2,ndimC,ndimD) - (nCMO1,nCMO2,nCMO3,ndimD)    SameLHS
  !CMO1A*CMO2B*CMO3D*CMO4C*(ndimA,ndimB,nOrbD,nOrbC) - (nCMO1,ndimB,nOrbD,nOrbC) - (nCMO1,nCMO2,ndimD,ndimC) - (nCMO1,nCMO2,nCMO3,ndimC)    SameRHS
  !CMO1B*CMO2A*CMO3D*CMO4C*(ndimB,ndimA,nOrbD,nOrbC) - (nCMO1,ndimA,nOrbD,nOrbC) - (nCMO1,nCMO2,ndimD,ndimC) - (nCMO1,nCMO2,nCMO3,ndimC)    SameLHS,SameRHS


  !CMO1C*CMO2D*CMO3A*CMO4B*(nOrbC,nOrbD,ndimA,ndimB) - (nCMO1,nOrbD,ndimA,ndimB) - (nCMO1,nCMO2,ndimA,ndimB) - (nCMO1,nCMO2,nCMO3,ndimB)    SameODs
  !CMO1D*CMO2C*CMO3A*CMO4B*(nOrbD,nOrbC,ndimA,ndimB) - (nCMO1,nOrbC,ndimA,ndimB) - (nCMO1,nCMO2,ndimA,ndimB) - (nCMO1,nCMO2,nCMO3,ndimB)    SameODs,SameRHS
  !CMO1C*CMO2D*CMO3B*CMO4A*(nOrbC,nOrbD,ndimB,ndimA) - (nCMO1,nOrbD,ndimB,ndimA) - (nCMO1,nCMO2,ndimB,ndimA) - (nCMO1,nCMO2,nCMO3,ndimA)    SameODs,SameLHS
  !CMO1D*CMO2C*CMO3B*CMO4A*(nOrbD,nOrbC,ndimB,ndimA) - (nCMO1,nOrbC,ndimB,ndimA) - (nCMO1,nCMO2,ndimB,ndimA) - (nCMO1,nCMO2,nCMO3,ndimA)    SameODs,SameLHS,SameRHS

  CALL Mem_ichor_dealloc(OutputA)
  deallocate(OutputA)

  call mem_ichor_dealloc(TmpArray1)
  deallocate(TmpArray1)
  call mem_ichor_dealloc(TmpArray2)
  deallocate(TmpArray2)    
  call mem_ichor_dealloc(IatomBPass) 
  deallocate(IatomBPass) 
  call mem_ichor_dealloc(IatomAPass) 
  deallocate(IatomAPass) 
  CALL Mem_ichor_dealloc(LocalIntPass1)
  deallocate(LocalIntPass1)
  CALL Mem_ichor_dealloc(LocalIntPass2)
  deallocate(LocalIntPass2)

  IF(PermuteRHSTypes)THEN
     allocate(OutputC(nCMO1,nCMO2,nCMO3*MAX(nDimD,nDimC)))
     CALL Mem_ichor_alloc(OutputC)

     !OutputC(nCMO1,nCMO2,nCMO3,ndimD)=OutputCD(nCMO1,nCMO2,ndimC,ndimD)*CMO3C(ndimC,nCMO3)
     CALL DistributeLink_MOtransformC(ndimC,ndimD,&
          & CMO3C,OutputCD,OutputC,lupri,nCMO1,nCMO2,nCMO3)
     
     !OutputStorage(nCMO1*nCMO2*nCMO3,nCMO4) = OutputStorage(nCMO1*nCMO2*nCMO3,nCMO4)
     ! + OutputC(nCMO1*nCMO2*nCMO3,nDimD) * CMO4D(nDimD,nCMO4)
     CALL DistributeLink_MOtransformD(ndimD,CMO4D,OutputC,OutputStorage,lupri,nCMO1,nCMO2,nCMO3,nCMO4)
     !  OR
     !call dgemm('N','N',nCMO1*nCMO2*nCMO3,nCMO4,nDimD,1.0E0_realk,OutputC,nCMO1*nCMO2*nCMO3,&
     !     & CMO4D,nDimD,1.0E0_realk,OutputStorage,nCMO1*nCMO2*nCMO3)

     !OutputC(nCMO1,nCMO2,nCMO3,ndimC)=OutputCD(nCMO1,nCMO2,ndimC,ndimD)*CMO3D(ndimD,nCMO3)
     CALL DistributeLink_MOtransformC2(ndimC,ndimD,&
          & CMO3D,OutputCD,OutputC,lupri,nCMO1,nCMO2,nCMO3)
     
     CALL Mem_ichor_dealloc(OutputCD)
     deallocate(OutputCD)

     !OutputStorage(nCMO1*nCMO2*nCMO3,nCMO4) = OutputStorage(nCMO1*nCMO2*nCMO3,nCMO4)
     ! + OutputC(nCMO1*nCMO2*nCMO3,nDimC) * CMO4C(nDimC,nCMO4)
     CALL DistributeLink_MOtransformD(ndimC,CMO4C,OutputC,OutputStorage,lupri,nCMO1,nCMO2,nCMO3,nCMO4)
     !  OR 
     !call dgemm('N','N',nCMO1*nCMO2*nCMO3,nCMO4,nDimC,1.0E0_realk,OutputC,nCMO1*nCMO2*nCMO3,&
     !     & CMO4C,nDimC,1.0E0_realk,OutputStorage,nCMO1*nCMO2*nCMO3)
  ELSE
     allocate(OutputC(nCMO1,nCMO2,nCMO3*nDimD))
     CALL Mem_ichor_alloc(OutputC)

     !OutputC(nCMO1,nCMO2,nCMO3,ndimD)=OutputCD(nCMO1,nCMO2,ndimC,ndimD)*CMO3C(ndimC,nVirt)
     CALL DistributeLink_MOtransformC(ndimC,ndimD,&
          & CMO3C,OutputCD,OutputC,lupri,nCMO1,nCMO2,nCMO3)
     
     CALL Mem_ichor_dealloc(OutputCD)
     deallocate(OutputCD)

     !OutputStorage(nCMO1*nCMO2*nCMO3,nCMO4) = OutputStorage(nCMO1*nCMO2*nCMO3,nCMO4)
     ! + OutputC(nCMO1*nCMO2*nCMO3,nDimD) * CMO4(nDimD,nCMO4)
     CALL DistributeLink_MOtransformD(ndimD,CMO4D,OutputC,OutputStorage,lupri,nCMO1,nCMO2,nCMO3,nCMO4)
     !  OR 
     !call dgemm('N','N',nCMO1*nCMO2*nCMO3,nCMO4,nDimD,1.0E0_realk,OutputC,nCMO1*nCMO2*nCMO3,&
     !     & CMO4D,nDimD,1.0E0_realk,OutputStorage,nCMO1*nCMO2*nCMO3)
  ENDIF

  CALL Mem_ichor_dealloc(OutputC)
  deallocate(OutputC)

end subroutine IchorTypeMOtransLoop

!scale ndimA*ndimB*nOrbD*nOrbC*nVirt
!IF room on accelerator this can be done on accelerator
subroutine DistributeLink_MOtransformA2(nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,&
     & CMO1B,LocalIntPass2,Output,lupri,nCMO1)
  implicit none
  integer,intent(in) :: nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD
  integer,intent(in) :: lupri,nCMO1
  real(realk),intent(in) :: CMO1B(ndimB,nCMO1)
  real(realk),intent(in) :: LocalIntPass2(ndimA,ndimB,nOrbC*nOrbD)
  real(realk),intent(inout) :: Output(nCMO1,ndimA,nOrbC*nOrbD)
  !
  integer :: iOrbCD,B,A,aVirt
  real(realk) :: TMP
!$OMP DO COLLAPSE(2) PRIVATE(iOrbCD,B,A,aVirt,TMP)
  DO iOrbCD = 1,nOrbD*nOrbC
     DO A = 1,ndimA !nAtomsB*nOrbB
        DO aVirt = 1,nCMO1
           TMP = 0.0E0_realk
           DO B = 1,ndimB !nAtomsA*nOrbA
              TMP = TMP + CMO1B(B,aVirt)*LocalIntPass2(A,B,iOrbCD)
           ENDDO
           Output(aVirt,A,iOrbCD) = TMP
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO
end subroutine DistributeLink_MOtransformA2

!scale ndimB*nOrbD*nOrbC*nOcc*nVirt
!IF room on accelerator this can be done on accelerator
subroutine DistributeLink_MOtransformB2(nOrbD,nOrbC,ndimA,ndimC,ndimD,&
     & CMO2A,Input,Output,lupri,nCMO1,nCMO2,iAtomC,iAtomD,PermuteRHS)
  implicit none
  integer,intent(in) :: nOrbD,nOrbC,ndimA,ndimC,ndimD
  integer,intent(in) :: lupri,nCMO1,nCMO2,iAtomC,iAtomD
  real(realk),intent(in) :: CMO2A(ndimA,nCMO2)
  real(realk),intent(in) :: Input(nCMO1,ndimA,nOrbC,nOrbD)
  real(realk),intent(inout) :: Output(nCMO1,nCMO2,nDimC,nDimD)
  logical,intent(in) :: PermuteRHS
  !
  integer :: iOrbC,iOrbD,iOcc,aVirt,A,C,D
  real(realk) :: TMP
  IF(PermuteRHS)THEN
     !$OMP DO COLLAPSE(3) PRIVATE(iOrbC,iOrbD,iOcc,aVirt,C,D,A,TMP)
     DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
           DO iOcc = 1,nCMO2
              D = iOrbD+(iAtomD-1)*nOrbD
              C = iOrbC+(iAtomC-1)*nOrbC
              !           already obtained a contribution DistributeLink_MOtransformB
              !           DO aVirt = 1,nCMO1
              !              Output(aVirt,iOcc,C,D) = 0.0E0_realk
              !           ENDDO
              DO A = 1,ndimA !nAtomsA*nOrbA
                 TMP = CMO2A(A,iOcc)
                 DO aVirt = 1,nCMO1
                    Output(aVirt,iOcc,C,D) = Output(aVirt,iOcc,C,D) + TMP*Input(aVirt,A,iOrbC,iOrbD)
                    Output(aVirt,iOcc,D,C) = Output(aVirt,iOcc,D,C) + TMP*Input(aVirt,A,iOrbC,iOrbD)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END DO
  ELSE
     !$OMP DO COLLAPSE(3) PRIVATE(iOrbC,iOrbD,iOcc,aVirt,C,D,A,TMP)
     DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
           DO iOcc = 1,nCMO2
              D = iOrbD+(iAtomD-1)*nOrbD
              C = iOrbC+(iAtomC-1)*nOrbC
              !           already obtained a contribution DistributeLink_MOtransformB
              !           DO aVirt = 1,nCMO1
              !              Output(aVirt,iOcc,C,D) = 0.0E0_realk
              !           ENDDO
              DO A = 1,ndimA !nAtomsA*nOrbA
                 TMP = CMO2A(A,iOcc)
                 DO aVirt = 1,nCMO1
                    Output(aVirt,iOcc,C,D) = Output(aVirt,iOcc,C,D) + TMP*Input(aVirt,A,iOrbC,iOrbD)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END DO
  ENDIF
end subroutine DistributeLink_MOtransformB2

!scale ndimA*ndimB*nOrbD*nOrbC*nVirt
!IF room on accelerator this can be done on accelerator
subroutine DistributeLink_MOtransformA(nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,&
     & CMO1,LocalIntPass2,Output,lupri,nCMO1)
  implicit none
  integer,intent(in) :: nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD
  integer,intent(in) :: lupri,nCMO1
  real(realk),intent(in) :: CMO1(ndimA,nCMO1)
  real(realk),intent(in) :: LocalIntPass2(ndimA,ndimB,nOrbC*nOrbD)
  real(realk),intent(inout) :: Output(nCMO1,ndimB,nOrbC*nOrbD)
  !
  integer :: iOrbCD,B,A,aVirt
  real(realk) :: TMP
!$OMP DO COLLAPSE(2) PRIVATE(iOrbCD,B,A,aVirt,TMP)
  DO iOrbCD = 1,nOrbD*nOrbC
     DO B = 1,ndimB !nAtomsB*nOrbB
        DO aVirt = 1,nCMO1
           TMP = 0.0E0_realk
           DO A = 1,ndimA !nAtomsA*nOrbA
              TMP = TMP + CMO1(A,aVirt)*LocalIntPass2(A,B,iOrbCD)
           ENDDO
           Output(aVirt,B,iOrbCD) = TMP
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO
end subroutine DistributeLink_MOtransformA

subroutine DistributeLink_MOtransformA1(nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD,&
     & CMO1,LocalIntPass2,Output,lupri,nCMO1)
  implicit none
  integer,intent(in) :: nOrbD,nOrbC,ndimA,ndimB,ndimC,ndimD
  integer,intent(in) :: lupri,nCMO1
  real(realk),intent(in) :: CMO1(nCMO1)
  real(realk),intent(in) :: LocalIntPass2(ndimB,nOrbC*nOrbD)
  real(realk),intent(inout) :: Output(nCMO1,ndimB,nOrbC*nOrbD)
  !
  integer :: iOrbCD,B,A,aVirt
  real(realk) :: TMP
!$OMP DO COLLAPSE(2) PRIVATE(iOrbCD,B,aVirt,TMP)
  DO iOrbCD = 1,nOrbD*nOrbC
     DO B = 1,ndimB !nAtomsB*nOrbB
        TMP = LocalIntPass2(B,iOrbCD)
        DO aVirt = 1,nCMO1
           Output(aVirt,B,iOrbCD) = CMO1(aVirt)*TMP
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO
end subroutine DistributeLink_MOtransformA1

!scale ndimB*nOrbD*nOrbC*nOcc*nVirt
!IF room on accelerator this can be done on accelerator
subroutine DistributeLink_MOtransformB(nOrbD,nOrbC,ndimB,ndimC,ndimD,&
     & CMO2,Input,Output,lupri,nCMO1,nCMO2,iAtomC,iAtomD,PermuteRHS)
  implicit none
  integer,intent(in) :: nOrbD,nOrbC,ndimB,ndimC,ndimD
  integer,intent(in) :: lupri,nCMO1,nCMO2,iAtomC,iAtomD
  real(realk),intent(in) :: CMO2(ndimB,nCMO2)
  real(realk),intent(in) :: Input(nCMO1,ndimB,nOrbC,nOrbD)
  real(realk),intent(inout) :: Output(nCMO1,nCMO2,nDimC,nDimD)
  logical,intent(in) :: PermuteRHS
  !
  integer :: iOrbC,iOrbD,iOcc,aVirt,B,C,D
  real(realk) :: TMP
  IF(PermuteRHS)THEN
!$OMP DO COLLAPSE(3) PRIVATE(iOrbC,iOrbD,iOcc,aVirt,C,D,B,TMP)
     DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
           DO iOcc = 1,nCMO2
              D = iOrbD+(iAtomD-1)*nOrbD
              C = iOrbC+(iAtomC-1)*nOrbC
              DO aVirt = 1,nCMO1
                 Output(aVirt,iOcc,C,D) = 0.0E0_realk
                 Output(aVirt,iOcc,D,C) = 0.0E0_realk
              ENDDO
              DO B = 1,ndimB !nAtomsB*nOrbB
                 TMP = CMO2(B,iOcc)
                 DO aVirt = 1,nCMO1
                    Output(aVirt,iOcc,C,D) = Output(aVirt,iOcc,C,D) + TMP*Input(aVirt,B,iOrbC,iOrbD)
                    Output(aVirt,iOcc,D,C) = Output(aVirt,iOcc,D,C) + TMP*Input(aVirt,B,iOrbC,iOrbD)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!$OMP END DO
  ELSE
!$OMP DO COLLAPSE(3) PRIVATE(iOrbC,iOrbD,iOcc,aVirt,C,D,B,TMP)
     DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
           DO iOcc = 1,nCMO2
              D = iOrbD+(iAtomD-1)*nOrbD
              C = iOrbC+(iAtomC-1)*nOrbC
              DO aVirt = 1,nCMO1
                 Output(aVirt,iOcc,C,D) = 0.0E0_realk
              ENDDO
              DO B = 1,ndimB !nAtomsB*nOrbB
                 TMP = CMO2(B,iOcc)
                 DO aVirt = 1,nCMO1
                    Output(aVirt,iOcc,C,D) = Output(aVirt,iOcc,C,D) + TMP*Input(aVirt,B,iOrbC,iOrbD)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!$OMP END DO
  ENDIF
end subroutine DistributeLink_MOtransformB

!scale nOrbD*nOrbC*nOcc*nVirt
subroutine CollectMOtransformCD(nOrbD,nOrbC,ndimC,ndimD,&
     & Input,Output,iAtomC,iAtomD,lupri,nCMO1,nCMO2)
  implicit none
  integer,intent(in) :: nOrbD,nOrbC,ndimC,ndimD,iAtomC,iAtomD
  integer,intent(in) :: lupri,nCMO1,nCMO2
  real(realk),intent(in) :: Input(nCMO1*nCMO2,nOrbC,nOrbD)
  real(realk),intent(inout) :: Output(nCMO1*nCMO2,ndimC,nDimD)
  !
  integer :: iOrbD,iOrbC,AI,C,D
!$OMP DO COLLAPSE(2) PRIVATE(iOrbD,iOrbC,AI)
  DO iOrbD = 1,nOrbD
     DO iOrbC = 1,nOrbC
        D = iOrbD+(iAtomD-1)*nOrbD
        C = iOrbC+(iAtomC-1)*nOrbC
        DO AI = 1,nCMO1*nCMO2
           Output(AI,C,D) = Input(AI,iOrbC,iOrbD)
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO
end subroutine CollectMOtransformCD

!scale nDimD*nDimC*nOcc*nVirt*nVirt
subroutine DistributeLink_MOtransformC(ndimC,ndimD,&
     & CMO3,Input,Output,lupri,nCMO1,nCMO2,nCMO3)
  implicit none
  integer,intent(in) :: ndimC,ndimD,lupri,nCMO1,nCMO2,nCMO3
  real(realk),intent(in) :: CMO3(ndimC,nCMO3)
  real(realk),intent(in) :: Input(nCMO1*nCMO2,ndimC,ndimD)
  real(realk),intent(inout) :: Output(nCMO1*nCMO2,nCMO3,ndimD)
  !
  integer :: D,bVirt,AI,C,nAI,nAI2,AIb
  real(realk) :: TMP
  logical :: MODAI
nAI = nCMO1*nCMO2
nAI2 = (nAI/64)*64
MODAI = (MOD(nAI,64).GT.0)
IF(nAI2.GT.0)THEN
   !$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) &   
   !$OMP PRIVATE(D,bVirt,AI,C,AIb,TMP) &
   !$OMP SHARED(ndimC,ndimD,nCMO1,nCMO2,nCMO3,Input,Output,nAI,nAI2,CMO3) 
   DO D = 1,ndimD
      DO AI = 1,nAI2,64
         DO bVirt = 1,nCMO3
            TMP = CMO3(1,bVirt)
            DO AIb = 0,63
               Output(AI+AIb,bVirt,D) = 0.0E0_realk
            ENDDO
            DO C = 1,ndimC
               TMP = CMO3(C,bVirt)
               DO AIb = 0,63
                  Output(AI+AIb,bVirt,D) = Output(AI+AIb,bVirt,D) + TMP*Input(AI+AIb,C,D)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
IF(MODAI)THEN
   !$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) &
   !$OMP PRIVATE(D,bVirt,C,TMP,AIb) &
   !$OMP SHARED(ndimC,ndimD,nCMO1,nCMO2,nCMO3,Input,Output,nAI,nAI2,CMO3) 
   DO D = 1,ndimD
      DO bVirt = 1,nCMO3
         TMP = CMO3(1,bVirt)
         DO AIb = nAI2+1,nAI
            Output(AIb,bVirt,D) = 0.0E0_realk
         ENDDO
         DO C = 1,ndimC
            TMP = CMO3(C,bVirt)
            DO AIb = nAI2+1,nAI
               Output(AIb,bVirt,D) = Output(AIb,bVirt,D) + TMP*Input(AIb,C,D)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
end subroutine DistributeLink_MOtransformC

subroutine DistributeLink_MOtransformC2(ndimC,ndimD,&
     & CMO3D,Input,Output,lupri,nCMO1,nCMO2,nCMO3)
  implicit none
  integer,intent(in) :: ndimC,ndimD,lupri,nCMO1,nCMO2,nCMO3
  real(realk),intent(in) :: CMO3D(ndimD,nCMO3)
  real(realk),intent(in) :: Input(nCMO1*nCMO2,ndimC,ndimD)
  real(realk),intent(inout) :: Output(nCMO1*nCMO2,nCMO3,ndimC)
  !
  integer :: D,bVirt,AI,C,nAI,nAI2,AIb
  real(realk) :: TMP
  logical :: MODAI
nAI = nCMO1*nCMO2
nAI2 = (nAI/64)*64
MODAI = (MOD(nAI,64).GT.0)
IF(nAI2.GT.0)THEN
   !$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) &   
   !$OMP PRIVATE(D,bVirt,AI,C,AIb,TMP) &
   !$OMP SHARED(ndimC,ndimD,nCMO1,nCMO2,nCMO3,Input,Output,nAI,nAI2,CMO3D) 
   DO C = 1,ndimC
      DO AI = 1,nAI2,64
         DO bVirt = 1,nCMO3
            DO AIb = 0,63
               Output(AI+AIb,bVirt,C) = 0.0E0_realk
            ENDDO
            DO D = 1,ndimD
               TMP = CMO3D(D,bVirt)
               DO AIb = 0,63
                  Output(AI+AIb,bVirt,C) = Output(AI+AIb,bVirt,C) + TMP*Input(AI+AIb,C,D)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
IF(MODAI)THEN
   !$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) &
   !$OMP PRIVATE(D,bVirt,C,TMP,AIb) &
   !$OMP SHARED(ndimC,ndimD,nCMO1,nCMO2,nCMO3,Input,Output,nAI,nAI2,CMO3D) 
   DO C = 1,ndimC
      DO bVirt = 1,nCMO3
         DO AIb = nAI2+1,nAI
            Output(AIb,bVirt,C) = 0.0E0_realk
         ENDDO
         DO D = 1,ndimD
            TMP = CMO3D(D,bVirt)
            DO AIb = nAI2+1,nAI
               Output(AIb,bVirt,C) = Output(AIb,bVirt,C) + TMP*Input(AIb,C,D)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
end subroutine DistributeLink_MOtransformC2

!scale nDimD*nOcc*nOcc*nVirt*nVirt called several times ntypes**4 times
!Normally it scales as 
!step1 : N*N*N*N*nOCC 
!step1 : N*N*N*nVirt*nOCC
!step1 : N*N*nOcc*nVirt*nOCC
!step4 : N*nVirt*nOcc*nVirt*nOCC
subroutine DistributeLink_MOtransformD(ndimD,&
     & CMO4,Input,Output,lupri,nCMO1,nCMO2,nCMO3,nCMO4)
  implicit none
  integer,intent(in) :: ndimD,lupri,nCMO1,nCMO2,nCMO3,nCMO4
  real(realk),intent(in) :: CMO4(ndimD,nCMO4)
  real(realk),intent(in) :: Input(nCMO1*nCMO2*nCMO3,ndimD)
  real(realk),intent(inout) :: Output(nCMO1*nCMO2*nCMO3,nCMO4)
  !
  integer :: jOcc,D,AIB,AIBb,nAIB,nAIB2 
  logical :: MODAIB 
  real(realk) :: TMP
!  This is the actual Output (nVirt,nOcc,nVirt,nOcc)
!   which was set to zero at beginning!
!     should not have
!     DO AIB = 1,nVirt*nOcc*nVirt  
!        Output(AIB,jOcc) = 0.0E0_realk
!     ENDDO

! I do not use AIB because that could easily require a 64 bit integer!
nAIB = nCMO1*nCMO2*nCMO3
IF(nCMO1*(nCMO2*(nCMO3*1_long)).GT.2147483640_long)call ichorQuit('int64 issue in MOtransformD',-1)
nAIB2 = (nAIB/64)*64
MODAIB = (MOD(nAIB,64).GT.0)
IF(nAIB2.GT.0)THEN
   !$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) &
   !$OMP PRIVATE(jOcc,D,AIB,AIBb,TMP) &
   !$OMP SHARED(ndimD,Input,Output,nAIB2,nCMO4,CMO4)
   DO AIB = 1,nAIB2,64
      DO jOcc = 1,nCMO4
         DO D = 1,ndimD
            TMP = CMO4(D,jOcc)
            DO AIBb = 0,63
               Output(AIB+AIBb,jOcc) = Output(AIB+AIBb,jOcc) + TMP*Input(AIB+AIBb,D)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
IF(MODAIB)THEN
   !$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) &
   !$OMP PRIVATE(jOcc,D,AIBb,TMP) &
   !$OMP SHARED(ndimD,Input,Output,nAIB2,nAIB,nCMO4,CMO4)
   DO jOcc = 1,nCMO4
      DO D = 1,ndimD
         TMP = CMO4(D,jOcc)
         DO AIBb = nAIB2+1,nAIB
            Output(AIBb,jOcc) = Output(AIBb,jOcc) + TMP*Input(AIBb,D)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
end subroutine DistributeLink_MOtransformD

subroutine DistributeLink_MOtransformDold(ndimD,&
     & CMO4,Input,Output,lupri,nCMO1,nCMO2,nCMO3,nCMO4)
  implicit none
  integer,intent(in) :: ndimD,lupri,nCMO1,nCMO2,nCMO3,nCMO4
  real(realk),intent(in) :: CMO4(ndimD,nCMO4)
  real(realk),intent(in) :: Input(nCMO1*nCMO2,nCMO3,ndimD)
  real(realk),intent(inout) :: Output(nCMO1*nCMO2,nCMO3,nCMO4)
  !
  integer :: jOcc,D,AI,B
  real(realk) :: TMP
  !$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) &
  !$OMP PRIVATE(jOcc,D,AI,B,TMP) &
  !$OMP SHARED(ndimD,Input,Output,nCMO1,nCMO2,nCMO3,nCMO4,CMO4)
  DO jOcc = 1,nCMO4
     DO D = 1,ndimD
        TMP = CMO4(D,jOcc)
        DO B = 1,nCMO3
           DO AI = 1,nCMO1*nCMO2
              Output(AI,B,jOcc) = Output(AI,B,jOcc) + TMP*Input(AI,B,D)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
end subroutine DistributeLink_MOtransformDold

END MODULE IchorEriLinkmod
