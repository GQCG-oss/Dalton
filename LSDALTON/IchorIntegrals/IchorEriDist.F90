!> @file
!> Contains the Distribution
MODULE IchorEriDistmodule
  use IchorprecisionModule
  use IchorCommonModule

CONTAINS
subroutine DetermineMaxPasses(nAtomsD,iBatchIndexOfTypeD,nAtomsC,nAtomsA,nAtomsB,&
     & iBatchIndexOfTypeC,iBatchIndexOfTypeA,nBatchB,nBatchA,iBatchIndexOfTypeB,&
     & TriangularRHSAtomLoop,CSscreen,TriangularLHSAtomLoop,TriangularODAtomLoop,&
     & noScreenCD2,BATCHGAB,THRESHOLD_CS,noScreenAB,BATCHGCD,nBatchC,nBatchD,MaxPasses)
  implicit none
  integer,intent(in) :: nAtomsD,iBatchIndexOfTypeD,nAtomsC
  integer,intent(in) :: nAtomsA,nAtomsB,iBatchIndexOfTypeC
  integer,intent(in) :: iBatchIndexOfTypeA,nBatchB,nBatchA
  integer,intent(in) :: iBatchIndexOfTypeB,nBatchC,nBatchD
  logical,intent(in) :: TriangularRHSAtomLoop,CSscreen
  logical,intent(in) :: TriangularLHSAtomLoop,TriangularODAtomLoop
  logical,intent(in) :: noScreenCD2(nAtomsC,nAtomsD)
  real(realk),intent(in) :: BATCHGAB(nBatchA,nBatchB),THRESHOLD_CS
  real(realk),intent(in) :: BATCHGCD(nBatchC,nBatchD)
  logical,intent(in) :: noScreenAB(natomsA,natomsB)
  integer,intent(inout) :: MaxPasses
  !
  !local variables
  integer :: iBatchD,IatomD,IatomC,nPasses
  real(realk) :: GABELM
  MaxPasses = 1
  DO IatomD = 1,nAtomsD
     GABELM = 0.0E0_realk
     iBatchD = iBatchIndexOfTypeD + IatomD
     DO IatomC = 1,nAtomsC
        IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
        IF(noScreenCD2(IatomC,IatomD))THEN
           nPasses = nAtomsA*nAtomsB
           IF(CSscreen)GABELM = BATCHGCD(iBatchIndexOfTypeC+IatomC,iBatchD)
           !output: IatomAPass,IatomBPass,nPasses
           CALL BUILD_noScreenRed(CSscreen,nAtomsA,nAtomsB,&
                & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
                & BATCHGAB,THRESHOLD_CS,GABELM,nPasses,&
                & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB)
           MaxPasses = MAX(MaxPasses,nPasses)
        ENDIF
     ENDDO
  ENDDO
end subroutine DetermineMaxPasses

SUBROUTINE BUILD_noScreenRed(CSscreen,nAtomsA,nAtomsB,&
     & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,BATCHGAB,&
     & THRESHOLD_CS,GABELM,nPasses,&
     & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,iAtomD,noScreenABin) 
  implicit none
  logical,intent(in) :: CSScreen,TriangularLHSAtomLoop,TriangularODAtomLoop
  integer,intent(in) :: nAtomsA,nAtomsB,iAtomC,iAtomD
  integer,intent(in) :: iBatchIndexOfTypeA,nBatchB,nBatchA
  integer,intent(in) :: iBatchIndexOfTypeB
  real(realk),intent(in) :: BATCHGAB(nBatchA,nBatchB),THRESHOLD_CS,GABELM
  logical,intent(in) :: noScreenABin(natomsA,natomsB)
  integer,intent(inout) :: nPasses
  !local variables
  integer :: iBatchA,IatomA,iBatchB,IatomB,iPass,IatomAstart,IatomBend
  iPass=0
!  IF(TriangularODAtomLoop)THEN
!     IatomAstart = iAtomC !Restrict AtomC =< AtomA
!  ELSE
     IatomAstart = 1     
!  ENDIF
  IatomBend = nAtomsB
  IF(CSScreen)THEN
     DO IatomA = IatomAstart,nAtomsA
        iBatchA = iBatchIndexOfTypeA + IatomA
        IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
        DO IatomB = 1,IatomBend
         IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
          IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
           IF(noScreenABin(IatomA,IatomB))THEN
            IF(GABELM*BATCHGAB(iBatchA,iBatchIndexOfTypeB + IatomB).GT.THRESHOLD_CS)THEN
               iPass = iPass + 1
            ENDIF
           ENDIF
          ENDIF
         ELSE
          IF(noScreenABin(IatomA,IatomB))THEN
           IF(GABELM*BATCHGAB(iBatchA,iBatchIndexOfTypeB + IatomB).GT.THRESHOLD_CS)THEN
              iPass = iPass + 1
           ENDIF
          ENDIF
         ENDIF
        ENDDO
     ENDDO
  ELSE
     DO IatomA = IatomAstart,nAtomsA
        IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
        DO IatomB = 1,IatomBend
         IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
          IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
           IF(noScreenABin(IatomA,IatomB))THEN
              iPass = iPass + 1
           ENDIF
          ENDIF
         ELSE
          IF(noScreenABin(IatomA,IatomB))THEN
             iPass = iPass + 1
          ENDIF
         ENDIF
        ENDDO
     ENDDO
  ENDIF
  nPasses = iPass
END SUBROUTINE BUILD_NOSCREENRed

SUBROUTINE BUILD_noScreen2(CSscreen,nAtomsA,nAtomsB,&
     & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,BATCHGAB,&
     & THRESHOLD_CS,GABELM,nPasses,IatomAPass,IatomBPass,MaxPasses,&
     & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,iAtomD,noScreenABin) 
  implicit none
  logical,intent(in) :: CSScreen,TriangularLHSAtomLoop,TriangularODAtomLoop
  integer,intent(in) :: nAtomsA,nAtomsB,iAtomC,iAtomD,MaxPasses
  integer,intent(in) :: iBatchIndexOfTypeA,nBatchB,nBatchA
  integer,intent(in) :: iBatchIndexOfTypeB
  real(realk),intent(in) :: BATCHGAB(nBatchA,nBatchB),THRESHOLD_CS,GABELM
  logical,intent(in) :: noScreenABin(natomsA,natomsB)
  integer,intent(inout) :: nPasses
  integer,intent(inout) :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  !local variables
  integer :: IatomAPass2(MaxPasses),IatomBPass2(MaxPasses),iPass2
  integer :: iBatchA,IatomA,iBatchB,IatomB,iPass,IatomAstart,IatomBend
  iPass=0
!  IF(TriangularODAtomLoop)THEN
!     IatomAstart = iAtomC !Restrict AtomC =< AtomA
!  ELSE
     IatomAstart = 1     
!  ENDIF
  IatomBend = nAtomsB
  IF(CSScreen)THEN
   DO IatomA = IatomAstart,nAtomsA
    iBatchA = iBatchIndexOfTypeA + IatomA
    IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
    DO IatomB = 1,IatomBend
     IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
      IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
       IF(noScreenABin(IatomA,IatomB))THEN
        IF(GABELM*BATCHGAB(iBatchA,iBatchIndexOfTypeB + IatomB).GT.THRESHOLD_CS)THEN
           iPass = iPass + 1
           IatomAPass2(iPass) = IatomA
           IatomBPass2(iPass) = IatomB
        ENDIF
       ENDIF
      ENDIF
     ELSE
      IF(noScreenABin(IatomA,IatomB))THEN
       IF(GABELM*BATCHGAB(iBatchA,iBatchIndexOfTypeB + IatomB).GT.THRESHOLD_CS)THEN
          iPass = iPass + 1
          IatomAPass2(iPass) = IatomA
          IatomBPass2(iPass) = IatomB
       ENDIF
      ENDIF
     ENDIF
    ENDDO
   ENDDO
  ELSE
   DO IatomA = IatomAstart,nAtomsA
    IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
    DO IatomB = 1,IatomBend
     IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
      IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
       IF(noScreenABin(IatomA,IatomB))THEN
          iPass = iPass + 1
          IatomAPass2(iPass) = IatomA
          IatomBPass2(iPass) = IatomB
       ENDIF
      ENDIF
     ELSE
      IF(noScreenABin(IatomA,IatomB))THEN
         iPass = iPass + 1
         IatomAPass2(iPass) = IatomA
         IatomBPass2(iPass) = IatomB
      ENDIF
     ENDIF
    ENDDO
   ENDDO
  ENDIF
  nPasses = iPass
  !sort IatomAPass,IatomBPass so that it is the same as 
  !Do iatomB = 1,nAtomsB
  ! Do iatomA = 1,nAtomsA   !which is how everything else is stored. 
  iPass2 = 1
  DO iPass = iPass,1,-1     
     IatomAPass(iPass2) = IatomAPass2(iPass)
     iAtomBPass(iPass2) = IatomBPass2(iPass)
     iPass2 = iPass2 + 1
  ENDDO
END SUBROUTINE BUILD_NOSCREEN2

subroutine TypeDistribution(PermuteLHSTypes,PermuteRHS,nOrbD,nOrbC,startC,startD,&
     & nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startOrbitalB,OutputDim1,OutputDim2,&
     & OutputDim3,OutputDim4,OutputStorage,LocalIntPass2,nOrbQ)
  implicit none
  logical,intent(in) :: PermuteLHSTypes,PermuteRHS
  integer,intent(in) :: nOrbD,nOrbC,startC,startD,nAtomsA,nAtomsB,nOrbA,nOrbB
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,nOrbQ        
  integer,intent(in) :: startOrbitalA(nAtomsA),startOrbitalB(nAtomsB)
  real(realk),intent(in) :: LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  !                                                                                                      
  integer :: iOrbQ,iOrbC,iOrbD,IatomB,i4,i3,startB

     IF(PermuteLHSTypes)THEN
      IF(PermuteRHS)THEN
!!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!!$OMP PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB) &
!!$OMP SHARED(startOrbitalB,nOrbD,nOrbC,nAtomsB,startD,&
!!$OMP        startC,nAtomsA,nOrbA,nOrbB,startOrbitalA,&
!!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
!!$OMP        OutputStorage,LocalIntPass2,nOrbQ)
!$OMP DO COLLAPSE(3) PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB)
       DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
         DO IatomB = 1,nAtomsB
          iOrbQ = iOrbC+(iOrbD-1)*nOrbC
          I4 = startD + iOrbD
          I3 = startC + iOrbC
          startB = startOrbitalB(iAtomB)
          CALL DIST1(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
               & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
               & OutputStorage,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iatomB)
         ENDDO
        ENDDO
       ENDDO
!$OMP END DO 
!!$OMP END PARALLEL DO
      ELSE  !PermuteLHSTypes NOT PermuteRHS
!!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!!$OMP PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB) &
!!$OMP SHARED(startOrbitalB,nOrbD,nOrbC,nAtomsB,startD,&
!!$OMP        startC,nAtomsA,nOrbA,nOrbB,startOrbitalA,&
!!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
!!$OMP        OutputStorage,LocalIntPass2,nOrbQ)
!$OMP DO COLLAPSE(3) PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB)
       DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
         DO IatomB = 1,nAtomsB
          iOrbQ = iOrbC+(iOrbD-1)*nOrbC
          I4 = startD + iOrbD
          I3 = startC + iOrbC
          startB = startOrbitalB(iAtomB)
          CALL DIST2(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
               & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
               & OutputStorage,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iatomB)
         ENDDO
        ENDDO
       ENDDO
!$OMP END DO 
!!$OMP END PARALLEL DO
      ENDIF
     ELSE !NOT PermuteLHSTypes
      IF(PermuteRHS)THEN
!!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!!$OMP PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB) &
!!$OMP SHARED(startOrbitalB,nOrbD,nOrbC,nAtomsB,startD,&
!!$OMP        startC,nAtomsA,nOrbA,nOrbB,startOrbitalA,&
!!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
!!$OMP        OutputStorage,LocalIntPass2,nOrbQ)
!$OMP DO COLLAPSE(3) PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB)
       DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
         DO IatomB = 1,nAtomsB
          iOrbQ = iOrbC+(iOrbD-1)*nOrbC
          I4 = startD + iOrbD
          I3 = startC + iOrbC
          startB = startOrbitalB(iAtomB)
          CALL DIST3(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
               & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
               & OutputStorage,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iatomB)
         ENDDO
        ENDDO
       ENDDO
!$OMP END DO 
!!$OMP END PARALLEL DO
      ELSE !NOT PermuteLHSTypes NOT PermuteRHS 
!!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!!$OMP PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB) &
!!$OMP SHARED(startOrbitalB,nOrbD,nOrbC,nAtomsB,startD,&
!!$OMP        startC,nAtomsA,nOrbA,nOrbB,startOrbitalA,&
!!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
!!$OMP        OutputStorage,LocalIntPass2,nOrbQ)
!$OMP DO COLLAPSE(3) PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB)
       DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
         DO IatomB = 1,nAtomsB
          iOrbQ = iOrbC+(iOrbD-1)*nOrbC
          I4 = startD + iOrbD
          I3 = startC + iOrbC
          startB = startOrbitalB(iAtomB)
          CALL DIST4(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
               & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
               & OutputStorage,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iatomB)
         ENDDO
        ENDDO
       ENDDO
!$OMP END DO 
!!$OMP END PARALLEL DO
      ENDIF
     ENDIF

   end subroutine TypeDistribution

!Reorder LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
!(including LHS permute) to LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
!this can be done on the accelerator
subroutine MainTriDistributetoLocalIntPass2CPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
     & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
     & MaxPasses,nPasses,TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1,&
     & LocalIntPass2,IatomAPass,iatomBPass,nContQ,nContP)
  implicit none
  integer,intent(in) :: TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,nContQ,nContP
  integer,intent(in) :: nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD
  integer,intent(in) :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
  logical,intent(in) :: TriangularLHSAtomLoop,Qsegmented,Psegmented
  real(realk),intent(in) :: LocalIntPass1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
  real(realk),intent(inout) :: LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
  integer,intent(in) :: IatomAPass(MaxPasses),iatomBPass(MaxPasses)
  IF(TriangularLHSAtomLoop)THEN
     IF(Qsegmented.AND.Psegmented)THEN
        IF(TotalAngmom.NE.0)THEN
           call TriDistributeCPUToLocalIntPassSeg(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                & nOrbCompD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,nPasses)
        ELSE !TotalAngmom=0
           call TriDistributeCPUToLocalIntPassSeg0000(LocalIntPass1,nAtomsA,nAtomsB,&
                & LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,nPasses)
        ENDIF
     ELSE
        IF(TotalAngmom.NE.0)THEN
           call TriDistributeCPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
                & MaxPasses,IatomAPass,iatomBPass,nPasses)
        ELSE !TotalAngmom=0
           call TriDistributeCPUToLocalIntPass0000(LocalIntPass1,&
                & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,&
                & IatomAPass,iatomBPass,nPasses)
        ENDIF
     ENDIF
  ELSE
     IF(Qsegmented.AND.Psegmented)THEN
        IF(TotalAngmom.NE.0)THEN
           call DistributeCPUToLocalIntPassSeg(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                & nOrbCompD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,&
                & nPasses)
        ELSE !TotalAngmom=0
           call DistributeCPUToLocalIntPassSeg0000(LocalIntPass1,nAtomsA,nAtomsB,&
                & LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,nPasses)
        ENDIF
     ELSE
        IF(TotalAngmom.NE.0)THEN
           call DistributeCPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
                & MaxPasses,IatomAPass,iatomBPass,nPasses)
        ELSE !TotalAngmom=0
           call DistributeCPUToLocalIntPass0000(LocalIntPass1,&
                & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,&
                & IatomAPass,iatomBPass,nPasses)
        ENDIF
     ENDIF
  ENDIF
end subroutine MainTriDistributetoLocalIntPass2CPU

!Reorder LocalIntPass(nContQ,nContP,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
!(including LHS permute) to LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
!this can be done on the accelerator
subroutine MainTriDistributetoLocalIntPass2GPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
     & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
     & MaxPasses,nPasses,TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1,&
     & LocalIntPass2,IatomAPass,iatomBPass,nContQ,nContP)
  implicit none
  integer,intent(in) :: TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,nContQ,nContP
  integer,intent(in) :: nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD
  integer,intent(in) :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
  logical,intent(in) :: TriangularLHSAtomLoop,Qsegmented,Psegmented
  real(realk),intent(in) :: LocalIntPass1(nContQ,nContP,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
  real(realk),intent(inout) :: LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
  integer,intent(in) :: IatomAPass(MaxPasses),iatomBPass(MaxPasses)
  IF(TriangularLHSAtomLoop)THEN
     IF(Qsegmented.AND.Psegmented)THEN
        IF(TotalAngmom.NE.0)THEN
           call TriDistributeGPUToLocalIntPassSeg(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                & nOrbCompD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,nPasses)
        ELSE !TotalAngmom=0
           call TriDistributeGPUToLocalIntPassSeg0000(LocalIntPass1,nAtomsA,nAtomsB,&
                & LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,nPasses)
        ENDIF
     ELSE
        IF(TotalAngmom.NE.0)THEN
           call TriDistributeGPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
                & MaxPasses,IatomAPass,iatomBPass,nPasses)
        ELSE !TotalAngmom=0
           call TriDistributeGPUToLocalIntPass0000(LocalIntPass1,&
                & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,&
                & IatomAPass,iatomBPass,nPasses)
        ENDIF
     ENDIF
  ELSE
     IF(Qsegmented.AND.Psegmented)THEN
        IF(TotalAngmom.NE.0)THEN
           call DistributeGPUToLocalIntPassSeg(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                & nOrbCompD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,&
                & nPasses)
        ELSE !TotalAngmom=0
           call DistributeGPUToLocalIntPassSeg0000(LocalIntPass1,nAtomsA,nAtomsB,&
                & LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,nPasses)
        ENDIF
     ELSE
        IF(TotalAngmom.NE.0)THEN
           call DistributeGPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
                & MaxPasses,IatomAPass,iatomBPass,nPasses)
        ELSE !TotalAngmom=0
           call DistributeGPUToLocalIntPass0000(LocalIntPass1,&
                & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,&
                & IatomAPass,iatomBPass,nPasses)
        ENDIF
     ENDIF
  ENDIF
end subroutine MainTriDistributetoLocalIntPass2GPU

subroutine DIST1(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
     & OutputStorage1,LocalIntPass1,nOrbQ,I3,I4,iOrbQ,iAtomB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB,I3,I4,iOrbQ
  integer,intent(in) :: startOrbitalA(nAtomsA),startB,iAtomB,nOrbQ
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4
  real(realk),intent(inout) :: OutputStorage1(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  real(realk),intent(in) :: LocalIntPass1(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  !local
  integer :: IatomA,startA,iOrbB,I2,iOrbA
  DO IatomA = 1,nAtomsA
   startA = startOrbitalA(iAtomA)
   DO iOrbB = 1,nOrbB
    I2 = startB + iOrbB
    DO iOrbA = 1,nOrbA
     OutputStorage1(iOrbA + startA,I2,I3,I4)=LocalIntPass1(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage1(iOrbA + startA,I2,I4,I3)=LocalIntPass1(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage1(I2,iOrbA + startA,I3,I4)=LocalIntPass1(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage1(I2,iOrbA + startA,I4,I3)=LocalIntPass1(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
    ENDDO
   ENDDO
  ENDDO
end subroutine DIST1

subroutine DIST2(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
     & OutputStorage2,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iAtomB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB,I3,I4,iOrbQ
  integer,intent(in) :: startOrbitalA(nAtomsA),startB,iAtomB,nOrbQ
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4
  real(realk),intent(inout) :: OutputStorage2(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  real(realk),intent(in) :: LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  !local
  integer :: IatomA,startA,iOrbB,I2,iOrbA
  DO IatomA = 1,nAtomsA
   startA = startOrbitalA(iAtomA)
   DO iOrbB = 1,nOrbB
    DO iOrbA = 1,nOrbA
     OutputStorage2(iOrbA + startA,iOrbB + startB,I3,I4)=LocalIntPass2(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage2(iOrbB + startB,iOrbA + startA,I3,I4)=LocalIntPass2(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
    ENDDO
   ENDDO
  ENDDO
end subroutine DIST2

subroutine DIST3(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
     & OutputStorage3,LocalIntPass3,nOrbQ,I3,I4,iOrbQ,iAtomB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB,I3,I4,iOrbQ
  integer,intent(in) :: startOrbitalA(nAtomsA),startB,iAtomB,nOrbQ
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4
  real(realk),intent(inout) :: OutputStorage3(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  real(realk),intent(in) :: LocalIntPass3(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  !local
  integer :: IatomA,startA,iOrbB,I2,iOrbA
  DO IatomA = 1,nAtomsA
   startA = startOrbitalA(iAtomA)
   DO iOrbB = 1,nOrbB
    DO iOrbA = 1,nOrbA
     OutputStorage3(iOrbA + startA,iOrbB + startB,I3,I4)=LocalIntPass3(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage3(iOrbA + startA,iOrbB + startB,I4,I3)=LocalIntPass3(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
    ENDDO
   ENDDO
  ENDDO
end subroutine DIST3

subroutine DIST4(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
     & OutputStorage4,LocalIntPass4,nOrbQ,I3,I4,iOrbQ,iAtomB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB,I3,I4,iOrbQ
  integer,intent(in) :: startOrbitalA(nAtomsA),startB,iAtomB,nOrbQ
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4
  real(realk),intent(inout) :: OutputStorage4(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  real(realk),intent(in) :: LocalIntPass4(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  !local
  integer :: IatomA,startA,iOrbB,I2,iOrbA
  DO IatomA = 1,nAtomsA
   startA = startOrbitalA(iAtomA)
   DO iOrbB = 1,nOrbB
    DO iOrbA = 1,nOrbA
     OutputStorage4(iOrbA + startA,iOrbB + startB,I3,I4)=LocalIntPass4(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
    ENDDO
   ENDDO
  ENDDO
end subroutine DIST4

subroutine DistributeIchor1(startA,startB,startC,startD,nOrbD,nOrbC,nOrbB,nOrbA,&
     & OutputStorage,LocalInt,Outputdim1,Outputdim2,Outputdim3,Outputdim4)
  implicit none
  integer,intent(in) :: startA,startB,startC,startD,nOrbA,nOrbB,nOrbC,nOrbD
  integer,intent(in) :: Outputdim1,Outputdim2,Outputdim3,Outputdim4
  real(realk),intent(in) :: LocalInt(nOrbA,nOrbB,nOrbC,nOrbD)
  real(realk),intent(inout) :: OutputStorage(Outputdim1,Outputdim2,Outputdim3,Outputdim4)
  !
  integer :: I2,I3,I4,iOrbA,iOrbB,iOrbC,iOrbD
  DO iOrbD = 1,nOrbD
   I4 = startD + iOrbD
   DO iOrbC = 1,nOrbC
    I3 = startC + iOrbC
    DO iOrbB = 1,nOrbB
     I2 = startB + iOrbB
     DO iOrbA = 1,nOrbA
      OutputStorage(iOrbA + startA,I2,I3,I4)=LocalInt(iOrbA,iOrbB,iOrbC,iOrbD)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine DistributeIchor1

subroutine IchorPermuteLHS1(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,LocalIntPass5,iAtomA,iAtomB,iOrbQ)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,iAtomA,iAtomB,iOrbQ
  integer,intent(in) :: nOrbA,nOrbB,nOrbC,nOrbD
  real(realk),intent(inout) :: LocalIntPass5(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC*nOrbD)  
  !local variables
  integer :: iOrbB,iOrbA
  DO iOrbB = 1,nOrbB
     DO iOrbA = 1,nOrbA
        LocalIntPass5(iOrbB,iAtomB,iOrbA,iAtomA,iOrbQ)=LocalIntPass5(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     ENDDO
  ENDDO
end subroutine IchorPermuteLHS1

!=====================================================================================================
!
!   Distribution routines to 
! LocalIntPass2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!
!=====================================================================================================
!
!      ======================
!      First the CPU routines - OpenMP
!      ======================
!
subroutine DistributeCPUToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
implicit none 
integer,intent(in)    :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
integer,intent(in)    :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
integer,intent(in)    :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
real(realk),intent(in)::LP1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB,nPasses)
real(realk),intent(inout)::LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!local variables
integer :: iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB
integer :: iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB

!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB,&
!$OMP         iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB) 
  DO IPass = 1,nPasses
   DO iContD = 1,nContD
    DO iContC = 1,nContC
     DO iAngD = 1,nOrbCompD
      DO iAngC = 1,nOrbCompC
       IatomB = IatomBPass(IPass)
       IatomA = IatomAPass(IPass)
       iContQ = iContC + (iContD-1)*nContC
       I4 = iAngD + (iContD-1)*nOrbCompD
       I3 = iAngC + (iContC-1)*nOrbCompC
       DO iContB = 1,nContB
        offsetB = (iContB-1)*nOrbCompB
        DO iContA = 1,nContA
         iContP = iContA+(iContB-1)*nContA
         offsetA = (iContA-1)*nOrbCompA
         DO iAngB = 1,nOrbCompB
          DO iAngA = 1,nOrbCompA
           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&
                & LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
end subroutine DistributeCPUToLocalIntPass

subroutine TriDistributeCPUToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
implicit none 
integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
integer,intent(in)        :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
real(realk),intent(in)   :: LP1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB,nPasses)
real(realk),intent(inout):: LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!local variables
integer :: iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,iContA,iContB
integer :: iContC,iContD,i4,i3,offsetA,offsetB,iAngB,iAngA

!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,iContA,iContB,iContC,&
!$OMP         iContD,i4,i3,offsetA,offsetB,iAngB,iAngA) 
  DO IPass = 1,nPasses
   DO iContP = 1,nContA*nContA
    DO iContQ = 1,nContC*nContD
     DO iAngD = 1,nOrbCompD
      DO iAngC = 1,nOrbCompC
       IatomB = IatomBPass(IPass)
       IatomA = IatomAPass(IPass)
       IF(IatomA.NE.IatomB)THEN
        !Ordering of Ipass is 
        !iPass = 0 
        !DO IatomA = 1,natomsA
        ! DO IatomB = 1,IatomBend
        !   iPass = iPass + 1
        ! ENDDO
        !ENDDO
        !Where IatomBend=IatomA for triangularLHSatomLoop
        !   or IatomBend=natomsB 
        
        iContA = iContP - ((iContP-1)/nContA)*nContA
        iContB = (iContP-1)/nContA+1
        iContC = iContQ - ((iContQ-1)/nContC)*nContC
        iContD = (iContQ-1)/nContC+1
        I4 = iAngD + (iContD-1)*nOrbCompD
        I3 = iAngC + (iContC-1)*nOrbCompC
        offsetB = (iContB-1)*nOrbCompA
        offsetA = (iContA-1)*nOrbCompA
        DO iAngB = 1,nOrbCompA
         DO iAngA = 1,nOrbCompA
          LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
          LP2(iAngB + offsetB,iatomB,iAngA + offsetA,iatomA,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
         ENDDO
        ENDDO
       ELSE
        iContA = iContP - ((iContP-1)/nContA)*nContA
        iContB = (iContP-1)/nContA+1
        iContC = iContQ - ((iContQ-1)/nContC)*nContC
        iContD = (iContQ-1)/nContC+1
        I4 = iAngD + (iContD-1)*nOrbCompD
        I3 = iAngC + (iContC-1)*nOrbCompC
        offsetB = (iContB-1)*nOrbCompA
        offsetA = (iContA-1)*nOrbCompA
!$OMP CRITICAL
        DO iAngB = 1,nOrbCompA
         DO iAngA = 1,nOrbCompA
          LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomA,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
          LP2(iAngB + offsetB,iatomA,iAngA + offsetA,iatomA,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
         ENDDO
        ENDDO
!$OMP END CRITICAL
       ENDIF
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
end subroutine TriDistributeCPUToLocalIntPass

subroutine TriDistributeCPUToLocalIntPassSeg0000(LP1,&
     & nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(MaxPasses)
  real(realk),intent(inout) :: LP2(nAtomsA,nAtomsB)
  !local variables
  integer :: iPass,iAtomA,iAtomB
!$OMP DO PRIVATE(iPass,iAtomA,iAtomB) 
  DO IPass = 1,nPasses
    IatomB = IatomBPass(IPass)
    IatomA = IatomAPass(IPass)    
    IF(IatomA.NE.IatomB)THEN
       LP2(iatomA,iatomB) = LP1(IPass)
       LP2(iatomB,iatomA) = LP1(IPass)
    ELSE
       LP2(iatomA,iatomA) = LP1(IPass)
    ENDIF
  ENDDO
!$OMP END DO
end subroutine TriDistributeCPUToLocalIntPassSeg0000

subroutine TriDistributeCPUToLocalIntPassSeg(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nOrbCompA,nOrbCompB,nOrbCompC*nOrbCompD,nPasses)
  real(realk),intent(inout) :: LP2(nOrbCompA,nAtomsA,nOrbCompB,nAtomsB,nOrbCompC*nOrbCompD)
  !local variables
  integer :: iPass,iAngQ,iAtomA,iAtomB,i4,i3,offsetA,offsetB,iAngB,iAngA

!$OMP DO COLLAPSE(2) &
!$OMP PRIVATE(iPass,iAngQ,iAtomA,iAtomB,i4,i3,offsetA,offsetB,iAngB,iAngA) 
  DO IPass = 1,nPasses
   DO iAngQ = 1,nOrbCompD*nOrbCompC
    IatomB = IatomBPass(IPass)
    IatomA = IatomAPass(IPass)    
    IF(IatomA.NE.IatomB)THEN
     DO iAngB = 1,nOrbCompA
      DO iAngA = 1,nOrbCompA
       LP2(iAngA,iatomA,iAngB,iatomB,iAngQ) = LP1(iAngA,iAngB,iAngQ,IPass)
       LP2(iAngB,iatomB,iAngA,iatomA,iAngQ) = LP1(iAngA,iAngB,iAngQ,IPass)
      ENDDO
     ENDDO
    ELSE
!$OMP CRITICAL
     DO iAngB = 1,nOrbCompA
      DO iAngA = 1,nOrbCompA
       LP2(iAngA,iatomA,iAngB,iatomA,iAngQ) = LP1(iAngA,iAngB,iAngQ,IPass)
       LP2(iAngB,iatomA,iAngA,iatomA,iAngQ) = LP1(iAngA,iAngB,iAngQ,IPass)
      ENDDO
     ENDDO
!$OMP END CRITICAL
    ENDIF
   ENDDO
  ENDDO
!$OMP END DO
end subroutine TriDistributeCPUToLocalIntPassSeg

subroutine TriDistributeCPUToLocalIntPass0000(LP1,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB
  integer,intent(in)        :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nContC*nContD,nContA,nContB,nPasses)
  real(realk),intent(inout) :: LP2(nContA,nAtomsA,nContB,nAtomsB,nContC*nContD)
  !local variables
  integer :: iPass,iContP,iContQ,iAtomA,iAtomB,iContA,iContB,iContC,iContD,i4,i3,offsetA,offsetB
!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iPass,iContP,iContQ,iAtomA,iAtomB,iContA,iContB, &
!$OMP         iContC,iContD,i4,i3,offsetA,offsetB) 
  DO IPass = 1,nPasses
   DO iContB = 1,nContA
    DO iContA = 1,nContA
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)    
     IF(IatomA.NE.IatomB)THEN
      DO iContQ = 1,nContC*nContD
       LP2(iContA,iatomA,iContB,iatomB,iContQ) = LP1(iContQ,iContA,iContB,IPass)
       LP2(iContB,iatomB,iContA,iatomA,iContQ) = LP1(iContQ,iContA,iContB,IPass)
      ENDDO
     ELSE
!$OMP CRITICAL
      DO iContQ = 1,nContC*nContD
       LP2(iContA,iatomA,iContB,iatomB,iContQ) = LP1(iContQ,iContA,iContB,IPass)
       LP2(iContB,iatomB,iContA,iatomA,iContQ) = LP1(iContQ,iContA,iContB,IPass)
      ENDDO
!$OMP END CRITICAL
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO

end subroutine TriDistributeCPUToLocalIntPass0000

subroutine DistributeCPUToLocalIntPass0000(LP1,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB,nContA,nContB,nContC,nContD,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nContC*nContD,nContA,nContB,nPasses)
  real(realk),intent(inout) :: LP2(nContA,nAtomsA,nContB,nAtomsB,nContC*nContD)
  !local variables
  integer :: iContA,iContB,iContQ,ipass,iAtomA,iAtomB
!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iContA,iContB,iContQ,ipass,iAtomA,iAtomB) 
  DO IPass = 1,nPasses
   DO iContB = 1,nContB
    DO iContA = 1,nContA
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)
     DO iContQ = 1,nContD*nContC
      LP2(iContA,iatomA,iContB,iatomB,iContQ) = LP1(iContQ,iContA,iContB,IPass)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
end subroutine DistributeCPUToLocalIntPass0000

subroutine DistributeCPUToLocalIntPassSeg0000(LP1,nAtomsA,nAtomsB,LP2,&
     & MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB,MaxPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses),nPasses
  real(realk),intent(in)    :: LP1(MaxPasses)
  real(realk),intent(inout) :: LP2(nAtomsA,nAtomsB)
  !local variables
  integer :: ipass,iAtomA,iAtomB
!$OMP DO PRIVATE(ipass,iAtomA,iAtomB) 
  DO IPass = 1,nPasses
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)
     LP2(iAtomA,iAtomB) = LP1(IPass)
  ENDDO
!$OMP END DO
end subroutine DistributeCPUToLocalIntPassSeg0000

subroutine DistributeCPUToLocalIntPassSeg(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,MaxPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses),nPasses
  real(realk),intent(in)    :: LP1(nOrbCompA,nOrbCompB,nOrbCompC*nOrbCompD,nPasses)
  real(realk),intent(inout) :: LP2(nOrbCompA,nAtomsA,nOrbCompB,nAtomsB,nOrbCompC*nOrbCompD)
  !local variables
  integer :: iAngA,iAngB,iAngQ,ipass,iAtomA,iAtomB
!$OMP DO COLLAPSE(3) &
!$OMP PRIVATE(iAngA,iAngB,iAngQ,ipass,iAtomA,iAtomB) 
  DO IPass = 1,nPasses
   DO iAngQ = 1,nOrbCompD*nOrbCompC
    DO iAngB = 1,nOrbCompB
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)
     DO iAngA = 1,nOrbCompA
      LP2(iAngA,iAtomA,iAngB,iAtomB,iAngQ) = LP1(iAngA,iAngB,iAngQ,IPass)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO

end subroutine DistributeCPUToLocalIntPassSeg
!
!      ====================
!      Now the GPU routines - OpenACC 
!      ====================
!
subroutine DistributeGPUToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
implicit none 
integer,intent(in)    :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
integer,intent(in)    :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
integer,intent(in)    :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
real(realk),intent(in)::LP1(nContC*nContD,nContA*nContB,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
real(realk),intent(inout)::LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!local variables
integer :: iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB
integer :: iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB
  DO IPass = 1,nPasses
   DO iContD = 1,nContD
    DO iContC = 1,nContC
     DO iAngD = 1,nOrbCompD
      DO iAngC = 1,nOrbCompC
       IatomB = IatomBPass(IPass)
       IatomA = IatomAPass(IPass)
       iContQ = iContC + (iContD-1)*nContC
       I4 = iAngD + (iContD-1)*nOrbCompD
       I3 = iAngC + (iContC-1)*nOrbCompC
       DO iContB = 1,nContB
        offsetB = (iContB-1)*nOrbCompB
        DO iContA = 1,nContA
         iContP = iContA+(iContB-1)*nContA
         offsetA = (iContA-1)*nOrbCompA
         DO iAngB = 1,nOrbCompB
          DO iAngA = 1,nOrbCompA
           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&
                & LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine DistributeGPUToLocalIntPass

subroutine TriDistributeGPUToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
implicit none 
integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
integer,intent(in)        :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
real(realk),intent(in)   :: LP1(nContC*nContD,nContA*nContB,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
real(realk),intent(inout):: LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!local variables
integer :: iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,iContA,iContB
integer :: iContC,iContD,i4,i3,offsetA,offsetB,iAngB,iAngA

  DO IPass = 1,nPasses
   DO iContP = 1,nContA*nContA
    DO iContQ = 1,nContC*nContD
     DO iAngD = 1,nOrbCompD
      DO iAngC = 1,nOrbCompC
       IatomB = IatomBPass(IPass)
       IatomA = IatomAPass(IPass)
       IF(IatomA.NE.IatomB)THEN
        !Ordering of Ipass is 
        !iPass = 0 
        !DO IatomA = 1,natomsA
        ! DO IatomB = 1,IatomBend
        !   iPass = iPass + 1
        ! ENDDO
        !ENDDO
        !Where IatomBend=IatomA for triangularLHSatomLoop
        !   or IatomBend=natomsB 
        
        iContA = iContP - ((iContP-1)/nContA)*nContA
        iContB = (iContP-1)/nContA+1
        iContC = iContQ - ((iContQ-1)/nContC)*nContC
        iContD = (iContQ-1)/nContC+1
        I4 = iAngD + (iContD-1)*nOrbCompD
        I3 = iAngC + (iContC-1)*nOrbCompC
        offsetB = (iContB-1)*nOrbCompA
        offsetA = (iContA-1)*nOrbCompA
        DO iAngB = 1,nOrbCompA
         DO iAngA = 1,nOrbCompA
          LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) = LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
          LP2(iAngB + offsetB,iatomB,iAngA + offsetA,iatomA,I3,I4) = LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
         ENDDO
        ENDDO
       ELSE
        iContA = iContP - ((iContP-1)/nContA)*nContA
        iContB = (iContP-1)/nContA+1
        iContC = iContQ - ((iContQ-1)/nContC)*nContC
        iContD = (iContQ-1)/nContC+1
        I4 = iAngD + (iContD-1)*nOrbCompD
        I3 = iAngC + (iContC-1)*nOrbCompC
        offsetB = (iContB-1)*nOrbCompA
        offsetA = (iContA-1)*nOrbCompA
        DO iAngB = 1,nOrbCompA
         DO iAngA = 1,nOrbCompA
          LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomA,I3,I4) = LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
          LP2(iAngB + offsetB,iatomA,iAngA + offsetA,iatomA,I3,I4) = LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
         ENDDO
        ENDDO
       ENDIF
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine TriDistributeGPUToLocalIntPass

subroutine TriDistributeGPUToLocalIntPassSeg0000(LP1,&
     & nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(MaxPasses)
  real(realk),intent(inout) :: LP2(nAtomsA,nAtomsB)
  !local variables
  integer :: iPass,iAtomA,iAtomB
  DO IPass = 1,nPasses
    IatomB = IatomBPass(IPass)
    IatomA = IatomAPass(IPass)    
    IF(IatomA.NE.IatomB)THEN
       LP2(iatomA,iatomB) = LP1(IPass)
       LP2(iatomB,iatomA) = LP1(IPass)
    ELSE
       LP2(iatomA,iatomA) = LP1(IPass)
    ENDIF
  ENDDO
end subroutine TriDistributeGPUToLocalIntPassSeg0000

subroutine TriDistributeGPUToLocalIntPassSeg(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nPasses,nOrbCompA,nOrbCompB,nOrbCompC*nOrbCompD)
  real(realk),intent(inout) :: LP2(nOrbCompA,nAtomsA,nOrbCompB,nAtomsB,nOrbCompC*nOrbCompD)
  !local variables
  integer :: iPass,iAngQ,iAtomA,iAtomB,i4,i3,offsetA,offsetB,iAngB,iAngA

  DO IPass = 1,nPasses
   DO iAngQ = 1,nOrbCompD*nOrbCompC
    IatomB = IatomBPass(IPass)
    IatomA = IatomAPass(IPass)    
    IF(IatomA.NE.IatomB)THEN
     DO iAngB = 1,nOrbCompA
      DO iAngA = 1,nOrbCompA
       LP2(iAngA,iatomA,iAngB,iatomB,iAngQ) = LP1(IPass,iAngA,iAngB,iAngQ)
       LP2(iAngB,iatomB,iAngA,iatomA,iAngQ) = LP1(IPass,iAngA,iAngB,iAngQ)
      ENDDO
     ENDDO
    ELSE
     DO iAngB = 1,nOrbCompA
      DO iAngA = 1,nOrbCompA
       LP2(iAngA,iatomA,iAngB,iatomA,iAngQ) = LP1(IPass,iAngA,iAngB,iAngQ)
       LP2(iAngB,iatomA,iAngA,iatomA,iAngQ) = LP1(IPass,iAngA,iAngB,iAngQ)
      ENDDO
     ENDDO
    ENDIF
   ENDDO
  ENDDO
end subroutine TriDistributeGPUToLocalIntPassSeg

subroutine TriDistributeGPUToLocalIntPass0000(LP1,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB
  integer,intent(in)        :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nContC*nContD,nContA,nContB,nPasses)
  real(realk),intent(inout) :: LP2(nContA,nAtomsA,nContB,nAtomsB,nContC*nContD)
  !local variables
  integer :: iPass,iContP,iContQ,iAtomA,iAtomB,iContA,iContB,iContC,iContD,i4,i3,offsetA,offsetB
  DO IPass = 1,nPasses
   DO iContB = 1,nContA
    DO iContA = 1,nContA
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)    
     IF(IatomA.NE.IatomB)THEN
      DO iContQ = 1,nContC*nContD
       LP2(iContA,iatomA,iContB,iatomB,iContQ) = LP1(iContQ,iContA,iContB,IPass)
       LP2(iContB,iatomB,iContA,iatomA,iContQ) = LP1(iContQ,iContA,iContB,IPass)
      ENDDO
     ELSE
      DO iContQ = 1,nContC*nContD
       LP2(iContA,iatomA,iContB,iatomB,iContQ) = LP1(iContQ,iContA,iContB,IPass)
       LP2(iContB,iatomB,iContA,iatomA,iContQ) = LP1(iContQ,iContA,iContB,IPass)
      ENDDO
     ENDIF
    ENDDO
   ENDDO
  ENDDO
end subroutine TriDistributeGPUToLocalIntPass0000

subroutine DistributeGPUToLocalIntPass0000(LP1,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB,nContA,nContB,nContC,nContD,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nContC*nContD,nContA,nContB,nPasses)
  real(realk),intent(inout) :: LP2(nContA,nAtomsA,nContB,nAtomsB,nContC*nContD)
  !local variables
  integer :: iContA,iContB,iContQ,ipass,iAtomA,iAtomB
  DO IPass = 1,nPasses
   DO iContB = 1,nContB
    DO iContA = 1,nContA
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)
     DO iContQ = 1,nContD*nContC
      LP2(iContA,iatomA,iContB,iatomB,iContQ) = LP1(iContQ,iContA,iContB,IPass)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine DistributeGPUToLocalIntPass0000

subroutine DistributeGPUToLocalIntPassSeg0000(LP1,nAtomsA,nAtomsB,LP2,&
     & MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB,MaxPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses),nPasses
  real(realk),intent(in)    :: LP1(nPasses)
  real(realk),intent(inout) :: LP2(nAtomsA,nAtomsB)
  !local variables
  integer :: ipass,iAtomA,iAtomB

  DO IPass = 1,nPasses
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)
     LP2(iAtomA,iAtomB) = LP1(IPass)
  ENDDO

end subroutine DistributeGPUToLocalIntPassSeg0000

subroutine DistributeGPUToLocalIntPassSeg(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,MaxPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses),nPasses
  real(realk),intent(in)    :: LP1(nPasses,nOrbCompA,nOrbCompB,nOrbCompC*nOrbCompD)
  real(realk),intent(inout) :: LP2(nOrbCompA,nAtomsA,nOrbCompB,nAtomsB,nOrbCompC*nOrbCompD)
  !local variables
  integer :: iAngA,iAngB,iAngQ,ipass,iAtomA,iAtomB
  DO IPass = 1,nPasses
   IatomB = IatomBPass(IPass)
   IatomA = IatomAPass(IPass)
   DO iAngQ = 1,nOrbCompD*nOrbCompC
    DO iAngB = 1,nOrbCompB
     DO iAngA = 1,nOrbCompA
      LP2(iAngA,iAtomA,iAngB,iAtomB,iAngQ) = LP1(IPass,iAngA,iAngB,iAngQ)
     ENDDO
    ENDDO
   ENDDO
  ENDDO

end subroutine DistributeGPUToLocalIntPassSeg
!==============================================================================
!
!   Done with Distribution routines to LocalIntPass2()
!
!==============================================================================

subroutine IchorPermuteLHS(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,LocalIntPass5)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  integer,intent(in) :: nOrbA,nOrbB,nOrbC,nOrbD
  real(realk),intent(inout) :: LocalIntPass5(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC*nOrbD)  
  !local variables
  integer :: iOrbQ,iAtomB,iOrbB,iatomA,iOrbA
!!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iOrbQ,iatomB,iOrbB,IatomA,&
!!$OMP iOrbA) FIRSTPRIVATE(nAtomsB,nAtomsA,nOrbA,nOrbB,nOrbC,&
!!$OMP nOrbD) SHARED(LocalIntPass5) SCHEDULE(DYNAMIC,1)
  DO iOrbQ = 1,nOrbC*nOrbD
   DO IatomB = 1,nAtomsB
    DO iOrbB = 1,nOrbB
     DO IatomA = 1,IatomB-1
      DO iOrbA = 1,nOrbA
       LocalIntPass5(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)=LocalIntPass5(iOrbB,iAtomB,iOrbA,iAtomA,iOrbQ)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!!$OMP END PARALLEL DO
end subroutine IchorPermuteLHS

subroutine IchorPermuteLHSSeg0000(nAtomsA,nAtomsB,LocalIntPass5)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  real(realk),intent(inout) :: LocalIntPass5(nAtomsA,nAtomsB)
  !local variables
  integer :: iAtomB,iatomA
!!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iatomB,&
!!$OMP IatomA) FIRSTPRIVATE(nAtomsB) SHARED(LocalIntPass5) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
     DO IatomA = 1,IatomB-1
        LocalIntPass5(iAtomA,iAtomB)=LocalIntPass5(iAtomB,iAtomA)
     ENDDO
  ENDDO
!!$OMP END PARALLEL DO
end subroutine IchorPermuteLHSSeg0000

subroutine IchorDistribute(nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,&
     & startC,startD,OutputStorage,dim1,dim2,dim3,dim4,&
     & LocalIntPass6,nOrbA,nOrbB,nOrbC,nOrbD,PermuteRHS,PermuteLHSTypes,&
     & TriangularODAtomLoop,lupri)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,startC,startD,dim1,dim2,dim3,dim4
  integer,intent(in) :: nOrbA,nOrbB,nOrbC,nOrbD,lupri
  integer,intent(in) :: startOrbitalA(nAtomsA),startOrbitalB(nAtomsB)
  logical,intent(in) :: PermuteRHS,PermuteLHSTypes,TriangularODAtomLoop
  real(realk),intent(inout) :: OutputStorage(dim1,dim2,dim3,dim4)
  real(realk),intent(in) :: LocalIntPass6(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
  !local variables
  integer :: iOrbD,i4,iOrbC,i3,iatomb,startb,iOrbB,i2,iatomA,startA,iOrbA
!!$OMP PARALLEL DEFAULT(none) PRIVATE(iOrbD,I4,iOrbC,I3,iatomB,startB,iOrbB,I2,IatomA,&
!!$OMP iOrbA,startA) FIRSTPRIVATE(nOrbD,nOrbC,nOrbB,nOrbA,nAtomsA,nAtomsB,PermuteLHSTypes,&
!!$OMP PermuteRHS) SHARED(startOrbitalB,startOrbitalA,startC,startD,LocalIntPass6,OutputStorage)
  IF(PermuteLHSTypes)THEN
   IF(PermuteRHS)THEN
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!!$OMP DO SCHEDULE(DYNAMIC,1)
       DO IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO iOrbB = 1,nOrbB
         I2 = startB + iOrbB
         DO IatomA = 1,nAtomsA
          startA = startOrbitalA(iAtomA)
          DO iOrbA = 1,nOrbA
           OutputStorage(iOrbA + startA,I2,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(iOrbA + startA,I2,I4,I3)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(I2,iOrbA + startA,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(I2,iOrbA + startA,I4,I3)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ELSE  !PermuteLHSTypes NOT PermuteRHS
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!!$OMP DO SCHEDULE(DYNAMIC,1)
       DO IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO iOrbB = 1,nOrbB
         I2 = startB + iOrbB
         DO IatomA = 1,nAtomsA
          startA = startOrbitalA(iAtomA)
          DO iOrbA = 1,nOrbA
           OutputStorage(iOrbA + startA,I2,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(I2,iOrbA + startA,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ENDIF
  ELSE !NOT PermuteLHSTypes
   IF(PermuteRHS)THEN
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!!$OMP DO SCHEDULE(DYNAMIC,1)
       DO IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO iOrbB = 1,nOrbB
         I2 = startB + iOrbB
         DO IatomA = 1,nAtomsA
          startA = startOrbitalA(iAtomA)
          DO iOrbA = 1,nOrbA
           OutputStorage(iOrbA + startA,I2,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(iOrbA + startA,I2,I4,I3)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ELSE !NOT PermuteLHSTypes NOT PermuteRHS 
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!!$OMP DO SCHEDULE(DYNAMIC,1)
       DO IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO iOrbB = 1,nOrbB
         I2 = startB + iOrbB
         DO IatomA = 1,nAtomsA
          startA = startOrbitalA(iAtomA)
          DO iOrbA = 1,nOrbA
           OutputStorage(iOrbA + startA,I2,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ENDIF
  ENDIF
!!$OMP END PARALLEL
end subroutine IchorDistribute

END MODULE IchorEriDistmodule
