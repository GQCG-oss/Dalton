!> @file
!> Contains some tools for the main Ichor integral drivers

MODULE IchorEriToolsmod
  use IchorprecisionMod
  use IchorCommonMod
  use IchorParametersMod

  logical,parameter :: UseSP = .TRUE.!SinglePrecision default double precision

CONTAINS
subroutine ObtainTypeInfo(nTypesD,ItypeDnon,OrderdListD,nAtomsOfTypeD,AngmomOfTypeD,nPrimOfTypeD,&
     & nContOfTypeD,ExtentOfTypeD,ItypeD,nAtomsD,AngmomD,nPrimD,nContD,extentD,nOrbCompD,nOrbD,&
     & nDimD,nCartOrbCompD,spherical)
  implicit none
  logical,intent(in) :: spherical
  integer,intent(in) :: nTypesD,ItypeDnon
  real(realk),intent(inout) :: extentD
  integer,intent(inout) :: nAtomsD,AngmomD,nPrimD,nContD,nCartOrbCompD
  integer,intent(inout) :: nOrbCompD,nOrbD,nDimD,ItypeD
  integer,intent(in) :: OrderdListD(nTypesD),nAtomsOfTypeD(nTypesD),AngmomOfTypeD(nTypesD)
  integer,intent(in) :: nPrimOfTypeD(nTypesD),nContOfTypeD(nTypesD)
  real(realk),intent(in) :: ExtentOfTypeD(nTypesD)
  
  ItypeD = OrderdListD(ItypeDnon)  
  nAtomsD = nAtomsOfTypeD(ItypeD)
  AngmomD = AngmomOfTypeD(ItypeD)
  nPrimD = nPrimOfTypeD(ItypeD)
  nContD = nContOfTypeD(ItypeD)
  extentD = ExtentOfTypeD(ItypeD)
  nCartOrbCompD = (AngmomD+1)*(AngmomD+2)/2
  IF (spherical) THEN
     nOrbCompD = 2*(AngmomD+1)-1
  ELSE
     nOrbCompD = nCartOrbCompD
  ENDIF
  nOrbD = nContD*nOrbCompD
  nDimD = nContD*nOrbCompD*nAtomsD
end subroutine ObtainTypeInfo

subroutine ObtainTypeInfoNoExtent(nTypesD,nAtomsOfTypeD,AngmomOfTypeD,nPrimOfTypeD,&
     & nContOfTypeD,ItypeD,nAtomsD,AngmomD,nPrimD,nContD,nOrbCompD,nOrbD,&
     & nDimD,nCartOrbCompD,spherical)
  implicit none
  logical,intent(in) :: spherical
  integer,intent(in) :: nTypesD,ItypeD
  integer,intent(inout) :: nAtomsD,AngmomD,nPrimD,nContD,nCartOrbCompD
  integer,intent(inout) :: nOrbCompD,nOrbD,nDimD
  integer,intent(in) :: nAtomsOfTypeD(nTypesD),AngmomOfTypeD(nTypesD)
  integer,intent(in) :: nPrimOfTypeD(nTypesD),nContOfTypeD(nTypesD)  
  nAtomsD = nAtomsOfTypeD(ItypeD)
  AngmomD = AngmomOfTypeD(ItypeD)
  nPrimD = nPrimOfTypeD(ItypeD)
  nContD = nContOfTypeD(ItypeD)
  nCartOrbCompD = (AngmomD+1)*(AngmomD+2)/2
  IF (spherical) THEN
     nOrbCompD = 2*(AngmomD+1)-1
  ELSE
     nOrbCompD = nCartOrbCompD
  ENDIF
  nOrbD = nContD*nOrbCompD
  nDimD = nContD*nOrbCompD*nAtomsD
end subroutine ObtainTypeInfoNoExtent

subroutine MakeCombCMO34(nDimC,nDimD,nCMO3,nCMO4,CMO3,CMO4,CombCMO34)
  implicit none
  integer,intent(in) :: nDimC,nDimD,nCMO3,nCMO4
  real(realk),intent(in) :: CMO3(nDimC,nCMO3),CMO4(nDimD,nCMO4)
  real(realk),intent(inout) :: CombCMO34(nDimC,nDimD,nCMO3,nCMO4)
  !
  integer :: J,B,ID,IC
  real(realk) :: TMP
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
  !$OMP PRIVATE(J,B,ID,IC,TMP) &
  !$OMP SHARED(nDimC,nDimD,nCMO3,nCMO4,CMO3,CMO4,CombCMO34)
  DO J=1,nCMO4
     DO B=1,nCMO3
        DO ID=1,ndimD
           TMP = CMO4(ID,J)
           DO IC=1,ndimC
              CombCMO34(IC,ID,B,J) = CMO3(IC,B)*TMP
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO 
end subroutine MakeCombCMO34

!Make from Screening AB (GAB a list of A atoms for given B atom)
!This list is restricted by symmetry to be AtomB =< AtomA
subroutine ConstructBraList(nBraList,BraList,nAtomsA,nAtomsB,THRESHOLD_CS,maxGABRHS,&
     & TriangularLHSAtomLoop,atomGAB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  integer,intent(inout) :: nBraList(nAtomsB),BraList(nAtomsA,nAtomsB)
  real(realk),intent(in) :: THRESHOLD_CS,maxGABRHS,atomGAB(nAtomsA,nAtomsB)
  logical,intent(in) :: TriangularLHSAtomLoop
  !
  integer :: IatomB,IatomA,I,IatomAstart
  real(realk) :: CRITERIA,Sorting(nAtomsA)
  CRITERIA = THRESHOLD_CS/maxGABRHS
  DO IatomB = 1,nAtomsB
     I = 0 
     IatomAstart = 1
     IF(TriangularLHSAtomLoop) IatomAstart = IatomB
     DO IatomA = IatomAstart,nAtomsA
        IF(atomGAB(IatomA,IatomB).GT.CRITERIA)THEN
           I = I + 1
           BraList(I,IatomB) = IatomA
           Sorting(I) = atomGAB(IatomA,IatomB)
        ENDIF
     ENDDO
     nBraList(IatomB) = I
     call IchorDecreasingSortingWindex(Sorting(1:I),nBraList(IatomB),BraList(1:I,IatomB))
  ENDDO
end subroutine ConstructBraList

!Make from Screening info a list of contributing B atoms for given D atom)
subroutine ConstructBraketList(nAtomsB,nAtomsD,nBraketList,BraketList,TriangularODAtomLoop,&
     & ReducedDmatBD,MaxGABVec,MaxGCDVec,THRESHOLD_CS)
  implicit none
  logical,intent(in) :: TriangularODAtomLoop
  integer,intent(in) :: nAtomsB,nAtomsD
  integer,intent(inout) :: nBraketList(nAtomsD),BraketList(nAtomsB,nAtomsD)
  real(realk),intent(in) :: ReducedDmatBD(nAtomsB,nAtomsD),MaxGABVec(nAtomsB)
  real(realk),intent(in) :: MaxGCDVec(nAtomsD),THRESHOLD_CS
  !
  integer :: iAtomD,startBatom,I,iAtomB
  real(realk) :: SORTING(nAtomsB)
  Do iAtomD=1,nAtomsD
     IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
        startBatom = iAtomD
     ELSE
        startBatom = 1
     ENDIF
     I=0 
     Do iAtomB = startBatom,nAtomsB
        IF(ReducedDmatBD(iAtomB,iAtomD)*MaxGABVec(iAtomB)*MaxGCDVec(iAtomD).GT.THRESHOLD_CS)THEN
           I=I+1
           BraketList(I,iAtomD) = iAtomB
           SORTING(I) = ReducedDmatBD(iAtomB,iAtomD)*MaxGCDVec(iAtomD)
        ENDIF
     Enddo
     nBraketList(iAtomD) = I
     call IchorDecreasingSortingWindex(Sorting(1:I),I,BraketList(1:I,IatomD))
  Enddo
end subroutine ConstructBraketList

subroutine FormActiveDmat(nDimB,nDimD,DmatBD,nAtomsB,nAtomsD,&
     & startOrbitalB,startOrbitalD,nOrbB,nOrbD,Dfull,nd1,nd2,nd3)
  implicit none
  integer,intent(in) :: nDimB,nDimD,nAtomsB,nAtomsD,nOrbB,nOrbD,nd1,nd2,nd3
  integer,intent(in) :: startOrbitalB(nAtomsB),startOrbitalD(nAtomsD)
  real(realk),intent(in) :: Dfull(nd1,nd2,nd3)
  real(realk),intent(inout) :: DmatBD(nDimB,nDimD,nd3)
  !
  integer :: iD3,IatomD,IatomB,IorbD,IorbB,startLocD,startLocB,startB,startD

  do iD3 = 1,nd3
   startLocD = 0
   do IatomD = 1,nAtomsD              
    startD = startOrbitalD(iAtomD)
    startLocB = 0
    do IatomB = 1,nAtomsB
     startB = startOrbitalB(iAtomB)
     DO IorbD = 1,nOrbD
      DO IorbB = 1,nOrbB
       DmatBD(startLocB + IorbB,startLocD + IorbD,iD3) = Dfull(startB+IorbB,startD+IorbD,iD3)
      ENDDO
     ENDDO
     startLocB = startLocB + nOrbB          
    enddo
    startLocD = startLocD + nOrbD          
   enddo
  enddo
end subroutine FormActiveDmat

subroutine addActiveKmatToOutput(nDimA,nDimB,nDimC,nDimD,nAtomsA,nAtomsB,nAtomsC,nAtomsD,&
     & nOrbA,nOrbB,nOrbC,nOrbD,startOrbitalA,startOrbitalB,startOrbitalC,startOrbitalD,&
     & OutputDim1,OutputDim2,OutputDim5,OutputStorage,KmatAC)
  implicit none
  integer,intent(in) ::  nDimA,nDimB,nDimC,nDimD,nAtomsA,nAtomsB,nAtomsC,nAtomsD
  integer,intent(in) ::  nOrbA,nOrbB,nOrbC,nOrbD
  integer,intent(in) ::  startOrbitalA(nAtomsA),startOrbitalB(nAtomsB)
  integer,intent(in) ::  startOrbitalC(nAtomsC),startOrbitalD(nAtomsD)
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim5
  real(realk) :: OutputStorage(OutputDim1,OutputDim2,OutputDim5)
  real(realk) :: KmatAC(nDimA,nDimC,OutputDim5)
  !
  integer :: idim5,IatomC,IatomA,IorbC,IorbA,startLocC,startLocA,startA,startC
  
  do idim5 = 1, OutputDim5
   startLocC = 0
   do IatomC = 1,nAtomsC
    startC = startOrbitalC(iAtomC)                 
    startLocA = 0
    do IatomA = 1,nAtomsA
     startA = startOrbitalA(iAtomA)                 
     DO IorbC = 1,nOrbC
      DO IorbA = 1,nOrbA
       OutputStorage(startA+IorbA,startC+IorbC,idim5) = KmatAC(startLocA+IorbA,startLocC+IorbC,idim5)
      ENDDO
     ENDDO
     startLocA = startLocA + nOrbA
    enddo
    startLocC = startLocC + nOrbC
   enddo
  enddo
end subroutine addActiveKmatToOutput

subroutine FormReducedDmat(DmatBD,ReducedDmatBD,nOrbB,nOrbD,nAtomsB,nAtomsD,nDmat)
  implicit none
  integer,intent(in) :: nOrbB,nOrbD,nAtomsB,nAtomsD,nDmat
  real(realk),intent(in) :: DmatBD(nOrbB,nAtomsB,nOrbD,nAtomsD,nDmat)
  real(realk),intent(inout) :: ReducedDmatBD(nAtomsB,nAtomsD)
  !
  integer :: IatomD,IatomB,D,B,idmat
  real(realk) :: TMP
  do IatomD=1,nAtomsD
     do IatomB=1,nAtomsB
        TMP = 0.0E0_realk
        do idmat = 1,nDmat
           do D=1,nOrbD
              do B=1,nOrbB
                 TMP = MAX(TMP,ABS(DmatBD(B,IatomB,D,IatomD,idmat)))
              enddo
           enddo
        enddo
        ReducedDmatBD(IatomB,IatomD) = TMP
     enddo
  enddo
end subroutine FormReducedDmat

subroutine Build_AtomGAB(nAtomsA,nAtomsB,nBatchA,nBatchB,ntypesA,ntypesB,ItypeA,ItypeB,&
     & BatchIndexOfTypeA,BatchIndexOfTypeB,BATCHGAB,AtomGAB,MaxGABvec,MaxAtomGabelmLHS)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,ItypeA,ItypeB,ntypesA,ntypesB
  integer,intent(in) :: nBatchA,nBatchB
  real(realk),intent(in) :: BATCHGAB(nBatchA,nBatchB)
  real(realk),intent(inout) :: AtomGAB(nAtomsA,nAtomsB)
  integer,intent(in) :: BatchIndexOfTypeA(ntypesA),BatchIndexOfTypeB(ntypesB)
  real(realk),intent(inout) :: maxGABvec(nAtomsB),MaxAtomGabelmLHS
  !
  integer :: IatomB,IatomA,iBatchA,iBatchB
  real(realk) :: TMP
  iBatchB = BatchIndexOfTypeB(ItypeB)
  iBatchA = BatchIndexOfTypeA(ItypeA)
  DO IatomB = 1,nAtomsB
     DO IatomA = 1,nAtomsA
        AtomGAB(IatomA,IatomB) = BATCHGAB(iBatchA+IatomA,iBatchB+IatomB)
     ENDDO
     maxGABvec(IAtomB) = MAXVAL(AtomGAB(:,IatomB))
  ENDDO
  MaxAtomGabelmLHS = MAXVAL(maxGABvec)
end subroutine Build_AtomGAB

!FIXME replace with improved sorting routine - bubble sort - bucket sort ,..
!possible put into buckets with 10-0, 0-0.1, 0.1-0.01,..... and sort each bucket using this
subroutine IchorDecreasingSortingWindex(Sorting,n,TrackIndex)
  implicit none
  integer :: n
  real(realk),intent(inout) :: Sorting(n)
  integer,intent(inout) :: TrackIndex(n)
  !
  real(realk) :: tmp
  integer :: tmp1,i
  logical :: swp
  swp=.true.
  do while (swp)
     swp=.false.
     do i=1,n-1
        if(Sorting(i) < Sorting(i+1)) then ! reverse order
           
           tmp = Sorting(i+1)
           Sorting(i+1) = Sorting(i)
           Sorting(i) = tmp
           
           tmp1 = TrackIndex(i+1)
           TrackIndex(i+1) = TrackIndex(i)
           TrackIndex(i) = tmp1
           
           swp=.true.
        end if
     end do
  end do
end subroutine IchorDecreasingSortingWindex

subroutine build_TYPESTRING(TYPESTRING,AngmomA,AngmomB,AngmomC,AngmomD,&
     & nTypesA,nTypesB,nTypesC,nTypesD,ItypeA,ItypeB,ItypeC,ItypeD)
  implicit none
  character(len=16),intent(inout)  :: TYPESTRING
  integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
  integer,intent(in) :: nTypesA,nTypesB,nTypesC,nTypesD
  integer,intent(in) :: ItypeA,ItypeB,ItypeC,ItypeD
  integer :: n  
  TYPESTRING(1:3) = 'ANG'
  WRITE(TYPESTRING(4:4),'(A1)') AngmomA
  WRITE(TYPESTRING(5:5),'(A1)') AngmomB
  WRITE(TYPESTRING(6:6),'(A1)') AngmomC
  WRITE(TYPESTRING(7:7),'(A1)') AngmomD
  TYPESTRING(8:8) = 'T'
  n=9
  call typeAspec(ItypeA,nTypesA,n,TYPESTRING)
  call typeAspec(ItypeB,nTypesB,n,TYPESTRING)
  call typeAspec(ItypeC,nTypesC,n,TYPESTRING)
  call typeAspec(ItypeD,nTypesD,n,TYPESTRING)
end subroutine build_TYPESTRING

subroutine typeAspec(ItypeA,nTypesA,n,TYPESTRING)
  integer,intent(in) :: ItypeA,nTypesA
  integer,intent(inout) :: n
  character(len=16),intent(inout)  :: TYPESTRING
  IF(nTypesA.LT.10)THEN
     WRITE(TYPESTRING(n:n),'(A1)') ItypeA
     n=n+1
  ELSEIF(nTypesA.LT.100)THEN
     IF(ItypeA.LT.10)THEN
        TYPESTRING(n:n) = '0'
        WRITE(TYPESTRING(n+1:n+1),'(A1)') ItypeA
     ELSE
        WRITE(TYPESTRING(n:n+1),'(A1)') ItypeA
     ENDIF
     n=n+2
  ELSE
     TYPESTRING(n:n+1) = '  '
     n=n+2
  ENDIF
end subroutine typeAspec

subroutine determinePermuteSym(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
implicit none
integer,intent(in) :: IchorPermuteSpec
logical,intent(inout) ::  SameLHSaos,SameRHSaos,SameODs
SELECT CASE(IchorPermuteSpec)
CASE(IchorPermuteTTT)
   SameLHSaos=.TRUE.;  SameRHSaos=.TRUE. ; SameODs = .TRUE.
CASE(IchorPermuteFFT)
   SameLHSaos=.FALSE.; SameRHSaos=.FALSE.; SameODs = .TRUE.
CASE(IchorPermuteTTF)
   SameLHSaos=.TRUE.;  SameRHSaos=.TRUE. ; SameODs = .FALSE.
CASE(IchorPermuteTFF)
   SameLHSaos=.TRUE.;  SameRHSaos=.FALSE.; SameODs = .FALSE.
CASE(IchorPermuteFTF)
   SameLHSaos=.FALSE.; SameRHSaos=.TRUE. ; SameODs = .FALSE.
CASE(IchorPermuteFFF)
   SameLHSaos=.FALSE.; SameRHSaos=.FALSE.; SameODs = .FALSE.
CASE DEFAULT
   call ichorquit('unknown case in determinePermuteSym',-1)
END SELECT
end subroutine determinePermuteSym

subroutine PermuteRHStypesSub(dim1dim2,dim3,dim4,OutputStorage,nAtomsC,nAtomsD,&
     & startOrbitalC,startOrbitalD,nOrbC,nOrbD,TriangularRHSAtomLoop,lupri)
  implicit none
  integer,intent(in) :: dim1dim2,dim3,dim4,nAtomsC,nAtomsD,nOrbC,nOrbD,lupri
  integer,intent(in) :: startOrbitalC(nAtomsC),startOrbitalD(nAtomsD)
  real(realk),intent(inout) :: OutputStorage(dim1dim2,dim3,dim4)
  logical,intent(in) :: TriangularRHSAtomLoop
  !
  integer :: I,C,D,startC,startD,IatomC,IatomD,A,B
  IF(TriangularRHSAtomLoop)THEN
!$OMP PARALLEL DEFAULT(none) PRIVATE(IatomD,startD,IatomC,startC,C,D,&
!$OMP I) SHARED(startOrbitalD,startOrbitalC,OutputStorage) FIRSTPRIVATE(nAtomsD,&
!$OMP nAtomsC,nOrbC,nOrbD,dim1dim2) 
   DO IatomC = 1,nAtomsC
    startC = startOrbitalC(iAtomC)
    DO IatomD = 1,IatomC
     startD = startOrbitalD(iAtomD)
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO D=startD+1,startD+nOrbD
      DO C=startC+1,startC+nOrbC
       DO I=1,dim1dim2
        OutputStorage(I,D,C) = OutputStorage(I,C,D)
       ENDDO
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
!$OMP END PARALLEL
  ELSE
!$OMP PARALLEL DEFAULT(none) PRIVATE(IatomD,startD,IatomC,startC,C,D,&
!$OMP I) SHARED(startOrbitalD,startOrbitalC,OutputStorage) FIRSTPRIVATE(nAtomsD,&
!$OMP nAtomsC,nOrbC,nOrbD,dim1dim2)
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
!$OMP DO SCHEDULE(DYNAMIC,1)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
     DO D=startD+1,startD+nOrbD
      DO C=startC+1,startC+nOrbC
       DO I=1,dim1dim2
        OutputStorage(I,D,C) = OutputStorage(I,C,D)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
!$OMP END DO NOWAIT
   ENDDO
!$OMP END PARALLEL
  ENDIF
end subroutine PermuteRHStypesSub

subroutine PermuteODtypesSub(dim1,dim2,dim3,dim4,OutputStorage,&
     & nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,nOrbA,nOrbB,&
     & nAtomsC,nAtomsD,startOrbitalC,startOrbitalD,nOrbC,nOrbD,&
     & TriangularODAtomLoop,TriangularRHSAtomLoop,&
     & TriangularLHSAtomLoop,SameLHSaos,SameRHSaos,lupri)
  implicit none
  integer,intent(in) :: dim1,dim2,dim3,dim4,lupri
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB
  integer,intent(in) :: nAtomsC,nAtomsD,nOrbC,nOrbD
  integer,intent(in) :: startOrbitalA(nAtomsA),startOrbitalB(nAtomsB)
  integer,intent(in) :: startOrbitalC(nAtomsC),startOrbitalD(nAtomsD)
  real(realk),intent(inout) :: OutputStorage(dim1,dim2,dim3,dim4)
  logical,intent(in) :: TriangularODAtomLoop,TriangularRHSAtomLoop
  logical,intent(in) :: TriangularLHSAtomLoop,SameLHSaos,SameRHSaos
  !
  integer :: A,B,C,D,startC,startD,IatomC,IatomD,IatomB,IatomA,startA,startB
  integer :: IatomBend
!$OMP PARALLEL DEFAULT(none) PRIVATE(IatomB,IatomBend,iatomD,IatomC,startC,startD,&
!$OMP IatomA,startA,startB,A,B,C,D) SHARED(startOrbitalD,startOrbitalC,startOrbitalA,&
!$OMP startOrbitalB,TriangularODAtomLoop,OutputStorage) FIRSTPRIVATE(nOrbA,nOrbB,nOrbC,&
!$OMP nOrbD,nAtomsC,nAtomsD,nAtomsB,nAtomsA,TriangularLHSAtomLoop,&
!$OMP TriangularRHSAtomLoop,SameRHSaos,SameLHSaos)
  IF(SameLHSaos.AND.SameRHSaos)THEN
   IatomBend = nAtomsB
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
     IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      IF(TriangularLHSAtomLoop)IatomBend = IatomA
      DO IatomB = 1,IatomBend       
       startB = startOrbitalB(iAtomB)
       IF(TriangularODAtomLoop)THEN
        IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
         DO D=startD+1,startD+nOrbD
          DO C=startC+1,startC+nOrbC
           DO B=startB+1,startB+nOrbB
            DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,B,A) = OutputStorage(A,B,C,D)
              OutputStorage(C,D,B,A) = OutputStorage(A,B,C,D) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ELSE
        DO D=startD+1,startD+nOrbD
         DO C=startC+1,startC+nOrbC
          DO B=startB+1,startB+nOrbB
           DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,B,A) = OutputStorage(A,B,C,D)
              OutputStorage(C,D,B,A) = OutputStorage(A,B,C,D)
           ENDDO !A
          ENDDO !B
         ENDDO !C
        ENDDO !D
       ENDIF
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
 ELSEIF(SameLHSaos)THEN
   IatomBend = nAtomsB
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      IF(TriangularLHSAtomLoop)IatomBend = IatomA
      DO IatomB = 1,IatomBend       
       startB = startOrbitalB(iAtomB)
       IF(TriangularODAtomLoop)THEN
        IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
         DO D=startD+1,startD+nOrbD
          DO C=startC+1,startC+nOrbC
           DO B=startB+1,startB+nOrbB
            DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(C,D,B,A) = OutputStorage(A,B,C,D)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ELSE
        DO D=startD+1,startD+nOrbD
         DO C=startC+1,startC+nOrbC
          DO B=startB+1,startB+nOrbB
           DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(C,D,B,A) = OutputStorage(A,B,C,D)
           ENDDO !A
          ENDDO !B
         ENDDO !C
        ENDDO !D
       ENDIF
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
 ELSEIF(SameRHSaos)THEN
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
     IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      DO IatomB = 1,nAtomsB
       startB = startOrbitalB(iAtomB)
       IF(TriangularODAtomLoop)THEN
        IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
         DO D=startD+1,startD+nOrbD
          DO C=startC+1,startC+nOrbC
           DO B=startB+1,startB+nOrbB
            DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,A,B) = OutputStorage(A,B,C,D)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ELSE
        DO D=startD+1,startD+nOrbD
         DO C=startC+1,startC+nOrbC
          DO B=startB+1,startB+nOrbB
           DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,A,B) = OutputStorage(A,B,C,D)
           ENDDO !A
          ENDDO !B
         ENDDO !C
        ENDDO !D
       ENDIF
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
 ELSE
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      DO IatomB = 1,nAtomsB
       startB = startOrbitalB(iAtomB)
       IF(TriangularODAtomLoop)THEN
        IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
         DO D=startD+1,startD+nOrbD
          DO C=startC+1,startC+nOrbC
           DO B=startB+1,startB+nOrbB
            DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ELSE
        DO D=startD+1,startD+nOrbD
         DO C=startC+1,startC+nOrbC
          DO B=startB+1,startB+nOrbB
           DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
           ENDDO !A
          ENDDO !B
         ENDDO !C
        ENDDO !D
       ENDIF
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
 ENDIF
!$OMP END PARALLEL
end subroutine PermuteODtypesSub

subroutine ExtractBatchGabFromFullGab(nBatchASmall,nBatchBSmall,&
     & BATCHGABSmall,nBatchABIG,nBatchBBIG,BATCHGABBig,&
     & startBatchA,endBatchA,startBatchB,endBatchB)
implicit none
integer,intent(in) :: nBatchASmall,nBatchBSmall
integer,intent(in) :: nBatchABIG,nBatchBBIG
integer,intent(in) :: startBatchA,endBatchA,startBatchB,endBatchB
real(realk),intent(in) :: BATCHGABBig(nBatchABIG,nBatchBBIG) 
real(realk),intent(inout) :: BATCHGABSmall(nBatchASmall,nBatchBSmall) 
!local
integer :: sA,sB,IBB,IA,IB
sA = startBatchA-1
sB = startBatchB-1
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(IB,IBB,&
!$OMP IA) SHARED(BATCHGABSmall,BATCHGABBig) FIRSTPRIVATE(sA,&
!$OMP sB,nBatchASmall,nBatchBSmall) SCHEDULE(DYNAMIC,3)
DO IB=1,nBatchBSmall
   IBB = sB+IB
   DO IA=1,nBatchASmall
      BATCHGABSmall(IA,IB) = BATCHGABBig(sA+IA,IBB)
   ENDDO
ENDDO 
!$OMP END PARALLEL DO 
end subroutine ExtractBatchGabFromFullGab

subroutine makeComCMO34(ComCMO34,nCMO3,nCMO4,nDimC,nDimD,&
     & CMO3,CMO4)
  implicit none
  integer,intent(in) :: nCMO3,nCMO4,nDimC,nDimD
  real(realk),intent(in) :: CMO3(nDimC,nCMO3)
  real(realk),intent(in) :: CMO4(nDimD,nCMO4)
  real(realk),intent(inout) :: ComCMO34(nDimC,nDimD,nCMO3,nCMO4)
  !
  integer :: IC,ID,J,B
  real(realk) :: TMP 
  !$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &   
  !$OMP PRIVATE(IC,ID,J,B,TMP) &
  !$OMP SHARED(ndimC,ndimD,nCMO3,nCMO4,CMO3,CMO4,ComCMO34) 
  DO J = 1,nCMO4
     DO B = 1,nCMO3
        DO ID = 1,nDimD
           TMP = CMO4(ID,J)
           DO IC = 1,nDimC
              ComCMO34(IC,ID,B,J) = CMO3(IC,B)*TMP
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
end subroutine makeComCMO34

subroutine FormActiveCmat(nDimB,nOcc,CmatB,nAtomsB,startOrbitalB,nOrbB,&
     & Cfull,nd1)
  implicit none
  integer,intent(in) :: nDimB,nAtomsB,nOrbB,nd1,nOcc
  integer,intent(in) :: startOrbitalB(nAtomsB)
  real(realk),intent(in) :: Cfull(nd1,nOcc)
  real(realk),intent(inout) :: CmatB(nDimB,nOcc)
  !local variables
  integer :: IatomB,IorbB,startLocB,startB,iOcc
  do IOcc= 1,nOcc
     startLocB = 0
     do IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO IorbB = 1,nOrbB
           CmatB(startLocB + IorbB,iOcc) = Cfull(startB+IorbB,iOcc)
        ENDDO
        startLocB = startLocB + nOrbB          
     ENDDO
  enddo
end subroutine FormActiveCmat

END MODULE IchorEriToolsmod
