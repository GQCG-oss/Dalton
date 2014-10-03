PROGRAM TUV
  use math
  use stringsMODULE
  implicit none
  logical :: nPrimLast
  integer,pointer :: TUVINDEX(:,:,:),TUVINDEXP(:,:,:)
  integer :: JMAX,J,JMAX1,JMAXP
  logical,pointer :: Enoscreen(:,:),EnoscreenS(:,:),zero(:)
  integer :: ijk1,ijk2,ijkcart,ijk,ijkcart1,ijkcart2,nTUV,ijkP
  integer :: iTUV,ilmP
  real(realk),pointer :: SCMAT1(:,:),SCMAT2(:,:),Spherical(:,:)
  logical :: sph1,sph2,sphericalTrans,sphericalGTO,newparam,ELSESTATEMENT
  integer :: LUMAIN,LUMOD1,LUMOD2,LUMOD3,LUMOD4,iparam,iparam2,nparam
  integer :: JMAX2,JMAX3,JMAX4,ituvP,jp,tp,up,vp,ntuvP,l1,l2,l12
  integer :: ip1,jp1,kp1,p1,ip2,jp2,kp2,p2,ijkcartP
  real(realk),pointer :: uniqeparam(:)
  character(len=15),pointer :: uniqeparamNAME(:)
  character(len=3) :: ARCSTRING
  integer :: GPUrun,nString
  logical :: DoOpenMP,DoOpenACC,CPU

  DO AngmomA = 0,2
   DO AngmomB = 0,2
    DO AngmomC = 0,2
     DO AngmomD = 0,2

      nOrbCompA = 2*AngmomA+1
      nOrbCompB = 2*AngmomB+1
      nOrbCompC = 2*AngmomC+1
      nOrbCompD = 2*AngmomD+1

      WRITE(ILUMOD,'(A)')''

      

      IF(Gen)THEN
         !(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
      ELSEIF(SegQ)THEN
         !(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContP,nPasses)
      ELSEIF(SegP)THEN
         !(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nPasses)
      ELSEIF(Seg.OR.Seg1Prim)THEN
         !(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nPasses)
      ENDIF
      ! =>   (nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD) 


      IF(Gen.OR.SegP)THEN
         WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE() PRIVATE(IPass,iContD,iContC,iAngD,iAngC)'
      ELSE

      ENDIF
      WRITE(ILUMOD,'(A)')'  DO IPass = 1,nPasses'
      IF(Gen.OR.SegP)THEN
      WRITE(ILUMOD,'(A)')'   DO iContD = 1,nContD'
      WRITE(ILUMOD,'(A)')'    DO iContC = 1,nContC'
      ENDIF
      IF(nOrbCompD.GT.1)THEN
      WRITE(ILUMOD,'(A)')'     DO iAngD = 1,',nOrbCompD
      ENDIF
      IF(nOrbCompC.GT.1)THEN
      WRITE(ILUMOD,'(A)')'      DO iAngC = 1,',nOrbCompC
      ENDIF
      WRITE(ILUMOD,'(A)')'       IatomB = IatomBPass(IPass)'
      WRITE(ILUMOD,'(A)')'       IatomA = IatomAPass(IPass)'
      IF(Gen.OR.SegP)THEN
      WRITE(ILUMOD,'(A)')'       iContQ = iContC + (iContD-1)*nContC'
      IF(nOrbCompD.GT.1)THEN
      WRITE(ILUMOD,'(A)')'       I4 = iAngD + (iContD-1)*',nOrbCompD
      ELSE
      WRITE(ILUMOD,'(A)')'       I4 = iContD'
      ENDIF
      IF(nOrbCompC.GT.1)THEN
      WRITE(ILUMOD,'(A)')'       I3 = iAngC + (iContC-1)*',nOrbCompC
      ELSE
      WRITE(ILUMOD,'(A)')'       I3 = iContC'
      ENDIF
      ENDIF !ENDIF(Gen.OR.SegP)THEN
      IF(Gen.OR.SegQ)THEN
      WRITE(ILUMOD,'(A)')'       DO iContB = 1,nContB'
      WRITE(ILUMOD,'(A)')'        offsetB = (iContB-1)*',nOrbCompB
      WRITE(ILUMOD,'(A)')'        DO iContA = 1,nContA'
      WRITE(ILUMOD,'(A)')'         iContP = iContA+(iContB-1)*nContA'
      WRITE(ILUMOD,'(A)')'         offsetA = (iContA-1)*',nOrbCompA
      ENDIF
      IF(nOrbCompA.GT.5.AND.nOrbCompB.GT.5)THEN
      WRITE(ILUMOD,'(A)')'         DO iAngB = 1,',nOrbCompB
      WRITE(ILUMOD,'(A)')'          DO iAngA = 1,',nOrbCompA
      IF(Gen)THEN
         WRITE(ILUMOD,'(A)')'           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&'
         WRITE(ILUMOD,'(A)')'                & LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)'
      IF(SegQ)THEN
         WRITE(ILUMOD,'(A)')'           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,iAngC,iAngD) =&'
         WRITE(ILUMOD,'(A)')'                & LP1(iAngA,iAngB,iAngC,iAngD,iContP,IPass)'
      IF(SegP)THEN
         WRITE(ILUMOD,'(A)')'           LP2(iAngA,iatomA,iAngB,iatomB,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,IPass)'
      IF(Seg.OR.Seg1Prim)THEN
         WRITE(ILUMOD,'(A)')'           LP2(iAngA,iatomA,iAngB,iatomB,iAngC,iAngD) = LP1(iAngA,iAngB,iAngC,iAngD,IPass)'
      ENDIF
      WRITE(ILUMOD,'(A)')'          ENDDO'
      WRITE(ILUMOD,'(A)')'         ENDDO'
      ELSEIF(nOrbCompA.GT.5)THEN
         !UNROLL nOrbCompB
      ELSEIF(nOrbCompB.GT.5)THEN
         !UNROLL nOrbCompA
      ELSE
         !UNROLL nOrbCompA AND UNROLL nOrbCompB
      ENDIF
      WRITE(ILUMOD,'(A)')'        ENDDO'
      WRITE(ILUMOD,'(A)')'       ENDDO'
      WRITE(ILUMOD,'(A)')'      ENDDO'
      WRITE(ILUMOD,'(A)')'     ENDDO'
      WRITE(ILUMOD,'(A)')'    ENDDO'
      WRITE(ILUMOD,'(A)')'   ENDDO'
      WRITE(ILUMOD,'(A)')'  ENDDO'


















              
              !(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)

              !(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD) 


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
     call TriDistributeCPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
          & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
          & MaxPasses,IatomAPass,iatomBPass,nPasses)
  ELSE
     call DistributeCPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
          & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
          & MaxPasses,IatomAPass,iatomBPass,nPasses)
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
     call TriDistributeGPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
          & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
          & MaxPasses,IatomAPass,iatomBPass,nPasses)
  ELSE
     call DistributeGPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
          & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
          & MaxPasses,IatomAPass,iatomBPass,nPasses)
  ENDIF
end subroutine MainTriDistributetoLocalIntPass2GPU

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





END PROGRAM

!contractecoeff_gen
