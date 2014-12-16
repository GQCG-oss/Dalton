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
  logical :: sph1,sph2,sphericalTrans,sphericalGTO,newparam,ELSESTATEMENT
  integer :: LUMAIN,LUMOD1,LUMOD2,LUMOD3,LUMOD4,iparam,iparam2,nparam
  integer :: JMAX2,JMAX3,JMAX4,ituvP,jp,tp,up,vp,ntuvP,l1,l2,l12
  integer :: ip1,jp1,kp1,p1,ip2,jp2,kp2,p2,ijkcartP
  character(len=15),pointer :: uniqeparamNAME(:)
  integer :: GPUrun,nString,I,LUFILE,iAngA,iAngB,iseg
  integer :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD
  integer :: AngmomA,AngmomB,ILUMOD,iseglabel,AngmomC,AngmomD
  logical :: DoOpenMP,DoOpenACC,CPU
  logical :: Gen,SegP,SegQ,Seg,UNROLLA,UNROLLB
  Character(len=48) :: FileName    
  Character(len=8)  :: SegLabel
  character(len=3) :: ARCSTRING
 ARCSTRING = 'CPU'
 DO iseg = 1,4
  Gen = .FALSE.; SegQ=.FALSE.; SegP=.FALSE.;Seg=.FALSE.
  IF(iseg.EQ.1)THEN
     Gen = .TRUE.      ; SegLabel = 'Gen     '; iSegLabel = 3
  ELSEIF(iseg.EQ.2)THEN
     SegQ = .TRUE.     ; SegLabel = 'SegQ    '; iSegLabel = 4
  ELSEIF(iseg.EQ.3)THEN
     SegP = .TRUE.     ; SegLabel = 'SegP    '; iSegLabel = 4
  ELSEIF(iseg.EQ.4)THEN
     Seg = .TRUE.      ; SegLabel = 'Seg     '; iSegLabel = 3
  ENDIF
  
  DO I =1,48
     FileName(I:I) = ' '
  ENDDO
  WRITE(FileName,'(3A)')'AGC_'//ARCSTRING//'_Distribute',SegLabel(1:iSegLabel),'.F90'
  print*,'FileName:',FileName
  LUFILE = 12
  open(unit = LUFILE, file=TRIM(FileName),status="unknown")
     
  ILUMOD = LUFILE
  WRITE(ILUMOD,'(2A)')'MODULE AGC_Distribute'//ARCSTRING,SegLabel(1:iSegLabel)
  WRITE(LUFILE,'(A)')' use IchorPrecisionMod'
  WRITE(LUFILE,'(A)')'  '
  WRITE(LUFILE,'(A)')' CONTAINS'
  DO AngmomA = 0,2
   DO AngmomB = 0,2
    DO AngmomC = 0,1
     DO AngmomD = 0,1

      nOrbCompA = 2*AngmomA+1
      nOrbCompB = 2*AngmomB+1
      nOrbCompC = 2*AngmomC+1
      nOrbCompD = 2*AngmomD+1

      IF(AngmomC.EQ.0.AND.AngmomD.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'subroutine DistributeCPUA',AngmomA,'B',AngmomB,'C0D0',SegLabel(1:iSegLabel),'(LP1,LP2,IatomAPass,IatomBPass,&'
      ELSEIF(AngmomC.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'subroutine DistributeCPUA',AngmomA,'B',AngmomB,'C0DX',SegLabel(1:iSegLabel),'(LP1,LP2,IatomAPass,IatomBPass,&'
      ELSEIF(AngmomD.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'subroutine DistributeCPUA',AngmomA,'B',AngmomB,'CXD0',SegLabel(1:iSegLabel),'(LP1,LP2,IatomAPass,IatomBPass,&'
      ELSE
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'subroutine DistributeCPUA',AngmomA,'B',AngmomB,'CXDX',SegLabel(1:iSegLabel),'(LP1,LP2,IatomAPass,IatomBPass,&'      
      ENDIF
      WRITE(ILUMOD,'(A)')'                    & nContA,nContB,nContC,nContD,nContQ,nContP,nPasses,&'
      WRITE(ILUMOD,'(A)')'                    & nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,MaxPasses,nOrbCompC,nOrbCompD)'
      WRITE(ILUMOD,'(A)')' implicit none'
      WRITE(ILUMOD,'(A)')' integer,intent(in) :: nContA,nContB,nContC,nContD,nOrbCompC'
      WRITE(ILUMOD,'(A)')' integer,intent(in) :: nAtomsB,nOrbC,nOrbD,MaxPasses,nOrbCompD'
      WRITE(ILUMOD,'(A)')' integer,intent(in) :: nContQ,nContP,nPasses,nOrbA,nAtomsA,nOrbB'
      WRITE(ILUMOD,'(A)')' integer,intent(in) :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)'

      !LP1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC,nContD,nContA,nContB,nPasses)
      call initString(1)
      call AddToString('real(realk),intent(in) :: LP1(')
      IF(nOrbCompA.GT.1)THEN
         call AddToString(nOrbCompA)
         call AddToString(',')
      ENDIF
      IF(nOrbCompB.GT.1)THEN
         call AddToString(nOrbCompB)
         call AddToString(',')
      ENDIF
      IF(nOrbCompC.GT.1)THEN
         call AddToString('nOrbCompC')
         call AddToString(',')
      ENDIF
      IF(nOrbCompD.GT.1)THEN
         call AddToString('nOrbCompD')
         call AddToString(',')
      ENDIF
      IF(Gen)THEN
         call AddToString('nContQ,nContP,nPasses)')
      ELSEIF(SegQ)THEN
         call AddToString('nContP,nPasses)')
      ELSEIF(SegP)THEN
         call AddToString('nContQ,nPasses)')
      ELSEIF(Seg)THEN
         call AddToString('nPasses)')
      ENDIF
      call writeString(ILUMOD)

      !LP2(nOrbCompA,nContA,nAtomsA,nOrbCompB,nContB,nAtomsB,nOrbCompC,nContC,nOrbCompD,nContD)
      call initString(1)
      call AddToString('real(realk),intent(inout) :: LP2(')
      IF(Seg.OR.SegP)THEN !nContA = 1
         IF(nOrbCompA.GT.1)THEN
            call AddToString(nOrbCompA)
            call AddToString(',')
         ENDIF
      ELSE
         IF(nOrbCompA.GT.1)THEN
            call AddToString('nOrbA,')
         ELSE
            call AddToString('nContA,')
         ENDIF
      ENDIF
      call AddToString('nAtomsA,')
      IF(Seg.OR.SegP)THEN !nContB = 1
         IF(nOrbCompB.GT.1)THEN
            call AddToString(nOrbCompB)
            call AddToString(',')
         ENDIF
      ELSE
         IF(nOrbCompB.GT.1)THEN
            call AddToString('nOrbB,')
         ELSE
            call AddToString('nContB,')
         ENDIF
      ENDIF
      call AddToString('nAtomsB')

      IF(Seg.OR.SegQ)THEN !nContC = 1, nContD = 1
         IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
            call AddToString(',nOrbCompC,nOrbCompD)')
         ELSEIF(nOrbCompC.GT.1)THEN
            call AddToString(',nOrbCompC)')
         ELSEIF(nOrbCompD.GT.1)THEN
            call AddToString(',nOrbCompD)')
         ELSE
            call AddToString(')')
         ENDIF
      ELSE
         IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
            call AddToString(',nOrbC,nOrbD)')
         ELSEIF(nOrbCompC.GT.1)THEN
            call AddToString(',nOrbC,nContD)')
         ELSEIF(nOrbCompD.GT.1)THEN
            call AddToString(',nContC,nOrbD)')
         ELSE
            call AddToString(',nContC,nContD)')
         ENDIF
      ENDIF
      call writeString(ILUMOD)
      WRITE(ILUMOD,'(A)')' !local variables'
      WRITE(ILUMOD,'(A)')' integer :: offsetB,offsetA,IPass,iAngA,iAngB,iAtomA,iAtomB,iContP'
      WRITE(ILUMOD,'(A)')' integer :: iContD,iContC,iAngD,iAngC,iContB,iContA,iContQ,I4,I3'
      IF(Gen)THEN !collapse pass and contc,cont or angD,angC
         WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iContD,iContC,iAngD,iAngC,iContB,iContA,&'
         WRITE(ILUMOD,'(A)')'!$OMP iContQ,I4,I3,offsetB,offsetA,&'
      ELSEIF(SegP)THEN!collapse pass and contc,cont or angD,angC
         WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iContD,iContC,iAngD,iAngC,iContQ,I4,I3,&'
      ELSEIF(SegQ)THEN!collapse pass and angD,angC
         IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iAngD,iAngC,iContB,iContA,offsetB,offsetA,&'
         ELSEIF(nOrbCompC.GT.1.OR.nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(2) PRIVATE(iAngD,iAngC,iContB,iContA,offsetB,offsetA,&'
         ELSE
            WRITE(ILUMOD,'(A)')'!$OMP DO PRIVATE(iAngD,iAngC,iContB,iContA,offsetB,offsetA,&'
         ENDIF
      ELSE !Seg Seg1Prim
         IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iAngD,iAngC,&'
         ELSEIF(nOrbCompC.GT.1.OR.nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(2) PRIVATE(iAngD,iAngC,&'
         ELSE
            WRITE(ILUMOD,'(A)')'!$OMP DO PRIVATE(iAngD,iAngC,&'
         ENDIF
      ENDIF
      IF(nOrbCompA.GT.5.AND.nOrbCompB.GT.5)THEN
         WRITE(ILUMOD,'(A)')'!$OMP IPass,iAngA,iAngB,iAtomA,iAtomB)'
      ELSEIF(nOrbCompA.GT.5)THEN
         WRITE(ILUMOD,'(A)')'!$OMP IPass,iAngA,iAtomA,iAtomB)'
      ELSEIF(nOrbCompB.GT.5)THEN
         WRITE(ILUMOD,'(A)')'!$OMP IPass,iAngB,iAtomA,iAtomB)'
      ELSE
         WRITE(ILUMOD,'(A)')'!$OMP IPass,iAtomA,iAtomB)'
      ENDIF
   
      WRITE(ILUMOD,'(A)')'  DO IPass = 1,nPasses'
      IF(Gen.OR.SegP)THEN
         WRITE(ILUMOD,'(A)')'   DO iContD = 1,nContD'
         WRITE(ILUMOD,'(A)')'    DO iContC = 1,nContC'
      ENDIF
      IF(nOrbCompD.GT.1)THEN
         WRITE(ILUMOD,'(A)')'     DO iAngD = 1,nOrbCompD'
      ENDIF
      IF(nOrbCompC.GT.1)THEN
         WRITE(ILUMOD,'(A)')'      DO iAngC = 1,nOrbCompC'
      ENDIF
      WRITE(ILUMOD,'(A)')'       IatomB = IatomBPass(IPass)'
      WRITE(ILUMOD,'(A)')'       IatomA = IatomAPass(IPass)'
      IF(Gen.OR.SegP)THEN
         WRITE(ILUMOD,'(A)')'       iContQ = iContC + (iContD-1)*nContC'
         IF(nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'       I4 = iAngD + (iContD-1)*nOrbCompD'
!         ELSE
!            WRITE(ILUMOD,'(A)')'       I4 = iContD'
         ENDIF
         IF(nOrbCompC.GT.1)THEN
            WRITE(ILUMOD,'(A)')'       I3 = iAngC + (iContC-1)*nOrbCompC'
!         ELSE
!            WRITE(ILUMOD,'(A)')'       I3 = iContC'
         ENDIF
      ENDIF !ENDIF(Gen.OR.SegP)THEN
      IF(Gen.OR.SegQ)THEN
         WRITE(ILUMOD,'(A)')'       DO iContB = 1,nContB'
         IF(nOrbCompB.GT.1)THEN
            WRITE(ILUMOD,'(A,I1)')'        offsetB = (iContB-1)*',nOrbCompB
         ENDIF
         WRITE(ILUMOD,'(A)')'        DO iContA = 1,nContA'
         WRITE(ILUMOD,'(A)')'         iContP = iContA+(iContB-1)*nContA'
         IF(nOrbCompA.GT.1)THEN
            WRITE(ILUMOD,'(A,I1)')'         offsetA = (iContA-1)*',nOrbCompA
         ENDIF
      ENDIF
      IF(nOrbCompA.GT.7.AND.nOrbCompB.GT.1)THEN
         UNROLLA = .FALSE.
         UNROLLB = .FALSE.
         WRITE(ILUMOD,'(A,I1)')'         DO iAngB = 1,',nOrbCompB
         WRITE(ILUMOD,'(A,I1)')'          DO iAngA = 1,',nOrbCompA
         CAll BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
         WRITE(ILUMOD,'(A)')'          ENDDO'
         WRITE(ILUMOD,'(A)')'         ENDDO'
      ELSEIF(nOrbCompA.GT.7)THEN
         UNROLLA = .FALSE.
         UNROLLB = .TRUE.
         WRITE(ILUMOD,'(A,I1)')'          DO iAngA = 1,',nOrbCompA
         DO iAngB = 1,nOrbCompB
            CAll BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
         ENDDO
         WRITE(ILUMOD,'(A)')'          ENDDO'
      ELSEIF(nOrbCompB.GT.1)THEN
         UNROLLA = .TRUE.
         UNROLLB = .FALSE.
         WRITE(ILUMOD,'(A,I1)')'          DO iAngB = 1,',nOrbCompB
         DO iAngA = 1,nOrbCompA
            CAll BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
         ENDDO
         WRITE(ILUMOD,'(A)')'          ENDDO'
      ELSE
         UNROLLA = .TRUE.
         UNROLLB = .TRUE.
         DO iAngB = 1,nOrbCompB
            DO iAngA = 1,nOrbCompA
               CAll BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
            ENDDO
         ENDDO
      ENDIF
      IF(Gen.OR.SegQ)THEN !ContA and ContB
         WRITE(ILUMOD,'(A)')'        ENDDO'
         WRITE(ILUMOD,'(A)')'       ENDDO'
      ENDIF
      !AngC and AngD
      IF(nOrbCompD.GT.1)THEN
         WRITE(ILUMOD,'(A)')'      ENDDO'
      ENDIF
      IF(nOrbCompC.GT.1)THEN
         WRITE(ILUMOD,'(A)')'     ENDDO'
      ENDIF
      IF(Gen.OR.SegP)THEN !ContC and ContD 
         WRITE(ILUMOD,'(A)')'    ENDDO'
         WRITE(ILUMOD,'(A)')'   ENDDO'
      ENDIF
      WRITE(ILUMOD,'(A)')'  ENDDO'
      WRITE(ILUMOD,'(A)')'!$OMP END DO'
      
      IF(AngmomC.EQ.0.AND.AngmomD.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'end subroutine DistributeCPUA',AngmomA,'B',AngmomB,'C0D0',SegLabel(1:iSegLabel)
      ELSEIF(AngmomC.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'end subroutine DistributeCPUA',AngmomA,'B',AngmomB,'C0DX',SegLabel(1:iSegLabel)
      ELSEIF(AngmomD.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'end subroutine DistributeCPUA',AngmomA,'B',AngmomB,'CXD0',SegLabel(1:iSegLabel)
      ELSE
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'end subroutine DistributeCPUA',AngmomA,'B',AngmomB,'CXDX',SegLabel(1:iSegLabel)
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  WRITE(ILUMOD,'(2A)')'END MODULE AGC_Distribute'//ARCSTRING,SegLabel(1:iSegLabel)
  close(unit = LUFILE)


 !=======================================================================================================================


  DO I =1,48
     FileName(I:I) = ' '
  ENDDO
  WRITE(FileName,'(3A)')'AGC_'//ARCSTRING//'_TriDistribute',SegLabel(1:iSegLabel),'.F90'
  print*,'FileName:',FileName
  LUFILE = 12
  open(unit = LUFILE, file=TRIM(FileName),status="unknown")
     
  ILUMOD = LUFILE
  WRITE(ILUMOD,'(2A)')'MODULE AGC_TriDistribute'//ARCSTRING,SegLabel(1:iSegLabel)
  WRITE(LUFILE,'(A)')' use IchorPrecisionMod'
  WRITE(LUFILE,'(A)')'  '
  WRITE(LUFILE,'(A)')' CONTAINS'
  DO AngmomA = 0,2
   DO AngmomB = 0,2
    DO AngmomC = 0,1
     DO AngmomD = 0,1

      nOrbCompA = 2*AngmomA+1
      nOrbCompB = 2*AngmomB+1
      nOrbCompC = 2*AngmomC+1
      nOrbCompD = 2*AngmomD+1

      IF(AngmomC.EQ.0.AND.AngmomD.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'subroutine TriDistributeCPUA',AngmomA,'B',AngmomB,'C0D0',SegLabel(1:iSegLabel),'(LP1,LP2,IatomAPass,IatomBPass,&'
      ELSEIF(AngmomC.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'subroutine TriDistributeCPUA',AngmomA,'B',AngmomB,'C0DX',SegLabel(1:iSegLabel),'(LP1,LP2,IatomAPass,IatomBPass,&'
      ELSEIF(AngmomD.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'subroutine TriDistributeCPUA',AngmomA,'B',AngmomB,'CXD0',SegLabel(1:iSegLabel),'(LP1,LP2,IatomAPass,IatomBPass,&'
      ELSE
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'subroutine TriDistributeCPUA',AngmomA,'B',AngmomB,'CXDX',SegLabel(1:iSegLabel),'(LP1,LP2,IatomAPass,IatomBPass,&'      
      ENDIF
      WRITE(ILUMOD,'(A)')'                    & nContA,nContB,nContC,nContD,nContQ,nContP,nPasses,&'
      WRITE(ILUMOD,'(A)')'                    & nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,MaxPasses,nOrbCompC,nOrbCompD)'
      WRITE(ILUMOD,'(A)')' implicit none'
      WRITE(ILUMOD,'(A)')' integer,intent(in) :: nContA,nContB,nContC,nContD,nOrbCompC'
      WRITE(ILUMOD,'(A)')' integer,intent(in) :: nAtomsB,nOrbC,nOrbD,MaxPasses,nOrbCompD'
      WRITE(ILUMOD,'(A)')' integer,intent(in) :: nContQ,nContP,nPasses,nOrbA,nAtomsA,nOrbB'
      WRITE(ILUMOD,'(A)')' integer,intent(in) :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)'

      !LP1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC,nContD,nContA,nContB,nPasses)
      call initString(1)
      call AddToString('real(realk),intent(in) :: LP1(')
      IF(nOrbCompA.GT.1)THEN
         call AddToString(nOrbCompA)
         call AddToString(',')
      ENDIF
      IF(nOrbCompB.GT.1)THEN
         call AddToString(nOrbCompB)
         call AddToString(',')
      ENDIF
      IF(nOrbCompC.GT.1)THEN
         call AddToString('nOrbCompC')
         call AddToString(',')
      ENDIF
      IF(nOrbCompD.GT.1)THEN
         call AddToString('nOrbCompD')
         call AddToString(',')
      ENDIF
      IF(Gen)THEN
         call AddToString('nContQ,nContP,nPasses)')
      ELSEIF(SegQ)THEN
         call AddToString('nContP,nPasses)')
      ELSEIF(SegP)THEN
         call AddToString('nContQ,nPasses)')
      ELSEIF(Seg)THEN
         call AddToString('nPasses)')
      ENDIF
      call writeString(ILUMOD)

      !LP2(nOrbCompA,nContA,nAtomsA,nOrbCompB,nContB,nAtomsB,nOrbCompC,nContC,nOrbCompD,nContD)
      call initString(1)
      call AddToString('real(realk),intent(inout) :: LP2(')
      IF(Seg.OR.SegP)THEN !nContA = 1
         IF(nOrbCompA.GT.1)THEN
            call AddToString(nOrbCompA)
            call AddToString(',')
         ENDIF
      ELSE
         IF(nOrbCompA.GT.1)THEN
            call AddToString('nOrbA,')
         ELSE
            call AddToString('nContA,')
         ENDIF
      ENDIF
      call AddToString('nAtomsA,')
      IF(Seg.OR.SegP)THEN !nContB = 1
         IF(nOrbCompB.GT.1)THEN
            call AddToString(nOrbCompB)
            call AddToString(',')
         ENDIF
      ELSE
         IF(nOrbCompB.GT.1)THEN
            call AddToString('nOrbB,')
         ELSE
            call AddToString('nContB,')
         ENDIF
      ENDIF
      call AddToString('nAtomsB')

      IF(Seg.OR.SegQ)THEN !nContC = 1, nContD = 1
         IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
            call AddToString(',nOrbCompC,nOrbCompD)')
         ELSEIF(nOrbCompC.GT.1)THEN
            call AddToString(',nOrbCompC)')
         ELSEIF(nOrbCompD.GT.1)THEN
            call AddToString(',nOrbCompD)')
         ELSE
            call AddToString(')')
         ENDIF
      ELSE
         IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
            call AddToString(',nOrbC,nOrbD)')
         ELSEIF(nOrbCompC.GT.1)THEN
            call AddToString(',nOrbC,nContD)')
         ELSEIF(nOrbCompD.GT.1)THEN
            call AddToString(',nContC,nOrbD)')
         ELSE
            call AddToString(',nContC,nContD)')
         ENDIF
      ENDIF
      call writeString(ILUMOD)
      WRITE(ILUMOD,'(A)')' !local variables'
      WRITE(ILUMOD,'(A)')' integer :: offsetB,offsetA,IPass,iAngA,iAngB,iAtomA,iAtomB,iContP'
      WRITE(ILUMOD,'(A)')' integer :: iContD,iContC,iAngD,iAngC,iContB,iContA,iContQ,I4,I3'
      IF(Gen)THEN !collapse pass and contc,cont or angD,angC
         WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iContD,iContC,iAngD,iAngC,iContB,iContA,&'
         WRITE(ILUMOD,'(A)')'!$OMP iContQ,I4,I3,offsetB,offsetA,&'
      ELSEIF(SegP)THEN!collapse pass and contc,cont or angD,angC
         WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iContD,iContC,iAngD,iAngC,iContQ,I4,I3,&'
      ELSEIF(SegQ)THEN!collapse pass and angD,angC
         IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iAngD,iAngC,iContB,iContA,offsetB,offsetA,&'
         ELSEIF(nOrbCompC.GT.1.OR.nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(2) PRIVATE(iAngD,iAngC,iContB,iContA,offsetB,offsetA,&'
         ELSE
            WRITE(ILUMOD,'(A)')'!$OMP DO PRIVATE(iAngD,iAngC,iContB,iContA,offsetB,offsetA,&'
         ENDIF
      ELSE !Seg Seg1Prim
         IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iAngD,iAngC,&'
         ELSEIF(nOrbCompC.GT.1.OR.nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'!$OMP DO COLLAPSE(2) PRIVATE(iAngD,iAngC,&'
         ELSE
            WRITE(ILUMOD,'(A)')'!$OMP DO PRIVATE(iAngD,iAngC,&'
         ENDIF
      ENDIF
      IF(nOrbCompA.GT.5.AND.nOrbCompB.GT.5)THEN
         WRITE(ILUMOD,'(A)')'!$OMP IPass,iAngA,iAngB,iAtomA,iAtomB)'
      ELSEIF(nOrbCompA.GT.5)THEN
         WRITE(ILUMOD,'(A)')'!$OMP IPass,iAngA,iAtomA,iAtomB)'
      ELSEIF(nOrbCompB.GT.5)THEN
         WRITE(ILUMOD,'(A)')'!$OMP IPass,iAngB,iAtomA,iAtomB)'
      ELSE
         WRITE(ILUMOD,'(A)')'!$OMP IPass,iAtomA,iAtomB)'
      ENDIF
   
      WRITE(ILUMOD,'(A)')'  DO IPass = 1,nPasses'
      IF(Gen.OR.SegP)THEN
         WRITE(ILUMOD,'(A)')'   DO iContD = 1,nContD'
         WRITE(ILUMOD,'(A)')'    DO iContC = 1,nContC'
      ENDIF
      IF(nOrbCompD.GT.1)THEN
         WRITE(ILUMOD,'(A)')'     DO iAngD = 1,nOrbCompD'
      ENDIF
      IF(nOrbCompC.GT.1)THEN
         WRITE(ILUMOD,'(A)')'      DO iAngC = 1,nOrbCompC'
      ENDIF
      WRITE(ILUMOD,'(A)')'       IatomB = IatomBPass(IPass)'
      WRITE(ILUMOD,'(A)')'       IatomA = IatomAPass(IPass)'
      IF(Gen.OR.SegP)THEN
         WRITE(ILUMOD,'(A)')'       iContQ = iContC + (iContD-1)*nContC'
         IF(nOrbCompD.GT.1)THEN
            WRITE(ILUMOD,'(A)')'       I4 = iAngD + (iContD-1)*nOrbCompD'
!         ELSE
!            WRITE(ILUMOD,'(A)')'       I4 = iContD'
         ENDIF
         IF(nOrbCompC.GT.1)THEN
            WRITE(ILUMOD,'(A)')'       I3 = iAngC + (iContC-1)*nOrbCompC'
!         ELSE
!            WRITE(ILUMOD,'(A)')'       I3 = iContC'
         ENDIF
      ENDIF !ENDIF(Gen.OR.SegP)THEN
      IF(Gen.OR.SegQ)THEN
         WRITE(ILUMOD,'(A)')'       DO iContB = 1,nContB'
         IF(nOrbCompB.GT.1)THEN
            WRITE(ILUMOD,'(A,I1)')'        offsetB = (iContB-1)*',nOrbCompB
         ENDIF
         WRITE(ILUMOD,'(A)')'        DO iContA = 1,nContA'
         WRITE(ILUMOD,'(A)')'         iContP = iContA+(iContB-1)*nContA'
         IF(nOrbCompA.GT.1)THEN
            WRITE(ILUMOD,'(A,I1)')'         offsetA = (iContA-1)*',nOrbCompA
         ENDIF
      ENDIF
      IF(nOrbCompA.GT.7.AND.nOrbCompB.GT.1)THEN
         UNROLLA = .FALSE.
         UNROLLB = .FALSE.
         WRITE(ILUMOD,'(A,I1)')'         DO iAngB = 1,',nOrbCompB
         WRITE(ILUMOD,'(A,I1)')'          DO iAngA = 1,',nOrbCompA
         CAll BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
         WRITE(ILUMOD,'(A)')'          ENDDO'
         WRITE(ILUMOD,'(A)')'         ENDDO'
      ELSEIF(nOrbCompA.GT.7)THEN
         UNROLLA = .FALSE.
         UNROLLB = .TRUE.
         WRITE(ILUMOD,'(A,I1)')'          DO iAngA = 1,',nOrbCompA
         DO iAngB = 1,nOrbCompB
            CAll BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
         ENDDO
         WRITE(ILUMOD,'(A)')'          ENDDO'
      ELSEIF(nOrbCompB.GT.1)THEN
         UNROLLA = .TRUE.
         UNROLLB = .FALSE.
         WRITE(ILUMOD,'(A,I1)')'          DO iAngB = 1,',nOrbCompB
         DO iAngA = 1,nOrbCompA
            CAll BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
         ENDDO
         WRITE(ILUMOD,'(A)')'          ENDDO'
      ELSE
         UNROLLA = .TRUE.
         UNROLLB = .TRUE.
         DO iAngB = 1,nOrbCompB
            DO iAngA = 1,nOrbCompA
               CAll BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
            ENDDO
         ENDDO
      ENDIF
      IF(Gen.OR.SegQ)THEN !ContA and ContB
         WRITE(ILUMOD,'(A)')'        ENDDO'
         WRITE(ILUMOD,'(A)')'       ENDDO'
      ENDIF
      !AngC and AngD
      IF(nOrbCompD.GT.1)THEN
         WRITE(ILUMOD,'(A)')'      ENDDO'
      ENDIF
      IF(nOrbCompC.GT.1)THEN
         WRITE(ILUMOD,'(A)')'     ENDDO'
      ENDIF
      IF(Gen.OR.SegP)THEN !ContC and ContD 
         WRITE(ILUMOD,'(A)')'    ENDDO'
         WRITE(ILUMOD,'(A)')'   ENDDO'
      ENDIF
      WRITE(ILUMOD,'(A)')'  ENDDO'
      WRITE(ILUMOD,'(A)')'!$OMP END DO'
      
      IF(AngmomC.EQ.0.AND.AngmomD.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'end subroutine TriDistributeCPUA',AngmomA,'B',AngmomB,'C0D0',SegLabel(1:iSegLabel)
      ELSEIF(AngmomC.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'end subroutine TriDistributeCPUA',AngmomA,'B',AngmomB,'C0DX',SegLabel(1:iSegLabel)
      ELSEIF(AngmomD.EQ.0)THEN
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'end subroutine TriDistributeCPUA',AngmomA,'B',AngmomB,'CXD0',SegLabel(1:iSegLabel)
      ELSE
         WRITE(ILUMOD,'(A,I1,A,I1,A,A,A)')'end subroutine TriDistributeCPUA',AngmomA,'B',AngmomB,'CXDX',SegLabel(1:iSegLabel)
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  WRITE(ILUMOD,'(2A)')'END MODULE AGC_TriDistribute'//ARCSTRING,SegLabel(1:iSegLabel)
  close(unit = LUFILE)

 ENDDO
CONTAINS

  subroutine BuildLine(Gen,Seg,SegP,SegQ,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,UNROLLA,UNROLLB,IangA,IangB)
    implicit none
    integer,intent(in) :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,IangA,IangB
    logical,intent(in) :: Gen,Seg,SegP,SegQ,UNROLLA,UNROLLB
    !
    logical :: SPLITLINE
    SPLITLINE = .FALSE.
    call initString(11)
    call AddToString('LP2(')
    IF(SegP.OR.Seg)THEN !nContA = 1 => offsetA = 0
       IF(nOrbCompA.GT.1)THEN
          IF(UNROLLA)THEN
             call AddToString(IangA)
             call AddToString(',')
          ELSE
             call AddToString('iAngA,')
          ENDIF
       ENDIF
    ELSE
       IF(nOrbCompA.GT.1)THEN
          IF(UNROLLA)THEN
             call AddToString(IangA)
             call AddToString(' + offsetA,')
          ELSE
             call AddToString('iAngA + offsetA,')
          ENDIF
          SPLITLINE = .TRUE.
       ELSE
          call AddToString('iContA,')
       ENDIF
    ENDIF
    call AddToString('iAtomA,')
    
    IF(SegP.OR.Seg)THEN !nContB = 1 => offsetB = 0
       IF(nOrbCompB.GT.1)THEN
          IF(UNROLLB)THEN
             call AddToString(IangB)
             call AddToString(',')
          ELSE
             call AddToString('iAngB,')
          ENDIF
       ENDIF
    ELSE
       IF(nOrbCompB.GT.1)THEN
          IF(UNROLLB)THEN
             call AddToString(IangB)
             call AddToString(' + offsetB,')
          ELSE
             call AddToString('iAngB + offsetB,')
          ENDIF
          SPLITLINE = .TRUE.
       ELSE
          call AddToString('iContB,')
       ENDIF
    ENDIF
    call AddToString('iAtomB')
    IF(Gen.OR.SegP)THEN !nContC and nContD > 1
       IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
          call AddToString(',I3,I4) =')
       ELSEIF(nOrbCompC.GT.1)THEN
          call AddToString(',I3,iContD) =')
       ELSEIF(nOrbCompD.GT.1)THEN
          call AddToString(',iContC,I4) =')
       ELSE
          call AddToString(',iContC,iContD) =')
       ENDIF
    ELSE
       IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
          call AddToString(',iAngC,iAngD) =')
       ELSEIF(nOrbCompC.GT.1)THEN
          call AddToString(',iAngC) =')
       ELSEIF(nOrbCompD.GT.1)THEN
          call AddToString(',iAngD) =')
       ELSE
          call AddToString(') =')
       ENDIF
    ENDIF
    IF(SPLITLINE)THEN
       call AddToString(' &')
       call writeString(ILUMOD)
       call initString(16)
       call AddToString('& LP1(')
    ELSE
       call AddToString(' LP1(')
    ENDIF
    IF(nOrbCompA.GT.1)THEN
       IF(UNROLLA)THEN
          call AddToString(IangA)
          call AddToString(',')
       ELSE
          call AddToString('iAngA,')
       ENDIF
    ENDIF
    IF(nOrbCompB.GT.1)THEN
       IF(UNROLLB)THEN
          call AddToString(IangB)
          call AddToString(',')
       ELSE
          call AddToString('iAngB,')
       ENDIF
    ENDIF
    IF(nOrbCompC.GT.1)THEN
       call AddToString('iAngC,')
    ENDIF
    IF(nOrbCompD.GT.1)THEN
       call AddToString('iAngD,')
    ENDIF
    IF(Gen.OR.SegP)THEN !nContQ > 1
       call AddToString('iContQ,')
    ENDIF
    IF(Gen.OR.SegQ)THEN !nContP > 1
       call AddToString('iContP,')
    ENDIF
    call AddToString('iPass)')
    call writeString(ILUMOD)
  end subroutine BuildLine
end PROGRAM TUV




!!$
!!$         IF(Gen)THEN
!!$            WRITE(ILUMOD,'(A)')'           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&'
!!$            WRITE(ILUMOD,'(A)')'                & LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)'
!!$         ELSEIF(SegQ)THEN
!!$            WRITE(ILUMOD,'(A)')'           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,iAngC,iAngD) =&'
!!$            WRITE(ILUMOD,'(A)')'                & LP1(iAngA,iAngB,iAngC,iAngD,iContP,IPass)'
!!$         ELSEIF(SegP)THEN
!!$            IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$               WRITE(ILUMOD,'(A)')'           LP2(iAngA,iatomA,iAngB,iatomB,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,IPass)'
!!$            ELSEIF(nOrbCompC.GT.1)THEN
!!$               WRITE(ILUMOD,'(A)')'           LP2(iAngA,iatomA,iAngB,iatomB,I3,iContD) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,IPass)'
!!$            ELSEIF(nOrbCompD.GT.1)THEN
!!$               WRITE(ILUMOD,'(A)')'           LP2(iAngA,iatomA,iAngB,iatomB,iContC,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,IPass)'
!!$            ELSE
!!$               WRITE(ILUMOD,'(A)')'           LP2(iAngA,iatomA,iAngB,iatomB,iContC,iContD) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,IPass)'
!!$            ENDIF
!!$         ELSEIF(Seg)THEN
!!$            WRITE(ILUMOD,'(A)')'           LP2(iAngA,iatomA,iAngB,iatomB,iAngC,iAngD) = LP1(iAngA,iAngB,iAngC,iAngD,IPass)'
!!$         ENDIF
!!$         WRITE(ILUMOD,'(A)')'          ENDDO'
!!$         WRITE(ILUMOD,'(A)')'         ENDDO'
!!$      ELSEIF(nOrbCompA.GT.5)THEN
!!$         !UNROLL nOrbCompB
!!$         IF(nOrbCompB.EQ.1)THEN
!!$            WRITE(ILUMOD,'(A,I1)')'          DO iAngA = 1,',nOrbCompA
!!$            IF(Gen)THEN
!!$               WRITE(ILUMOD,'(A,I1,A)')'           LP2(iAngA + offsetA,iatomA,iContB,iatomB,I3,I4) =&'
!!$               WRITE(ILUMOD,'(A,I1,A)')'                & LP1(iAngA,iAngC,iAngD,iContQ,iContP,IPass)'
!!$            ELSEIF(SegQ)THEN
!!$               WRITE(ILUMOD,'(A,I1,A)')'           LP2(iAngA + offsetA,iatomA,iContB,iatomB,iAngC,iAngD) =&'
!!$               WRITE(ILUMOD,'(A,I1,A)')'                & LP1(iAngA,iAngC,iAngD,iContP,IPass)'
!!$            ELSEIF(SegP)THEN
!!$               IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,iContB,iatomB,I3,I4) = LP1(iAngA,iAngC,iAngD,iContQ,IPass)'
!!$               ELSEIF(nOrbCompC.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,iContB,iatomB,I3,iContD) = LP1(iAngA,iAngC,iAngD,iContQ,IPass)'
!!$               ELSEIF(nOrbCompD.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,iContB,iatomB,iContC,I4) = LP1(iAngA,iAngC,iAngD,iContQ,IPass)'
!!$               ELSE
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,iContB,iatomB,iContC,iContD) = LP1(iAngA,iAngC,iAngD,iContQ,IPass)'
!!$               ENDIF
!!$            ELSEIF(Seg)THEN
!!$               WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,iContB,iatomB,iAngC,iAngD) = LP1(iAngA,iAngC,iAngD,IPass)'
!!$            ENDIF
!!$            WRITE(ILUMOD,'(A)')'          ENDDO'
!!$         ELSE
!!$            WRITE(ILUMOD,'(A,I1)')'          DO iAngA = 1,',nOrbCompA
!!$            IF(Gen)THEN
!!$               WRITE(ILUMOD,'(A,I1,A)')'           LP2(iAngA + offsetA,iatomA,',iAngB,' + offsetB,iatomB,I3,I4) =&'
!!$               WRITE(ILUMOD,'(A,I1,A)')'                & LP1(iAngA,',iAngB,',iAngC,iAngD,iContQ,iContP,IPass)'
!!$            ELSEIF(SegQ)THEN
!!$               WRITE(ILUMOD,'(A,I1,A)')'           LP2(iAngA + offsetA,iatomA,',iAngB,' + offsetB,iatomB,iAngC,iAngD) =&'
!!$               WRITE(ILUMOD,'(A,I1,A)')'                & LP1(iAngA,',iAngB,',iAngC,iAngD,iContP,IPass)'
!!$            ELSEIF(SegP)THEN
!!$               IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,',iAngB,',iatomB,I3,I4) = LP1(iAngA,',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$               ELSEIF(nOrbCompC.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,',iAngB,',iatomB,I3,iContD) = LP1(iAngA,',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$               ELSEIF(nOrbCompD.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,',iAngB,',iatomB,iContC,I4) = LP1(iAngA,',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$               ELSE
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,',iAngB,',iatomB,iContC,iContD) = LP1(iAngA,',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$               ENDIF
!!$            ELSEIF(Seg)THEN
!!$               WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iAngA,iatomA,',iAngB,',iatomB,iAngC,iAngD) = LP1(iAngA,',iAngB,',iAngC,iAngD,IPass)'
!!$            ENDIF
!!$            WRITE(ILUMOD,'(A)')'          ENDDO'
!!$         ENDIF
!!$      ELSEIF(nOrbCompB.GT.5)THEN
!!$         !UNROLL nOrbCompA
!!$         IF(nOrbCompA.EQ.1)THEN
!!$            WRITE(ILUMOD,'(A,I1)')'         DO iAngB = 1,',nOrbCompB
!!$            IF(Gen)THEN
!!$               WRITE(ILUMOD,'(A,I1,A)')'           LP2(iContA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&'
!!$               WRITE(ILUMOD,'(A,I1,A)')'                & LP1(iAngB,iAngC,iAngD,iContQ,iContP,IPass)'
!!$            ELSEIF(SegQ)THEN
!!$               WRITE(ILUMOD,'(A,I1,A)')'           LP2(iContA,iatomA,iAngB + offsetB,iatomB,iAngC,iAngD) =&'
!!$               WRITE(ILUMOD,'(A,I1,A)')'                & LP1(iAngB,iAngC,iAngD,iContP,IPass)'
!!$            ELSEIF(SegP)THEN
!!$               IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iatomA,iAngB,iatomB,I3,I4) = LP1(iAngB,iAngC,iAngD,iContQ,IPass)'
!!$               ELSEIF(nOrbCompC.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iatomA,iAngB,iatomB,I3,iContD) = LP1(iAngB,iAngC,iAngD,iContQ,IPass)'
!!$               ELSEIF(nOrbCompD.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iatomA,iAngB,iatomB,iContC,I4) = LP1(iAngB,iAngC,iAngD,iContQ,IPass)'
!!$               ELSE
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iatomA,iAngB,iatomB,iContC,iContD) = LP1(iAngB,iAngC,iAngD,iContQ,IPass)'
!!$               ENDIF
!!$            ELSEIF(Seg)THEN
!!$               WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iatomA,iAngB,iatomB,iAngC,iAngD) = LP1(iAngB,iAngC,iAngD,IPass)'
!!$            ENDIF
!!$            WRITE(ILUMOD,'(A)')'         ENDDO'
!!$         ELSE
!!$            DO iAngA = 1,nOrbCompA
!!$               WRITE(ILUMOD,'(A,I1)')'         DO iAngB = 1,',nOrbCompB
!!$               IF(Gen)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A)')'           LP2(',iAngA,' + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&'
!!$                  WRITE(ILUMOD,'(A,I1,A)')'                & LP1(',iAngA,',iAngB,iAngC,iAngD,iContQ,iContP,IPass)'
!!$               ELSEIF(SegQ)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A)')'           LP2(',iAngA,' + offsetA,iatomA,iAngB + offsetB,iatomB,iAngC,iAngD) =&'
!!$                  WRITE(ILUMOD,'(A,I1,A)')'                & LP1(',iAngA,',iAngB,iAngC,iAngD,iContP,IPass)'
!!$               ELSEIF(SegP)THEN
!!$                  IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iAngB,iatomB,I3,I4) = LP1(',iAngA,',iAngB,iAngC,iAngD,iContQ,IPass)'
!!$                  ELSEIF(nOrbCompC.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iAngB,iatomB,I3,iContD) = LP1(',iAngA,',iAngB,iAngC,iAngD,iContQ,IPass)'
!!$                  ELSEIF(nOrbCompD.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iAngB,iatomB,iContC,I4) = LP1(',iAngA,',iAngB,iAngC,iAngD,iContQ,IPass)'
!!$                  ELSE
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iAngB,iatomB,iContC,iContD) = LP1(',iAngA,',iAngB,iAngC,iAngD,iContQ,IPass)'
!!$                  ENDIF
!!$               ELSEIF(Seg)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iAngB,iatomB,iAngC,iAngD) = LP1(',iAngA,',iAngB,iAngC,iAngD,IPass)'
!!$               ENDIF
!!$               WRITE(ILUMOD,'(A)')'         ENDDO'
!!$            ENDDO
!!$         ENDIF
!!$      ELSE
!!$         IF(nOrbCompA.EQ.1.AND.nOrbCompB.EQ.1)THEN
!!$            !UNROLL nOrbCompA AND UNROLL nOrbCompB
!!$            IF(Gen)THEN
!!$               WRITE(ILUMOD,'(A)')'           LP2(iContA,iatomA,iContB,iatomB,I3,I4) = LP1(iAngC,iAngD,iContQ,iContP,IPass)'
!!$            ELSEIF(SegQ)THEN
!!$               WRITE(ILUMOD,'(A)')'           LP2(iContA,iatomA,iContB,iatomB,iAngC,iAngD) = LP1(iAngC,iAngD,iContP,IPass)'
!!$            ELSEIF(SegP)THEN
!!$               IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A)')'           LP2(iatomA,iatomB,I3,I4) = LP1(iAngC,iAngD,iContQ,IPass)'
!!$               ELSEIF(nOrbCompC.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A)')'           LP2(iatomA,iatomB,I3,iContD) = LP1(iAngC,iAngD,iContQ,IPass)'
!!$               ELSEIF(nOrbCompD.GT.1)THEN
!!$                  WRITE(ILUMOD,'(A)')'           LP2(iatomA,iatomB,iContC,I4) = LP1(iAngC,iAngD,iContQ,IPass)'
!!$               ELSE
!!$                  WRITE(ILUMOD,'(A)')'           LP2(iatomA,iatomB,iContC,iContD) = LP1(iAngC,iAngD,iContQ,IPass)'
!!$               ENDIF
!!$            ELSEIF(Seg)THEN
!!$               WRITE(ILUMOD,'(A)')'           LP2(iatomA,iatomB,iAngC,iAngD) = LP1(iAngC,iAngD,IPass)'
!!$            ENDIF
!!$         ELSEIF(nOrbCompA.EQ.1)THEN
!!$            !UNROLL nOrbCompA AND UNROLL nOrbCompB
!!$            DO iAngB = 1,nOrbCompB
!!$               IF(Gen)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iContA,iatomA,',iAngB,' + offsetB,iatomB,I3,I4) =&'
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'                & LP1(',iAngB,',iAngC,iAngD,iContQ,iContP,IPass)'
!!$               ELSEIF(SegQ)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(iContA,iatomA,',iAngB,' + offsetB,iatomB,iAngC,iAngD) =&'
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'                & LP1(',iAngB,',iAngC,iAngD,iContP,IPass)'
!!$               ELSEIF(SegP)THEN
!!$                  IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(iatomA,',iAngB,',iatomB,I3,I4) = LP1(',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$                  ELSEIF(nOrbCompC.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(iatomA,',iAngB,',iatomB,I3,iContD) = LP1(',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$                  ELSEIF(nOrbCompD.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(iatomA,',iAngB,',iatomB,iContC,I4) = LP1(',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$                  ELSE
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(iatomA,',iAngB,',iatomB,iContC,iContD) = LP1(',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$                  ENDIF
!!$               ELSEIF(Seg)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(iatomA,',iAngB,',iatomB,iAngC,iAngD) = LP1(',iAngB,',iAngC,iAngD,IPass)'
!!$               ENDIF
!!$            ENDDO
!!$         ELSEIF(nOrbCompB.EQ.1)THEN
!!$            !UNROLL nOrbCompA AND UNROLL nOrbCompB
!!$            DO iAngA = 1,nOrbCompA
!!$               IF(Gen)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,' + offsetA,iatomA,iContB,iatomB,I3,I4) =&'
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'                & LP1(',iAngA,',iAngC,iAngD,iContQ,iContP,IPass)'
!!$               ELSEIF(SegQ)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,' + offsetA,iatomA,iContB,iatomB,iAngC,iAngD) =&'
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A)')'                & LP1(',iAngA,',iAngC,iAngD,iContP,IPass)'
!!$               ELSEIF(SegP)THEN
!!$                  IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iatomB,I3,I4) = LP1(',iAngA,',iAngC,iAngD,iContQ,IPass)'
!!$                  ELSEIF(nOrbCompC.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iatomB,I3,iContD) = LP1(',iAngA,',iAngC,iAngD,iContQ,IPass)'
!!$                  ELSEIF(nOrbCompD.GT.1)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iatomB,iContC,I4) = LP1(',iAngA,',iAngC,iAngD,iContQ,IPass)'
!!$                  ELSE
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iatomB,iContC,iContD) = LP1(',iAngA,',iAngC,iAngD,iContQ,IPass)'
!!$                  ENDIF
!!$               ELSEIF(Seg)THEN
!!$                  WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,iatomB,iAngC,iAngD) = LP1(',iAngA,',iAngC,iAngD,IPass)'
!!$               ENDIF
!!$            ENDDO
!!$         ELSE
!!$            !UNROLL nOrbCompA AND UNROLL nOrbCompB
!!$            DO iAngB = 1,nOrbCompB
!!$               DO iAngA = 1,nOrbCompA
!!$                  IF(Gen)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,' + offsetA,iatomA,',iAngB,' + offsetB,iatomB,I3,I4) =&'
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A)')'                & LP1(',iAngA,',',iAngB,',iAngC,iAngD,iContQ,iContP,IPass)'
!!$                  ELSEIF(SegQ)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A)')'           LP2(',iAngA,' + offsetA,iatomA,',iAngB,' + offsetB,iatomB,iAngC,iAngD) =&'
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A)')'                & LP1(',iAngA,',',iAngB,',iAngC,iAngD,iContP,IPass)'
!!$                  ELSEIF(SegP)THEN
!!$                     IF(nOrbCompC.GT.1.AND.nOrbCompD.GT.1)THEN
!!$                        WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,',iAngB,',iatomB,I3,I4) = LP1(',iAngA,',',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$                     ELSEIF(nOrbCompC.GT.1)THEN
!!$                        WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,',iAngB,',iatomB,I3,iContD) = LP1(',iAngA,',',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$                     ELSEIF(nOrbCompD.GT.1)THEN
!!$                        WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,',iAngB,',iatomB,iContC,I4) = LP1(',iAngA,',',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$                     ELSE
!!$                        WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,',iAngB,',iatomB,iContC,iContD) = LP1(',iAngA,',',iAngB,',iAngC,iAngD,iContQ,IPass)'
!!$                     ENDIF
!!$                  ELSEIF(Seg)THEN
!!$                     WRITE(ILUMOD,'(A,I1,A,I1,A,I1,A,I1,A)')'           LP2(',iAngA,',iatomA,',iAngB,',iatomB,iAngC,iAngD) = LP1(',iAngA,',',iAngB,',iAngC,iAngD,IPass)'
!!$                  ENDIF
!!$               ENDDO
!!$            ENDDO







!!$
!!$
!!$
!!$
!!$
!!$              
!!$              !(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
!!$
!!$              !(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD) 
!!$
!!$
!!$!Reorder LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
!!$!(including LHS permute) to LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
!!$!this can be done on the accelerator
!!$subroutine MainTriDistributetoLocalIntPass2CPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
!!$     & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
!!$     & MaxPasses,nPasses,TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1,&
!!$     & LocalIntPass2,IatomAPass,iatomBPass,nContQ,nContP)
!!$  implicit none
!!$  integer,intent(in) :: TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,nContQ,nContP
!!$  integer,intent(in) :: nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD
!!$  integer,intent(in) :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
!!$  logical,intent(in) :: TriangularLHSAtomLoop,Qsegmented,Psegmented
!!$  real(realk),intent(in) :: LocalIntPass1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
!!$  real(realk),intent(inout) :: LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
!!$  integer,intent(in) :: IatomAPass(MaxPasses),iatomBPass(MaxPasses)
!!$  IF(TriangularLHSAtomLoop)THEN
!!$     call TriDistributeCPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
!!$          & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
!!$          & MaxPasses,IatomAPass,iatomBPass,nPasses)
!!$  ELSE
!!$     call DistributeCPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
!!$          & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
!!$          & MaxPasses,IatomAPass,iatomBPass,nPasses)
!!$  ENDIF
!!$end subroutine MainTriDistributetoLocalIntPass2CPU
!!$
!!$!Reorder LocalIntPass(nContQ,nContP,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
!!$!(including LHS permute) to LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
!!$!this can be done on the accelerator
!!$subroutine MainTriDistributetoLocalIntPass2GPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
!!$     & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
!!$     & MaxPasses,nPasses,TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1,&
!!$     & LocalIntPass2,IatomAPass,iatomBPass,nContQ,nContP)
!!$  implicit none
!!$  integer,intent(in) :: TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,nContQ,nContP
!!$  integer,intent(in) :: nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD
!!$  integer,intent(in) :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
!!$  logical,intent(in) :: TriangularLHSAtomLoop,Qsegmented,Psegmented
!!$  real(realk),intent(in) :: LocalIntPass1(nContQ,nContP,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
!!$  real(realk),intent(inout) :: LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
!!$  integer,intent(in) :: IatomAPass(MaxPasses),iatomBPass(MaxPasses)
!!$  IF(TriangularLHSAtomLoop)THEN
!!$     call TriDistributeGPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
!!$          & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
!!$          & MaxPasses,IatomAPass,iatomBPass,nPasses)
!!$  ELSE
!!$     call DistributeGPUToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
!!$          & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
!!$          & MaxPasses,IatomAPass,iatomBPass,nPasses)
!!$  ENDIF
!!$end subroutine MainTriDistributetoLocalIntPass2GPU
!!$
!!$subroutine DistributeCPUToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
!!$     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
!!$implicit none 
!!$integer,intent(in)    :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
!!$integer,intent(in)    :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
!!$integer,intent(in)    :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
!!$real(realk),intent(in)::LP1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB,nPasses)
!!$real(realk),intent(inout)::LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!!$!local variables
!!$integer :: iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB
!!$integer :: iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB
!!$
!!$!$OMP DO COLLAPSE(3) &
!!$!$OMP PRIVATE(iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB,&
!!$!$OMP         iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB) 
!!$  DO IPass = 1,nPasses
!!$   DO iContD = 1,nContD
!!$    DO iContC = 1,nContC
!!$     DO iAngD = 1,nOrbCompD
!!$      DO iAngC = 1,nOrbCompC
!!$       IatomB = IatomBPass(IPass)
!!$       IatomA = IatomAPass(IPass)
!!$       iContQ = iContC + (iContD-1)*nContC
!!$       I4 = iAngD + (iContD-1)*nOrbCompD
!!$       I3 = iAngC + (iContC-1)*nOrbCompC
!!$       DO iContB = 1,nContB
!!$        offsetB = (iContB-1)*nOrbCompB
!!$        DO iContA = 1,nContA
!!$         iContP = iContA+(iContB-1)*nContA
!!$         offsetA = (iContA-1)*nOrbCompA
!!$         DO iAngB = 1,nOrbCompB
!!$          DO iAngA = 1,nOrbCompA
!!$           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&
!!$                & LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
!!$          ENDDO
!!$         ENDDO
!!$        ENDDO
!!$       ENDDO
!!$      ENDDO
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$  ENDDO
!!$!$OMP END DO
!!$end subroutine DistributeCPUToLocalIntPass
!!$
!!$subroutine TriDistributeCPUToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
!!$     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
!!$implicit none 
!!$integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
!!$integer,intent(in)        :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
!!$integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
!!$real(realk),intent(in)   :: LP1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB,nPasses)
!!$real(realk),intent(inout):: LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!!$!local variables
!!$integer :: iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,iContA,iContB
!!$integer :: iContC,iContD,i4,i3,offsetA,offsetB,iAngB,iAngA
!!$
!!$!$OMP DO COLLAPSE(3) &
!!$!$OMP PRIVATE(iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,iContA,iContB,iContC,&
!!$!$OMP         iContD,i4,i3,offsetA,offsetB,iAngB,iAngA) 
!!$  DO IPass = 1,nPasses
!!$   DO iContP = 1,nContA*nContA
!!$    DO iContQ = 1,nContC*nContD
!!$     DO iAngD = 1,nOrbCompD
!!$      DO iAngC = 1,nOrbCompC
!!$       IatomB = IatomBPass(IPass)
!!$       IatomA = IatomAPass(IPass)
!!$       IF(IatomA.NE.IatomB)THEN
!!$        !Ordering of Ipass is 
!!$        !iPass = 0 
!!$        !DO IatomA = 1,natomsA
!!$        ! DO IatomB = 1,IatomBend
!!$        !   iPass = iPass + 1
!!$        ! ENDDO
!!$        !ENDDO
!!$        !Where IatomBend=IatomA for triangularLHSatomLoop
!!$        !   or IatomBend=natomsB 
!!$        
!!$        iContA = iContP - ((iContP-1)/nContA)*nContA
!!$        iContB = (iContP-1)/nContA+1
!!$        iContC = iContQ - ((iContQ-1)/nContC)*nContC
!!$        iContD = (iContQ-1)/nContC+1
!!$        I4 = iAngD + (iContD-1)*nOrbCompD
!!$        I3 = iAngC + (iContC-1)*nOrbCompC
!!$        offsetB = (iContB-1)*nOrbCompA
!!$        offsetA = (iContA-1)*nOrbCompA
!!$        DO iAngB = 1,nOrbCompA
!!$         DO iAngA = 1,nOrbCompA
!!$          LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
!!$          LP2(iAngB + offsetB,iatomB,iAngA + offsetA,iatomA,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
!!$         ENDDO
!!$        ENDDO
!!$       ELSE
!!$        iContA = iContP - ((iContP-1)/nContA)*nContA
!!$        iContB = (iContP-1)/nContA+1
!!$        iContC = iContQ - ((iContQ-1)/nContC)*nContC
!!$        iContD = (iContQ-1)/nContC+1
!!$        I4 = iAngD + (iContD-1)*nOrbCompD
!!$        I3 = iAngC + (iContC-1)*nOrbCompC
!!$        offsetB = (iContB-1)*nOrbCompA
!!$        offsetA = (iContA-1)*nOrbCompA
!!$!$OMP CRITICAL
!!$        DO iAngB = 1,nOrbCompA
!!$         DO iAngA = 1,nOrbCompA
!!$          LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomA,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
!!$          LP2(iAngB + offsetB,iatomA,iAngA + offsetA,iatomA,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
!!$         ENDDO
!!$        ENDDO
!!$!$OMP END CRITICAL
!!$       ENDIF
!!$      ENDDO
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$  ENDDO
!!$!$OMP END DO
!!$end subroutine TriDistributeCPUToLocalIntPass
!!$!
!!$!      ====================
!!$!      Now the GPU routines - OpenACC 
!!$!      ====================
!!$!
!!$subroutine DistributeGPUToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
!!$     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
!!$implicit none 
!!$integer,intent(in)    :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
!!$integer,intent(in)    :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
!!$integer,intent(in)    :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
!!$real(realk),intent(in)::LP1(nContC*nContD,nContA*nContB,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
!!$real(realk),intent(inout)::LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!!$!local variables
!!$integer :: iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB
!!$integer :: iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB
!!$  DO IPass = 1,nPasses
!!$   DO iContD = 1,nContD
!!$    DO iContC = 1,nContC
!!$     DO iAngD = 1,nOrbCompD
!!$      DO iAngC = 1,nOrbCompC
!!$       IatomB = IatomBPass(IPass)
!!$       IatomA = IatomAPass(IPass)
!!$       iContQ = iContC + (iContD-1)*nContC
!!$       I4 = iAngD + (iContD-1)*nOrbCompD
!!$       I3 = iAngC + (iContC-1)*nOrbCompC
!!$       DO iContB = 1,nContB
!!$        offsetB = (iContB-1)*nOrbCompB
!!$        DO iContA = 1,nContA
!!$         iContP = iContA+(iContB-1)*nContA
!!$         offsetA = (iContA-1)*nOrbCompA
!!$         DO iAngB = 1,nOrbCompB
!!$          DO iAngA = 1,nOrbCompA
!!$           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&
!!$                & LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
!!$          ENDDO
!!$         ENDDO
!!$        ENDDO
!!$       ENDDO
!!$      ENDDO
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$  ENDDO
!!$end subroutine DistributeGPUToLocalIntPass
!!$
!!$subroutine TriDistributeGPUToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
!!$     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
!!$implicit none 
!!$integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
!!$integer,intent(in)        :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
!!$integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
!!$real(realk),intent(in)   :: LP1(nContC*nContD,nContA*nContB,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
!!$real(realk),intent(inout):: LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!!$!local variables
!!$integer :: iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,iContA,iContB
!!$integer :: iContC,iContD,i4,i3,offsetA,offsetB,iAngB,iAngA
!!$
!!$  DO IPass = 1,nPasses
!!$   DO iContP = 1,nContA*nContA
!!$    DO iContQ = 1,nContC*nContD
!!$     DO iAngD = 1,nOrbCompD
!!$      DO iAngC = 1,nOrbCompC
!!$       IatomB = IatomBPass(IPass)
!!$       IatomA = IatomAPass(IPass)
!!$       IF(IatomA.NE.IatomB)THEN
!!$        !Ordering of Ipass is 
!!$        !iPass = 0 
!!$        !DO IatomA = 1,natomsA
!!$        ! DO IatomB = 1,IatomBend
!!$        !   iPass = iPass + 1
!!$        ! ENDDO
!!$        !ENDDO
!!$        !Where IatomBend=IatomA for triangularLHSatomLoop
!!$        !   or IatomBend=natomsB 
!!$        
!!$        iContA = iContP - ((iContP-1)/nContA)*nContA
!!$        iContB = (iContP-1)/nContA+1
!!$        iContC = iContQ - ((iContQ-1)/nContC)*nContC
!!$        iContD = (iContQ-1)/nContC+1
!!$        I4 = iAngD + (iContD-1)*nOrbCompD
!!$        I3 = iAngC + (iContC-1)*nOrbCompC
!!$        offsetB = (iContB-1)*nOrbCompA
!!$        offsetA = (iContA-1)*nOrbCompA
!!$        DO iAngB = 1,nOrbCompA
!!$         DO iAngA = 1,nOrbCompA
!!$          LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) = LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
!!$          LP2(iAngB + offsetB,iatomB,iAngA + offsetA,iatomA,I3,I4) = LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
!!$         ENDDO
!!$        ENDDO
!!$       ELSE
!!$        iContA = iContP - ((iContP-1)/nContA)*nContA
!!$        iContB = (iContP-1)/nContA+1
!!$        iContC = iContQ - ((iContQ-1)/nContC)*nContC
!!$        iContD = (iContQ-1)/nContC+1
!!$        I4 = iAngD + (iContD-1)*nOrbCompD
!!$        I3 = iAngC + (iContC-1)*nOrbCompC
!!$        offsetB = (iContB-1)*nOrbCompA
!!$        offsetA = (iContA-1)*nOrbCompA
!!$        DO iAngB = 1,nOrbCompA
!!$         DO iAngA = 1,nOrbCompA
!!$          LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomA,I3,I4) = LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
!!$          LP2(iAngB + offsetB,iatomA,iAngA + offsetA,iatomA,I3,I4) = LP1(iContQ,iContP,IPass,iAngA,iAngB,iAngC,iAngD)
!!$         ENDDO
!!$        ENDDO
!!$       ENDIF
!!$      ENDDO
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$  ENDDO
!!$end subroutine TriDistributeGPUToLocalIntPass
!!$
!!$
!!$
!!$
!!$
!!$END PROGRAM
!!$
!!$!contractecoeff_gen
