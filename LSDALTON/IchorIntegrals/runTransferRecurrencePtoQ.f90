MODULE TESTMODULE
  use stringsMODULE
  logical,save :: nPrimLast
  integer,save :: nLines,nFiles
  logical,save :: CPU
CONTAINS
subroutine PASSsub
  IMPLICIT NONE
  INTEGER :: JMAX,nTUV,nTUVprev,ituvP,J,Tp,Up,Vp,N,N2,ntuvP,ituv,C
  INTEGER :: nTUVQ,nTUVTMPP,nTUVTMPQ,JPQ,JP,JQ,nTUVTMP,nTUVTMPprev
  integer :: tq,uq,vq,ituvq,ituvqminus1x,ituvqminus1y,ituvqminus1z
  integer :: iTUVplus1x,iTUVplus1y,iTUVplus1z,nTUVTMPQs,I,nTUVTMPP2,nTUVP2
  Integer :: MaxAngmomQP,nTUVplus,JTMQ,ntuvprev2,ntuvprev3,iTUVPminus1
  !    logical :: CREATED(-2:8,-2:8,-2:8)
  logical,pointer :: CREATED(:,:,:),ituvpplus1LEnTUVParray(:)
  logical,pointer :: ituvpminus1LEnTUVParray(:)
  logical :: TREC,UREC,VREC,TREC2,UREC2,VREC2,DoneCartDir(3),ituvpplus1LEnTUVP
  logical :: ituvpminus1LEnTUVP
  integer,pointer :: TUVINDEX(:,:,:)
  integer,pointer :: TINDEX(:)
  integer,pointer :: UINDEX(:)
  integer,pointer :: VINDEX(:)
  integer,pointer :: JINDEX(:)
  integer :: nTUVLIST,nTUVLISTactual,CARTDIR,iTUVQminus1,iTUVQminus2,Tqminus1
  integer,pointer :: TwoTermTUVLIST(:)
  integer :: iTUVPplus1,nLength
  integer :: LUFILE,GPUrun
  Character(len=52) :: FileName    
  Character(len=44) :: FromPrimLabel, ToPrimLabel
  Character(len=1)  :: FromLabel,ToLabel,FromExpLabel,ToExpLabel,SIGN2
  Character(len=8)  :: SegLabel
  character(len=3) :: ARCSTRING
  integer :: iseg,ifile,iseglabel
  logical :: Gen,SegQ,Segp,Seg,Seg1Prim,LOOPUNROLL,DoOpenMP,DoOpenACC
  logical :: Collapse,WRITETHING
  integer :: LUFILEMOD
  integer,pointer :: IfacX(:,:),TUVindexX(:,:)
  logical,pointer :: UniqeTUVindexX(:,:),UniqeIfacX(:,:)

  !Make Module file containing 
  !       TUVindexX1,TUVindexX2,TUVindexX3,IfacX1,IfacX2,IfacX3
  DO I =1,LEN(FileName)
     FileName(I:I) = ' '
  ENDDO
  print*,'WRITE FileName1:'
  WRITE(FileName,'(A)')'AGC_CPU_TransferRecurrenceParam.F90'
  print*,'FileName:',FileName
  LUFILEMOD = 31
  open(unit = LUFILEMOD, file=TRIM(FileName),status="unknown")
  WRITE(LUFILEMOD,'(A)')'MODULE AGC_CPU_OBS_TRParamMod'
  WRITE(LUFILEMOD,'(A)')'  '

  MaxAngmomQP = 8
  nTUVTMPP=(MaxAngmomQP+1)*(MaxAngmomQP+2)*(MaxAngmomQP+3)/6
  nTUVTMPP2=(MaxAngmomQP+1)*(MaxAngmomQP+2)*(MaxAngmomQP+3)/6
  allocate(UniqeTUVindexX(nTUVTMPP,3))
  allocate(UniqeIfacX(nTUVTMPP2,3))
  UniqeTUVindexX = .TRUE.
  UniqeIfacX = .TRUE.

  DO GPUrun = 1,2
    CPU = .TRUE.
    IF(GPUrun.EQ.2)CPU = .FALSE.
    nPrimLAST = .FALSE.
    IF(CPU)nPrimLAST = .TRUE.
    DoOpenMP = .FALSE.
    DoOpenACC = .FALSE.
    IF(CPU)DoOpenMP = .TRUE.
    IF(.NOT.CPU)DoOpenACC = .TRUE.
    COLLAPSE=.TRUE.
    IF(CPU)THEN
       ARCSTRING = 'CPU'
    ELSE
       ARCSTRING = 'GPU'
    ENDIF
    DO iseg = 1,5
       DO ifile = 1,4
          nFiles = 1
          IF(ifile.EQ.1)THEN
             FromLabel = 'A'; ToLabel = 'C'; FromExpLabel = 'B'; ToExpLabel = 'D'
             FromPrimLabel = 'iPrimB = (iPrimP-1)/nPrimA+1                '; SIGN2= '+'  
             ToPrimLabel   = 'iPrimD = (iPrimQ-1)/nPrimC+1                '
          ELSEIF(ifile.EQ.2)THEN
             FromLabel = 'B'; ToLabel = 'C'; FromExpLabel = 'A'; ToExpLabel = 'D' 
             FromPrimLabel = 'iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA'; SIGN2= '+'  
             ToPrimLabel   = 'iPrimD = (iPrimQ-1)/nPrimC+1                '
          ELSEIF(ifile.EQ.3)THEN
             FromLabel = 'A'; ToLabel = 'D'; FromExpLabel = 'B'; ToExpLabel = 'C' 
             FromPrimLabel = 'iPrimB = (iPrimP-1)/nPrimA+1                '; SIGN2= '-'  
             ToPrimLabel   = 'iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC'
          ELSEIF(ifile.EQ.4)THEN
             FromLabel = 'B'; ToLabel = 'D'; FromExpLabel = 'A'; ToExpLabel = 'C' 
             FromPrimLabel = 'iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA'; SIGN2= '-'  
             ToPrimLabel   = 'iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC'
          ENDIF
          Gen = .FALSE.; SegQ=.FALSE.; SegP=.FALSE.;Seg=.FALSE.;Seg1Prim=.FALSE.
          IF(iseg.EQ.1)THEN
             Gen = .TRUE.      ; SegLabel = 'Gen     '; iSegLabel = 3
          ELSEIF(iseg.EQ.2)THEN
             SegQ = .TRUE.     ; SegLabel = 'SegQ    '; iSegLabel = 4
          ELSEIF(iseg.EQ.3)THEN
             SegP = .TRUE.     ; SegLabel = 'SegP    '; iSegLabel = 4
          ELSEIF(iseg.EQ.4)THEN
             Seg = .TRUE.      ; SegLabel = 'Seg     '; iSegLabel = 3
          ELSEIF(iseg.EQ.5)THEN
             Seg1Prim = .TRUE. ; SegLabel = 'Seg1Prim'; iSegLabel = 8
          ENDIF
          DO I =1,LEN(FileName)
             FileName(I:I) = ' '
          ENDDO
          print*,'WRITE FileName1:'
          WRITE(FileName,'(7A,I1,A)')'AGC_',ARCSTRING,'_TransferRecurrence',FromLabel,'to',ToLabel,SegLabel(1:iSegLabel),nFiles,'.F90'
          print*,'FileName:',FileName
          LUFILE = 1
          nLines = 0 
          open(unit = LUFILE, file=TRIM(FileName),status="unknown")
          WRITE(LUFILE,'(7A,I1)')'MODULE AGC_',ARCSTRING,'_OBS_TRMOD',FromLabel,'to',ToLabel,SegLabel(1:iSegLabel),nFiles
          WRITE(LUFILE,'(A)')' use IchorPrecisionMod'
          WRITE(LUFILE,'(A)')'  '
          WRITE(LUFILE,'(A)')' CONTAINS'
          nLines = nLines + 4

          MaxAngmomQP = 8
          IF((ifile.EQ.2))MaxAngmomQP = 6 !BtoC (PDDP)
          IF((ifile.EQ.4))MaxAngmomQP = 6 !BtoD (PDPD)
          IF((ifile.EQ.3))MaxAngmomQP = 7 !AtoD (DDPD)
          DO JMAX=2,MaxAngmomQP
             DO JP = 1, JMAX
                JQ = JMAX-JP
                IF(JQ.GT.JP)CYCLE
                IF(JQ.EQ.0)CYCLE
                IF(JQ.GT.4)CYCLE
                IF((ifile.EQ.2).AND.JP.GT.3)CYCLE
                IF((ifile.EQ.4).AND.JP.GT.3)CYCLE
                IF((ifile.EQ.2.OR.ifile.EQ.4).AND.JQ.GT.3)CYCLE
                IF((ifile.EQ.3).AND.JQ.GT.3)CYCLE
                IF(JP.GT.4)CYCLE
                IF(JMAX.LT.5)THEN
                   LOOPUNROLL = .TRUE.
                ELSE
                   LOOPUNROLL = .FALSE.
                ENDIF
                nTUVP = (JP+1)*(JP+2)*(JP+3)/6   
                nTUVQ = (JQ+1)*(JQ+2)*(JQ+3)/6   
                JPQ = JP + JQ
                nTUV = (JMAX+1)*(JMAX+2)*(JMAX+3)/6   
                nTUVPLUS = (JMAX+2)*(JMAX+3)*(JMAX+4)/6   
                allocate(TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
                TUVINDEX = 0
                allocate(TINDEX(nTUVPLUS))
                allocate(UINDEX(nTUVPLUS))
                allocate(VINDEX(nTUVPLUS))
                allocate(JINDEX(nTUVPLUS))

                nTUVprev3 = (JMAX-2)*(JMAX-1)*(JMAX)/6
                nTUVprev2 = (JMAX-1)*(JMAX)*(JMAX+1)/6
                nTUVprev = (JMAX)*(JMAX+1)*(JMAX+2)/6
                ituv = 0 
                TUVINDEX = 0
                DO J = 0, JMAX+1
                   DO Tp=J,0,-1       
                      DO Up=J-Tp,0,-1
                         Vp=J-Tp-Up
                         ituv = ituv + 1 
                         TUVINDEX(Tp,Up,Vp) = ituv
                         TINDEX(iTUV) = Tp
                         UINDEX(iTUV) = Up
                         VINDEX(iTUV) = Vp
                         JINDEX(iTUV) = J
                      ENDDO
                   ENDDO
                ENDDO
!                print*,'nLines',nLines,'nLines.GT.1000',nLines.GT.1000
                IF(nLines.GT.1000)THEN
                   WRITE(LUFILE,'(7A,I1)')'END MODULE AGC_',ARCSTRING,'_OBS_TRMOD',FromLabel,'to',ToLabel,SegLabel(1:iSegLabel),nFiles
                   DO I =1,LEN(FileName)
                      FileName(I:I) = ' '
                   ENDDO
                   close(unit = LUFILE)

                   nFiles = nFiles + 1 
!                   print*,'WRITE FileName2:'
                   WRITE(FileName,'(7A,I1,A)')'AGC_',ARCSTRING,'_TransferRecurrence',FromLabel,'to',ToLabel,SegLabel(1:iSegLabel),nFiles,'.F90'
!                   print*,'FileName:',FileName
                   LUFILE = 1
                   nLines = 0 
                   open(unit = LUFILE, file=TRIM(FileName),status="unknown")
                   WRITE(LUFILE,'(7A,I1)')'MODULE AGC_',ARCSTRING,'_OBS_TRMOD',FromLabel,'to',ToLabel,SegLabel(1:iSegLabel),nFiles
                   WRITE(LUFILE,'(A)')' use IchorPrecisionMod'
                   WRITE(LUFILE,'(A)')'  '
                   WRITE(LUFILE,'(A)')' CONTAINS'
                   nLines=0
                ENDIF

                call initString(1)          
                call AddToString('subroutine TransferRecurrence'//ARCSTRING//'P')
                call AddToString(JP)
                call AddToString('Q')
                call AddToString(JQ)
                call AddToString(FromLabel)
                call AddToString('to')
                call AddToString(ToLabel)
                call AddToString(SegLabel(1:iSegLabel))
                call AddToString('(nPasses,nPrimP,nPrimQ,reducedExponents,&')
                call writeString(LUFILE)
                nLines = nLines + 1

                call initString(9)          
                call AddToString('& Pexp,Qexp,Pdistance12,Qdistance12,')
                call AddToString(FromExpLabel)
                call AddToString('exp,')
                call AddToString(ToExpLabel)
                call AddToString('exp,nPrimA,nPrimB,nPrimC,nPrimD,&')
                call writeString(LUFILE)
                nLines = nLines + 1
                call initString(9)          
                call AddToString('& MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2')
                IF(DoOpenACC)THEN
                   call AddToString(',iASync)')
                ELSE
                   call AddToString(')')
                ENDIF
                call writeString(LUFILE)
                nLines = nLines + 1
                IF(.NOT.LOOPUNROLL)THEN
                   WRITE(LUFILE,'(A,A,A)')'  use AGC_',ARCSTRING,'_OBS_TRParamMod'
                ENDIF
                WRITE(LUFILE,'(A)')'  implicit none'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses'
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)'
                ELSE
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)'
                ENDIF
                WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
                nLines = nLines + 5
                call initString(2)          
                call AddToString('real(realk),intent(in) :: ')
                call AddToString(FromExpLabel)
                IF(.NOT.Seg1Prim)THEN
                   call AddToString('exp(nPrim')
                   call AddToString(FromExpLabel)
                ELSE
                   call AddToString('exp(1')
                ENDIF
                call AddToString('),')
                call AddToString(ToExpLabel)
                IF(.NOT.Seg1Prim)THEN
                   call AddToString('exp(nPrim')
                   call AddToString(ToExpLabel)
                ELSE
                   call AddToString('exp(1')
                ENDIF
                call AddToString(')')
                call writeString(LUFILE)
                nLines = nLines + 1
                !          WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Bexp(nPrimB),Dexp(nPrimD)'
                IF(nPrimLast)THEN
                   IF(.NOT.Seg1Prim)THEN
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPrimQ*nPrimP*nPasses)'
                      ELSE
                         WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPrimQ,nPrimP,nPasses)'
                      ENDIF
                   ELSE
                      WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPasses)'
                   ENDIF
                ELSE
                   IF(.NOT.Seg1Prim)THEN
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(nPrimQ*nPrimP*nPasses,',nTUV,')'
                      ELSE
                         WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,',nTUV,')'
                      ENDIF
                   ELSE
                      WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(nPasses,',nTUV,')'
                   ENDIF
                ENDIF
                nLines = nLines + 1
                IF(COLLAPSE)THEN
                   IF(nPrimLast)THEN
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPrimQ*nPrimP*nPasses)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPrimP*nPasses)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPrimQ*nPasses)'
                      ELSEIF(Seg.OR.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPasses)'
                      ENDIF
                   ELSE
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(nPrimQ*nPrimP*nPasses,',nTUVP,',',nTUVQ,')'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(nPrimP*nPasses,',nTUVP,',',nTUVQ,')'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,',nTUVP,',',nTUVQ,')'
                      ELSEIF(Seg.OR.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(nPasses,',nTUVP,',',nTUVQ,')'
                      ENDIF
                   ENDIF
                ELSE
                   IF(nPrimLast)THEN
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPrimQ*nPrimP*nPasses)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPrimP,nPasses)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPrimQ,nPasses)'
                      ELSEIF(Seg.OR.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPasses)'
                      ENDIF
                   ELSE
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(nPrimQ*nPrimP*nPasses,',nTUVP,',',nTUVQ,')'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(nPrimP,nPasses,',nTUVP,',',nTUVQ,')'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(nPrimQ,nPasses,',nTUVP,',',nTUVQ,')'
                      ELSEIF(Seg.OR.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(nPasses,',nTUVP,',',nTUVQ,')'
                      ENDIF
                   ENDIF
                ENDIF
                nLines = nLines + 1
                IF(nPrimLast)THEN
                   WRITE(LUFILE,'(A)')'!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)'
                ELSE
                   WRITE(LUFILE,'(A)')'!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)'
                ENDIF
                nLines = nLines + 1
                IF(DoOpenACC)WRITE(LUFILE,'(A)')'  integer(kind=acckind),intent(in) :: iASync'
                WRITE(LUFILE,'(A)')'  !Local variables'
                !             WRITE(LUFILE,'(A)')'  real(realk) :: Pexpfac,PREF'
                !             WRITE(LUFILE,'(A,i4,A)')'  real(realk) :: TwoTerms(',MAX(1,nTUVprev2-nTUVprev3),')'
                WRITE(LUFILE,'(A,I3,A,I3,A)')'  real(realk) :: Tmp0(',nTUVP,',',nTUVQ,')'
                DO JTMQ=1,JQ-1
                   nTUVTMPQs=(JTMQ)*(JTMQ+1)*(JTMQ+2)/6
                   nTUVTMPQ=(JTMQ+1)*(JTMQ+2)*(JTMQ+3)/6
                   nTUVTMPP=(JPQ-JTMQ+1)*(JPQ-JTMQ+2)*(JPQ-JTMQ+3)/6
                   if(JTMQ.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMQ,'(',nTUVP+1,':',nTUVTMPP,',',nTUVTMPQs+1,':',nTUVTMPQ,')'
                   else
                      WRITE(LUFILE,'(A,I2,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMQ,'(',nTUVP+1,':',nTUVTMPP,',',nTUVTMPQs+1,':',nTUVTMPQ,')'
                   endif
                   nLines = nLines + 1
                ENDDO
                WRITE(LUFILE,'(A)')'!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering'
                WRITE(LUFILE,'(A)')'  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB'
                nLines = nLines + 2
                IF(.NOT.Seg1Prim)THEN
                   call initString(2)          
                   call AddToString('integer :: iPrim')
                   call AddToString(FromExpLabel)
                   call AddToString(',iPrim')
                   call AddToString(ToExpLabel)
                   call writeString(LUFILE)
                   nLines = nLines + 1
                   !          WRITE(LUFILE,'(A)')'  integer :: iPrimB,iPrimD'
                ENDIF
                WRITE(LUFILE,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk'
                WRITE(LUFILE,'(A)')'  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP'
                call initString(2)          
                call AddToString('real(realk) :: exp')
                call AddToString(FromExpLabel)
                call AddToString('X,exp')
                call AddToString(FromExpLabel)
                call AddToString('Y,exp')
                call AddToString(FromExpLabel)
                call AddToString('Z')
                call writeString(LUFILE); nLines = nLines + 1
!                WRITE(LUFILE,'(A)')'  real(realk) :: expBX,expBY,expBZ'
                WRITE(LUFILE,'(A)')'  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq'
                nLines = nLines + 1
                !==========================================================================================================
                !         Build the TUVindexX
                !==========================================================================================================
                allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
                CREATED  = .FALSE.
                CREATED(0,0,0) = .TRUE.
                nTUVTMPP=(JPQ-1+1)*(JPQ-1+2)*(JPQ-1+3)/6                         
                allocate(TUVindexX(nTUVTMPP,3))
                DO JTMQ=1,1!JQ
                   DoneCartDir = .FALSE.
                   DO Tq=JTMQ,0,-1       
                      DO Uq=JTMQ-Tq,0,-1
                         Vq=JTMQ-Tq-Uq  
                         CALL DETERMINE_CARTDIR(CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1,Tq,Uq,Vq,CREATED,JMAX,TUVINDEX)
                         nTUVTMPP=(JPQ-JTMQ+1)*(JPQ-JTMQ+2)*(JPQ-JTMQ+3)/6
                         IF(.NOT.DoneCartDir(CARTDIR))THEN
                            WRITETHING=.FALSE.
                         !   print*,'A LOOPUNROLL',LOOPUNROLL
                            IF(.NOT.LOOPUNROLL)THEN
                               WRITETHING = .FALSE.
                          !     print*,'A UniqeTUVindexX(nTUVTMPP,CARTDIR)',UniqeTUVindexX(nTUVTMPP,CARTDIR),nTUVTMPP,CARTDIR
                               IF(UniqeTUVindexX(nTUVTMPP,CARTDIR))THEN
                           !       print*,'A GPU  Uni TUVindex'
                                  call initString(2)          
                                  call AddToString('integer,parameter, dimension(')
                                  call AddToString(nTUVTMPP)
                                  call AddToString(') :: TUVindex')
                                  !                      call AddToString(JTMQ)
                                  call AddToString('X')
                                  call AddToString(CARTDIR)
                                  call AddToString('_')
                                  call AddToString(nTUVTMPP)
                                  UniqeTUVindexX(nTUVTMPP,CARTDIR) = .FALSE.
                                  call AddToString(' = (/ ')
                                  WRITETHING = .TRUE.
                            !   ELSE
                            !      print*,'A GPU  NOT Uni TUVindex'
                               ENDIF
                            !   print*,'A WRITETHING',WRITETHING
                            ENDIF
                            nLength = 10
                            do iTUVP = 1,nTUVTMPP!nTUVP
                               nLength = nLength + 1
                               Tp = Tindex(iTUVp) 
                               Up = Uindex(iTUVp) 
                               Vp = Vindex(iTUVp)
                               IF(CARTDIR.EQ.1)THEN
                                  iTUVPplus1 = TUVINDEX(Tp+1,Up,Vp)
                               ELSEIF(CARTDIR.EQ.2)THEN
                                  iTUVPplus1 = TUVINDEX(Tp,Up+1,Vp)
                               ELSEIF(CARTDIR.EQ.3)THEN
                                  iTUVPplus1 = TUVINDEX(Tp,Up,Vp+1)
                               ENDIF
                               IF(.NOT.LOOPUNROLL)THEN
                                  !print*,'A2 Add to string',ituvpplus1
                                  IF(WRITETHING)call AddToString(ituvpplus1)
                               ENDIF
                               TUVindexX(iTUVP,CARTDIR) = ituvpplus1
                               IF(.NOT.LOOPUNROLL)THEN
                                  IF(iTUVP.NE.nTUVTMPP)THEN
                                     IF(WRITETHING)call AddToString(',')
                                  ENDIF
                                  IF(nLength.EQ.17)THEN
                                     nLength = 0
                                     IF(WRITETHING)THEN
                                        call AddToString('&')
                                        !print*,'A1 WRITE TO LUMOD'
                                        call writeString(LUFILEMOD); !nLines = nLines + 1
                                        call initString(5)          
                                        call AddToString('     & ')
                                     ENDIF
                                  ENDIF
                               ENDIF
                            enddo
                            IF(.NOT.LOOPUNROLL)THEN
                               IF(WRITETHING)THEN
                                  call AddToString(' /)')
                                  !print*,'A2 WRITE TO LUMOD'
                                  call writeString(LUFILEMOD);! nLines = nLines + 1                   
                               ENDIF
                            ENDIF
                            DoneCartDir(CARTDIR) = .TRUE.
                         ENDIF
                         CREATED(Tq,Uq,Vq) = .TRUE.
                      ENDDO
                   ENDDO
                ENDDO
                !==========================================================================================================
                !         Build the IfacX integer          
                !==========================================================================================================
                CREATED  = .FALSE.
                CREATED(0,0,0) = .TRUE.
                nTUVTMPP2=(JPQ-1)*(JPQ-1+1)*(JPQ-1+2)/6                         
                allocate(IfacX(nTUVTMPP2,3))
                DO JTMQ=1,1!JQ
                   DoneCartDir = .FALSE.
                   DO Tq=JTMQ,0,-1       
                      DO Uq=JTMQ-Tq,0,-1
                         Vq=JTMQ-Tq-Uq  
                         CALL DETERMINE_CARTDIR(CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1,Tq,Uq,Vq,CREATED,JMAX,TUVINDEX)
                         nTUVTMPP2=(JPQ-JTMQ)*(JPQ-JTMQ+1)*(JPQ-JTMQ+2)/6
                         IF(.NOT.DoneCartDir(CARTDIR))THEN
                            WRITETHING=.FALSE.
                          !  print*,'B LOOPUNROLL',LOOPUNROLL
                            IF(.NOT.LOOPUNROLL)THEN
                               WRITETHING=.FALSE.
                           !    print*,'A UniqeIfacX(nTUVTMPP,CARTDIR)',UniqeIfacX(nTUVTMPP2,CARTDIR),nTUVTMPP2,CARTDIR
                               IF(UniqeIfacX(nTUVTMPP2,CARTDIR))THEN
                                  call initString(2)          
                                  call AddToString('integer,parameter, dimension(')
                                  call AddToString(nTUVTMPP2)
                                  call AddToString(') :: Ifac')
                                  call AddToString('X')
                                  call AddToString(CARTDIR)
                                  call AddToString('_')
                                  call AddToString(nTUVTMPP2)
                                  UniqeIfacX(nTUVTMPP2,CARTDIR) = .FALSE.
                                  call AddToString(' = (/ ')
                                  WRITETHING = .TRUE.
                               ELSE
                                  WRITETHING = .FALSE.
                               ENDIF
                            ENDIF                            
                            !print*,'B WRITETHING',WRITETHING
                            nLength = 10
                            nTUVTMPP2=(JPQ-JTMQ)*(JPQ-JTMQ+1)*(JPQ-JTMQ+2)/6
                            do iTUVPminus1 = 1,nTUVTMPP2
                               nLength = nLength + 1
                               Tp = Tindex(iTUVpminus1) 
                               Up = Uindex(iTUVpminus1) 
                               Vp = Vindex(iTUVpminus1)
                               IF(CARTDIR.EQ.1)THEN
                                  I = Tp+1
                               ELSEIF(CARTDIR.EQ.2)THEN
                                  I = Up+1
                               ELSEIF(CARTDIR.EQ.3)THEN
                                  I = Vp+1
                               ENDIF
                               IF(.NOT.LOOPUNROLL)THEN
                                  !print*,'I',I
                                  IF(WRITETHING)call AddToString(I)
                               ENDIF
                               IfacX(iTUVPminus1,CARTDIR) = I
                               IF(.NOT.LOOPUNROLL)THEN
                                  IF(iTUVPminus1.NE.nTUVTMPP2)THEN
                                     IF(WRITETHING)call AddToString(',')
                                  ENDIF
                                  IF(nLength.EQ.17)THEN
                                     nLength = 0
                                     IF(WRITETHING)THEN
                                        call AddToString('&')
                                        !print*,'B1 WRITE LUMOD1'
                                        call writeString(LUFILEMOD)!; nLines = nLines + 1
                                        call initString(5)          
                                        call AddToString('     & ')
                                     ENDIF
                                  ENDIF
                               ENDIF
                            enddo
                            IF(.NOT.LOOPUNROLL)THEN
                               IF(WRITETHING)THEN
                                  call AddToString(' /)')
                                  !print*,'B2 WRITE LUMOD2'
                                  call writeString(LUFILEMOD)!; nLines = nLines + 1                   
                               ENDIF
                            ENDIF
                            DoneCartDir(CARTDIR) = .TRUE.
                            CREATED(Tq,Uq,Vq) = .TRUE.
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
                deallocate(CREATED)
               ! print*,'DONE'

                !OPENMP
                IF(DoOpenMP)THEN
!                   WRITE(LUFILE,'(A)')'!$OMP PARALLEL DO DEFAULT(none)&'
                   WRITE(LUFILE,'(A)')'!$OMP DO &'
                   WRITE(LUFILE,'(A)')'!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&'
                   IF(Seg)THEN
                      WRITE(LUFILE,'(A)')'!$OMP         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&'
                   ELSE
                      WRITE(LUFILE,'(A)')'!$OMP         iP,iPrimQ,iPrimP,iPassP,&'
                   ENDIF
                   WRITE(LUFILE,'(A,A,A,A,A,A,A)')'!$OMP         exp',FromExpLabel,'X,exp',FromExpLabel,'Y,exp',FromExpLabel,'Z,&'
                   IF(.NOT.seg1Prim)THEN
                      WRITE(LUFILE,'(A,A,A,A,A,A,A)')'!$OMP         iPrim',FromExpLabel,',iPrim',ToExpLabel,',&'
                   ENDIF
                   WRITE(LUFILE,'(A)')'!$OMP         Tmp0,&'
                   DO JTMQ=1,JQ-1
                      if(JTMQ.LT.10)THEN
                         WRITE(LUFILE,'(A,I1,A)')'!$OMP         Tmp',JTMQ,',&'
                      else
                         WRITE(LUFILE,'(A,I2,A)')'!$OMP         Tmp',JTMQ,',&'
                      endif
                   ENDDO                   
                   WRITE(LUFILE,'(A)')'!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)'
!                   WRITE(LUFILE,'(A)')'!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &'
!                   WRITE(LUFILE,'(A)')'!$OMP SHARED(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&'
!                   WRITE(LUFILE,'(A)')'!$OMP        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&'
!           WRITE(LUFILE,'(A,A,A,A,A)')'!$OMP       ',FromExpLabel,'exp,',ToExpLabel,'exp,&'
!                   WRITE(LUFILE,'(A)')'!$OMP        IatomApass,IatomBpass,Aux2,Aux)'
                   nLines = nLines + 6
                ENDIF
                IF(DoOpenACC)THEN
                   WRITE(LUFILE,'(A)')'!$ACC PARALLEL LOOP &'
                   WRITE(LUFILE,'(A)')'!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&'
                   IF(Seg)THEN
                      WRITE(LUFILE,'(A)')'!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&'
                   ELSE
                      WRITE(LUFILE,'(A)')'!$ACC         iP,iPrimQ,iPrimP,iPassP,&'
                   ENDIF
                   WRITE(LUFILE,'(A,A,A,A,A,A,A)')'!$ACC         exp',FromExpLabel,'X,exp',FromExpLabel,'Y,exp',FromExpLabel,'Z,&'
                   IF(.NOT.seg1Prim)THEN
                      WRITE(LUFILE,'(A,A,A,A,A,A,A)')'!$ACC         iPrim',FromExpLabel,',iPrim',ToExpLabel,',&'
                   ENDIF
                   WRITE(LUFILE,'(A)')'!$ACC         Tmp0,&'
                   DO JTMQ=1,JQ-1
                      if(JTMQ.LT.10)THEN
                         WRITE(LUFILE,'(A,I1,A)')'!$ACC         Tmp',JTMQ,',&'
                      else
                         WRITE(LUFILE,'(A,I2,A)')'!$ACC         Tmp',JTMQ,',&'
                      endif
                   ENDDO                   
                   WRITE(LUFILE,'(A)')'!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &'
                   WRITE(LUFILE,'(A)')'!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&'
                   IF(.NOT.LOOPUNROLL)THEN
                      call initString(0)          
                      call AddToString('!$ACC        TUVindexX1_')
                      call AddToString(nTUVTMPP)
                      call AddToString(',TUVindexX2_')
                      call AddToString(nTUVTMPP)
                      call AddToString(',TUVindexX3_')
                      call AddToString(nTUVTMPP)
                      call AddToString(', &')
                      call writeString(LUFILE)
                      call initString(0)          
                      call AddToString('!$ACC        IfacX1_')
                      call AddToString(nTUVTMPP2)
                      call AddToString(',IfacX2_')
                      call AddToString(nTUVTMPP2)
                      call AddToString(',IfacX3_')
                      call AddToString(nTUVTMPP2)
                      call AddToString(', &')
                      call writeString(LUFILE)
                   ENDIF
                   WRITE(LUFILE,'(A)')'!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&'
                   WRITE(LUFILE,'(A,A,A,A,A)')'!$ACC       ',FromExpLabel,'exp,',ToExpLabel,'exp,&'
                   WRITE(LUFILE,'(A)')'!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)'
                   nLines = nLines + 6
                ENDIF

                !START LOOP
                IF(COLLAPSE)THEN
                   IF(Gen)THEN
                      WRITE(LUFILE,'(A)') '  DO iP = 1,nPrimQ*nPrimP*nPasses'
                      WRITE(LUFILE,'(A)') '   iPrimQ = mod(IP-1,nPrimQ)+1'
                      WRITE(LUFILE,'(A)') '   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1'
                      WRITE(LUFILE,'(A)') '   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1'
                   ELSEIF(SegP)THEN
                      WRITE(LUFILE,'(A)') '  DO iP = 1,nPrimQ*nPasses'
                      !init AUX
                      IF(COLLAPSE)THEN
                         IF(nPrimLast)THEN
                            WRITE(LUFILE,'(A,I3)')'   DO iTUVQ=1,',nTUVQ
                            WRITE(LUFILE,'(A,I3)')'    DO iTUVP=1,',nTUVP
                            WRITE(LUFILE,'(A)')   '     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk'
                            WRITE(LUFILE,'(A)')   '    ENDDO'
                            WRITE(LUFILE,'(A)')   '   ENDDO'
                         ELSE
                            WRITE(LUFILE,'(A,I3)')'   DO iTUVQ=1,',nTUVQ
                            WRITE(LUFILE,'(A,I3)')'    DO iTUVP=1,',nTUVP
                            WRITE(LUFILE,'(A)')   '     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk'
                            WRITE(LUFILE,'(A)')   '    ENDDO'
                            WRITE(LUFILE,'(A)')   '   ENDDO'
                         ENDIF
                         nLines = nLines + 5
                      ELSE
                         !no need to init
                      ENDIF
                      WRITE(LUFILE,'(A)') '   DO iPrimP=1, nPrimP'
                      WRITE(LUFILE,'(A)') '    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ'
                      WRITE(LUFILE,'(A)') '    iPassP = (iP-1)/nPrimQ + 1'
                   ELSEIF(SegQ)THEN
                      WRITE(LUFILE,'(A)') '  DO iP = 1,nPrimP*nPasses'
                      IF(COLLAPSE)THEN
                         IF(nPrimLast)THEN
                            WRITE(LUFILE,'(A,I3)')'   DO iTUVQ=1,',nTUVQ
                            WRITE(LUFILE,'(A,I3)')'    DO iTUVP=1,',nTUVP
                            WRITE(LUFILE,'(A)')   '     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk'
                            WRITE(LUFILE,'(A)')   '    ENDDO'
                            WRITE(LUFILE,'(A)')   '   ENDDO'
                         ELSE
                            WRITE(LUFILE,'(A,I3)')'   DO iTUVQ=1,',nTUVQ
                            WRITE(LUFILE,'(A,I3)')'    DO iTUVP=1,',nTUVP
                            WRITE(LUFILE,'(A)')   '     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk'
                            WRITE(LUFILE,'(A)')   '    ENDDO'
                            WRITE(LUFILE,'(A)')   '   ENDDO'
                         ENDIF
                         nLines = nLines + 5
                      ELSE
                         !no need to init
                      ENDIF
                      WRITE(LUFILE,'(A)') '   DO iPrimQ=1, nPrimQ'
                      WRITE(LUFILE,'(A)') '    iPrimP = iP - ((iP-1)/nPrimP)*nPrimP'
                      WRITE(LUFILE,'(A)') '    iPassP = (iP-1)/nPrimP + 1'
                   ELSEIF(Seg)THEN
                      WRITE(LUFILE,'(A)') '  DO iP = 1,nPasses'
                      IF(COLLAPSE)THEN
                         IF(nPrimLast)THEN
                            WRITE(LUFILE,'(A,I3)')'   DO iTUVQ=1,',nTUVQ
                            WRITE(LUFILE,'(A,I3)')'    DO iTUVP=1,',nTUVP
                            WRITE(LUFILE,'(A)')   '     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk'
                            WRITE(LUFILE,'(A)')   '    ENDDO'
                            WRITE(LUFILE,'(A)')   '   ENDDO'
                         ELSE
                            WRITE(LUFILE,'(A,I3)')'   DO iTUVQ=1,',nTUVQ
                            WRITE(LUFILE,'(A,I3)')'    DO iTUVP=1,',nTUVP
                            WRITE(LUFILE,'(A)')   '     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk'
                            WRITE(LUFILE,'(A)')   '    ENDDO'
                            WRITE(LUFILE,'(A)')   '   ENDDO'
                         ENDIF
                         nLines = nLines + 5
                      ELSE
                         !no need to init
                      ENDIF
                      WRITE(LUFILE,'(A)') '   DO iPrimQP=1,nPrimQ*nPrimP'
                      WRITE(LUFILE,'(A)') '    iPrimQ = iPrimQP - ((iPrimQP-1)/nPrimQ)*nPrimQ'
                      WRITE(LUFILE,'(A)') '    iPrimP = (iPrimQP-1)/nPrimQ + 1'       
                      WRITE(LUFILE,'(A)') '    iPassP = iP'
                   ELSEIF(Seg1Prim)THEN
                      WRITE(LUFILE,'(A)') '  DO iP = 1,nPasses'
                      WRITE(LUFILE,'(A)') '   iPassP = iP'
                      WRITE(LUFILE,'(A)') '   iPrimP=1'
                      WRITE(LUFILE,'(A)') '   iPrimQ=1'
                   ENDIF
                   nLines = nLines + 5
                   WRITE(LUFILE,'(A)')'   Xcd = Qdistance12(1)'
                   WRITE(LUFILE,'(A)')'   Ycd = Qdistance12(2)'
                   WRITE(LUFILE,'(A)')'   Zcd = Qdistance12(3)'
                   WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                   WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                   WRITE(LUFILE,'(A)')'   Xab = Pdistance12(1,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Yab = Pdistance12(2,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Zab = Pdistance12(3,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    expP = Pexp(iPrimP)'
                   nLines = nLines + 9
                   IF(FromExpLabel.EQ.'A')THEN
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,A)')'    ',FromPrimLabel
                         !          WRITE(LUFILE,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'          
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'X = -',FromExpLabel,'exp(iPrim',FromExpLabel,')*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Y = -',FromExpLabel,'exp(iPrim',FromExpLabel,')*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Z = -',FromExpLabel,'exp(iPrim',FromExpLabel,')*Zab'
                      ELSE
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'X = -',FromExpLabel,'exp(1)*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Y = -',FromExpLabel,'exp(1)*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Z = -',FromExpLabel,'exp(1)*Zab'
                      ENDIF
                   ELSE
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,A)')'    ',FromPrimLabel
                         !          WRITE(LUFILE,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'          
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'X = ',FromExpLabel,'exp(iPrim',FromExpLabel,')*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Y = ',FromExpLabel,'exp(iPrim',FromExpLabel,')*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Z = ',FromExpLabel,'exp(iPrim',FromExpLabel,')*Zab'
                      ELSE
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'X = ',FromExpLabel,'exp(1)*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Y = ',FromExpLabel,'exp(1)*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Z = ',FromExpLabel,'exp(1)*Zab'
                      ENDIF
                   ENDIF
                   nLines = nLines + 3
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A,A)')'     ',ToPrimLabel
                   ENDIF
                ELSE ! no collapse
                   WRITE(LUFILE,'(A)')'  DO iPassP = 1,nPasses'                
                   WRITE(LUFILE,'(A)')'   Xcd = Qdistance12(1)'
                   WRITE(LUFILE,'(A)')'   Ycd = Qdistance12(2)'
                   WRITE(LUFILE,'(A)')'   Zcd = Qdistance12(3)'
                   WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                   WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                   WRITE(LUFILE,'(A)')'   Xab = Pdistance12(1,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Yab = Pdistance12(2,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Zab = Pdistance12(3,iAtomA,iAtomB)'
                   nLines = nLines + 9
                   IF(Seg)THEN
                      WRITE(LUFILE,'(A,I3)')'   DO iTUVQ=1,',nTUVQ
                      WRITE(LUFILE,'(A,I3)')'    DO iTUVP=1,',nTUVP
                      IF(nPrimLast)THEN
                         WRITE(LUFILE,'(A)')'     Aux2(iTUVP,iTUVQ,iPassP) = 0.0E0_realk'
                      ELSE
                         WRITE(LUFILE,'(A)')'     Aux2(iPassP,iTUVP,iTUVQ) = 0.0E0_realk'
                      ENDIF
                      WRITE(LUFILE,'(A)')'    ENDDO'
                      WRITE(LUFILE,'(A)')'   ENDDO'
                   ENDIF
                   nLines = nLines + 5
                   IF(SegP)THEN
                      WRITE(LUFILE,'(A)')'   DO iPrimQ = 1,nPrimQ'
                      WRITE(LUFILE,'(A,I3)')'    DO iTUVQ=1,',nTUVQ
                      WRITE(LUFILE,'(A,I3)')'     DO iTUVP=1,',nTUVP
                      IF(nPrimLast)THEN
                         WRITE(LUFILE,'(A)')'      Aux2(iTUVP,iTUVQ,iPrimQ,iPassP) = 0.0E0_realk'
                      ELSE
                         WRITE(LUFILE,'(A)')'      Aux2(iPassP,iTUVP,iTUVQ,iPrimQ) = 0.0E0_realk'
                      ENDIF
                      WRITE(LUFILE,'(A)')'     ENDDO'
                      WRITE(LUFILE,'(A)')'    ENDDO'
                      WRITE(LUFILE,'(A)')'   ENDDO'
                      nLines = nLines + 8
                   ENDIF
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'   IP = (iPassP-1)*nPrimQ*nPrimP'
                      WRITE(LUFILE,'(A)')'   DO iPrimP=1, nPrimP'
                      nLines = nLines + 2
                   ELSE
                      WRITE(LUFILE,'(A)')'   IP = iPassP'
                      WRITE(LUFILE,'(A)')'   iPrimP=1'
                      nLines = nLines + 2
                   ENDIF
                   IF(SegQ)THEN
                      WRITE(LUFILE,'(A,I3)')'    DO iTUVQ=1,',nTUVQ
                      WRITE(LUFILE,'(A,I3)')'     DO iTUVP=1,',nTUVP
                      IF(nPrimLast)THEN
                         WRITE(LUFILE,'(A)')'      Aux2(iTUVP,iTUVQ,iPrimP,iPassP) = 0.0E0_realk'
                      ELSE
                         WRITE(LUFILE,'(A)')'      Aux2(iPrimP,iPassP,iTUVP,iTUVQ) = 0.0E0_realk'
                      ENDIF
                      WRITE(LUFILE,'(A)')'     ENDDO'
                      WRITE(LUFILE,'(A)')'    ENDDO'
                      nLines = nLines + 6
                   ENDIF
                   WRITE(LUFILE,'(A)')'    expP = Pexp(iPrimP)'
                   nLines = nLines + 1
                   IF(FromExpLabel.EQ.'A')THEN
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,A)')'    ',FromPrimLabel
                         !          WRITE(LUFILE,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'          
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'X = -',FromExpLabel,'exp(iPrim',FromExpLabel,')*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Y = -',FromExpLabel,'exp(iPrim',FromExpLabel,')*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Z = -',FromExpLabel,'exp(iPrim',FromExpLabel,')*Zab'
                      ELSE
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'X = -',FromExpLabel,'exp(1)*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Y = -',FromExpLabel,'exp(1)*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Z = -',FromExpLabel,'exp(1)*Zab'
                      ENDIF
                      nLines = nLines + 3
                   ELSE
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,A)')'    ',FromPrimLabel
                         !          WRITE(LUFILE,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'          
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'X = ',FromExpLabel,'exp(iPrim',FromExpLabel,')*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Y = ',FromExpLabel,'exp(iPrim',FromExpLabel,')*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Z = ',FromExpLabel,'exp(iPrim',FromExpLabel,')*Zab'
                      ELSE
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'X = ',FromExpLabel,'exp(1)*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Y = ',FromExpLabel,'exp(1)*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',FromExpLabel,'Z = ',FromExpLabel,'exp(1)*Zab'
                      ENDIF
                      nLines = nLines + 3
                   ENDIF
                   !          WRITE(LUFILE,'(A)')'    expBX = Bexp(iPrimB)*Zab'
                   !          WRITE(LUFILE,'(A)')'    expBY = Bexp(iPrimB)*Zab'
                   !          WRITE(LUFILE,'(A)')'    expBZ = Bexp(iPrimB)*Zab'                   
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    DO iPrimQ=1, nPrimQ'
                   ELSE
                      WRITE(LUFILE,'(A)')'    iPrimQ=1'
                   ENDIF
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A,A)')'     ',ToPrimLabel
                      !          WRITE(LUFILE,'(A)')'     iPrimD = (iPrimQ-1)/nPrimC+1'
                      WRITE(LUFILE,'(A)')'     IP = IP + 1'
                   ENDIF
                   nLines = nLines + 1
                ENDIF !COLLAPSE

                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'     invexpQ = D1/Qexp(iPrimQ)'
                ELSE
                   WRITE(LUFILE,'(A)')'     invexpQ = D1/Qexp(1)'
                ENDIF
                WRITE(LUFILE,'(A)')'     inv2expQ = D05*invexpQ'
                nLines = nLines + 2
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(8A)')'     facX = -(exp',FromExpLabel,'X',SIGN2,ToExpLabel,'exp(iPrim',ToExpLabel,')*Xcd)*invexpQ'
                   WRITE(LUFILE,'(8A)')'     facY = -(exp',FromExpLabel,'Y',SIGN2,ToExpLabel,'exp(iPrim',ToExpLabel,')*Ycd)*invexpQ'
                   WRITE(LUFILE,'(8A)')'     facZ = -(exp',FromExpLabel,'Z',SIGN2,ToExpLabel,'exp(iPrim',ToExpLabel,')*Zcd)*invexpQ'
                ELSE
                   WRITE(LUFILE,'(8A)')'     facX = -(exp',FromExpLabel,'X',SIGN2,ToExpLabel,'exp(1)*Xcd)*invexpQ'
                   WRITE(LUFILE,'(8A)')'     facY = -(exp',FromExpLabel,'Y',SIGN2,ToExpLabel,'exp(1)*Ycd)*invexpQ'
                   WRITE(LUFILE,'(8A)')'     facZ = -(exp',FromExpLabel,'Z',SIGN2,ToExpLabel,'exp(1)*Zcd)*invexpQ'
                ENDIF
                nLines = nLines + 3
                !          WRITE(LUFILE,'(A)')'     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpQ'
                !          WRITE(LUFILE,'(A)')'     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpQ'
                !          WRITE(LUFILE,'(A)')'     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpQ'
                WRITE(LUFILE,'(A)')'     pinvq = -expP*invexpQ'
                nLines = nLines + 1
                CALL SUBROUTINE_MAIN(LUFILE,JMAX,JP,JQ,nTUVP,nTUVQ,JPQ,nTUV,nTUVPLUS,TUVINDEX,&
                     & TINDEX,UINDEX,VINDEX,JINDEX,nTUVprev3,nTUVprev2,nTUVprev,IfacX,TUVindexX,LOOPUNROLL,&
                     & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPP,nTUVTMPP2)               

                IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ'
                WRITE(LUFILE,'(A,I3)')'     DO iTUVQ=1,',nTUVQ
                IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ'
                WRITE(LUFILE,'(A,I3)')'      DO iTUVP=1,',nTUVP
                IF(COLLAPSE)THEN
                   IF(nPrimLast)THEN
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,IP) = Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(Seg)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,IP) = Tmp0(iTUVP,iTUVQ)'
                      ENDIF
                   ELSE
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(Seg)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)'
                      ENDIF
                   ENDIF
                ELSE
                   IF(nPrimLast)THEN
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,IP) = Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iPrimP,iPassP) = Aux2(iTUVP,iTUVQ,iPrimP,iPassP) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iPrimQ,iPassP) = Aux2(iTUVP,iTUVQ,iPrimQ,iPassP) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(Seg)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iPassP) = Aux2(iTUVP,iTUVQ,iPassP) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iPassP) = Tmp0(iTUVP,iTUVQ)'
                      ENDIF
                   ELSE
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iPrimP,iPassP,iTUVP,iTUVQ) = Aux2(iPrimP,iPassP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iPrimQ,iPassP,iTUVP,iTUVQ) = Aux2(iPrimQ,iPassP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(Seg)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iPassP,iTUVP,iTUVQ) = Aux2(iPassP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)'
                      ELSEIF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iPassP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)'
                      ENDIF
                   ENDIF
                ENDIF
                nLines = nLines + 1
                WRITE(LUFILE,'(A)')   '      ENDDO'
                WRITE(LUFILE,'(A)')   '     ENDDO'
                nLines = nLines + 2
                IF(COLLAPSE)THEN
                   IF(Gen)THEN
                      WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPrimQ*nPrimP*nPasses'
                   ELSEIF(SegQ)THEN
                      WRITE(LUFILE,'(A)') '   ENDDO !iPrimP=1, nPrimP'       
                      WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPrimQ*nPasses'
                   ELSEIF(SegP)THEN
                      WRITE(LUFILE,'(A)') '   ENDDO !iPrimQ=1, nPrimQ'
                      WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPrimP*nPasses'
                   ELSEIF(Seg)THEN
                      WRITE(LUFILE,'(A)') '   ENDDO !iPrimQP = 1,nPrimQ*nPrimP'
                      WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPasses'
                   ELSEIF(Seg1Prim)THEN
                      WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPasses'
                   ENDIF
                ELSE
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    ENDDO'
                      WRITE(LUFILE,'(A)')'   ENDDO'
                   ENDIF
                   WRITE(LUFILE,'(A)')'  ENDDO'
                ENDIF
                nLines = nLines + 2
                IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
!                IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'

                call initString(1)          
                call AddToString('end subroutine TransferRecurrence'//ARCSTRING//'P')
                call AddToString(JP)
                call AddToString('Q')
                call AddToString(JQ)
                call AddToString(FromLabel)
                call AddToString('to')
                call AddToString(ToLabel)
                call AddToString(SegLabel(1:iSegLabel))
                call writeString(LUFILE); nLines = nLines + 1

!!$          IF(JP.LT.10)THEN
!!$             IF(JQ.LT.10)THEN
!!$                WRITE(LUFILE,'(A,I1,A,I1,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC'
!!$             ELSE
!!$                WRITE(LUFILE,'(A,I1,A,I2,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC'
!!$             ENDIF
!!$          ELSE
!!$             IF(JQ.LT.10)THEN
!!$                WRITE(LUFILE,'(A,I2,A,I1,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC'
!!$             ELSE
!!$                WRITE(LUFILE,'(A,I2,A,I2,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC'
!!$             ENDIF
!!$          ENDIF
                deallocate(TUVINDEX)
                deallocate(TINDEX)
                deallocate(UINDEX)
                deallocate(VINDEX)
                deallocate(JINDEX)
             enddo
          enddo
          WRITE(LUFILE,'(A)')'end module'
          close(unit = LUFILE)
       ENDDO
    ENDDO
 ENDDO
 WRITE(LUFILEMOD,'(7A)')'END MODULE AGC_CPU_OBS_TRParamMod'
 close(unit = LUFILEMOD)
END subroutine PASSsub

  subroutine SUBROUTINE_MAIN(LUFILE,JMAX,JP,JQ,nTUVP,nTUVQ,JPQ,nTUV,nTUVPLUS,TUVINDEX,&
               & TINDEX,UINDEX,VINDEX,JINDEX,nTUVprev3,nTUVprev2,nTUVprev,IfacX,&
               & TUVindexX,LOOPUNROLL,Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
    implicit none
    INTEGER,intent(in) :: LUFILE,JMAX,JP,JQ,nTUVP,nTUVQ,JPQ,nTUV,nTUVPLUS,nTUVTMPPX,nTUVTMPP2X
    integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    integer :: TINDEX(nTUVPLUS),IfacX(:,:),TUVindexX(:,:)
    integer :: UINDEX(nTUVPLUS)
    integer :: VINDEX(nTUVPLUS)
    integer :: JINDEX(nTUVPLUS)
    logical :: LOOPUNROLL,Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC
    INTEGER,intent(in) :: nTUVprev3,nTUVprev2,nTUVprev
    !local
    INTEGER :: ituvP,J,Tp,Up,Vp,N,N2,ituv,C
    INTEGER :: nTUVTMPP,nTUVTMPQ,nTUVTMP,nTUVTMPprev
    integer :: tq,uq,vq,ituvq,ituvqminus1x,ituvqminus1y,ituvqminus1z
    integer :: iTUVplus1x,iTUVplus1y,iTUVplus1z,nTUVTMPQs,I,nTUVTMPP2,nTUVP2
    Integer :: MaxAngmomQP,JTMQ,iTUVPminus1
    !    logical :: CREATED(-2:8,-2:8,-2:8)
    logical,pointer :: CREATED(:,:,:),ituvpplus1LEnTUVParray(:)
    logical,pointer :: ituvpminus1LEnTUVParray(:)
    logical :: TREC,UREC,VREC,TREC2,UREC2,VREC2,DoneCartDir(3),ituvpplus1LEnTUVP
    logical :: ituvpminus1LEnTUVP
    integer :: nTUVLIST,nTUVLISTactual,CARTDIR,iTUVQminus1,iTUVQminus2,Tqminus1
    integer,pointer :: TwoTermTUVLIST(:)
    integer :: iTUVPplus1,nLength
    allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
    CREATED  = .FALSE.
    CREATED(0,0,0) = .TRUE.
    
    DO JTMQ=0,JQ
       !           print*,'!JTMQ = ',JTMQ
       WRITE(LUFILE,'(A,I2)')   ' ! Building for Angular momentum Jq =',JTMQ
       nLines = nLines + 1
       IF(JTMQ.EQ.0)THEN
          IF(COLLAPSE)THEN
             IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ';nLines = nLines + 1
             IF(Gen.OR.Seg1Prim)THEN
                WRITE(LUFILE,'(A,I3)')'     DO iTUVP=1,',nTUVP;nLines = nLines + 1
                IF(nPrimLast)THEN
                   WRITE(LUFILE,'(A)')   '      Tmp0(iTUVP,1) = Aux(iTUVP,IP)';nLines = nLines + 1
                ELSE
                   WRITE(LUFILE,'(A)')   '      Tmp0(iTUVP,1) = Aux(IP,iTUVP)';nLines = nLines + 1
                ENDIF
                WRITE(LUFILE,'(A)')   '     ENDDO';nLines = nLines + 1
             ELSE
                WRITE(LUFILE,'(A,I3)')'     DO iTUVP=1,',nTUVP;nLines = nLines + 1
                IF(nPrimLast)THEN
                   WRITE(LUFILE,'(A)')   '      Tmp0(iTUVP,1) = Aux(iTUVP,iPrimQ,iPrimP,iPassP)';nLines = nLines + 1
                ELSE
                   WRITE(LUFILE,'(A)')   '      Tmp0(iTUVP,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVP)';nLines = nLines + 1
                ENDIF
                WRITE(LUFILE,'(A)')   '     ENDDO';nLines = nLines + 1
             ENDIF
          ELSE
             IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ';nLines = nLines + 1
             WRITE(LUFILE,'(A,I3)')'     DO iTUVP=1,',nTUVP;nLines = nLines + 1
             IF(nPrimLast)THEN
                WRITE(LUFILE,'(A)')   '      Tmp0(iTUVP,1) = Aux(iTUVP,IP)';nLines = nLines + 1
             ELSE
                WRITE(LUFILE,'(A)')   '      Tmp0(iTUVP,1) = Aux(IP,iTUVP)';nLines = nLines + 1
             ENDIF
             WRITE(LUFILE,'(A)')   '     ENDDO';nLines = nLines + 1
          ENDIF
          CYCLE
       ELSE
          !==========================================================================
          !            Build the     Theta(i,0,k,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k-1,0) + (k-1)/(2q)*Theta(i,0,k-2,0)
          !==========================================================================
          nTUVTMPQ=(JTMQ+1)*(JTMQ+2)*(JTMQ+3)/6
          nTUVTMPP=(JPQ-JTMQ+1)*(JPQ-JTMQ+2)*(JPQ-JTMQ+3)/6
          !slowly increase JQ=JTMQ with 1 
          nTUVTMPQ=(JTMQ+1)*(JTMQ+2)*(JTMQ+3)/6
          DO Tq=JTMQ,0,-1       
             DO Uq=JTMQ-Tq,0,-1
                Vq=JTMQ-Tq-Uq  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1,Tq,Uq,Vq,CREATED,JMAX,TUVINDEX)
                IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ';nLines = nLines + 1
                WRITE(LUFILE,'(A,I3)')'     do iTUVP = 1,',nTUVP;nLines = nLines + 1
                !Theta(i,0,k,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k-1,0) + (k-1)/(2q)*Theta(i,0,k-2,0)
                CALL LOOPRECURRENCE1(iTUVQ,Tqminus1,iTUVQminus2,iTUVQminus1,CARTDIR,LUFILE,&
                     & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,nTUVTMPPX,nTUVTMPP2X)
                WRITE(LUFILE,'(A,I3)')'     enddo';nLines = nLines + 1
                
                !                   CREATED(Tq,Uq,Vq) = .TRUE.
             ENDDO
          ENDDO
          
          DO Tq=JTMQ,0,-1       
             DO Uq=JTMQ-Tq,0,-1
                Vq=JTMQ-Tq-Uq  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1,Tq,Uq,Vq,CREATED,JMAX,TUVINDEX)
                
                IF(nTUVTMPP.GE.nTUVP+1)THEN
                   IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ';nLines = nLines + 1
                   WRITE(LUFILE,'(A,I3,A,I3)')'     do iTUVP = ',nTUVP+1,',',nTUVTMPP !place in tmp array 
                   nLines = nLines + 1
                   !Theta(i,0,k,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k-1,0) + (k-1)/(2q)*Theta(i,0,k-2,0)
                   CALL LOOPRECURRENCE2(iTUVQ,Tqminus1,iTUVQminus2,iTUVQminus1,CARTDIR,JTMQ,LUFILE,&
                        & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,nTUVTMPPX,nTUVTMPP2X)
                   WRITE(LUFILE,'(A,I3,A,I3)')'     enddo';nLines = nLines + 1
                ENDIF
                
                !                   CREATED(Tq,Uq,Vq) = .TRUE.
             ENDDO
          ENDDO
          !==========================================================================
          !            Build the     Theta(i,0,k,0) += i/(2q)*Theta(i-1,0,k-1,0) 
          !==========================================================================
          DO Tq=JTMQ,0,-1       
             DO Uq=JTMQ-Tq,0,-1
                Vq=JTMQ-Tq-Uq  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1,Tq,Uq,Vq,CREATED,JMAX,TUVINDEX)
                nTUVP2 = (JP+0)*(JP+1)*(JP+2)/6   
                iTUVPminus1 = 1
                ituvpminus1LEnTUVP = .TRUE. 
              !  print*,'WRITERECURRENCE51'
                CALL WRITERECURRENCE5(CARTDIR,Tq,Uq,Vq,1,nTUVP2,&
                     & nTUVP2,TUVINDEX,JTMQ,JMAX,LUFILE,iTUVPminus1,.TRUE.,ituvpminus1LEnTUVP,IfacX,TUVindexX,LOOPUNROLL,&
                     & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
              !  print*,'DONE WRITERECURRENCE51'
                
                !                   do iTUVP = 1,nTUVP
                !                      Tp = Tindex(iTUVp) 
                !                      Up = Uindex(iTUVp) 
                !                      Vp = Vindex(iTUVp)
                !                      !Theta(i,0,k,0) += i/(2q)*Theta(i-1,0,k-1,0) 
                !                      CALL WRITERECURRENCE1(CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE)
                !                   enddo                   
             ENDDO
          ENDDO
          
          DO Tq=JTMQ,0,-1       
             DO Uq=JTMQ-Tq,0,-1
                Vq=JTMQ-Tq-Uq  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1,Tq,Uq,Vq,CREATED,JMAX,TUVINDEX)
                IF(nTUVTMPP.GE.nTUVP+1)THEN
                   nTUVP2 = (JP+0)*(JP+1)*(JP+2)/6   
                   nTUVTMPP2=(JPQ-JTMQ)*(JPQ-JTMQ+1)*(JPQ-JTMQ+2)/6
                   IF(iTUVQminus1.EQ.1)THEN
                      iTUVPminus1 = 1
                      ituvpminus1LEnTUVP = .TRUE. !dummy
                !print*,'WRITERECURRENCE52'
                      CALL WRITERECURRENCE5(CARTDIR,Tq,Uq,Vq,nTUVP2+1,nTUVTMPP2,&
                           & nTUVP2,TUVINDEX,JTMQ,JMAX,LUFILE,iTUVPminus1,.FALSE.,ituvpminus1LEnTUVP,IfacX,TUVindexX,LOOPUNROLL,&
                           & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE52'
                   ELSE
                      IF(nTUVTMPP2.LE.nTUVP)THEN
                         ituvpminus1LEnTUVP = .TRUE.
                !print*,'WRITERECURRENCE53'
                         CALL WRITERECURRENCE5(CARTDIR,Tq,Uq,Vq,nTUVP2+1,nTUVTMPP2,&
                              & nTUVP2,TUVINDEX,JTMQ,JMAX,LUFILE,iTUVPminus1,.FALSE.,ituvpminus1LEnTUVP,IfacX,TUVindexX,LOOPUNROLL,&
                              Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE53'
                      ELSE
                !print*,'WRITERECURRENCE54'
                         ituvpminus1LEnTUVP = .TRUE.
                         CALL WRITERECURRENCE5(CARTDIR,Tq,Uq,Vq,nTUVP2+1,nTUVP,&
                              & nTUVP2,TUVINDEX,JTMQ,JMAX,LUFILE,iTUVPminus1,.FALSE.,ituvpminus1LEnTUVP,IfacX,TUVindexX,LOOPUNROLL,&
                              & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE54'

                         ituvpminus1LEnTUVP = .FALSE.
                !print*,'WRITERECURRENCE55'
                         CALL WRITERECURRENCE5(CARTDIR,Tq,Uq,Vq,nTUVP+1,nTUVTMPP2,&
                              & nTUVP2,TUVINDEX,JTMQ,JMAX,LUFILE,iTUVPminus1,.FALSE.,ituvpminus1LEnTUVP,IfacX,TUVindexX,LOOPUNROLL,&
                              & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE55'
                      ENDIF
                   ENDIF
                   !                      do iTUVP = nTUVP+1,nTUVTMPP
                   !                         Tp = Tindex(iTUVp) 
                   !                         Up = Uindex(iTUVp) 
                   !                         Vp = Vindex(iTUVp)
                   !                         !Theta(i,0,k,0) += i/(2q)*Theta(i-1,0,k-1,0)
                   !                         CALL WRITERECURRENCE1(CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE)
                   !                      enddo
                   
                ENDIF
                
             ENDDO
          ENDDO
          !==========================================================================
          !            Build the     Theta(i,0,k,0) += - p/q*Theta(i+1,0,k-1,0) 
          !==========================================================================
          DO Tq=JTMQ,0,-1       
             DO Uq=JTMQ-Tq,0,-1
                Vq=JTMQ-Tq-Uq  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1,Tq,Uq,Vq,CREATED,JMAX,TUVINDEX)
                
                !Theta(i,0,k,0) += - p/q*Theta(i+1,0,k-1,0) 
                
                IF(iTUVQminus1.EQ.1)THEN
                   ituvpplus1LEnTUVP = .TRUE. !dummy argument 
                !print*,'WRITERECURRENCE4'
                   CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,1,ituvpplus1LEnTUVP,nTUVP,nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                        & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                ELSE
                   allocate(ituvpplus1LEnTUVParray(nTUVP))
                   DO iTUVP = 1, nTUVP
                      Tp = Tindex(iTUVp) 
                      Up = Uindex(iTUVp) 
                      Vp = Vindex(iTUVp)
                      IF(CARTDIR.EQ.1)iTUVPplus1 = TUVINDEX(Tp+1,Up,Vp)
                      IF(CARTDIR.EQ.2)iTUVPplus1 = TUVINDEX(Tp,Up+1,Vp)
                      IF(CARTDIR.EQ.3)iTUVPplus1 = TUVINDEX(Tp,Up,Vp+1)
                      ituvpplus1LEnTUVParray(iTUVP) = ituvpplus1.LE.nTUVP
                   ENDDO
                   IF(ALL(ituvpplus1LEnTUVParray))THEN
                      !all true
                      ituvpplus1LEnTUVP = .TRUE.
                !print*,'WRITERECURRENCE4'
                      CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,1,ituvpplus1LEnTUVP,nTUVP,&
                           & nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                           & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                   ELSEIF(COUNT(ituvpplus1LEnTUVParray).EQ.0)THEN
                      !all false
                      ituvpplus1LEnTUVP = .FALSE.
                !print*,'WRITERECURRENCE4'
                      CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,1,ituvpplus1LEnTUVP,nTUVP,&
                           & nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                           & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                   ELSE                      
                      IF(ALL(ituvpplus1LEnTUVParray(1:COUNT(ituvpplus1LEnTUVParray))))THEN
                         !all T are sequential in the beginning
                         ituvpplus1LEnTUVP = .TRUE.
                !print*,'WRITERECURRENCE4'
                         CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,1,ituvpplus1LEnTUVP,COUNT(ituvpplus1LEnTUVParray),&
                              & nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                              & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                         !all F are sequential at the end
                         ituvpplus1LEnTUVP = .FALSE.
                !print*,'WRITERECURRENCE4'
                         CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,COUNT(ituvpplus1LEnTUVParray)+1,ituvpplus1LEnTUVP,nTUVP,&
                              & nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                              & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                      ELSE
                         print*,'ituvpplus1LEnTUVParray',ituvpplus1LEnTUVParray
                         stop 'ERROR'
                      ENDIF
                   ENDIF
                   deallocate(ituvpplus1LEnTUVParray)
                ENDIF
                
                !                   do iTUVP = 1,nTUVP
                !                      Tp = Tindex(iTUVp) 
                !                      Up = Uindex(iTUVp) 
                !                      Vp = Vindex(iTUVp)
                !                      !Theta(i,0,k,0) += - p/q*Theta(i+1,0,k-1,0) 
                !                      CALL WRITERECURRENCE2(CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE)
                !                   enddo
                
                IF(nTUVTMPP.GE.nTUVP+1)THEN
                   IF(iTUVQminus1.EQ.1)THEN
                      ituvpplus1LEnTUVP = .TRUE. !dummy argument 
                !print*,'WRITERECURRENCE4'
                      CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,nTUVP+1,ituvpplus1LEnTUVP,nTUVTMPP,nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                           & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                   ELSE
                      !Theta(i,0,k,0) += - p/q*Theta(i+1,0,k-1,0) 
                      allocate(ituvpplus1LEnTUVParray(nTUVP+1:nTUVTMPP))
                      DO iTUVP = nTUVP+1, nTUVTMPP
                         Tp = Tindex(iTUVp) 
                         Up = Uindex(iTUVp) 
                         Vp = Vindex(iTUVp)
                         IF(CARTDIR.EQ.1)iTUVPplus1 = TUVINDEX(Tp+1,Up,Vp)
                         IF(CARTDIR.EQ.2)iTUVPplus1 = TUVINDEX(Tp,Up+1,Vp)
                         IF(CARTDIR.EQ.3)iTUVPplus1 = TUVINDEX(Tp,Up,Vp+1)
                         ituvpplus1LEnTUVParray(iTUVP) = ituvpplus1.LE.nTUVP
                      ENDDO
                      IF(ALL(ituvpplus1LEnTUVParray))THEN
                         !all true
                         ituvpplus1LEnTUVP = .TRUE.
                !print*,'DONE WRITERECURRENCE4'
                         CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,nTUVP+1,ituvpplus1LEnTUVP,nTUVTMPP,&
                              & nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                              & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'WRITERECURRENCE4'
                      ELSEIF(COUNT(ituvpplus1LEnTUVParray).EQ.0)THEN
                         !all false
                         ituvpplus1LEnTUVP = .FALSE.
                !print*,'WRITERECURRENCE4'
                         CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,nTUVP+1,ituvpplus1LEnTUVP,nTUVTMPP,&
                              & nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                              & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                      ELSE
                         IF(ALL(ituvpplus1LEnTUVParray(nTUVP+1:nTUVP+COUNT(ituvpplus1LEnTUVParray))))THEN
                            !all T are sequential in the beginning
                            ituvpplus1LEnTUVP = .TRUE.
                !print*,'WRITERECURRENCE4'
                            CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,nTUVP+1,ituvpplus1LEnTUVP,nTUVP+COUNT(ituvpplus1LEnTUVParray),&
                                 & nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                                 & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                            !all F are sequential at the end
                            ituvpplus1LEnTUVP = .FALSE.
                !print*,'WRITERECURRENCE4'
                            CALL WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,nTUVP+COUNT(ituvpplus1LEnTUVParray)+1,ituvpplus1LEnTUVP,nTUVTMPP,&
                                 & nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE,TUVindexX,LOOPUNROLL,&
                                 & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPPX,nTUVTMPP2X)
                !print*,'DONE WRITERECURRENCE4'
                         ELSE
                            print*,'ituvpplus1LEnTUVParray',ituvpplus1LEnTUVParray
                            stop 'ERROR2'
                         ENDIF
                      ENDIF
                      deallocate(ituvpplus1LEnTUVParray)
                   ENDIF
                   !                      do iTUVP = nTUVP+1,nTUVTMPP
                   !                         Tp = Tindex(iTUVp) 
                   !                         Up = Uindex(iTUVp) 
                   !                         Vp = Vindex(iTUVp)
                   !                         !Theta(i,0,k,0) = i/(2q)*Theta(i-1,0,k-1,0) - p/q*Theta(i+1,0,k-1,0) 
                   !                         CALL WRITERECURRENCE2(CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,TUVINDEX,JTMQ,JMAX,LUFILE)
                   !                      enddo
                ENDIF
                
                CREATED(Tq,Uq,Vq) = .TRUE.
             ENDDO
          ENDDO
       ENDIF
    ENDDO
  end subroutine SUBROUTINE_MAIN

  subroutine DETERMINE_CARTDIR(CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1,Tq,Uq,Vq,CREATED,JMAX,TUVINDEX)
    implicit none
    integer,intent(inout) :: CARTDIR,iTUVQ,iTUVQminus1,iTUVQminus2,Tqminus1
    integer,intent(in) :: Tq,Uq,Vq,JMAX
    logical,intent(in) :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    !
    integer :: N,N2
    logical :: Trec,Urec,Vrec,Trec2,Urec2,Vrec2
    iTUVQ = TUVINDEX(Tq,Uq,Vq)
    !Theta(i,0,k,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k-1,0) + i/(2q)*Theta(i-1,0,k-1,0) + (k-1)/(2q)*Theta(i,0,k-2,0) - p/q*Theta(i+1,0,k-1,0) 
    
    !choose to use X, Y or Z 
    !how can the (Tp,Up,Vp,Tq,Uq,Vq) be built from lower aka (210) can be built from (110) and (200)
    

    
    !how can the (Tp,Up,Vp) be built from lower aka (210) can be built from (110) and (200)
    !However building (210) from (110) using the X recurrence requires both (110) and (010)   (210) = (110) + 1*(010)
    !building (210) from (200) using the Y recurrence requires both (200)   (210) = (200) + 0*(2-10)
    TREC = CREATED(Tq-1,Uq,Vq)     !Example(210): test if (110) is build  (TRUE) 
    UREC = CREATED(Tq,Uq-1,Vq)     !Example(210): test if (200) is build  (TRUE) 
    VREC = CREATED(Tq,Uq,Vq-1)     !Example(210): test if (21-1) is build (FALSE) 
    N=0
    IF(TREC)N=N+1                  !Example(210): N=1
    IF(UREC)N=N+1                  !Example(210): N=2
    IF(VREC)N=N+1                  !Example(210): N=2
    IF(N.EQ.1)THEN
       !only one possible way to construct it like (200) can only be constructed from (100)
       IF(TREC)THEN
          CARTDIR = 1; iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq); iTUVQminus2 = TUVINDEX(Tq-2,Uq,Vq); Tqminus1=Tq-1
       ELSEIF(UREC)THEN
          CARTDIR = 2; iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq); iTUVQminus2 = TUVINDEX(Tq,Uq-2,Vq); Tqminus1=Uq-1
       ELSEIF(VREC)THEN
          CARTDIR = 3; iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1); iTUVQminus2 = TUVINDEX(Tq,Uq,Vq-2); Tqminus1=Vq-1          
       ELSE
          STOP 'TK1'
       ENDIF
    ELSE
       !several ways to construct it. For instance the (210) Example
       TREC2 = CREATED(Tq-2,Uq,Vq)    !Example(210): test if (010) is build (TRUE) 
       UREC2 = CREATED(Tq,Uq-2,Vq)    !Example(210): test if (2-10) is build (FALSE) 
       VREC2 = CREATED(Tq,Uq,Vq-2)    !Example(210): test if (21-1) is build (FALSE) 
       N2=0
       IF(TREC2)N2=N2+1               !Example(210): N2 = 1
       IF(UREC2)N2=N2+1               !Example(210): N2 = 1
       IF(VREC2)N2=N2+1               !Example(210): N2 = 1
       IF(N2.LT.N)THEN                !Example(210): N2 = 1 .LT. N = 2 => TRUE
          !3 term recurrence possible for one or more the possibilities 
          IF(.NOT.(TREC.AND.TREC2).AND.TREC)THEN       !Example(210): TREC.AND.TREC2 = .TRUE. => FALSE
             CARTDIR = 1; iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq); iTUVQminus2 = TUVINDEX(Tq-2,Uq,Vq); Tqminus1=Tq-1             
          ELSEIF(.NOT.(UREC.AND.UREC2).AND.UREC)THEN   !Example(210): UREC.AND.UREC2 = .FALSE. => TRUE
             CARTDIR = 2; iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq); iTUVQminus2 = TUVINDEX(Tq,Uq-2,Vq); Tqminus1=Uq-1          
          ELSEIF(.NOT.(VREC.AND.VREC2).AND.VREC)THEN
             CARTDIR = 3; iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1); iTUVQminus2 = TUVINDEX(Tq,Uq,Vq-2); Tqminus1=Vq-1                     
          ENDIF
       ELSE
          !all requires 4 term recurrence possible for one or more the possibilities 
          IF(TREC)THEN
             CARTDIR = 1; iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq); iTUVQminus2 = TUVINDEX(Tq-2,Uq,Vq); Tqminus1=Tq-1             
          ELSEIF(UREC)THEN
             CARTDIR = 2; iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq); iTUVQminus2 = TUVINDEX(Tq,Uq-2,Vq); Tqminus1=Uq-1          
          ELSEIF(VREC)THEN
             CARTDIR = 3; iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1); iTUVQminus2 = TUVINDEX(Tq,Uq,Vq-2); Tqminus1=Vq-1          
          ELSE
             STOP 'TK2'
          ENDIF
       ENDIF
    ENDIF
  END subroutine DETERMINE_CARTDIR

  SUBROUTINE LOOPRECURRENCE1(iTUVQ,Tqminus1,iTUVQminus2,iTUVQminus1,CARTDIR,LUPRI,&
       & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,nTUVTMPP,nTUVTMPP2)
    implicit none
    integer :: iTUVQ,Tqminus1,iTUVQminus2,iTUVQminus1,CARTDIR,lupri,nTUVTMPP,nTUVTMPP2
    character(len=132) :: STRING 
    integer :: iString
    character(len=4) :: DIRECTIONSTRING
    logical :: Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE
    STRING(1:8) = '        '
    iSTRING = 9
    IF(CARTDIR.EQ.1)THEN
       DIRECTIONSTRING = 'facX'
    ELSEIF(CARTDIR.EQ.2)THEN
       DIRECTIONSTRING = 'facY'
    ELSE
       DIRECTIONSTRING = 'facZ'
    ENDIF
    call initString(6)
    call AddToString('Tmp0(iTUVP,')
    call AddToString(iTUVQ)
    call AddToString(') = ')
    call AddToString(DIRECTIONSTRING)
    call AddToString('*')
!    WRITE(STRING(iSTRING:iSTRING+26),'(A11,I4,A7,A4,A1)') 'Aux2(iTUVP,',iTUVQ,',IP) = ',DIRECTIONSTRING,'*'
!    iString = iSTRING+27
    IF(iTUVQminus1.EQ.1)THEN
       IF(COLLAPSE)THEN
          IF(Gen.OR.Seg1Prim)THEN
             IF(nPrimLast)THEN
                call AddToString('Aux(iTUVP,IP)')             
             ELSE
                call AddToString('Aux(IP,iTUVP)')             
             ENDIF
          ELSE
             IF(nPrimLast)THEN
                call AddToString('Aux(iTUVP,iPrimQ,iPrimP,iPassP)')
             ELSE
                call AddToString('Aux(iPrimQ,iPrimP,iPassP,iTUVP)')
             ENDIF
          ENDIF
       ELSE
          IF(nPrimLast)THEN
             call AddToString('Aux(iTUVP,IP)')
          ELSE
             call AddToString('Aux(IP,iTUVP)')
          ENDIF
       ENDIF
!       WRITE(STRING(iSTRING:iSTRING+12),'(A13)') 'Aux(iTUVP,IP)'
!       iString = iSTRING+13
    ELSE
       call AddToString('Tmp0(iTUVP,')
       call AddToString(iTUVQminus1)
       call AddToString(')')
!       WRITE(STRING(iSTRING:iSTRING+18),'(A11,I4,A4)') 'Aux2(iTUVP,',iTUVQminus1,',IP)'
!       iString = iSTRING+19
    ENDIF
    !possibly add term 3   
    IF(Tqminus1.EQ.1.AND.iTUVQminus2.GT.0)THEN
       IF(iTUVQminus2.EQ.1)THEN
          IF(nPrimLast)THEN
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('+ inv2expQ*Aux(iTUVP,IP)')
                ELSE
                   call AddToString('+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)')
                ENDIF
             ELSE
                call AddToString('+ inv2expQ*Aux(iTUVP,IP)')
             ENDIF
          ELSE
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('+ inv2expQ*Aux(IP,iTUVP)')
                ELSE
                   call AddToString('+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)')
                ENDIF
             ELSE
                call AddToString('+ inv2expQ*Aux(IP,iTUVP)')
             ENDIF
          ENDIF
!          WRITE(STRING(iSTRING:iSTRING+23),'(A24)') '+ inv2expQ*Aux(iTUVP,IP)'
!          iString = iSTRING+24
       ELSE
          call AddToString('+ inv2expQ*Tmp0(iTUVP,')
          call AddToString(iTUVQminus2)
          call AddToString(')')
!          WRITE(STRING(iSTRING:iSTRING+29),'(A22,I4,A4)') '+ inv2expQ*Aux2(iTUVP,',iTUVQminus2,',IP)'
!          iString = iSTRING+30
       ENDIF
    ELSEIF(Tqminus1.GT.0.AND.iTUVQminus2.GT.0)THEN
       IF(iTUVQminus2.EQ.1)THEN
          call AddToString('+')
          call AddToString(Tqminus1)
          IF(nPrimLast)THEN
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('*inv2expQ*Aux(iTUVP,IP)')
                ELSE
                   call AddToString('*inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)')
                ENDIF
             ELSE
                call AddToString('*inv2expQ*Aux(iTUVP,IP)')
             ENDIF
          ELSE
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('*inv2expQ*Aux(IP,iTUVP)')
                ELSE
                   call AddToString('*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)')
                ENDIF
             ELSE
                call AddToString('*inv2expQ*Aux(IP,iTUVP)')
             ENDIF
          ENDIF
!          WRITE(STRING(iSTRING:iSTRING+24),'(A1,I1,A23)') '+',Tqminus1,'*inv2expQ*Aux(iTUVP,IP)'
!          iString = iSTRING+25
       ELSE
          call AddToString('+')
          call AddToString(Tqminus1)
          call AddToString('*inv2expQ*Tmp0(iTUVP,')
          call AddToString(iTUVQminus2)
          call AddToString(')')
!          call AddToString(',IP)')
!          WRITE(STRING(iSTRING:iSTRING+30),'(A1,I1,A21,I4,A4)') '+',Tqminus1,'*inv2expQ*Aux2(iTUVP,',iTUVQminus2,',IP)'
!          iString = iSTRING+31
       ENDIF
    ENDIF
    call writeString(LUPRI);nLines = nLines + 1
    !Final step write the string
!    WRITE(LUFILE,'(A)') STRING(1:iSTRING-1)
  END SUBROUTINE LOOPRECURRENCE1
  
  SUBROUTINE LOOPRECURRENCE2(iTUVQ,Tqminus1,iTUVQminus2,iTUVQminus1,CARTDIR,JTMQ,LUPRI,&
       & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,nTUVTMPP,nTUVTMPP2)
    implicit none
    integer :: iTUVQ,Tqminus1,iTUVQminus2,iTUVQminus1,CARTDIR,JTMQ,LUPRI,nTUVTMPP,nTUVTMPP2
    character(len=132) :: STRING 
    integer :: iString
    character(len=4) :: DIRECTIONSTRING
    logical :: Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE
    STRING(1:8) = '        '
    iSTRING = 9
    IF(CARTDIR.EQ.1)THEN
       DIRECTIONSTRING = 'facX'
    ELSEIF(CARTDIR.EQ.2)THEN
       DIRECTIONSTRING = 'facY'
    ELSE
       DIRECTIONSTRING = 'facZ'
    ENDIF
    call initString(6)

    call AddToString('Tmp')
    call AddToString(JTMQ)
    call AddToString('(iTUVP,')
    call AddToString(iTUVQ)
    call AddToString(') = ')
    call AddToString(DIRECTIONSTRING)
    call AddToString('*')
!    IF(JTMQ.LT.10)THEN
!       WRITE(STRING(iSTRING:iSTRING+23),'(A3,I1,A7,I4,A4,A4,A1)') 'Tmp',JTMQ,'(iTUVP,',iTUVQ,') = ',DIRECTIONSTRING,'*'
!       iString = iSTRING+24
!    ELSE
!       WRITE(STRING(iSTRING:iSTRING+24),'(A3,I2,A7,I4,A4,A4,A1)') 'Tmp',JTMQ,'(iTUVP,',iTUVQ,') = ',DIRECTIONSTRING,'*'
!       iString = iSTRING+25
!    ENDIF
    IF(iTUVQminus1.EQ.1)THEN
       IF(nPrimLast)THEN
          IF(COLLAPSE)THEN
             IF(Gen.OR.Seg1Prim)THEN
                call AddToString('Aux(iTUVP,IP)')
             ELSE
                call AddToString('Aux(iTUVP,iPrimQ,iPrimP,iPassP)')
             ENDIF
          ELSE
             call AddToString('Aux(iTUVP,IP)')
          ENDIF
       ELSE
          IF(COLLAPSE)THEN
             IF(Gen.OR.Seg1Prim)THEN
                call AddToString('Aux(IP,iTUVP)')
             ELSE
                call AddToString('Aux(iPrimQ,iPrimP,iPassP,iTUVP)')
             ENDIF
          ELSE
             call AddToString('Aux(IP,iTUVP)')
          ENDIF
       ENDIF
!       WRITE(STRING(iSTRING:iSTRING+13),'(A14)') 'Aux(iTUVP,IP)'
!       iString = iSTRING+14
    ELSE
       call AddToString('Tmp')
       call AddToString(JTMQ-1)
       call AddToString('(iTUVP,')
       call AddToString(iTUVQminus1)
       call AddToString(')')
!       IF(JTMQ-1.LT.10)THEN
!          WRITE(STRING(iSTRING:iSTRING+15),'(A3,I1,A7,I4,A1)') 'Tmp',JTMQ-1,'(iTUVP,',iTUVQminus1,')'
!          iString = iSTRING+16
!       ELSE
!          WRITE(STRING(iSTRING:iSTRING+16),'(A3,I2,A7,I4,A1)') 'Tmp',JTMQ-1,'(iTUVP,',iTUVQminus1,')'
!          iString = iSTRING+17
!       ENDIF
    ENDIF
    !possibly add term 3   
    IF(Tqminus1.EQ.1.AND.iTUVQminus2.GT.0)THEN
       IF(iTUVQminus2.EQ.1)THEN
          IF(nPrimLast)THEN
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('+ inv2expQ*Aux(iTUVP,IP)') 
                ELSE
                   call AddToString('+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)') 
                ENDIF
             ELSE             
                call AddToString('+ inv2expQ*Aux(iTUVP,IP)') 
             ENDIF
          ELSE
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('+ inv2expQ*Aux(IP,iTUVP)') 
                ELSE
                   call AddToString('+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)') 
                ENDIF
             ELSE             
                call AddToString('+ inv2expQ*Aux(IP,iTUVP)') 
             ENDIF
          ENDIF
!           WRITE(STRING(iSTRING:iSTRING+23),'(A24)') '+ inv2expQ*Aux(iTUVP,IP)'
!          iString = iSTRING+24
       ELSE
          call AddToString('+ inv2expQ*Tmp') 
          call AddToString(JTMQ-2) 
          call AddToString('(iTUVP,') 
          call AddToString(iTUVQminus2) 
          call AddToString(')') 
!          IF(JTMQ-1.LT.10)THEN
!             WRITE(STRING(iSTRING:iSTRING+26),'(A14,I1,A7,I4,A1)') '+ inv2expQ*Tmp',JTMQ-2,'(iTUVP,',iTUVQminus2,')'
!             iString = iSTRING+27
!          ELSE
!             WRITE(STRING(iSTRING:iSTRING+27),'(A14,I2,A7,I4,A1)') '+ inv2expQ*Tmp',JTMQ-2,'(iTUVP,',iTUVQminus2,')'
!             iString = iSTRING+28
!          ENDIF
       ENDIF
    ELSEIF(Tqminus1.GT.0.AND.iTUVQminus2.GT.0)THEN
       IF(iTUVQminus2.EQ.1)THEN
          call AddToString('+') 
          call AddToString(Tqminus1) 
          IF(nPrimLast)THEN
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('*inv2expQ*Aux(iTUVP,IP)') 
                ELSE
                   call AddToString('*inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)') 
                ENDIF
             ELSE             
                call AddToString('*inv2expQ*Aux(iTUVP,IP)') 
             ENDIF
          ELSE
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('*inv2expQ*Aux(IP,iTUVP)') 
                ELSE
                   call AddToString('*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)') 
                ENDIF
             ELSE             
                call AddToString('*inv2expQ*Aux(IP,iTUVP)') 
             ENDIF
          ENDIF
!          WRITE(STRING(iSTRING:iSTRING+24),'(A1,I1,A23)') '+',Tqminus1,'*inv2expQ*Aux(iTUVP,IP)'
!          iString = iSTRING+25
       ELSE
          call AddToString('+') 
          call AddToString(Tqminus1) 
          call AddToString('*inv2expQ*Tmp') 
          call AddToString(JTMQ-2) 
          call AddToString('(iTUVP,') 
          call AddToString(iTUVQminus2) 
          call AddToString(')') 
!          IF(JTMQ-1.LT.10)THEN
!             WRITE(STRING(iSTRING:iSTRING+27),'(A1,I1,A13,I1,A7,I4,A1)') '+',Tqminus1,'*inv2expQ*Tmp',JTMQ-2,'(iTUVP,',iTUVQminus2,')'
!             iString = iSTRING+28
!          ELSE
!             WRITE(STRING(iSTRING:iSTRING+28),'(A1,I1,A13,I2,A7,I4,A1)') '+',Tqminus1,'*inv2expQ*Tmp',JTMQ-2,'(iTUVP,',iTUVQminus2,')'
!             iString = iSTRING+29
!          ENDIF
       ENDIF
    ENDIF
    !Final step write the string
    call writeString(LUPRI);nLines = nLines + 1
!    WRITE(LUFILE,'(A)') STRING(1:iSTRING-1)
  END SUBROUTINE LOOPRECURRENCE2

  subroutine WRITERECURRENCE1(CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,TUVINDEX,JTMQ,JMAX,LUPRI,&
       & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE)
    implicit none
    !A to C
    !Theta(i,0,k,0) = i/(2q)*Theta(i-1,0,k-1,0) - p/q*Theta(i+1,0,k-1,0) 
    integer,intent(in) :: CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,JTMQ,JMAX
    integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),LUPRI
    logical :: Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE
    !
    integer :: iTUVQ,iTUVP,iTUVPplus1,iTUVPminus1,I,iTUVQminus1
    character(len=132) :: STRING 
    integer :: iString
    iTUVP = TUVINDEX(Tp,Up,Vp)
    iTUVQ = TUVINDEX(Tq,Uq,Vq)
    IF(CARTDIR.EQ.1)THEN
       !X direction
       iTUVPplus1 = TUVINDEX(Tp+1,Up,Vp)
       iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp)
       I = Tp
       iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq)
    ELSEIF(CARTDIR.EQ.2)THEN
       !Y direction
       iTUVPplus1 = TUVINDEX(Tp,Up+1,Vp)
       iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp)
       I = Up
       iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq)                        
    ELSEIF(CARTDIR.EQ.3)THEN
       !Z direction
       iTUVPplus1 = TUVINDEX(Tp,Up,Vp+1)
       iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1)
       I = Vp
       iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1)
    ENDIF
    !Result(iTUVP,iTUVQ) = I * Term2(iTUVPminus1,iTUVQminus1) + Term4(iTUVPplus1,iTUVQminus1)

    !step 1 add blanks
    IF(I.GT.0)THEN
       call initString(5)
       !step 2 determine where to put the result
       IF(iTUVP.LE.nTUVP)THEN
          call AddToString('Tmp0(')
       ELSE
          call AddToString('Tmp')
          call AddToString(JTMQ)
          call AddToString('(')
       ENDIF
       call AddToString(iTUVP)
       call AddToString(',')
       call AddToString(iTUVQ)
       call AddToString(') = ')

       IF(iTUVP.LE.nTUVP)THEN
          call AddToString('Tmp0(')
       ELSE
          call AddToString('Tmp')
          call AddToString(JTMQ)
          call AddToString('(')
       ENDIF
       call AddToString(iTUVP)
       call AddToString(',')
       call AddToString(iTUVQ)
       call AddToString(') ')

       !=====================
       !step 4 determine if the second term: 
       !  I*inv2expQ*Aux(',ituvpminus1x,',IP) 
       !  I*inv2expQ*Tmp',JTMQ-1,'(',ituvpminus1x,',',iTUVQminus1x,')
       !should be included and if it should use Aux or Tmp 
       IF(I.GT.2)THEN
          call AddToString('+ ')
          call AddToString(I)
          call AddToString('*inv2expQ*')
       ELSEIF(I.EQ.2)THEN
          call AddToString('+ invexpQ*')
       ELSEIF(I.EQ.1)THEN
          call AddToString('+ inv2expQ*')
       ELSE
          !do not include this term 
       ENDIF
       IF(iTUVQminus1.EQ.1)THEN
          call AddToString('Aux(')
          IF(nPrimLast)THEN
             call AddToString(ituvpminus1)
             call AddToString(',')
          ENDIF
          IF(COLLAPSE)THEN
             IF(Gen.OR.Seg1Prim)THEN
                call AddToString('IP')
             ELSE
                !iPrimQ,iPrimP,iPassP
                call AddToString('IPrimQ,iPrimP,iPassP')
             ENDIF
          ELSE             
             call AddToString('IP')
          ENDIF
          IF(nPrimLast)THEN
             call AddToString(') ')
          ELSE
             call AddToString(',')             
             call AddToString(ituvpminus1)
             call AddToString(') ')
          ENDIF
       ELSE
          IF(ituvpminus1.LE.nTUVP)THEN
             call AddToString('Tmp0(')
          ELSE
             call AddToString('Tmp')
             call AddToString(JTMQ-1)
             call AddToString('(')
          ENDIF
          call AddToString(ituvpminus1)
          call AddToString(',')
          IF(ituvpminus1.LE.nTUVP)THEN
             call AddToString(iTUVQminus1)
             call AddToString(') ')
          ELSE
             call AddToString(iTUVQminus1)
             call AddToString(') ')
          ENDIF
       ENDIF
!       call AddToString(' ! NNNNNNNNNNNNNNNN ')
       call writeString(LUPRI);nLines = nLines + 1
    ELSE
       !do not include this term       
    ENDIF
  END subroutine WRITERECURRENCE1

  subroutine WRITERECURRENCE5(CARTDIR,Tq,Uq,Vq,iTUVstart,nTUVP_tmp,&
       & nTUVP,TUVINDEX,JTMQ,JMAX,LUPRI,iTUVPminus1,iTUVPLEnTUVP,&
       & ituvpminus1LEnTUVP,IfacX,TUVindexX,LOOPUNROLL,&
       & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPP,nTUVTMPP2)
    implicit none
    !A to C
    !Theta(i,0,k,0) = i/(2q)*Theta(i-1,0,k-1,0) - p/q*Theta(i+1,0,k-1,0) 
    integer,intent(in) :: CARTDIR,Tq,Uq,Vq,nTUVP,JTMQ,JMAX,iTUVPminus1,nTUVP_tmp
    integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),LUPRI,iTUVstart
    integer,intent(in) :: IfacX(:,:),TUVindexX(:,:),nTUVTMPP,nTUVTMPP2
    logical,intent(in) :: iTUVPLEnTUVP,ituvpminus1LEnTUVP,LOOPUNROLL
    logical,intent(in) :: Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC
    !
    integer :: iTUVQ,iTUVP,iTUVPplus1,I,iTUVQminus1,ituvpminus1x
    character(len=132) :: STRING 
    integer :: iString
    iTUVQ = TUVINDEX(Tq,Uq,Vq)
    IF(CARTDIR.EQ.1)THEN
       !X direction
       iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq)
    ELSEIF(CARTDIR.EQ.2)THEN
       !Y direction
       iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq)                        
    ELSEIF(CARTDIR.EQ.3)THEN
       !Z direction
       iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1)
    ENDIF
    !Result(iTUVP,iTUVQ) = I * Term2(iTUVPminus1,iTUVQminus1) + Term4(iTUVPplus1,iTUVQminus1)
    DO ituvpminus1x = iTUVstart,nTUVP_tmp
       IF(LOOPUNROLL.OR.ituvpminus1x.EQ.iTUVstart)THEN
          IF(LOOPUNROLL)THEN
             iTUVP = TUVindexX(ituvpminus1x,CARTDIR)
             call initString(5)          
          ELSE
             IF(DoOpenACC)WRITE(LUPRI,'(A)')'!$ACC LOOP SEQ'
             call initString(5)          
             call AddToString('do ituvpminus1 = ')
             call AddToString(iTUVstart)
             call AddToString(',')
             call AddToString(nTUVP_tmp)
             call writeString(LUPRI);nLines = nLines + 1
             
             call initString(6)          
             call AddToString('iTUVP = TUVindex')
             call AddToString('X')
             call AddToString(CARTDIR)
             call AddToString('_')
             call AddToString(nTUVTMPP)
             call AddToString('(ituvpminus1)')
             call writeString(LUPRI);nLines = nLines + 1
             call initString(6)          
          ENDIF

          IF(iTUVPLEnTUVP)THEN
             call AddToString('Tmp0(')
          ELSE
             call AddToString('Tmp')
             call AddToString(JTMQ)
             call AddToString('(')
          ENDIF
          IF(LOOPUNROLL)THEN
             call AddToString(iTUVP)
          ELSE
             call AddToString('iTUVP')
          ENDIF
          call AddToString(',')
          call AddToString(iTUVQ)
          call AddToString(') = ')

          IF(iTUVPLEnTUVP)THEN
             call AddToString('Tmp0(')
          ELSE
             call AddToString('Tmp')
             call AddToString(JTMQ)
             call AddToString('(')
          ENDIF
          IF(LOOPUNROLL)THEN
             call AddToString(iTUVP)
          ELSE
             call AddToString('iTUVP')
          ENDIF
          call AddToString(',')
          call AddToString(iTUVQ)
          call AddToString(') ')

          IF(LOOPUNROLL)THEN
             IF(IfacX(ituvpminus1x,CARTDIR).EQ.1)THEN
                call AddToString('+ inv2expQ*')
             ELSE
                call AddToString('+ ')
                call AddToString(IfacX(ituvpminus1x,CARTDIR))
                call AddToString('*inv2expQ*')
             ENDIF
          ELSE
                call AddToString('+ IfacX')
                call AddToString(CARTDIR)
                call AddToString('_')
                call AddToString(nTUVTMPP2)
                call AddToString('(ituvpminus1)*inv2expQ*')             
          ENDIF
          IF(iTUVQminus1.EQ.1)THEN
             call AddToString('Aux(')
             IF(nPrimLast)THEN
                IF(LOOPUNROLL)THEN
                   call AddToString(ituvpminus1x)
                ELSE
                   call AddToString('ituvpminus1')
                ENDIF
                call AddToString(',')
             ENDIF
             IF(COLLAPSE)THEN
                IF(Gen.OR.Seg1Prim)THEN
                   call AddToString('IP')
                ELSE
                   !iPrimQ,iPrimP,iPassP
                   call AddToString('iPrimQ,iPrimP,iPassP')
                ENDIF
             ELSE             
                call AddToString('IP')
             ENDIF             
             IF(nPrimLast)THEN
                call AddToString(') ')
             ELSE
                call AddToString(',')
                IF(LOOPUNROLL)THEN
                   call AddToString(ituvpminus1x)
                ELSE
                   call AddToString('ituvpminus1')
                ENDIF
                call AddToString(') ')                
             ENDIF
          ELSE
             IF(ituvpminus1LEnTUVP)THEN
                call AddToString('Tmp0(')
             ELSE
                call AddToString('Tmp')
                call AddToString(JTMQ-1)
                call AddToString('(')
             ENDIF
             IF(LOOPUNROLL)THEN
                call AddToString(ituvpminus1x)
             ELSE
                call AddToString('ituvpminus1')
             ENDIF
             call AddToString(',')
             call AddToString(iTUVQminus1)
             call AddToString(') ')
          ENDIF
          call writeString(LUPRI);nLines = nLines + 1
       ENDIF
    ENDDO
    IF(.NOT.LOOPUNROLL)THEN    
       call initString(5)          
       call AddToString('enddo')
       call writeString(LUPRI);nLines = nLines + 1
    ENDIF

  END subroutine WRITERECURRENCE5

subroutine WRITERECURRENCE2(CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,TUVINDEX,JTMQ,JMAX,LUPRI,&
       & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE)
  implicit none
  !A to C
  !Theta(i,0,k,0) = i/(2q)*Theta(i-1,0,k-1,0) - p/q*Theta(i+1,0,k-1,0) 
  integer,intent(in) :: CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,JTMQ,JMAX
  integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),LUPRI
  logical,intent(in) :: Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE
  !
  integer :: iTUVQ,iTUVP,iTUVPplus1,iTUVPminus1,I,iTUVQminus1
  character(len=132) :: STRING 
  integer :: iString
  iTUVP = TUVINDEX(Tp,Up,Vp)
  iTUVQ = TUVINDEX(Tq,Uq,Vq)
  IF(CARTDIR.EQ.1)THEN
     !X direction
     iTUVPplus1 = TUVINDEX(Tp+1,Up,Vp)
     iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp)
     I = Tp
     iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq)
  ELSEIF(CARTDIR.EQ.2)THEN
     !Y direction
     iTUVPplus1 = TUVINDEX(Tp,Up+1,Vp)
     iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp)
     I = Up
     iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq)                        
  ELSEIF(CARTDIR.EQ.3)THEN
     !Z direction
     iTUVPplus1 = TUVINDEX(Tp,Up,Vp+1)
     iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1)
     I = Vp
     iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1)
  ENDIF
  !Result(iTUVP,iTUVQ) = I * Term2(iTUVPminus1,iTUVQminus1) + Term4(iTUVPplus1,iTUVQminus1)

  !step 1 add blanks
  call initString(5)
  !step 2 determine where to put the 
  IF(iTUVP.LE.nTUVP)THEN
     call AddToString('Tmp0(')
!     call AddToString('Aux2(')
  ELSE
     call AddToString('Tmp')
     call AddToString(JTMQ)
     call AddToString('(')
  ENDIF
  call AddToString(iTUVP)
  call AddToString(',')
  call AddToString(iTUVQ)
!  IF(iTUVP.LE.nTUVP)THEN
!     call AddToString(',IP) = ')
!  ELSE
     call AddToString(') = ')
!  ENDIF
!====================

  IF(iTUVP.LE.nTUVP)THEN
     call AddToString('Tmp0(')
  ELSE
     call AddToString('Tmp')
     call AddToString(JTMQ)
     call AddToString('(')
  ENDIF
  call AddToString(iTUVP)
  call AddToString(',')
  call AddToString(iTUVQ)
!  IF(iTUVP.LE.nTUVP)THEN
!     call AddToString(',IP) ')
!  ELSE
     call AddToString(') ')
!  ENDIF

!=====================
  !step 6 determine if the third term: 
  !  pinvq*Aux(',ituvpplus1x,',IP)'
  !  pinvq*Tmp',JTMQ-1,'(',ituvpplus1x,',',iTUVQminus1x,')'
  !should be included and if it should use Aux or Tmp 
  call AddToString('+ pinvq*')
  IF(iTUVQminus1.EQ.1)THEN
     call AddToString('Aux(')
     IF(nPrimLast)THEN
        call AddToString(ituvpplus1)
        call AddToString(',')
     ENDIF
     IF(COLLAPSE)THEN
        IF(Gen.OR.Seg1Prim)THEN
           call AddToString('IP')
        ELSE
           !iPrimQ,iPrimP,iPassP
           call AddToString('iPrimQ,iPrimP,iPassP')
        ENDIF
     ELSE             
        call AddToString('IP')
     ENDIF
     IF(nPrimLast)THEN
        call AddToString(') ')
     ELSE
        call AddToString(',')
        call AddToString(ituvpplus1)
        call AddToString(') ')
     ENDIF
  ELSE
     IF(ituvpplus1.LE.nTUVP)THEN
        call AddToString('Tmp0(')
     ELSE
        call AddToString('Tmp')
        call AddToString(JTMQ-1)
        call AddToString('(')
     ENDIF
     call AddToString(ituvpplus1)
     call AddToString(',')
     call AddToString(iTUVQminus1)
     call AddToString(') ')
  ENDIF
  !Final step write the string
  call writeString(LUPRI);nLines = nLines + 1
!  WRITE(LUFILE,'(A)') STRING(1:iSTRING-1)
END subroutine WRITERECURRENCE2

subroutine WRITERECURRENCE4(CARTDIR,Tq,Uq,Vq,iTUVstart,ituvpplus1LEnTUVP,nTUVP_tmp,&
     & nTUVP,TUVINDEX,JTMQ,JMAX,LUPRI,TUVindexX,LOOPUNROLL,&
     & Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE,DoOpenACC,nTUVTMPP,nTUVTMPP2)
  implicit none
  !A to C
  !Theta(i,0,k,0) = i/(2q)*Theta(i-1,0,k-1,0) - p/q*Theta(i+1,0,k-1,0) 
  integer,intent(in) :: CARTDIR,Tq,Uq,Vq,nTUVP,JTMQ,JMAX,iTUVstart,nTUVP_tmp
  integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),LUPRI,TUVindexX(:,:)
  integer,intent(in) :: nTUVTMPP,nTUVTMPP2
  logical,intent(in) :: ituvpplus1LEnTUVP,LOOPUNROLL,DoOpenACC
  logical,intent(in) :: Gen,SegQ,SegP,Seg,Seg1Prim,COLLAPSE
  !
  integer :: iTUVQ,iTUVP,iTUVPplus1,iTUVPminus1,I,iTUVQminus1
  character(len=132) :: STRING 
  integer :: iString
  iTUVQ = TUVINDEX(Tq,Uq,Vq)
  IF(CARTDIR.EQ.1)THEN
     !X direction
     iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq)
  ELSEIF(CARTDIR.EQ.2)THEN
     !Y direction
     iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq)                        
  ELSEIF(CARTDIR.EQ.3)THEN
     !Z direction
     iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1)
  ENDIF
  !Result(iTUVP,iTUVQ) = I * Term2(iTUVPminus1,iTUVQminus1) + Term4(iTUVPplus1,iTUVQminus1)
  DO ituvP = iTUVstart,nTUVP_tmp
     IF(LOOPUNROLL.OR.ituvP .EQ. iTUVstart)THEN        
        IF(LOOPUNROLL)THEN
           iTUVPplus1 = TUVindexX(ituvP,CARTDIR)
           call initString(5)          
        ELSE
           call initString(5)                 
           IF(DoOpenACC)WRITE(LUPRI,'(A)')'!$ACC LOOP SEQ'     
           call AddToString('do iTUVP = ')
           call AddToString(iTUVstart)
           call AddToString(',')
           call AddToString(nTUVP_tmp)
           call writeString(LUPRI);nLines = nLines + 1
        
           call initString(6)          
           call AddToString('iTUVplus1 = TUVindexX')
           call AddToString(CARTDIR)
           call AddToString('_')
           call AddToString(nTUVTMPP)
           call AddToString('(iTUVP)')
           call writeString(LUPRI);nLines = nLines + 1
           call initString(6)          
        ENDIF

        IF(iTUVstart.LE.nTUVP)THEN
           call AddToString('Tmp0(')
        ELSE
           call AddToString('Tmp')
           call AddToString(JTMQ)
           call AddToString('(')
        ENDIF
        IF(LOOPUNROLL)THEN
           call AddToString(iTUVP)
        ELSE
           call AddToString('iTUVP')
        ENDIF
        call AddToString(',')
        call AddToString(iTUVQ)
        call AddToString(')')
        call AddToString(' = ')
        IF(iTUVstart.LE.nTUVP)THEN
           call AddToString('Tmp0(')
        ELSE
           call AddToString('Tmp')
           call AddToString(JTMQ)
           call AddToString('(')
        ENDIF
        IF(LOOPUNROLL)THEN
           call AddToString(iTUVP)
        ELSE
           call AddToString('iTUVP')
        ENDIF
        call AddToString(',')
        call AddToString(iTUVQ)
        call AddToString(')')
        call AddToString(' + ')
        call AddToString('pinvq*')
        IF(iTUVQminus1.EQ.1)THEN
           call AddToString('Aux(')
           IF(nPrimLast)THEN
              IF(LOOPUNROLL)THEN
                 call AddToString(iTUVPplus1)
              ELSE
                 call AddToString('iTUVplus1')
              ENDIF
              call AddToString(',')
           ENDIF
           IF(COLLAPSE)THEN
              IF(Gen.OR.Seg1Prim)THEN
                 call AddToString('IP')
              ELSE
                 !iPrimQ,iPrimP,iPassP
                 call AddToString('iPrimQ,iPrimP,iPassP')
              ENDIF
           ELSE             
              call AddToString('IP')
           ENDIF
           IF(nPrimLast)THEN
              call AddToString(')')
           ELSE
              call AddToString(',')
              IF(LOOPUNROLL)THEN
                 call AddToString(iTUVPplus1)
              ELSE
                 call AddToString('iTUVplus1')
              ENDIF
              call AddToString(')')
           ENDIF
        ELSE
           IF(ituvpplus1LEnTUVP)THEN
              call AddToString('Tmp0')
           ELSE
              call AddToString('Tmp')
              call AddToString(JTMQ-1)
           ENDIF
           call AddToString('(')
           IF(LOOPUNROLL)THEN
              call AddToString(iTUVPplus1)
           ELSE
              call AddToString('iTUVplus1')
           ENDIF
           call AddToString(',')
           call AddToString(iTUVQminus1)
           call AddToString(')')
        ENDIF
        call writeString(LUPRI);nLines = nLines + 1
     ENDIF
  ENDDO
!  WRITE(LUFILE,'(A,I3,A,I3,A)')'      Tmp0(iTUVP,',iTUVQ,') = Tmp0(iTUVP,',iTUVQ,') + pinvq*Aux(iTUVplus1,IP) '
  IF(.NOT.LOOPUNROLL)THEN
     WRITE(LUPRI,'(A)')'     enddo'
  ENDIF
END subroutine WRITERECURRENCE4

end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
