MODULE TESTMODULE
  use stringsMODULE
  logical,save :: nPrimLast

CONTAINS
subroutine PASSsub
  IMPLICIT NONE
  INTEGER :: JMAX,nTUV,nTUVprev,ituvP,J,Tp,Up,Vp,N,N2,ntuvP,ituv,C
  INTEGER :: nTUVQ,nTUVTMPP,nTUVTMPQ,JPQ,JP,JQ,nTUVTMP,nTUVTMPprev
  integer :: tq,uq,vq,ituvq,ituvpminus1x,ituvpminus1y,ituvpminus1z
  integer :: iTUVplus1x,iTUVplus1y,iTUVplus1z,nTUVTMPPs,I,nTUVTMPQ2,nTUVP2
  Integer :: MaxAngmomQP,nTUVplus,JTMP,ntuvprev2,ntuvprev3,iTUVQminus1
  !    logical :: CREATED(-2:8,-2:8,-2:8)
  logical,pointer :: CREATED(:,:,:),ituvqplus1LEnTUVQarray(:)
  logical,pointer :: ituvpminus1LEnTUVParray(:)
  logical :: TREC,UREC,VREC,TREC2,UREC2,VREC2,DoneCartDir(3),ituvqplus1LEnTUVQ
  logical :: ituvqminus1LEnTUVQ
  integer,pointer :: TUVINDEX(:,:,:)
  integer,pointer :: TINDEX(:)
  integer,pointer :: UINDEX(:)
  integer,pointer :: VINDEX(:)
  integer,pointer :: JINDEX(:)
  integer :: nTUVLIST,nTUVLISTactual,CARTDIR,iTUVPminus1,iTUVPminus2,Tpminus1
  integer,pointer :: TwoTermTUVLIST(:)
  integer :: iTUVQplus1,nLength
  integer :: LUFILE,LUSPECIAL,LUFILE1,GPUrun
  Character(len=51) :: FileName    
  Character(len=44) :: FromPrimLabel, ToPrimLabel
  Character(len=1)  :: FromLabel,ToLabel,FromExpLabel,ToExpLabel,SIGN2
  Character(len=8)  :: SegLabel
  character(len=3) :: ARCSTRING
  character(len=20) :: PrimLabelAux
  integer :: iPrimLabelAux
  integer :: iseg,ifile,iseglabel
  logical :: Gen,SegQ,Segp,Seg,Seg1Prim,LOOPUNROLL,DoOpenMP,DoOpenACC
  logical :: Collapse,CPU
  integer,pointer :: IfacX(:,:),TUVindexX(:,:)
  !    LUSPECIAL = 2
  !    open(unit = LUSPECIAL, file="runNewTRQPBasic.F90",status="unknown")
  !    WRITE(LUSPECIAL,'(A)')'MODULE AGC_OBS_TRANSFERRECURRENCEMODBASIC'
  !    WRITE(LUSPECIAL,'(A)')' use IchorPrecisionModule'
  !    WRITE(LUSPECIAL,'(A)')'  '
  !    WRITE(LUSPECIAL,'(A)')' CONTAINS'
  DO GPUrun = 1,2
    CPU = .TRUE.
    IF(GPUrun.EQ.2)CPU = .FALSE.
    nPrimLAST = .FALSE.
    IF(CPU)nPrimLAST = .TRUE.
    DoOpenMP = .FALSE.
    DoOpenACC = .FALSE.
    COLLAPSE=.TRUE.
    IF(CPU)DoOpenMP = .TRUE.
    IF(.NOT.CPU)DoOpenACC = .TRUE.
    IF(CPU)THEN
       ARCSTRING = 'CPU'
    ELSE
       ARCSTRING = 'GPU'
    ENDIF
    DO iseg = 1,5
       DO ifile = 1,4
          IF(ifile.EQ.1)THEN
             FromLabel = 'C'; ToLabel = 'A'; FromExpLabel = 'D'; ToExpLabel = 'B'
             FromPrimLabel = 'iPrimD = (iPrimQ-1)/nPrimC+1                '; SIGN2= '+'  
             ToPrimLabel   = 'iPrimB = (iPrimP-1)/nPrimA+1                '
          ELSEIF(ifile.EQ.2)THEN
             FromLabel = 'C'; ToLabel = 'B'; FromExpLabel = 'D'; ToExpLabel = 'A' 
             FromPrimLabel = 'iPrimD = (iPrimQ-1)/nPrimC+1                '; SIGN2= '+'  
             ToPrimLabel   = 'iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA'
          ELSEIF(ifile.EQ.3)THEN
             FromLabel = 'D'; ToLabel = 'A'; FromExpLabel = 'C'; ToExpLabel = 'B' 
             FromPrimLabel = 'iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC'; SIGN2= '-'  
             ToPrimLabel   = 'iPrimB = (iPrimP-1)/nPrimA+1                '
          ELSEIF(ifile.EQ.4)THEN
             FromLabel = 'D'; ToLabel = 'B'; FromExpLabel = 'C'; ToExpLabel = 'A' 
             FromPrimLabel = 'iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC'; SIGN2= '-'  
             ToPrimLabel   = 'iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA'
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
          DO I = 1,48
             FileName(I:I) = ' '
          ENDDO
          WRITE(FileName,'(8A)')'runNewTransferRecurrenceQP',FromLabel,'to',ToLabel,SegLabel(1:iSegLabel),'output',ARCSTRING,'.F90'
          print*,'FileName:',FileName
          LUFILE1 = 1
          open(unit = LUFILE1, file=TRIM(FileName),status="unknown")

          WRITE(LUFILE1,'(7A)')'MODULE AGC_',ARCSTRING,'_OBS_TRMOD',FromLabel,'to',ToLabel,SegLabel(1:iSegLabel)
          WRITE(LUFILE1,'(A)')' use IchorPrecisionModule'
          WRITE(LUFILE1,'(A)')'  '
          WRITE(LUFILE1,'(A)')' CONTAINS'
          MaxAngmomQP = 8-1 !PDDD highest possible 

          DO JMAX=2,MaxAngmomQP
             IF((ifile.EQ.3.OR.ifile.EQ.4).AND.JMAX.GT.MaxAngmomQP-1)CYCLE !PDPD higest possible
             DO JP = 1, JMAX
                JQ = JMAX-JP
                IF(JP.GE.JQ)CYCLE
                IF(JP.EQ.0)CYCLE
                IF(JQ.GT.4)CYCLE
                IF((ifile.EQ.3.OR.ifile.EQ.4).AND.JQ.GT.4-1)CYCLE !PDPD higest possible
                IF(JP.GT.4)CYCLE
                IF(JMAX.LT.5)THEN
!                   LUFILE = LUSPECIAL
                   LUFILE = LUFILE1
                   LOOPUNROLL = .TRUE.
                ELSE
                   LUFILE = LUFILE1
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

                call initString(9)          
                call AddToString('& Pexp,Qexp,Pdistance12,Qdistance12,')
                call AddToString(FromExpLabel)
                call AddToString('exp,')
                call AddToString(ToExpLabel)
                call AddToString('exp,nPrimA,nPrimB,nPrimC,nPrimD,&')
                call writeString(LUFILE)
                call initString(9)          
                call AddToString('& MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)')
                call writeString(LUFILE)

                WRITE(LUFILE,'(A)')'  implicit none'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses'
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)'
                ELSE
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)'
                ENDIF
                WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
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
                !          WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Bexp(nPrimB),Dexp(nPrimD)'
                IF(.NOT.Seg1Prim)THEN
                   IF(COLLAPSE)THEN
                      IF(Gen)THEN
                         IF(nPrimLast)THEN
                            WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPrimQ*nPrimP*nPasses)'
                         ELSE
                            WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(nPrimQ*nPrimP*nPasses,',nTUV,')'
                         ENDIF
                         PrimLabelAux = 'iP'
                         iPrimLabelAux = 2
                      ELSE
                         IF(nPrimLast)THEN
                            WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPrimQ,nPrimP,nPasses)'
                         ELSE
                            WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,',nTUV,')'
                         ENDIF
                         PrimLabelAux = 'iPrimQ,iPrimP,iPassP'
                         iPrimLabelAux = 20
                      ENDIF
                   ELSE
                      IF(nPrimLast)THEN
                         WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPrimQ*nPrimP*nPasses)'
                      ELSE
                         WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(nPrimQ*nPrimP*nPasses,',nTUV,')'
                      ENDIF
                      PrimLabelAux = 'iP'
                      iPrimLabelAux = 2
                   ENDIF
                ELSE
                   IF(nPrimLast)THEN
                      WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPasses)'
                   ELSE
                      WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: Aux(nPasses,',nTUV,')'
                   ENDIF
                   PrimLabelAux = 'iP'
                   iPrimLabelAux = 2
                ENDIF
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
                IF(nPrimLast)THEN
                   WRITE(LUFILE,'(A)')'!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)'
                ELSE
                   WRITE(LUFILE,'(A)')'!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)'
                ENDIF
                WRITE(LUFILE,'(A)')'  !Local variables'
                !             WRITE(LUFILE,'(A)')'  real(realk) :: Pexpfac,PREF'
                !             WRITE(LUFILE,'(A,i4,A)')'  real(realk) :: TwoTerms(',MAX(1,nTUVprev2-nTUVprev3),')'
                WRITE(LUFILE,'(A,I3,A,I3,A)')'  real(realk) :: Tmp0(',nTUVQ,',',nTUVP,')'
                WRITE(LUFILE,'(A)')'! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2'
                DO JTMP=1,JP-1
                   nTUVTMPPs=(JTMP)*(JTMP+1)*(JTMP+2)/6
                   nTUVTMPP=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
                   nTUVTMPQ=(JPQ-JTMP+1)*(JPQ-JTMP+2)*(JPQ-JTMP+3)/6
                   if(JTMP.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVQ+1,':',nTUVTMPQ,',',nTUVTMPPs+1,':',nTUVTMPP,')'
                   else
                      WRITE(LUFILE,'(A,I2,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVQ+1,':',nTUVTMPQ,',',nTUVTMPPs+1,':',nTUVTMPP,')'
                   endif
                ENDDO
                WRITE(LUFILE,'(A)')'!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines '
                WRITE(LUFILE,'(A)')'  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB'
                IF(.NOT.Seg1Prim)THEN
                   call initString(2)          
                   call AddToString('integer :: iPrim')
                   call AddToString(FromExpLabel)
                   call AddToString(',iPrim')
                   call AddToString(ToExpLabel)
                   call writeString(LUFILE)
                   !          WRITE(LUFILE,'(A)')'  integer :: iPrimB,iPrimD'
                ENDIF
                WRITE(LUFILE,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk'
                WRITE(LUFILE,'(A)')'  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP'
                call initString(2)          
                call AddToString('real(realk) :: exp')
                call AddToString(ToExpLabel)
                call AddToString('X,exp')
                call AddToString(ToExpLabel)
                call AddToString('Y,exp')
                call AddToString(ToExpLabel)
                call AddToString('Z')
                call writeString(LUFILE)                
!                WRITE(LUFILE,'(A)')'  real(realk) :: expBX,expBY,expBZ'
                WRITE(LUFILE,'(A)')'  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp'
                !==========================================================================================================
                !         Build the TUVindexX
                !==========================================================================================================
                allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
                CREATED  = .FALSE.
                CREATED(0,0,0) = .TRUE.
                nTUVTMPQ=(JPQ-1+1)*(JPQ-1+2)*(JPQ-1+3)/6                         
                allocate(TUVindexX(nTUVTMPQ,3))
                DO JTMP=1,1!JP
                   DoneCartDir = .FALSE.
                   DO Tp=JTMP,0,-1       
                      DO Up=JTMP-Tp,0,-1
                         Vp=JTMP-Tp-Up  
                         CALL DETERMINE_CARTDIR(CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1,Tp,Up,Vp,CREATED,JMAX,TUVINDEX)
                         nTUVTMPQ=(JPQ-JTMP+1)*(JPQ-JTMP+2)*(JPQ-JTMP+3)/6                         
                         IF(.NOT.DoneCartDir(CARTDIR))THEN
                            IF(.NOT.LOOPUNROLL)THEN
                               call initString(2)          
                               call AddToString('!CARTDIR = ')
                               call AddToString(CARTDIR)
                               call writeString(LUFILE)
                               call initString(2)          
                               call AddToString('integer,parameter, dimension(')
                               call AddToString(nTUVTMPQ)
                               call AddToString(') :: TUVindex')
                               !                      call AddToString(JTMQ)
                               call AddToString('X')
                               call AddToString(CARTDIR)
                               call AddToString(' = (/ ')
                            ENDIF
                            nLength = 10
                            do iTUVQ = 1,nTUVTMPQ!nTUVQ
                               nLength = nLength + 1
                               Tq = Tindex(iTUVq) 
                               Uq = Uindex(iTUVq) 
                               Vq = Vindex(iTUVq)
                               IF(CARTDIR.EQ.1)THEN
                                  iTUVQplus1 = TUVINDEX(Tq+1,Uq,Vq)
                               ELSEIF(CARTDIR.EQ.2)THEN
                                  iTUVQplus1 = TUVINDEX(Tq,Uq+1,Vq)
                               ELSEIF(CARTDIR.EQ.3)THEN
                                  iTUVQplus1 = TUVINDEX(Tq,Uq,Vq+1)
                               ENDIF
                               IF(.NOT.LOOPUNROLL)THEN
                                  call AddToString(ituvqplus1)
                               ENDIF
                               TUVindexX(iTUVQ,CARTDIR) = ituvqplus1
                               IF(.NOT.LOOPUNROLL)THEN
                                  IF(iTUVQ.NE.nTUVTMPQ) call AddToString(',')
                                  IF(nLength.EQ.17)THEN
                                     nLength = 0
                                     call AddToString('&')
                                     call writeString(LUFILE)
                                     call initString(5)          
                                     call AddToString('     & ')
                                  ENDIF
                               ENDIF
                            enddo
                            IF(.NOT.LOOPUNROLL)THEN
                               call AddToString(' /)')
                               call writeString(LUFILE)                   
                            ENDIF
                            DoneCartDir(CARTDIR) = .TRUE.
                         ENDIF
                         CREATED(Tp,Up,Vp) = .TRUE.
                      ENDDO
                   ENDDO
                ENDDO
                !==========================================================================================================
                !         Build the IfacX integer          
                !==========================================================================================================
                CREATED  = .FALSE.
                CREATED(0,0,0) = .TRUE.
                nTUVTMPQ2=(JPQ-1)*(JPQ-1+1)*(JPQ-1+2)/6                         
                allocate(IfacX(nTUVTMPQ2,3))
                DO JTMP=1,1!JP
                   DoneCartDir = .FALSE.
                   DO Tp=JTMP,0,-1       
                      DO Up=JTMP-Tp,0,-1
                         Vp=JTMP-Tp-Up  
                         CALL DETERMINE_CARTDIR(CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1,Tp,Up,Vp,CREATED,JMAX,TUVINDEX)
                         nTUVTMPQ2=(JPQ-JTMP)*(JPQ-JTMP+1)*(JPQ-JTMP+2)/6
                         IF(.NOT.DoneCartDir(CARTDIR))THEN
                            IF(.NOT.LOOPUNROLL)THEN
                               call initString(2)          
                               call AddToString('!CARTDIR = ')
                               call AddToString(CARTDIR)
                               call writeString(LUFILE)
                               call initString(2)          
                               call AddToString('integer,parameter, dimension(')
                               call AddToString(nTUVTMPQ2)
                               call AddToString(') :: Ifac')
                               call AddToString('X')
                               call AddToString(CARTDIR)
                               call AddToString(' = (/ ')
                            ENDIF
                            nLength = 10
                            nTUVTMPQ2=(JPQ-JTMP)*(JPQ-JTMP+1)*(JPQ-JTMP+2)/6
                            do iTUVQminus1 = 1,nTUVTMPQ2
                               nLength = nLength + 1
                               Tq = Tindex(iTUVqminus1) 
                               Uq = Uindex(iTUVqminus1) 
                               Vq = Vindex(iTUVqminus1)
                               IF(CARTDIR.EQ.1)THEN
                                  I = Tq+1
                               ELSEIF(CARTDIR.EQ.2)THEN
                                  I = Uq+1
                               ELSEIF(CARTDIR.EQ.3)THEN
                                  I = Vq+1
                               ENDIF
                               IF(.NOT.LOOPUNROLL)THEN
                                  call AddToString(I)
                               ENDIF
                               IfacX(iTUVQminus1,CARTDIR) = I                                
                               IF(.NOT.LOOPUNROLL)THEN
                                  IF(iTUVQminus1.NE.nTUVTMPQ2) call AddToString(',')
                                  IF(nLength.EQ.17)THEN
                                     nLength = 0
                                     call AddToString('&')
                                     call writeString(LUFILE)
                                     call initString(5)          
                                     call AddToString('     & ')
                                  ENDIF
                               ENDIF
                            enddo
                            IF(.NOT.LOOPUNROLL)THEN
                               call AddToString(' /)')
                               call writeString(LUFILE)
                            ENDIF
                            DoneCartDir(CARTDIR) = .TRUE.
                         ENDIF
                         CREATED(Tp,Up,Vp) = .TRUE.
                      ENDDO
                   ENDDO
                ENDDO
                deallocate(CREATED)

                !init AUX
                IF(COLLAPSE)THEN
                   IF(SegP.OR.SegQ.OR.Seg)THEN
!                      IF(Doopenmp)WRITE(LUFILE,'(A)')'!$OMP PARALLEL DO DEFAULT(shared) COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ) '
                      IF(Doopenmp)WRITE(LUFILE,'(A)')'!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)'
                      IF(SegP)THEN
                         IF(DoopenACC)WRITE(LUFILE,'(A)')'!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,Aux2)'
                         WRITE(LUFILE,'(A)')'  DO iP = 1,nPrimQ*nPasses'
                      ENDIF
                      IF(SegQ)THEN
                         IF(DoopenACC)WRITE(LUFILE,'(A)')'!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimP,nPasses,Aux2)'
                         WRITE(LUFILE,'(A)')'  DO iP = 1,nPrimP*nPasses'
                      ENDIF
                      IF(Seg)THEN
                         IF(DoopenACC)WRITE(LUFILE,'(A)')'!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPasses,Aux2)'
                         WRITE(LUFILE,'(A)') '  DO iP = 1,nPasses'
                      ENDIF
                      WRITE(LUFILE,'(A,I3)')'   DO iTUVQ=1,',nTUVQ
                      WRITE(LUFILE,'(A,I3)')'    DO iTUVP=1,',nTUVP
                      IF(nPrimLast)THEN
                         WRITE(LUFILE,'(A)')   '     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk'
                      ELSE
                         WRITE(LUFILE,'(A)')   '     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk'
                      ENDIF
                      WRITE(LUFILE,'(A)')   '    ENDDO'
                      WRITE(LUFILE,'(A)')   '   ENDDO'
                      WRITE(LUFILE,'(A)')   '  ENDDO'
                      IF(Doopenmp)WRITE(LUFILE,'(A)')'!$OMP END DO'
!                      IF(Doopenmp)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'
                   ELSE
                      !no need to init
                   ENDIF
                ENDIF

                !OPENMP
                IF(Doopenmp)THEN
                   WRITE(LUFILE,'(A)')'!$OMP DO &'
!                   WRITE(LUFILE,'(A)')'!$OMP PARALLEL DO DEFAULT(none) &'
                   WRITE(LUFILE,'(A)')'!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&'
                   IF(Seg)THEN
                      WRITE(LUFILE,'(A)')'!$OMP         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&'
                   ELSE
                      WRITE(LUFILE,'(A)')'!$OMP         iP,iPrimQ,iPrimP,iPassP,&'
                   ENDIF
                   WRITE(LUFILE,'(A,A,A,A,A,A,A)')'!$OMP         exp',ToExpLabel,'X,exp',ToExpLabel,'Y,exp',ToExpLabel,'Z,&'
                   IF(.NOT.seg1Prim)THEN
                      WRITE(LUFILE,'(A,A,A,A,A,A,A)')'!$OMP         iPrim',FromExpLabel,',iPrim',ToExpLabel,',&'
                   ENDIF
                   WRITE(LUFILE,'(A)')'!$OMP         Tmp0,&'
                   DO JTMP=1,JP-1
                      if(JTMP.LT.10)THEN
                         WRITE(LUFILE,'(A,I1,A)')'!$OMP         Tmp',JTMP,',&'
                      else
                         WRITE(LUFILE,'(A,I2,A)')'!$OMP         Tmp',JTMP,',&'
                      endif
                   ENDDO
                   WRITE(LUFILE,'(A)')'!$OMP         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) '
!                   WRITE(LUFILE,'(A)')'!$OMP         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &'
!                   WRITE(LUFILE,'(A)')'!$OMP SHARED(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&'
!                   WRITE(LUFILE,'(A,A,A,A,A)')'!$OMP        ',ToExpLabel,'exp,',FromExpLabel,'exp,&'
!                   WRITE(LUFILE,'(A)')'!$OMP        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux)'
                ENDIF
                !OpenACC
                IF(DoopenACC)THEN
                   WRITE(LUFILE,'(A)')'!$ACC PARALLEL LOOP &'
                   WRITE(LUFILE,'(A)')'!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&'
                   IF(Seg)THEN
                      WRITE(LUFILE,'(A)')'!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&'
                   ELSE
                      WRITE(LUFILE,'(A)')'!$ACC         iP,iPrimQ,iPrimP,iPassP,&'
                   ENDIF
                   WRITE(LUFILE,'(A,A,A,A,A,A,A)')'!$ACC         exp',ToExpLabel,'X,exp',ToExpLabel,'Y,exp',ToExpLabel,'Z,&'
                   IF(.NOT.seg1Prim)THEN
                      WRITE(LUFILE,'(A,A,A,A,A,A,A)')'!$ACC         iPrim',FromExpLabel,',iPrim',ToExpLabel,',&'
                   ENDIF
                   WRITE(LUFILE,'(A)')'!$ACC         Tmp0,&'
                   DO JTMP=1,JP-1
                      if(JTMP.LT.10)THEN
                         WRITE(LUFILE,'(A,I1,A)')'!$ACC         Tmp',JTMP,',&'
                      else
                         WRITE(LUFILE,'(A,I2,A)')'!$ACC         Tmp',JTMP,',&'
                      endif
                   ENDDO
                   WRITE(LUFILE,'(A)')'!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &'
                   WRITE(LUFILE,'(A)')'!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&'
                   WRITE(LUFILE,'(A,A,A,A,A)')'!$ACC        ',ToExpLabel,'exp,',FromExpLabel,'exp,&'
                   WRITE(LUFILE,'(A)')'!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux)'
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
                      WRITE(LUFILE,'(A)') '   DO iPrimP=1, nPrimP'
                      WRITE(LUFILE,'(A)') '    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ'
                      WRITE(LUFILE,'(A)') '    iPassP = (iP-1)/nPrimQ + 1'
                   ELSEIF(SegQ)THEN
                      WRITE(LUFILE,'(A)') '  DO iP = 1,nPrimP*nPasses'
                      WRITE(LUFILE,'(A)') '   DO iPrimQ=1, nPrimQ'
                      WRITE(LUFILE,'(A)') '    iPrimP = iP - ((iP-1)/nPrimP)*nPrimP'
                      WRITE(LUFILE,'(A)') '    iPassP = (iP-1)/nPrimP + 1'
                   ELSEIF(Seg)THEN
                      WRITE(LUFILE,'(A)') '  DO iP = 1,nPasses'
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
                   WRITE(LUFILE,'(A)')'   Xcd = Qdistance12(1)'
                   WRITE(LUFILE,'(A)')'   Ycd = Qdistance12(2)'
                   WRITE(LUFILE,'(A)')'   Zcd = Qdistance12(3)'
                   WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                   WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                   WRITE(LUFILE,'(A)')'   Xab = Pdistance12(1,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Yab = Pdistance12(2,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Zab = Pdistance12(3,iAtomA,iAtomB)'
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    expP = Pexp(iPrimP)'
                      WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(iPrimP)'
                   ELSE
                      WRITE(LUFILE,'(A)')'    expP = Pexp(1)'
                      WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(1)'
                   ENDIF
                   WRITE(LUFILE,'(A)')'    inv2expP = D05*invexpP'
                   IF(ToExpLabel.EQ.'A')THEN
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,A)')'    ',ToPrimLabel
                         !          WRITE(LUFILE,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'          
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'X = -',ToExpLabel,'exp(iPrim',ToExpLabel,')*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Y = -',ToExpLabel,'exp(iPrim',ToExpLabel,')*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Z = -',ToExpLabel,'exp(iPrim',ToExpLabel,')*Zab'
                      ELSE
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'X = -',ToExpLabel,'exp(1)*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Y = -',ToExpLabel,'exp(1)*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Z = -',ToExpLabel,'exp(1)*Zab'
                      ENDIF
                   ELSE !B
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,A)')'    ',ToPrimLabel
                         !          WRITE(LUFILE,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'          
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'X = ',ToExpLabel,'exp(iPrim',ToExpLabel,')*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Y = ',ToExpLabel,'exp(iPrim',ToExpLabel,')*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Z = ',ToExpLabel,'exp(iPrim',ToExpLabel,')*Zab'
                      ELSE
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'X = ',ToExpLabel,'exp(1)*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Y = ',ToExpLabel,'exp(1)*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Z = ',ToExpLabel,'exp(1)*Zab'
                      ENDIF
                   ENDIF
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A,A)')'     ',FromPrimLabel
                   ENDIF
                ELSE !no collapse
                   WRITE(LUFILE,'(A)')'  Xcd = Qdistance12(1)'
                   WRITE(LUFILE,'(A)')'  Ycd = Qdistance12(2)'
                   WRITE(LUFILE,'(A)')'  Zcd = Qdistance12(3)'
                   WRITE(LUFILE,'(A)')'  DO iPassP = 1,nPasses'
                   WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                   WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                   WRITE(LUFILE,'(A)')'   Xab = Pdistance12(1,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Yab = Pdistance12(2,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Zab = Pdistance12(3,iAtomA,iAtomB)'
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
                   IF(SegP)THEN
                      WRITE(LUFILE,'(A)')'   DO iPrimQ = 1,nPrimQ'
                      WRITE(LUFILE,'(A,I3)')'    DO iTUVQ=1,',nTUVQ
                      WRITE(LUFILE,'(A,I3)')'     DO iTUVP=1,',nTUVP
                      IF(nPrimLast)THEN
                         WRITE(LUFILE,'(A)')'      Aux2(iTUVP,iTUVQ,iPrimQ,iPassP) = 0.0E0_realk'
                      ELSE
                         WRITE(LUFILE,'(A)')'      Aux2(iPrimQ,iPassP,iTUVP,iTUVQ) = 0.0E0_realk'
                      ENDIF
                      WRITE(LUFILE,'(A)')'     ENDDO'
                      WRITE(LUFILE,'(A)')'    ENDDO'
                      WRITE(LUFILE,'(A)')'   ENDDO'
                   ENDIF
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'   IP = (iPassP-1)*nPrimQ*nPrimP'
                      WRITE(LUFILE,'(A)')'   DO iPrimP=1, nPrimP'
                   ELSE
                      WRITE(LUFILE,'(A)')'   IP = iPassP'
                      WRITE(LUFILE,'(A)')'   iPrimP=1'
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
                   ENDIF
                   WRITE(LUFILE,'(A)')'    expP = Pexp(iPrimP)'
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(iPrimP)'
                   ELSE
                      WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(1)'
                   ENDIF
                   WRITE(LUFILE,'(A)')'    inv2expP = D05*invexpP'
                   IF(ToExpLabel.EQ.'A')THEN
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,A)')'    ',ToPrimLabel
                         !          WRITE(LUFILE,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'          
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'X = -',ToExpLabel,'exp(iPrim',ToExpLabel,')*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Y = -',ToExpLabel,'exp(iPrim',ToExpLabel,')*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Z = -',ToExpLabel,'exp(iPrim',ToExpLabel,')*Zab'
                      ELSE
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'X = -',ToExpLabel,'exp(1)*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Y = -',ToExpLabel,'exp(1)*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Z = -',ToExpLabel,'exp(1)*Zab'
                      ENDIF
                   ELSE !B
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A,A)')'    ',ToPrimLabel
                         !          WRITE(LUFILE,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'          
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'X = ',ToExpLabel,'exp(iPrim',ToExpLabel,')*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Y = ',ToExpLabel,'exp(iPrim',ToExpLabel,')*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Z = ',ToExpLabel,'exp(iPrim',ToExpLabel,')*Zab'
                      ELSE
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'X = ',ToExpLabel,'exp(1)*Xab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Y = ',ToExpLabel,'exp(1)*Yab'
                         WRITE(LUFILE,'(7A)')'    exp',ToExpLabel,'Z = ',ToExpLabel,'exp(1)*Zab'
                      ENDIF
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
                      WRITE(LUFILE,'(A,A)')'     ',FromPrimLabel
                      !          WRITE(LUFILE,'(A)')'     iPrimD = (iPrimQ-1)/nPrimC+1'
                      WRITE(LUFILE,'(A)')'     IP = IP + 1'
                   ENDIF
                ENDIF !COLLAPSE
                   
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(8A)')'     facX = -(exp',ToExpLabel,'X',SIGN2,FromExpLabel,'exp(iPrim',FromExpLabel,')*Xcd)*invexpP'
                   WRITE(LUFILE,'(8A)')'     facY = -(exp',ToExpLabel,'Y',SIGN2,FromExpLabel,'exp(iPrim',FromExpLabel,')*Ycd)*invexpP'
                   WRITE(LUFILE,'(8A)')'     facZ = -(exp',ToExpLabel,'Z',SIGN2,FromExpLabel,'exp(iPrim',FromExpLabel,')*Zcd)*invexpP'
                ELSE
                   WRITE(LUFILE,'(8A)')'     facX = -(exp',ToExpLabel,'X',SIGN2,FromExpLabel,'exp(1)*Xcd)*invexpP'
                   WRITE(LUFILE,'(8A)')'     facY = -(exp',ToExpLabel,'Y',SIGN2,FromExpLabel,'exp(1)*Ycd)*invexpP'
                   WRITE(LUFILE,'(8A)')'     facZ = -(exp',ToExpLabel,'Z',SIGN2,FromExpLabel,'exp(1)*Zcd)*invexpP'
                ENDIF
                !          WRITE(LUFILE,'(A)')'     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpQ'
                !          WRITE(LUFILE,'(A)')'     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpQ'
                !          WRITE(LUFILE,'(A)')'     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpQ'
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'     qinvp = -Qexp(iPrimQ)*invexpP'
                ELSE
                   WRITE(LUFILE,'(A)')'     qinvp = -Qexp(1)*invexpP'
                ENDIF

                CALL SUBROUTINE_MAIN(LUFILE,JMAX,JP,JQ,nTUVP,nTUVQ,JPQ,nTUV,nTUVPLUS,TUVINDEX,&
                     & TINDEX,UINDEX,VINDEX,JINDEX,nTUVprev3,nTUVprev2,nTUVprev,IfacX,TUVindexX,LOOPUNROLL,&
                     & PrimLabelAux,iPrimLabelAux,DoOpenACC)

                WRITE(LUFILE,'(A,I3)')'!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. '
                WRITE(LUFILE,'(A,I3)')'!    Hopefully Tmp0 is small enough that it can be in cache. '
                IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ'
                WRITE(LUFILE,'(A,I3)')'     DO iTUVQ=1,',nTUVQ
                IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ'
                WRITE(LUFILE,'(A,I3)')'      DO iTUVP=1,',nTUVP
                IF(nPrimLast)THEN
                   IF(COLLAPSE)THEN
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,IP) = Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(Seg)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iP) = Tmp0(iTUVQ,iTUVP)'
                      ENDIF
                   ELSE
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,IP) = Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iPrimP,iPassP) = Aux2(iTUVP,iTUVQ,iPrimP,iPassP) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iPrimQ,iPassP) = Aux2(iTUVP,iTUVQ,iPrimQ,iPassP) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(Seg)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iPassP) = Aux2(iTUVP,iTUVQ,iPassP) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iTUVP,iTUVQ,iPassP) = Tmp0(iTUVQ,iTUVP)'
                      ENDIF
                   ENDIF
                ELSE
                   IF(COLLAPSE)THEN
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(Seg)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)'
                      ENDIF
                   ELSE
                      IF(Gen)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(SegQ)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iPrimP,iPassP,iTUVP,iTUVQ) = Aux2(iPrimP,iPassP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(SegP)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iPrimQ,iPassP,iTUVP,iTUVQ) = Aux2(iPrimQ,iPassP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(Seg)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iPassP,iTUVP,iTUVQ) = Aux2(iPassP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)'
                      ELSEIF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')   '        Aux2(iPassP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)'
                      ENDIF
                   ENDIF
                ENDIF
                WRITE(LUFILE,'(A)')   '      ENDDO'
                WRITE(LUFILE,'(A)')   '     ENDDO'
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
                IF(Doopenmp)WRITE(LUFILE,'(A)')'!$OMP END DO'
!                IF(Doopenmp)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'

                call initString(1)          
                call AddToString('end subroutine TransferRecurrence'//ARCSTRING//'P')
                call AddToString(JP)
                call AddToString('Q')
                call AddToString(JQ)
                call AddToString(FromLabel)
                call AddToString('to')
                call AddToString(ToLabel)
                call AddToString(SegLabel(1:iSegLabel))
                call writeString(LUFILE)

                deallocate(TUVINDEX)
                deallocate(TINDEX)
                deallocate(UINDEX)
                deallocate(VINDEX)
                deallocate(JINDEX)
             enddo
          enddo
          WRITE(LUFILE1,'(A)')'end module'
          close(unit = LUFILE1)
       ENDDO
    ENDDO
 ENDDO
!    WRITE(LUFILE1,'(A)')'end module'
!    close(unit = LUSPECIAL)
END subroutine PASSsub

  subroutine SUBROUTINE_MAIN(LUFILE,JMAX,JP,JQ,nTUVP,nTUVQ,JPQ,nTUV,nTUVPLUS,TUVINDEX,&
               & TINDEX,UINDEX,VINDEX,JINDEX,nTUVprev3,nTUVprev2,nTUVprev,IfacX,TUVindexX,&
               & LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
    implicit none
    INTEGER,intent(in) :: LUFILE,JMAX,JP,JQ,nTUVP,nTUVQ,JPQ,nTUV,nTUVPLUS
    integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    integer :: TINDEX(nTUVPLUS),IfacX(:,:),TUVindexX(:,:)
    integer :: UINDEX(nTUVPLUS)
    integer :: VINDEX(nTUVPLUS)
    integer :: JINDEX(nTUVPLUS)
    logical :: LOOPUNROLL,DoOpenACC
    INTEGER,intent(in) :: nTUVprev3,nTUVprev2,nTUVprev
    character(len=20),intent(in) :: PrimLabelAux
    integer,intent(in) :: iPrimLabelAux
    !local
    INTEGER :: ituvP,J,Tp,Up,Vp,N,N2,ituv,C
    INTEGER :: nTUVTMPP,nTUVTMPQ,nTUVTMP,nTUVTMPprev
    integer :: tq,uq,vq,ituvq,ituvqminus1x,ituvqminus1y,ituvqminus1z
    integer :: iTUVplus1x,iTUVplus1y,iTUVplus1z,nTUVTMPQs,I,nTUVTMPQ2,nTUVQ2
    Integer :: MaxAngmomQP,JTMP,iTUVQminus1
    !    logical :: CREATED(-2:8,-2:8,-2:8)
    logical,pointer :: CREATED(:,:,:),ituvqplus1LEnTUVQarray(:)
    logical,pointer :: ituvqminus1LEnTUVQarray(:)
    logical :: TREC,UREC,VREC,TREC2,UREC2,VREC2,DoneCartDir(3),ituvqplus1LEnTUVQ
    logical :: ituvqminus1LEnTUVQ
    integer :: nTUVLIST,nTUVLISTactual,CARTDIR,iTUVPminus1,iTUVPminus2,Tpminus1
    integer,pointer :: TwoTermTUVLIST(:)
    integer :: iTUVQplus1,nLength

    allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
    CREATED  = .FALSE.
    CREATED(0,0,0) = .TRUE.
    
    DO JTMP=0,JP
       !           print*,'!JTMP = ',JTMP
       WRITE(LUFILE,'(A,I2)')   ' ! Building for Angular momentum Jp =',JTMP
       IF(JTMP.EQ.0)THEN
          IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ'
          WRITE(LUFILE,'(A,I3)')'     DO iTUVQ=1,',nTUVQ          
          IF(nPrimLast)THEN
             WRITE(LUFILE,'(A,A,A)')   '      Tmp0(iTUVQ,1) = Aux(iTUVQ,',PrimLabelAux(1:iPrimLabelAux),')'
          ELSE
             WRITE(LUFILE,'(A,A,A)')   '      Tmp0(iTUVQ,1) = Aux(',PrimLabelAux(1:iPrimLabelAux),',iTUVQ)'
          ENDIF
          WRITE(LUFILE,'(A)')   '     ENDDO'
          CYCLE
       ELSE
          !==========================================================================
          !            Build the     Theta(i,0,k,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k-1,0) + (k-1)/(2q)*Theta(i,0,k-2,0)
          !==========================================================================
          nTUVTMPP=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
          nTUVTMPQ=(JPQ-JTMP+1)*(JPQ-JTMP+2)*(JPQ-JTMP+3)/6
          !increase JP=JTMP with 1 
          nTUVTMPP=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
          DO Tp=JTMP,0,-1       
             DO Up=JTMP-Tp,0,-1
                Vp=JTMP-Tp-Up  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1,Tp,Up,Vp,CREATED,JMAX,TUVINDEX)
                IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ'
                WRITE(LUFILE,'(A,I3)')'     do iTUVQ = 1,',nTUVQ
                !Theta(i,0,k,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k-1,0) + (k-1)/(2q)*Theta(i,0,k-2,0)
                CALL LOOPRECURRENCE1(iTUVP,Tpminus1,iTUVPminus2,iTUVPminus1,CARTDIR,LUFILE,PrimLabelAux,iPrimLabelAux)
                WRITE(LUFILE,'(A,I3)')'     enddo'
                !                   CREATED(Tq,Uq,Vq) = .TRUE.
             ENDDO
          ENDDO
          
          DO Tp=JTMP,0,-1       
             DO Up=JTMP-Tp,0,-1
                Vp=JTMP-Tp-Up  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1,Tp,Up,Vp,CREATED,JMAX,TUVINDEX)
                
                IF(nTUVTMPQ.GE.nTUVQ+1)THEN
                   IF(DoOpenACC)WRITE(LUFILE,'(A)')'!$ACC LOOP SEQ'
                   WRITE(LUFILE,'(A,I3,A,I3)')'     do iTUVQ = ',nTUVQ+1,',',nTUVTMPQ !place in tmp array 
                   !Theta(i,0,k,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k-1,0) + (k-1)/(2q)*Theta(i,0,k-2,0)
                   CALL LOOPRECURRENCE2(iTUVP,Tpminus1,iTUVPminus2,iTUVPminus1,CARTDIR,JTMP,LUFILE,PrimLabelAux,iPrimLabelAux)
                   WRITE(LUFILE,'(A,I3,A,I3)')'     enddo'
                ENDIF
                
                !                   CREATED(Tq,Uq,Vq) = .TRUE.
             ENDDO
          ENDDO
          !==========================================================================
          !            Build the     Theta(i,0,k,0) += i/(2q)*Theta(i-1,0,k-1,0) 
          !==========================================================================
          DO Tp=JTMP,0,-1       
             DO Up=JTMP-Tp,0,-1
                Vp=JTMP-Tp-Up  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1,Tp,Up,Vp,CREATED,JMAX,TUVINDEX)
                nTUVQ2 = (JQ+0)*(JQ+1)*(JQ+2)/6   
                iTUVQminus1 = 1
                ituvQminus1LEnTUVQ = .TRUE. 
                CALL WRITERECURRENCE5(CARTDIR,Tp,Up,Vp,1,nTUVQ2,&
                     & nTUVQ2,TUVINDEX,JTMP,JMAX,LUFILE,iTUVQminus1,.TRUE.,ituvQminus1LEnTUVQ,IfacX,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
             ENDDO
          ENDDO

          DO Tp=JTMP,0,-1       
             DO Up=JTMP-Tp,0,-1
                Vp=JTMP-Tp-Up  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1,Tp,Up,Vp,CREATED,JMAX,TUVINDEX)
                IF(nTUVTMPQ.GE.nTUVQ+1)THEN
                   nTUVQ2 = (JQ+0)*(JQ+1)*(JQ+2)/6   
                   nTUVTMPQ2=(JPQ-JTMP)*(JPQ-JTMP+1)*(JPQ-JTMP+2)/6
                   IF(iTUVPminus1.EQ.1)THEN
                      iTUVQminus1 = 1
                      ituvQminus1LEnTUVQ = .TRUE. !dummy
                      CALL WRITERECURRENCE5(CARTDIR,Tp,Up,Vp,nTUVQ2+1,nTUVTMPQ2,&
                           & nTUVQ2,TUVINDEX,JTMP,JMAX,LUFILE,iTUVQminus1,.FALSE.,ituvQminus1LEnTUVQ,IfacX,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                   ELSE
                      IF(nTUVTMPQ2.LE.nTUVQ)THEN
                         ituvqminus1LEnTUVQ = .TRUE.
                         CALL WRITERECURRENCE5(CARTDIR,Tp,Up,Vp,nTUVQ2+1,nTUVTMPQ2,&
                              & nTUVQ2,TUVINDEX,JTMP,JMAX,LUFILE,iTUVQminus1,.FALSE.,ituvQminus1LEnTUVQ,IfacX,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                      ELSE
                         ituvqminus1LEnTUVQ = .TRUE.
                         CALL WRITERECURRENCE5(CARTDIR,Tp,Up,Vp,nTUVQ2+1,nTUVQ,&
                              & nTUVQ2,TUVINDEX,JTMP,JMAX,LUFILE,iTUVQminus1,.FALSE.,ituvQminus1LEnTUVQ,IfacX,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                         ituvqminus1LEnTUVQ = .FALSE.
                         CALL WRITERECURRENCE5(CARTDIR,Tp,Up,Vp,nTUVQ+1,nTUVTMPQ2,&
                              & nTUVQ2,TUVINDEX,JTMP,JMAX,LUFILE,iTUVQminus1,.FALSE.,ituvqminus1LEnTUVQ,IfacX,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                      ENDIF
                   ENDIF
                ENDIF                
             ENDDO
          ENDDO
          !==========================================================================
          !            Build the     Theta(i,0,k,0) += - p/q*Theta(i+1,0,k-1,0) 
          !==========================================================================
          DO Tp=JTMP,0,-1       
             DO Up=JTMP-Tp,0,-1
                Vp=JTMP-Tp-Up  
                CALL DETERMINE_CARTDIR(CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1,Tp,Up,Vp,CREATED,JMAX,TUVINDEX)
                
                !Theta(i,0,k,0) += - p/q*Theta(i+1,0,k-1,0) 
                
                IF(iTUVPminus1.EQ.1)THEN
                   ituvQplus1LEnTUVQ = .TRUE. !dummy argument 
                   CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,1,ituvQplus1LEnTUVQ,nTUVQ,nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                ELSE
                   allocate(ituvqplus1LEnTUVQarray(nTUVQ))
                   DO iTUVQ = 1, nTUVQ
                      Tq = Tindex(iTUVQ) 
                      Uq = Uindex(iTUVQ) 
                      Vq = Vindex(iTUVQ)
                      IF(CARTDIR.EQ.1)iTUVQplus1 = TUVINDEX(Tq+1,Uq,Vq)
                      IF(CARTDIR.EQ.2)iTUVQplus1 = TUVINDEX(Tq,Uq+1,Vq)
                      IF(CARTDIR.EQ.3)iTUVQplus1 = TUVINDEX(Tq,Uq,Vq+1)
                      ituvqplus1LEnTUVQarray(iTUVQ) = ituvqplus1.LE.nTUVQ
                   ENDDO
                   IF(ALL(ituvqplus1LEnTUVQarray))THEN
                      !all true
                      ituvqplus1LEnTUVQ = .TRUE.
                      CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,1,ituvqplus1LEnTUVQ,nTUVQ,&
                           & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                   ELSEIF(COUNT(ituvqplus1LEnTUVQarray).EQ.0)THEN
                      !all false
                      ituvqplus1LEnTUVQ = .FALSE.
                      CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,1,ituvqplus1LEnTUVQ,nTUVQ,&
                           & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                   ELSE                      
                      IF(ALL(ituvqplus1LEnTUVQarray(1:COUNT(ituvqplus1LEnTUVQarray))))THEN
                         !all T are sequential in the beginning
                         ituvqplus1LEnTUVQ = .TRUE.
                         CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,1,ituvqplus1LEnTUVQ,COUNT(ituvqplus1LEnTUVQarray),&
                              & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                         !all F are sequential at the end
                         ituvqplus1LEnTUVQ = .FALSE.
                         CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,COUNT(ituvqplus1LEnTUVQarray)+1,ituvqplus1LEnTUVQ,nTUVQ,&
                              & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                      ELSE
                         print*,'ituvqplus1LEnTUVQarray',ituvqplus1LEnTUVQarray
                         stop 'ERROR'
                      ENDIF
                   ENDIF
                   deallocate(ituvqplus1LEnTUVQarray)
                ENDIF

                IF(nTUVTMPQ.GE.nTUVQ+1)THEN
                   IF(iTUVPminus1.EQ.1)THEN
                      ituvQplus1LEnTUVQ = .TRUE. !dummy argument 
                      CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,nTUVQ+1,ituvqplus1LEnTUVQ,nTUVTMPQ,&
                           & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                   ELSE
                      !Theta(i,0,k,0) += - p/q*Theta(i+1,0,k-1,0) 
                      allocate(ituvqplus1LEnTUVQarray(nTUVQ+1:nTUVTMPQ))
                      DO iTUVQ = nTUVQ+1, nTUVTMPQ
                         Tq = Tindex(iTUVQ) 
                         Uq = Uindex(iTUVQ) 
                         Vq = Vindex(iTUVQ)
                         IF(CARTDIR.EQ.1)iTUVQplus1 = TUVINDEX(Tq+1,Uq,Vq)
                         IF(CARTDIR.EQ.2)iTUVQplus1 = TUVINDEX(Tq,Uq+1,Vq)
                         IF(CARTDIR.EQ.3)iTUVQplus1 = TUVINDEX(Tq,Uq,Vq+1)
                         ituvqplus1LEnTUVQarray(iTUVQ) = ituvqplus1.LE.nTUVQ
                      ENDDO
                      IF(ALL(ituvqplus1LEnTUVQarray))THEN
                         !all true
                         ituvqplus1LEnTUVQ = .TRUE.
                         CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,nTUVQ+1,ituvqplus1LEnTUVQ,nTUVTMPQ,&
                              & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                      ELSEIF(COUNT(ituvqplus1LEnTUVQarray).EQ.0)THEN
                         !all false
                         ituvqplus1LEnTUVQ = .FALSE.
                         CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,nTUVQ+1,ituvqplus1LEnTUVQ,nTUVTMPQ,&
                              & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                      ELSE
                         IF(ALL(ituvqplus1LEnTUVQarray(nTUVQ+1:nTUVQ+COUNT(ituvqplus1LEnTUVQarray))))THEN
                            !all T are sequential in the beginning
                            ituvqplus1LEnTUVQ = .TRUE.
                            CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,nTUVQ+1,ituvqplus1LEnTUVQ,nTUVQ+COUNT(ituvqplus1LEnTUVQarray),&
                                 & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                            !all F are sequential at the end
                            ituvqplus1LEnTUVQ = .FALSE.
                            CALL WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,nTUVQ+COUNT(ituvqplus1LEnTUVQarray)+1,ituvqplus1LEnTUVQ,nTUVTMPQ,&
                                 & nTUVQ,TUVINDEX,JTMP,JMAX,LUFILE,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
                         ELSE
                            print*,'ituvpplus1LEnTUVParray',ituvqplus1LEnTUVQarray
                            stop 'ERROR2'
                         ENDIF
                      ENDIF
                      deallocate(ituvqplus1LEnTUVQarray)
                   ENDIF
                ENDIF
                
                CREATED(Tp,Up,Vp) = .TRUE.
             ENDDO
          ENDDO
       ENDIF
    ENDDO
  end subroutine SUBROUTINE_MAIN

  subroutine DETERMINE_CARTDIR(CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1,Tp,Up,Vp,CREATED,JMAX,TUVINDEX)
    implicit none
    integer,intent(inout) :: CARTDIR,iTUVP,iTUVPminus1,iTUVPminus2,Tpminus1
    integer,intent(in) :: Tp,Up,Vp,JMAX
    logical,intent(in) :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    !
    integer :: N,N2
    logical :: Trec,Urec,Vrec,Trec2,Urec2,Vrec2
    iTUVP = TUVINDEX(Tp,Up,Vp)
    !Theta(i,0,k,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k-1,0) + i/(2q)*Theta(i-1,0,k-1,0) + (k-1)/(2q)*Theta(i,0,k-2,0) - p/q*Theta(i+1,0,k-1,0) 
    
    !choose to use X, Y or Z 
    !how can the (Tp,Up,Vp,Tq,Uq,Vq) be built from lower aka (210) can be built from (110) and (200)
    

    
    !how can the (Tp,Up,Vp) be built from lower aka (210) can be built from (110) and (200)
    !However building (210) from (110) using the X recurrence requires both (110) and (010)   (210) = (110) + 1*(010)
    !building (210) from (200) using the Y recurrence requires both (200)   (210) = (200) + 0*(2-10)
    TREC = CREATED(Tp-1,Up,Vp)     !Example(210): test if (110) is build  (TRUE) 
    UREC = CREATED(Tp,Up-1,Vp)     !Example(210): test if (200) is build  (TRUE) 
    VREC = CREATED(Tp,Up,Vp-1)     !Example(210): test if (21-1) is build (FALSE) 
    N=0
    IF(TREC)N=N+1                  !Example(210): N=1
    IF(UREC)N=N+1                  !Example(210): N=2
    IF(VREC)N=N+1                  !Example(210): N=2
    IF(N.EQ.1)THEN
       !only one possible way to construct it like (200) can only be constructed from (100)
       IF(TREC)THEN
          CARTDIR = 1; iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp); iTUVPminus2 = TUVINDEX(Tp-2,Up,Vp); Tpminus1=Tp-1
       ELSEIF(UREC)THEN
          CARTDIR = 2; iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp); iTUVPminus2 = TUVINDEX(Tp,Up-2,Vp); Tpminus1=Up-1
       ELSEIF(VREC)THEN
          CARTDIR = 3; iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1); iTUVPminus2 = TUVINDEX(Tp,Up,Vp-2); Tpminus1=Vp-1          
       ELSE
          STOP 'TK1'
       ENDIF
    ELSE
       !several ways to construct it. For instance the (210) Example
       TREC2 = CREATED(Tp-2,Up,Vp)    !Example(210): test if (010) is build (TRUE) 
       UREC2 = CREATED(Tp,Up-2,Vp)    !Example(210): test if (2-10) is build (FALSE) 
       VREC2 = CREATED(Tp,Up,Vp-2)    !Example(210): test if (21-1) is build (FALSE) 
       N2=0
       IF(TREC2)N2=N2+1               !Example(210): N2 = 1
       IF(UREC2)N2=N2+1               !Example(210): N2 = 1
       IF(VREC2)N2=N2+1               !Example(210): N2 = 1
       IF(N2.LT.N)THEN                !Example(210): N2 = 1 .LT. N = 2 => TRUE
          !3 term recurrence possible for one or more the possibilities 
          IF(.NOT.(TREC.AND.TREC2).AND.TREC)THEN       !Example(210): TREC.AND.TREC2 = .TRUE. => FALSE
             CARTDIR = 1; iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp); iTUVPminus2 = TUVINDEX(Tp-2,Up,Vp); Tpminus1=Tp-1             
          ELSEIF(.NOT.(UREC.AND.UREC2).AND.UREC)THEN   !Example(210): UREC.AND.UREC2 = .FALSE. => TRUE
             CARTDIR = 2; iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp); iTUVPminus2 = TUVINDEX(Tp,Up-2,Vp); Tpminus1=Up-1          
          ELSEIF(.NOT.(VREC.AND.VREC2).AND.VREC)THEN
             CARTDIR = 3; iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1); iTUVPminus2 = TUVINDEX(Tp,Up,Vp-2); Tpminus1=Vp-1                     
          ENDIF
       ELSE
          !all requires 4 term recurrence possible for one or more the possibilities 
          IF(TREC)THEN
             CARTDIR = 1; iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp); iTUVPminus2 = TUVINDEX(Tp-2,Up,Vp); Tpminus1=Tp-1             
          ELSEIF(UREC)THEN
             CARTDIR = 2; iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp); iTUVPminus2 = TUVINDEX(Tp,Up-2,Vp); Tpminus1=Up-1          
          ELSEIF(VREC)THEN
             CARTDIR = 3; iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1); iTUVPminus2 = TUVINDEX(Tp,Up,Vp-2); Tpminus1=Vp-1          
          ELSE
             STOP 'TK2'
          ENDIF
       ENDIF
    ENDIF
  END subroutine DETERMINE_CARTDIR

  SUBROUTINE LOOPRECURRENCE1(iTUVP,Tpminus1,iTUVPminus2,iTUVPminus1,CARTDIR,&
       & LUPRI,PrimLabelAux,iPrimLabelAux)
    implicit none
    integer :: iTUVP,Tpminus1,iTUVPminus2,iTUVPminus1,CARTDIR,lupri
    character(len=132) :: STRING 
    integer :: iString
    character(len=4) :: DIRECTIONSTRING
    character(len=20),intent(in) :: PrimLabelAux
    integer,intent(in) :: iPrimLabelAux
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
    call AddToString('Tmp0(iTUVQ,')
    call AddToString(iTUVP)
    call AddToString(') = ')
    call AddToString(DIRECTIONSTRING)
    call AddToString('*')
    IF(iTUVPminus1.EQ.1)THEN
       IF(nPrimLast)THEN
          call AddToString('Aux(iTUVQ,')
          call AddToString(PrimLabelAux(1:iPrimLabelAux))
          call AddToString(')')
       ELSE
          call AddToString('Aux(')
          call AddToString(PrimLabelAux(1:iPrimLabelAux))
          call AddToString(',iTUVQ)')
       ENDIF
    ELSE
       call AddToString('Tmp0(iTUVQ,')
       call AddToString(iTUVPminus1)
       call AddToString(')')
    ENDIF
    !possibly add term 3   
    IF(Tpminus1.EQ.1.AND.iTUVPminus2.GT.0)THEN
       IF(iTUVPminus2.EQ.1)THEN
          IF(nPrimLast)THEN
             call AddToString('+ inv2expP*Aux(iTUVQ,')
             call AddToString(PrimLabelAux(1:iPrimLabelAux))
             call AddToString(')')
          ELSE
             call AddToString('+ inv2expP*Aux(')
             call AddToString(PrimLabelAux(1:iPrimLabelAux))
             call AddToString(',iTUVQ)')
          ENDIF
       ELSE
          call AddToString('+ inv2expP*Tmp0(iTUVQ,')
          call AddToString(iTUVPminus2)
          call AddToString(')')
       ENDIF
    ELSEIF(Tpminus1.GT.0.AND.iTUVPminus2.GT.0)THEN
       IF(iTUVPminus2.EQ.1)THEN
          call AddToString('+')
          call AddToString(Tpminus1)
          IF(nPrimLast)THEN
             call AddToString('*inv2expP*Aux(iTUVQ,')
             call AddToString(PrimLabelAux(1:iPrimLabelAux))
             call AddToString(')')
          ELSE
             call AddToString('*inv2expP*Aux(')
             call AddToString(PrimLabelAux(1:iPrimLabelAux))
             call AddToString(',iTUVQ)')
          ENDIF
       ELSE
          call AddToString('+')
          call AddToString(Tpminus1)
          call AddToString('*inv2expP*Tmp0(iTUVQ,')
          call AddToString(iTUVPminus2)
          call AddToString(')')
       ENDIF
    ENDIF
    call writeString(LUPRI)
  END SUBROUTINE LOOPRECURRENCE1
  
  SUBROUTINE LOOPRECURRENCE2(iTUVP,Tpminus1,iTUVPminus2,iTUVPminus1,CARTDIR,JTMP,&
       & LUPRI,PrimLabelAux,iPrimLabelAux)
    implicit none
    integer :: iTUVP,Tpminus1,iTUVPminus2,iTUVPminus1,CARTDIR,JTMP,LUPRI
    character(len=132) :: STRING 
    integer :: iString
    character(len=4) :: DIRECTIONSTRING
    character(len=20),intent(in) :: PrimLabelAux
    integer,intent(in) :: iPrimLabelAux
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
    call AddToString(JTMP)
    call AddToString('(iTUVQ,')
    call AddToString(iTUVP)
    call AddToString(') = ')
    call AddToString(DIRECTIONSTRING)
    call AddToString('*')
    IF(iTUVPminus1.EQ.1)THEN
       IF(nPrimLast)THEN
          call AddToString('Aux(iTUVQ,')
          call AddToString(PrimLabelAux(1:iPrimLabelAux))
          call AddToString(')')
       ELSE
          call AddToString('Aux(')
          call AddToString(PrimLabelAux(1:iPrimLabelAux))
          call AddToString(',iTUVQ)')
       ENDIF
    ELSE
       call AddToString('Tmp')
       call AddToString(JTMP-1)
       call AddToString('(iTUVQ,')
       call AddToString(iTUVPminus1)
       call AddToString(')')
    ENDIF
    !possibly add term 3   
    IF(Tpminus1.EQ.1.AND.iTUVPminus2.GT.0)THEN
       IF(iTUVPminus2.EQ.1)THEN
          IF(nPrimLast)THEN
             call AddToString('+ inv2expP*Aux(iTUVQ,')
             call AddToString(PrimLabelAux(1:iPrimLabelAux))
             call AddToString(')') 
          ELSE
             call AddToString('+ inv2expP*Aux(')
             call AddToString(PrimLabelAux(1:iPrimLabelAux))
             call AddToString(',iTUVQ)') 
          ENDIF
       ELSE
          call AddToString('+ inv2expP*Tmp') 
          call AddToString(JTMP-2) 
          call AddToString('(iTUVQ,') 
          call AddToString(iTUVPminus2) 
          call AddToString(')') 
       ENDIF
    ELSEIF(Tpminus1.GT.0.AND.iTUVPminus2.GT.0)THEN
       IF(iTUVPminus2.EQ.1)THEN
          call AddToString('+') 
          call AddToString(Tpminus1) 
          IF(nPrimLast)THEN
             call AddToString('*inv2expP*Aux(iTUVQ,') 
             call AddToString(PrimLabelAux(1:iPrimLabelAux))
             call AddToString(')') 
          ELSE
             call AddToString('*inv2expP*Aux(') 
             call AddToString(PrimLabelAux(1:iPrimLabelAux))
             call AddToString(',iTUVQ)') 
          ENDIF
       ELSE
          call AddToString('+') 
          call AddToString(Tpminus1) 
          call AddToString('*inv2expP*Tmp') 
          call AddToString(JTMP-2) 
          call AddToString('(iTUVQ,') 
          call AddToString(iTUVPminus2) 
          call AddToString(')') 
       ENDIF
    ENDIF
    !Final step write the string
    call writeString(LUPRI)
!    WRITE(LUFILE,'(A)') STRING(1:iSTRING-1)
  END SUBROUTINE LOOPRECURRENCE2

  subroutine WRITERECURRENCE5(CARTDIR,Tp,Up,Vp,iTUVstart,nTUVQ_tmp,&
       & nTUVQ,TUVINDEX,JTMP,JMAX,LUPRI,iTUVQminus1,iTUVQLEnTUVQ,&
       & ituvqminus1LEnTUVQ,IfacX,TUVindexX,LOOPUNROLL,&
       & PrimLabelAux,iPrimLabelAux,DoOpenACC)
    implicit none
    !A to C
    !Theta(i,0,k,0) = i/(2p)*Theta(i-1,0,k-1,0) 
    integer,intent(in) :: CARTDIR,Tp,Up,Vp,nTUVQ,JTMP,JMAX,iTUVQminus1,nTUVQ_tmp
    integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),LUPRI,iTUVstart
    integer,intent(in) :: IfacX(:,:),TUVindexX(:,:)
    logical,intent(in) :: iTUVQLEnTUVQ,ituvqminus1LEnTUVQ,LOOPUNROLL,DoOpenACC
    character(len=20),intent(in) :: PrimLabelAux
    integer,intent(in) :: iPrimLabelAux
    !
    integer :: iTUVQ,iTUVP,iTUVPplus1,I,iTUVPminus1
    character(len=132) :: STRING 
    integer :: ituvqminus1x,iString
    iTUVP = TUVINDEX(Tp,Up,Vp)
    IF(CARTDIR.EQ.1)THEN
       !X direction
       iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp)
    ELSEIF(CARTDIR.EQ.2)THEN
       !Y direction
       iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp)                        
    ELSEIF(CARTDIR.EQ.3)THEN
       !Z direction
       iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1)
    ENDIF
    !Result(iTUVP,iTUVQ) = I * Term2(iTUVPminus1,iTUVQminus1) + Term4(iTUVPplus1,iTUVQminus1)
    DO ituvqminus1x = iTUVstart,nTUVQ_tmp
       IF(LOOPUNROLL.OR.ituvqminus1x.EQ.iTUVstart)THEN
          IF(LOOPUNROLL)THEN
             iTUVQ = TUVindexX(ituvqminus1x,CARTDIR)
             call initString(5)          
          ELSE
             call initString(5)          
             IF(DoOpenACC)WRITE(LUPRI,'(A)')'!$ACC LOOP SEQ'
             call AddToString('do ituvqminus1 = ')
             call AddToString(iTUVstart)
             call AddToString(',')
             call AddToString(nTUVQ_tmp)
             call writeString(LUPRI)
             
             call initString(6)          
             call AddToString('iTUVQ = TUVindex')
             call AddToString('X')
             call AddToString(CARTDIR)
             call AddToString('(ituvqminus1)')
             call writeString(LUPRI)
             call initString(6)          
          ENDIF
       
          IF(iTUVQLEnTUVQ)THEN
             call AddToString('Tmp0(')
          ELSE
             call AddToString('Tmp')
             call AddToString(JTMP)
             call AddToString('(')
          ENDIF
          IF(LOOPUNROLL)THEN
             call AddToString(iTUVQ)
          ELSEIF(ituvqminus1x.EQ.iTUVstart)THEN
             call AddToString('iTUVQ')
          ENDIF
          call AddToString(',')
          call AddToString(iTUVP)
          call AddToString(') = ')
          
          IF(iTUVQLEnTUVQ)THEN
             call AddToString('Tmp0(')
          ELSE
             call AddToString('Tmp')
             call AddToString(JTMP)
             call AddToString('(')
          ENDIF
          
          IF(LOOPUNROLL)THEN
             call AddToString(iTUVQ)
          ELSEIF(ituvqminus1x.EQ.iTUVstart)THEN
             call AddToString('iTUVQ')
          ENDIF
          call AddToString(',')
          call AddToString(iTUVP)
          call AddToString(') ')
          
          IF(LOOPUNROLL)THEN
             IF(IfacX(ituvqminus1x,CARTDIR).EQ.1)THEN
                call AddToString('+ inv2expP*')
             ELSE
                call AddToString('+ ')
                call AddToString(IfacX(ituvqminus1x,CARTDIR))
                call AddToString('*inv2expP*')
             ENDIF
          ELSEIF(ituvqminus1x.EQ.iTUVstart)THEN
             call AddToString('+ IfacX')
             call AddToString(CARTDIR)
             call AddToString('(ituvqminus1)*inv2expP*')
          ENDIF
          
          IF(iTUVPminus1.EQ.1)THEN
             call AddToString('Aux(')

             IF(.NOT.nPrimLast)THEN
                call AddToString(PrimLabelAux(1:iPrimLabelAux))
                call AddToString(',')
             ENDIF
             IF(LOOPUNROLL)THEN
                call AddToString(ituvqminus1x)
             ELSEIF(ituvqminus1x.EQ.iTUVstart)THEN
                call AddToString('ituvqminus1')
             ENDIF
             IF(nPrimLast)THEN
                call AddToString(',')
                call AddToString(PrimLabelAux(1:iPrimLabelAux))
             ENDIF
             call AddToString(') ') 
          ELSE
             IF(ituvQminus1LEnTUVQ)THEN
                call AddToString('Tmp0(')
             ELSE
                call AddToString('Tmp')
                call AddToString(JTMP-1)
                call AddToString('(')
             ENDIF
             IF(LOOPUNROLL)THEN
                call AddToString(ituvqminus1x)
             ELSEIF(ituvqminus1x.EQ.iTUVstart)THEN
                call AddToString('ituvqminus1')
             ENDIF
             call AddToString(',')
             call AddToString(iTUVPminus1)
             call AddToString(') ')
          ENDIF
          IF(LOOPUNROLL)THEN
             call writeString(LUPRI)
          ELSEIF(ituvqminus1x.EQ.iTUVstart)THEN
             call writeString(LUPRI)
          ENDIF
       ENDIF
    ENDDO
    IF(.NOT.LOOPUNROLL)THEN
       call initString(5)          
       call AddToString('enddo')
       call writeString(LUPRI)
    ENDIF
  END subroutine WRITERECURRENCE5

subroutine WRITERECURRENCE4(CARTDIR,Tp,Up,Vp,iTUVstart,ituvqplus1LEnTUVQ,nTUVQ_tmp,&
     & nTUVQ,TUVINDEX,JTMP,JMAX,LUPRI,TUVindexX,LOOPUNROLL,PrimLabelAux,iPrimLabelAux,DoOpenACC)
  implicit none
  !A to C
  !Theta(i,0,k,0) = i/(2q)*Theta(i-1,0,k-1,0) - p/q*Theta(i+1,0,k-1,0) 
  integer,intent(in) :: CARTDIR,Tp,Up,Vp,nTUVQ,JTMP,JMAX,iTUVstart,nTUVQ_tmp
  integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),LUPRI,TUVindexX(:,:)
  logical,intent(in) :: ituvqplus1LEnTUVQ,LOOPUNROLL,DoOpenACC
  character(len=20),intent(in) :: PrimLabelAux
  integer,intent(in) :: iPrimLabelAux
  !
  integer :: iTUVQ,iTUVP,iTUVQplus1,iTUVQminus1,I,iTUVPminus1
  character(len=132) :: STRING 
  integer :: iString
  iTUVP = TUVINDEX(Tp,Up,Vp)
  IF(CARTDIR.EQ.1)THEN
     !X direction
     iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp)
  ELSEIF(CARTDIR.EQ.2)THEN
     !Y direction
     iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp)                        
  ELSEIF(CARTDIR.EQ.3)THEN
     !Z direction
     iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1)
  ENDIF
  !Result(iTUVP,iTUVQ) = I * Term2(iTUVPminus1,iTUVQminus1) + Term4(iTUVPplus1,iTUVQminus1)

  !determine if part of TUVindexX(iTUVstart:nTUVQ_tmp,CARTDIR) is sequential - def as more than 8.
  call DetermineSequentialRegion(iTUVstart,nTUVQ_tmp,CARTDIR,TUVindexX)

  DO ituvQ = iTUVstart,nTUVQ_tmp
     IF(LOOPUNROLL)THEN
        iTUVQplus1 = TUVindexX(ituvQ,CARTDIR)
        call initString(5)          
     ELSEIF(ituvQ.EQ.iTUVstart)THEN
        call initString(5)          
        IF(DoOpenACC)WRITE(LUPRI,'(A)')'!$ACC LOOP SEQ'
        call AddToString('do iTUVQ = ')
        call AddToString(iTUVstart)
        call AddToString(',')
        call AddToString(nTUVQ_tmp)
        call writeString(LUPRI)

        call initString(6)          
        call AddToString('iTUVplus1 = TUVindexX')
        call AddToString(CARTDIR)
        call AddToString('(iTUVQ)')
        call writeString(LUPRI)
        call initString(6)          
     ENDIF

     IF(LOOPUNROLL.OR.ituvQ.EQ.iTUVstart)THEN
        IF(iTUVstart.LE.nTUVQ)THEN
           call AddToString('Tmp0(')
        ELSE
           call AddToString('Tmp')
           call AddToString(JTMP)
           call AddToString('(')
        ENDIF
     ENDIF

     IF(LOOPUNROLL)THEN
        call AddToString(iTUVQ)
     ELSEIF(ituvQ.EQ.iTUVstart)THEN
        call AddToString('iTUVQ')
     ENDIF
     IF(LOOPUNROLL.OR.ituvQ.EQ.iTUVstart)THEN
        call AddToString(',')
        call AddToString(iTUVP)
        call AddToString(')')
        call AddToString(' = ')
        IF(iTUVstart.LE.nTUVQ)THEN
           call AddToString('Tmp0(')
        ELSE
           call AddToString('Tmp')
           call AddToString(JTMP)
           call AddToString('(')
        ENDIF
     ENDIF
     IF(LOOPUNROLL)THEN
        call AddToString(iTUVQ)
     ELSEIF(ituvQ.EQ.iTUVstart)THEN
        call AddToString('iTUVQ')
     ENDIF
     IF(LOOPUNROLL.OR.ituvQ.EQ.iTUVstart)THEN
        call AddToString(',')
        call AddToString(iTUVP)
        call AddToString(')')
        call AddToString(' + ')
        call AddToString('qinvp*')
        IF(iTUVPminus1.EQ.1)THEN
           call AddToString('Aux(')
           IF(.NOT.nPrimLast)THEN
              call AddToString(PrimLabelAux(1:iPrimLabelAux))
              call AddToString(',') 
           ENDIF
           IF(LOOPUNROLL)THEN
              call AddToString(iTUVQplus1)
           ELSEIF(ituvQ.EQ.iTUVstart)THEN
              call AddToString('iTUVplus1')
           ENDIF
!           call AddToString(',IP)')
           IF(nPrimLast)THEN
              call AddToString(',') 
              call AddToString(PrimLabelAux(1:iPrimLabelAux))
           ENDIF
           call AddToString(')') 
        ELSE
           IF(ituvqplus1LEnTUVQ)THEN
              call AddToString('Tmp0')
           ELSE
              call AddToString('Tmp')
              call AddToString(JTMP-1)
           ENDIF
           call AddToString('(')
           IF(LOOPUNROLL)THEN
              call AddToString(iTUVQplus1)
           ELSEIF(ituvQ.EQ.iTUVstart)THEN
              call AddToString('iTUVplus1')
           ENDIF
           call AddToString(',')
           call AddToString(iTUVPminus1)
           call AddToString(')')
        ENDIF
        call writeString(LUPRI)
     ENDIF
  enddo
  IF(.NOT.LOOPUNROLL)THEN
     WRITE(LUPRI,'(A)')'     enddo'
  ENDIF
END subroutine WRITERECURRENCE4

subroutine DetermineSequentialRegion(iTUVstart,nTUVQ_tmp,CARTDIR,TUVindexX)
  implicit none
  integer :: iTUVstart,nTUVQ_tmp,CARTDIR,TUVindexX(:,:)
  !Local
  integer :: nSegMax,nseg,iOld,ituvQ,iTUVQplus1,NumberSeqRegions
  logical :: InsidesequentialRegion
  nSegMax = 0
  nSeg = 0
  NumberSeqRegions = 0
  iOlD = -1
  InsidesequentialRegion = .FALSE.
  DO ituvQ = iTUVstart,nTUVQ_tmp
     iTUVQplus1 = TUVindexX(ituvQ,CARTDIR)
     IF(InsidesequentialRegion)THEN
        IF(iTUVQplus1.EQ.iOLD+1)THEN
           nSeg = nSeg + 1
        ELSE
           nSegMax = MAX(nSeg,nSegMax)
           nSeg = 0 
           NumberSeqRegions = NumberSeqRegions + 1
           InsidesequentialRegion = .FALSE.
        ENDIF
     ELSE
        IF(iTUVQplus1.EQ.iOLD+1)THEN
           nSeg = nSeg + 1
           InsidesequentialRegion = .TRUE.
        ENDIF
     ENDIF
     iOld = iTUVQplus1
  ENDDO
  IF(nSegMax.GT.5)THEN
!     print*,'nSegMax',nSegMax
!     print*,'ARRAY: ',TUVindexX(iTUVstart:nTUVQ_tmp,CARTDIR)
!     print*,'NumberSeqRegions',NumberSeqRegions
  ENDIF
  


end subroutine DetermineSequentialRegion


!!$  subroutine WRITERECURRENCE1(CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,TUVINDEX,JTMQ,JMAX,LUPRI)
!!$    implicit none
!!$    !A to C
!!$    !Theta(i,0,k,0) = i/(2q)*Theta(i-1,0,k-1,0) - p/q*Theta(i+1,0,k-1,0) 
!!$    integer,intent(in) :: CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,JTMQ,JMAX
!!$    integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),LUPRI
!!$    !
!!$    integer :: iTUVQ,iTUVP,iTUVPplus1,iTUVPminus1,I,iTUVQminus1
!!$    character(len=132) :: STRING 
!!$    integer :: iString
!!$    iTUVP = TUVINDEX(Tp,Up,Vp)
!!$    iTUVQ = TUVINDEX(Tq,Uq,Vq)
!!$    IF(CARTDIR.EQ.1)THEN
!!$       !X direction
!!$       iTUVPplus1 = TUVINDEX(Tp+1,Up,Vp)
!!$       iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp)
!!$       I = Tp
!!$       iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq)
!!$    ELSEIF(CARTDIR.EQ.2)THEN
!!$       !Y direction
!!$       iTUVPplus1 = TUVINDEX(Tp,Up+1,Vp)
!!$       iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp)
!!$       I = Up
!!$       iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq)                        
!!$    ELSEIF(CARTDIR.EQ.3)THEN
!!$       !Z direction
!!$       iTUVPplus1 = TUVINDEX(Tp,Up,Vp+1)
!!$       iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1)
!!$       I = Vp
!!$       iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1)
!!$    ENDIF
!!$    !Result(iTUVP,iTUVQ) = I * Term2(iTUVPminus1,iTUVQminus1) + Term4(iTUVPplus1,iTUVQminus1)
!!$
!!$    !step 1 add blanks
!!$    IF(I.GT.0)THEN
!!$       call initString(5)
!!$       !step 2 determine where to put the result
!!$       IF(iTUVP.LE.nTUVP)THEN
!!$          call AddToString('Tmp0(')
!!$       ELSE
!!$          call AddToString('Tmp')
!!$          call AddToString(JTMQ)
!!$          call AddToString('(')
!!$       ENDIF
!!$       call AddToString(iTUVP)
!!$       call AddToString(',')
!!$       call AddToString(iTUVQ)
!!$       call AddToString(') = ')
!!$
!!$       IF(iTUVP.LE.nTUVP)THEN
!!$          call AddToString('Tmp0(')
!!$       ELSE
!!$          call AddToString('Tmp')
!!$          call AddToString(JTMQ)
!!$          call AddToString('(')
!!$       ENDIF
!!$       call AddToString(iTUVP)
!!$       call AddToString(',')
!!$       call AddToString(iTUVQ)
!!$       call AddToString(') ')
!!$
!!$       !=====================
!!$       !step 4 determine if the second term: 
!!$       !  I*inv2expQ*Aux(',ituvpminus1x,',IP) 
!!$       !  I*inv2expQ*Tmp',JTMQ-1,'(',ituvpminus1x,',',iTUVQminus1x,')
!!$       !should be included and if it should use Aux or Tmp 
!!$       IF(I.GT.2)THEN
!!$          call AddToString('+ ')
!!$          call AddToString(I)
!!$          call AddToString('*inv2expQ*')
!!$       ELSEIF(I.EQ.2)THEN
!!$          call AddToString('+ invexpQ*')
!!$       ELSEIF(I.EQ.1)THEN
!!$          call AddToString('+ inv2expQ*')
!!$       ELSE
!!$          !do not include this term 
!!$       ENDIF
!!$
!!$       IF(iTUVQminus1.EQ.1)THEN
!!$          call AddToString('Aux(')
!!$       ELSE
!!$          IF(ituvpminus1.LE.nTUVP)THEN
!!$             call AddToString('Tmp0(')
!!$          ELSE
!!$             call AddToString('Tmp')
!!$             call AddToString(JTMQ-1)
!!$             call AddToString('(')
!!$          ENDIF
!!$       ENDIF
!!$       call AddToString(ituvpminus1)
!!$       call AddToString(',')
!!$       IF(iTUVQminus1.EQ.1)THEN
!!$          call AddToString('IP) ')
!!$       ELSE
!!$          IF(ituvpminus1.LE.nTUVP)THEN
!!$             call AddToString(iTUVQminus1)
!!$             call AddToString(') ')
!!$          ELSE
!!$             call AddToString(iTUVQminus1)
!!$             call AddToString(') ')
!!$          ENDIF
!!$       ENDIF
!!$       call writeString(LUPRI)
!!$    ELSE
!!$       !do not include this term       
!!$    ENDIF
!!$  END subroutine WRITERECURRENCE1
!!$
!!$subroutine WRITERECURRENCE2(CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,TUVINDEX,JTMQ,JMAX,LUPRI)
!!$  implicit none
!!$  !A to C
!!$  !Theta(i,0,k,0) = i/(2q)*Theta(i-1,0,k-1,0) - p/q*Theta(i+1,0,k-1,0) 
!!$  integer,intent(in) :: CARTDIR,Tq,Uq,Vq,Tp,Up,Vp,nTUVP,JTMQ,JMAX
!!$  integer,intent(in) :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),LUPRI
!!$  !
!!$  integer :: iTUVQ,iTUVP,iTUVPplus1,iTUVPminus1,I,iTUVQminus1
!!$  character(len=132) :: STRING 
!!$  integer :: iString
!!$  iTUVP = TUVINDEX(Tp,Up,Vp)
!!$  iTUVQ = TUVINDEX(Tq,Uq,Vq)
!!$  IF(CARTDIR.EQ.1)THEN
!!$     !X direction
!!$     iTUVPplus1 = TUVINDEX(Tp+1,Up,Vp)
!!$     iTUVPminus1 = TUVINDEX(Tp-1,Up,Vp)
!!$     I = Tp
!!$     iTUVQminus1 = TUVINDEX(Tq-1,Uq,Vq)
!!$  ELSEIF(CARTDIR.EQ.2)THEN
!!$     !Y direction
!!$     iTUVPplus1 = TUVINDEX(Tp,Up+1,Vp)
!!$     iTUVPminus1 = TUVINDEX(Tp,Up-1,Vp)
!!$     I = Up
!!$     iTUVQminus1 = TUVINDEX(Tq,Uq-1,Vq)                        
!!$  ELSEIF(CARTDIR.EQ.3)THEN
!!$     !Z direction
!!$     iTUVPplus1 = TUVINDEX(Tp,Up,Vp+1)
!!$     iTUVPminus1 = TUVINDEX(Tp,Up,Vp-1)
!!$     I = Vp
!!$     iTUVQminus1 = TUVINDEX(Tq,Uq,Vq-1)
!!$  ENDIF
!!$  !Result(iTUVP,iTUVQ) = I * Term2(iTUVPminus1,iTUVQminus1) + Term4(iTUVPplus1,iTUVQminus1)
!!$
!!$  !step 1 add blanks
!!$  call initString(5)
!!$  !step 2 determine where to put the 
!!$  IF(iTUVP.LE.nTUVP)THEN
!!$     call AddToString('Tmp0(')
!!$!     call AddToString('Aux2(')
!!$  ELSE
!!$     call AddToString('Tmp')
!!$     call AddToString(JTMQ)
!!$     call AddToString('(')
!!$  ENDIF
!!$  call AddToString(iTUVP)
!!$  call AddToString(',')
!!$  call AddToString(iTUVQ)
!!$!  IF(iTUVP.LE.nTUVP)THEN
!!$!     call AddToString(',IP) = ')
!!$!  ELSE
!!$     call AddToString(') = ')
!!$!  ENDIF
!!$!====================
!!$
!!$  IF(iTUVP.LE.nTUVP)THEN
!!$     call AddToString('Tmp0(')
!!$  ELSE
!!$     call AddToString('Tmp')
!!$     call AddToString(JTMQ)
!!$     call AddToString('(')
!!$  ENDIF
!!$  call AddToString(iTUVP)
!!$  call AddToString(',')
!!$  call AddToString(iTUVQ)
!!$!  IF(iTUVP.LE.nTUVP)THEN
!!$!     call AddToString(',IP) ')
!!$!  ELSE
!!$     call AddToString(') ')
!!$!  ENDIF
!!$
!!$!=====================
!!$  !step 6 determine if the third term: 
!!$  !  pinvq*Aux(',ituvpplus1x,',IP)'
!!$  !  pinvq*Tmp',JTMQ-1,'(',ituvpplus1x,',',iTUVQminus1x,')'
!!$  !should be included and if it should use Aux or Tmp 
!!$  call AddToString('+ qinvp*')
!!$  IF(iTUVQminus1.EQ.1)THEN
!!$     call AddToString('Aux(')
!!$  ELSE
!!$     IF(ituvpplus1.LE.nTUVP)THEN
!!$        call AddToString('Tmp0(')
!!$     ELSE
!!$        call AddToString('Tmp')
!!$        call AddToString(JTMQ-1)
!!$        call AddToString('(')
!!$     ENDIF
!!$  ENDIF
!!$  call AddToString(ituvpplus1)
!!$  call AddToString(',')
!!$  IF(iTUVQminus1.EQ.1)THEN
!!$     call AddToString('IP) ')
!!$  ELSE
!!$     call AddToString(iTUVQminus1)
!!$!     IF(ituvpplus1.LE.nTUVP)THEN
!!$!        call AddToString(',IP) ')
!!$!     ELSE
!!$        call AddToString(') ')
!!$!     ENDIF
!!$  ENDIF
!!$  !Final step write the string
!!$  call writeString(LUPRI)
!!$!  WRITE(LUFILE,'(A)') STRING(1:iSTRING-1)
!!$END subroutine WRITERECURRENCE2

end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
