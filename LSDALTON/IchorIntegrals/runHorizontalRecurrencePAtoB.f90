MODULE TESTMODULE
CONTAINS
  subroutine PASSsub
    IMPLICIT NONE
    INTEGER :: JMAX,nTUV,nTUVprev,ituvP,J,Tp,Up,Vp,N,N2,ntuvP,ituv,C
    INTEGER :: nTUVQ,nTUVTMPP,nTUVTMPQ,JPQ,JP,JQ,nTUVTMP,nTUVTMPprev
    integer :: tq,uq,vq,ituvq,ituvqminus1x,ituvqminus1y,ituvqminus1z
    integer :: iTUVplus1x,iTUVplus1y,iTUVplus1z,nTUVTMPBprev
    integer :: MaxAngmomP, NTUVMAX, AngmomA, AngmomB, NTUVAspec, NTUVBspec
    Integer :: NTUVAstart, NTUVBstart,Jb,Jab,nTUVA,nTUVB,tb,ub,vb
    Integer :: MaxAngmomQP,nTUVplus,JTMP,ntuvprev2,ntuvprev3
    Integer :: nTUVTMPB,nTUVTMPA
    !    logical :: CREATED(-2:8,-2:8,-2:8)
    logical,pointer :: CREATED(:,:,:)
    logical :: TREC,UREC,VREC,TREC2,UREC2,VREC2
    integer,pointer :: TUVINDEX(:,:,:)
    integer,pointer :: TINDEX(:)
    integer,pointer :: UINDEX(:)
    integer,pointer :: VINDEX(:)
    integer,pointer :: JINDEX(:)
    integer :: nTUVLIST,nTUVLISTactual
    integer,pointer :: TwoTermTUVLIST(:)
    character(len=3) :: ARCSTRING
    integer :: GPUrun,MaxAngmomSingle,I,LUFILE
    logical :: DoOpenMP,DoOpenACC,CPU
    Character(len=48) :: FileName    

DO GPUrun = 1,2
    CPU = .TRUE.
    IF(GPUrun.EQ.2)CPU = .FALSE.
    DoOpenMP = .FALSE.
    DoOpenACC = .FALSE.
    IF(CPU)DoOpenMP = .TRUE.
    IF(.NOT.CPU)DoOpenACC = .TRUE.
    IF(CPU)THEN
       ARCSTRING = 'CPU'
    ELSE
       ARCSTRING = 'GPU'
    ENDIF

    DO I =1,48
       FileName(I:I) = ' '
    ENDDO
    WRITE(FileName,'(4A)')'runHorizontalRecurrence'//ARCSTRING//'LHSModAtoB.F90'
    print*,'FileName:',FileName
    LUFILE = 12
    open(unit = LUFILE, file=TRIM(FileName),status="unknown")

    WRITE(LUFILE,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceLHSModAtoB'
    WRITE(LUFILE,'(A)')' use IchorPrecisionModule'
    WRITE(LUFILE,'(A)')'  '
    WRITE(LUFILE,'(A)')' CONTAINS'


    MaxAngmomP = 4 !D currently
    MaxAngmomSingle = 2 !D currently
!    IF(GPUrun.EQ.2)WRITE(LUFILE,'(A)')'#ifdef VAR_OPENACC'
    DO JMAX=0,MaxAngmomP

       nTUV = (JMAX+1)*(JMAX+2)*(JMAX+3)/6   
       nTUVPLUS = (JMAX+2)*(JMAX+3)*(JMAX+4)/6   
       allocate(TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
       TUVINDEX = 0
       allocate(TINDEX(nTUVPLUS))
       allocate(UINDEX(nTUVPLUS))
       allocate(VINDEX(nTUVPLUS))
       allocate(JINDEX(nTUVPLUS))

       ituv = 0 
       TUVINDEX = 0
       DO J = 0, JMAX
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

       JP = JMAX          
       nTUVP = (JP+1)*(JP+2)*(JP+3)/6   
       NTUVMAX = (JP+1)*(JP+2)*(JP+3)/6
       DO AngmomA = 0,JP
          AngmomB = JP - AngmomA
          IF(AngmomB.GT.AngmomA)CYCLE
          IF(AngmomA.GT.MaxAngmomSingle)CYCLE
          IF(AngmomB.GT.MaxAngmomSingle)CYCLE

          NTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
          NTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6

          NTUVAspec = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6 - &
               &         (AngmomA)*(AngmomA+1)*(AngmomA+2)/6
          NTUVBspec = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6 - &
               &         (AngmomB)*(AngmomB+1)*(AngmomB+2)/6
          NTUVAstart = (AngmomA)*(AngmomA+1)*(AngmomA+2)/6
          NTUVBstart = (AngmomB)*(AngmomB+1)*(AngmomB+2)/6
          
          IF(AngmomA.EQ.0.AND.AngmomB.EQ.0)CYCLE
          WRITE(LUFILE,'(A)')''
          IF(JP.LT.10)THEN
             WRITE(LUFILE,'(A,I1,A,I1,A,I1,A)')'subroutine HorizontalRR_'//ARCSTRING//'_LHS_P',JP,'A',AngmomA,'B',AngmomB,'AtoB(nContQP,nPasses,nTUVQ,&'
          ELSE
             WRITE(LUFILE,'(A,I2,A,I1,A,I1,A)')'subroutine HorizontalRR_'//ARCSTRING//'_LHS_P',JP,'A',AngmomA,'B',AngmomB,'AtoB(nContQP,nPasses,nTUVQ,&'
          ENDIF
          WRITE(LUFILE,'(A)')'         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)'
          WRITE(LUFILE,'(A)')'  implicit none'
          WRITE(LUFILE,'(A)')'  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB'
          WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)'
          WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: AuxCont(',nTUVP,',nTUVQ*nContQP*nPasses)'
          WRITE(LUFILE,'(A)')'  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
          !             WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(',nTUVAspec,',',nTUVBspec,',nTUVQ,nContPasses)'
          IF(nTUVBstart+1.EQ.1.AND.nTUVB.EQ.1)THEN
             WRITE(LUFILE,'(A,I5,A1,I5,A)')'  real(realk),intent(inout) :: ThetaP(',NTUVAstart+1,':',nTUVA,',1,nTUVQ*nContQP*nPasses)'
          ELSE
             WRITE(LUFILE,'(A,I5,A1,I5,A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(',NTUVAstart+1,':',nTUVA,',',nTUVBstart+1,':',nTUVB,',nTUVQ*nContQP*nPasses)'
          ENDIF
          WRITE(LUFILE,'(A)')'  !Local variables'
          IF(ANGMOMB.NE.0)THEN
             WRITE(LUFILE,'(A)')'  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB'
             WRITE(LUFILE,'(A)')'  real(realk) :: Xab,Yab,Zab'
          ELSE
             WRITE(LUFILE,'(A)')'  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB'
          ENDIF
          allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
          CREATED  = .FALSE.
          CREATED(0,0,0) = .TRUE.
          JAB = ANGMOMA + ANGMOMB
          JB = ANGMOMB
          DO JTMP=1,JB-1
             nTUVTMPBprev=(JTMP)*(JTMP+1)*(JTMP+2)/6
             nTUVTMPB=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
             nTUVTMPA=(Jab-JTMP+1)*(Jab-JTMP+2)*(Jab-JTMP+3)/6
             if(JTMP.LT.10)THEN
                WRITE(LUFILE,'(A,I1,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVAstart+1,':',nTUVTMPA,',',nTUVTMPBprev+1,':',nTUVTMPB,')'
             else
                WRITE(LUFILE,'(A,I2,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVAstart+1,':',nTUVTMPA,',',nTUVTMPBprev+1,':',nTUVTMPB,')'
             endif
          ENDDO
          IF(JB.GT.1)THEN
             WRITE(LUFILE,'(A)')'!  real(realk) :: Tmp(nTUVA,nTUVB) ordering'
          ENDIF
          IF(DoOpenMP)THEN
!             WRITE(LUFILE,'(A)')'!$OMP PARALLEL DO DEFAULT(none) &'
             WRITE(LUFILE,'(A)')'!$OMP DO &'
             WRITE(LUFILE,'(A)')'!$OMP PRIVATE(iP,&'
             IF(JB.NE.0)THEN
                DO JTMP=1,JB-1
                   if(JTMP.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A)')'!$OMP         Tmp',JTMP,',&'
                   else
                      WRITE(LUFILE,'(A,I2,A)')'!$OMP         Tmp',JTMP,',&'
                   endif
                ENDDO
                WRITE(LUFILE,'(A)')'!$OMP         iPassP,iTUVA,iAtomA,iAtomB,Xab,Yab,Zab) '
             ELSE
                WRITE(LUFILE,'(A)')'!$OMP         iTUVA) '
             ENDIF
!             WRITE(LUFILE,'(A)')'!$OMP SHARED(nTUVQ,nContQP,nPasses,&'
!             IF(JB.NE.0)THEN
!                WRITE(LUFILE,'(A)')'!$OMP         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP)'
!             ELSE
!                WRITE(LUFILE,'(A)')'!$OMP         AuxCont,ThetaP)'
!             ENDIF
          ENDIF
          IF(DoOpenACC)THEN
             WRITE(LUFILE,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUFILE,'(A)')'!$ACC PRIVATE(iP,&'
             IF(JB.NE.0)THEN
                DO JTMP=1,JB-1
                   if(JTMP.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A)')'!$ACC         Tmp',JTMP,',&'
                   else
                      WRITE(LUFILE,'(A,I2,A)')'!$ACC         Tmp',JTMP,',&'
                   endif
                ENDDO
                WRITE(LUFILE,'(A)')'!$ACC         iPassP,iTUVA,iAtomA,iAtomB,Xab,Yab,Zab) &'
             ELSE
                WRITE(LUFILE,'(A)')'!$ACC         iTUVA) &'
             ENDIF
             WRITE(LUFILE,'(A)')'!$ACC PRESENT(nTUVQ,nContQP,nPasses,&'
             IF(JB.NE.0)THEN
                WRITE(LUFILE,'(A)')'!$ACC         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP)'
             ELSE
                WRITE(LUFILE,'(A)')'!$ACC         AuxCont,ThetaP)'
             ENDIF
          ENDIF
             !          WRITE(LUFILE,'(A)')'  DO iP = 1,nContQP*nPasses'
!          WRITE(LUFILE,'(A)')'   DO iTUVQ = 1,nTUVQ'
!             WRITE(LUFILE,'(A)')'    iPassP = (iP-1)/nContQP+1'
!DO i123 = 1,n1*n2*n3'
!i1 = mod(I123-1,n1)+1'
!i2 = mod((I123-(mod(I123-1,nPrim1)+1))/n1,n2)+1'
!i3 = (I123-1)/(n1*n2) + 1'
          WRITE(LUFILE,'(A)')'  DO iP = 1,nTUVQ*nContQP*nPasses'
          IF(JB.NE.0)THEN
             WRITE(LUFILE,'(A)')'!    iTUVQ = mod(IP-1,nTUVQ)+1'
             WRITE(LUFILE,'(A)')'!    iContQP = mod((IP-(mod(IP-1,nTUVQ)+1))/nTUVQ,nContQP)+1'
             WRITE(LUFILE,'(A)')'    iPassP = (IP-1)/(nTUVQ*nContQP) + 1'
             WRITE(LUFILE,'(A)')'    iAtomA = iAtomApass(iPassP)'
             WRITE(LUFILE,'(A)')'    iAtomB = iAtomBpass(iPassP)'
             WRITE(LUFILE,'(A)')'    Xab = Pdistance12(1,iAtomA,iAtomB)'
             WRITE(LUFILE,'(A)')'    Yab = Pdistance12(2,iAtomA,iAtomB)'
             WRITE(LUFILE,'(A)')'    Zab = Pdistance12(3,iAtomA,iAtomB)'
          ENDIF

          DO JTMP=0,JB
             !           print*,'!JTMP = ',JTMP
             IF(JTMP.EQ.0)THEN
                !nTUVTMPP=(Jab-JTMP)*(Jab-JTMP+1)*(Jab-JTMP+2)/6
                IF(nTUVBstart+1.EQ.1)THEN
                   WRITE(LUFILE,'(A,I3,A,I3)')'     DO iTUVA=',NTUVAstart+1,',',nTUVA
                   WRITE(LUFILE,'(A)')   '        ThetaP(iTUVA,1,IP) = AuxCont(iTUVA,IP)'
                   WRITE(LUFILE,'(A)')   '     ENDDO'
                ENDIF
                CYCLE
             ELSE
                nTUVTMPB=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
                nTUVTMPA=(JAB-JTMP+1)*(JAB-JTMP+2)*(JAB-JTMP+3)/6
                !slowly increase JB=JTMP with 1 
                nTUVTMPB=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
                DO Tb=JTMP,0,-1       
                   DO Ub=JTMP-Tb,0,-1
                      Vb=JTMP-Tb-Ub                
                      !                   print*,'!Tq,Uq,Vq = ',Tq,Uq,Vq
                      !step 1.
                      !how can the (Tp,Up,Vp) be built from lower
                      TREC = CREATED(Tb-1,Ub,Vb)
                      UREC = CREATED(Tb,Ub-1,Vb)
                      VREC = CREATED(Tb,Ub,Vb-1)
                      N=0
                      IF(TREC)N=N+1
                      IF(UREC)N=N+1
                      IF(VREC)N=N+1
                      IF(TREC)THEN
                         !                         print*,'!A TRECURRENCE iTUV',iTUV
                         call TRECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB,lufile)
                      ELSEIF(UREC)THEN
                         !                         print*,'!A URECURRENCE iTUV',iTUV
                         call URECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB,lufile)
                      ELSEIF(VREC)THEN
                         !                         print*,'!A VRECURRENCE iTUV',iTUV
                         call VRECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB,lufile)
                      ELSE
                         STOP 'TK1'
                      ENDIF
                      CREATED(Tb,Ub,Vb) = .TRUE.
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
          deallocate(CREATED)
!          WRITE(LUFILE,'(A)')'   ENDDO'
          WRITE(LUFILE,'(A)')'  ENDDO'
!          IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'
          IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
          IF(JP.LT.10)THEN
             WRITE(LUFILE,'(A,I1,A,I1,A,I1,A)')'end subroutine HorizontalRR_'//ARCSTRING//'_LHS_P',JP,'A',AngmomA,'B',AngmomB,'AtoB'
          ELSE
             WRITE(LUFILE,'(A,I2,A,I1,A,I1,A)')'end subroutine HorizontalRR_'//ARCSTRING//'_LHS_P',JP,'A',AngmomA,'B',AngmomB,'AtoB'
          ENDIF
       enddo
       deallocate(TUVINDEX)
       deallocate(TINDEX)
       deallocate(UINDEX)
       deallocate(VINDEX)
       deallocate(JINDEX)
    enddo
!    IF(GPUrun.EQ.2)WRITE(LUFILE,'(A)')'#endif'
    WRITE(LUFILE,'(A)')'end module'
    close(unit = LUFILE)
 enddo
END subroutine PASSsub
  
  subroutine TRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ,lufile)
    implicit none
    integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVQminus1x,iTUVQminus1y
    integer :: iTUVQminus1z,nTUVP,JQ,TMQ,TMP,JP,nTUVASTART,nTUVQSTART,nTUVQ
    integer :: iTUVQminus2x, iTUVP,nTUVTMPP,iTUVminus1x,iTUVplus1x,JMAX
    integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    integer :: TINDEX(:)
    integer :: UINDEX(:)
    integer :: VINDEX(:)
    integer :: JINDEX(:)
    logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    character(len=132) :: STRING 
    integer :: iString,lufile
    !Tb = Tq
    iTUVQ = TUVINDEX(Tq,Uq,Vq)
    iTUVQminus1x = TUVINDEX(Tq-1,Uq,Vq)
    TMQ = Tq-1
    do iTUVP = nTUVASTART+1,nTUVTMPP
       Tp = Tindex(iTUVp) 
       Up = Uindex(iTUVp) 
       Vp = Vindex(iTUVp)
       TMP = Tp
       iTUVplus1x = TUVindex(Tp+1,Up,Vp)
       Jp = Jindex(iTUVp)   

       call XYZRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,&
            & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,nTUVQSTART,nTUVQ,&
            & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Xab',lufile)

    ENDDO

  end subroutine TRECURRENCE

  subroutine XYZRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,&
       & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,nTUVQSTART,nTUVQ,&
       & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,DIRSTRING,LUFILE)
    implicit none
    !step 1 add blanks
    integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVQminus1x,iTUVQminus1y
    integer :: iTUVQminus1z,nTUVP,JQ,TMQ,TMP,JP,nTUVASTART,nTUVQSTART,nTUVQ
    integer :: iTUVQminus2x, iTUVP,nTUVTMPP,iTUVminus1x,iTUVplus1x,JMAX
    integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    integer :: TINDEX(:)
    integer :: UINDEX(:)
    integer :: VINDEX(:)
    integer :: JINDEX(:)
    logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    character(len=132) :: STRING 
    character(len=3) :: DIRSTRING
    integer :: iString,LUFILE
    STRING(1:5) = '     '
    iSTRING = 6
    !step 2 determine where to put the 
    IF(iTUVP.LE.nTUVP.AND.((iTUVQ.GE.nTUVQSTART+1).AND.(iTUVQ.LE.nTUVQ)))THEN
       WRITE(STRING(iSTRING:iSTRING+6),'(A7)') 'ThetaP('
       iString = iSTRING+7
    ELSE
       IF(JTMP.LT.10)THEN
          WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMP,'('
          iString = iSTRING+5
       ELSE
          WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A)') 'Tmp',JTMP,'('
          iString = iSTRING+6
       ENDIF
    ENDIF
    IF(iTUVP.LT.100)THEN
       WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') iTUVP,','
       iString = iSTRING+3
    ELSEIF(iTUVP.LT.1000)THEN
       WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') iTUVP,','
       iString = iSTRING+4
    ELSEIF(iTUVP.LT.10000)THEN
       WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') iTUVP,','
       iString = iSTRING+5
    ELSE
       STOP 'Recurrent iTUVP'
    ENDIF
    IF(iTUVQ.LT.100)THEN
       WRITE(STRING(iSTRING:iSTRING+2),'(I2)') iTUVQ
       iString = iSTRING+2
    ELSEIF(iTUVQ.LT.1000)THEN
       WRITE(STRING(iSTRING:iSTRING+3),'(I3)') iTUVQ
       iString = iSTRING+3
    ELSEIF(iTUVQ.LT.10000)THEN
       WRITE(STRING(iSTRING:iSTRING+4),'(I4)') iTUVQ
       iString = iSTRING+4
    ELSE
       STOP 'Recurrent iTUVQ'
    ENDIF
    IF(iTUVP.LE.nTUVP.AND.((iTUVQ.GE.nTUVQSTART+1).AND.(iTUVQ.LE.nTUVQ)))THEN
       WRITE(STRING(iSTRING:iSTRING+12),'(A13)') ',IP) = '
       iString = iSTRING+13
    ELSE
       WRITE(STRING(iSTRING:iSTRING+3),'(A4)') ') = '
       iString = iSTRING+4
    ENDIF
    !step 3.  the first term: i+1
    
    IF(iTUVQminus1x.EQ.1)THEN
       !Aux(',iTUVp,',IP)
       WRITE(STRING(iSTRING:iSTRING+7),'(A8)') 'AuxCont('
       iString = iSTRING+8
    ELSE
       !          IF(iTUVplus1x.LE.nTUVP)THEN
       IF(iTUVplus1x.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+6),'(A7)') 'ThetaP('
          iString = iSTRING+7
       ELSE
          IF(JTMP-1.LT.10)THEN
             WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMP-1,'('
             iString = iSTRING+5
          ELSE
             WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMP-1,'('
             iString = iSTRING+6
          ENDIF
       ENDIF
    ENDIF
    IF(iTUVplus1x.LT.100)THEN
       WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVplus1x,','
       iString = iSTRING+3
    ELSEIF(iTUVplus1x.LT.1000)THEN
       WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVplus1x,','
       iString = iSTRING+4
    ELSEIF(iTUVplus1x.LT.10000)THEN
       WRITE(STRING(iSTRING:iSTRING+4),'(I4,A1)') iTUVplus1x,','
       iString = iSTRING+5
    ELSE
       STOP 'Recurrent iTUVplus1x'
    ENDIF
    IF(iTUVQminus1x.EQ.1)THEN
       !Aux(',iTUVp,',IP)
       WRITE(STRING(iSTRING:iSTRING+9),'(A10)') 'IP) '
       iString = iSTRING+10
    ELSE
       IF(iTUVplus1x.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+13),'(I2,A12)') iTUVQminus1x,',IP) '
          iString = iSTRING+14
       ELSE
          IF(iTUVQminus1x.LT.100)THEN
             WRITE(STRING(iSTRING:iSTRING+3),'(I2,A2)') iTUVQminus1x,') '
             iString = iSTRING+4
          ELSEIF(iTUVQminus1x.LT.1000)THEN
             WRITE(STRING(iSTRING:iSTRING+4),'(I3,A2)') iTUVQminus1x,') '
             iString = iSTRING+5
          ELSE
             STOP 'Recurrent iTUVPminus'
          ENDIF
       ENDIF
    ENDIF
    !step 4: the second term: X*Theta(i,j,k,l) (p,q-1)
    !either  Xab*AuxCont( 2,iTUVQ,IP) 
    !        Xab*Tmp1( 5, 2) 
    IF(iTUVQminus1x.EQ.1)THEN
       WRITE(STRING(iSTRING:iSTRING+13),'(A2,A3,A9)') '+ ',DIRSTRING,'*AuxCont('
       iString = iSTRING+14
    ELSE
       IF(iTUVP.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
!          WRITE(STRING(iSTRING:iSTRING+12),'(A2,A3,A8)') '+ ',DIRSTRING,'*ThetaP('
!          iString = iSTRING+13
          STOP 'ERROR '
       ELSE
          IF(JTMP-1.LT.10)THEN
             WRITE(STRING(iSTRING:iSTRING+10),'(A2,A3,A4,I1,A1)') '+ ',DIRSTRING,'*Tmp',JTMP-1,'('
             iString = iSTRING+11
          ELSE
             WRITE(STRING(iSTRING:iSTRING+11),'(A2,A3,A4,I2,A1)') '+ ',DIRSTRING,'*Tmp',JTMP-1,'('
             iString = iSTRING+12
          ENDIF
       ENDIF
    ENDIF
    IF(iTUVP.LT.100)THEN
       WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') iTUVP,','
       iString = iSTRING+3
    ELSEIF(iTUVP.LT.1000)THEN
       WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') iTUVP,','
       iString = iSTRING+4
    ELSEIF(iTUVP.LT.10000)THEN
       WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') iTUVP,','
       iString = iSTRING+5
    ELSE
       STOP 'Recurrent B iTUVP'
    ENDIF
    IF(iTUVQminus1x.EQ.1)THEN
       !Aux(',iTUVp,',IP)
       WRITE(STRING(iSTRING:iSTRING+9),'(A10)') 'IP) '
       iString = iSTRING+10
    ELSE
       IF(iTUVP.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+12),'(I2,A11)') iTUVQminus1x,',IP) '
          iString = iSTRING+13
       ELSE
          IF(iTUVQminus1x.LT.100)THEN
             WRITE(STRING(iSTRING:iSTRING+3),'(I2,A2)') iTUVQminus1x,') '
             iString = iSTRING+4
          ELSEIF(iTUVQminus1x.LT.1000)THEN
             WRITE(STRING(iSTRING:iSTRING+4),'(I3,A2)') iTUVQminus1x,') '
             iString = iSTRING+5
          ELSE
             STOP 'Recurrent iTUVPminus'
          ENDIF
       ENDIF
    ENDIF
    !Final step write the string
    WRITE(LUFILE,'(A)') STRING(1:iSTRING-1)
    
  end subroutine XYZRECURRENCE

subroutine URECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ,lufile)
    implicit none
    integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVQminus1x,iTUVQminus1y
    integer :: iTUVQminus1z,nTUVP,JQ,TMQ,TMP,JP,nTUVASTART,nTUVQSTART,nTUVQ
    integer :: iTUVQminus2x, iTUVP,nTUVTMPP,iTUVminus1x,iTUVplus1x,JMAX
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: TINDEX(:)
integer :: UINDEX(:)
integer :: VINDEX(:)
integer :: JINDEX(:)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
character(len=132) :: STRING 
integer :: iString,lufile
iTUVQ = TUVINDEX(Tq,Uq,Vq)
iTUVQminus1x = TUVINDEX(Tq,Uq-1,Vq)
iTUVQminus2x = TUVINDEX(Tq,Uq-2,Vq)
TMQ = Uq-1
do iTUVP = nTUVASTART+1,nTUVTMPP
   Tp = Tindex(iTUVp) 
   Up = Uindex(iTUVp) 
   Vp = Vindex(iTUVp)
   TMP = Up
   iTUVplus1x = TUVindex(Tp,Up+1,Vp)
   iTUVminus1x = TUVindex(Tp,Up-1,Vp)
   Jp = Jindex(iTUVp)   

   call XYZRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,&
        & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,nTUVQSTART,nTUVQ,&
        & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Yab',lufile)
ENDDO

end subroutine URECURRENCE

subroutine VRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ,lufile)
    implicit none
    integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVQminus1x,iTUVQminus1y
    integer :: iTUVQminus1z,nTUVP,JQ,TMQ,TMP,JP,nTUVASTART,nTUVQSTART,nTUVQ
    integer :: iTUVQminus2x, iTUVP,nTUVTMPP,iTUVminus1x,iTUVplus1x,JMAX
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: TINDEX(:)
integer :: UINDEX(:)
integer :: VINDEX(:)
integer :: JINDEX(:)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
character(len=132) :: STRING 
integer :: iString,lufile
iTUVQ = TUVINDEX(Tq,Uq,Vq)
iTUVQminus1x = TUVINDEX(Tq,Uq,Vq-1)
iTUVQminus2x = TUVINDEX(Tq,Uq,Vq-2)
TMQ = Vq-1
do iTUVP = nTUVASTART+1,nTUVTMPP
   Tp = Tindex(iTUVp) 
   Up = Uindex(iTUVp) 
   Vp = Vindex(iTUVp)
   TMP = Vp
   iTUVplus1x = TUVindex(Tp,Up,Vp+1)
   iTUVminus1x = TUVindex(Tp,Up,Vp-1)
   Jp = Jindex(iTUVp)   

   call XYZRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,&
        & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,nTUVQSTART,nTUVQ,&
        & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Zab',lufile)
ENDDO

end subroutine VRECURRENCE


end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
