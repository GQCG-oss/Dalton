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
    Integer :: nTUVTMPB,nTUVTMPA,nTUVTMPAprev,Ta,Ua,Va,JA
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
    logical :: CPU,DoOpenMP,DoOpenACC
    integer :: GPUrun,MaxAngmomSingle,I,LUFILE
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
       WRITE(FileName,'(4A)')'runHorizontalRecurrence'//ARCSTRING//'LHSModBtoA.F90'
       print*,'FileName:',FileName
       LUFILE = 12
       open(unit = LUFILE, file=TRIM(FileName),status="unknown")
       
       WRITE(LUFILE,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceLHSModBtoA'
       WRITE(LUFILE,'(A)')' use IchorPrecisionModule'
       WRITE(LUFILE,'(A)')'  '
       WRITE(LUFILE,'(A)')' CONTAINS'

!    IF(GPUrun.EQ.2)WRITE(LUFILE,'(A)')'#ifdef VAR_OPENACC'
    MaxAngmomP = 6
    MaxAngmomSingle = 2 !D functions currently


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
          IF(AngmomB.LE.AngmomA)CYCLE
          IF(AngmomA.GT.MaxAngmomSingle)CYCLE
          IF(AngmomB.GT.MaxAngmomSingle)CYCLE
          IF(AngmomA.EQ.0.AND.AngmomB.EQ.0)CYCLE
          NTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
          NTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6

          NTUVAspec = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6 - &
               &         (AngmomA)*(AngmomA+1)*(AngmomA+2)/6
          NTUVBspec = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6 - &
               &         (AngmomB)*(AngmomB+1)*(AngmomB+2)/6
          NTUVAstart = (AngmomA)*(AngmomA+1)*(AngmomA+2)/6
          NTUVBstart = (AngmomB)*(AngmomB+1)*(AngmomB+2)/6

          WRITE(LUFILE,'(A)')''
          IF(JP.LT.10)THEN
             WRITE(LUFILE,'(A,I1,A,I1,A,I1,A)')'subroutine HorizontalRR_'//ARCSTRING//'_LHS_P',JP,'A',AngmomA,'B',AngmomB,'BtoA(nContQP,nPasses,nTUVQ,&'
          ELSE
             WRITE(LUFILE,'(A,I2,A,I1,A,I1,A)')'subroutine HorizontalRR_'//ARCSTRING//'_LHS_P',JP,'A',AngmomA,'B',AngmomB,'BtoA(nContQP,nPasses,nTUVQ,&'
          ENDIF
          WRITE(LUFILE,'(A)')'         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)'
          WRITE(LUFILE,'(A)')'  implicit none'
          WRITE(LUFILE,'(A)')'  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB'
          WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)'
          WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: AuxCont(',nTUVP,',nTUVQ*nContQP*nPasses)'
          WRITE(LUFILE,'(A)')'  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
          !             WRITE(LUFILE,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(',nTUVAspec,',',nTUVBspec,',nTUVQ,nContPasses)'
          WRITE(LUFILE,'(A,I5,A1,I5,A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(',NTUVAstart+1,':',nTUVA,',',nTUVBstart+1,':',nTUVB,',nTUVQ*nContQP*nPasses)'
          WRITE(LUFILE,'(A)')'  !Local variables'
          IF(ANGMOMA.NE.0)THEN
             WRITE(LUFILE,'(A)')'  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB'
             WRITE(LUFILE,'(A)')'  real(realk) :: Xab,Yab,Zab'
          ELSE
             WRITE(LUFILE,'(A)')'  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB'
          ENDIF
          allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
          CREATED  = .FALSE.
          CREATED(0,0,0) = .TRUE.
          JAB = ANGMOMA + ANGMOMB
          JA = ANGMOMA
          JB = ANGMOMB
          DO JTMP=1,JA-1
             nTUVTMPAprev=(JTMP)*(JTMP+1)*(JTMP+2)/6
             nTUVTMPA=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
             nTUVTMPB=(Jab-JTMP+1)*(Jab-JTMP+2)*(Jab-JTMP+3)/6
             nTUVTMPBprev=(Jab-JTMP)*(Jab-JTMP+1)*(Jab-JTMP+2)/6
             if(JTMP.LT.10)THEN
                WRITE(LUFILE,'(A,I1,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVTMPAprev+1,':',nTUVTMPA,',',nTUVBstart+1,':',nTUVTMPB,')'
             else
                WRITE(LUFILE,'(A,I2,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVTMPAprev+1,':',nTUVTMPA,',',nTUVBstart+1,':',nTUVTMPB,')'
             endif
          ENDDO
          IF(JA.GT.1)THEN
             WRITE(LUFILE,'(A)')'!  real(realk) :: Tmp(nTUVA,nTUVB) ordering'
          ENDIF
!          WRITE(LUFILE,'(A)')'  DO iP = 1,nContQP*nPasses'
!          WRITE(LUFILE,'(A)')'   DO iTUVQ = 1,nTUVQ'
          IF(DoOpenMP)THEN
!             WRITE(LUFILE,'(A)')'!$OMP PARALLEL DO DEFAULT(none) &'
             WRITE(LUFILE,'(A)')'!$OMP DO &'
             WRITE(LUFILE,'(A)')'!$OMP PRIVATE(iP,&'
             IF(JA.NE.0)THEN
                DO JTMP=1,JA-1
                   if(JTMP.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A)')'!$OMP         Tmp',JTMP,',&'
                   else
                      WRITE(LUFILE,'(A,I2,A)')'!$OMP         Tmp',JTMP,',&'
                   endif
                ENDDO
                WRITE(LUFILE,'(A)')'!$OMP         iPassP,iTUVB,iAtomA,iAtomB,Xab,Yab,Zab) '
             ELSE
                WRITE(LUFILE,'(A)')'!$OMP         iTUVB) '
             ENDIF
!             WRITE(LUFILE,'(A)')'!$OMP SHARED(nTUVQ,nContQP,nPasses,&'
!             IF(JA.NE.0)THEN
!                WRITE(LUFILE,'(A)')'!$OMP         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP)'
!             ELSE
!                WRITE(LUFILE,'(A)')'!$OMP         AuxCont,ThetaP)'
!             ENDIF
          ENDIF
          IF(DoOpenACC)THEN
             WRITE(LUFILE,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUFILE,'(A)')'!$ACC PRIVATE(iP,&'
             IF(JA.NE.0)THEN
                DO JTMP=1,JA-1
                   if(JTMP.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A)')'!$ACC         Tmp',JTMP,',&'
                   else
                      WRITE(LUFILE,'(A,I2,A)')'!$ACC         Tmp',JTMP,',&'
                   endif
                ENDDO
                WRITE(LUFILE,'(A)')'!$ACC         iPassP,iTUVB,iAtomA,iAtomB,Xab,Yab,Zab) &'
             ELSE
                WRITE(LUFILE,'(A)')'!$ACC         iTUVB) &'
             ENDIF
             WRITE(LUFILE,'(A)')'!$ACC PRESENT(nTUVQ,nContQP,nPasses,&'
             IF(JA.NE.0)THEN
                WRITE(LUFILE,'(A)')'!$ACC         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP)'
             ELSE
                WRITE(LUFILE,'(A)')'!$ACC         AuxCont,ThetaP)'
             ENDIF
          ENDIF
          WRITE(LUFILE,'(A)')'  DO iP = 1,nTUVQ*nContQP*nPasses'
          IF(JA.NE.0)THEN
             WRITE(LUFILE,'(A)')'   iPassP = (iP-1)/(nTUVQ*nContQP)+1'
             WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
             WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
             WRITE(LUFILE,'(A)')'   Xab = -Pdistance12(1,iAtomA,iAtomB)'
             WRITE(LUFILE,'(A)')'   Yab = -Pdistance12(2,iAtomA,iAtomB)'
             WRITE(LUFILE,'(A)')'   Zab = -Pdistance12(3,iAtomA,iAtomB)'
          ENDIF


          DO JTMP=0,JA
             !           print*,'!JTMP = ',JTMP
             IF(JTMP.EQ.0)THEN
                
                IF(nTUVAstart+1.EQ.1)THEN
                   WRITE(LUFILE,'(A,I3,A,I3)')'     DO iTUVB=',NTUVBstart+1,',',nTUVB
                   WRITE(LUFILE,'(A)')   '        ThetaP(1,iTUVB,iP) = AuxCont(iTUVB,iP)'
                   WRITE(LUFILE,'(A)')   '     ENDDO'
                ENDIF
                CYCLE
             ELSE
                nTUVTMPA=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
                nTUVTMPB=(JAB-JTMP+1)*(JAB-JTMP+2)*(JAB-JTMP+3)/6
                !slowly increase JB=JTMP with 1 
                nTUVTMPA=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
                DO Ta=JTMP,0,-1       
                   DO Ua=JTMP-Ta,0,-1
                      Va=JTMP-Ta-Ua                
                      !                   print*,'!Tq,Uq,Vq = ',Tq,Uq,Vq
                      !step 1.
                      !how can the (Tp,Up,Vp) be built from lower
                      TREC = CREATED(Ta-1,Ua,Va)
                      UREC = CREATED(Ta,Ua-1,Va)
                      VREC = CREATED(Ta,Ua,Va-1)
                      N=0
                      IF(TREC)N=N+1
                      IF(UREC)N=N+1
                      IF(VREC)N=N+1
                      IF(TREC)THEN
                         !                         print*,'!A TRECURRENCE iTUV',iTUV
                         call TRECURRENCE(Ta,Ua,Va,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPB,JMAX,JA,nTUVB,nTUVBSTART,nTUVAstart,nTUVA,LUFILE)
                      ELSEIF(UREC)THEN
                         !                         print*,'!A URECURRENCE iTUV',iTUV
                         call URECURRENCE(Ta,Ua,Va,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPB,JMAX,JA,nTUVB,nTUVBSTART,nTUVAstart,nTUVA,LUFILE)
                      ELSEIF(VREC)THEN
                         !                         print*,'!A VRECURRENCE iTUV',iTUV
                         call VRECURRENCE(Ta,Ua,Va,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPB,JMAX,JA,nTUVB,nTUVBSTART,nTUVAstart,nTUVA,LUFILE)
                      ELSE
                         STOP 'TK1'
                      ENDIF
                      CREATED(Ta,Ua,Va) = .TRUE.
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
          deallocate(CREATED)
!          WRITE(LUFILE,'(A)')'   ENDDO'
          WRITE(LUFILE,'(A)')'  ENDDO'
!          IF(DoOpenMP) WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'
          IF(DoOpenMP) WRITE(LUFILE,'(A)')'!$OMP END DO'
          IF(JP.LT.10)THEN
             WRITE(LUFILE,'(A,I1,A,I1,A,I1,A)')'end subroutine HorizontalRR_'//ARCSTRING//'_LHS_P',JP,'A',AngmomA,'B',AngmomB,'BtoA'
          ELSE
             WRITE(LUFILE,'(A,I2,A,I1,A,I1,A)')'end subroutine HorizontalRR_'//ARCSTRING//'_LHS_P',JP,'A',AngmomA,'B',AngmomB,'BtoA'
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
  
  subroutine TRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVPSTART,&
       & nTUVQSTART,nTUVQ,LUFILE)
    implicit none
    integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVQminus1x,iTUVQminus1y
    integer :: iTUVQminus1z,nTUVP,JQ,TMQ,TMP,JP,nTUVPSTART,nTUVQSTART,nTUVQ
    integer :: iTUVQminus2x, iTUVP,nTUVTMPP,iTUVminus1x,iTUVplus1x,JMAX
    integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    integer :: TINDEX(:)
    integer :: UINDEX(:)
    integer :: VINDEX(:)
    integer :: JINDEX(:)
    logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
    character(len=132) :: STRING 
    integer :: iString,LUFILE
    !Ta = Tq  iTUVA = iTUVC
    iTUVQ = TUVINDEX(Tq,Uq,Vq)
    iTUVQminus1x = TUVINDEX(Tq-1,Uq,Vq)
    TMQ = Tq-1
    do iTUVP = nTUVPSTART+1,nTUVTMPP
       Tp = Tindex(iTUVp) 
       Up = Uindex(iTUVp) 
       Vp = Vindex(iTUVp)
       TMP = Tp
       iTUVplus1x = TUVindex(Tp+1,Up,Vp)
       Jp = Jindex(iTUVp)   
       
       call XYZRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,&
            & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVPSTART,nTUVQSTART,nTUVQ,&
            & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Xab',LUFILE)

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
    !Q = A so A first A have the 
    !iTUVQ
    !iTUVQminus1x
    !P = B (the one with all the angmom)
    !iTUVP
    !iTUVPplus1x
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
    IF(iTUVQ.LT.100)THEN
       WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') iTUVQ,','
       iString = iSTRING+3
    ELSEIF(iTUVQ.LT.1000)THEN
       WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') iTUVQ,','
       iString = iSTRING+4
    ELSEIF(iTUVQ.LT.10000)THEN
       WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') iTUVQ,','
       iString = iSTRING+5
    ELSE
       STOP 'Recurrent iTUVQ'
    ENDIF
    IF(iTUVP.LT.100)THEN
       WRITE(STRING(iSTRING:iSTRING+2),'(I2)') iTUVP
       iString = iSTRING+2
    ELSEIF(iTUVP.LT.1000)THEN
       WRITE(STRING(iSTRING:iSTRING+3),'(I3)') iTUVP
       iString = iSTRING+3
    ELSEIF(iTUVP.LT.10000)THEN
       WRITE(STRING(iSTRING:iSTRING+4),'(I4)') iTUVP
       iString = iSTRING+4
    ELSE
       STOP 'Recurrent iTUVP'
    ENDIF
    IF(iTUVP.LE.nTUVP.AND.((iTUVQ.GE.nTUVQSTART+1).AND.(iTUVQ.LE.nTUVQ)))THEN
       WRITE(STRING(iSTRING:iSTRING+12),'(A13)') ',iP) = '
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
    IF(iTUVQminus1x.NE.1)THEN
       IF(iTUVQminus1x.LT.100)THEN
          WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVQminus1x,','
          iString = iSTRING+3
       ELSEIF(iTUVQminus1x.LT.1000)THEN
          WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVQminus1x,','
          iString = iSTRING+4
       ELSEIF(iTUVQminus1x.LT.10000)THEN
          WRITE(STRING(iSTRING:iSTRING+4),'(I4,A1)') iTUVQminus1x,','
          iString = iSTRING+5
       ELSE
          STOP 'Recurrent iTUVQplus1x'
       ENDIF
    ENDIF
    IF(iTUVQminus1x.EQ.1)THEN
       IF(iTUVplus1x.LT.100)THEN
          WRITE(STRING(iSTRING:iSTRING+11),'(I2,A10)') iTUVplus1x,',iP)'
          iString = iSTRING+12
       ELSEIF(iTUVplus1x.LT.1000)THEN
          WRITE(STRING(iSTRING:iSTRING+12),'(I3,A10)') iTUVplus1x,',iP)'
          iString = iSTRING+13
       ELSEIF(iTUVplus1x.LT.10000)THEN
          WRITE(STRING(iSTRING:iSTRING+13),'(I4,A10)') iTUVplus1x,',iP)'
          iString = iSTRING+14
       ELSE
          STOP 'Recurrent iTUVplus1x'
       ENDIF
    ELSE
       IF(iTUVplus1x.LT.100)THEN
          WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVplus1x,')'
          iString = iSTRING+3
       ELSEIF(iTUVplus1x.LT.1000)THEN
          WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVplus1x,')'
          iString = iSTRING+4
       ELSEIF(iTUVplus1x.LT.10000)THEN
          WRITE(STRING(iSTRING:iSTRING+4),'(I4,A1)') iTUVplus1x,')'
          iString = iSTRING+5
       ELSE
          STOP 'Recurrent iTUVplus1x'
       ENDIF
    ENDIF
    !step 4: the second term: X*Theta(i,j,k,l) (q-1,p)
    IF(iTUVQminus1x.EQ.1)THEN
       WRITE(STRING(iSTRING:iSTRING+13),'(A2,A3,A9)') '+ ',DIRSTRING,'*AuxCont('
       iString = iSTRING+14

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
       WRITE(STRING(iSTRING:iSTRING+9),'(A10)') 'ip) '
       iString = iSTRING+10

    ELSE
       IF(JTMP-1.LT.10)THEN
          WRITE(STRING(iSTRING:iSTRING+10),'(A2,A3,A4,I1,A1)') '+ ',DIRSTRING,'*Tmp',JTMP-1,'('
          iString = iSTRING+11
       ELSE
          WRITE(STRING(iSTRING:iSTRING+11),'(A2,A3,A4,I2,A1)') '+ ',DIRSTRING,'*Tmp',JTMP-1,'('
          iString = iSTRING+12
       ENDIF
       IF(iTUVQminus1x.LT.100)THEN
          WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVQminus1x,','
          iString = iSTRING+3
       ELSEIF(iTUVQminus1x.LT.1000)THEN
          WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVQminus1x,','
          iString = iSTRING+4
       ELSE
          STOP 'Recurrent iTUVPminus'
       ENDIF
       IF(iTUVP.LT.100)THEN
          WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') iTUVP,')'
          iString = iSTRING+3
       ELSEIF(iTUVP.LT.1000)THEN
          WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') iTUVP,')'
          iString = iSTRING+4
       ELSEIF(iTUVP.LT.10000)THEN
          WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') iTUVP,')'
          iString = iSTRING+5
       ELSE
          STOP 'Recurrent B iTUVP'
       ENDIF
    ENDIF
    !Final step write the string
    WRITE(LUFILE,'(A)') STRING(1:iSTRING-1)
    
  end subroutine XYZRECURRENCE

subroutine URECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ,LUFILE)
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
integer :: iString,LUFILE
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
            & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Yab',LUFILE)

ENDDO

end subroutine URECURRENCE

subroutine VRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ,LUFILE)
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
integer :: iString,LUFILE
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
            & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Zab',LUFILE)

ENDDO

end subroutine VRECURRENCE


end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
