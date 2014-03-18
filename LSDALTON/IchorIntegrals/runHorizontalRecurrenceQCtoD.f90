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
    integer :: GPUrun,MaxAngmomSingle
    logical :: DoOpenMP,DoOpenACC,CPU
    WRITE(*,'(A)')'MODULE AGC_OBS_HorizontalRecurrenceRHSModCtoD'
    WRITE(*,'(A)')' use IchorPrecisionModule'
    WRITE(*,'(A)')'  '
    WRITE(*,'(A)')' CONTAINS'
    MaxAngmomP = 4      ! currently D highest possibel XXDD
    MaxAngmomSingle = 2 ! currently D
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
    IF(GPUrun.EQ.2)WRITE(*,'(A)')'#ifdef VAR_OPENACC'
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
          IF(JP.EQ.0)THEN
!!$             WRITE(*,'(A)')' '
!!$             WRITE(*,'(A)')'!Unnecesarry as this is a simpel copy'
!!$             WRITE(*,'(A)')'!Transfer angmom from C to D'
!!$             WRITE(*,'(A)')'!subroutine HorizontalRR_RHS_Q0C0D0CtoD(nContPQ,nPasses,nlmP,&'
!!$             WRITE(*,'(A)')'!         & Qdistance12,ThetaP2,ThetaP,lupri)'
!!$             WRITE(*,'(A)')'!  implicit none'
!!$             WRITE(*,'(A)')'!  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri'
!!$             WRITE(*,'(A)')'!  real(realk),intent(in) :: Qdistance12(3)'
!!$             WRITE(*,'(A)')'!  real(realk),intent(in) :: ThetaP2(nlmP,1,nContPQ*nPasses)'
!!$             WRITE(*,'(A)')'!  real(realk),intent(inout) :: ThetaP(nlmP, 1,1,nContPQ*nPasses)'
!!$             WRITE(*,'(A)')'!  !Local variables'
!!$             WRITE(*,'(A)')'!  integer :: iP,ilmP'
!!$             WRITE(*,'(A)')'!  DO iP = 1,nPasses*nContPQ'
!!$             WRITE(*,'(A)')'!     DO ilmP = 1,nlmP'
!!$             WRITE(*,'(A)')'!        ThetaP(ilmP,1,1,IP) = ThetaP2(ilmP,1,IP)'
!!$             WRITE(*,'(A)')'!     ENDDO'
!!$             WRITE(*,'(A)')'!  ENDDO'
!!$             WRITE(*,'(A)')'!end subroutine HorizontalRR_RHS_Q0C0D0CtoD'
             CYCLE
          ENDIF
          WRITE(*,'(A)')''
          WRITE(*,'(A)')'!Transfer angmom from C to D'
          IF(JP.LT.10)THEN
             WRITE(*,'(A,I1,A,I1,A,I1,A)')'subroutine HorizontalRR_'//ARCSTRING//'_RHS_Q',JP,'C',AngmomA,'D',AngmomB,'CtoD(nContPQ,nPasses,nlmP,&'
          ELSE
             WRITE(*,'(A,I2,A,I1,A,I1,A)')'subroutine HorizontalRR_'//ARCSTRING//'_RHS_Q',JP,'C',AngmomA,'D',AngmomB,'CtoD(nContPQ,nPasses,nlmP,&'
          ENDIF
          WRITE(*,'(A)')'         & Qdistance12,ThetaP2,ThetaP,lupri)'
          WRITE(*,'(A)')'  implicit none'
          WRITE(*,'(A)')'  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri'
          WRITE(*,'(A)')'  real(realk),intent(in) :: Qdistance12(3)'
          WRITE(*,'(A,I5,A)')'  real(realk),intent(in) :: ThetaP2(nlmP,',nTUVP,',nContPQ*nPasses)'
          WRITE(*,'(A,I5,A1,I5,A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(nlmP,',NTUVAstart+1,':',nTUVA,',',nTUVBstart+1,':',nTUVB,',nContPQ*nPasses)'
          WRITE(*,'(A)')'  !Local variables'
          IF(ANGMOMB.NE.0)THEN
             WRITE(*,'(A)')'  integer :: iP,iC,iPassQ,ilmP,iTUVC'
             WRITE(*,'(A)')'  real(realk) :: Xcd,Ycd,Zcd'
          ELSE
             WRITE(*,'(A)')'  integer :: iP,ilmP,iTUVC'
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
                WRITE(*,'(A,I1,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVAstart+1,':',nTUVTMPA,',',nTUVTMPBprev+1,':',nTUVTMPB,')'
             else
                WRITE(*,'(A,I2,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVAstart+1,':',nTUVTMPA,',',nTUVTMPBprev+1,':',nTUVTMPB,')'
             endif
          ENDDO
          WRITE(*,'(A)')'!  real(realk) :: Tmp(nTUVA,nTUVB) ordering'
          IF(DoOpenMP)THEN
             WRITE(*,'(A)')'!$OMP PARALLEL DO DEFAULT(none) &'
             IF(JB.NE.0)THEN
                WRITE(*,'(A)')'!$OMP PRIVATE(iP,&'
                DO JTMP=1,JB-1
                   if(JTMP.LT.10)THEN
                      WRITE(*,'(A,I1,A)')'!$OMP         Tmp',JTMP,',&'
                   else
                      WRITE(*,'(A,I2,A)')'!$OMP         Tmp',JTMP,',&'
                   endif
                ENDDO
                WRITE(*,'(A)')'!$OMP         iTUVC,ilmP,Xcd,Ycd,Zcd) &'
                WRITE(*,'(A)')'!$OMP SHARED(nlmP,nContPQ,nPasses,Qdistance12,ThetaP,ThetaP2)'
             ELSE
                WRITE(*,'(A)')'!$OMP PRIVATE(iP,iTUVC,ilmP) &'
                WRITE(*,'(A)')'!$OMP SHARED(nlmP,nContPQ,nPasses,ThetaP,ThetaP2)'
             ENDIF
          ENDIF
          IF(DoOpenACC)THEN
             WRITE(*,'(A)')'!$ACC PARALLEL LOOP &'
             IF(JB.NE.0)THEN
                WRITE(*,'(A)')'!$ACC PRIVATE(iP,&'
                DO JTMP=1,JB-1
                   if(JTMP.LT.10)THEN
                      WRITE(*,'(A,I1,A)')'!$ACC         Tmp',JTMP,',&'
                   else
                      WRITE(*,'(A,I2,A)')'!$ACC         Tmp',JTMP,',&'
                   endif
                ENDDO
                WRITE(*,'(A)')'!$ACC         iTUVC,ilmP,Xcd,Ycd,Zcd) &'
                WRITE(*,'(A)')'!$ACC PRESENT(nlmP,nContPQ,nPasses,Qdistance12,ThetaP,ThetaP2)'
             ELSE
                WRITE(*,'(A)')'!$ACC PRIVATE(iP,iTUVC,ilmP) &'
                WRITE(*,'(A)')'!$ACC PRESENT(nlmP,nContPQ,nPasses,ThetaP,ThetaP2)'
             ENDIF
          ENDIF
          WRITE(*,'(A)')'  DO iP = 1,nContPQ*nPasses'
          IF(JB.NE.0)THEN
             WRITE(*,'(A)')'   Xcd = Qdistance12(1)'
             WRITE(*,'(A)')'   Ycd = Qdistance12(2)'
             WRITE(*,'(A)')'   Zcd = Qdistance12(3)'
          ENDIF
          DO JTMP=0,JB
             !           print*,'!JTMP = ',JTMP
             IF(JTMP.EQ.0)THEN
                nTUVTMPP=(Jab-JTMP)*(Jab-JTMP+1)*(Jab-JTMP+2)/6
                IF(nTUVBstart+1.EQ.1)THEN
                   WRITE(*,'(A,I3,A,I3)')'    DO iTUVC=',NTUVAstart+1,',',nTUVA
                   WRITE(*,'(A)')'     DO ilmP = 1,nlmP'
                   WRITE(*,'(A)')   '        ThetaP(ilmP,iTUVC,1,IP) = ThetaP2(ilmP,iTUVC,IP)'
                   WRITE(*,'(A)')   '     ENDDO'
                   WRITE(*,'(A)')   '    ENDDO'
                   IF(JB.GT.0)THEN
                      WRITE(*,'(A)')   '    DO ilmP = 1,nlmP'
                   ENDIF
                ELSE
                   IF(JB.GT.0)THEN
                      WRITE(*,'(A)')   '    DO ilmP = 1,nlmP'
                   ENDIF
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
                         call TRECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB)
                      ELSEIF(UREC)THEN
                         !                         print*,'!A URECURRENCE iTUV',iTUV
                         call URECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB)
                      ELSEIF(VREC)THEN
                         !                         print*,'!A VRECURRENCE iTUV',iTUV
                         call VRECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB)
                      ELSE
                         STOP 'TK1'
                      ENDIF
                      CREATED(Tb,Ub,Vb) = .TRUE.
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
          deallocate(CREATED)
          IF(JB.GT.0)THEN
             WRITE(*,'(A)')'    ENDDO'
          ENDIF
          WRITE(*,'(A)')'  ENDDO'
          IF(DoOpenMP)WRITE(*,'(A)')'!$OMP END PARALLEL DO'
          IF(JP.LT.10)THEN
             WRITE(*,'(A,I1,A,I1,A,I1,A)')'end subroutine HorizontalRR_'//ARCSTRING//'_RHS_Q',JP,'C',AngmomA,'D',AngmomB,'CtoD'
          ELSE
             WRITE(*,'(A,I2,A,I1,A,I1,A)')'end subroutine HorizontalRR_'//ARCSTRING//'_RHS_Q',JP,'C',AngmomA,'D',AngmomB,'CtoD'
          ENDIF
       enddo
       deallocate(TUVINDEX)
       deallocate(TINDEX)
       deallocate(UINDEX)
       deallocate(VINDEX)
       deallocate(JINDEX)
    enddo
    IF(GPUrun.EQ.2)WRITE(*,'(A)')'#endif'
 enddo
    WRITE(*,'(A)')'end module'
  END subroutine PASSsub
  
  subroutine TRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ)
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
    integer :: iString
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
            & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
            & nTUVQSTART,nTUVQ,iTUVQ,iTUVQminus1x,TMQ,iTUVP,&
            & Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Xcd')
    ENDDO

  end subroutine TRECURRENCE

  subroutine XYZRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,&
       & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ,iTUVQ,iTUVQminus1x,TMQ,iTUVP,&
       & Tp,Up,Vp,TMP,iTUVplus1x,Jp,DIRSTRING)
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
    character(len=3) :: DIRSTRING
    integer :: iString
    !step 1 add blanks
    STRING(1:5) = '     '
    iSTRING = 6
    !step 2 determine where to put the 
    IF(iTUVP.LE.nTUVP.AND.((iTUVQ.GE.nTUVQSTART+1).AND.(iTUVQ.LE.nTUVQ)))THEN
       WRITE(STRING(iSTRING:iSTRING+11),'(A12)') 'ThetaP(ilmP,'
       iString = iSTRING+12
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
       WRITE(STRING(iSTRING:iSTRING+6),'(A7)') ',IP) = '
       iString = iSTRING+7
    ELSE
       WRITE(STRING(iSTRING:iSTRING+3),'(A4)') ') = '
       iString = iSTRING+4
    ENDIF
    !step 3.  the first term: i+1 (p+1,q-1)
    
    IF(iTUVQminus1x.EQ.1)THEN
       !Aux(',iTUVp,',IP)
       WRITE(STRING(iSTRING:iSTRING+12),'(A13)') 'ThetaP2(ilmP,'
       iString = iSTRING+13
    ELSE
       !          IF(iTUVplus1x.LE.nTUVP)THEN
       IF(iTUVplus1x.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+11),'(A12)') 'ThetaP(ilmP,'
          iString = iSTRING+12
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
       WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'IP) '
       iString = iSTRING+4
    ELSE
       IF(iTUVplus1x.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+6),'(I2,A5)') iTUVQminus1x,',IP) '
          iString = iSTRING+7
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
    !step 4: the second term: X*Theta(i,j,k,l)  (p,q-1)
    IF(iTUVQminus1x.EQ.1)THEN
       WRITE(STRING(iSTRING:iSTRING+18),'(A2,A3,A14)') '+ ',DIRSTRING,'*ThetaP2(ilmP,'
       iString = iSTRING+19
    ELSE
       IF(iTUVP.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+17),'(A2,A3,A13)') '+ ',DIRSTRING,'*ThetaP(ilmP,'
          iString = iSTRING+18
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
       WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'IP) '
       iString = iSTRING+4
    ELSE
       IF(iTUVP.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+6),'(I2,A5)') iTUVQminus1x,',IP) '
          iString = iSTRING+7
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
    WRITE(*,'(A)') STRING(1:iSTRING-1)
  end subroutine XYZRECURRENCE

  subroutine URECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ)
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
    integer :: iString
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
            & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
            & nTUVQSTART,nTUVQ,iTUVQ,iTUVQminus1x,TMQ,iTUVP,&
            & Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Ycd')
    ENDDO

  end subroutine URECURRENCE

  subroutine VRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ)
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
    integer :: iString
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
            & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
            & nTUVQSTART,nTUVQ,iTUVQ,iTUVQminus1x,TMQ,iTUVP,&
            & Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Zcd')
    ENDDO

  end subroutine VRECURRENCE


end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
