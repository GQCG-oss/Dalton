MODULE TESTMODULE
  use stringsMODULE
 logical,save :: nPrimLast
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
    integer :: GPUrun,MaxAngmomSingle,LUFILE,I
    logical :: DoOpenMP,DoOpenACC,CPU   
    Character(len=48) :: FileName    
    MaxAngmomP = 4      ! currently D highest possibel XXDD
    MaxAngmomSingle = 2 ! currently D
DO GPUrun = 1,2
    CPU = .TRUE.
    IF(GPUrun.EQ.2)CPU = .FALSE.
    nPrimLAST = .FALSE.
    IF(CPU)nPrimLAST = .TRUE.
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
    WRITE(FileName,'(4A)')'runHorizontalRecurrence'//ARCSTRING//'RHSModCtoD.F90'

    print*,'FileName:',FileName
    LUFILE = 12
    open(unit = LUFILE, file=TRIM(FileName),status="unknown")

    WRITE(LUFILE,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceRHSModCtoD'
    WRITE(LUFILE,'(A)')' use IchorPrecisionModule'
    WRITE(LUFILE,'(A)')'  '
    WRITE(LUFILE,'(A)')' CONTAINS'


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
          IF(JP.EQ.0)THEN
!!$             WRITE(LUFILE,'(A)')' '
!!$             WRITE(LUFILE,'(A)')'!Unnecesarry as this is a simpel copy'
!!$             WRITE(LUFILE,'(A)')'!Transfer angmom from C to D'
!!$             WRITE(LUFILE,'(A)')'!subroutine HorizontalRR_RHS_Q0C0D0CtoD(nContPQ,nPasses,nlmP,&'
!!$             WRITE(LUFILE,'(A)')'!         & Qdistance12,ThetaP2,ThetaP,lupri)'
!!$             WRITE(LUFILE,'(A)')'!  implicit none'
!!$             WRITE(LUFILE,'(A)')'!  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri'
!!$             WRITE(LUFILE,'(A)')'!  real(realk),intent(in) :: Qdistance12(3)'
!!$             WRITE(LUFILE,'(A)')'!  real(realk),intent(in) :: ThetaP2(nlmP,1,nContPQ*nPasses)'
!!$             WRITE(LUFILE,'(A)')'!  real(realk),intent(inout) :: ThetaP(nlmP, 1,1,nContPQ*nPasses)'
!!$             WRITE(LUFILE,'(A)')'!  !Local variables'
!!$             WRITE(LUFILE,'(A)')'!  integer :: iP,ilmP'
!!$             WRITE(LUFILE,'(A)')'!  DO iP = 1,nPasses*nContPQ'
!!$             WRITE(LUFILE,'(A)')'!     DO ilmP = 1,nlmP'
!!$             WRITE(LUFILE,'(A)')'!        ThetaP(ilmP,1,1,IP) = ThetaP2(ilmP,1,IP)'
!!$             WRITE(LUFILE,'(A)')'!     ENDDO'
!!$             WRITE(LUFILE,'(A)')'!  ENDDO'
!!$             WRITE(LUFILE,'(A)')'!end subroutine HorizontalRR_RHS_Q0C0D0CtoD'
             CYCLE
          ENDIF
          WRITE(LUFILE,'(A)')''
          WRITE(LUFILE,'(A)')'!Transfer angmom from C to D'
          IF(JP.LT.10)THEN
             WRITE(LUFILE,'(A,I1,A,I1,A,I1,A)')'subroutine HorizontalRR_'//ARCSTRING//'_RHS_Q',JP,'C',AngmomA,'D',AngmomB,'CtoD(nContPQ,nPasses,nlmP,&'
          ELSE
             WRITE(LUFILE,'(A,I2,A,I1,A,I1,A)')'subroutine HorizontalRR_'//ARCSTRING//'_RHS_Q',JP,'C',AngmomA,'D',AngmomB,'CtoD(nContPQ,nPasses,nlmP,&'
          ENDIF
          WRITE(LUFILE,'(A)')'         & Qdistance12,ThetaP2,ThetaP,lupri)'
          WRITE(LUFILE,'(A)')'  implicit none'
          WRITE(LUFILE,'(A)')'  integer,intent(in) :: nContPQ,nPasses,nlmP,lupri'
          WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Qdistance12(3)'
          IF(nPrimLAST)THEN
             WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: ThetaP2(nlmP,',nTUVP,',nContPQ*nPasses)'
             WRITE(LUFILE,'(A,I5,A1,I5,A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(nlmP,',NTUVAstart+1,':',nTUVA,',',nTUVBstart+1,':',nTUVB,',nContPQ*nPasses)'
          ELSE
             WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(in) :: ThetaP2(nContPQ*nPasses,nlmP,',nTUVP,')'
             WRITE(LUFILE,'(A,I5,A1,I5,A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(nContPQ*nPasses,nlmP,',NTUVAstart+1,':',nTUVA,',',nTUVBstart+1,':',nTUVB,')'
          ENDIF
          WRITE(LUFILE,'(A)')'  !Local variables'
          IF(ANGMOMB.NE.0)THEN
             WRITE(LUFILE,'(A)')'  integer :: iP,iC,iPassQ,ilmP,iTUVC'
             WRITE(LUFILE,'(A)')'  real(realk) :: Xcd,Ycd,Zcd'
          ELSE
             WRITE(LUFILE,'(A)')'  integer :: iP,ilmP,iTUVC'
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
          WRITE(LUFILE,'(A)')'!  real(realk) :: Tmp(nTUVA,nTUVB) ordering'
          IF(DoOpenMP)THEN
!             WRITE(LUFILE,'(A)')'!$OMP PARALLEL DO DEFAULT(none) &'
             WRITE(LUFILE,'(A)')'!$OMP DO &'
             IF(JB.NE.0)THEN
                WRITE(LUFILE,'(A)')'!$OMP PRIVATE(iP,&'
                DO JTMP=1,JB-1
                   if(JTMP.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A)')'!$OMP         Tmp',JTMP,',&'
                   else
                      WRITE(LUFILE,'(A,I2,A)')'!$OMP         Tmp',JTMP,',&'
                   endif
                ENDDO
                WRITE(LUFILE,'(A)')'!$OMP         iTUVC,ilmP,Xcd,Ycd,Zcd) '
!                WRITE(LUFILE,'(A)')'!$OMP SHARED(nlmP,nContPQ,nPasses,Qdistance12,ThetaP,ThetaP2)'
             ELSE
                WRITE(LUFILE,'(A)')'!$OMP PRIVATE(iP,iTUVC,ilmP) '
!                WRITE(LUFILE,'(A)')'!$OMP SHARED(nlmP,nContPQ,nPasses,ThetaP,ThetaP2)'
             ENDIF
          ENDIF
          IF(DoOpenACC)THEN
             WRITE(LUFILE,'(A)')'!$ACC PARALLEL LOOP &'
             IF(JB.NE.0)THEN
                WRITE(LUFILE,'(A)')'!$ACC PRIVATE(iP,&'
                DO JTMP=1,JB-1
                   if(JTMP.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A)')'!$ACC         Tmp',JTMP,',&'
                   else
                      WRITE(LUFILE,'(A,I2,A)')'!$ACC         Tmp',JTMP,',&'
                   endif
                ENDDO
                WRITE(LUFILE,'(A)')'!$ACC         iTUVC,ilmP,Xcd,Ycd,Zcd) &'
                WRITE(LUFILE,'(A)')'!$ACC PRESENT(nContPQ,nPasses,Qdistance12,ThetaP,ThetaP2)'
             ELSE
                WRITE(LUFILE,'(A)')'!$ACC PRIVATE(iP,iTUVC,ilmP) &'
                WRITE(LUFILE,'(A)')'!$ACC PRESENT(nContPQ,nPasses,ThetaP,ThetaP2)'
             ENDIF
          ENDIF
          WRITE(LUFILE,'(A)')'  DO iP = 1,nContPQ*nPasses'
          IF(JB.NE.0)THEN
             WRITE(LUFILE,'(A)')'   Xcd = Qdistance12(1)'
             WRITE(LUFILE,'(A)')'   Ycd = Qdistance12(2)'
             WRITE(LUFILE,'(A)')'   Zcd = Qdistance12(3)'
          ENDIF
          DO JTMP=0,JB
             !           print*,'!JTMP = ',JTMP
             IF(JTMP.EQ.0)THEN
                nTUVTMPP=(Jab-JTMP)*(Jab-JTMP+1)*(Jab-JTMP+2)/6
                IF(nTUVBstart+1.EQ.1)THEN
                   WRITE(LUFILE,'(A,I3,A,I3)')'    DO iTUVC=',NTUVAstart+1,',',nTUVA
                   WRITE(LUFILE,'(A)')'     DO ilmP = 1,nlmP'
                   IF(nPrimLAST)THEN
                      WRITE(LUFILE,'(A)')   '        ThetaP(ilmP,iTUVC,1,IP) = ThetaP2(ilmP,iTUVC,IP)'
                   ELSE
                      WRITE(LUFILE,'(A)')   '        ThetaP(IP,ilmP,iTUVC,1) = ThetaP2(IP,ilmP,iTUVC)'
                   ENDIF
                   WRITE(LUFILE,'(A)')   '     ENDDO'
                   WRITE(LUFILE,'(A)')   '    ENDDO'
                   IF(JB.GT.0)THEN
                      WRITE(LUFILE,'(A)')   '    DO ilmP = 1,nlmP'
                   ENDIF
                ELSE
                   IF(JB.GT.0)THEN
                      WRITE(LUFILE,'(A)')   '    DO ilmP = 1,nlmP'
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
                         call TRECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB,LUFILE)
                      ELSEIF(UREC)THEN
                         !                         print*,'!A URECURRENCE iTUV',iTUV
                         call URECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB,LUFILE)
                      ELSEIF(VREC)THEN
                         !                         print*,'!A VRECURRENCE iTUV',iTUV
                         call VRECURRENCE(Tb,Ub,Vb,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPA,JMAX,JB,nTUVA,nTUVASTART,nTUVBstart,nTUVB,LUFILE)
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
             WRITE(LUFILE,'(A)')'    ENDDO'
          ENDIF
          WRITE(LUFILE,'(A)')'  ENDDO'
!          IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'
          IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
          IF(JP.LT.10)THEN
             WRITE(LUFILE,'(A,I1,A,I1,A,I1,A)')'end subroutine HorizontalRR_'//ARCSTRING//'_RHS_Q',JP,'C',AngmomA,'D',AngmomB,'CtoD'
          ELSE
             WRITE(LUFILE,'(A,I2,A,I1,A,I1,A)')'end subroutine HorizontalRR_'//ARCSTRING//'_RHS_Q',JP,'C',AngmomA,'D',AngmomB,'CtoD'
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
            & Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Xcd',LUFILE)
    ENDDO

  end subroutine TRECURRENCE

  subroutine XYZRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,&
       & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
       & nTUVQSTART,nTUVQ,iTUVQ,iTUVQminus1x,TMQ,iTUVP,&
       & Tp,Up,Vp,TMP,iTUVplus1x,Jp,DIRSTRING,LUFILE)
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
    integer :: iString,LUFILE
    !step 1 add blanks
    call initString(5)      
    !step 2 determine where to put the 
    IF(iTUVP.LE.nTUVP.AND.((iTUVQ.GE.nTUVQSTART+1).AND.(iTUVQ.LE.nTUVQ)))THEN
       IF(nPrimLAST)THEN
          call AddToString('ThetaP(ilmP,')
       ELSE
          call AddToString('ThetaP(IP,ilmP,')
       ENDIF
    ELSE
       call AddToString('Tmp')
       call AddToString(JTMP)
       call AddToString('(')
    ENDIF
    call AddToString(iTUVP)
    call AddToString(',')
    call AddToString(iTUVQ)
    IF(iTUVP.LE.nTUVP.AND.((iTUVQ.GE.nTUVQSTART+1).AND.(iTUVQ.LE.nTUVQ)))THEN
       IF(nPrimLAST)THEN
          call AddToString(',IP) = ')
       ELSE
          call AddToString(') = ')
       ENDIF
    ELSE
       call AddToString(') = ')
    ENDIF
    !step 3.  the first term: i+1 (p+1,q-1)
    
    IF(iTUVQminus1x.EQ.1)THEN
       IF(nPrimLAST)THEN
          call AddToString('ThetaP2(ilmP,')
       ELSE
          call AddToString('ThetaP2(IP,ilmP,')
       ENDIF
    ELSE
       IF(iTUVplus1x.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          IF(nPrimLAST)THEN
             call AddToString('ThetaP(ilmP,')
          ELSE
             call AddToString('ThetaP(IP,ilmP,')
          ENDIF
       ELSE
          call AddToString('Tmp')
          call AddToString(JTMP-1)
          call AddToString('(')
       ENDIF
    ENDIF
    call AddToString(iTUVplus1x)
    IF(iTUVQminus1x.EQ.1)THEN
       IF(nPrimLAST)THEN
          call AddToString(',')
          call AddToString('IP) ')
       ELSE
          call AddToString(') ')
       ENDIF
    ELSE
       IF(iTUVplus1x.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          call AddToString(iTUVQminus1x)
          IF(nPrimLAST)THEN
             call AddToString(',IP) ')
          ELSE
             call AddToString(') ')
          ENDIF
       ELSE
          call AddToString(',')
          call AddToString(iTUVQminus1x)
          call AddToString(') ')
       ENDIF
    ENDIF
    !step 4: the second term: X*Theta(i,j,k,l)  (p,q-1)
    IF(iTUVQminus1x.EQ.1)THEN
       call AddToString('+ ')
       call AddToString(DIRSTRING)
       IF(nPrimLAST)THEN
          call AddToString('*ThetaP2(ilmP,')
       ELSE
          call AddToString('*ThetaP2(iP,ilmP,')
       ENDIF
    ELSE
       IF(iTUVP.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          call AddToString('+ ')
          call AddToString(DIRSTRING)
          IF(nPrimLAST)THEN
             call AddToString('*ThetaP(ilmP,')
          ELSE
             call AddToString('*ThetaP(iP,ilmP,')
          ENDIF
       ELSE
          call AddToString('+ ')
          call AddToString(DIRSTRING)
          call AddToString('*Tmp')
          call AddToString(JTMP-1)
          call AddToString('(')
       ENDIF
    ENDIF
    call AddToString(iTUVP)
    IF(iTUVQminus1x.EQ.1)THEN
       IF(nPrimLAST)THEN
          call AddToString(',IP) ')
       ELSE
          call AddToString(') ')
       ENDIF
    ELSE
       IF(iTUVP.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          call AddToString(iTUVQminus1x)
          IF(nPrimLAST)THEN
             call AddToString(',IP) ')
          ELSE
             call AddToString(') ')
          ENDIF
       ELSE
          call AddToString(',')
          call AddToString(iTUVQminus1x)
          call AddToString(') ')
       ENDIF
    ENDIF
    !Final step write the string
    call writeString(LUFILE)

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
            & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
            & nTUVQSTART,nTUVQ,iTUVQ,iTUVQminus1x,TMQ,iTUVP,&
            & Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Ycd',LUFILE)
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
            & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,&
            & nTUVQSTART,nTUVQ,iTUVQ,iTUVQminus1x,TMQ,iTUVP,&
            & Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Zcd',LUFILE)
    ENDDO

  end subroutine VRECURRENCE


end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
