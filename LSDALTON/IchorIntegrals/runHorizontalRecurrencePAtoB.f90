MODULE TESTMODULE
CONTAINS
!TODO you can easily change
!  real(realk),intent(in) :: AuxCont(    1,nTUVQ,nContPasses)
!  real(realk),intent(inout) :: ThetaP(    1:    1,    1:    1,nTUVQ,nContPasses)
!  !Local variables
!  integer :: iP,iTUVQ,iTUVA
!  DO iP = 1,nContPasses
!   DO iTUVQ = 1,nTUVQ
!     DO iTUVA=  1,  1
!        ThetaP(iTUVA,1,iTUVQ,IP) = AuxCont(iTUVA,iTUVQ,IP)

!  to :
!  real(realk),intent(in) :: AuxCont(    1,nTUVQ*nContPasses)
!  real(realk),intent(inout) :: ThetaP(    1:    1,    1:    1,nTUVQ*nContPasses)
!  !Local variables
!  integer :: iP,iTUVQ,iTUVA
!  DO iP = 1,nContPasses*nTUVQ
!     DO iTUVA=  1,  1
!        ThetaP(iTUVA,1,IP) = AuxCont(iTUVA,IP)

!better OpenMP loadbalance 


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

    WRITE(*,'(A)')'MODULE AGC_OBS_HorizontalRecurrenceLHSModAtoB'
    WRITE(*,'(A)')' use IchorPrecisionModule'
    WRITE(*,'(A)')'  '
    WRITE(*,'(A)')' CONTAINS'
    MaxAngmomP = 6


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
          IF(AngmomA.GT.3)CYCLE
          IF(AngmomB.GT.3)CYCLE

          NTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
          NTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6

          NTUVAspec = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6 - &
               &         (AngmomA)*(AngmomA+1)*(AngmomA+2)/6
          NTUVBspec = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6 - &
               &         (AngmomB)*(AngmomB+1)*(AngmomB+2)/6
          NTUVAstart = (AngmomA)*(AngmomA+1)*(AngmomA+2)/6
          NTUVBstart = (AngmomB)*(AngmomB+1)*(AngmomB+2)/6

          WRITE(*,'(A)')''
          IF(JP.LT.10)THEN
             WRITE(*,'(A,I1,A,I1,A,I1,A)')'subroutine HorizontalRR_LHS_P',JP,'A',AngmomA,'B',AngmomB,'AtoB(nContPasses,nTUVQ,&'
          ELSE
             WRITE(*,'(A,I2,A,I1,A,I1,A)')'subroutine HorizontalRR_LHS_P',JP,'A',AngmomA,'B',AngmomB,'AtoB(nContPasses,nTUVQ,&'
          ENDIF
          WRITE(*,'(A)')'         & Pdistance12,AuxCont,ThetaP,lupri)'
          WRITE(*,'(A)')'  implicit none'
          WRITE(*,'(A)')'  integer,intent(in) :: nContPasses,nTUVQ,lupri'
          WRITE(*,'(A)')'  real(realk),intent(in) :: Pdistance12(3)'
          WRITE(*,'(A,I5,A)')'  real(realk),intent(in) :: AuxCont(',nTUVP,',nTUVQ,nContPasses)'
          !             WRITE(*,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(',nTUVAspec,',',nTUVBspec,',nTUVQ,nContPasses)'
          IF(nTUVBstart+1.EQ.1.AND.nTUVB.EQ.1)THEN
             WRITE(*,'(A,I5,A1,I5,A)')'  real(realk),intent(inout) :: ThetaP(',NTUVAstart+1,':',nTUVA,',1,nTUVQ,nContPasses)'
          ELSE
             WRITE(*,'(A,I5,A1,I5,A,I5,A,I5,A)')'  real(realk),intent(inout) :: ThetaP(',NTUVAstart+1,':',nTUVA,',',nTUVBstart+1,':',nTUVB,',nTUVQ,nContPasses)'
          ENDIF
          WRITE(*,'(A)')'  !Local variables'
          IF(ANGMOMB.NE.0)THEN
             WRITE(*,'(A)')'  integer :: iP,iTUVQ,iTUVA'
             WRITE(*,'(A)')'  real(realk) :: Xab,Yab,Zab'
          ELSE
             WRITE(*,'(A)')'  integer :: iP,iTUVQ,iTUVA'
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
          IF(JB.GT.1)THEN
             WRITE(*,'(A)')'!  real(realk) :: Tmp(nTUVA,nTUVB) ordering'
          ENDIF
          IF(JB.NE.0)THEN
             WRITE(*,'(A)')'  Xab = Pdistance12(1)'
             WRITE(*,'(A)')'  Yab = Pdistance12(2)'
             WRITE(*,'(A)')'  Zab = Pdistance12(3)'
          ENDIF
          WRITE(*,'(A)')'  DO iP = 1,nContPasses'
          WRITE(*,'(A)')'   DO iTUVQ = 1,nTUVQ'


          DO JTMP=0,JB
             !           print*,'!JTMP = ',JTMP
             IF(JTMP.EQ.0)THEN
                !nTUVTMPP=(Jab-JTMP)*(Jab-JTMP+1)*(Jab-JTMP+2)/6
                IF(nTUVBstart+1.EQ.1)THEN
                   WRITE(*,'(A,I3,A,I3)')'     DO iTUVA=',NTUVAstart+1,',',nTUVA
                   WRITE(*,'(A)')   '        ThetaP(iTUVA,1,iTUVQ,IP) = AuxCont(iTUVA,iTUVQ,IP)'
                   WRITE(*,'(A)')   '     ENDDO'
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
          WRITE(*,'(A)')'   ENDDO'
          WRITE(*,'(A)')'  ENDDO'
          IF(JP.LT.10)THEN
             WRITE(*,'(A,I1,A,I1,A,I1,A)')'end subroutine HorizontalRR_LHS_P',JP,'A',AngmomA,'B',AngmomB,'AtoB'
          ELSE
             WRITE(*,'(A,I2,A,I1,A,I1,A)')'end subroutine HorizontalRR_LHS_P',JP,'A',AngmomA,'B',AngmomB,'AtoB'
          ENDIF
       enddo
       deallocate(TUVINDEX)
       deallocate(TINDEX)
       deallocate(UINDEX)
       deallocate(VINDEX)
       deallocate(JINDEX)
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
            & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,nTUVQSTART,nTUVQ,&
            & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Xab')

    ENDDO

  end subroutine TRECURRENCE

  subroutine XYZRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,&
       & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,nTUVQSTART,nTUVQ,&
       & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,DIRSTRING)
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
    integer :: iString
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
       WRITE(STRING(iSTRING:iSTRING+12),'(A13)') ',iTUVQ,IP) = '
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
       WRITE(STRING(iSTRING:iSTRING+9),'(A10)') 'iTUVQ,IP) '
       iString = iSTRING+10
    ELSE
       IF(iTUVplus1x.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+13),'(I2,A12)') iTUVQminus1x,',iTUVQ,IP) '
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
       WRITE(STRING(iSTRING:iSTRING+9),'(A10)') 'iTUVQ,IP) '
       iString = iSTRING+10
    ELSE
       IF(iTUVP.LE.nTUVP.AND.((iTUVQminus1x.GE.nTUVQSTART+1).AND.(iTUVQminus1x.LE.nTUVQ)))THEN
          WRITE(STRING(iSTRING:iSTRING+12),'(I2,A11)') iTUVQminus1x,',iTUVQ,IP) '
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
        & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,nTUVQSTART,nTUVQ,&
        & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Yab')
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
        & CREATED,JTMP,nTUVTMPP,JMAX,JQ,nTUVP,nTUVASTART,nTUVQSTART,nTUVQ,&
        & iTUVQ,iTUVQminus1x,TMQ,iTUVP,Tp,Up,Vp,TMP,iTUVplus1x,Jp,'Zab')
ENDDO

end subroutine VRECURRENCE


end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
