MODULE TESTMODULE
CONTAINS
  subroutine PASSsub
    IMPLICIT NONE
    INTEGER :: JMAX,nTUV,nTUVprev,ituvP,J,Tp,Up,Vp,N,N2,ntuvP,ituv,C
    INTEGER :: nTUVQ,nTUVTMPP,nTUVTMPQ,JPQ,JP,JQ,nTUVTMP,nTUVTMPprev
    integer :: tq,uq,vq,ituvq,ituvqminus1x,ituvqminus1y,ituvqminus1z
    integer :: iTUVplus1x,iTUVplus1y,iTUVplus1z,nTUVTMPQs
    Integer :: MaxAngmomQP,nTUVplus,JTMQ,ntuvprev2,ntuvprev3
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

    WRITE(*,'(A)')'MODULE AGC_OBS_TRANSFERRECURRENCEMODAtoC'
    WRITE(*,'(A)')' use IchorPrecisionModule'
    WRITE(*,'(A)')'  '
    WRITE(*,'(A)')' CONTAINS'
    MaxAngmomQP = 8

    DO JMAX=2,MaxAngmomQP
       DO JP = 1, JMAX
          JQ = JMAX-JP
          IF(JQ.GT.JP)CYCLE
          IF(JQ.EQ.0)CYCLE
          IF(JQ.GT.4)CYCLE
          IF(JP.GT.4)CYCLE
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



          WRITE(*,'(A)')''
          IF(JP.LT.10)THEN
             IF(JQ.LT.10)THEN
                WRITE(*,'(A,I1,A,I1,A)')'subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(*,'(A,I1,A,I2,A)')'subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ELSE
             IF(JQ.LT.10)THEN
                WRITE(*,'(A,I2,A,I1,A)')'subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(*,'(A,I2,A,I2,A)')'subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ENDIF

          WRITE(*,'(A)')'         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,Aux,Aux2)'
          WRITE(*,'(A)')'  implicit none'
          WRITE(*,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD'
          WRITE(*,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)'
          WRITE(*,'(A)')'  real(realk),intent(in) :: Pdistance12(3),Qdistance12(3,nPasses)'
          WRITE(*,'(A)')'  real(realk),intent(in) :: Bexp(nPrimB),Dexp(nPrimD)'
          WRITE(*,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPrimQ*nPrimP*nPasses)'
          WRITE(*,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPrimQ*nPrimP*nPasses)'
          WRITE(*,'(A)')'!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ*nPrimP*nPasses)'
          WRITE(*,'(A)')'  !Local variables'
          !             WRITE(*,'(A)')'  real(realk) :: Pexpfac,PREF'
          !             WRITE(*,'(A,i4,A)')'  real(realk) :: TwoTerms(',MAX(1,nTUVprev2-nTUVprev3),')'
          DO JTMQ=1,JQ-1
             nTUVTMPQs=(JTMQ)*(JTMQ+1)*(JTMQ+2)/6
             nTUVTMPQ=(JTMQ+1)*(JTMQ+2)*(JTMQ+3)/6
             nTUVTMPP=(JPQ-JTMQ+1)*(JPQ-JTMQ+2)*(JPQ-JTMQ+3)/6
             if(JTMQ.LT.10)THEN
                WRITE(*,'(A,I1,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMQ,'(',nTUVP+1,':',nTUVTMPP,',',nTUVTMPQs+1,':',nTUVTMPQ,')'
             else
                WRITE(*,'(A,I2,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMQ,'(',nTUVP+1,':',nTUVTMPP,',',nTUVTMPQs+1,':',nTUVTMPQ,')'

             endif
          ENDDO
          WRITE(*,'(A)')'!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering'
          WRITE(*,'(A)')'  integer :: iPassQ,iPrimP,iPrimQ,IP,iTUVP,iPrimB,iPrimD'
          WRITE(*,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk'
          WRITE(*,'(A)')'  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,expBX,expBY,expBZ'
          WRITE(*,'(A)')'  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq'
          WRITE(*,'(A)')'  Xab = Pdistance12(1)'
          WRITE(*,'(A)')'  Yab = Pdistance12(2)'
          WRITE(*,'(A)')'  Zab = Pdistance12(3)'
          WRITE(*,'(A)')'  DO iPassQ = 1,nPasses'
          WRITE(*,'(A)')'   Xcd = Qdistance12(1,iPassQ)'
          WRITE(*,'(A)')'   Ycd = Qdistance12(2,iPassQ)'
          WRITE(*,'(A)')'   Zcd = Qdistance12(3,iPassQ)'
          WRITE(*,'(A)')'   IP = (iPassQ-1)*nPrimQ*nPrimP'
          WRITE(*,'(A)')'   DO iPrimP=1, nPrimP'
          WRITE(*,'(A)')'    expP = Pexp(iPrimP)'
          WRITE(*,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'
          WRITE(*,'(A)')'    expBX = Bexp(iPrimB)*Xab'
          WRITE(*,'(A)')'    expBY = Bexp(iPrimB)*Yab'
          WRITE(*,'(A)')'    expBZ = Bexp(iPrimB)*Zab'
          WRITE(*,'(A)')'    DO iPrimQ=1, nPrimQ'
          WRITE(*,'(A)')'     iPrimD = (iPrimQ-1)/nPrimC+1'
          WRITE(*,'(A)')'     IP = IP + 1'
          WRITE(*,'(A)')'     invexpQ = D1/Qexp(iPrimQ)'
          WRITE(*,'(A)')'     inv2expQ = D05*invexpQ'
          WRITE(*,'(A)')'     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpQ'
          WRITE(*,'(A)')'     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpQ'
          WRITE(*,'(A)')'     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpQ'
          WRITE(*,'(A)')'     pinvq = -expP*invexpQ'
          allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
          CREATED  = .FALSE.
          CREATED(0,0,0) = .TRUE.

          DO JTMQ=0,JQ
!           print*,'!JTMQ = ',JTMQ
           IF(JTMQ.EQ.0)THEN
              nTUVTMPP=(JPQ-JTMQ)*(JPQ-JTMQ+1)*(JPQ-JTMQ+2)/6
              WRITE(*,'(A,I3)')'     DO iTUVP=1,',nTUVP
              WRITE(*,'(A)')   '        Aux2(iTUVP,1,IP) = Aux(iTUVp,IP)'
              WRITE(*,'(A)')   '     ENDDO'
              CYCLE
           ELSE
             nTUVTMPQ=(JTMQ+1)*(JTMQ+2)*(JTMQ+3)/6
             nTUVTMPP=(JPQ-JTMQ+1)*(JPQ-JTMQ+2)*(JPQ-JTMQ+3)/6
             !slowly increase JQ=JTMQ with 1 
             nTUVTMPQ=(JTMQ+1)*(JTMQ+2)*(JTMQ+3)/6
             DO Tq=JTMQ,0,-1       
                DO Uq=JTMQ-Tq,0,-1
                   Vq=JTMQ-Tq-Uq                
!                   print*,'!Tq,Uq,Vq = ',Tq,Uq,Vq
                   !step 1.
                   !how can the (Tp,Up,Vp) be built from lower
                   TREC = CREATED(Tq-1,Uq,Vq)
                   UREC = CREATED(Tq,Uq-1,Vq)
                   VREC = CREATED(Tq,Uq,Vq-1)

                   N=0
                   IF(TREC)N=N+1
                   IF(UREC)N=N+1
                   IF(VREC)N=N+1
                   IF(N.EQ.1)THEN
                      !only one possible way to construct it
                      IF(TREC)THEN
!                         print*,'!A TRECURRENCE iTUV',iTUV
                         call TRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                      ELSEIF(UREC)THEN
!                         print*,'!A URECURRENCE iTUV',iTUV
                         call URECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                      ELSEIF(VREC)THEN
!                         print*,'!A VRECURRENCE iTUV',iTUV
                         call VRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                      ELSE
                         STOP 'TK1'
                      ENDIF
                   ELSE
                      !several ways to construct it
                      TREC2 = CREATED(Tq-2,Uq,Vq)
                      UREC2 = CREATED(Tq,Uq-2,Vq)
                      VREC2 = CREATED(Tq,Uq,Vq-2)
                      N2=0
                      IF(TREC2)N2=N2+1
                      IF(UREC2)N2=N2+1
                      IF(VREC2)N2=N2+1
                      IF(N2.LT.N)THEN
                         !two term recurrence possible for one or more the possibilities 
                         !we chose the one term possibility
                         IF(.NOT.(TREC.AND.TREC2).AND.TREC)THEN
!                            print*,'!B TRECURRENCE iTUV',iTUV
                            call TRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                         ELSEIF(.NOT.(UREC.AND.UREC2).AND.UREC)THEN
!                            print*,'!B URECURRENCE iTUV',iTUV
                            call URECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                         ELSEIF(.NOT.(VREC.AND.VREC2).AND.VREC)THEN
!                            print*,'!B VRECURRENCE iTUV',iTUV
                            call VRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                         ENDIF
                      ELSE
                         !chose one of the possibilities
                         IF(TREC)THEN
!                            print*,'!C TRECURRENCE iTUV',iTUV
                            call TRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                         ELSEIF(UREC)THEN
!                            print*,'!C URECURRENCE iTUV',iTUV
                            call URECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                         ELSEIF(VREC)THEN
!                            print*,'!D VRECURRENCE iTUV',iTUV
                            call VRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
                         ELSE
                            STOP 'TK2'
                         ENDIF
                      ENDIF
                   ENDIF
                   CREATED(Tq,Uq,Vq) = .TRUE.
                   
                ENDDO
             ENDDO
           ENDIF
          ENDDO
!!$                deallocate(TwoTermTUVLIST)
          !WRITE(*,'(A,I4)')'           ENDDO'
!       ENDDO
          WRITE(*,'(A)')'    ENDDO'
          WRITE(*,'(A)')'   ENDDO'
          WRITE(*,'(A)')'  ENDDO'
          IF(JP.LT.10)THEN
             IF(JQ.LT.10)THEN
                WRITE(*,'(A,I1,A,I1,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC'
             ELSE
                WRITE(*,'(A,I1,A,I2,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC'
             ENDIF
          ELSE
             IF(JQ.LT.10)THEN
                WRITE(*,'(A,I2,A,I1,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC'
             ELSE
                WRITE(*,'(A,I2,A,I2,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'AtoC'
             ENDIF
          ENDIF
          deallocate(TUVINDEX)
          deallocate(TINDEX)
          deallocate(UINDEX)
          deallocate(VINDEX)
          deallocate(JINDEX)
       enddo
    enddo
    WRITE(*,'(A)')'end module'
  END subroutine PASSsub
  

subroutine TRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
implicit none
integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMQ,iTUVQ,iTUVQminus1x,iTUVQminus1y,iTUVQminus1z,nTUVP
integer :: iTUVQminus2x, iTUVP,nTUVTMPP,iTUVpminus1x,iTUVpplus1x,JMAX,JQ,TMQ,TMP,JP
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: TINDEX(:)
integer :: UINDEX(:)
integer :: VINDEX(:)
integer :: JINDEX(:)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
character(len=132) :: STRING 
integer :: iString
iTUVQ = TUVINDEX(Tq,Uq,Vq)
iTUVQminus1x = TUVINDEX(Tq-1,Uq,Vq)
iTUVQminus2x = TUVINDEX(Tq-2,Uq,Vq)
TMQ = Tq-1
do iTUVP = 1,nTUVTMPP
   Tp = Tindex(iTUVp) 
   Up = Uindex(iTUVp) 
   Vp = Vindex(iTUVp)
   TMP = Tp
   ituvpplus1x = TUVindex(Tp+1,Up,Vp)
   ituvpminus1x = TUVindex(Tp-1,Up,Vp)
   Jp = Jindex(iTUVp)   

   !A to C
   !Theta(i,0,k+1,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k,0) - p/q*Theta(i+1,0,k,0) + k/(2q)*Theta(i,0,k-1,0) + i/(2q)*Theta(i-1,0,k,0) 

   !step 1 add blanks

   STRING(1:5) = '     '
   iSTRING = 6
   !step 2 determine where to put the 
   IF(iTUVP.LE.nTUVP)THEN
      WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
      iString = iSTRING+5
   ELSE
      IF(JTMQ.LT.10)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ,'('
         iString = iSTRING+5
      ELSE
         WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A)') 'Tmp',JTMQ,'('
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
   IF(iTUVP.LE.nTUVP)THEN
      WRITE(STRING(iSTRING:iSTRING+6),'(A7)') ',IP) = '
      iString = iSTRING+7
   ELSE
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)') ') = '
      iString = iSTRING+4
   ENDIF
   !step 3 determine if the first term: facX*THETA should be included 
   !and if it should use Aux or Tmp 

   IF(iTUVQminus1x.EQ.1)THEN
      !Aux(',iTUVp,',IP)
      WRITE(STRING(iSTRING:iSTRING+8),'(A9)') 'facX*Aux('
      iString = iSTRING+9
   ELSE
      IF(iTUVP.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+9),'(A10)') 'facX*Aux2('
         iString = iSTRING+10
      ELSE
         IF(JTMQ-1.LT.10)THEN
            WRITE(STRING(iSTRING:iSTRING+9),'(A8,I1,A1)') 'facX*Tmp',JTMQ-1,'('
            iString = iSTRING+10
         ELSE
            WRITE(STRING(iSTRING:iSTRING+10),'(A8,I2,A1)') 'facX*Tmp',JTMQ-1,'('
            iString = iSTRING+11
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
      IF(iTUVP.LE.nTUVP)THEN
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
   !step 4 determine if the second term: 
   !  TMP*inv2expQ*Aux(',ituvpminus1x,',IP) 
   !  TMp*inv2expQ*Tmp',JTMQ-1,'(',ituvpminus1x,',',iTUVQminus1x,')
   !should be included and if it should use Aux or Tmp 
   IF(TMp.GT.2)THEN
      IF(TMp.GE.10)THEN
         WRITE(STRING(iSTRING:iSTRING+13),'(A2,I2,A10)') '+ ',TMP,'*inv2expQ*'
         iString = iSTRING+14
      ELSE
         WRITE(STRING(iSTRING:iSTRING+12),'(A2,I1,A10)') '+ ',TMP,'*inv2expQ*'
         iString = iSTRING+13
      ENDIF
   ELSEIF(TMp.EQ.2)THEN
      WRITE(STRING(iSTRING:iSTRING+9),'(A10)') '+ invexpQ*'
      iString = iSTRING+10
   ELSEIF(TMp.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+10),'(A11)') '+ inv2expQ*'
      iString = iSTRING+11
   ELSE
      !do not include this term 
   ENDIF
   
   IF(TMp.GT.0)THEN
      IF(iTUVQminus1x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
         iString = iSTRING+4
      ELSE
         IF(ituvpminus1x.LE.nTUVP)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
            iString = iSTRING+5
         ELSE
            IF(JTMQ-1.LT.10)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-1,'('
               iString = iSTRING+5
            ELSE
               WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMQ-1,'('
               iString = iSTRING+6
            ENDIF
         ENDIF
      ENDIF
      IF(ituvpminus1x.LT.100)THEN
         WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') ituvpminus1x,','
         iString = iSTRING+3
      ELSEIF(ituvpminus1x.LT.1000)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') ituvpminus1x,','
         iString = iSTRING+4
      ELSEIF(ituvpminus1x.LT.10000)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') ituvpminus1x,','
         iString = iSTRING+5
      ELSE
         STOP 'Recurrent ituvpminus1x'
      ENDIF
      IF(iTUVQminus1x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
         iString = iSTRING+4
      ELSE
         IF(ituvpminus1x.LE.nTUVP)THEN
            IF(iTUVQminus1x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+6),'(I2,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+7
            ELSEIF(iTUVQminus1x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+7),'(I3,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+8
            ELSEIF(iTUVQminus1x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+8),'(I4,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+9
            ELSE
               STOP 'Recurrent iTUVQminus1x'
            ENDIF
         ELSE
            IF(iTUVQminus1x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+3),'(I2,A2)') iTUVQminus1x,') '
               iString = iSTRING+4
            ELSEIF(iTUVQminus1x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(I3,A2)') iTUVQminus1x,') '
               iString = iSTRING+5
            ELSEIF(iTUVQminus1x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+5),'(I4,A2)') iTUVQminus1x,') '
               iString = iSTRING+6
            ELSE
               STOP 'Recurrent iTUVQminus1x'
            ENDIF
         ENDIF
      ENDIF
   ELSE
      !do not include this term       
   ENDIF
   !step 5 determine if the third term: 
   !  Tmq*inv2expQ*Aux(',iTUVp,',IP)
   !  Tmq*inv2expQ*Tmp',JTMQ-2,'(',iTUVp,',',iTUVQminus2x,')
   !should be included and if it should use Aux or Tmp 
   IF(TMq.GT.2)THEN
      IF(TMq.GE.10)THEN
         WRITE(STRING(iSTRING:iSTRING+13),'(A2,I2,A10)') '+ ',TMQ,'*inv2expQ*'
         iString = iSTRING+14
      ELSE
         WRITE(STRING(iSTRING:iSTRING+12),'(A2,I1,A10)') '+ ',TMQ,'*inv2expQ*'
         iString = iSTRING+13
      ENDIF
   ELSEIF(TMQ.EQ.2)THEN
      WRITE(STRING(iSTRING:iSTRING+9),'(A10)') '+ invexpQ*'
      iString = iSTRING+10
   ELSEIF(TMQ.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+10),'(A11)') '+ inv2expQ*'
      iString = iSTRING+11
   ELSE
      !do not include this term 
   ENDIF

   IF(TMq.GT.0)THEN
      IF(iTUVQminus2x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
         iString = iSTRING+4
      ELSE
         IF(iTUVP.LE.nTUVP)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
            iString = iSTRING+5
         ELSE
            IF(JTMQ-1.LT.10)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-2,'('
               iString = iSTRING+5
            ELSE
               WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMQ-2,'('
               iString = iSTRING+6
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
      IF(iTUVQminus2x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
         iString = iSTRING+4
      ELSE
         IF(iTUVP.LE.nTUVP)THEN
            IF(iTUVQminus2x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+6),'(I2,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+7
            ELSEIF(iTUVQminus2x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+7),'(I3,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+8
            ELSEIF(iTUVQminus2x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+8),'(I4,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+9
            ELSE
               STOP 'Recurrent iTUVQminus2x'
            ENDIF
         ELSE
            IF(iTUVQminus2x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+3),'(I2,A2)') iTUVQminus2x,') '
               iString = iSTRING+4
            ELSEIF(iTUVQminus2x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(I3,A2)') iTUVQminus2x,') '
               iString = iSTRING+5
            ELSEIF(iTUVQminus2x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+5),'(I4,A2)') iTUVQminus2x,') '
               iString = iSTRING+6
            ELSE
               STOP 'Recurrent iTUVQminus2x'
            ENDIF
         ENDIF
      ENDIF
   ELSE
      !do not include this term       
   ENDIF
   !step 6 determine if the third term: 
   !  pinvq*Aux(',ituvpplus1x,',IP)'
   !  pinvq*Tmp',JTMQ-1,'(',ituvpplus1x,',',iTUVQminus1x,')'
   !should be included and if it should use Aux or Tmp 
   WRITE(STRING(iSTRING:iSTRING+7),'(A8)') '+ pinvq*'
   iString = iSTRING+8
   IF(iTUVQminus1x.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
      iString = iSTRING+4
   ELSE
      IF(ituvpplus1x.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
         iString = iSTRING+5
      ELSE
         IF(JTMQ-1.LT.10)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-1,'('
            iString = iSTRING+5
         ELSE
            WRITE(STRING(iSTRING:iSTRING+4),'(A3,I2,A1)') 'Tmp',JTMQ-1,'('
            iString = iSTRING+5
         ENDIF
      ENDIF
   ENDIF
   IF(ituvpplus1x.LT.100)THEN
      WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') ituvpplus1x,','
      iString = iSTRING+3
   ELSEIF(ituvpplus1x.LT.1000)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') ituvpplus1x,','
      iString = iSTRING+4
   ELSEIF(ituvpplus1x.LT.10000)THEN
      WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') ituvpplus1x,','
      iString = iSTRING+5
   ELSE
      STOP 'Recurrent ituvpplus1x'
   ENDIF
   IF(iTUVQminus1x.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
      iString = iSTRING+4
   ELSE
      IF(iTUVQminus1x.LT.100)THEN
         WRITE(STRING(iSTRING:iSTRING+1),'(I2)') iTUVQminus1x
         iString = iSTRING+2
      ELSEIF(iTUVQminus1x.LT.1000)THEN
         WRITE(STRING(iSTRING:iSTRING+2),'(I3)') iTUVQminus1x
         iString = iSTRING+3
      ELSEIF(iTUVQminus1x.LT.10000)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(I4)') iTUVQminus1x
         iString = iSTRING+4
      ELSE
         STOP 'Recurrent iTUVQminus2x'
      ENDIF
      IF(ituvpplus1x.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A5)') ',IP) '
         iString = iSTRING+5
      ELSE
         WRITE(STRING(iSTRING:iSTRING+1),'(A2)') ') '
         iString = iSTRING+2
      ENDIF
   ENDIF
   !Final step write the string
   WRITE(*,'(A)') STRING(1:iSTRING-1)
ENDDO

end subroutine TRECURRENCE

subroutine URECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
implicit none
integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMQ,iTUVQ,iTUVQminus1x,iTUVQminus1y,iTUVQminus1z,nTUVP
integer :: iTUVQminus2x, iTUVP,nTUVTMPP,ituvpminus1x,ituvpplus1x,JMAX,JQ,TMQ,TMP,JP
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
do iTUVP = 1,nTUVTMPP
   Tp = Tindex(iTUVp) 
   Up = Uindex(iTUVp) 
   Vp = Vindex(iTUVp)
   TMP = Up
   ituvpplus1x = TUVindex(Tp,Up+1,Vp)
   ituvpminus1x = TUVindex(Tp,Up-1,Vp)
   Jp = Jindex(iTUVp)   

   !step 1 add blanks

   STRING(1:5) = '     '
   iSTRING = 6
   !step 2 determine where to put the 
   IF(iTUVP.LE.nTUVP)THEN
      WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
      iString = iSTRING+5
   ELSE
      IF(JTMQ.LT.10)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ,'('
         iString = iSTRING+5
      ELSE
         WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A)') 'Tmp',JTMQ,'('
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
   IF(iTUVP.LE.nTUVP)THEN
      WRITE(STRING(iSTRING:iSTRING+6),'(A7)') ',IP) = '
      iString = iSTRING+7
   ELSE
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)') ') = '
      iString = iSTRING+4
   ENDIF
   !step 3 determine if the first term: facX*THETA should be included 
   !and if it should use Aux or Tmp 

   IF(iTUVQminus1x.EQ.1)THEN
      !Aux(',iTUVp,',IP)
      WRITE(STRING(iSTRING:iSTRING+8),'(A9)') 'facY*Aux('
      iString = iSTRING+9
   ELSE
      IF(iTUVP.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+9),'(A10)') 'facY*Aux2('
         iString = iSTRING+10
      ELSE
         IF(JTMQ-1.LT.10)THEN
            WRITE(STRING(iSTRING:iSTRING+9),'(A8,I1,A1)') 'facY*Tmp',JTMQ-1,'('
            iString = iSTRING+10
         ELSE
            WRITE(STRING(iSTRING:iSTRING+10),'(A8,I2,A1)') 'facY*Tmp',JTMQ-1,'('
            iString = iSTRING+11
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
      IF(iTUVP.LE.nTUVP)THEN
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
   !step 4 determine if the second term: 
   !  TMP*inv2expQ*Aux(',ituvpminus1x,',IP) 
   !  TMp*inv2expQ*Tmp',JTMQ-1,'(',ituvpminus1x,',',iTUVQminus1x,')
   !should be included and if it should use Aux or Tmp 
   IF(TMp.GT.2)THEN
      IF(TMp.GE.10)THEN
         WRITE(STRING(iSTRING:iSTRING+13),'(A2,I2,A10)') '+ ',TMP,'*inv2expQ*'
         iString = iSTRING+14
      ELSE
         WRITE(STRING(iSTRING:iSTRING+12),'(A2,I1,A10)') '+ ',TMP,'*inv2expQ*'
         iString = iSTRING+13
      ENDIF
   ELSEIF(TMp.EQ.2)THEN
      WRITE(STRING(iSTRING:iSTRING+9),'(A10)') '+ invexpQ*'
      iString = iSTRING+10
   ELSEIF(TMp.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+10),'(A11)') '+ inv2expQ*'
      iString = iSTRING+11
   ELSE
      !do not include this term 
   ENDIF
   
   IF(TMp.GT.0)THEN
      IF(iTUVQminus1x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
         iString = iSTRING+4
      ELSE
         IF(ituvpminus1x.LE.nTUVP)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
            iString = iSTRING+5
         ELSE
            IF(JTMQ-1.LT.10)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-1,'('
               iString = iSTRING+5
            ELSE
               WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMQ-1,'('
               iString = iSTRING+6
            ENDIF
         ENDIF
      ENDIF
      IF(ituvpminus1x.LT.100)THEN
         WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') ituvpminus1x,','
         iString = iSTRING+3
      ELSEIF(ituvpminus1x.LT.1000)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') ituvpminus1x,','
         iString = iSTRING+4
      ELSEIF(ituvpminus1x.LT.10000)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') ituvpminus1x,','
         iString = iSTRING+5
      ELSE
         STOP 'Recurrent ituvpminus1x'
      ENDIF
      IF(iTUVQminus1x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
         iString = iSTRING+4
      ELSE
         IF(ituvpminus1x.LE.nTUVP)THEN
            IF(iTUVQminus1x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+6),'(I2,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+7
            ELSEIF(iTUVQminus1x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+7),'(I3,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+8
            ELSEIF(iTUVQminus1x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+8),'(I4,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+9
            ELSE
               STOP 'Recurrent iTUVQminus1x'
            ENDIF
         ELSE
            IF(iTUVQminus1x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+3),'(I2,A2)') iTUVQminus1x,') '
               iString = iSTRING+4
            ELSEIF(iTUVQminus1x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(I3,A2)') iTUVQminus1x,') '
               iString = iSTRING+5
            ELSEIF(iTUVQminus1x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+5),'(I4,A2)') iTUVQminus1x,') '
               iString = iSTRING+6
            ELSE
               STOP 'Recurrent iTUVQminus1x'
            ENDIF
         ENDIF
      ENDIF
   ELSE
      !do not include this term       
   ENDIF
   !step 5 determine if the third term: 
   !  Tmq*inv2expQ*Aux(',iTUVp,',IP)
   !  Tmq*inv2expQ*Tmp',JTMQ-2,'(',iTUVp,',',iTUVQminus2x,')
   !should be included and if it should use Aux or Tmp 
   IF(TMq.GT.2)THEN
      IF(TMq.GE.10)THEN
         WRITE(STRING(iSTRING:iSTRING+13),'(A2,I2,A10)') '+ ',TMQ,'*inv2expQ*'
         iString = iSTRING+14
      ELSE
         WRITE(STRING(iSTRING:iSTRING+12),'(A2,I1,A10)') '+ ',TMQ,'*inv2expQ*'
         iString = iSTRING+13
      ENDIF
   ELSEIF(TMQ.EQ.2)THEN
      WRITE(STRING(iSTRING:iSTRING+9),'(A10)') '+ invexpQ*'
      iString = iSTRING+10
   ELSEIF(TMQ.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+10),'(A11)') '+ inv2expQ*'
      iString = iSTRING+11
   ELSE
      !do not include this term 
   ENDIF

   IF(TMq.GT.0)THEN
      IF(iTUVQminus2x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
         iString = iSTRING+4
      ELSE
         IF(iTUVP.LE.nTUVP)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
            iString = iSTRING+5
         ELSE
            IF(JTMQ-1.LT.10)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-2,'('
               iString = iSTRING+5
            ELSE
               WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMQ-2,'('
               iString = iSTRING+6
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
      IF(iTUVQminus2x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
         iString = iSTRING+4
      ELSE
         IF(iTUVP.LE.nTUVP)THEN
            IF(iTUVQminus2x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+6),'(I2,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+7
            ELSEIF(iTUVQminus2x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+7),'(I3,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+8
            ELSEIF(iTUVQminus2x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+8),'(I4,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+9
            ELSE
               STOP 'Recurrent iTUVQminus2x'
            ENDIF
         ELSE
            IF(iTUVQminus2x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+3),'(I2,A2)') iTUVQminus2x,') '
               iString = iSTRING+4
            ELSEIF(iTUVQminus2x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(I3,A2)') iTUVQminus2x,') '
               iString = iSTRING+5
            ELSEIF(iTUVQminus2x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+5),'(I4,A2)') iTUVQminus2x,') '
               iString = iSTRING+6
            ELSE
               STOP 'Recurrent iTUVQminus2x'
            ENDIF
         ENDIF
      ENDIF
   ELSE
      !do not include this term       
   ENDIF
   !step 6 determine if the third term: 
   !  pinvq*Aux(',ituvpplus1x,',IP)'
   !  pinvq*Tmp',JTMQ-1,'(',ituvpplus1x,',',iTUVQminus1x,')'
   !should be included and if it should use Aux or Tmp 
   WRITE(STRING(iSTRING:iSTRING+7),'(A8)') '+ pinvq*'
   iString = iSTRING+8
   IF(iTUVQminus1x.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
      iString = iSTRING+4
   ELSE
      IF(ituvpplus1x.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
         iString = iSTRING+5
      ELSE
         IF(JTMQ-1.LT.10)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-1,'('
            iString = iSTRING+5
         ELSE
            WRITE(STRING(iSTRING:iSTRING+4),'(A3,I2,A1)') 'Tmp',JTMQ-1,'('
            iString = iSTRING+5
         ENDIF
      ENDIF
   ENDIF
   IF(ituvpplus1x.LT.100)THEN
      WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') ituvpplus1x,','
      iString = iSTRING+3
   ELSEIF(ituvpplus1x.LT.1000)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') ituvpplus1x,','
      iString = iSTRING+4
   ELSEIF(ituvpplus1x.LT.10000)THEN
      WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') ituvpplus1x,','
      iString = iSTRING+5
   ELSE
      STOP 'Recurrent ituvpplus1x'
   ENDIF
   IF(iTUVQminus1x.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
      iString = iSTRING+4
   ELSE
      IF(iTUVQminus1x.LT.100)THEN
         WRITE(STRING(iSTRING:iSTRING+1),'(I2)') iTUVQminus1x
         iString = iSTRING+2
      ELSEIF(iTUVQminus1x.LT.1000)THEN
         WRITE(STRING(iSTRING:iSTRING+2),'(I3)') iTUVQminus1x
         iString = iSTRING+3
      ELSEIF(iTUVQminus1x.LT.10000)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(I4)') iTUVQminus1x
         iString = iSTRING+4
      ELSE
         STOP 'Recurrent iTUVQminus2x'
      ENDIF
      IF(ituvpplus1x.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A5)') ',IP) '
         iString = iSTRING+5
      ELSE
         WRITE(STRING(iSTRING:iSTRING+1),'(A2)') ') '
         iString = iSTRING+2
      ENDIF
   ENDIF
   !Final step write the string
   WRITE(*,'(A)') STRING(1:iSTRING-1)
ENDDO

end subroutine URECURRENCE

subroutine VRECURRENCE(Tq,Uq,Vq,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMQ,nTUVTMPP,JMAX,JQ,nTUVP)
implicit none
integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMQ,iTUVQ,iTUVQminus1x,iTUVQminus1y,iTUVQminus1z
integer :: iTUVQminus2x, iTUVP,nTUVTMPP,ituvpminus1x,ituvpplus1x,JMAX,JQ,TMQ,TMP,JP,nTUVP
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
do iTUVP = 1,nTUVTMPP
   Tp = Tindex(iTUVp) 
   Up = Uindex(iTUVp) 
   Vp = Vindex(iTUVp)
   TMP = Vp
   ituvpplus1x = TUVindex(Tp,Up,Vp+1)
   ituvpminus1x = TUVindex(Tp,Up,Vp-1)
   Jp = Jindex(iTUVp)   

   !step 1 add blanks

   STRING(1:5) = '     '
   iSTRING = 6
   !step 2 determine where to put the 
   IF(iTUVP.LE.nTUVP)THEN
      WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
      iString = iSTRING+5
   ELSE
      IF(JTMQ.LT.10)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ,'('
         iString = iSTRING+5
      ELSE
         WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A)') 'Tmp',JTMQ,'('
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
   IF(iTUVP.LE.nTUVP)THEN
      WRITE(STRING(iSTRING:iSTRING+6),'(A7)') ',IP) = '
      iString = iSTRING+7
   ELSE
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)') ') = '
      iString = iSTRING+4
   ENDIF
   !step 3 determine if the first term: facX*THETA should be included 
   !and if it should use Aux or Tmp 

   IF(iTUVQminus1x.EQ.1)THEN
      !Aux(',iTUVp,',IP)
      WRITE(STRING(iSTRING:iSTRING+8),'(A9)') 'facZ*Aux('
      iString = iSTRING+9
   ELSE
      IF(iTUVP.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+9),'(A10)') 'facZ*Aux2('
         iString = iSTRING+10
      ELSE
         IF(JTMQ-1.LT.10)THEN
            WRITE(STRING(iSTRING:iSTRING+9),'(A8,I1,A1)') 'facZ*Tmp',JTMQ-1,'('
            iString = iSTRING+10
         ELSE
            WRITE(STRING(iSTRING:iSTRING+10),'(A8,I2,A1)') 'facZ*Tmp',JTMQ-1,'('
            iString = iSTRING+11
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
      IF(iTUVP.LE.nTUVP)THEN
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
   !step 4 determine if the second term: 
   !  TMP*inv2expQ*Aux(',ituvpminus1x,',IP) 
   !  TMp*inv2expQ*Tmp',JTMQ-1,'(',ituvpminus1x,',',iTUVQminus1x,')
   !should be included and if it should use Aux or Tmp 
   IF(TMp.GT.2)THEN
      IF(TMp.GE.10)THEN
         WRITE(STRING(iSTRING:iSTRING+13),'(A2,I2,A10)') '+ ',TMP,'*inv2expQ*'
         iString = iSTRING+14
      ELSE
         WRITE(STRING(iSTRING:iSTRING+12),'(A2,I1,A10)') '+ ',TMP,'*inv2expQ*'
         iString = iSTRING+13
      ENDIF
   ELSEIF(TMp.EQ.2)THEN
      WRITE(STRING(iSTRING:iSTRING+9),'(A10)') '+ invexpQ*'
      iString = iSTRING+10
   ELSEIF(TMp.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+10),'(A11)') '+ inv2expQ*'
      iString = iSTRING+11
   ELSE
      !do not include this term 
   ENDIF
   
   IF(TMp.GT.0)THEN
      IF(iTUVQminus1x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
         iString = iSTRING+4
      ELSE
         IF(ituvpminus1x.LE.nTUVP)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
            iString = iSTRING+5
         ELSE
            IF(JTMQ-1.LT.10)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-1,'('
               iString = iSTRING+5
            ELSE
               WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMQ-1,'('
               iString = iSTRING+6
            ENDIF
         ENDIF
      ENDIF
      IF(ituvpminus1x.LT.100)THEN
         WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') ituvpminus1x,','
         iString = iSTRING+3
      ELSEIF(ituvpminus1x.LT.1000)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') ituvpminus1x,','
         iString = iSTRING+4
      ELSEIF(ituvpminus1x.LT.10000)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') ituvpminus1x,','
         iString = iSTRING+5
      ELSE
         STOP 'Recurrent ituvpminus1x'
      ENDIF
      IF(iTUVQminus1x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
         iString = iSTRING+4
      ELSE
         IF(ituvpminus1x.LE.nTUVP)THEN
            IF(iTUVQminus1x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+6),'(I2,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+7
            ELSEIF(iTUVQminus1x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+7),'(I3,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+8
            ELSEIF(iTUVQminus1x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+8),'(I4,A5)') iTUVQminus1x,',IP) '
               iString = iSTRING+9
            ELSE
               STOP 'Recurrent iTUVQminus1x'
            ENDIF
         ELSE
            IF(iTUVQminus1x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+3),'(I2,A2)') iTUVQminus1x,') '
               iString = iSTRING+4
            ELSEIF(iTUVQminus1x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(I3,A2)') iTUVQminus1x,') '
               iString = iSTRING+5
            ELSEIF(iTUVQminus1x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+5),'(I4,A2)') iTUVQminus1x,') '
               iString = iSTRING+6
            ELSE
               STOP 'Recurrent iTUVQminus1x'
            ENDIF
         ENDIF
      ENDIF
   ELSE
      !do not include this term       
   ENDIF
   !step 5 determine if the third term: 
   !  Tmq*inv2expQ*Aux(',iTUVp,',IP)
   !  Tmq*inv2expQ*Tmp',JTMQ-2,'(',iTUVp,',',iTUVQminus2x,')
   !should be included and if it should use Aux or Tmp 
   IF(TMq.GT.2)THEN
      IF(TMq.GE.10)THEN
         WRITE(STRING(iSTRING:iSTRING+13),'(A2,I2,A10)') '+ ',TMQ,'*inv2expQ*'
         iString = iSTRING+14
      ELSE
         WRITE(STRING(iSTRING:iSTRING+12),'(A2,I1,A10)') '+ ',TMQ,'*inv2expQ*'
         iString = iSTRING+13
      ENDIF
   ELSEIF(TMQ.EQ.2)THEN
      WRITE(STRING(iSTRING:iSTRING+9),'(A10)') '+ invexpQ*'
      iString = iSTRING+10
   ELSEIF(TMQ.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+10),'(A11)') '+ inv2expQ*'
      iString = iSTRING+11
   ELSE
      !do not include this term 
   ENDIF

   IF(TMq.GT.0)THEN
      IF(iTUVQminus2x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
         iString = iSTRING+4
      ELSE
         IF(iTUVP.LE.nTUVP)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
            iString = iSTRING+5
         ELSE
            IF(JTMQ-1.LT.10)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-2,'('
               iString = iSTRING+5
            ELSE
               WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMQ-2,'('
               iString = iSTRING+6
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
      IF(iTUVQminus2x.EQ.1)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
         iString = iSTRING+4
      ELSE
         IF(iTUVP.LE.nTUVP)THEN
            IF(iTUVQminus2x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+6),'(I2,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+7
            ELSEIF(iTUVQminus2x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+7),'(I3,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+8
            ELSEIF(iTUVQminus2x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+8),'(I4,A5)') iTUVQminus2x,',IP) '
               iString = iSTRING+9
            ELSE
               STOP 'Recurrent iTUVQminus2x'
            ENDIF
         ELSE
            IF(iTUVQminus2x.LT.100)THEN
               WRITE(STRING(iSTRING:iSTRING+3),'(I2,A2)') iTUVQminus2x,') '
               iString = iSTRING+4
            ELSEIF(iTUVQminus2x.LT.1000)THEN
               WRITE(STRING(iSTRING:iSTRING+4),'(I3,A2)') iTUVQminus2x,') '
               iString = iSTRING+5
            ELSEIF(iTUVQminus2x.LT.10000)THEN
               WRITE(STRING(iSTRING:iSTRING+5),'(I4,A2)') iTUVQminus2x,') '
               iString = iSTRING+6
            ELSE
               STOP 'Recurrent iTUVQminus2x'
            ENDIF
         ENDIF
      ENDIF
   ELSE
      !do not include this term       
   ENDIF
   !step 6 determine if the third term: 
   !  pinvq*Aux(',ituvpplus1x,',IP)'
   !  pinvq*Tmp',JTMQ-1,'(',ituvpplus1x,',',iTUVQminus1x,')'
   !should be included and if it should use Aux or Tmp 
   WRITE(STRING(iSTRING:iSTRING+7),'(A8)') '+ pinvq*'
   iString = iSTRING+8
   IF(iTUVQminus1x.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
      iString = iSTRING+4
   ELSE
      IF(ituvpplus1x.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
         iString = iSTRING+5
      ELSE
         IF(JTMQ-1.LT.10)THEN
            WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMQ-1,'('
            iString = iSTRING+5
         ELSE
            WRITE(STRING(iSTRING:iSTRING+4),'(A3,I2,A1)') 'Tmp',JTMQ-1,'('
            iString = iSTRING+5
         ENDIF
      ENDIF
   ENDIF
   IF(ituvpplus1x.LT.100)THEN
      WRITE(STRING(iSTRING:iSTRING+2),'(I2,A)') ituvpplus1x,','
      iString = iSTRING+3
   ELSEIF(ituvpplus1x.LT.1000)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(I3,A)') ituvpplus1x,','
      iString = iSTRING+4
   ELSEIF(ituvpplus1x.LT.10000)THEN
      WRITE(STRING(iSTRING:iSTRING+4),'(I4,A)') ituvpplus1x,','
      iString = iSTRING+5
   ELSE
      STOP 'Recurrent ituvpplus1x'
   ENDIF
   IF(iTUVQminus1x.EQ.1)THEN
      WRITE(STRING(iSTRING:iSTRING+3),'(A4)')'IP) '
      iString = iSTRING+4
   ELSE
      IF(iTUVQminus1x.LT.100)THEN
         WRITE(STRING(iSTRING:iSTRING+1),'(I2)') iTUVQminus1x
         iString = iSTRING+2
      ELSEIF(iTUVQminus1x.LT.1000)THEN
         WRITE(STRING(iSTRING:iSTRING+2),'(I3)') iTUVQminus1x
         iString = iSTRING+3
      ELSEIF(iTUVQminus1x.LT.10000)THEN
         WRITE(STRING(iSTRING:iSTRING+3),'(I4)') iTUVQminus1x
         iString = iSTRING+4
      ELSE
         STOP 'Recurrent iTUVQminus2x'
      ENDIF
      IF(ituvpplus1x.LE.nTUVP)THEN
         WRITE(STRING(iSTRING:iSTRING+4),'(A5)') ',IP) '
         iString = iSTRING+5
      ELSE
         WRITE(STRING(iSTRING:iSTRING+1),'(A2)') ') '
         iString = iSTRING+2
      ENDIF
   ENDIF
   !Final step write the string
   WRITE(*,'(A)') STRING(1:iSTRING-1)
ENDDO

end subroutine VRECURRENCE

end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
