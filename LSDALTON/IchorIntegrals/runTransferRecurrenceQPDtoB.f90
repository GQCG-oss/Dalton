MODULE TESTMODULE
CONTAINS
  subroutine PASSsub
    IMPLICIT NONE
    INTEGER :: JMAX,nTUV,nTUVprev,ituvP,J,Tp,Up,Vp,N,N2,ntuvP,ituv,C
    INTEGER :: nTUVQ,nTUVTMPP,nTUVTMPQ,JPQ,JP,JQ,nTUVTMP,nTUVTMPprev
    integer :: tq,uq,vq,ituvq,ituvqminus1x,ituvqminus1y,ituvqminus1z
    integer :: iTUVplus1x,iTUVplus1y,iTUVplus1z,nTUVTMPPs
    Integer :: MaxAngmomQP,nTUVplus,JTMP,ntuvprev2,ntuvprev3
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

    WRITE(*,'(A)')'MODULE AGC_OBS_TRANSFERRECURRENCEMODDtoA'
    WRITE(*,'(A)')' use IchorPrecisionModule'
    WRITE(*,'(A)')'  '
    WRITE(*,'(A)')' CONTAINS'
    MaxAngmomQP = 8

    DO JMAX=2,MaxAngmomQP
       DO JP = 1, JMAX
          JQ = JMAX-JP
          IF(JP.GT.JQ)CYCLE
          IF(JP.EQ.0)CYCLE
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

          !The Equation D to A looks like
          !D to A
          !Theta(i+1,0,0,l) = -(b*X_{ab}-cX_{cd})/p Theta(i,0,0,l) - q/p*Theta(i,0,0,l+1) + i/(2p)*Theta(i-1,0,0,l) + l/(2p)*Theta(i,0,0,l-1)
          !C to A
          !Theta(i+1,0,k,0) = -(b*X_{ab}+dX_{cd})/p Theta(i,0,k,0) - q/p*Theta(i,0,k+1,0) + i/(2p)*Theta(i-1,0,k,0) + k/(2p)*Theta(i,0,k-1,0)
          !A to C
          !Theta(i,0,k+1,0) = -(b*X_{ab}+dX_{cd})/q Theta(i,0,k,0) - p/q*Theta(i+1,0,k,0) + k/(2q)*Theta(i,0,k-1,0) + i/(2q)*Theta(i-1,0,k,0) 
          !A to D
          !Theta(i,0,0,l+1) = -(b*X_{ab}-cX_{cd})/q Theta(i,0,0,l) - p/q*Theta(i+1,0,0,l) + l/(2q)*Theta(i,0,0,l-1) + i/(2q)*Theta(i-1,0,0,l) 

          !The Equation D to B looks like
          !D to B
          !Theta(0,j+1,0,l) = -(-a*X_{ab}-cX_{cd})/p Theta(0,j,0,l) - q/p*Theta(0,j,0,l+1) + i/(2p)*Theta(0,j-1,0,l) + l/(2p)*Theta(i,0,0,l-1)
          !C to B
          !Theta(i+1,0,k,0) = -(-a*X_{ab}+dX_{cd})/p Theta(0,j,k,0) - q/p*Theta(0,j,k+1,0) + i/(2p)*Theta(0,j-1,k,0) + k/(2p)*Theta(i,0,k-1,0)
          !B to C
          !Theta(i,0,k+1,0) = -(-a*X_{ab}+dX_{cd})/q Theta(0,j,k,0) - p/q*Theta(0,j+1,k,0) + k/(2q)*Theta(0,j,k-1,0) + i/(2q)*Theta(i-1,0,k,0) 
          !B to D
          !Theta(i,0,0,l+1) = -(-a*X_{ab}-cX_{cd})/q Theta(0,j,0,l) - p/q*Theta(0,j+1,0,l) + l/(2q)*Theta(0,j,0,l-1) + i/(2q)*Theta(i-1,0,0,l) 

          WRITE(*,'(A)')''
          IF(JP.LT.10)THEN
             IF(JQ.LT.10)THEN
                WRITE(*,'(A,I1,A,I1,A)')'subroutine TransferRecurrenceP',JP,'Q',JQ,'DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(*,'(A,I1,A,I2,A)')'subroutine TransferRecurrenceP',JP,'Q',JQ,'DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ELSE
             IF(JQ.LT.10)THEN
                WRITE(*,'(A,I2,A,I1,A)')'subroutine TransferRecurrenceP',JP,'Q',JQ,'DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(*,'(A,I2,A,I2,A)')'subroutine TransferRecurrenceP',JP,'Q',JQ,'DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ENDIF

          WRITE(*,'(A)')'         & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,Aux,Aux2)'
          WRITE(*,'(A)')'  implicit none'
          WRITE(*,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD'
          WRITE(*,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)'
          WRITE(*,'(A)')'  real(realk),intent(in) :: Pdistance12(3),Qdistance12(3,nPasses)'
          WRITE(*,'(A)')'  real(realk),intent(in) :: Aexp(nPrimA),Cexp(nPrimC)'
          WRITE(*,'(A,I5,A)')'  real(realk),intent(in) :: Aux(',nTUV,',nPrimQ*nPrimP*nPasses)'
          WRITE(*,'(A,I5,A,I5,A)')'  real(realk),intent(inout) :: Aux2(',nTUVP,',',nTUVQ,',nPrimQ*nPrimP*nPasses)'
          WRITE(*,'(A)')'!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ*nPrimP*nPasses)'
          WRITE(*,'(A)')'  !Local variables'
          !             WRITE(*,'(A)')'  real(realk) :: Pexpfac,PREF'
          !             WRITE(*,'(A,i4,A)')'  real(realk) :: TwoTerms(',MAX(1,nTUVprev2-nTUVprev3),')'
          DO JTMP=1,JP-1
             nTUVTMPPs=(JTMP)*(JTMP+1)*(JTMP+2)/6
             nTUVTMPP=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
             nTUVTMPQ=(JPQ-JTMP+1)*(JPQ-JTMP+2)*(JPQ-JTMP+3)/6
             if(JTMP.LT.10)THEN
                WRITE(*,'(A,I1,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVTMPPs+1,':',nTUVTMPP,',',nTUVQ+1,':',nTUVTMPQ,')'
             else
                WRITE(*,'(A,I2,A,I3,A,I3,A,I3,A,I3,A)')'  real(realk) :: Tmp',JTMP,'(',nTUVTMPPs+1,':',nTUVTMPP,',',nTUVQ+1,':',nTUVTMPQ,')'

             endif
          ENDDO
          WRITE(*,'(A)')'!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering'
          WRITE(*,'(A)')'  integer :: iPassQ,iPrimP,iPrimQ,IP,iTUVQ,iPrimA,iPrimC'
          WRITE(*,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk'
          WRITE(*,'(A)')'  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,expAX,expAY,expAZ'
          WRITE(*,'(A)')'  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp'
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
          WRITE(*,'(A)')'    invexpP = D1/expP'
!          WRITE(*,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'
          WRITE(*,'(A)')'    iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA'
          WRITE(*,'(A)')'    expAX = -Aexp(iPrimA)*Xab'
          WRITE(*,'(A)')'    expAY = -Aexp(iPrimA)*Yab'
          WRITE(*,'(A)')'    expAZ = -Aexp(iPrimA)*Zab'
          WRITE(*,'(A)')'    DO iPrimQ=1, nPrimQ'
!          WRITE(*,'(A)')'     iPrimD = (iPrimQ-1)/nPrimC+1'
          WRITE(*,'(A)')'     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC'
          WRITE(*,'(A)')'     IP = IP + 1'          
          WRITE(*,'(A)')'     qinvp = -Qexp(iPrimQ)*invexpP'
          WRITE(*,'(A)')'     inv2expP = D05*invexpP'
          WRITE(*,'(A)')'     facX = -(expAX-Cexp(iPrimC)*Xcd)*invexpP'
          WRITE(*,'(A)')'     facY = -(expAY-Cexp(iPrimC)*Ycd)*invexpP'
          WRITE(*,'(A)')'     facZ = -(expAZ-Cexp(iPrimC)*Zcd)*invexpP'


          allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
          CREATED  = .FALSE.
          CREATED(0,0,0) = .TRUE.

          DO JTMP=0,JP
           IF(JTMP.EQ.0)THEN
              WRITE(*,'(A,I3)')'     DO iTUVQ=1,',nTUVQ
              WRITE(*,'(A)')   '        Aux2(1,iTUVQ,IP) = Aux(iTUVQ,IP)'
              WRITE(*,'(A)')   '     ENDDO'
              CYCLE
           ELSE
             nTUVTMPP=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
             nTUVTMPQ=(JPQ-JTMP+1)*(JPQ-JTMP+2)*(JPQ-JTMP+3)/6
             !slowly increase JP=JTMP with 1 
             nTUVTMPP=(JTMP+1)*(JTMP+2)*(JTMP+3)/6
             DO Tp=JTMP,0,-1       
                DO Up=JTMP-Tp,0,-1
                   Vp=JTMP-Tp-Up                
!                   print*,'!Tq,Uq,Vq = ',Tq,Uq,Vq
                   !step 1.
                   !how can the (Tp,Up,Vp) be built from lower
                   TREC = CREATED(Tp-1,Up,Vp)
                   UREC = CREATED(Tp,Up-1,Vp)
                   VREC = CREATED(Tp,Up,Vp-1)

                   N=0
                   IF(TREC)N=N+1
                   IF(UREC)N=N+1
                   IF(VREC)N=N+1
                   IF(N.EQ.1)THEN
                      !only one possible way to construct it
                      IF(TREC)THEN
!                         print*,'!A TRECURRENCE iTUV',iTUV
                         call TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                      ELSEIF(UREC)THEN
!                         print*,'!A URECURRENCE iTUV',iTUV
                         call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                      ELSEIF(VREC)THEN
!                         print*,'!A VRECURRENCE iTUV',iTUV
                         call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                      ELSE
                         STOP 'TK1'
                      ENDIF
                   ELSE
                      !several ways to construct it
                      TREC2 = CREATED(Tp-2,Up,Vp)
                      UREC2 = CREATED(Tp,Up-2,Vp)
                      VREC2 = CREATED(Tp,Up,Vp-2)
                      N2=0
                      IF(TREC2)N2=N2+1
                      IF(UREC2)N2=N2+1
                      IF(VREC2)N2=N2+1
                      IF(N2.LT.N)THEN
                         !two term recurrence possible for one or more the possibilities 
                         !we chose the one term possibility
                         IF(.NOT.(TREC.AND.TREC2).AND.TREC)THEN
!                            print*,'!B TRECURRENCE iTUV',iTUV
                            call TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                         ELSEIF(.NOT.(UREC.AND.UREC2).AND.UREC)THEN
!                            print*,'!B URECURRENCE iTUV',iTUV
                            call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                         ELSEIF(.NOT.(VREC.AND.VREC2).AND.VREC)THEN
!                            print*,'!B VRECURRENCE iTUV',iTUV
                            call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                         ENDIF
                      ELSE
                         !chose one of the possibilities
                         IF(TREC)THEN
!                            print*,'!C TRECURRENCE iTUV',iTUV
                            call TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                         ELSEIF(UREC)THEN
!                            print*,'!C URECURRENCE iTUV',iTUV
                            call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                         ELSEIF(VREC)THEN
!                            print*,'!D VRECURRENCE iTUV',iTUV
                            call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
                         ELSE
                            STOP 'TK2'
                         ENDIF
                      ENDIF
                   ENDIF
                   CREATED(Tp,Up,Vp) = .TRUE.
                   
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
                WRITE(*,'(A,I1,A,I1,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'DtoA'
             ELSE
                WRITE(*,'(A,I1,A,I2,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'DtoA'
             ENDIF
          ELSE
             IF(JQ.LT.10)THEN
                WRITE(*,'(A,I2,A,I1,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'DtoA'
             ELSE
                WRITE(*,'(A,I2,A,I2,A)')'end subroutine TransferRecurrenceP',JP,'Q',JQ,'DtoA'
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
  

subroutine TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
implicit none
integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVPminus1x,nTUVQ
integer :: iTUVPminus2x, iTUVP,nTUVTMPQ,iTUVQminus1x,iTUVQplus1x,JMAX,JQ,TMQ,TMP,JP
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: TINDEX(:)
integer :: UINDEX(:)
integer :: VINDEX(:)
integer :: JINDEX(:)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
character(len=132) :: STRING 
integer :: iString
iTUVP = TUVINDEX(Tp,Up,Vp)
iTUVPminus1x = TUVINDEX(Tp-1,Up,Vp)
iTUVPminus2x = TUVINDEX(Tp-2,Up,Vp)
TMP = Tp-1
do iTUVQ = 1,nTUVTMPQ
   Tq = Tindex(iTUVq) 
   Uq = Uindex(iTUVq) 
   Vq = Vindex(iTUVq)
   TMQ = Tq
   iTUVQplus1x = TUVindex(Tq+1,Uq,Vq)
   iTUVQminus1x = TUVindex(Tq-1,Uq,Vq)
   Jq = Jindex(iTUVq)   

   CALL XYZRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ,&
        & 'facX',iTUVP,iTUVPminus1x,iTUVPminus2x,TMP,iTUVQ,Tq,Uq,Vq,TMQ,iTUVQplus1x,iTUVQminus1x,Jq)

ENDDO

end subroutine TRECURRENCE

SUBROUTINE XYZRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ,&
     & DIRECTIONSTRING,iTUVP,iTUVPminus1x,iTUVPminus2x,TMP,iTUVQ,Tq,Uq,Vq,TMQ,iTUVQplus1x,iTUVQminus1x,Jq)
  implicit none
  character(len=4) :: DIRECTIONSTRING
  integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVPminus1x
  integer :: iTUVPminus2x, iTUVP,nTUVTMPQ,iTUVQminus1x,iTUVQplus1x,JMAX,JQ,TMQ,TMP,JP,nTUVQ
  integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
  integer :: TINDEX(:)
  integer :: UINDEX(:)
  integer :: VINDEX(:)
  integer :: JINDEX(:)
  logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
  character(len=132) :: STRING 
  integer :: iString
  
  !D to A
  !Theta(i+1,0,0,l) = -(b*X_{ab}-cX_{cd})/p Theta(i,0,0,l) - q/p*Theta(i,0,0,l+1) + i/(2p)*Theta(i-1,0,0,l) + l/(2p)*Theta(i,0,0,l-1)

  !step 1 add blanks
!  print*,'!TEST1A1'
  STRING(1:5) = '     '
  iSTRING = 6
  !step 2 determine where to put the 
  IF(iTUVQ.LE.nTUVQ)THEN
     WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
     iString = iSTRING+5
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
!  print*,'!TEST1A2'
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
  
!  print*,'!TEST1A3'
  IF(iTUVQ.LE.nTUVQ)THEN
     WRITE(STRING(iSTRING:iSTRING+6),'(A7)') ',IP) = '
     iString = iSTRING+7
  ELSE
     WRITE(STRING(iSTRING:iSTRING+3),'(A4)') ') = '
     iString = iSTRING+4
  ENDIF
!  print*,'!TEST1A4'
  !step 3 determine if the first term: facX*THETA should be included 
  !and if it should use Aux or Tmp 
  
  IF(iTUVPminus1x.EQ.1)THEN
     !Aux(',iTUVp,',IP)
     WRITE(STRING(iSTRING:iSTRING+8),'(A4,A5)') DIRECTIONSTRING,'*Aux('
     iString = iSTRING+9
  ELSE
     IF(iTUVQ.LE.nTUVQ)THEN
        WRITE(STRING(iSTRING:iSTRING+9),'(A4,A6)') DIRECTIONSTRING,'*Aux2('
        iString = iSTRING+10
     ELSE
        IF(JTMP-1.LT.10)THEN
           WRITE(STRING(iSTRING:iSTRING+9),'(A4,A4,I1,A1)') DIRECTIONSTRING,'*Tmp',JTMP-1,'('
           iString = iSTRING+10
        ELSE
           WRITE(STRING(iSTRING:iSTRING+10),'(A4,A4,I2,A1)') DIRECTIONSTRING,'*Tmp',JTMP-1,'('
           iString = iSTRING+11
        ENDIF
     ENDIF
  ENDIF
!  print*,'!TEST1A5'
  IF(iTUVPminus1x.EQ.1)THEN
     IF(iTUVQ.LT.100)THEN
        WRITE(STRING(iSTRING:iSTRING+5),'(I2,A4)') iTUVQ,',IP)'
        iString = iSTRING+6
     ELSEIF(iTUVQ.LT.1000)THEN
        WRITE(STRING(iSTRING:iSTRING+6),'(I3,A4)') iTUVQ,',IP)'
        iString = iSTRING+7
     ELSEIF(iTUVQ.LT.10000)THEN
        WRITE(STRING(iSTRING:iSTRING+7),'(I4,A4)') iTUVQ,',IP)'
        iString = iSTRING+8
     ELSE
        STOP 'Recurrent B iTUVP'
     ENDIF
  ELSE
     IF(iTUVPminus1x.LT.100)THEN
        WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVPminus1x
        iString = iSTRING+3
     ELSEIF(iTUVPminus1x.LT.1000)THEN
        WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVPminus1x
        iString = iSTRING+4
     ELSE
        STOP 'Recurrent iTUVPminus'
     ENDIF
     IF(iTUVQ.LE.nTUVQ)THEN
        IF(iTUVQ.LT.100)THEN
           WRITE(STRING(iSTRING:iSTRING+6),'(A1,I2,A4)')',',iTUVQ,',IP)'
           iString = iSTRING+7
        ELSEIF(iTUVQ.LT.1000)THEN
           WRITE(STRING(iSTRING:iSTRING+7),'(A1,I3,A4)') ',',iTUVQ,',IP)'
           iString = iSTRING+8
        ELSEIF(iTUVQ.LT.10000)THEN
           WRITE(STRING(iSTRING:iSTRING+8),'(A1,I4,A4)') ',',iTUVQ,',IP)'
           iString = iSTRING+9
        ELSE
           STOP 'Recurrent B iTUVP'
        ENDIF
     ELSE
        IF(iTUVQ.LT.100)THEN
           WRITE(STRING(iSTRING:iSTRING+3),'(A1,I2,A1)')',',iTUVQ,')'
           iString = iSTRING+4
        ELSEIF(iTUVQ.LT.1000)THEN
           WRITE(STRING(iSTRING:iSTRING+4),'(A1,I3,A1)') ',',iTUVQ,')'
           iString = iSTRING+5
        ELSEIF(iTUVQ.LT.10000)THEN
           WRITE(STRING(iSTRING:iSTRING+5),'(A1,I4,A1)') ',',iTUVQ,')'
           iString = iSTRING+6
        ELSE
           STOP 'Recurrent B iTUVP'
        ENDIF
     ENDIF
  ENDIF
  !print*,'!TEST1A6'
  !step 4 determine if the second term: 
  !  TMP*inv2expQ*Aux(',iTUVminus1x,',IP) 
  !  TMp*inv2expQ*Tmp',JTMP-1,'(',iTUVminus1x,',',iTUVQminus1x,')
  !should be included and if it should use Aux or Tmp 
  IF(TMq.GT.2)THEN
     IF(TMq.GE.10)THEN
        WRITE(STRING(iSTRING:iSTRING+13),'(A2,I2,A10)') '+ ',TMQ,'*inv2expP*'
        iString = iSTRING+14
     ELSE
        WRITE(STRING(iSTRING:iSTRING+12),'(A2,I1,A10)') '+ ',TMQ,'*inv2expP*'
        iString = iSTRING+13
     ENDIF
  ELSEIF(TMq.EQ.2)THEN
     WRITE(STRING(iSTRING:iSTRING+9),'(A10)') '+ invexpP*'
     iString = iSTRING+10
  ELSEIF(TMq.EQ.1)THEN
     WRITE(STRING(iSTRING:iSTRING+10),'(A11)') '+ inv2expP*'
     iString = iSTRING+11
  ELSE
     !do not include this term 
  ENDIF
  
  !print*,'!TEST1A71'
  IF(TMq.GT.0)THEN
     IF(iTUVPminus1x.EQ.1)THEN
  !print*,'!TEST1A72'
        WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
        iString = iSTRING+4
     ELSE
        IF(iTUVQminus1x.LE.nTUVQ)THEN
  !print*,'!TEST1A73'
           WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
           iString = iSTRING+5
        ELSE
           IF(JTMP-1.LT.10)THEN
  !print*,'!TEST1A74'
              WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMP-1,'('
              iString = iSTRING+5
           ELSE
  !print*,'!TEST1A75'
              WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMP-1,'('
              iString = iSTRING+6
           ENDIF
        ENDIF
     ENDIF
     IF(iTUVPminus1x.LT.100)THEN
  !print*,'!TEST1A76'
        WRITE(STRING(iSTRING:iSTRING+2),'(I2)') iTUVPminus1x
        iString = iSTRING+3
     ELSEIF(iTUVPminus1x.LT.1000)THEN
  !print*,'!TEST1A77'
        WRITE(STRING(iSTRING:iSTRING+3),'(I3)') iTUVPminus1x
        iString = iSTRING+4
     ELSEIF(iTUVPminus1x.LT.10000)THEN
  !print*,'!TEST1A78'
        WRITE(STRING(iSTRING:iSTRING+4),'(I4)') iTUVPminus1x
        iString = iSTRING+5
     ELSE
        STOP 'Recurrent iTUVQminus1x'
     ENDIF
     IF(iTUVPminus1x.EQ.1)THEN
        WRITE(STRING(iSTRING:iSTRING+4),'(A5)')',IP) '
        iString = iSTRING+5
     ELSE
        IF(iTUVQminus1x.LE.nTUVQ)THEN
           IF(iTUVQminus1x.LT.100)THEN
              WRITE(STRING(iSTRING:iSTRING+6),'(A1,I2,A4)') ',',iTUVQminus1x,',IP)'
              iString = iSTRING+7
           ELSEIF(iTUVQminus1x.LT.1000)THEN
              WRITE(STRING(iSTRING:iSTRING+7),'(A1,I3,A4)') ',',iTUVQminus1x,',IP)'
              iString = iSTRING+8
           ELSEIF(iTUVQminus1x.LT.10000)THEN
              WRITE(STRING(iSTRING:iSTRING+8),'(A1,I4,A4)') ',',iTUVQminus1x,',IP)'
              iString = iSTRING+9
           ELSE
              STOP 'Recurrent iTUVminus1x'
           ENDIF
        ELSE
           IF(iTUVQminus1x.LT.100)THEN
              WRITE(STRING(iSTRING:iSTRING+3),'(A1,I2,A1)') ',',iTUVQminus1x,')'
              iString = iSTRING+4
           ELSEIF(iTUVQminus1x.LT.1000)THEN
              WRITE(STRING(iSTRING:iSTRING+4),'(A1,I3,A1)') ',',iTUVQminus1x,')'
              iString = iSTRING+5
           ELSEIF(iTUVQminus1x.LT.10000)THEN
              WRITE(STRING(iSTRING:iSTRING+5),'(A1,I4,A1)') ',',iTUVQminus1x,')'
              iString = iSTRING+6
           ELSE
              STOP 'Recurrent iTUVminus1x'
           ENDIF           
        ENDIF
     ENDIF
  ELSE
     !do not include this term       
  ENDIF
  !print*,'!TEST1A8'
  !step 5 determine if the third term: 
  !  Tmq*inv2expQ*Aux(',iTUVp,',IP)
  !  Tmq*inv2expQ*Tmp',JTMP-2,'(',iTUVp,',',iTUVQminus2x,')
  !should be included and if it should use Aux or Tmp 
  IF(TMp.GT.2)THEN
     IF(TMp.GE.10)THEN
        WRITE(STRING(iSTRING:iSTRING+13),'(A2,I2,A10)') '+ ',TMP,'*inv2expP*'
        iString = iSTRING+14
     ELSE
        WRITE(STRING(iSTRING:iSTRING+12),'(A2,I1,A10)') '+ ',TMP,'*inv2expP*'
        iString = iSTRING+13
     ENDIF
  ELSEIF(TMP.EQ.2)THEN
     WRITE(STRING(iSTRING:iSTRING+9),'(A10)') '+ invexpP*'
     iString = iSTRING+10
  ELSEIF(TMP.EQ.1)THEN
     WRITE(STRING(iSTRING:iSTRING+10),'(A11)') '+ inv2expP*'
     iString = iSTRING+11
  ELSE
     !do not include this term 
  ENDIF
  
  !print*,'!TEST1A9'
  IF(TMp.GT.0)THEN
     IF(iTUVPminus2x.EQ.1)THEN
        WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
        iString = iSTRING+4
     ELSE
        IF(iTUVQ.LE.nTUVQ)THEN
           WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
           iString = iSTRING+5
        ELSE
           IF(JTMP-1.LT.10)THEN
              WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMP-2,'('
              iString = iSTRING+5
           ELSE
              WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMP-2,'('
              iString = iSTRING+6
           ENDIF
        ENDIF
     ENDIF
     IF(iTUVPminus2x.LT.100)THEN
        WRITE(STRING(iSTRING:iSTRING+1),'(I2)') iTUVPminus2x
        iString = iSTRING+2
     ELSEIF(iTUVPminus2x.LT.1000)THEN
        WRITE(STRING(iSTRING:iSTRING+2),'(I3)') iTUVPminus2x
        iString = iSTRING+3
     ELSEIF(iTUVPminus2x.LT.10000)THEN
        WRITE(STRING(iSTRING:iSTRING+3),'(I4)') iTUVPminus2x
        iString = iSTRING+4
     ELSE
        STOP 'Recurrent iTUVPminus2x'
     ENDIF
     IF(iTUVPminus2x.EQ.1)THEN
        WRITE(STRING(iSTRING:iSTRING+4),'(A5)')',IP) '
        iString = iSTRING+5
     ELSE
        IF(iTUVQ.LE.nTUVQ)THEN
           IF(iTUVQ.LT.100)THEN
              WRITE(STRING(iSTRING:iSTRING+6),'(A1,I2,A4)') ',',iTUVQ,',IP)'
              iString = iSTRING+7
           ELSEIF(iTUVQ.LT.1000)THEN
              WRITE(STRING(iSTRING:iSTRING+7),'(A1,I3,A4)') ',',iTUVQ,',IP)'
              iString = iSTRING+8
           ELSEIF(iTUVQ.LT.10000)THEN
              WRITE(STRING(iSTRING:iSTRING+8),'(A1,I4,A4)') ',',iTUVQ,',IP)'
              iString = iSTRING+9
           ELSE
              STOP 'Recurrent B iTUVP'
           ENDIF
        ELSE
           IF(iTUVQ.LT.100)THEN
              WRITE(STRING(iSTRING:iSTRING+3),'(A1,I2,A1)') ',',iTUVQ,')'
              iString = iSTRING+4
           ELSEIF(iTUVQ.LT.1000)THEN
              WRITE(STRING(iSTRING:iSTRING+4),'(A1,I3,A1)') ',',iTUVQ,')'
              iString = iSTRING+5
           ELSEIF(iTUVQ.LT.10000)THEN
              WRITE(STRING(iSTRING:iSTRING+5),'(A1,I4,A1)') ',',iTUVQ,')'
              iString = iSTRING+6
           ELSE
              STOP 'Recurrent B iTUVP'
           ENDIF
        ENDIF
     ENDIF
  ELSE
     !do not include this term       
  ENDIF
  !print*,'!TEST1A10'
  !step 6 determine if the third term: 
  !  qinvp*Aux(',iTUVplus1x,',IP)'
  !  qinvp*Tmp',JTMP-1,'(',iTUVplus1x,',',iTUVQminus1x,')'
  !should be included and if it should use Aux or Tmp 
  !print*,'!TEST1'
  WRITE(STRING(iSTRING:iSTRING+7),'(A8)') '+ qinvp*'
  iString = iSTRING+8
  IF(iTUVPminus1x.EQ.1)THEN
  !print*,'!TEST2'
     WRITE(STRING(iSTRING:iSTRING+3),'(A4)') 'Aux('
     iString = iSTRING+4
  ELSE
     IF(iTUVQplus1x.LE.nTUVQ)THEN
  !print*,'!TEST3'
        WRITE(STRING(iSTRING:iSTRING+4),'(A5)') 'Aux2('
        iString = iSTRING+5
     ELSE
        IF(JTMP-1.LT.10)THEN
  !print*,'!TEST4'
           WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMP-1,'('
           iString = iSTRING+5
        ELSE
  !print*,'!TEST5'
           WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMP-1,'('
           iString = iSTRING+6
        ENDIF
     ENDIF
  ENDIF
  IF(iTUVPminus1x.NE.1)THEN
     IF(iTUVPminus1x.LT.100)THEN
        WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVPminus1x,','
        iString = iSTRING+3
     ELSEIF(iTUVPminus1x.LT.1000)THEN
        WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVPminus1x,','
        iString = iSTRING+4
     ELSEIF(iTUVPminus1x.LT.10000)THEN
        WRITE(STRING(iSTRING:iSTRING+4),'(I4,A1)') iTUVPminus1x,','
        iString = iSTRING+5
     ELSE
        STOP 'Recurrent iTUVQminus2x'
     ENDIF
  ENDIF
  IF(iTUVQplus1x.LT.100)THEN
     !print*,'!TEST6'
     WRITE(STRING(iSTRING:iSTRING+1),'(I2)') iTUVQplus1x
     iString = iSTRING+2
  ELSEIF(iTUVQplus1x.LT.1000)THEN
     !print*,'!TEST7'
     WRITE(STRING(iSTRING:iSTRING+2),'(I3)') iTUVQplus1x
     iString = iSTRING+3
  ELSEIF(iTUVQplus1x.LT.10000)THEN
     !print*,'!TEST8'
     WRITE(STRING(iSTRING:iSTRING+3),'(I4)') iTUVQplus1x
     iString = iSTRING+4
  ELSE
     STOP 'Recurrent iTUVplus1x'
  ENDIF
  IF(iTUVPminus1x.EQ.1)THEN
  !print*,'!TEST9'
     WRITE(STRING(iSTRING:iSTRING+3),'(A4)')',IP)'
     iString = iSTRING+4
  ELSE
     IF(iTUVQplus1x.LE.nTUVQ)THEN
        WRITE(STRING(iSTRING:iSTRING+4),'(A5)') ',IP) '
        iString = iSTRING+5
     ELSE
        WRITE(STRING(iSTRING:iSTRING+1),'(A2)') ') '
        iString = iSTRING+2
     ENDIF
  ENDIF
  !Final step write the string
  WRITE(*,'(A)') STRING(1:iSTRING-1)
end SUBROUTINE XYZRECURRENCE


subroutine URECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
implicit none
integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVPminus1x,nTUVQ
integer :: iTUVPminus2x, iTUVP,nTUVTMPQ,iTUVQminus1x,iTUVQplus1x,JMAX,JQ,TMQ,TMP,JP
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: TINDEX(:)
integer :: UINDEX(:)
integer :: VINDEX(:)
integer :: JINDEX(:)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
character(len=132) :: STRING 
integer :: iString
iTUVP = TUVINDEX(Tp,Up,Vp)
iTUVPminus1x = TUVINDEX(Tp,Up-1,Vp)
iTUVPminus2x = TUVINDEX(Tp,Up-2,Vp)
TMP = Up-1
do iTUVQ = 1,nTUVTMPQ
   Tq = Tindex(iTUVq) 
   Uq = Uindex(iTUVq) 
   Vq = Vindex(iTUVq)
   TMQ = Uq
   iTUVQplus1x = TUVindex(Tq,Uq+1,Vq)
   iTUVQminus1x = TUVindex(Tq,Uq-1,Vq)
   Jq = Jindex(iTUVq)   

   CALL XYZRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ,&
        & 'facY',iTUVP,iTUVPminus1x,iTUVPminus2x,TMP,iTUVQ,Tq,Uq,Vq,TMQ,iTUVQplus1x,iTUVQminus1x,Jq)

ENDDO

end subroutine URECURRENCE

subroutine VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ)
implicit none
integer :: Tq,Uq,Vq,Tp,Up,Vp,J,JTMP,iTUVQ,iTUVPminus1x,nTUVQ
integer :: iTUVPminus2x, iTUVP,nTUVTMPQ,iTUVQminus1x,iTUVQplus1x,JMAX,JQ,TMQ,TMP,JP
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: TINDEX(:)
integer :: UINDEX(:)
integer :: VINDEX(:)
integer :: JINDEX(:)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
character(len=132) :: STRING 
integer :: iString
iTUVP = TUVINDEX(Tp,Up,Vp)
iTUVPminus1x = TUVINDEX(Tp,Up,Vp-1)
iTUVPminus2x = TUVINDEX(Tp,Up,Vp-2)
TMP = Vp-1
do iTUVQ = 1,nTUVTMPQ
   Tq = Tindex(iTUVq) 
   Uq = Uindex(iTUVq) 
   Vq = Vindex(iTUVq)
   TMQ = Vq
   iTUVQplus1x = TUVindex(Tq,Uq,Vq+1)
   iTUVQminus1x = TUVindex(Tq,Uq,Vq-1)
   Jq = Jindex(iTUVq)   
   CALL XYZRECURRENCE(Tp,Up,Vp,J,TUVINDEX,TINDEX,UINDEX,VINDEX,JINDEX,CREATED,JTMP,nTUVTMPQ,JMAX,JP,nTUVQ,&
        & 'facZ',iTUVP,iTUVPminus1x,iTUVPminus2x,TMP,iTUVQ,Tq,Uq,Vq,TMQ,iTUVQplus1x,iTUVQminus1x,Jq)
ENDDO

end subroutine VRECURRENCE

end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
