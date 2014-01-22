MODULE TESTMODULE
  use TESTMODULEQXYZ
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

          WRITE(*,'(A)')'         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,Aux,Aux2)'
          WRITE(*,'(A)')'  implicit none'
          WRITE(*,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD'
          WRITE(*,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)'
          WRITE(*,'(A)')'  real(realk),intent(in) :: Pdistance12(3),Qdistance12(3,nPasses)'
          WRITE(*,'(A)')'  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)'
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
          WRITE(*,'(A)')'  integer :: iPassQ,iPrimP,iPrimQ,IP,iTUVQ,iPrimB,iPrimC'
          WRITE(*,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk'
          WRITE(*,'(A)')'  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,expBX,expBY,expBZ'
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
          WRITE(*,'(A)')'    iPrimB = (iPrimP-1)/nPrimA+1'
          WRITE(*,'(A)')'    expBX = Bexp(iPrimB)*Xab'
          WRITE(*,'(A)')'    expBY = Bexp(iPrimB)*Yab'
          WRITE(*,'(A)')'    expBZ = Bexp(iPrimB)*Zab'
          WRITE(*,'(A)')'    DO iPrimQ=1, nPrimQ'
!          WRITE(*,'(A)')'     iPrimD = (iPrimQ-1)/nPrimC+1'
          WRITE(*,'(A)')'     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC'
          WRITE(*,'(A)')'     IP = IP + 1'          
          WRITE(*,'(A)')'     qinvp = -Qexp(iPrimQ)*invexpP'
          WRITE(*,'(A)')'     inv2expP = D05*invexpP'
          WRITE(*,'(A)')'     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpP'
          WRITE(*,'(A)')'     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpP'
          WRITE(*,'(A)')'     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpP'


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
