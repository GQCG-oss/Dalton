MODULE TESTMODULEQXYZ
CONTAINS
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
           IF(iTUVQminus1x.LT.100)THEN
              WRITE(STRING(iSTRING:iSTRING+5),'(I2,A4)') iTUVQminus1x,',IP)'
              iString = iSTRING+6
           ELSEIF(iTUVQminus1x.LT.1000)THEN
              WRITE(STRING(iSTRING:iSTRING+6),'(I3,A4)') iTUVQminus1x,',IP)'
              iString = iSTRING+7
           ELSEIF(iTUVQminus1x.LT.10000)THEN
              WRITE(STRING(iSTRING:iSTRING+7),'(I4,A4)') iTUVQminus1x,',IP)'
              iString = iSTRING+8
           ELSE
              STOP 'Recurrent iTUVminus1x'
           ENDIF
!!$        ELSE
!!$           IF(JTMP-1.LT.10)THEN
!!$              !print*,'!TEST1A74'
!!$              WRITE(STRING(iSTRING:iSTRING+4),'(A3,I1,A1)') 'Tmp',JTMP-1,'('
!!$              iString = iSTRING+5
!!$           ELSE
!!$              !print*,'!TEST1A75'
!!$              WRITE(STRING(iSTRING:iSTRING+5),'(A3,I2,A1)') 'Tmp',JTMP-1,'('
!!$              iString = iSTRING+6
!!$           ENDIF
!!$           IF(iTUVQminus1x.LT.100)THEN
!!$              WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVQminus1x,')'
!!$              iString = iSTRING+3
!!$           ELSEIF(iTUVQminus1x.LT.1000)THEN
!!$              WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVQminus1x,')'
!!$              iString = iSTRING+4
!!$           ELSEIF(iTUVQminus1x.LT.10000)THEN
!!$              WRITE(STRING(iSTRING:iSTRING+4),'(I4,A1)') iTUVQminus1x,')'
!!$              iString = iSTRING+5
!!$           ELSE
!!$              STOP 'Recurrent iTUVminus1x'
!!$           ENDIF
!!$        ENDIF 
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
        IF(iTUVPminus1x.LT.100)THEN
           !print*,'!TEST1A76'
           WRITE(STRING(iSTRING:iSTRING+1),'(I2)') iTUVPminus1x
           iString = iSTRING+2
        ELSEIF(iTUVPminus1x.LT.1000)THEN
           !print*,'!TEST1A77'
           WRITE(STRING(iSTRING:iSTRING+2),'(I3)') iTUVPminus1x
           iString = iSTRING+3
        ELSEIF(iTUVPminus1x.LT.10000)THEN
           !print*,'!TEST1A78'
           WRITE(STRING(iSTRING:iSTRING+3),'(I4)') iTUVPminus1x
           iString = iSTRING+4
        ELSE
           STOP 'Recurrent iTUVQminus1x'
        ENDIF
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
        IF(iTUVPminus2x.LT.100)THEN
           WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVPminus2x,','
           iString = iSTRING+3
        ELSEIF(iTUVPminus2x.LT.1000)THEN
           WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVPminus2x,','
           iString = iSTRING+4
        ELSEIF(iTUVPminus2x.LT.10000)THEN
           WRITE(STRING(iSTRING:iSTRING+4),'(I4,A1)') iTUVPminus2x,','
           iString = iSTRING+5
        ELSE
           STOP 'Recurrent iTUVPminus2x'
        ENDIF
     ENDIF

     IF(iTUVPminus2x.NE.1.AND.iTUVQ.GT.nTUVQ)THEN !tmp
        IF(iTUVQ.LT.100)THEN
           WRITE(STRING(iSTRING:iSTRING+2),'(I2,A1)') iTUVQ,')'
           iString = iSTRING+3
        ELSEIF(iTUVQ.LT.1000)THEN
           WRITE(STRING(iSTRING:iSTRING+3),'(I3,A1)') iTUVQ,')'
           iString = iSTRING+4
        ELSEIF(iTUVQ.LT.10000)THEN
           WRITE(STRING(iSTRING:iSTRING+4),'(I4,A1)') iTUVQ,')'
           iString = iSTRING+5
        ELSE
           STOP 'Recurrent B iTUVP'
        ENDIF
     ELSE !aux or Aux2
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

end MODULE TESTMODULEQXYZ
