PROGRAM TUV
  use math
  use stringsMODULE
  implicit none
  logical :: nPrimLast
  integer,pointer :: TUVINDEX(:,:,:),TUVINDEXP(:,:,:)
  integer :: JMAX,J,JMAX1,JMAXP
  logical,pointer :: Enoscreen(:,:),EnoscreenS(:,:),zero(:)
  integer :: ijk1,ijk2,ijkcart,ijk,ijkcart1,ijkcart2,nTUV,ijkP
  integer :: iTUV,ilmP
  real(realk),pointer :: SCMAT1(:,:),SCMAT2(:,:),Spherical(:,:)
  logical :: sph1,sph2,sphericalTrans,sphericalGTO,newparam,ELSESTATEMENT
  integer :: LUMAIN,LUMOD1,LUMOD2,LUMOD3,LUMOD4,iparam,iparam2,nparam
  integer :: JMAX2,JMAX3,JMAX4,ituvP,jp,tp,up,vp,ntuvP,l1,l2,l12
  integer :: ip1,jp1,kp1,p1,ip2,jp2,kp2,p2,ijkcartP
  real(realk),pointer :: uniqeparam(:)
  character(len=15),pointer :: uniqeparamNAME(:)
    character(len=3) :: ARCSTRING
    integer :: GPUrun,nString
    logical :: DoOpenMP,DoOpenACC,CPU
  !buildtuvindex
  sphericalGTO = .TRUE.
  LUMOD3=3
  ! 0 1 2 3 4
  ! S P D F 
  !
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

  open(unit = LUMOD3, file="AutoGenCoderunSphContractOBS1_"//ARCSTRING//"_new.F90",status="unknown")
  WRITE(LUMOD3,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_Sphcontract1Mod'
  WRITE(LUMOD3,'(A)')'!Automatic Generated Code (AGC) by runSphContractOBS1.f90 in tools directory'
  WRITE(LUMOD3,'(A)')'use IchorPrecisionMod'

  WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')' CONTAINS'
  WRITE(LUMOD3,'(A)')'  '
!  IF(GPUrun.EQ.2)WRITE(LUMOD3,'(A)')'#ifdef VAR_OPENACC'

  JMAX1=2
  JMAX2=2
  do JMAXP = 0,JMAX1+JMAX2
     allocate(TUVINDEXP(0:JMAXP,0:JMAXP,0:JMAXP))
     ituvP = 0 
     DO JP = 0, JMAXP
        DO Tp=JP,0,-1       
           DO Up=JP-Tp,0,-1
              Vp=JP-Tp-Up
              ituvP = ituvP + 1 
              TUVINDEXP(Tp,Up,Vp) = ituvP
           ENDDO
        ENDDO
     ENDDO
     ntuvP = ituvP
     ntuv = ntuvp
     IF(nTUVP.NE.(JMAXP+1)*(JMAXP+2)*(JMAXP+3)/6)STOP 'ERROR nTUVQ'

     do l1 = MIN(JMAXP,JMAX1),0,-1 
        l2 = JMAXP-l1          !            J2=3
        IF(l2.GT.JMAX2)CYCLE
        l12 = l1 + l2

        ijkcart1 = (l1 + 1)*(l1 + 2)/2
        ijkcart2 = (l2 + 1)*(l2 + 2)/2
        ijkcart = ijkcart1*ijkcart2
        ijkcartP = ijkcart
        ijk1 = 2*l1 + 1
        ijk2 = 2*l2 + 1
        ijk = ijk1*ijk2
        ijkP = ijk

        Sph1 = sphericalGTO.AND.(l1.GT. 1)
        Sph2 = sphericalGTO.AND.(l2.GT. 1)
        SphericalTrans = Sph1.OR.Sph2
        !spherical Transform
        !ijkcart, ijk depend on l1,l2 SPHMAT depend on it
        !construct SPHMAT
        allocate(SCMAT1(ijk1,ijkcart1))
        allocate(SCMAT2(ijk2,ijkcart2))
        allocate(Spherical(ijkcart,ijk))
        !CSMAT1(ijkcart1,ijk1)
        call Sph_to_Cart_matrix(L1,SCMAT1,ijk1,ijkcart1,6,0)
        !CSMAT2(ijkcart2,ijk2)
        call Sph_to_Cart_matrix(L2,SCMAT2,ijk2,ijkcart2,6,0)
        call Buildsphericaltransformation(Spherical,SCMAT1,SCMAT2,&
             & ijk1,ijk2,ijkcart1,ijkcart2)
        !print*,'spherical'
        !call output(Spherical,1,ijkcart,1,ijk,ijkcart,ijk,1,6)
!        do ilmP=1,ijk
!           do ijkP = 1,ijkcart
!              WRITE(LUMOD3,'(A,i3,A,i3,A,ES16.8)')&
!                   &'!    SPHERICALMAT(',ijkP,',',ilmP,') = ',Spherical(ijkP,ilmP)
!           enddo
!        enddo

        IF(SphericalTrans.AND.(l12.LT.5.OR.(l12.EQ.5.AND.l1.EQ.3)) )THEN
           call initString(1)
           call AddToString('subroutine SphericalContractOBS1_')
           call AddToString(ARCSTRING)
           call AddToString('_maxAngP')
           call AddToString(l1+l2)
           call AddToString('_maxAngA')
           call AddToString(l1)
           call AddToString('(ijkQcart,nContPasses,IN,OUT')
           IF(DoOpenACC)THEN
              call AddToString(',iASync)')
           ELSE
              call AddToString(')')
           ENDIF
           call writeString(LUMOD3)

           WRITE(LUMOD3,'(A)')'  implicit none'
           WRITE(LUMOD3,'(A)')'  integer,intent(in)        :: ijkQcart,nContPasses'
           IF(nPrimLast)THEN
              WRITE(LUMOD3,'(A,I3,A)')'  real(realk),intent(in)    :: IN(',ijkcartP,',ijkQcart*nContPasses)'
              WRITE(LUMOD3,'(A,I3,A)')'  real(realk),intent(inout) :: OUT(',ijkP,',ijkQcart*nContPasses)'
           ELSE
              WRITE(LUMOD3,'(A,I3,A)')'  real(realk),intent(in)    :: IN(nContPasses,',ijkcartP,',ijkQcart)'
              WRITE(LUMOD3,'(A,I3,A)')'  real(realk),intent(inout) :: OUT(nContPasses,',ijkP,',ijkQcart)'
           ENDIF
           IF(DoOpenACC)WRITE(LUMOD3,'(A)')'  integer(kind=acckind),intent(in) :: iASync'
           WRITE(LUMOD3,'(A)')'  integer :: iP,ijkQ'
           
              iparam = 0
              do ijkP = 1,ijkcart
                 do ilmP=1,ijk
                    IF(ABS(Spherical(ijkP,ilmP)).GT.1.0E-8_realk)THEN
                       iparam = iparam + 1
                    ENDIF
                 enddo
              enddo
              allocate(uniqeparam(iparam))
              allocate(uniqeparamNAME(iparam))
              iparam = 0
              do ijkP = 1,ijkcart
                 do ilmP=1,ijk
                    IF(ABS(Spherical(ijkP,ilmP)).GT.1.0E-8_realk)THEN
                       IF(iparam.EQ.0)THEN
                          iparam = 1
                          uniqeparam(1) = Spherical(ijkP,ilmP)
                          IF(ijkP.LT.10.AND.ilmP.LT.10)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i1,A1,i1,A6)')'SPHMAT',ijkP,'_',ilmP,'      ' !15
                          ELSEIF(ijkP.LT.100.AND.ilmP.LT.10)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i2,A1,i1,A5)')'SPHMAT',ijkP,'_',ilmP,'     '               
                          ELSEIF(ijkP.LT.1000.AND.ilmP.LT.10)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i3,A1,i1,A4)')'SPHMAT',ijkP,'_',ilmP,'    '               
                          ELSEIF(ijkP.LT.10000.AND.ilmP.LT.10)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i4,A1,i1,A3)')'SPHMAT',ijkP,'_',ilmP,'   '               
                          ELSEIF(ijkP.LT.10.AND.ilmP.LT.100)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i1,A1,i2,A5)')'SPHMAT',ijkP,'_',ilmP,'     ' !15
                          ELSEIF(ijkP.LT.100.AND.ilmP.LT.100)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i2,A1,i2,A4)')'SPHMAT',ijkP,'_',ilmP,'    '               
                          ELSEIF(ijkP.LT.1000.AND.ilmP.LT.100)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i3,A1,i2,A3)')'SPHMAT',ijkP,'_',ilmP,'   '               
                          ELSEIF(ijkP.LT.10000.AND.ilmP.LT.100)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i4,A1,i2,A2)')'SPHMAT',ijkP,'_',ilmP,'  '               
                          ELSEIF(ijkP.LT.10.AND.ilmP.LT.1000)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i1,A1,i3,A4)')'SPHMAT',ijkP,'_',ilmP,'    ' !15
                          ELSEIF(ijkP.LT.100.AND.ilmP.LT.1000)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i2,A1,i3,A3)')'SPHMAT',ijkP,'_',ilmP,'   '               
                          ELSEIF(ijkP.LT.1000.AND.ilmP.LT.1000)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i3,A1,i3,A2)')'SPHMAT',ijkP,'_',ilmP,'  '               
                          ELSEIF(ijkP.LT.10000.AND.ilmP.LT.1000)THEN
                             WRITE(uniqeparamNAME(1),'(A6,i4,A1,i3,A1)')'SPHMAT',ijkP,'_',ilmP,' '               
                          ELSE
                             STOP 'ERROR '
                          ENDIF
                       ELSE
                          newparam = .TRUE.
                          do iparam2 = 1,iparam
                             IF(ABS(Spherical(ijkP,ilmP)-uniqeparam(iparam2)).LT.1.0E-13_realk)THEN
                                newparam = .FALSE.
                             ENDIF
                          enddo
                          IF(newparam)THEN
                             iparam = iparam + 1
                             uniqeparam(iparam) = Spherical(ijkP,ilmP)
                             IF(ijkP.LT.10.AND.ilmP.LT.10)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i1,A1,i1,A6)')'SPHMAT',ijkP,'_',ilmP,'      ' !15
                             ELSEIF(ijkP.LT.100.AND.ilmP.LT.10)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i2,A1,i1,A5)')'SPHMAT',ijkP,'_',ilmP,'     '               
                             ELSEIF(ijkP.LT.1000.AND.ilmP.LT.10)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i3,A1,i1,A4)')'SPHMAT',ijkP,'_',ilmP,'    '               
                             ELSEIF(ijkP.LT.10000.AND.ilmP.LT.10)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i4,A1,i1,A3)')'SPHMAT',ijkP,'_',ilmP,'   '               
                             ELSEIF(ijkP.LT.10.AND.ilmP.LT.100)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i1,A1,i2,A5)')'SPHMAT',ijkP,'_',ilmP,'     ' !15
                             ELSEIF(ijkP.LT.100.AND.ilmP.LT.100)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i2,A1,i2,A4)')'SPHMAT',ijkP,'_',ilmP,'    '               
                             ELSEIF(ijkP.LT.1000.AND.ilmP.LT.100)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i3,A1,i2,A3)')'SPHMAT',ijkP,'_',ilmP,'   '               
                             ELSEIF(ijkP.LT.10000.AND.ilmP.LT.100)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i4,A1,i2,A2)')'SPHMAT',ijkP,'_',ilmP,'  '               
                             ELSEIF(ijkP.LT.10.AND.ilmP.LT.1000)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i1,A1,i3,A4)')'SPHMAT',ijkP,'_',ilmP,'    ' !15
                             ELSEIF(ijkP.LT.100.AND.ilmP.LT.1000)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i2,A1,i3,A3)')'SPHMAT',ijkP,'_',ilmP,'   '               
                             ELSEIF(ijkP.LT.1000.AND.ilmP.LT.1000)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i3,A1,i3,A2)')'SPHMAT',ijkP,'_',ilmP,'  '               
                             ELSEIF(ijkP.LT.10000.AND.ilmP.LT.1000)THEN
                                WRITE(uniqeparamNAME(iparam),'(A6,i4,A1,i3,A1)')'SPHMAT',ijkP,'_',ilmP,' '               
                             ELSE
                                STOP 'ERROR '
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDIF
                 enddo
              enddo
              do iparam2 = 1,iparam
                 WRITE(LUMOD3,'(A,A15,A,ES26.16,A)')'  real(realk),parameter :: ',uniqeparamNAME(iparam2),'=',uniqeparam(iparam2),'_realk'
              enddo
              nparam = iparam

              allocate(zero(ijk))
              zero = .TRUE.

              !              WRITE(LUMOD3,'(A,I3,A,I3,A)')'  real(realk),intent(in)    :: IN(',ijkcartP,',ijkQcart,nPasses)'
              !              WRITE(LUMOD3,'(A,I3,A,I3,A)')'  real(realk),intent(inout) :: OUT(',ijkP,',ijkQcart,nPasses)'

!              IF(DoOpenMP) WRITE(LUMOD3,'(A)')'!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iP) SHARED(nContPasses,ijkQcart,IN,OUT)'
              IF(nPrimLast)THEN
                 IF(DoOpenMP) WRITE(LUMOD3,'(A)')'!$OMP DO PRIVATE(iP)'
                 IF(DoOpenACC) WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(IN,OUT) ASYNC(iASync)'
              ELSE
                 IF(DoOpenMP) WRITE(LUMOD3,'(A)')'!$OMP DO PRIVATE(iP,ijkQ)'
                 IF(DoOpenACC) WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP PRIVATE(iP,ijkQ) PRESENT(IN,OUT) ASYNC(iASync)'
              ENDIF
              IF(nPrimLast)THEN
                 WRITE(LUMOD3,'(A)')'  DO iP=1,ijkQcart*nContPasses'
              ELSE
                 WRITE(LUMOD3,'(A)')'  DO iP=1,nContPasses'
                 WRITE(LUMOD3,'(A)')'   DO ijkQ=1,ijkQcart'
              ENDIF
              do ilmP=1,ijk
                 call initString(4)
                 nString = 0 
                 do ijkP = 1,ijkcart
                    IF(ABS(Spherical(ijkP,ilmP)).GT.1.0E-8)THEN
                       IF(zero(ilmP))THEN 
                          IF(ABS(Spherical(ijkP,ilmP)-1.0E0_realk).GT.1.0E-10)THEN
                             iparam = -1
                             do iparam2 = 1,nparam
                                IF(ABS(Spherical(ijkP,ilmP)-uniqeparam(iparam2)).LT.1.0E-13_realk)THEN
                                   iparam = iparam2
                                ENDIF
                             enddo
                             IF(nPrimLast)THEN
                                call AddToString('OUT(')
                                call AddToString(ilmP)
                                call AddToString(',iP) = IN(')
                                call AddToString(ijkP)
                                call AddToString(',iP)*')
                             ELSE
                                call AddToString('OUT(iP,')
                                call AddToString(ilmP)
                                call AddToString(',ijkQ) = IN(iP,')
                                call AddToString(ijkP)
                                call AddToString(',ijkQ)*')
                             ENDIF
                             call AddToString(TRIM(uniqeparamNAME(iparam)))
                             nString = 4 + 4 + 3 + 10 + 3 + 5 + 12 

!                             WRITE(LUMOD3,'(A,i3,A,i3,A,A15)')&
!                                  &'    OUT(',ilmP,',iP) = IN(',ijkP,',iP)*',uniqeparamNAME(iparam)
                          ELSE
                             IF(nPrimLast)THEN
                                call AddToString('OUT(')
                                call AddToString(ilmP)
                                call AddToString(',iP) = IN(')
                                call AddToString(ijkP)
                                call AddToString(',iP)')
                             ELSE
                                call AddToString('OUT(iP,')
                                call AddToString(ilmP)
                                call AddToString(',ijkQ) = IN(iP,')
                                call AddToString(ijkP)
                                call AddToString(',ijkQ)')
                             ENDIF
                             nString = 4 + 4 + 3 + 10 + 3 + 4
!                             WRITE(LUMOD3,'(A,i3,A,i3,A)')&
!                                  &'    OUT(',ilmP,',iP) = IN(',ijkP,',iP)'
                          ENDIF
                          zero(ilmP) = .FALSE.
                       ELSE
                          IF(ABS(Spherical(ijkP,ilmP)-1.0E0_realk).GT.1.0E-10)THEN
                             IF(nString.GT.104)THEN
                                call AddToString(' &')
                                call writeString(LUMOD3)
                                call initString(15)
                                call AddToString('&')
                                nString = 16
                             ENDIF
                             iparam = 0
                             do iparam2 = 1,nparam
                                IF(ABS(Spherical(ijkP,ilmP)-uniqeparam(iparam2)).LT.1.0E-13_realk)THEN
                                   iparam = iparam2
                                ENDIF
                             enddo
                             IF(nPrimLast)THEN
                                call AddToString(' + IN(')
                                call AddToString(ijkP)
                                call AddToString(',iP)*')
                             ELSE
                                call AddToString(' + IN(iP,')
                                call AddToString(ijkP)
                                call AddToString(',ijkQ)*')
                             ENDIF
                             call AddToString(TRIM(uniqeparamNAME(iparam)))
                             nString = nString + 6 + 3 + 5 + 12

!                             WRITE(LUMOD3,'(A,i3,A,i3,A,i3,A,A15)')&
!                                  &'    OUT(',ilmP,',iP) = OUT(',ilmP,',iP) + IN(',ijkP,',iP)*',uniqeparamNAME(iparam)
                          ELSE
                             IF(nString.GT.118)THEN
                                call AddToString(' &')
                                call writeString(LUMOD3)
                                call initString(15)
                                call AddToString('&')
                                nString = 16
                             ENDIF
                             IF(nPrimLast)THEN
                                call AddToString(' + IN(')
                                call AddToString(ijkP)
                                call AddToString(',iP)')
                             ELSE
                                call AddToString(' + IN(iP,')
                                call AddToString(ijkP)
                                call AddToString(',ijkQ)')
                             ENDIF
                             nString = nString + 6 + 3 + 5
!                             WRITE(LUMOD3,'(A,i3,A,i3,A,i3,A)')&
!                                  &'    OUT(',ilmP,',iP) = OUT(',ilmP,',iP) + IN(',ijkP,',iP)'
                          ENDIF
                       ENDIF
                    ENDIF
                 enddo
                 call writeString(LUMOD3)
              enddo
              IF(nPrimLast)THEN
                 WRITE(LUMOD3,'(A)')'  ENDDO'
              ELSE
                 WRITE(LUMOD3,'(A)')'   ENDDO'
                 WRITE(LUMOD3,'(A)')'  ENDDO'
              ENDIF
!              IF(DoOpenMP) WRITE(LUMOD3,'(A)')'!$OMP END PARALLEL DO'
              IF(DoOpenMP) WRITE(LUMOD3,'(A)')'!$OMP END DO'
              IF(l12.LT.10)THEN
                 WRITE(LUMOD3,'(A,I1,A,I1,A)')'end subroutine SphericalContractOBS1_'//ARCSTRING//'_maxAngP',l1+l2,'_maxAngA',l1,' '
              ELSE
                 WRITE(LUMOD3,'(A,I2,A,I1,A)')'end subroutine SphericalContractOBS1_'//ARCSTRING//'_maxAngP',l1+l2,'_maxAngA',l1,' '
              ENDIF
              WRITE(LUMOD3,'(A)')'  '
              WRITE(LUMOD3,'(A)')'  '
           !endif
        ENDIF !SphericalTrans
     enddo
  enddo
!    IF(GPUrun.EQ.2)WRITE(LUMOD3,'(A)')'#endif'
  WRITE(LUMOD3,'(A)')'END MODULE AGC_'//ARCSTRING//'_OBS_Sphcontract1Mod'
  close(unit = LUMOD3)
  enddo


END PROGRAM

!contractecoeff_gen
