!> @file
!> Contains the main Thermite integral drivers

!> \brief Main Thermite drivers for the calculation of integrals 
!> based on the McMurchie-Davidson scheme, using (Turbo-)hermite functions. (PCCP 2007 9, 4771) 
!> \author T. Kjaergaard and S. Reine
!> \date 2008 
MODULE IchorEriCoulombintegralMod
use precision
use ThermiteMem_module
private 
public :: ICI_SSSS
CONTAINS
  subroutine ICI_SSSS(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimPQ,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & PQorder,CDAB)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimPQ,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented,PQorder
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP),qcent(3*nPrimQ*MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimPQ),reducedExponents(nPrimPQ)
    real(realk),pointer :: squaredDistance(:)
    !
!    real(realk),target :: TMPWORK((2+2*nPasses)*nPrimP*nPrimQ)
    real(realk),pointer :: RJ000(:),RE(:)
    integer :: nPrimPassQ,i,j
    IF(.NOT.PQorder)THEN
       !ordering is (PrimQ,nPrimP)
       !contract P first as done in ContractEcoeffPQ_seg_SSSS
       call lsquit('PrimQ,nPrimP',-1)
    ENDIF

    nPrimPassQ = nPrimQ*nPasses 
!    call mem_alloc(reducedExponents,nPrimPQ)
!    call mem_alloc(integralPrefactor,nPrimPQ)
    !build reducedExponents(nPrimP,nPrimQ),integralPrefactor(nPrimP,nPrimQ)
!    call build_exp_redexp_intprefactor_SSSS(nPrimQ,nPrimP,pexp,qexp,&
!         & reducedExponents,integralPrefactor)
    !build squaredDistance(nPrimP,nPrimQ,nPasses)
    call mem_alloc(squaredDistance,nPasses*nPrimPQ)
    call build_squaredDistance_SSSS(nPrimPassQ,nPrimP,Qcent,Pcent,&
         & squaredDistance)
    !build RJ000(nPrimP,nPrimQ,nPasses)
    call mem_alloc(RJ000,nPasses*nPrimPQ)
    call buildRJ000_SSSS(nPasses,nPrimPQ,nTABFJW1,nTABFJW2,reducedExponents,squaredDistance,&
         & TABFJW,RJ000)
!    call mem_dealloc(reducedExponents)
    call mem_dealloc(squaredDistance)
    !scale RJ000
    call scal_RJ000_SSSS(nPrimPQ,nPasses,RJ000,integralPrefactor)
    IF (INTPRINT .GE. 10) THEN
       WRITE(lupri,*)'Output from R000'
       DO I=1,NPrimP*nPrimQ
          WRITE(LUPRI,'(2X,A6,I4,A2,ES16.8)')'RJ000(',I,')=',RJ000(I)
       ENDDO
    END IF
!    call mem_dealloc(integralPrefactor)
    !RE(nPrimP,nPasses)
    IF(.NOT.Qsegmented.AND..NOT.Psegmented)THEN
       call lsquit('Psegment error1',-1)
       IF(nPrimP*nContQ*nPasses.LE.nContP*nPrimQ*nPasses)THEN
          call ContractEcoeffPQ_genQ1_SSSS(nPrimP,nPrimQ,nPasses,RJ000,QpreExpFac)
!          call ContractBasisPQ_genQ1_SSSS(nPrimP,nPrimQ,nPasses,RJ000,QpreExpFac)
       ELSE
          call ContractEcoeffPQ_genP1_SSSS(nPrimP,nPrimQ,nPasses,RJ000,PpreExpFac)
!          call ContractBasisPQ_genP1_SSSS(nPrimP,nPrimQ,nPasses,RJ000,PpreExpFac)
       ENDIF

    ELSEIF(.NOT.Qsegmented.AND.Psegmented)THEN
       call lsquit('Psegment error2',-1)
       !P first
       !RJ000(nPrimP,nPrimQ,nPasses) => ER(nPrimQ,nPasses)
       !ER(nPrimQ,nPasses) => ERE(nPrimQ,nPasses)
       !ERE(nPrimQ,nPasses) => EREC(nContQ,nPasses)
    ELSEIF(Qsegmented.AND..NOT.Psegmented)THEN
       call lsquit('Psegment error3',-1)
       !Q first
       call mem_alloc(RE,nPrimP*nPasses)
       call ContractEcoeffQ_seg_SSSS(nPrimP,nPrimQ,nPasses,RJ000,QpreExpFac,RE)
       call mem_dealloc(RJ000)
       !RE(nPrimP,nPasses) => ERE(nPrimP,nPasses)
       !ERE(nPrimP,nPasses) => CERE(nContQ,nPasses)
    ELSE !both segmented
!       call mem_alloc(RE,nPrimP*nPasses)
!       call ContractEcoeffQ_seg_SSSS(nPrimP,nPrimQ,nPasses,RJ000,QpreExpFac,RE)
!       call mem_dealloc(RJ000)
!       call ContractEcoeffP_seg_SSSS(nPrimP,nPasses,RE,PpreExpFac,CDAB)
!       call mem_dealloc(RE)
       call ContractEcoeffPQ_seg_SSSS(nPrimP,nPrimQ,nPasses,MaxPasses,&
            & RJ000,PpreExpFac,QpreExpFac,CDAB)
    ENDIF

  end subroutine ICI_SSSS

  subroutine ContractEcoeffPQ_seg_SSSS(nPrimP,nPrimQ,nPasses,MaxPasses,&
       & RJ000,PpreExpFac,QpreExpFac,CDAB)
    implicit none
    integer,intent(in) :: nPrimP,nPrimQ,nPasses,MaxPasses
    real(realk),intent(in) :: RJ000(nPrimP*nPrimQ*nPasses),QpreExpFac(nPrimQ,MaxPasses)
    real(realk),intent(in) :: PpreExpFac(nPrimP)
    real(realk),intent(inout) :: CDAB(nPasses)
    !
    integer :: iPrimQ,iPrimP,offset,iPassQ,offset2
    real(realk) :: tmp,TMPQ
    DO iPassQ = 1,nPasses
       tmp = 0.0E0_realk
       offset2 = (iPassQ-1)*nPrimP*nPrimQ
       DO iPrimQ =1,nPrimQ
          TMPQ = QpreExpFac(iPrimQ,iPassQ)
          offset = (iPrimQ-1)*nPrimP+offset2
          DO iPrimP = 1, nPrimP
             tmp = tmp + PpreExpFac(iPrimP)*RJ000(iPrimP+offset)*TMPQ
          ENDDO
       ENDDO
       CDAB(iPassQ) = tmp
    ENDDO
  end subroutine ContractEcoeffPQ_seg_SSSS

  subroutine ContractEcoeffPQ_genQ1_SSSS(nPrimP,nPrimQ,nPasses,RJ000,QpreExpFac)
    implicit none
    integer,intent(in) :: nPrimP,nPrimQ,nPasses
    real(realk),intent(inout) :: RJ000(nPrimP,nPrimQ*nPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*nPasses)
    !
    integer :: iPrimQ,iPrimP
    DO iPrimQ = 1,nPasses*nPrimQ
       DO iPrimP = 1, nPrimP
          RJ000(iPrimP,iPrimQ) = RJ000(iPrimP,iPrimQ)*QpreExpFac(iPrimQ)
       ENDDO
    ENDDO
  end subroutine ContractEcoeffPQ_genQ1_SSSS

  subroutine ContractEcoeffPQ_genP1_SSSS(nPrimP,nPrimQ,nPasses,RJ000,PpreExpFac)
    implicit none
    integer,intent(in) :: nPrimP,nPrimQ,nPasses
    real(realk),intent(inout) :: RJ000(nPrimP,nPrimQ*nPasses)
    real(realk),intent(in) :: PpreExpFac(nPrimP)
    !
    integer :: iPrimQ,iPrimP
    DO iPrimQ =1,nPrimQ*nPasses
       DO iPrimP = 1, nPrimP
          RJ000(iPrimP,iPrimQ) = RJ000(iPrimP,iPrimQ)*PpreExpFac(iPrimP)
       ENDDO
    ENDDO
  end subroutine ContractEcoeffPQ_genP1_SSSS

!!$  subroutine build_exp_redexp_intprefactor_SSSS(nPrimQ,nPrimP,pexp,qexp,&
!!$       & reducedExponents,integralPrefactor) !exponents&
!!$    implicit none
!!$    integer,intent(in) :: nPrimQ,nPrimP
!!$    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
!!$    !      real(realk),intent(inout) :: exponents(nPrimP,nPrimQ)
!!$    real(realk),intent(inout) :: integralPrefactor(nPrimP,nPrimQ)
!!$    real(realk),intent(inout) :: reducedExponents(nPrimP,nPrimQ)
!!$    !
!!$    real(realk) :: p,q,p_q
!!$    Real(realk), parameter :: PIFAC = 34.986836655249725E0_realk !Two*PI**TwoHalf      
!!$    integer :: iPrimQ,iPrimP
!!$    DO iPrimQ = 1, nPrimQ
!!$       q  = qexp(iPrimQ)
!!$       DO iPrimP=1, nPrimP
!!$          p  = pexp(iPrimP)
!!$          p_q = p + q
!!$          !(nPrimP,nPrimQ)
!!$          !            exponents(iPrimP,iPrimQ) = p_q  
!!$          reducedExponents(iPrimP,iPrimQ) = p*q/p_q
!!$          integralPrefactor(iPrimP,iPrimQ) = PIFAC/(p*q*SQRT(p_q))
!!$       ENDDO
!!$    ENDDO
!!$  end subroutine build_exp_redexp_intprefactor_SSSS

  subroutine build_squaredDistance_SSSS(nPrimPassQ,nPrimP,Qcent,Pcent,squaredDistance)
    implicit none
    integer,intent(in) :: nPrimPassQ,nPrimP
    real(realk),intent(in) :: Qcent(3,nPrimPassQ),Pcent(3,nPrimP)
    real(realk),intent(inout) :: squaredDistance(nPrimP,nPrimPassQ)
    !
    integer :: iPrimPassQ,iPrimP
    real(realk) :: qx,qy,qz,pqx,pqy,pqz
    DO iPrimPassQ=1, nPrimPassQ
       qx = -Qcent(1,iPrimPassQ)
       qy = -Qcent(2,iPrimPassQ)
       qz = -Qcent(3,iPrimPassQ)
       DO iPrimP=1, nPrimP
          pqx = Pcent(1,iPrimP) + qx
          pqy = Pcent(2,iPrimP) + qy
          pqz = Pcent(3,iPrimP) + qz          
          !(nPrimP,nPrimQ,nPasses)
          squaredDistance(iPrimP,iPrimPassQ) = pqx*pqx+pqy*pqy+pqz*pqz
       ENDDO
    ENDDO
  end subroutine build_squaredDistance_SSSS

  subroutine buildRJ000_SSSS(nPasses,nPrimPQ,nTABFJW1,nTABFJW2,reducedExponents,squaredDistance,&
       & TABFJW,RJ000)
    implicit none
    integer,intent(in) :: nPasses,nPrimPQ,nTABFJW1,nTABFJW2
    real(realk),intent(in) :: reducedExponents(nPrimPQ),squaredDistance(nPrimPQ,nPasses)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    real(realk),intent(inout) :: RJ000(nPrimPQ,nPasses)
    !
    integer :: ipq,ipassq,ipnt
    Real(realk), parameter :: PI=3.14159265358979323846E0_realk,HALF=0.5E0_realk
    REAL(REALK), PARAMETER :: D2=2.E0_realk,D100=100E0_realk,D4=4.0E0_realk
    REAL(REALK), PARAMETER :: D36=36.E0_realk,D1=1E0_realk,D12 = 12E0_realk
    REAL(REALK), PARAMETER :: SMALL = 1E-15_realk,TENTH = 0.01E0_realk
    REAL(REALK), PARAMETER :: COEF2 = 0.5E0_realk, COEF3=-D1/6.0E0_realk
    REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
    REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
    REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
    REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
    REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
    REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2,PID4 = PI/D4, PID4I = D4/PI
    real(realk) :: WVAL,R,RWVAL,REXPW,GVAL,WDIFF,W2,W3
    DO iPassq=1, nPasses
       DO ipq=1, nPrimPQ
          !      (nPrimP,nPrimQ)      *     (nPrimP,nPrimQ,nPasses)
          WVAL = reducedExponents(ipq)*squaredDistance(ipq,iPassQ)
          IF (ABS(WVAL) .LT. SMALL) THEN
             RJ000(ipq,iPassQ) = D1
          ELSE IF (WVAL .LT. D12) THEN
             IPNT = NINT(D100*WVAL)
             WDIFF = WVAL - TENTH*IPNT
             W2    = WDIFF*WDIFF
             W3    = W2*WDIFF
             W2    = W2*COEF2
             W3    = W3*COEF3
             R = TABFJW(0,IPNT)
             R = R -TABFJW(1,IPNT)*WDIFF
             R = R + TABFJW(2,IPNT)*W2
             R = R + TABFJW(3,IPNT)*W3
             RJ000(ipq,iPassQ) = R
             !  12 < WVAL <= (2J+36) 
          ELSE IF (WVAL.LE.D36) THEN
             REXPW = HALF*EXP(-WVAL)
             RWVAL = D1/WVAL
             GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
             RJ000(ipq,iPassQ) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
             !  (2J+36) < WVAL 
          ELSE
             RWVAL = PID4/WVAL
             RJ000(ipq,iPassQ) = SQRT(RWVAL)
             RWVAL = RWVAL*PID4I
          END IF
       ENDDO
    ENDDO
  end subroutine buildRJ000_SSSS

  subroutine scal_RJ000_SSSS(nPrimPQ,nPasses,RJ000,integralPrefactor)
    implicit none
    integer :: nPrimPQ,nPasses
    real(realk) :: RJ000(nPrimPQ,nPasses),integralPrefactor(nPrimPQ)
    !
    integer :: iPassq,ipq
    ! Scaling (nPrimP,nPrimQ,nPasses)
    DO iPassq=1, nPasses
       DO ipq=1, nPrimPQ
          RJ000(ipq,iPassQ) = integralPrefactor(Ipq)*RJ000(ipq,iPassQ)
       ENDDO
    ENDDO
  end subroutine scal_RJ000_SSSS

  subroutine build_pref_full_SSSS(nPrimQ,nPrimC,nPrimD,Dexp,Cexp,QpreExpFac,pref)
    implicit none
    integer     :: nPrimC,nPrimD,nPrimQ
    real(realk),intent(in) :: Dexp(nPrimD),Cexp(nPrimC),QpreExpFac(nPrimC,nPrimD)
    real(realk),intent(inout) :: pref(nPrimC,nPrimD)
    !
    integer :: iPrimD,iPrimC
    real(realk) :: DexpTmp
    REAL(REALK), PARAMETER :: D4=4.0E0_realk      
    DO iPrimD=1,nPrimD
       DexpTmp = D4*Dexp(iPrimD)
       DO iPrimC=1,nPrimC
          pref(iPrimC,iPrimD) = QpreExpFac(iPrimC,iPrimD)*DexpTmp*Cexp(iPrimC)
       ENDDO
    ENDDO
  end subroutine build_pref_full_SSSS

  subroutine build_pref_screen_SSSS(nPrimQ,nPrimC,nPrimD,Dexp,Cexp,Qiprim1,Qiprim2,QpreExpFac,pref)
    implicit none
    integer :: nPrimQ,nPrimC,nPrimD
    integer :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)
    real(realk) :: Dexp(nPrimD),Cexp(nPrimC),pref(nPrimQ),QpreExpFac(nPrimQ)
    !
    integer :: iPrimD,iPrimC,iPrimQ
    REAL(REALK), PARAMETER :: D4=4.0E0_realk      
    !primitive screening 
    DO iPrimQ=1,nPrimQ
       iPrimC = Qiprim1(iPrimQ)
       iPrimD = Qiprim2(iPrimQ)
       pref(iPrimQ) = QpreExpFac(iPrimQ)*D4*Cexp(iPrimC)*Dexp(iPrimD)
    ENDDO
  end subroutine build_pref_screen_SSSS

  subroutine ContractEcoeffQ_seg_SSSS(nPrimP,nPrimQ,nPasses,RJ000,pref,RE)
    implicit none
    integer,intent(in) :: nPrimP,nPrimQ,nPasses
    real(realk),intent(in) :: RJ000(nPrimP,nPrimQ,nPasses),pref(nPrimQ,nPasses)
    real(realk),intent(inout) :: RE(nPrimP,nPasses)
    !
    integer :: iPrimP,iPassQ,iPrimQ
    real(realk) :: tmp

    DO iPassQ =1,nPasses
       DO iPrimP = 1, nPrimP
          tmp = 0.0E0_realk
          DO iPrimQ =1,nPrimQ
             tmp = tmp + RJ000(iPrimP,iPrimQ,iPassQ)*pref(iPrimQ,iPassQ)
          ENDDO
          RE(iPrimP,iPassQ) = tmp
       ENDDO
    ENDDO
  end subroutine ContractEcoeffQ_seg_SSSS

  subroutine ContractEcoeffQ_gen_SSSS(nPrimP,nPrimQ,nPasses,RJ000,pref,RE)
    implicit none
    integer,intent(in) :: nPrimP,nPrimQ,nPasses
    real(realk),intent(in) :: RJ000(nPrimP,nPrimQ,nPasses),pref(nPrimQ)
    real(realk),intent(inout) :: RE(nPrimP,nPrimQ,nPasses)
    !
    integer :: iPrimP,iPassQ,iPrimQ
    real(realk) :: tmp

    DO iPassQ =1,nPasses
       DO iPrimP = 1, nPrimP
          DO iPrimQ =1,nPrimQ
             RE(iPrimP,iPrimQ,iPassQ) = RJ000(iPrimP,iPrimQ,iPassQ)*pref(iPrimQ)
          ENDDO
       ENDDO
    ENDDO
  end subroutine ContractEcoeffQ_gen_SSSS

  subroutine ContractEcoeffP_seg_SSSS(nPrimP,nPasses,RE,pref,CDAB)
    implicit none
    integer,intent(in) :: nPrimP,nPasses
    real(realk),intent(in) :: RE(nPrimP,nPasses),pref(nPrimP)
    real(realk) :: CDAB(nPasses)
    !
    integer :: iPassQ,iPrimP
    real(realk) :: tmp
    DO iPassQ=1,nPasses
       tmp = 0.0E0_realk
       DO iPrimP = 1, nPrimP
          tmp = tmp + RE(iPrimP,iPassQ)*pref(iPrimP)
       ENDDO
       CDAB(iPassQ) = tmp
    ENDDO
  end subroutine ContractEcoeffP_seg_SSSS

  SUBROUTINE ConstructContraction_screen_ssss(CC,nPrim,nCont1,nCont2,&
       & nPrim1,nPrim2,CC1,CC2,Piprim1,Piprim2,LUPRI,IPRINT)
    Implicit none
    Integer       :: nPrim,nCont1,nCont2,nPrim1,nPrim2,LUPRI,IPRINT
    Real(realk),intent(in)    :: CC1(nPrim1,nCont1)
    Real(realk),intent(in)    :: CC2(nPrim2,nCont2)
    integer,intent(in)        :: Piprim1(nPrim),Piprim2(nPrim)
    Real(realk),intent(inout) :: CC(nPrim,nCont1,nCont2)
    !
    Integer       :: iP,i1,i2,iC1,iC2,nP1,nP2
    real(realk)   :: TMP
    DO iP=1,nPrim
       i1 = Piprim1(iP)
       i2 = Piprim2(iP)
       DO iC2=1,nCont2
          TMP = CC2(i2,iC2) 
          DO iC1=1,nCont1
             CC(iP,iC1,iC2)=CC1(i1,iC1)*TMP
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE ConstructContraction_screen_ssss

  SUBROUTINE ConstructContraction_full_ssss(CC,nPrim,nCont1,nCont2,&
       & nPrim1,nPrim2,CC1,CC2,LUPRI,IPRINT)
    Implicit none
    Integer       :: nPrim,nCont1,nCont2,nPrim1,nPrim2,LUPRI,IPRINT
    Real(realk),intent(in)    :: CC1(nPrim1,nCont1)
    Real(realk),intent(in)    :: CC2(nPrim2,nCont2)
    Real(realk),intent(inout) :: CC(nPrim1,nPrim2,nCont1,nCont2)
    !
    Integer       :: iP,i1,i2,iC1,iC2,nP1,nP2
    real(realk)   :: TMP
    DO iC2=1,nCont2
       DO iC1=1,nCont1
          DO i2=1,nPrim2
             TMP = CC2(i2,iC2) 
             DO i1=1,nPrim1
                CC(i1,i2,iC1,iC2)=CC1(i1,iC1)*TMP
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE ConstructContraction_full_ssss


END MODULE IchorEriCoulombintegralMod
