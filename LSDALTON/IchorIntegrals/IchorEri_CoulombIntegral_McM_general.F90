!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriCoulombintegralGeneralMod
use precision
use IchorCommonModule
use ThermiteMem_module
use AGCContractEcoeffPSegSubsMod
use AGCDirectContractEcoeffSegSubsMod1
use AGCDirectContractEcoeffSegSubsMod2
use AGCDirectContractEcoeffSegSubsMod3
use AGCDirectContractEcoeffSegSubsMod4
use AGCDirectContractEcoeffSegSubsMod5
use AGCDirectContractEcoeffSegSubsMod6
use AGCDirectContractEcoeffSegSubsMod7
use AGCDirectContractEcoeffSegSubsMod8
use AGCSphPcontractSegSubsMod
use AGCSphQcontractSegSubsMod
private 
public :: IchorCoulombIntegral_McM_general
CONTAINS
  subroutine IchorCoulombIntegral_McM_general(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimPQ,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimPQ,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimPQ)
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimPQ)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses)
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)
    real(realk),pointer :: squaredDistance(:),Rpq(:)
    !
!    real(realk),target :: TMPWORK((2+2*nPasses)*nPrimP*nPrimQ)
    integer :: nPrimPassQ,AngmomP,AngmomQ,AngmomPQ,nTUV,nTUVQ
    logical :: PQorder
    integer :: ijk1,ijk2,ijkPcart,ijk1s,ijk2s,ijkPsph
    integer :: ijk3,ijk4,ijkQcart,ijk3s,ijk4s,ijkQsph,nPassesP,nTUVP
    logical :: Sph1,Sph2,sphericalGTO,Sph3,Sph4,SphericalTransP,SphericalTransQ
    integer :: AngmomID
    real(realk),pointer :: RJ000(:),WTUV(:),RE(:),EcoeffN(:),ERE(:),SERE(:),SERES(:)
    integer :: TUV,J,T,U,V,I
    logical :: RHS
    sphericalGTO = .TRUE.
    !build from old IchorEri_CoulombIntegral_general.f90 in
    !/home/tkjaer/DaltonDevelopment/ExplicitIntegrals/LSint
    IF(PQorder)THEN
       call ichorquit('PQorder McM general expect to get QP ordering',-1)
    ENDIF
    !TODO ipassP!
    nPrimPassQ = nPrimQ*nPasses 
    call mem_alloc(squaredDistance,nPasses*nPrimPQ)
    call mem_alloc(Rpq,3*nPasses*nPrimPQ)
    !builds squaredDistance(nPrimQ,nPrimP,nPassQ)
    !builds Rpq(3,nPrimQ,nPrimP,nPassQ)
    call build_QP_squaredDistance_and_Rpq(nPrimQ,nPasses,nPrimP,Qcent,Pcent,&
         & squaredDistance,Rpq)

    AngmomP = AngmomA+AngmomB
    AngmomQ = AngmomC+AngmomD
    AngmomPQ  = AngmomP + AngmomQ
    call mem_alloc(RJ000,(AngmomPQ+1)*nPasses*nPrimPQ)
    call buildRJ000_general(nPasses,nPrimPQ,nTABFJW1,nTABFJW2,reducedExponents,&
         & squaredDistance,TABFJW,RJ000,AngmomPQ,integralPrefactor)
    !builds RJ000(0:AngmomPQ,nPrimQ,nPrimP,nPasses)

    IF (INTPRINT .GE. 10) THEN
       WRITE(lupri,*)'Output from W000'
       DO I=1,nPrimQ*nPrimP*nPasses
          DO J=0,AngmomPQ
             WRITE(LUPRI,'(2X,A6,I4,A1,I4,A2,ES16.8)')'RJ000(',J,',',I,')=',RJ000(1+J+(I-1)*(AngmomPQ+1))
          ENDDO
       ENDDO
    END IF

    nTUV=(AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6
    call mem_alloc(WTUV,nTUV*nPrimPQ*nPasses)
    IF(AngmomPQ.LT.5)then
       IF (AngmomPQ.EQ. 0) THEN
          call IchorwtuvRecurrenceJMIN0JMAX0(WTUV,RJ000,nPrimPQ*nPasses)
       ELSEIF (AngmomPQ.EQ. 1) THEN
          call IchorwtuvRecurrenceJMIN0JMAX1(WTUV,RJ000,Rpq,nPrimPQ*nPasses)
       ELSEIF (AngmomPQ.EQ. 2) THEN
          call IchorwtuvRecurrenceJMIN0JMAX2(WTUV,RJ000,Rpq,nPrimPQ*nPasses)
       ELSEIF (AngmomPQ.EQ. 3) THEN
          call IchorwtuvRecurrenceJMIN0JMAX3(WTUV,RJ000,Rpq,nPrimPQ*nPasses)
       ELSEIF (AngmomPQ.EQ. 4) THEN
          call IchorwtuvRecurrenceJMIN0JMAX4(WTUV,RJ000,Rpq,nPrimPQ*nPasses)
       ENDIF
    ELSE
       !generate more using RECURRENCE.F90 in tools
       call ichorquit('AngmomPQ .GT. 4',-1)
    ENDIF
    !builds WTUV(nTUV,nPrimP,nPrimQ,nPassesQ)
    IF (IntPrint .GE. 25) THEN
       WRITE(LUPRI,*)'Output from WTUVrecurrence'
       WRITE (LUPRI,'(2X,A,I10)') 'JMAX  ', AngmomPQ
       WRITE (LUPRI,'(2X,A,I10)') 'NPrim ', nPrimPQ*nPasses
       WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', nTUV
       WRITE(LUPRI,*)'Hermite integrals S(t,u,v)'
       TUV=0
       DO J = 0, AngmomPQ
          DO T = J,0,-1
             DO U = J-T,0,-1
                V=J-T-U
                TUV=TUV+1
                WRITE (LUPRI,'(2X,A2,I3,A1,I3,A1,I3,A1,2X,5ES16.8/,(18X,5ES16.8))')&
                     & 'W(',T,',',U,',',V,')', &
                     &(WTUV(TUV + (I-1)*nTUV),I=1,nPrimPQ*nPasses)
                WRITE (LUPRI,*) ' '
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    ijk1 = (AngmomA + 1)*(AngmomA + 2)/2
    ijk2 = (AngmomB + 1)*(AngmomB + 2)/2
    ijk3 = (AngmomC + 1)*(AngmomC + 2)/2
    ijk4 = (AngmomD + 1)*(AngmomD + 2)/2
    ijkPcart = ijk1*ijk2
    ijkQcart = ijk3*ijk4
    ijk1s = 2*AngmomA + 1
    ijk2s = 2*AngmomB + 1
    ijk3s = 2*AngmomC + 1
    ijk4s = 2*AngmomD + 1
    ijkPsph = ijk1s*ijk2s
    ijkQsph = ijk3s*ijk4s
    
    nPassesP = 1

    nTUVP=(AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
    nTUVQ=(AngmomQ+1)*(AngmomQ+2)*(AngmomQ+3)/6

    Sph1 = sphericalGTO.AND.(AngmomA.GT. 1)
    Sph2 = sphericalGTO.AND.(AngmomB.GT. 1)
    Sph3 = sphericalGTO.AND.(AngmomC.GT. 1)
    Sph4 = sphericalGTO.AND.(AngmomD.GT. 1)
    SphericalTransP = Sph1.OR.Sph2
    SphericalTransQ = Sph3.OR.Sph4

    call mem_alloc(EcoeffN,nTUVQ*ijkQcart*nPrimQ*nPasses)

    !WARNING ONLY ZERO FOR DEBUGGING PURPOSES
    CALL LS_DZERO(EcoeffN,nTUVQ*ijkQcart*nPrimQ*nPasses)
    
    RHS=.TRUE.
    call Ichorbuild_Ecoeff(nPrimQ,nPrimC,nPrimD,AngmomQ,AngmomC,AngmomD,nTUVQ,&
         & ijkQcart,Cexp,Dexp,EcoeffN,Qdistance12,Qpreexpfac,nPasses,RHS,&
         & intprint,lupri)
    !builds Ecoeff(nTUVQ,ijkQcart,nPrimQ,nPasses)
    IF (IntPrint .GE. 25) THEN
       call printEcoeff(EcoeffN,nTUVQ,ijkQcart,nPrimQ,nPasses,lupri)
    ENDIF
    IF(Psegmented.AND.Qsegmented)THEN
     !RE(nTUVP,ijkQcart,nPrimP,nPasses)
     call mem_alloc(RE,nTUVP*ijkQcart*nPrimP*nPasses)

     !WARNING ONLY ZERO FOR DEBUGGING PURPOSES
     CALL LS_DZERO(RE,nTUVP*ijkQcart*nPrimP*nPasses)

     !Max AngmomA:AngmomB:AngmomC:AngmomD  0:9     
     AngmomID = AngmomC + AngmomD*10 + AngmomP*100
     SELECT CASE(AngmomID)
     CASE(0) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ0_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(100) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ0_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(200) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ0_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(300) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ0_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(400) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ0_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(1) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ1_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(101) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ1_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(201) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ1_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(301) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ1_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(401) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ1_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(10) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ1_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(110) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ1_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(210) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ1_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(310) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ1_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(410) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ1_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(2) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ2_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(102) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ2_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(202) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ2_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(302) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ2_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(402) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ2_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(11) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ2_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(111) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ2_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(211) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ2_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(311) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ2_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(411) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ2_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(20) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ2_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(120) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ2_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(220) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ2_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(320) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ2_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(420) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ2_maxAngC0(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(12) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ3_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(112) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ3_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(212) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ3_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(312) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ3_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(412) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ3_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(21) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ3_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(121) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ3_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(221) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ3_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(321) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ3_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(421) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ3_maxAngC1(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(22) 
        call DirectContractEcoeffN_seg_maxAngP0_maxAngQ4_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(122) 
        call DirectContractEcoeffN_seg_maxAngP1_maxAngQ4_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(222) 
        call DirectContractEcoeffN_seg_maxAngP2_maxAngQ4_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(322) 
        call DirectContractEcoeffN_seg_maxAngP3_maxAngQ4_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE(422) 
        call DirectContractEcoeffN_seg_maxAngP4_maxAngQ4_maxAngC2(&
             & nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,RE)
     CASE DEFAULT        
        call ichorquit('DirectContractEcoeffN',-1)
     END SELECT
     !
     !RE(nTUVP,ijkQcart,nPrimP,nPasses)
!IFDEF
    IF (IntPrint .GE. 25) THEN
       call PrintIchorTensorRE(RE,nTUVP,ijkQcart,nPrimP,nPasses,lupri)
    ENDIF
    call mem_alloc(EcoeffN,nTUVP*ijkPcart*nPrimP)
    !WARNING ONLY ZERO FOR DEBUGGING PURPOSES
    CALL LS_DZERO(EcoeffN,nTUVP*ijkPcart*nPrimP)
    RHS=.FALSE.
    call Ichorbuild_Ecoeff(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
         & ijkPcart,Aexp,Bexp,EcoeffN,Pdistance12,Ppreexpfac,nPassesP,RHS,&
         & intprint,lupri)
    !builds Ecoeff(nTUVP,ijkPcart,nPrimP,1)

!IFDEF
    IF (IntPrint .GE. 25) THEN
       call printEcoeff(EcoeffN,nTUVP,ijkPcart,nPrimP,1,lupri)
    ENDIF
    call mem_alloc(ERE,ijkPcart*ijkQcart*nPasses)
    !WARNING ONLY ZERO FOR DEBUGGING PURPOSES
    CALL LS_DZERO(ERE,ijkPcart*ijkQcart*nPasses)
     

     AngmomID = AngmomA + AngmomP*10
     SELECT CASE(AngmomID)
     CASE(0) 
        call ContractPEcoeffN_seg_maxAngP0_maxAngA0(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE(10) 
        call ContractPEcoeffN_seg_maxAngP1_maxAngA0(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE(11) 
        call ContractPEcoeffN_seg_maxAngP1_maxAngA1(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE(20) 
        call ContractPEcoeffN_seg_maxAngP2_maxAngA0(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE(21) 
        call ContractPEcoeffN_seg_maxAngP2_maxAngA1(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE(22) 
        call ContractPEcoeffN_seg_maxAngP2_maxAngA2(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE(31) 
        call ContractPEcoeffN_seg_maxAngP3_maxAngA1(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE(32) 
        call ContractPEcoeffN_seg_maxAngP3_maxAngA2(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE(42) 
        call ContractPEcoeffN_seg_maxAngP4_maxAngA2(&
             & nPrimP,nPasses,ijkQcart,Ecoeffn,RE,ERE)
     CASE DEFAULT        
        call ichorquit('ContractPEcoeffN',-1)
     END SELECT
     !builds ERE(ijkPcart,ijkQcart,nPasses)
!IFDEF
     IF (IntPrint .GE. 25) THEN
        call PrintIchorTensorERE(ERE,ijkPcart,ijkQcart,nPasses,lupri)
     ENDIF
     call mem_dealloc(RE)
     
     
     IF(SphericalTransP)THEN
        !Spherical Transformation of ERE(ijkPcart,ijkQcart,nPasses)
        !to ERE(ijkPsph,ijkQcart,nPasses)
        call mem_alloc(SERE,ijkPsph*ijkQcart*nPasses)
        !WARNING ONLY ZERO FOR DEBUGGING PURPOSES
        CALL LS_DZERO(SERE,ijkPsph*ijkQcart*nPasses)

        !we should contract the one with the highest angular momentum first - on the 
        !other hand the P is consequtive in memory so maybe it is best to do that one first. 

        AngmomID = AngmomA + AngmomP*10
        SELECT CASE(AngmomID)
!        CASE(0) 
!           call SphericalContractP_maxAngP0_maxAngA0(ijkQcart,nPasses,ERE,SERE)
!        CASE(10) 
!           call SphericalContractP_maxAngP1_maxAngA0(ijkQcart,nPasses,ERE,SERE)
!        CASE(11) 
!           call SphericalContractP_maxAngP1_maxAngA1(ijkQcart,nPasses,ERE,SERE)
        CASE(20) 
           call SphericalContractP_maxAngP2_maxAngA0(ijkQcart,nPasses,ERE,SERE)
!        CASE(21) 
!           call SphericalContractP_maxAngP2_maxAngA1(ijkQcart,nPasses,ERE,SERE)
        CASE(22) 
           call SphericalContractP_maxAngP2_maxAngA2(ijkQcart,nPasses,ERE,SERE)
! TODO: SphericalContractP_maxAngP2_maxAngA0 and SphericalContractP_maxAngP2_maxAngA2
!       is identical 

        CASE(31) 
           call SphericalContractP_maxAngP3_maxAngA1(ijkQcart,nPasses,ERE,SERE)
        CASE(32) 
           call SphericalContractP_maxAngP3_maxAngA2(ijkQcart,nPasses,ERE,SERE)
        CASE(42) 
           call SphericalContractP_maxAngP4_maxAngA2(ijkQcart,nPasses,ERE,SERE)
        CASE DEFAULT        
           call ichorquit('ContractSphP',-1)
        END SELECT
        call mem_dealloc(ERE)        
     ELSE
        SERE => ERE
     ENDIF

     IF(SphericalTransQ)THEN
        !Spherical Transformation of SERE(ijkPsph,ijkQcart,nPasses)
        !to SERES(ijkPsph,ijkQsph,nPasses)
        call mem_alloc(SERES,ijkPsph*ijkQsph*nPasses)
        !WARNING ONLY ZERO FOR DEBUGGING PURPOSES
        CALL LS_DZERO(SERES,ijkPsph*ijkQsph*nPasses)

        AngmomID = AngmomC + AngmomQ*10
        SELECT CASE(AngmomID)
!        CASE(0) 
!           call SphericalContractQ_maxAngQ0_maxAngC0(ijkPsph,nPasses,SERE,SERES)
!        CASE(10) 
!           call SphericalContractQ_maxAngQ1_maxAngC0(ijkPsph,nPasses,SERE,SERES)
!        CASE(11) 
!           call SphericalContractQ_maxAngQ1_maxAngC1(ijkPsph,nPasses,SERE,SERES)
        CASE(20) 
           call SphericalContractQ_maxAngQ2_maxAngC0(ijkPsph,nPasses,SERE,SERES)
!        CASE(21) 
!           call SphericalContractQ_maxAngQ2_maxAngC1(ijkPsph,nPasses,SERE,SERES)
        CASE(22) 
           call SphericalContractQ_maxAngQ2_maxAngC2(ijkPsph,nPasses,SERE,SERES)
        CASE(31) 
           call SphericalContractQ_maxAngQ3_maxAngC1(ijkPsph,nPasses,SERE,SERES)
        CASE(32) 
           call SphericalContractQ_maxAngQ3_maxAngC2(ijkPsph,nPasses,SERE,SERES)
        CASE(42) 
           call SphericalContractQ_maxAngQ4_maxAngC2(ijkPsph,nPasses,SERE,SERES)
        CASE DEFAULT        
           call ichorquit('ContractSphQ',-1)
        END SELECT
        call mem_dealloc(ERE)        
     ELSE
        SERES => SERE
     ENDIF

     !Now CDAB
     print*,'ijkPcart',ijkPcart
     print*,'ijkQcart',ijkQcart
     print*,'ijkPsph',ijkPsph
     print*,'ijkQsph',ijkQsph
     print*,'dim',ijkQsph*ijkPsph*nPasses
     call PQreorderSeg(ijkPsph,ijkQsph,nPasses,CDAB,SERES)

    ELSEIF(Psegmented)THEN

    ELSEIF(Qsegmented)THEN

    ELSE
!!$       !WTUV(nTUV,nPrimP,nPrimQ,nPassesQ)*EcoeffN(nTUVQ,ijkQcart,nPrimQ,nPasses)
!!$       !OUTPUT(ijkQcart,nTUVP,nPrimP,nPrimQ,nPasses)
!!$       call mem_alloc(RE,ijkQcart*nPrimP*nPrimQ*nPasses)
!!$       IF (AngmomPQ.EQ. 0) THEN
!!$          call DirectContractEcoeffN_maxAngP0_maxAngQ0_maxAngC0(nPrimP,nPrimQ,&
!!$               & nPasses,EcoeffN,WTUV,RE)
!!$       ENDIF
!!$       call mem_dealloc(EcoeffN)
!!$
!!$       IF(.NOT.Qsegmented)THEN
!!$          !ContInt(ijkQcart,nTUVP,nPrimP,nContQ,nPasses) = CC*PrimInt(ijkQcart,nTUVP,nPrimP,nPrimQ,nPasses)
!!$          call ichorquit('ContQ',-1)
!!$       ENDIF
!!$
!!$       !SphericalTransform
!!$       !ContInt(ijkQsph,nTUVP,nPrimP,nContQ,nPasses) = CC*ContInt(ijkQcart,nTUVP,nPrimP,nContQ,nPasses)
!!$       
!!$       nPassesP = 1
!!$       nTUVP=(AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
!!$       ijk1 = (AngmomA + 1)*(AngmomA + 2)/2
!!$       ijk2 = (AngmomB + 1)*(AngmomB + 2)/2
!!$       ijkPcart = ijk1*ijk2
!!$       call mem_alloc(EcoeffN,nTUVP*ijkPcart*nPrimP)
!!$       call Ichorbuild_Ecoeff(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
!!$            & ijkPcart,Aexp,Bexp,EcoeffN,Pdistance12,Ppreexpfac,nPassesP,lupri)
!!$       !builds Ecoeff(nTUVP,ijkPcart,nPrimP,1)
!!$       
!!$       call mem_alloc(ERECS,ijkQsph*nPrimP*nContQ*nPasses)
!!$       !OUT(ijkQsph,ijkPcart,nPrimP,nContQ,nPasses) =+ IN(ijkQsph,nTUVP,nPrimP,nContQ,nPasses)*Etuv(nTUVP,ijkPcart,nPrimP)
!!$       IF(Psegmented)THEN
!!$          call contractEcoeffSeg(RE,ERECS,EcoeffN,&
!!$               & ijkQsph,ijkPcart,nTUVP,nPrimP,nContQ,nPasses)
!!$       ELSE
!!$          call contractEcoeffGen(RE,ERECS,EcoeffN,&
!!$               & ijkQsph,ijkPcart,nTUVP,nPrimP,nContQ,nPasses)
!!$       ENDIF
!!$       
!!$       IF(Psegmented)THEN
!!$          
!!$       ELSE
!!$          !ContInt(ijkQcart,nTUVP,nPrimP,nContQ,nPasses) = CC*PrimInt(ijkQcart,nTUVP,nPrimP,nPrimQ,nPasses)
!!$          call ichorquit('ContQ',-1)
!!$       ENDIF
    ENDIF
  end subroutine IchorCoulombIntegral_McM_general

  subroutine PQreorderSeg(ijkPsph,ijkQsph,nPasses,CDAB,IN)
    implicit none
    integer,intent(in) :: ijkPsph,ijkQsph,nPasses
    real(realk),intent(in) :: IN(ijkPsph,ijkQsph,nPasses)
    real(realk),intent(inout) :: CDAB(ijkQsph,ijkPsph,nPasses)
!
    integer :: ijkQ,ijkP,iPass
    do iPass = 1,nPasses
     do ijkP=1,ijkPsph
      do ijkQ=1,ijkPsph
       CDAB(ijkQ,ijkP,iPass) = IN(ijkP,ijkQ,iPass)
      enddo
     enddo
    enddo
  end subroutine PQreorderSeg

  subroutine build_QP_squaredDistance_and_Rpq(nPrimQ,nPassQ,nPrimP,Qcent,Pcent,squaredDistance,Rpq)
    implicit none
    integer,intent(in) :: nPrimQ,nPassQ,nPrimP
    real(realk),intent(in) :: Qcent(3,nPrimQ,nPassQ),Pcent(3,nPrimP)
    real(realk),intent(inout) :: squaredDistance(nPrimQ,nPrimP,nPassQ)
    real(realk),intent(inout) :: Rpq(3,nPrimQ,nPrimP,nPassQ)
    !
    integer :: iPrimQ,iPassQ,iPrimP
    real(realk) :: px,py,pz,pqx,pqy,pqz
    DO iPassQ=1, nPassQ
       DO iPrimP=1, nPrimP
          px = Pcent(1,iPrimP)
          py = Pcent(2,iPrimP)
          pz = Pcent(3,iPrimP)
          DO iPrimQ=1, nPrimQ
             pqx = px - Qcent(1,iPrimQ,nPassQ)
             pqy = py - Qcent(2,iPrimQ,nPassQ)
             pqz = pz - Qcent(3,iPrimQ,nPassQ)
             rPQ(1,iPrimQ,iPrimP,iPassQ) = pqx
             rPQ(2,iPrimQ,iPrimP,iPassQ) = pqy
             rPQ(3,iPrimQ,iPrimP,iPassQ) = pqz
             squaredDistance(iPrimQ,iPrimP,iPassQ) = pqx*pqx+pqy*pqy+pqz*pqz
          ENDDO
       ENDDO
    ENDDO
  end subroutine build_QP_squaredDistance_and_Rpq

  SUBROUTINE buildRJ000_general(nPasses,nPrimPQ,nTABFJW1,nTABFJW2,reducedExponents,&
       & squaredDistance,TABFJW,RJ000,JMAX,integralPrefactor)
    IMPLICIT NONE
    INTEGER,intent(in)         :: nPrimPQ,Jmax,nTABFJW1,nTABFJW2,nPasses
    REAL(REALK),intent(in)     :: reducedExponents(nPrimPQ)
    REAL(REALK),intent(in)     :: squaredDistance(nPrimPQ,nPasses)
    REAL(REALK),intent(in)     :: integralPrefactor(nprimPQ)
    REAL(REALK),intent(in)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrimPQ,nPasses)
    !
    REAL(REALK)     :: D2JP36,WVAL
    REAL(REALK),PARAMETER :: HALF =0.5E0_realk,D1=1E0_realk,D2 = 2E0_realk, D4 = 4E0_realk, D100=100E0_realk
    Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk
    REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
    REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk
    Integer :: IPNT,J
    Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
    REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
    REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
    REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
    REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
    Real(realk), parameter :: PI=3.14159265358979323846E0_realk
    REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
    REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
    REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
    Real(realk) :: W2,W3,R
    REAL(REALK), PARAMETER :: SMALL = 1E-15_realk
    integer :: iPassQ,iPQ

    D2JP36 = 2*JMAX + 36
    DO iPassq=1, nPasses
       DO ipq=1, nPrimPQ
          !(nPrimP,nPrimQ)*(nPrimP,nPrimQ,nPasses)
          WVAL = reducedExponents(ipq)*squaredDistance(ipq,iPassQ)
          !  0 < WVAL < 0.000001
          IF (ABS(WVAL) .LT. SMALL) THEN         
             RJ000(0,ipq,ipassq) = D1
             DO J=1,JMAX
                RJ000(J,ipq,ipassq)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
             ENDDO
             !  0 < WVAL < 12 
          ELSE IF (WVAL .LT. D12) THEN
             IPNT = NINT(D100*WVAL)
             WDIFF = WVAL - TENTH*IPNT
             W2    = WDIFF*WDIFF
             W3    = W2*WDIFF
             W2    = W2*COEF2
             W3    = W3*COEF3
             DO J=0,JMAX
                R = TABFJW(J,IPNT)
                R = R -TABFJW(J+1,IPNT)*WDIFF
                R = R + TABFJW(J+2,IPNT)*W2
                R = R + TABFJW(J+3,IPNT)*W3
                RJ000(J,ipq,ipassq) = R
             ENDDO
             !  12 < WVAL <= (2J+36) 
          ELSE IF (WVAL.LE.D2JP36) THEN
             REXPW = HALF*EXP(-WVAL)
             RWVAL = D1/WVAL
             GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
             RJ000(0,ipq,ipassq) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
             DO J=1,JMAX
                RJ000(J,ipq,ipassq) = RWVAL*((J - HALF)*RJ000(J-1,ipq,ipassq)-REXPW)
             ENDDO
             !  (2J+36) < WVAL 
          ELSE
             RWVAL = PID4/WVAL
             RJ000(0,ipq,ipassq) = SQRT(RWVAL)
             RWVAL = RWVAL*PID4I
             DO J = 1, JMAX
                RJ000(J,ipq,ipassq) = RWVAL*(J - HALF)*RJ000(J-1,ipq,ipassq)
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    ! Scaling
    DO iPassq=1, nPasses
       DO ipq=1, nPrimPQ
          PREF = integralPrefactor(ipq)
          RJ000(0,ipq,iPassq) = PREF*RJ000(0,ipq,iPassq)
!          IF (jmax.GT. 0) THEN
             D2MALPHA = -D2*reducedExponents(ipq)
             DO j=1,jmax
                PREF = PREF*D2MALPHA
                RJ000(J,ipq,iPassq) = PREF*RJ000(J,ipq,iPassq)
             ENDDO
!          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE buildRJ000_general

  SUBROUTINE IchorwtuvRecurrenceJMIN0JMAX0(WTUV,WJ000,nPrim)
    implicit none
    Integer,intent(in)           :: nPrim
    Real(realk),intent(in)       :: WJ000(0:0,nPrim)
    Real(realk),intent(inout)    :: WTUV(1,nPrim)
    !
    integer :: n
    DO n=1,nPrim
       WTUV(1,n) = WJ000(0,n)
    ENDDO
  END SUBROUTINE ICHORWTUVRECURRENCEJMIN0JMAX0

SUBROUTINE ichorwtuvrecurrenceJMIN0JMAX1(WTUV,WJ000,Rpq,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:1,nPrim),Rpq(3,nPrim)
  Real(realk),intent(inout)    :: WTUV(4,nPrim)
  !
  integer :: n
  DO n=1,nPrim
     WTUV(1,n) = WJ000(0,n)
     WTUV(2,n) = Rpq(1,n)*WJ000(1,n)
     WTUV(3,n) = Rpq(2,n)*WJ000(1,n)
     WTUV(4,n) = Rpq(3,n)*WJ000(1,n)
  ENDDO
END SUBROUTINE ICHORWTUVRECURRENCEJMIN0JMAX1

SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX2(WTUV,WJ000,Rpq,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:2,nPrim),Rpq(3,nPrim)
  Real(realk),intent(inout)    :: WTUV(10,nPrim)
  !
  integer :: n
  !000
  DO n=1,nPrim
     WTUV(1,n) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(2,n) = Rpq(1,n)*WJ000(1,n)
     WTUV(3,n) = Rpq(2,n)*WJ000(1,n)
     WTUV(4,n) = Rpq(3,n)*WJ000(1,n)
     WTUV(5,n) = WJ000(1,n) + Rpq(1,n)*Rpq(1,n)*WJ000(2,n)
     WTUV(8,n) = WJ000(1,n) + Rpq(2,n)*Rpq(2,n)*WJ000(2,n)
     WTUV(10,n)= WJ000(1,n) + Rpq(3,n)*Rpq(3,n)*WJ000(2,n)
     WTUV(6,n) = Rpq(1,n)*Rpq(2,n)*WJ000(2,n)
     WTUV(7,n) = Rpq(1,n)*Rpq(3,n)*WJ000(2,n)
     WTUV(9,n) = Rpq(2,n)*Rpq(3,n)*WJ000(2,n)
  ENDDO
end SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX2

SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX3(WTUV,WJ000,Rpq,nPrim)
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:3,nPrim),Rpq(3,nPrim)
  Real(realk),intent(inout)    :: WTUV(20,nPrim)
  !
  Real(realk) :: WTUVTEMP(2:4,2)
  Real(realk) :: WTUVTEMP25,WTUVTEMP28
  Real(realk) :: WTUVTEMP210,WTUVTEMP29
  integer :: n
  DO n=1,nPrim
     WTUV(1,n) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(2,n)     = Rpq(1,n)*WJ000(1,n)
     WTUVTEMP(2,1) = Rpq(1,n)*WJ000(2,n)
     WTUVTEMP(2,2) = Rpq(1,n)*WJ000(3,n)
     WTUV(3,n)     = Rpq(2,n)*WJ000(1,n)
     WTUVTEMP(3,1) = Rpq(2,n)*WJ000(2,n)
     WTUVTEMP(3,2) = Rpq(2,n)*WJ000(3,n)
     WTUV(4,n)     = Rpq(3,n)*WJ000(1,n)
     WTUVTEMP(4,1) = Rpq(3,n)*WJ000(2,n)
     WTUVTEMP(4,2) = Rpq(3,n)*WJ000(3,n)
     WTUV(6,n) = Rpq(1,n)*WTUVTEMP(3,1)
     WTUV(7,n) = Rpq(1,n)*WTUVTEMP(4,1)
     WTUV(9,n)  = Rpq(2,n)*WTUVTEMP(4,1)
     WTUVTEMP29 = Rpq(2,n)*WTUVTEMP(4,2)
     WTUV(5,n)  = WJ000(1,n) + Rpq(1,n)*WTUVTEMP(2,1)
     WTUV(8,n)  = WJ000(1,n) + Rpq(2,n)*WTUVTEMP(3,1)
     WTUV(10,n) = WJ000(1,n) + Rpq(3,n)*WTUVTEMP(4,1)
     WTUVTEMP25 = WJ000(2,n) + Rpq(1,n)*WTUVTEMP(2,2)
     WTUVTEMP28 = WJ000(2,n) + Rpq(2,n)*WTUVTEMP(3,2)
     WTUVTEMP210 = WJ000(2,n) + Rpq(3,n)*WTUVTEMP(4,2)
     WTUV(14,n)= Rpq(1,n)*WTUVTEMP28
     WTUV(15,n)= Rpq(1,n)*WTUVTEMP29
     WTUV(16,n)= Rpq(1,n)*WTUVTEMP210
     WTUV(12,n)= Rpq(2,n)*WTUVTEMP25
     WTUV(19,n)= Rpq(2,n)*WTUVTEMP210
     WTUV(13,n)= Rpq(3,n)*WTUVTEMP25
     WTUV(18,n)= Rpq(3,n)*WTUVTEMP28
     WTUV(11,n)= 2*WTUVTEMP(2,1) + Rpq(1,n)*WTUVTEMP25
     WTUV(17,n)= 2*WTUVTEMP(3,1) + Rpq(2,n)*WTUVTEMP28
     WTUV(20,n)= 2*WTUVTEMP(4,1) + Rpq(3,n)*WTUVTEMP210
  ENDDO
END SUBROUTINE ichorWTUVRECURRENCEJMIN0JMAX3

SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX4(WTUV,WJ000,Rpq,nPrim)
  use memory_handling
  implicit none
  Integer,intent(in)           :: nPrim
  Real(realk),intent(in)       :: WJ000(0:4,nPrim),Rpq(3,nPrim)
  Real(realk),intent(inout)    :: WTUV(35,nPrim)
  !
  Integer     :: n,M
  Real(realk) :: WTUVTEMP(2:4,3)
  Real(realk) :: WTUVTEMP2(5:10,2)
  Real(realk) :: WTUVTEMP3(11:20)
  !000
  DO n=1,nPrim
     WTUV(1,n) = WJ000(0,n)
  ENDDO
  !it can reuse registers with this collection of this in 1 loop
  DO n=1,nPrim
     WTUV(2,n) = Rpq(1,n)*WJ000(1,n)
     WTUV(3,n) = Rpq(2,n)*WJ000(1,n)
     WTUV(4,n) = Rpq(3,n)*WJ000(1,n)
     WTUVTEMP(2,1) = Rpq(1,n)*WJ000(2,n)
     WTUVTEMP(3,1) = Rpq(2,n)*WJ000(2,n)
     WTUVTEMP(4,1) = Rpq(3,n)*WJ000(2,n)     
     WTUVTEMP(2,2) = Rpq(1,n)*WJ000(3,n)
     WTUVTEMP(3,2) = Rpq(2,n)*WJ000(3,n)
     WTUVTEMP(4,2) = Rpq(3,n)*WJ000(3,n)     
     WTUVTEMP(2,3) = Rpq(1,n)*WJ000(4,n)
     WTUVTEMP(3,3) = Rpq(2,n)*WJ000(4,n)
     WTUVTEMP(4,3) = Rpq(3,n)*WJ000(4,n)
     WTUV(5,n) = WJ000(1,n) + Rpq(1,n)*WTUVTEMP(2,1)
     WTUV(6,n) =              Rpq(1,n)*WTUVTEMP(3,1)
     WTUV(7,n) =              Rpq(1,n)*WTUVTEMP(4,1)
     WTUV(8,n) = WJ000(1,n) + Rpq(2,n)*WTUVTEMP(3,1)
     WTUV(9,n) =              Rpq(2,n)*WTUVTEMP(4,1)
     WTUV(10,n)= WJ000(1,n) + Rpq(3,n)*WTUVTEMP(4,1)
     WTUVTEMP2(5,1) = WJ000(2,n) + Rpq(1,n)*WTUVTEMP(2,2)
     WTUVTEMP2(6,1) =              Rpq(1,n)*WTUVTEMP(3,2)
     WTUVTEMP2(7,1) =              Rpq(1,n)*WTUVTEMP(4,2)
     WTUVTEMP2(8,1) = WJ000(2,n) + Rpq(2,n)*WTUVTEMP(3,2)
     WTUVTEMP2(9,1) =              Rpq(2,n)*WTUVTEMP(4,2)
     WTUVTEMP2(10,1)= WJ000(2,n) + Rpq(3,n)*WTUVTEMP(4,2)     
     WTUVTEMP2(5,2) = WJ000(3,n) + Rpq(1,n)*WTUVTEMP(2,3)
!     WTUVTEMP2(6,2) =              Rpq(1,n)*WTUVTEMP(3,3)
!     WTUVTEMP2(7,2) =              Rpq(1,n)*WTUVTEMP(4,3)
     WTUVTEMP2(8,2) = WJ000(3,n) + Rpq(2,n)*WTUVTEMP(3,3)
     WTUVTEMP2(9,2) =              Rpq(2,n)*WTUVTEMP(4,3)
     WTUVTEMP2(10,2)= WJ000(3,n) + Rpq(3,n)*WTUVTEMP(4,3)
     WTUV(11,n)= 2*WTUVTEMP(2,1) + Rpq(1,n)*WTUVTEMP2(5,1)
     WTUV(12,n)=  Rpq(2,n)*WTUVTEMP2(5,1)
     WTUV(13,n)= Rpq(3,n)*WTUVTEMP2(5,1)
     WTUV(14,n)= Rpq(1,n)*WTUVTEMP2(8,1)
     WTUV(15,n)= Rpq(1,n)*WTUVTEMP2(9,1)
     WTUV(16,n)= Rpq(1,n)*WTUVTEMP2(10,1)
     WTUV(17,n)= 2*WTUVTEMP(3,1) + Rpq(2,n)*WTUVTEMP2(8,1)
     WTUV(18,n)= Rpq(3,n)*WTUVTEMP2(8,1)
     WTUV(19,n)= Rpq(2,n)*WTUVTEMP2(10,1)
     WTUV(20,n)= 2*WTUVTEMP(4,1) + Rpq(3,n)*WTUVTEMP2(10,1)
     WTUVTEMP3(11)= 2*WTUVTEMP(2,2) + Rpq(1,n)*WTUVTEMP2(5,2)
    ! WTUVTEMP3(12)=  Rpq(2,n)*WTUVTEMP2(5,2)
     WTUVTEMP3(13)= Rpq(3,n)*WTUVTEMP2(5,2)
     WTUVTEMP3(14)= Rpq(1,n)*WTUVTEMP2(8,2)
    ! WTUVTEMP3(15)= Rpq(1,n)*WTUVTEMP2(9,2)
     WTUVTEMP3(16)= Rpq(1,n)*WTUVTEMP2(10,2)
     WTUVTEMP3(17)= 2*WTUVTEMP(3,2) + Rpq(2,n)*WTUVTEMP2(8,2)
     WTUVTEMP3(18)= Rpq(3,n)*WTUVTEMP2(8,2)
     WTUVTEMP3(19)= Rpq(2,n)*WTUVTEMP2(10,2)
     WTUVTEMP3(20)= 2*WTUVTEMP(4,2) + Rpq(3,n)*WTUVTEMP2(10,2)
     WTUV(21,n)= 3*WTUVTEMP2(5,1) + Rpq(1,n)*WTUVTEMP3(11)
     WTUV(22,n)= Rpq(2,n)*WTUVTEMP3(11)
     WTUV(23,n)= Rpq(3,n)*WTUVTEMP3(11)
     WTUV(24,n)= WTUVTEMP2(8,1) + Rpq(1,n)*WTUVTEMP3(14)
     WTUV(25,n)= Rpq(2,n)*WTUVTEMP3(13)
     WTUV(26,n)= WTUVTEMP2(10,1) + Rpq(1,n)*WTUVTEMP3(16)
     WTUV(27,n)= Rpq(1,n)*WTUVTEMP3(17)
     WTUV(28,n)= Rpq(1,n)*WTUVTEMP3(18)
     WTUV(29,n)= Rpq(1,n)*WTUVTEMP3(19)
     WTUV(30,n)= Rpq(1,n)*WTUVTEMP3(20)
     WTUV(31,n)= 3*WTUVTEMP2(8,1) + Rpq(2,n)*WTUVTEMP3(17) 
     WTUV(32,n)= Rpq(3,n)*WTUVTEMP3(17)
     WTUV(33,n)= WTUVTEMP2(10,1) + Rpq(2,n)*WTUVTEMP3(19)
     WTUV(34,n)= Rpq(2,n)*WTUVTEMP3(20)
     WTUV(35,n)= 3*WTUVTEMP2(10,1) + Rpq(3,n)*WTUVTEMP3(20)
  ENDDO
END SUBROUTINE ichorwtuvRecurrenceJMIN0JMAX4

Subroutine Ichorbuild_Ecoeff(nPrimP,nPrimA,nPrimB,maxAngP,maxAng1,maxAng2,nTUV,&
     & ijk,Aexp,Bexp,EcoeffN,Pdistance12,pref,nPass,RHS,intprint,lupri)
  implicit none
  integer,intent(in)        :: nPrimP,nPrimA,nPrimB,lupri,nTUV,ijk
  integer,intent(in)        :: maxAngP,maxAng1,maxAng2,nPass,intprint
  real(realk),intent(in)    :: Aexp(nPrimA),Bexp(nPrimB),Pdistance12(3,nPass)
  real(realk),intent(in)    :: pref(nPrimP)
  real(realk),intent(inout) :: Ecoeffn(nTUV,ijk,nPrimP,nPass)
  logical,intent(in)        :: RHS
  !
  real(realk) :: PINV(nPrimP),PA(nPrimP,nPass),PB(nPrimP,nPass)
  real(realk) :: APINV(nPrimP),BPINV(nPrimP)
  real(realk) :: HPINV(nPrimP),TWOA(nPrimP)
  real(realk) :: TWOB(nPrimP),dist
  real(realk) :: ETIJ(0:maxAngP,0:maxAng1,0:maxAng2,3,nPrimP,nPass)
  real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk
  !
  integer :: K,iPass,iPrimP
  call build_auxiliary(nPrimA,nPrimB,Aexp,Bexp,PINV,APINV,BPINV,HPINV,TWOA,TWOB)
  !X DIST
  do iPass=1,nPass
     dist=-Pdistance12(1,ipass)
     DO iPrimP=1,nPrimP
        PA(iPrimP,iPass) = BPINV(iPrimP)*dist
     ENDDO
     dist=Pdistance12(1,ipass)
     DO iPrimP=1,nPrimP
        PB(iPrimP,iPass) =  APINV(iPrimP)*dist
     ENDDO
  enddo
  CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPass,maxAng1,maxAng2,&
       & PA,PB,HPINV,TWOA,TWOB,1,LUPRI)  
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPass,maxAng1,maxAng2,1,LUPRI)
  ENDIF
  !Y DIST
  do iPass=1,nPass
     dist=-Pdistance12(2,ipass)
     DO iPrimP=1,nPrimP
        PA(iPrimP,iPass) = BPINV(iPrimP)*dist
     ENDDO
     dist=Pdistance12(2,ipass)
     DO iPrimP=1,nPrimP
        PB(iPrimP,iPass) =  APINV(iPrimP)*dist
     ENDDO
  enddo
  CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPass,maxAng1,maxAng2,&
       & PA,PB,HPINV,TWOA,TWOB,2,LUPRI)
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPass,maxAng1,maxAng2,2,LUPRI)
  ENDIF
  !Z DIST
  do iPass=1,nPass
     dist=-Pdistance12(3,ipass)
     DO iPrimP=1,nPrimP
        PA(iPrimP,iPass) = BPINV(iPrimP)*dist
     ENDDO
     dist=Pdistance12(3,ipass)
     DO iPrimP=1,nPrimP
        PB(iPrimP,iPass) =  APINV(iPrimP)*dist
     ENDDO
  enddo
  CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPass,maxAng1,maxAng2,&
       & PA,PB,HPINV,TWOA,TWOB,3,LUPRI)
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPass,maxAng1,maxAng2,3,LUPRI)
  ENDIF
  IF(maxAngP.EQ.0)THEN
     IF(maxAng1.EQ.0)THEN
        call ICHOR_Ecoeffn_maxAngP0_maxAngA0(nPrimP*nPass,Ecoeffn,ETIJ,pref)
     ENDIF
  ELSE
     call ICHOR_Ecoeffn_general(nPrimP*nPass,nTUV,ijk,maxAngP,maxAng1,maxAng2,Ecoeffn,ETIJ,pref,RHS)
  ENDIF
end Subroutine Ichorbuild_Ecoeff
  
SUBROUTINE PRINT_ETIJ(ETIJ,nPrimP,nPass,MAXI,MAXJ,X,LUPRI)
implicit none
INTEGER          ::  MAXI,MAXJ,X,nPrimP,nPass,LUPRI
CHARACTER(len=4) ::  WORD
REAL(REALK)      ::  ETIJ(0:MAXI+MAXJ,0:MAXI,0:MAXJ,3,nPrimP,nPass)
!
INTEGER          ::  I,J,T,K,IP

WRITE(LUPRI,*)'Output from Icor ETIJ'
WRITE(LUPRI,*)'----------------------'
WRITE (LUPRI,'(2X,A,2I5,I7,I7)') 'MAXI,MAXJ,nPrim,nPass:', MAXI, MAXJ,nPrimP,nPass
IF (X .EQ. 1) WORD = 'EX00'
IF (X .EQ. 2) WORD = 'EY00'
IF (X .EQ. 3) WORD = 'EZ00'
DO IP = 1, nPass
   DO I = 0, MAXI 
      DO J = 0, MAXJ 
         DO T = 0, I + J
            WRITE (LUPRI,'(/,2X,A4,A1,I1,A1,I1,A1,I1,A7,I2,A1/)')&
                 &              WORD, '(', I, ',', J, ', ',  T, ',iPass=',iP,')'
            WRITE (LUPRI,'(1X,6ES16.8)') (ETIJ(T,I,J,X,K,IP),K=1,nPrimP)
         END DO
      END DO
   END DO
ENDDO
END SUBROUTINE PRINT_ETIJ

subroutine build_auxiliary(nPrimA,nPrimB,Aexp,Bexp,PINV,APINV,BPINV,HPINV,TWOA,TWOB)
    implicit none
    integer,intent(in) :: nPrimA,nPrimB
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)
    real(realk),intent(inout) :: APINV(nPrimA,nPrimB),BPINV(nPrimA,nPrimB)
    real(realk),intent(inout) :: HPINV(nPrimA,nPrimB),TWOA(nPrimA,nPrimB)
    real(realk),intent(inout) :: TWOB(nPrimA,nPrimB),PINV(nPrimA,nPrimB)
    real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk,D1=1.0E0_realk
    !
    real(realk) :: bexpo
    integer :: J,I
    DO J=1,nPrimB
       bexpo = Bexp(J)
       DO I=1,nPrimA
          PINV(I,J)  = D1/(Aexp(I)+bexpo)   ! 1/p
          HPINV(I,J) = DHALF   * PINV(I,J)  ! 1/2p
          BPINV(I,J) = bexpo   * PINV(I,J)  ! b/p
          APINV(I,J) = Aexp(I) * PINV(I,J)  ! a/p
       ENDDO
       bexpo = D2*bexpo
       DO I=1,nPrimA
          TWOB(I,J)  = bexpo                ! 2b
       ENDDO
    ENDDO
    DO J=1,nPrimB
       DO I=1,nPrimA
          TWOA(I,J)  = D2*Aexp(I)           ! 2a
       ENDDO
    ENDDO
end subroutine build_auxiliary

SUBROUTINE ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPass,MAXI,MAXJ,PA,PB,&
     & HPINV,TWOA,TWOB,X,LUPRI)
implicit none
INTEGER,intent(in)        :: nPrimP,nPass,MAXI,MAXJ,LUPRI,X
REAL(REALK),intent(in)    :: PA(nPrimP,nPass), PB(nPrimP,nPass)
REAL(REALK),intent(inout) :: ETIJ(0:MAXI+MAXJ,0:MAXI,0:MAXJ,3,nPrimP,nPass)
REAL(REALK),intent(in)    :: HPINV(nPrimP)
REAL(REALK)               :: TWOA(nPrimP),TWOB(nPrimP)
!
REAL(REALK),PARAMETER     :: D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER                   :: MAXIJ,iPass,K
MAXIJ = MAXI + MAXJ
IF(maxIJ.EQ.0 )THEN
   IF(maxI.EQ.0 )THEN
      !maxI.EQ.0.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.1 )THEN
   IF(maxI.EQ.1 )THEN
      !maxI.EQ.1.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.2 )THEN
   IF(maxI.EQ.2 )THEN
      !maxI.EQ.2.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.3 )THEN
   IF(maxI.EQ.3 )THEN
      !maxI.EQ.3.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.2)THEN
      !maxI.EQ.2.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.3
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.4 )THEN
   IF(maxI.EQ.4 )THEN
      !maxI.EQ.4.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,4,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,3,0,X,K,ipass)+ ETIJ(1,3,0,X,K,ipass)&
                 & - 3*ETIJ( 0,2,0,X,K,ipass)/TWOA(K)
            ETIJ(3,4,0,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,4,0,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,4,0,X,K,ipass) = HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass) - 3*ETIJ( 1,2,0,X,K,ipass)/TWOA(K)
            ETIJ(2,4,0,X,K,ipass) = HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass) - 3*ETIJ( 2,2,0,X,K,ipass)/TWOA(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.3)THEN
      !maxI.EQ.3.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,3,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,0,X,K,ipass)+ ETIJ( 1,3,0,X,K,ipass)
            ETIJ(3,3,1,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,3,1,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,3,1,X,K,ipass)= HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass)! - 0*ETIJ(1,3,*,X,K,ipass)/TWOB(K)
            ETIJ(2,3,1,X,K,ipass)= HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass)! - 0*ETIJ(2,3,*,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.2)THEN
      !maxI.EQ.2.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)

            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ(0,1,1,X,K,ipass)+ ETIJ(1,1,1,X,K,ipass)&
                                & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass) = HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                                & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                                & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,2,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,1,X,K,ipass)+ ETIJ( 1,2,1,X,K,ipass)&
                                & -1*ETIJ(0,2,0,X,K,ipass)/TWOB(K)
            ETIJ(3,2,2,X,K,ipass) = HPINV(K)*ETIJ(2,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(4,2,2,X,K,ipass) = HPINV(K)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(1,2,2,X,K,ipass) = HPINV(K)*ETIJ(0,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,1,X,K,ipass)&
                                & + 2*ETIJ(2,2,1,X,K,ipass) - 1*ETIJ(1,2,0,X,K,ipass)/TWOB(K)
            ETIJ(2,2,2,X,K,ipass) = HPINV(K)*ETIJ(1,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,1,X,K,ipass)&
                                & + 3*ETIJ(3,2,1,X,K,ipass) - 1*ETIJ(2,2,0,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.3
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,1,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,2,X,K,ipass)+ ETIJ( 1,1,2,X,K,ipass)&
                 & -2*ETIJ(0,1,1,X,K,ipass)/TWOB(K)
            ETIJ(3,1,3,X,K,ipass) = HPINV(K)*ETIJ(2,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(4,1,3,X,K,ipass) = HPINV(K)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(1,1,3,X,K,ipass)= HPINV(K)*ETIJ(0,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,2,X,K,ipass)&
                 & + 2*ETIJ(2,1,2,X,K,ipass) - 2*ETIJ(1,1,1,X,K,ipass)/TWOB(K)
            ETIJ(2,1,3,X,K,ipass)= HPINV(K)*ETIJ(1,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,2,X,K,ipass)&
                 & + 3*ETIJ(3,1,2,X,K,ipass) - 2*ETIJ(2,1,1,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.4
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,0,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,3,X,K,ipass)+ ETIJ( 1,0,3,X,K,ipass)&
                 & -3*ETIJ(0,0,2,X,K,ipass)/TWOB(K)
            ETIJ(3,0,4,X,K,ipass) = HPINV(K)*ETIJ(2,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(4,0,4,X,K,ipass) = HPINV(K)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(1,0,4,X,K,ipass)= HPINV(K)*ETIJ(0,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,3,X,K,ipass)&
                 & + 2*ETIJ(2,0,3,X,K,ipass) - 3*ETIJ(1,0,2,X,K,ipass)/TWOB(K)
            ETIJ(2,0,4,X,K,ipass)= HPINV(K)*ETIJ(1,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,3,X,K,ipass)&
                 & + 3*ETIJ(3,0,3,X,K,ipass) - 3*ETIJ(2,0,2,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.5 )THEN
   IF(maxI.EQ.5 )THEN
      !maxI.EQ.5.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,4,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,3,0,X,K,ipass)+ ETIJ(1,3,0,X,K,ipass)&
                 & - 3*ETIJ( 0,2,0,X,K,ipass)/TWOA(K)
            ETIJ(3,4,0,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,4,0,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,4,0,X,K,ipass) = HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass) - 3*ETIJ( 1,2,0,X,K,ipass)/TWOA(K)
            ETIJ(2,4,0,X,K,ipass) = HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass) - 3*ETIJ( 2,2,0,X,K,ipass)/TWOA(K)
            ETIJ(0,5,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,4,0,X,K,ipass)+ ETIJ(1,4,0,X,K,ipass)&
                 & - 4*ETIJ( 0,3,0,X,K,ipass)/TWOA(K)
            ETIJ(4,5,0,X,K,ipass) = HPINV(K)*ETIJ(3,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(5,5,0,X,K,ipass) = HPINV(K)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(1,5,0,X,K,ipass) = HPINV(K)*ETIJ(0,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,4,0,X,K,ipass)&
                 & + 2*ETIJ(2,4,0,X,K,ipass) - 4*ETIJ( 1,3,0,X,K,ipass)/TWOA(K)
            ETIJ(2,5,0,X,K,ipass) = HPINV(K)*ETIJ(1,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,4,0,X,K,ipass)&
                 & + 3*ETIJ(3,4,0,X,K,ipass) - 4*ETIJ( 2,3,0,X,K,ipass)/TWOA(K)
            ETIJ(3,5,0,X,K,ipass) = HPINV(K)*ETIJ(2,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,4,0,X,K,ipass)&
                 & + 4*ETIJ(4,4,0,X,K,ipass) - 4*ETIJ( 3,3,0,X,K,ipass)/TWOA(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.4)THEN
      !maxI.EQ.4.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,3,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,0,X,K,ipass)+ ETIJ( 1,3,0,X,K,ipass)
            ETIJ(3,3,1,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,3,1,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,3,1,X,K,ipass)= HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass)! - 0*ETIJ(1,3,*,X,K,ipass)/TWOB(K)
            ETIJ(2,3,1,X,K,ipass)= HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass)! - 0*ETIJ(2,3,*,X,K,ipass)/TWOB(K)
            ETIJ(0,4,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,3,0,X,K,ipass)+ ETIJ(1,3,0,X,K,ipass)&
                 & - 3*ETIJ( 0,2,0,X,K,ipass)/TWOA(K)
            ETIJ(3,4,0,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,4,0,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,4,0,X,K,ipass) = HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass) - 3*ETIJ( 1,2,0,X,K,ipass)/TWOA(K)
            ETIJ(2,4,0,X,K,ipass) = HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass) - 3*ETIJ( 2,2,0,X,K,ipass)/TWOA(K)
            ETIJ(0,4,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,4,0,X,K,ipass)+ ETIJ( 1,4,0,X,K,ipass)
            ETIJ(4,4,1,X,K,ipass) = HPINV(K)*ETIJ(3,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(5,4,1,X,K,ipass) = HPINV(K)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(1,4,1,X,K,ipass)= HPINV(K)*ETIJ(0,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,4,0,X,K,ipass)&
                 & + 2*ETIJ(2,4,0,X,K,ipass)! - 0*ETIJ(1,4,*,X,K,ipass)/TWOB(K)
            ETIJ(2,4,1,X,K,ipass)= HPINV(K)*ETIJ(1,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,4,0,X,K,ipass)&
                 & + 3*ETIJ(3,4,0,X,K,ipass)! - 0*ETIJ(2,4,*,X,K,ipass)/TWOB(K)
            ETIJ(3,4,1,X,K,ipass)= HPINV(K)*ETIJ(2,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,4,0,X,K,ipass)&
                 & + 4*ETIJ(4,4,0,X,K,ipass)! - 0*ETIJ(3,4,*,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.3)THEN
      !maxI.EQ.3.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,2,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,1,X,K,ipass)+ ETIJ( 1,2,1,X,K,ipass)&
                 & -1*ETIJ(0,2,0,X,K,ipass)/TWOB(K)
            ETIJ(3,2,2,X,K,ipass) = HPINV(K)*ETIJ(2,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(4,2,2,X,K,ipass) = HPINV(K)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(1,2,2,X,K,ipass)= HPINV(K)*ETIJ(0,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,1,X,K,ipass)&
                 & + 2*ETIJ(2,2,1,X,K,ipass) - 1*ETIJ(1,2,0,X,K,ipass)/TWOB(K)
            ETIJ(2,2,2,X,K,ipass)= HPINV(K)*ETIJ(1,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,1,X,K,ipass)&
                 & + 3*ETIJ(3,2,1,X,K,ipass) - 1*ETIJ(2,2,0,X,K,ipass)/TWOB(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,3,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,0,X,K,ipass)+ ETIJ( 1,3,0,X,K,ipass)
            ETIJ(3,3,1,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,3,1,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,3,1,X,K,ipass)= HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass)! - 0*ETIJ(1,3,*,X,K,ipass)/TWOB(K)
            ETIJ(2,3,1,X,K,ipass)= HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass)! - 0*ETIJ(2,3,*,X,K,ipass)/TWOB(K)
            ETIJ(0,3,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,1,X,K,ipass)+ ETIJ( 1,3,1,X,K,ipass)&
                 & -1*ETIJ(0,3,0,X,K,ipass)/TWOB(K)
            ETIJ(4,3,2,X,K,ipass) = HPINV(K)*ETIJ(3,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(4,3,1,X,K,ipass)
            ETIJ(5,3,2,X,K,ipass) = HPINV(K)*ETIJ(4,3,1,X,K,ipass)
            ETIJ(1,3,2,X,K,ipass)= HPINV(K)*ETIJ(0,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,1,X,K,ipass)&
                 & + 2*ETIJ(2,3,1,X,K,ipass) - 1*ETIJ(1,3,0,X,K,ipass)/TWOB(K)
            ETIJ(2,3,2,X,K,ipass)= HPINV(K)*ETIJ(1,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,1,X,K,ipass)&
                 & + 3*ETIJ(3,3,1,X,K,ipass) - 1*ETIJ(2,3,0,X,K,ipass)/TWOB(K)
            ETIJ(3,3,2,X,K,ipass)= HPINV(K)*ETIJ(2,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,1,X,K,ipass)&
                 & + 4*ETIJ(4,3,1,X,K,ipass) - 1*ETIJ(3,3,0,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.2)THEN
      !maxI.EQ.2.AND.maxJ.EQ.3
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,1,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,2,X,K,ipass)+ ETIJ( 1,1,2,X,K,ipass)&
                 & -2*ETIJ(0,1,1,X,K,ipass)/TWOB(K)
            ETIJ(3,1,3,X,K,ipass) = HPINV(K)*ETIJ(2,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(4,1,3,X,K,ipass) = HPINV(K)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(1,1,3,X,K,ipass)= HPINV(K)*ETIJ(0,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,2,X,K,ipass)&
                 & + 2*ETIJ(2,1,2,X,K,ipass) - 2*ETIJ(1,1,1,X,K,ipass)/TWOB(K)
            ETIJ(2,1,3,X,K,ipass)= HPINV(K)*ETIJ(1,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,2,X,K,ipass)&
                 & + 3*ETIJ(3,1,2,X,K,ipass) - 2*ETIJ(2,1,1,X,K,ipass)/TWOB(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,2,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,1,X,K,ipass)+ ETIJ( 1,2,1,X,K,ipass)&
                 & -1*ETIJ(0,2,0,X,K,ipass)/TWOB(K)
            ETIJ(3,2,2,X,K,ipass) = HPINV(K)*ETIJ(2,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(4,2,2,X,K,ipass) = HPINV(K)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(1,2,2,X,K,ipass)= HPINV(K)*ETIJ(0,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,1,X,K,ipass)&
                 & + 2*ETIJ(2,2,1,X,K,ipass) - 1*ETIJ(1,2,0,X,K,ipass)/TWOB(K)
            ETIJ(2,2,2,X,K,ipass)= HPINV(K)*ETIJ(1,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,1,X,K,ipass)&
                 & + 3*ETIJ(3,2,1,X,K,ipass) - 1*ETIJ(2,2,0,X,K,ipass)/TWOB(K)
            ETIJ(0,2,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,2,X,K,ipass)+ ETIJ( 1,2,2,X,K,ipass)&
                 & -2*ETIJ(0,2,1,X,K,ipass)/TWOB(K)
            ETIJ(4,2,3,X,K,ipass) = HPINV(K)*ETIJ(3,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(4,2,2,X,K,ipass)
            ETIJ(5,2,3,X,K,ipass) = HPINV(K)*ETIJ(4,2,2,X,K,ipass)
            ETIJ(1,2,3,X,K,ipass)= HPINV(K)*ETIJ(0,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,2,X,K,ipass)&
                 & + 2*ETIJ(2,2,2,X,K,ipass) - 2*ETIJ(1,2,1,X,K,ipass)/TWOB(K)
            ETIJ(2,2,3,X,K,ipass)= HPINV(K)*ETIJ(1,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,2,X,K,ipass)&
                 & + 3*ETIJ(3,2,2,X,K,ipass) - 2*ETIJ(2,2,1,X,K,ipass)/TWOB(K)
            ETIJ(3,2,3,X,K,ipass)= HPINV(K)*ETIJ(2,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,2,X,K,ipass)&
                 & + 4*ETIJ(4,2,2,X,K,ipass) - 2*ETIJ(3,2,1,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.4
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,0,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,3,X,K,ipass)+ ETIJ( 1,0,3,X,K,ipass)&
                 & -3*ETIJ(0,0,2,X,K,ipass)/TWOB(K)
            ETIJ(3,0,4,X,K,ipass) = HPINV(K)*ETIJ(2,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(4,0,4,X,K,ipass) = HPINV(K)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(1,0,4,X,K,ipass)= HPINV(K)*ETIJ(0,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,3,X,K,ipass)&
                 & + 2*ETIJ(2,0,3,X,K,ipass) - 3*ETIJ(1,0,2,X,K,ipass)/TWOB(K)
            ETIJ(2,0,4,X,K,ipass)= HPINV(K)*ETIJ(1,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,3,X,K,ipass)&
                 & + 3*ETIJ(3,0,3,X,K,ipass) - 3*ETIJ(2,0,2,X,K,ipass)/TWOB(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,1,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,2,X,K,ipass)+ ETIJ( 1,1,2,X,K,ipass)&
                 & -2*ETIJ(0,1,1,X,K,ipass)/TWOB(K)
            ETIJ(3,1,3,X,K,ipass) = HPINV(K)*ETIJ(2,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(4,1,3,X,K,ipass) = HPINV(K)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(1,1,3,X,K,ipass)= HPINV(K)*ETIJ(0,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,2,X,K,ipass)&
                 & + 2*ETIJ(2,1,2,X,K,ipass) - 2*ETIJ(1,1,1,X,K,ipass)/TWOB(K)
            ETIJ(2,1,3,X,K,ipass)= HPINV(K)*ETIJ(1,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,2,X,K,ipass)&
                 & + 3*ETIJ(3,1,2,X,K,ipass) - 2*ETIJ(2,1,1,X,K,ipass)/TWOB(K)
            ETIJ(0,1,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,3,X,K,ipass)+ ETIJ( 1,1,3,X,K,ipass)&
                 & -3*ETIJ(0,1,2,X,K,ipass)/TWOB(K)
            ETIJ(4,1,4,X,K,ipass) = HPINV(K)*ETIJ(3,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(4,1,3,X,K,ipass)
            ETIJ(5,1,4,X,K,ipass) = HPINV(K)*ETIJ(4,1,3,X,K,ipass)
            ETIJ(1,1,4,X,K,ipass)= HPINV(K)*ETIJ(0,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,3,X,K,ipass)&
                 & + 2*ETIJ(2,1,3,X,K,ipass) - 3*ETIJ(1,1,2,X,K,ipass)/TWOB(K)
            ETIJ(2,1,4,X,K,ipass)= HPINV(K)*ETIJ(1,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,3,X,K,ipass)&
                 & + 3*ETIJ(3,1,3,X,K,ipass) - 3*ETIJ(2,1,2,X,K,ipass)/TWOB(K)
            ETIJ(3,1,4,X,K,ipass)= HPINV(K)*ETIJ(2,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,3,X,K,ipass)&
                 & + 4*ETIJ(4,1,3,X,K,ipass) - 3*ETIJ(3,1,2,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.5
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,0,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,3,X,K,ipass)+ ETIJ( 1,0,3,X,K,ipass)&
                 & -3*ETIJ(0,0,2,X,K,ipass)/TWOB(K)
            ETIJ(3,0,4,X,K,ipass) = HPINV(K)*ETIJ(2,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(4,0,4,X,K,ipass) = HPINV(K)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(1,0,4,X,K,ipass)= HPINV(K)*ETIJ(0,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,3,X,K,ipass)&
                 & + 2*ETIJ(2,0,3,X,K,ipass) - 3*ETIJ(1,0,2,X,K,ipass)/TWOB(K)
            ETIJ(2,0,4,X,K,ipass)= HPINV(K)*ETIJ(1,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,3,X,K,ipass)&
                 & + 3*ETIJ(3,0,3,X,K,ipass) - 3*ETIJ(2,0,2,X,K,ipass)/TWOB(K)
            ETIJ(0,0,5,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,4,X,K,ipass)+ ETIJ( 1,0,4,X,K,ipass)&
                 & -4*ETIJ(0,0,3,X,K,ipass)/TWOB(K)
            ETIJ(4,0,5,X,K,ipass) = HPINV(K)*ETIJ(3,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(4,0,4,X,K,ipass)
            ETIJ(5,0,5,X,K,ipass) = HPINV(K)*ETIJ(4,0,4,X,K,ipass)
            ETIJ(1,0,5,X,K,ipass)= HPINV(K)*ETIJ(0,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,4,X,K,ipass)&
                 & + 2*ETIJ(2,0,4,X,K,ipass) - 4*ETIJ(1,0,3,X,K,ipass)/TWOB(K)
            ETIJ(2,0,5,X,K,ipass)= HPINV(K)*ETIJ(1,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,4,X,K,ipass)&
                 & + 3*ETIJ(3,0,4,X,K,ipass) - 4*ETIJ(2,0,3,X,K,ipass)/TWOB(K)
            ETIJ(3,0,5,X,K,ipass)= HPINV(K)*ETIJ(2,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,4,X,K,ipass)&
                 & + 4*ETIJ(4,0,4,X,K,ipass) - 4*ETIJ(3,0,3,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.6 )THEN
   IF(maxI.EQ.6 )THEN
      !maxI.EQ.6.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,4,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,3,0,X,K,ipass)+ ETIJ(1,3,0,X,K,ipass)&
                 & - 3*ETIJ( 0,2,0,X,K,ipass)/TWOA(K)
            ETIJ(3,4,0,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,4,0,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,4,0,X,K,ipass) = HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass) - 3*ETIJ( 1,2,0,X,K,ipass)/TWOA(K)
            ETIJ(2,4,0,X,K,ipass) = HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass) - 3*ETIJ( 2,2,0,X,K,ipass)/TWOA(K)
            ETIJ(0,5,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,4,0,X,K,ipass)+ ETIJ(1,4,0,X,K,ipass)&
                 & - 4*ETIJ( 0,3,0,X,K,ipass)/TWOA(K)
            ETIJ(4,5,0,X,K,ipass) = HPINV(K)*ETIJ(3,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(5,5,0,X,K,ipass) = HPINV(K)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(1,5,0,X,K,ipass) = HPINV(K)*ETIJ(0,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,4,0,X,K,ipass)&
                 & + 2*ETIJ(2,4,0,X,K,ipass) - 4*ETIJ( 1,3,0,X,K,ipass)/TWOA(K)
            ETIJ(2,5,0,X,K,ipass) = HPINV(K)*ETIJ(1,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,4,0,X,K,ipass)&
                 & + 3*ETIJ(3,4,0,X,K,ipass) - 4*ETIJ( 2,3,0,X,K,ipass)/TWOA(K)
            ETIJ(3,5,0,X,K,ipass) = HPINV(K)*ETIJ(2,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,4,0,X,K,ipass)&
                 & + 4*ETIJ(4,4,0,X,K,ipass) - 4*ETIJ( 3,3,0,X,K,ipass)/TWOA(K)
            ETIJ(0,6,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,5,0,X,K,ipass)+ ETIJ(1,5,0,X,K,ipass)&
                 & - 5*ETIJ( 0,4,0,X,K,ipass)/TWOA(K)
            ETIJ(5,6,0,X,K,ipass) = HPINV(K)*ETIJ(4,5,0,X,K,ipass)+ PA(K,ipass)*ETIJ(5,5,0,X,K,ipass)
            ETIJ(6,6,0,X,K,ipass) = HPINV(K)*ETIJ(5,5,0,X,K,ipass)
            ETIJ(1,6,0,X,K,ipass) = HPINV(K)*ETIJ(0,5,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,5,0,X,K,ipass)&
                 & + 2*ETIJ(2,5,0,X,K,ipass) - 5*ETIJ( 1,4,0,X,K,ipass)/TWOA(K)
            ETIJ(2,6,0,X,K,ipass) = HPINV(K)*ETIJ(1,5,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,5,0,X,K,ipass)&
                 & + 3*ETIJ(3,5,0,X,K,ipass) - 5*ETIJ( 2,4,0,X,K,ipass)/TWOA(K)
            ETIJ(3,6,0,X,K,ipass) = HPINV(K)*ETIJ(2,5,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,5,0,X,K,ipass)&
                 & + 4*ETIJ(4,5,0,X,K,ipass) - 5*ETIJ( 3,4,0,X,K,ipass)/TWOA(K)
            ETIJ(4,6,0,X,K,ipass) = HPINV(K)*ETIJ(3,5,0,X,K,ipass)+ PA(K,ipass)*ETIJ(4,5,0,X,K,ipass)&
                 & + 5*ETIJ(5,5,0,X,K,ipass) - 5*ETIJ( 4,4,0,X,K,ipass)/TWOA(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.5)THEN
      !maxI.EQ.5.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,3,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,0,X,K,ipass)+ ETIJ( 1,3,0,X,K,ipass)
            ETIJ(3,3,1,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,3,1,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,3,1,X,K,ipass)= HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass)! - 0*ETIJ(1,3,*,X,K,ipass)/TWOB(K)
            ETIJ(2,3,1,X,K,ipass)= HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass)! - 0*ETIJ(2,3,*,X,K,ipass)/TWOB(K)
            ETIJ(0,4,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,3,0,X,K,ipass)+ ETIJ(1,3,0,X,K,ipass)&
                 & - 3*ETIJ( 0,2,0,X,K,ipass)/TWOA(K)
            ETIJ(3,4,0,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,4,0,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,4,0,X,K,ipass) = HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass) - 3*ETIJ( 1,2,0,X,K,ipass)/TWOA(K)
            ETIJ(2,4,0,X,K,ipass) = HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass) - 3*ETIJ( 2,2,0,X,K,ipass)/TWOA(K)
            ETIJ(0,4,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,4,0,X,K,ipass)+ ETIJ( 1,4,0,X,K,ipass)
            ETIJ(4,4,1,X,K,ipass) = HPINV(K)*ETIJ(3,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(5,4,1,X,K,ipass) = HPINV(K)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(1,4,1,X,K,ipass)= HPINV(K)*ETIJ(0,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,4,0,X,K,ipass)&
                 & + 2*ETIJ(2,4,0,X,K,ipass)! - 0*ETIJ(1,4,*,X,K,ipass)/TWOB(K)
            ETIJ(2,4,1,X,K,ipass)= HPINV(K)*ETIJ(1,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,4,0,X,K,ipass)&
                 & + 3*ETIJ(3,4,0,X,K,ipass)! - 0*ETIJ(2,4,*,X,K,ipass)/TWOB(K)
            ETIJ(3,4,1,X,K,ipass)= HPINV(K)*ETIJ(2,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,4,0,X,K,ipass)&
                 & + 4*ETIJ(4,4,0,X,K,ipass)! - 0*ETIJ(3,4,*,X,K,ipass)/TWOB(K)
            ETIJ(0,5,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,4,0,X,K,ipass)+ ETIJ(1,4,0,X,K,ipass)&
                 & - 4*ETIJ( 0,3,0,X,K,ipass)/TWOA(K)
            ETIJ(4,5,0,X,K,ipass) = HPINV(K)*ETIJ(3,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(5,5,0,X,K,ipass) = HPINV(K)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(1,5,0,X,K,ipass) = HPINV(K)*ETIJ(0,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,4,0,X,K,ipass)&
                 & + 2*ETIJ(2,4,0,X,K,ipass) - 4*ETIJ( 1,3,0,X,K,ipass)/TWOA(K)
            ETIJ(2,5,0,X,K,ipass) = HPINV(K)*ETIJ(1,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,4,0,X,K,ipass)&
                 & + 3*ETIJ(3,4,0,X,K,ipass) - 4*ETIJ( 2,3,0,X,K,ipass)/TWOA(K)
            ETIJ(3,5,0,X,K,ipass) = HPINV(K)*ETIJ(2,4,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,4,0,X,K,ipass)&
                 & + 4*ETIJ(4,4,0,X,K,ipass) - 4*ETIJ( 3,3,0,X,K,ipass)/TWOA(K)
            ETIJ(0,5,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,5,0,X,K,ipass)+ ETIJ( 1,5,0,X,K,ipass)
            ETIJ(5,5,1,X,K,ipass) = HPINV(K)*ETIJ(4,5,0,X,K,ipass)+ PB(K,ipass)*ETIJ(5,5,0,X,K,ipass)
            ETIJ(6,5,1,X,K,ipass) = HPINV(K)*ETIJ(5,5,0,X,K,ipass)
            ETIJ(1,5,1,X,K,ipass)= HPINV(K)*ETIJ(0,5,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,5,0,X,K,ipass)&
                 & + 2*ETIJ(2,5,0,X,K,ipass)! - 0*ETIJ(1,5,*,X,K,ipass)/TWOB(K)
            ETIJ(2,5,1,X,K,ipass)= HPINV(K)*ETIJ(1,5,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,5,0,X,K,ipass)&
                 & + 3*ETIJ(3,5,0,X,K,ipass)! - 0*ETIJ(2,5,*,X,K,ipass)/TWOB(K)
            ETIJ(3,5,1,X,K,ipass)= HPINV(K)*ETIJ(2,5,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,5,0,X,K,ipass)&
                 & + 4*ETIJ(4,5,0,X,K,ipass)! - 0*ETIJ(3,5,*,X,K,ipass)/TWOB(K)
            ETIJ(4,5,1,X,K,ipass)= HPINV(K)*ETIJ(3,5,0,X,K,ipass)+ PB(K,ipass)*ETIJ(4,5,0,X,K,ipass)&
                 & + 5*ETIJ(5,5,0,X,K,ipass)! - 0*ETIJ(4,5,*,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.4)THEN
      !maxI.EQ.4.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,2,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,1,X,K,ipass)+ ETIJ( 1,2,1,X,K,ipass)&
                 & -1*ETIJ(0,2,0,X,K,ipass)/TWOB(K)
            ETIJ(3,2,2,X,K,ipass) = HPINV(K)*ETIJ(2,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(4,2,2,X,K,ipass) = HPINV(K)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(1,2,2,X,K,ipass)= HPINV(K)*ETIJ(0,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,1,X,K,ipass)&
                 & + 2*ETIJ(2,2,1,X,K,ipass) - 1*ETIJ(1,2,0,X,K,ipass)/TWOB(K)
            ETIJ(2,2,2,X,K,ipass)= HPINV(K)*ETIJ(1,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,1,X,K,ipass)&
                 & + 3*ETIJ(3,2,1,X,K,ipass) - 1*ETIJ(2,2,0,X,K,ipass)/TWOB(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,3,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,0,X,K,ipass)+ ETIJ( 1,3,0,X,K,ipass)
            ETIJ(3,3,1,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,3,1,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,3,1,X,K,ipass)= HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass)! - 0*ETIJ(1,3,*,X,K,ipass)/TWOB(K)
            ETIJ(2,3,1,X,K,ipass)= HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass)! - 0*ETIJ(2,3,*,X,K,ipass)/TWOB(K)
            ETIJ(0,3,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,1,X,K,ipass)+ ETIJ( 1,3,1,X,K,ipass)&
                 & -1*ETIJ(0,3,0,X,K,ipass)/TWOB(K)
            ETIJ(4,3,2,X,K,ipass) = HPINV(K)*ETIJ(3,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(4,3,1,X,K,ipass)
            ETIJ(5,3,2,X,K,ipass) = HPINV(K)*ETIJ(4,3,1,X,K,ipass)
            ETIJ(1,3,2,X,K,ipass)= HPINV(K)*ETIJ(0,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,1,X,K,ipass)&
                 & + 2*ETIJ(2,3,1,X,K,ipass) - 1*ETIJ(1,3,0,X,K,ipass)/TWOB(K)
            ETIJ(2,3,2,X,K,ipass)= HPINV(K)*ETIJ(1,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,1,X,K,ipass)&
                 & + 3*ETIJ(3,3,1,X,K,ipass) - 1*ETIJ(2,3,0,X,K,ipass)/TWOB(K)
            ETIJ(3,3,2,X,K,ipass)= HPINV(K)*ETIJ(2,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,1,X,K,ipass)&
                 & + 4*ETIJ(4,3,1,X,K,ipass) - 1*ETIJ(3,3,0,X,K,ipass)/TWOB(K)
            ETIJ(0,4,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,3,0,X,K,ipass)+ ETIJ(1,3,0,X,K,ipass)&
                 & - 3*ETIJ( 0,2,0,X,K,ipass)/TWOA(K)
            ETIJ(3,4,0,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,4,0,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,4,0,X,K,ipass) = HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass) - 3*ETIJ( 1,2,0,X,K,ipass)/TWOA(K)
            ETIJ(2,4,0,X,K,ipass) = HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass) - 3*ETIJ( 2,2,0,X,K,ipass)/TWOA(K)
            ETIJ(0,4,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,4,0,X,K,ipass)+ ETIJ( 1,4,0,X,K,ipass)
            ETIJ(4,4,1,X,K,ipass) = HPINV(K)*ETIJ(3,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(5,4,1,X,K,ipass) = HPINV(K)*ETIJ(4,4,0,X,K,ipass)
            ETIJ(1,4,1,X,K,ipass)= HPINV(K)*ETIJ(0,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,4,0,X,K,ipass)&
                 & + 2*ETIJ(2,4,0,X,K,ipass)! - 0*ETIJ(1,4,*,X,K,ipass)/TWOB(K)
            ETIJ(2,4,1,X,K,ipass)= HPINV(K)*ETIJ(1,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,4,0,X,K,ipass)&
                 & + 3*ETIJ(3,4,0,X,K,ipass)! - 0*ETIJ(2,4,*,X,K,ipass)/TWOB(K)
            ETIJ(3,4,1,X,K,ipass)= HPINV(K)*ETIJ(2,4,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,4,0,X,K,ipass)&
                 & + 4*ETIJ(4,4,0,X,K,ipass)! - 0*ETIJ(3,4,*,X,K,ipass)/TWOB(K)
            ETIJ(0,4,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,4,1,X,K,ipass)+ ETIJ( 1,4,1,X,K,ipass)&
                 & -1*ETIJ(0,4,0,X,K,ipass)/TWOB(K)
            ETIJ(5,4,2,X,K,ipass) = HPINV(K)*ETIJ(4,4,1,X,K,ipass)+ PB(K,ipass)*ETIJ(5,4,1,X,K,ipass)
            ETIJ(6,4,2,X,K,ipass) = HPINV(K)*ETIJ(5,4,1,X,K,ipass)
            ETIJ(1,4,2,X,K,ipass)= HPINV(K)*ETIJ(0,4,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,4,1,X,K,ipass)&
                 & + 2*ETIJ(2,4,1,X,K,ipass) - 1*ETIJ(1,4,0,X,K,ipass)/TWOB(K)
            ETIJ(2,4,2,X,K,ipass)= HPINV(K)*ETIJ(1,4,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,4,1,X,K,ipass)&
                 & + 3*ETIJ(3,4,1,X,K,ipass) - 1*ETIJ(2,4,0,X,K,ipass)/TWOB(K)
            ETIJ(3,4,2,X,K,ipass)= HPINV(K)*ETIJ(2,4,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,4,1,X,K,ipass)&
                 & + 4*ETIJ(4,4,1,X,K,ipass) - 1*ETIJ(3,4,0,X,K,ipass)/TWOB(K)
            ETIJ(4,4,2,X,K,ipass)= HPINV(K)*ETIJ(3,4,1,X,K,ipass)+ PB(K,ipass)*ETIJ(4,4,1,X,K,ipass)&
                 & + 5*ETIJ(5,4,1,X,K,ipass) - 1*ETIJ(4,4,0,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.3)THEN
      !maxI.EQ.3.AND.maxJ.EQ.3
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,1,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,2,X,K,ipass)+ ETIJ( 1,1,2,X,K,ipass)&
                 & -2*ETIJ(0,1,1,X,K,ipass)/TWOB(K)
            ETIJ(3,1,3,X,K,ipass) = HPINV(K)*ETIJ(2,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(4,1,3,X,K,ipass) = HPINV(K)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(1,1,3,X,K,ipass)= HPINV(K)*ETIJ(0,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,2,X,K,ipass)&
                 & + 2*ETIJ(2,1,2,X,K,ipass) - 2*ETIJ(1,1,1,X,K,ipass)/TWOB(K)
            ETIJ(2,1,3,X,K,ipass)= HPINV(K)*ETIJ(1,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,2,X,K,ipass)&
                 & + 3*ETIJ(3,1,2,X,K,ipass) - 2*ETIJ(2,1,1,X,K,ipass)/TWOB(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,2,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,1,X,K,ipass)+ ETIJ( 1,2,1,X,K,ipass)&
                 & -1*ETIJ(0,2,0,X,K,ipass)/TWOB(K)
            ETIJ(3,2,2,X,K,ipass) = HPINV(K)*ETIJ(2,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(4,2,2,X,K,ipass) = HPINV(K)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(1,2,2,X,K,ipass)= HPINV(K)*ETIJ(0,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,1,X,K,ipass)&
                 & + 2*ETIJ(2,2,1,X,K,ipass) - 1*ETIJ(1,2,0,X,K,ipass)/TWOB(K)
            ETIJ(2,2,2,X,K,ipass)= HPINV(K)*ETIJ(1,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,1,X,K,ipass)&
                 & + 3*ETIJ(3,2,1,X,K,ipass) - 1*ETIJ(2,2,0,X,K,ipass)/TWOB(K)
            ETIJ(0,2,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,2,X,K,ipass)+ ETIJ( 1,2,2,X,K,ipass)&
                 & -2*ETIJ(0,2,1,X,K,ipass)/TWOB(K)
            ETIJ(4,2,3,X,K,ipass) = HPINV(K)*ETIJ(3,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(4,2,2,X,K,ipass)
            ETIJ(5,2,3,X,K,ipass) = HPINV(K)*ETIJ(4,2,2,X,K,ipass)
            ETIJ(1,2,3,X,K,ipass)= HPINV(K)*ETIJ(0,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,2,X,K,ipass)&
                 & + 2*ETIJ(2,2,2,X,K,ipass) - 2*ETIJ(1,2,1,X,K,ipass)/TWOB(K)
            ETIJ(2,2,3,X,K,ipass)= HPINV(K)*ETIJ(1,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,2,X,K,ipass)&
                 & + 3*ETIJ(3,2,2,X,K,ipass) - 2*ETIJ(2,2,1,X,K,ipass)/TWOB(K)
            ETIJ(3,2,3,X,K,ipass)= HPINV(K)*ETIJ(2,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,2,X,K,ipass)&
                 & + 4*ETIJ(4,2,2,X,K,ipass) - 2*ETIJ(3,2,1,X,K,ipass)/TWOB(K)
            ETIJ(0,3,0,X,K,ipass) = PA(K,ipass)*ETIJ(0,2,0,X,K,ipass)+ ETIJ(1,2,0,X,K,ipass)&
                 & - 2*ETIJ( 0,1,0,X,K,ipass)/TWOA(K)
            ETIJ(2,3,0,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,3,0,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,3,0,X,K,ipass) = HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PA(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass) - 2*ETIJ( 1,1,0,X,K,ipass)/TWOA(K)
            ETIJ(0,3,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,0,X,K,ipass)+ ETIJ( 1,3,0,X,K,ipass)
            ETIJ(3,3,1,X,K,ipass) = HPINV(K)*ETIJ(2,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(4,3,1,X,K,ipass) = HPINV(K)*ETIJ(3,3,0,X,K,ipass)
            ETIJ(1,3,1,X,K,ipass)= HPINV(K)*ETIJ(0,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,0,X,K,ipass)&
                 & + 2*ETIJ(2,3,0,X,K,ipass)! - 0*ETIJ(1,3,*,X,K,ipass)/TWOB(K)
            ETIJ(2,3,1,X,K,ipass)= HPINV(K)*ETIJ(1,3,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,0,X,K,ipass)&
                 & + 3*ETIJ(3,3,0,X,K,ipass)! - 0*ETIJ(2,3,*,X,K,ipass)/TWOB(K)
            ETIJ(0,3,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,1,X,K,ipass)+ ETIJ( 1,3,1,X,K,ipass)&
                 & -1*ETIJ(0,3,0,X,K,ipass)/TWOB(K)
            ETIJ(4,3,2,X,K,ipass) = HPINV(K)*ETIJ(3,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(4,3,1,X,K,ipass)
            ETIJ(5,3,2,X,K,ipass) = HPINV(K)*ETIJ(4,3,1,X,K,ipass)
            ETIJ(1,3,2,X,K,ipass)= HPINV(K)*ETIJ(0,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,1,X,K,ipass)&
                 & + 2*ETIJ(2,3,1,X,K,ipass) - 1*ETIJ(1,3,0,X,K,ipass)/TWOB(K)
            ETIJ(2,3,2,X,K,ipass)= HPINV(K)*ETIJ(1,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,1,X,K,ipass)&
                 & + 3*ETIJ(3,3,1,X,K,ipass) - 1*ETIJ(2,3,0,X,K,ipass)/TWOB(K)
            ETIJ(3,3,2,X,K,ipass)= HPINV(K)*ETIJ(2,3,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,1,X,K,ipass)&
                 & + 4*ETIJ(4,3,1,X,K,ipass) - 1*ETIJ(3,3,0,X,K,ipass)/TWOB(K)
            ETIJ(0,3,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,3,2,X,K,ipass)+ ETIJ( 1,3,2,X,K,ipass)&
                 & -2*ETIJ(0,3,1,X,K,ipass)/TWOB(K)
            ETIJ(5,3,3,X,K,ipass) = HPINV(K)*ETIJ(4,3,2,X,K,ipass)+ PB(K,ipass)*ETIJ(5,3,2,X,K,ipass)
            ETIJ(6,3,3,X,K,ipass) = HPINV(K)*ETIJ(5,3,2,X,K,ipass)
            ETIJ(1,3,3,X,K,ipass)= HPINV(K)*ETIJ(0,3,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,3,2,X,K,ipass)&
                 & + 2*ETIJ(2,3,2,X,K,ipass) - 2*ETIJ(1,3,1,X,K,ipass)/TWOB(K)
            ETIJ(2,3,3,X,K,ipass)= HPINV(K)*ETIJ(1,3,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,3,2,X,K,ipass)&
                 & + 3*ETIJ(3,3,2,X,K,ipass) - 2*ETIJ(2,3,1,X,K,ipass)/TWOB(K)
            ETIJ(3,3,3,X,K,ipass)= HPINV(K)*ETIJ(2,3,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,3,2,X,K,ipass)&
                 & + 4*ETIJ(4,3,2,X,K,ipass) - 2*ETIJ(3,3,1,X,K,ipass)/TWOB(K)
            ETIJ(4,3,3,X,K,ipass)= HPINV(K)*ETIJ(3,3,2,X,K,ipass)+ PB(K,ipass)*ETIJ(4,3,2,X,K,ipass)&
                 & + 5*ETIJ(5,3,2,X,K,ipass) - 2*ETIJ(4,3,1,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.2)THEN
      !maxI.EQ.2.AND.maxJ.EQ.4
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,0,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,3,X,K,ipass)+ ETIJ( 1,0,3,X,K,ipass)&
                 & -3*ETIJ(0,0,2,X,K,ipass)/TWOB(K)
            ETIJ(3,0,4,X,K,ipass) = HPINV(K)*ETIJ(2,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(4,0,4,X,K,ipass) = HPINV(K)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(1,0,4,X,K,ipass)= HPINV(K)*ETIJ(0,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,3,X,K,ipass)&
                 & + 2*ETIJ(2,0,3,X,K,ipass) - 3*ETIJ(1,0,2,X,K,ipass)/TWOB(K)
            ETIJ(2,0,4,X,K,ipass)= HPINV(K)*ETIJ(1,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,3,X,K,ipass)&
                 & + 3*ETIJ(3,0,3,X,K,ipass) - 3*ETIJ(2,0,2,X,K,ipass)/TWOB(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,1,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,2,X,K,ipass)+ ETIJ( 1,1,2,X,K,ipass)&
                 & -2*ETIJ(0,1,1,X,K,ipass)/TWOB(K)
            ETIJ(3,1,3,X,K,ipass) = HPINV(K)*ETIJ(2,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(4,1,3,X,K,ipass) = HPINV(K)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(1,1,3,X,K,ipass)= HPINV(K)*ETIJ(0,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,2,X,K,ipass)&
                 & + 2*ETIJ(2,1,2,X,K,ipass) - 2*ETIJ(1,1,1,X,K,ipass)/TWOB(K)
            ETIJ(2,1,3,X,K,ipass)= HPINV(K)*ETIJ(1,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,2,X,K,ipass)&
                 & + 3*ETIJ(3,1,2,X,K,ipass) - 2*ETIJ(2,1,1,X,K,ipass)/TWOB(K)
            ETIJ(0,1,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,3,X,K,ipass)+ ETIJ( 1,1,3,X,K,ipass)&
                 & -3*ETIJ(0,1,2,X,K,ipass)/TWOB(K)
            ETIJ(4,1,4,X,K,ipass) = HPINV(K)*ETIJ(3,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(4,1,3,X,K,ipass)
            ETIJ(5,1,4,X,K,ipass) = HPINV(K)*ETIJ(4,1,3,X,K,ipass)
            ETIJ(1,1,4,X,K,ipass)= HPINV(K)*ETIJ(0,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,3,X,K,ipass)&
                 & + 2*ETIJ(2,1,3,X,K,ipass) - 3*ETIJ(1,1,2,X,K,ipass)/TWOB(K)
            ETIJ(2,1,4,X,K,ipass)= HPINV(K)*ETIJ(1,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,3,X,K,ipass)&
                 & + 3*ETIJ(3,1,3,X,K,ipass) - 3*ETIJ(2,1,2,X,K,ipass)/TWOB(K)
            ETIJ(3,1,4,X,K,ipass)= HPINV(K)*ETIJ(2,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,3,X,K,ipass)&
                 & + 4*ETIJ(4,1,3,X,K,ipass) - 3*ETIJ(3,1,2,X,K,ipass)/TWOB(K)
            ETIJ(0,2,0,X,K,ipass) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(1,2,0,X,K,ipass) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(2,2,0,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,2,1,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,0,X,K,ipass)+ ETIJ( 1,2,0,X,K,ipass)
            ETIJ(2,2,1,X,K,ipass) = HPINV(K)*ETIJ(1,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(3,2,1,X,K,ipass) = HPINV(K)*ETIJ(2,2,0,X,K,ipass)
            ETIJ(1,2,1,X,K,ipass)= HPINV(K)*ETIJ(0,2,0,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,0,X,K,ipass)&
                 & + 2*ETIJ(2,2,0,X,K,ipass)! - 0*ETIJ(1,2,*,X,K,ipass)/TWOB(K)
            ETIJ(0,2,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,1,X,K,ipass)+ ETIJ( 1,2,1,X,K,ipass)&
                 & -1*ETIJ(0,2,0,X,K,ipass)/TWOB(K)
            ETIJ(3,2,2,X,K,ipass) = HPINV(K)*ETIJ(2,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(4,2,2,X,K,ipass) = HPINV(K)*ETIJ(3,2,1,X,K,ipass)
            ETIJ(1,2,2,X,K,ipass)= HPINV(K)*ETIJ(0,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,1,X,K,ipass)&
                 & + 2*ETIJ(2,2,1,X,K,ipass) - 1*ETIJ(1,2,0,X,K,ipass)/TWOB(K)
            ETIJ(2,2,2,X,K,ipass)= HPINV(K)*ETIJ(1,2,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,1,X,K,ipass)&
                 & + 3*ETIJ(3,2,1,X,K,ipass) - 1*ETIJ(2,2,0,X,K,ipass)/TWOB(K)
            ETIJ(0,2,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,2,X,K,ipass)+ ETIJ( 1,2,2,X,K,ipass)&
                 & -2*ETIJ(0,2,1,X,K,ipass)/TWOB(K)
            ETIJ(4,2,3,X,K,ipass) = HPINV(K)*ETIJ(3,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(4,2,2,X,K,ipass)
            ETIJ(5,2,3,X,K,ipass) = HPINV(K)*ETIJ(4,2,2,X,K,ipass)
            ETIJ(1,2,3,X,K,ipass)= HPINV(K)*ETIJ(0,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,2,X,K,ipass)&
                 & + 2*ETIJ(2,2,2,X,K,ipass) - 2*ETIJ(1,2,1,X,K,ipass)/TWOB(K)
            ETIJ(2,2,3,X,K,ipass)= HPINV(K)*ETIJ(1,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,2,X,K,ipass)&
                 & + 3*ETIJ(3,2,2,X,K,ipass) - 2*ETIJ(2,2,1,X,K,ipass)/TWOB(K)
            ETIJ(3,2,3,X,K,ipass)= HPINV(K)*ETIJ(2,2,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,2,X,K,ipass)&
                 & + 4*ETIJ(4,2,2,X,K,ipass) - 2*ETIJ(3,2,1,X,K,ipass)/TWOB(K)
            ETIJ(0,2,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,2,3,X,K,ipass)+ ETIJ( 1,2,3,X,K,ipass)&
                 & -3*ETIJ(0,2,2,X,K,ipass)/TWOB(K)
            ETIJ(5,2,4,X,K,ipass) = HPINV(K)*ETIJ(4,2,3,X,K,ipass)+ PB(K,ipass)*ETIJ(5,2,3,X,K,ipass)
            ETIJ(6,2,4,X,K,ipass) = HPINV(K)*ETIJ(5,2,3,X,K,ipass)
            ETIJ(1,2,4,X,K,ipass)= HPINV(K)*ETIJ(0,2,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,2,3,X,K,ipass)&
                 & + 2*ETIJ(2,2,3,X,K,ipass) - 3*ETIJ(1,2,2,X,K,ipass)/TWOB(K)
            ETIJ(2,2,4,X,K,ipass)= HPINV(K)*ETIJ(1,2,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,2,3,X,K,ipass)&
                 & + 3*ETIJ(3,2,3,X,K,ipass) - 3*ETIJ(2,2,2,X,K,ipass)/TWOB(K)
            ETIJ(3,2,4,X,K,ipass)= HPINV(K)*ETIJ(2,2,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,2,3,X,K,ipass)&
                 & + 4*ETIJ(4,2,3,X,K,ipass) - 3*ETIJ(3,2,2,X,K,ipass)/TWOB(K)
            ETIJ(4,2,4,X,K,ipass)= HPINV(K)*ETIJ(3,2,3,X,K,ipass)+ PB(K,ipass)*ETIJ(4,2,3,X,K,ipass)&
                 & + 5*ETIJ(5,2,3,X,K,ipass) - 3*ETIJ(4,2,2,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.5
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,0,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,3,X,K,ipass)+ ETIJ( 1,0,3,X,K,ipass)&
                 & -3*ETIJ(0,0,2,X,K,ipass)/TWOB(K)
            ETIJ(3,0,4,X,K,ipass) = HPINV(K)*ETIJ(2,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(4,0,4,X,K,ipass) = HPINV(K)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(1,0,4,X,K,ipass)= HPINV(K)*ETIJ(0,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,3,X,K,ipass)&
                 & + 2*ETIJ(2,0,3,X,K,ipass) - 3*ETIJ(1,0,2,X,K,ipass)/TWOB(K)
            ETIJ(2,0,4,X,K,ipass)= HPINV(K)*ETIJ(1,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,3,X,K,ipass)&
                 & + 3*ETIJ(3,0,3,X,K,ipass) - 3*ETIJ(2,0,2,X,K,ipass)/TWOB(K)
            ETIJ(0,0,5,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,4,X,K,ipass)+ ETIJ( 1,0,4,X,K,ipass)&
                 & -4*ETIJ(0,0,3,X,K,ipass)/TWOB(K)
            ETIJ(4,0,5,X,K,ipass) = HPINV(K)*ETIJ(3,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(4,0,4,X,K,ipass)
            ETIJ(5,0,5,X,K,ipass) = HPINV(K)*ETIJ(4,0,4,X,K,ipass)
            ETIJ(1,0,5,X,K,ipass)= HPINV(K)*ETIJ(0,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,4,X,K,ipass)&
                 & + 2*ETIJ(2,0,4,X,K,ipass) - 4*ETIJ(1,0,3,X,K,ipass)/TWOB(K)
            ETIJ(2,0,5,X,K,ipass)= HPINV(K)*ETIJ(1,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,4,X,K,ipass)&
                 & + 3*ETIJ(3,0,4,X,K,ipass) - 4*ETIJ(2,0,3,X,K,ipass)/TWOB(K)
            ETIJ(3,0,5,X,K,ipass)= HPINV(K)*ETIJ(2,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,4,X,K,ipass)&
                 & + 4*ETIJ(4,0,4,X,K,ipass) - 4*ETIJ(3,0,3,X,K,ipass)/TWOB(K)
            ETIJ(0,1,0,X,K,ipass) = PA(K,ipass)
            ETIJ(1,1,0,X,K,ipass) = HPINV(K)
            ETIJ(0,1,1,X,K,ipass) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(1,1,1,X,K,ipass) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(2,1,1,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,1,2,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,1,X,K,ipass)+ ETIJ( 1,1,1,X,K,ipass)&
                 & -1*ETIJ(0,1,0,X,K,ipass)/TWOB(K)
            ETIJ(2,1,2,X,K,ipass) = HPINV(K)*ETIJ(1,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(3,1,2,X,K,ipass) = HPINV(K)*ETIJ(2,1,1,X,K,ipass)
            ETIJ(1,1,2,X,K,ipass)= HPINV(K)*ETIJ(0,1,1,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,1,X,K,ipass)&
                 & + 2*ETIJ(2,1,1,X,K,ipass) - 1*ETIJ(1,1,0,X,K,ipass)/TWOB(K)
            ETIJ(0,1,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,2,X,K,ipass)+ ETIJ( 1,1,2,X,K,ipass)&
                 & -2*ETIJ(0,1,1,X,K,ipass)/TWOB(K)
            ETIJ(3,1,3,X,K,ipass) = HPINV(K)*ETIJ(2,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(4,1,3,X,K,ipass) = HPINV(K)*ETIJ(3,1,2,X,K,ipass)
            ETIJ(1,1,3,X,K,ipass)= HPINV(K)*ETIJ(0,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,2,X,K,ipass)&
                 & + 2*ETIJ(2,1,2,X,K,ipass) - 2*ETIJ(1,1,1,X,K,ipass)/TWOB(K)
            ETIJ(2,1,3,X,K,ipass)= HPINV(K)*ETIJ(1,1,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,2,X,K,ipass)&
                 & + 3*ETIJ(3,1,2,X,K,ipass) - 2*ETIJ(2,1,1,X,K,ipass)/TWOB(K)
            ETIJ(0,1,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,3,X,K,ipass)+ ETIJ( 1,1,3,X,K,ipass)&
                 & -3*ETIJ(0,1,2,X,K,ipass)/TWOB(K)
            ETIJ(4,1,4,X,K,ipass) = HPINV(K)*ETIJ(3,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(4,1,3,X,K,ipass)
            ETIJ(5,1,4,X,K,ipass) = HPINV(K)*ETIJ(4,1,3,X,K,ipass)
            ETIJ(1,1,4,X,K,ipass)= HPINV(K)*ETIJ(0,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,3,X,K,ipass)&
                 & + 2*ETIJ(2,1,3,X,K,ipass) - 3*ETIJ(1,1,2,X,K,ipass)/TWOB(K)
            ETIJ(2,1,4,X,K,ipass)= HPINV(K)*ETIJ(1,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,3,X,K,ipass)&
                 & + 3*ETIJ(3,1,3,X,K,ipass) - 3*ETIJ(2,1,2,X,K,ipass)/TWOB(K)
            ETIJ(3,1,4,X,K,ipass)= HPINV(K)*ETIJ(2,1,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,3,X,K,ipass)&
                 & + 4*ETIJ(4,1,3,X,K,ipass) - 3*ETIJ(3,1,2,X,K,ipass)/TWOB(K)
            ETIJ(0,1,5,X,K,ipass) = PB(K,ipass)*ETIJ( 0,1,4,X,K,ipass)+ ETIJ( 1,1,4,X,K,ipass)&
                 & -4*ETIJ(0,1,3,X,K,ipass)/TWOB(K)
            ETIJ(5,1,5,X,K,ipass) = HPINV(K)*ETIJ(4,1,4,X,K,ipass)+ PB(K,ipass)*ETIJ(5,1,4,X,K,ipass)
            ETIJ(6,1,5,X,K,ipass) = HPINV(K)*ETIJ(5,1,4,X,K,ipass)
            ETIJ(1,1,5,X,K,ipass)= HPINV(K)*ETIJ(0,1,4,X,K,ipass)+ PB(K,ipass)*ETIJ(1,1,4,X,K,ipass)&
                 & + 2*ETIJ(2,1,4,X,K,ipass) - 4*ETIJ(1,1,3,X,K,ipass)/TWOB(K)
            ETIJ(2,1,5,X,K,ipass)= HPINV(K)*ETIJ(1,1,4,X,K,ipass)+ PB(K,ipass)*ETIJ(2,1,4,X,K,ipass)&
                 & + 3*ETIJ(3,1,4,X,K,ipass) - 4*ETIJ(2,1,3,X,K,ipass)/TWOB(K)
            ETIJ(3,1,5,X,K,ipass)= HPINV(K)*ETIJ(2,1,4,X,K,ipass)+ PB(K,ipass)*ETIJ(3,1,4,X,K,ipass)&
                 & + 4*ETIJ(4,1,4,X,K,ipass) - 4*ETIJ(3,1,3,X,K,ipass)/TWOB(K)
            ETIJ(4,1,5,X,K,ipass)= HPINV(K)*ETIJ(3,1,4,X,K,ipass)+ PB(K,ipass)*ETIJ(4,1,4,X,K,ipass)&
                 & + 5*ETIJ(5,1,4,X,K,ipass) - 4*ETIJ(4,1,3,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.6
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(0,0,0,X,K,ipass) = D1
            ETIJ(0,0,1,X,K,ipass) = PB(K,ipass)    
            ETIJ(1,0,1,X,K,ipass) = HPINV(K) 
            ETIJ(0,0,2,X,K,ipass) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(1,0,2,X,K,ipass) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(2,0,2,X,K,ipass) = HPINV(K)*HPINV(K)
            ETIJ(0,0,3,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,2,X,K,ipass)+ ETIJ( 1,0,2,X,K,ipass)&
                 & -2*ETIJ(0,0,1,X,K,ipass)/TWOB(K)
            ETIJ(2,0,3,X,K,ipass) = HPINV(K)*ETIJ(1,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(3,0,3,X,K,ipass) = HPINV(K)*ETIJ(2,0,2,X,K,ipass)
            ETIJ(1,0,3,X,K,ipass)= HPINV(K)*ETIJ(0,0,2,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,2,X,K,ipass)&
                 & + 2*ETIJ(2,0,2,X,K,ipass) - 2*ETIJ(1,0,1,X,K,ipass)/TWOB(K)
            ETIJ(0,0,4,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,3,X,K,ipass)+ ETIJ( 1,0,3,X,K,ipass)&
                 & -3*ETIJ(0,0,2,X,K,ipass)/TWOB(K)
            ETIJ(3,0,4,X,K,ipass) = HPINV(K)*ETIJ(2,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(4,0,4,X,K,ipass) = HPINV(K)*ETIJ(3,0,3,X,K,ipass)
            ETIJ(1,0,4,X,K,ipass)= HPINV(K)*ETIJ(0,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,3,X,K,ipass)&
                 & + 2*ETIJ(2,0,3,X,K,ipass) - 3*ETIJ(1,0,2,X,K,ipass)/TWOB(K)
            ETIJ(2,0,4,X,K,ipass)= HPINV(K)*ETIJ(1,0,3,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,3,X,K,ipass)&
                 & + 3*ETIJ(3,0,3,X,K,ipass) - 3*ETIJ(2,0,2,X,K,ipass)/TWOB(K)
            ETIJ(0,0,5,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,4,X,K,ipass)+ ETIJ( 1,0,4,X,K,ipass)&
                 & -4*ETIJ(0,0,3,X,K,ipass)/TWOB(K)
            ETIJ(4,0,5,X,K,ipass) = HPINV(K)*ETIJ(3,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(4,0,4,X,K,ipass)
            ETIJ(5,0,5,X,K,ipass) = HPINV(K)*ETIJ(4,0,4,X,K,ipass)
            ETIJ(1,0,5,X,K,ipass)= HPINV(K)*ETIJ(0,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,4,X,K,ipass)&
                 & + 2*ETIJ(2,0,4,X,K,ipass) - 4*ETIJ(1,0,3,X,K,ipass)/TWOB(K)
            ETIJ(2,0,5,X,K,ipass)= HPINV(K)*ETIJ(1,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,4,X,K,ipass)&
                 & + 3*ETIJ(3,0,4,X,K,ipass) - 4*ETIJ(2,0,3,X,K,ipass)/TWOB(K)
            ETIJ(3,0,5,X,K,ipass)= HPINV(K)*ETIJ(2,0,4,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,4,X,K,ipass)&
                 & + 4*ETIJ(4,0,4,X,K,ipass) - 4*ETIJ(3,0,3,X,K,ipass)/TWOB(K)
            ETIJ(0,0,6,X,K,ipass) = PB(K,ipass)*ETIJ( 0,0,5,X,K,ipass)+ ETIJ( 1,0,5,X,K,ipass)&
                 & -5*ETIJ(0,0,4,X,K,ipass)/TWOB(K)
            ETIJ(5,0,6,X,K,ipass) = HPINV(K)*ETIJ(4,0,5,X,K,ipass)+ PB(K,ipass)*ETIJ(5,0,5,X,K,ipass)
            ETIJ(6,0,6,X,K,ipass) = HPINV(K)*ETIJ(5,0,5,X,K,ipass)
            ETIJ(1,0,6,X,K,ipass)= HPINV(K)*ETIJ(0,0,5,X,K,ipass)+ PB(K,ipass)*ETIJ(1,0,5,X,K,ipass)&
                 & + 2*ETIJ(2,0,5,X,K,ipass) - 5*ETIJ(1,0,4,X,K,ipass)/TWOB(K)
            ETIJ(2,0,6,X,K,ipass)= HPINV(K)*ETIJ(1,0,5,X,K,ipass)+ PB(K,ipass)*ETIJ(2,0,5,X,K,ipass)&
                 & + 3*ETIJ(3,0,5,X,K,ipass) - 5*ETIJ(2,0,4,X,K,ipass)/TWOB(K)
            ETIJ(3,0,6,X,K,ipass)= HPINV(K)*ETIJ(2,0,5,X,K,ipass)+ PB(K,ipass)*ETIJ(3,0,5,X,K,ipass)&
                 & + 4*ETIJ(4,0,5,X,K,ipass) - 5*ETIJ(3,0,4,X,K,ipass)/TWOB(K)
            ETIJ(4,0,6,X,K,ipass)= HPINV(K)*ETIJ(3,0,5,X,K,ipass)+ PB(K,ipass)*ETIJ(4,0,5,X,K,ipass)&
                 & + 5*ETIJ(5,0,5,X,K,ipass) - 5*ETIJ(4,0,4,X,K,ipass)/TWOB(K)
         ENDDO
      ENDDO
   ENDIF
ELSE
   CALL ichorquit('EXTEND ICHOR_HERM_ECOEFFS (MAIN2 IN etij)',-1)
ENDIF
END SUBROUTINE ICHOR_HERM_ECOEFFS

subroutine ICHOR_Ecoeffn_maxAngP0_maxAngA0(nPrimPnPass,Ecoeffn,ETIJ,pref)
  implicit none
  integer,intent(in) :: nPrimPnPass
  real(realk),intent(in) :: ETIJ(0:0,0:0,0:0,3,nprimPnPass),pref(nPrimPnPass)
  real(realk),intent(inout) :: Ecoeffn(1,  1,nPrimPnPass)
!  
  integer :: i
  DO i = 1, nPrimPnPass
   EcoeffN(1,1,i) = ETIJ(0,0,0,1,i)*ETIJ(0,0,0,2,i)*ETIJ(0,0,0,3,i)*pref(i)
  ENDDO
End subroutine ICHOR_Ecoeffn_maxAngP0_maxAngA0
  
subroutine ICHOR_Ecoeffn_general(nPrimPnPass,nTUV,ijk,JMAX,l1,l2,EcoeffN,ETIJ,pref,RHS)!,Enoscreen)
  implicit none
  integer,intent(in) :: nPrimPnPass,nTUV,ijk,JMAX
  real(realk),intent(in) :: ETIJ(0:JMAX,0:l1,0:l2,3,nPrimPnPass)
  real(realk),intent(in) :: pref(nPrimPnPass)
  real(realk),intent(inout) :: EcoeffN(nTUV,ijk,nPrimPnPass)
  logical :: RHS
!  logical,intent(inout) :: Enoscreen(nTUV,ijk)
  !
  integer :: TUVINDEX(0:JMAX,0:JMAX,0:JMAX)
  integer :: l2,l1,P2,iP2,JP2,kp2,P1,jp1,ip1,kp1,tp,up,vp,ituvp,J,ijkQ,i
  real(realk),parameter :: D1=1.0E0_realk,signQ=-1.0E0_realk
  real(realk) :: sign,signijk
!  Enoscreen = .FALSE.
  sign = D1
  IF(RHS) sign = signQ

  ituvP = 0
  DO J = 0, JMAX
     DO Tp=J,0,-1
        DO Up=J-Tp,0,-1
           Vp=J-Tp-Up
           ituvP = ituvP + 1 
           TUVINDEX(Tp,Up,Vp) = ituvP
        ENDDO
     ENDDO
  ENDDO

  ijkQ=0
  DO P2 = 0,l2
   DO iP2=P2,0,-1
    DO jP2=P2-iP2,0,-1
     kP2=P2-iP2-jP2
     DO P1 = 0,l1
      DO iP1=P1,0,-1
       DO jP1=P1-iP1,0,-1
        kP1=P1-iP1-jP1
        IF(iP1+jP1+kP1+iP2+jP2+kP2 .EQ. l1+l2) then
         ijkQ=ijkQ+1
         DO tP=0,iP1+iP2
          DO uP=0,jP1+jP2
           DO vP=0,kP1+kP2
            IF(MOD(tP+uP+vP,2).EQ. 0)THEN
               signijk = D1
            ELSE
               signijk = sign
            ENDIF
            ituvP = TUVINDEX(tp,up,vp)
            DO i = 1, nPrimPnPass
             EcoeffN(ituvP,ijkQ,i) = signijk*ETIJ(tp,iP1,iP2,1,i)*ETIJ(up,jP1,jP2,2,i)*ETIJ(vp,kP1,kP2,3,i)*pref(i)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine ICHOR_Ecoeffn_general

!!$subroutine DirectContractEcoeffN_maxAngP0_maxAngQ0_maxAngC0(nPrimP,nPrimQ,nPasses,Ecoeffn,WTUV,OUT)
!!$ implicit none
!!$ integer,intent(in) :: nPrimQ,nPasses,nPrimP
!!$ real(realk),intent(in) :: Ecoeffn(1,1,nPrimQ,nPasses)
!!$ real(realk),intent(in) :: WTUV(1,nPrimP,nPrimQ,nPasses)
!!$ real(realk),intent(inout) :: OUT(1,nPrimP,nPrimQ,nPasses)
!!$ !
!!$ integer :: iPasses,iPrimQ,iPrimP
!!$ DO iPasses=1,nPasses
!!$  DO iPrimQ=1,nPrimQ
!!$   DO iPrimP=1,nPrimP
!!$    OUT(1,iPrimP,iPrimQ,iPasses) = Ecoeffn(1,1,iPrimQ,iPasses)*WTUV(1,iPrimP,iPrimQ,iPasses)
!!$   ENDDO
!!$  ENDDO
!!$ ENDDO
!!$END SUBROUTINE DirectContractEcoeffN_maxAngP0_maxAngQ0_maxAngC0
!!$
!!$SUBROUTINE contractEcoeffSeg(IntegralIN,IntegralOUT,Etuv,&
!!$     & ijkQsph,ijkPcart,nTUVP,nPrimP,nContQ,nPasses)
!!$implicit none
!!$Integer,intent(IN) :: nOrb,nTUV,nAng,nPrim
!!$Real(realk),intent(IN)  :: IN(ijkQsph,nTUVP,nPrimP,nContQ,nPasses)
!!$Real(realk),intent(IN)  :: Etuv(nTUVP,ijkPcart,nPrim)
!!$Real(realk),intent(OUT) :: OUT(ijkQsph,ijkPcart,nPrimP,nContQ,nPasses)
!!$!
!!$Integer      :: iPrimP,iTUV,iAng,iOrb
!!$Real(realk) :: tmp
!!$!
!!$!iContQ,i
!!$DO iPasses = 1,nPasses
!!$ DO iContQ = 1,nContQ
!!$  DO iijkPcart = 1,ijkPcart
!!$   DO iPrimP = 1,nPrim
!!$    iTUVP=1
!!$     DO iijkQsph = 1,ijkQsph
!!$      OUT(iijkQsph,iijkPcart,iPrimP,iContQ,iPasses) = &
!!$          & IN(iijkQsph,iTUVP,iPrimP,iContQ,iPasses)*Etuv(iTUVP,iijkPcart,iPrimP)
!!$     ENDDO
!!$    DO iTUVP=2,nTUVP
!!$     DO iijkQsph = 1,ijkQsph
!!$      OUT(iijkQsph,iijkPcart,iPrimP,iContQ,iPasses) = OUT(iijkQsph,iijkPcart,iPrimP,iContQ,iPasses) +&
!!$          & IN(iijkQsph,iTUVP,iPrimP,iContQ,iPasses)*Etuv(iTUVP,iijkPcart,iPrimP)
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO  
!!$  ENDDO
!!$ ENDDO
!!$ENDDO
!!$!
!!$END SUBROUTINE contractEcoeffSeg
!!$
!!$SUBROUTINE contractEcoeffGen(IntegralIN,IntegralOUT,Etuv,&
!!$     & ijkQsph,ijkPcart,nTUVP,nPrimP,nContQ,nPasses)
!!$implicit none
!!$Integer,intent(IN) :: nOrb,nTUV,nAng,nPrim
!!$Real(realk),intent(IN)  :: IN(ijkQsph,nTUVP,nPrimP,nContQ,nPasses)
!!$Real(realk),intent(IN)  :: Etuv(nTUVP,ijkPcart,nPrim)
!!$Real(realk),intent(OUT) :: OUT(ijkQsph,ijkPcart,nContQ,nPasses)
!!$!
!!$Integer      :: iPrimP,iTUV,iAng,iOrb
!!$Real(realk) :: tmp
!!$!
!!$iContQ,i
!!$DO iPasses = 1,nPasses
!!$ DO iContQ = 1,nContQ
!!$  DO iijkPcart = 1,ijkPcart
!!$   DO iijkQsph = 1,ijkQsph
!!$    tmp = 0.0E0_realk
!!$    DO iPrimP = 1,nPrim
!!$     DO iTUVP=1,nTUVP
!!$        tmp = tmp + IN(iijkQsph,iTUVP,iPrimP,iContQ,iPasses)*Etuv(iTUVP,iijkPcart,iPrimP)
!!$     ENDDO
!!$    ENDDO
!!$    OUT(iijkQsph,iijkPcart,iContQ,iPasses) = tmp
!!$   ENDDO  
!!$  ENDDO
!!$ ENDDO
!!$ENDDO
!!$!
!!$END SUBROUTINE contractEcoeffGen
subroutine PrintIchorTensorRE(RE,nTUVP,ijkQcart,nPrimP,nPasses,lupri)
implicit none
integer,intent(in) :: nTUVP,ijkQcart,nPrimP,nPasses,lupri
real(realk),intent(in) :: RE(nTUVP,ijkQcart,nPrimP,nPasses)
integer :: iTUV,iPass,iP,ijkQ
iTUV=0
WRITE(lupri,'(A)')'Print the RTUV tensor contracted with Ecoeff segmented: RE(nTUVP,ijkQcart,nPrimP,nPasses)'
WRITE(lupri,'(A,I7)')'nTUVP:   ',nTUVP
WRITE(lupri,'(A,I7)')'ijkQcart:',ijkQcart
WRITE(lupri,'(A,I7)')'nPrimP:  ',nPrimP
WRITE(lupri,'(A,I7)')'nPasses: ',nPasses
do iPass = 1,nPasses
 do iP = 1,nPrimP
  do ijkQ = 1,ijkQcart
   WRITE(lupri,'(A4,1X,A6,1X,A5,1X,A12)')'ijkQ','iPrimP','iPass','iTUV=1,nTUVP'
   WRITE(LUPRI,'(I4,I7,I6,5ES18.9/,(17X,5ES18.9))') ijkQ,iP,iPass,(RE(iTUV,ijkQ,iP,iPass),iTUV=1,nTUVP)
  enddo
 enddo
enddo
end subroutine PrintIchorTensorRE

subroutine PrintIchorTensorERE(ERE,ijkPcart,ijkQcart,nPasses,lupri)
implicit none
integer,intent(in) :: ijkQcart,ijkPcart,nPasses,lupri
real(realk),intent(in) :: ERE(ijkPcart,ijkQcart,nPasses)
!
integer :: iPass,iP,ijkQ,ijkP
WRITE(lupri,'(A)')'Print the RTUV tensor contracted with Ecoeff segmented: RE(nTUVP,ijkQcart,nPrimP,nPasses)'
WRITE(lupri,'(A,I7)')'ijkPcart:',ijkPcart
WRITE(lupri,'(A,I7)')'ijkQcart:',ijkQcart
WRITE(lupri,'(A,I7)')'nPasses: ',nPasses
do iPass = 1,nPasses
 do ijkQ = 1,ijkQcart
  WRITE(lupri,'(A4,1X,A6,1X,A15)')'ijkQ','iPrimP','ijkP=1,ijkPcart'
  WRITE(LUPRI,'(I4,I7,5ES18.9/,(11X,5ES18.9))') ijkQ,iPass,(ERE(ijkP,ijkQ,iPass),ijkP=1,ijkPcart)
 enddo
enddo
end subroutine PrintIchorTensorERE

subroutine printEcoeff(EcoeffN,nTUVQ,ijkQcart,nPrimP,nPass,lupri)
implicit none
integer     :: nTUVQ,ijkQcart,nPrimP,nPass,lupri
real(realk) :: EcoeffN(nTUVQ,ijkQcart,nPrimP,nPass)
!
integer :: i,j,k,ip
k=0
WRITE(LUPRI,'(A)') 'EcoeffN'
WRITE(lupri,'(A,I7)')'nTUVQ:    ',nTUVQ
WRITE(lupri,'(A,I7)')'ijkQcart: ',ijkQcart
WRITE(lupri,'(A,I7)')'nPrimP:   ',nPrimP
WRITE(lupri,'(A,I7)')'nPass:    ',nPass
WRITE(LUPRI,'(1X,A6,1X,A6,1X,A5,12X,A14 )') 'prim  ','tuv   ','iPass ','ijk=1,ijkQcart'
DO ip=1,nPass
 DO i=1,nPrimP
   DO j=1,nTUVQ
      WRITE(LUPRI,'(1X,I4,2X,I4,1X,I3,5ES18.9/,(15X,5ES18.9))') i,j,ip,(EcoeffN(j,k,i,ip),k=1,ijkQcart)
   ENDDO
 ENDDO
ENDDO

end subroutine printEcoeff

end MODULE IchorEriCoulombintegralGeneralMod
