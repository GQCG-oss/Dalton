!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriCoulombintegralCPUMcMGeneralMod
use IchorPrecisionModule
use IchorMemory
use IchorCommonModule
use IchorEriCoulombintegralCPUMcMGeneralEcoeffMod, only: &
     & Ichorbuild_Ecoeff_RHS,Ichorbuild_Ecoeff_LHS, printEcoeff
use IchorEriCoulombintegralCPUMcMGeneralWTUVMod

public 
!build from old IchorEri_CoulombIntegral_general.f90 in
!/home/tkjaer/DaltonDevelopment/ExplicitIntegrals/LSint

! different paths depending on nCont, ...
! Segmentation !!!
! Spherical Ecoefficients

type vector
real(realk),pointer :: elms(:)
end type vector
type(vector),allocatable :: SPH_MAT(:)
integer :: nTmpArray3,nTmpArray4
real(realk),allocatable :: TmpArray3(:),TmpArray4(:)

CONTAINS
  subroutine DetermineSizeTmpArray34(nTUVQ,nCartOrbCompQ,nPrimQ,nTUVP,nCartOrbCompP,nPrimP,nPasses,&
       & AngmomA,AngmomB,AngmomC,AngmomD,AngmomP,AngmomQ,AngmomPQ)
    implicit none
    integer,intent(in) :: nTUVQ,nCartOrbCompQ,nPrimQ,nTUVP,nCartOrbCompP,nPrimP,nPasses
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD,AngmomP,AngmomQ,AngmomPQ
    !TmpArray3 used Ecoeff and Rpq
    nTmpArray3 = MAX(nTUVQ*nCartOrbCompQ*nPrimQ,nTUVP*nCartOrbCompP*nPrimP*nPasses,nPrimQ*nPrimP*nPasses*3)
    !TmpArray4 used for RJ000 and ETIJ (tmp array in building Ecoeff)
    nTmpArray4 = MAX(nPasses*nPrimQ*nPrimP*(AngmomPQ+1),nPrimQ*(AngmomQ+1)*(AngmomC+1)*(AngmomD+1)*3)
    nTmpArray4 = MAX(nTmpArray4,nPasses*nPrimP*(AngmomP+1)*(AngmomA+1)*(AngmomB+1)*3)
  end subroutine DetermineSizeTmpArray34

  subroutine ICI_CPU_McM_general(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,&
       & CDAB,nLocalIntPass,Acenter,Bcenter,Ccenter,Dcenter,&
       & nAtomsA,nAtomsB,SphericalGTO,TmpArray1,TMParray1maxsize,TmpArray2,&
       & TMParray2maxsize,IatomAPass,iatomBPass)
    implicit none
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri,nLocalIntPass
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD,nAtomsA,nAtomsB
    integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
    integer,intent(in) :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD
    integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
    integer,intent(in) :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented,sphericalGTO,PQorder
    real(realk),intent(in) :: Acenter(3,nAtomsA),Bcenter(3,nAtomsB),Ccenter(3),Dcenter(3)
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP*nAtomsA*nAtomsB) !pcent(3,nPrimP,nPasses)
    real(realk),intent(in) :: qcent(3*nPrimQ)                 !qcent(3,nPrimQ)
    real(realk),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP*nAtomsA*nAtomsB)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(nLocalIntPass)    
    real(realk),intent(in) :: integralPrefactor(nPrimQP)!integralPrefactor(nPrimQ,nPrimP)    
    real(realk),intent(in) :: reducedExponents(nPrimQP) !reducedExponents(nPrimQ,nPrimP)
    real(realk),intent(in) :: Qdistance12(3)
    real(realk),intent(in) :: Pdistance12(3*nAtomsA*nAtomsB) !full set nAtomsA*nAtomsB
    !local variables
!    real(realk),allocatable :: squaredDistance(:),Rpq(:)
    integer :: AngmomP,AngmomQ,AngmomPQ
    logical :: RHS,TMP1
!    integer :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
!    integer :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD
!    integer :: nCartOrbCompP,nOrbCompP,nCartOrbCompQ,nOrbCompQ
    logical :: Sph1,Sph2,Sph3,Sph4,SphericalTransP,SphericalTransQ
    integer :: AngmomID,TUV,J,T,U,V,I
    integer,parameter :: nPassQ = 1
    !FIXME this should be called at a much higher place

    AngmomP = AngmomA+AngmomB     
    AngmomQ = AngmomC+AngmomD
    AngmomPQ  = AngmomP + AngmomQ
    Sph1 = sphericalGTO.AND.(AngmomA.GT. 1)
    Sph2 = sphericalGTO.AND.(AngmomB.GT. 1)
    Sph3 = sphericalGTO.AND.(AngmomC.GT. 1)
    Sph4 = sphericalGTO.AND.(AngmomD.GT. 1)
    SphericalTransP = Sph1.OR.Sph2
    SphericalTransQ = Sph3.OR.Sph4

#ifdef VAR_DEBUGICHOR
    IF(nTmpArray3.LT.nPrimQ*nPrimP*nPasses*3)call ichorquit('IchorTmp0G1',-1)
#endif
    !Rpq(nPrimQ,nPrimP,nPasses,3)
    call build_Rpq(nPrimQ,nPasses,nPrimP,Qcent,Pcent,TmpArray3,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB)

    IF(PQorder) call lsquit('PQorder CPU_McM general expect to get QP ordering',-1)

    !
    !      builds RJ000(0:AngmomPQ,nPrimQ,nPrimP,nPasses) Store in TmpArray4
    !
#ifdef VAR_DEBUGICHOR
    IF(nTmpArray4.LT.nPasses*nPrimQ*nPrimP*(AngmomPQ+1))call ichorquit('IchorTmp1G1',-1)
#endif
    call buildRJ000_general(nPasses,nPrimQ,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&
         & TABFJW,TmpArray4,AngmomPQ,integralPrefactor,TmpArray3)
    IF (INTPRINT .GE. 10) CALL PrintRJ000(TmpArray4,AngmomPQ,nPrimQ*nPrimP,nPasses,lupri)    
    !
    !     Build WTUV(nPrimQ,nPrimP,nPasses,nTUV) RJ000 = TmpArray4, Rpq = TmpArray3
    !
    TMP1 = .TRUE.
    IF (AngmomPQ.EQ. 0) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimQP*nPasses*ntuv)call ichorquit('IchorTmp1G3A',-1)
#endif
       call IchorwtuvRecurrenceJMIN0JMAX0(TmpArray4,TmpArray1,nPrimQP*nPasses)
    ELSEIF (AngmomPQ.EQ. 1) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimQP*nPasses*ntuv)call ichorquit('IchorTmp1G3B',-1)
#endif
       call IchorwtuvRecurrenceJMIN0JMAX1(TmpArray4,TmpArray1,TmpArray3,nPrimQP*nPasses)
    ELSEIF (AngmomPQ.EQ. 2) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimQP*nPasses*ntuv)call ichorquit('IchorTmp1G3C',-1)
#endif
       call IchorwtuvRecurrenceJMIN0JMAX2(TmpArray4,TmpArray1,TmpArray3,nPrimQP*nPasses)
    ELSEIF (AngmomPQ.EQ. 3) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimQP*nPasses*ntuv)call ichorquit('IchorTmp1G3D',-1)
#endif
       call IchorwtuvRecurrenceJMIN0JMAX3(TmpArray4,TmpArray1,TmpArray3,nPrimQP*nPasses)
    ELSE !AngmomPQ > 3
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimQP*nPasses*ntuv)call ichorquit('IchorTmp1G3a',-1)
       IF(TMParray2maxsize.LT.nPrimQP*nPasses*ntuv)call ichorquit('IchorTmp1G3b',-1)
#endif
       J = AngmomPQ-3
       call IchorwtuvRecurrenceJMIN0JMAX3J(TmpArray4,TmpArray1,TmpArray3,nPrimQP*nPasses,AngmomPQ,J)
       TMP1 = .TRUE.
       !J minimum goes from 0 to 0. For  AngmomPQ=6 J takes the values 2,1,0 
       DO j=AngmomPQ-4,0,-1 
          IF(TMP1)THEN
             call IchorwtuvRecurrenceCurrent(TmpArray1,TmpArray2,J,AngmomPQ,nPrimQ,nPrimP,nPasses,&
                  & ntuv,TmpArray3,TmpArray4)
          ELSE
             call IchorwtuvRecurrenceCurrent(TmpArray2,TmpArray1,J,AngmomPQ,nPrimQ,nPrimP,nPasses,&
                  & ntuv,TmpArray3,TmpArray4)
          ENDIF
          TMP1 = .NOT.TMP1
       ENDDO
    ENDIF    
    IF (IntPrint .GE. 25)THEN
       IF(TMP1)THEN
          call PrintWTUV(TmpArray1,AngmomPQ,nPrimQP,nPasses,nTUV,lupri)
       ELSE
          call PrintWTUV(TmpArray2,AngmomPQ,nPrimQP,nPasses,nTUV,lupri)
       ENDIF
    ENDIF

    !builds Ecoeff(nPrimQ,nPasses,nTUVQ,nCartOrbCompQ)
    !FOR NOW THIS IS NOT OpenMP parallized - unclear what the best method is. 
    call Ichorbuild_Ecoeff_RHS(nPrimQ,nPrimC,nPrimD,AngmomQ,AngmomC,AngmomD,nTUVQ,&
         & nCartOrbCompQ,Cexp,Dexp,TmpArray3,Qdistance12,Qpreexpfac,intprint,lupri,TmpArray4)

    IF (IntPrint .GE. 25)call printEcoeff(TmpArray3,nTUVQ,nCartOrbCompQ,nPrimQ,nPassQ,lupri)

    IF(TMP1)THEN !current intermediate WTUV reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nPrimQP*nPasses*ntuvP*nCartOrbCompQ)call ichorquit('IchorTmp1G3Q1',-1)
#endif
       !builds RE(nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ)
       call DirectcontractEQgen(TmpArray1,TmpArray2,nPrimQP,nPrimP*nPasses,nPrimQ,nTUV,ntuvP,ntuvQ,&
            & TmpArray3,nCartOrbCompQ,AngmomA,AngmomB,AngmomC,AngmomD,AngmomPQ)
       IF (IntPrint .GE. 25)call PrintIchorTensorRE(TmpArray2,nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ,lupri)
    ELSE !current intermediate WTUV reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimQP*nPasses*ntuvP*nCartOrbCompQ)call ichorquit('IchorTmp1G3Q2',-1)
#endif
       !builds RE(nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ)
       call DirectcontractEQgen(TmpArray2,TmpArray1,nPrimQP,nPrimP*nPasses,nPrimQ,nTUV,ntuvP,ntuvQ,&
            & TmpArray3,nCartOrbCompQ,AngmomA,AngmomB,AngmomC,AngmomD,AngmomPQ)
       IF (IntPrint .GE. 25)call PrintIchorTensorRE(TmpArray1,nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ,lupri)
    ENDIF
    TMP1 = .NOT.TMP1
    
    IF(TMP1)THEN !current intermediate RE(nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ) reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContQ*nPrimP*nPasses*nTUVP*nCartOrbCompQ)call ichorquit('IchorTmp2G4',-1)
#endif
       !Build REC(nContQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ)
       IF(Qsegmented)THEN
          call contractBasisSegQ(TmpArray1,TmpArray2,nPrimC,nPrimD,nPrimP*nPasses*nTUVP*nCartOrbCompQ) 
       ELSE
          call contractBasisGenQ(TmpArray1,TmpArray2,CCC,DCC,nPrimC,nPrimD,nContC,nContD,nPrimP*nPasses*nTUVP*nCartOrbCompQ) 
       ENDIF
       IF (IntPrint .GE. 25)call PrintIchorTensorREC(TmpArray2,nContQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ,lupri)
    ELSE !current intermediate RE(nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ) reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContQ*nPrimP*nPasses*nTUVP*nCartOrbCompQ)call ichorquit('IchorTmp2G4',-1)
#endif
       !Build REC(nContQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ)
       IF(Qsegmented)THEN
          call contractBasisSegQ(TmpArray2,TmpArray1,nPrimC,nPrimD,nPrimP*nPasses*nTUVP*nCartOrbCompQ) 
       ELSE
          call contractBasisGenQ(TmpArray2,TmpArray1,CCC,DCC,nPrimC,nPrimD,nContC,nContD,nPrimP*nPasses*nTUVP*nCartOrbCompQ) 
       ENDIF
       IF (IntPrint .GE. 25)call PrintIchorTensorREC(TmpArray1,nContQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ,lupri)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransQ)THEN
       IF(TMP1)THEN !current intermediate REC(nContQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ) reside in TmpArray1
#ifdef VAR_DEBUGICHOR
          IF(TMParray2maxsize.LT.nContQ*nPrimP*nPasses*nTUVP*nOrbCompQ)call ichorquit('IchorTmp1G5',-1)
#endif
          !Build RECS(nContQ,nPrimP,nPasses,nTUVP,nOrbCompQ)
          IF(Sph3.AND.Sph4)THEN
             call SphericalTransformGenQCD(TmpArray1,TmpArray2,SPH_MAT(AngmomC)%elms,SPH_MAT(AngmomD)%elms,&
                  & nCartOrbCompC,nOrbCompC,nCartOrbCompD,nOrbCompD,nContQ*nPrimP*nPasses*nTUVP)
          ELSEIF(Sph3)THEN
             call SphericalTransformGenQC(TmpArray1,TmpArray2,SPH_MAT(AngmomC)%elms,&
                  & nCartOrbCompC,nOrbCompC,nOrbCompD,nContQ*nPrimP*nPasses*nTUVP)
          ELSE
             call SphericalTransformGenQD(TmpArray1,TmpArray2,SPH_MAT(AngmomD)%elms,&
                  & nOrbCompC,nCartOrbCompD,nOrbCompD,nContQ*nPrimP*nPasses*nTUVP)
          ENDIF
          IF (IntPrint .GE. 25) call PrintIchorTensorRECS(TmpArray2,nContQ,nPrimP,nPasses,nTUVP,nOrbCompQ,lupri)
       ELSE !current intermediate REC(nContQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ) reside in TmpArray2
#ifdef VAR_DEBUGICHOR
          IF(TMParray1maxsize.LT.nContQ*nPrimP*nPasses*nTUVP*nOrbCompQ)call ichorquit('IchorTmp1G5',-1)
#endif
          !Build RECS(nContQ,nPrimP,nPasses,nTUVP,nOrbCompQ)
          IF(Sph3.AND.Sph4)THEN
             call SphericalTransformGenQCD(TmpArray2,TmpArray1,SPH_MAT(AngmomC)%elms,SPH_MAT(AngmomD)%elms,&
                  & nCartOrbCompC,nOrbCompC,nCartOrbCompD,nOrbCompD,nContQ*nPrimP*nPasses*nTUVP)
          ELSEIF(Sph3)THEN
             call SphericalTransformGenQC(TmpArray2,TmpArray1,SPH_MAT(AngmomC)%elms,&
                  & nCartOrbCompC,nOrbCompC,nOrbCompD,nContQ*nPrimP*nPasses*nTUVP)
          ELSE
             call SphericalTransformGenQD(TmpArray2,TmpArray1,SPH_MAT(AngmomD)%elms,&
                  & nOrbCompC,nCartOrbCompD,nOrbCompD,nContQ*nPrimP*nPasses*nTUVP)
          ENDIF
          IF (IntPrint .GE. 25) call PrintIchorTensorRECS(TmpArray1,nContQ,nPrimP,nPasses,nTUVP,nOrbCompQ,lupri)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

    !builds Ecoeff(nPrimP,nPasses,nTUVP,nCartOrbCompP)
    !currently not OpenMP parallel 
    call Ichorbuild_Ecoeff_LHS(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
         & nCartOrbCompP,Aexp,Bexp,TmpArray3,Pdistance12,Ppreexpfac,nPasses,&
         & nAtomsA,nAtomsB,IatomApass,IatomBpass,MaxPasses,intprint,lupri,TmpArray4)

    IF (IntPrint .GE. 25)call printEcoeff(TmpArray3,nTUVP,nCartOrbCompP,nPrimP,nPasses,lupri)
    
    !builds ERECS(nContQ,nPrimP,nPasses,nCartOrbCompP,nOrbCompQ)
    IF(TMP1)THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContQ*nPrimP*nPasses*nCartOrbCompP*nOrbCompQ)call ichorquit('IchorTmp2G6A',-1)
#endif
       call contractEcoeffGenP(TmpArray1,TmpArray2,TmpArray3,nOrbCompQ,nCartOrbCompP,nTUVP,nContQ,nPrimP*nPasses)

       IF (IntPrint .GE. 25)call PrintIchorTensorERECS(TmpArray2,nContQ,nPrimP,nPasses,nCartOrbCompP,nOrbCompQ,lupri)
    ELSE
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContQ*nPrimP*nPasses*nCartOrbCompP*nOrbCompQ)call ichorquit('IchorTmp1G6B',-1)
#endif
       call contractEcoeffGenP(TmpArray2,TmpArray1,TmpArray3,nOrbCompQ,nCartOrbCompP,nTUVP,nContQ,nPrimP*nPasses)
       IF (IntPrint .GE. 25)call PrintIchorTensorERECS(TmpArray1,nContQ,nPrimP,nPasses,nCartOrbCompP,nOrbCompQ,lupri)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(.NOT.TMP1)THEN !current intermediate reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContQ*nContP*nPasses*nCartOrbCompP*nOrbCompQ)call ichorquit('IchorTmp1G7',-1)
#endif
       !builds CERECS(nContQ,nContP,nPasses,nCartOrbCompP,nOrbCompQ)
       IF(Psegmented)THEN
          call contractBasisSegP(TmpArray2,TmpArray1,nContQ,nPasses*nCartOrbCompP*nOrbCompQ,nPrimA,nPrimB) 
       ELSE
          call contractBasisGenP(TmpArray2,TmpArray1,ACC,BCC,nContQ,nPasses*nCartOrbCompP*nOrbCompQ,nPrimA,nPrimB,nContA,nContB) 
       ENDIF
       IF (IntPrint .GE. 25) call PrintIchorTensorCERECS(TmpArray1,nContQ,nContP,nPasses,nCartOrbCompP,nOrbCompQ,lupri)
    ELSE !current intermediate reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContQ*nContP*nPasses*nCartOrbCompP*nOrbCompQ)call ichorquit('IchorTmp2G7',-1)
#endif
       !builds CERECS(nContQ,nContP,nPasses,nCartOrbCompP,nOrbCompQ)
       IF(Psegmented)THEN
          call contractBasisSegP(TmpArray1,TmpArray2,nContQ,nPasses*nCartOrbCompP*nOrbCompQ,nPrimA,nPrimB) 
       ELSE
          call contractBasisGenP(TmpArray1,TmpArray2,ACC,BCC,nContQ,nPasses*nCartOrbCompP*nOrbCompQ,nPrimA,nPrimB,nContA,nContB) 
       ENDIF
       IF (IntPrint .GE. 25) call PrintIchorTensorCERECS(TmpArray2,nContQ,nContP,nPasses,nCartOrbCompP,nOrbCompQ,lupri)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransP)THEN
       IF(TMP1)THEN !current intermediate reside in TmpArray1
#ifdef VAR_DEBUGICHOR
          IF(TMParray2maxsize.LT.nContQ*nContP*nPasses*nOrbCompP*nOrbCompQ)call ichorquit('IchorTmp1G8',-1)
#endif
          !builds SCEREC(nContQ,nContP,nPasses,nOrbCompP,nOrbCompQ)
          IF(Sph1.AND.Sph2)THEN
             call SphericalTransformGenPAB(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,&
                  & SPH_MAT(AngmomB)%elms,nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,&
                  & nContQ*nContP*nPasses,nOrbCompQ)
          ELSEIF(Sph1)THEN
             call SphericalTransformGenPA(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContQ*nContP*nPasses,nOrbCompQ)
          ELSE
             call SphericalTransformGenPB(TmpArray1,TmpArray2,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContQ*nContP*nPasses,nOrbCompQ)
          ENDIF
          IF(IntPrint.GE.25)call PrintIchorTensorSCERECS(TmpArray2,&
               & nContQ,nContP,nPasses,nOrbCompP,nOrbCompQ,lupri)
       ELSE !current intermediate reside in TmpArray2
#ifdef VAR_DEBUGICHOR
          IF(TMParray1maxsize.LT.nContQ*nContP*nPasses*nOrbCompP*nOrbCompQ)call ichorquit('IchorTmp2G8',-1)
#endif
          !builds SCEREC(nContQ,nContP,nPasses,nOrbCompP,nOrbCompQ)
          IF(Sph1.AND.Sph2)THEN
             call SphericalTransformGenPAB(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,&
                  & SPH_MAT(AngmomB)%elms,nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,&
                  & nContQ*nContP*nPasses,nOrbCompQ)
          ELSEIF(Sph1)THEN
             call SphericalTransformGenPA(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContQ*nContP*nPasses,nOrbCompQ)
          ELSE
             call SphericalTransformGenPB(TmpArray2,TmpArray1,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContQ*nContP*nPasses,nOrbCompQ)
          ENDIF          
          IF(IntPrint.GE.25)call PrintIchorTensorSCERECS(TmpArray1,&
               & nContQ,nContP,nPasses,nOrbCompP,nOrbCompQ,lupri)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

    !reorder From SCERECS(nContQ,nContP,nPasses,nOrbCompP,nOrbCompQ) to 
    !LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
    IF(.NOT.TMP1)THEN
       call reorderABCD(TmpArray2,CDAB,nContQ*nContP*nPasses,nOrbCompP*nOrbCompQ)
    ELSE
       call reorderABCD(TmpArray1,CDAB,nContQ*nContP*nPasses,nOrbCompP*nOrbCompQ)
    ENDIF
  end subroutine ICI_CPU_McM_general

  subroutine ICI_CPU_McM_general_size(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP,Psegmented,Qsegmented)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP
    integer,intent(in) :: nContP,nContQ,nContA,nContB,nContC,nContD
    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD,nPrimQP,nContQP
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    logical,intent(in) :: Qsegmented,Psegmented
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    !local variables
    integer :: AngmomP,AngmomQ,AngmomPQ,nTUV,nTUVQ
    integer :: ijk1,ijk2,nCartOrbCompP,ijk1s,ijk2s,nOrbCompP
    integer :: ijk3,ijk4,nCartOrbCompQ,ijk3s,ijk4s,nOrbCompQ,nTUVP,j
    logical :: Sph1,Sph2,Sph3,Sph4,SphericalTransP,SphericalTransQ
    logical :: TMP1

    TMParray1maxsize = 1
    TMParray2maxsize = 1

    AngmomP = AngmomA+AngmomB     
    AngmomQ = AngmomC+AngmomD
    AngmomPQ  = AngmomP + AngmomQ
    nTUV=(AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6
    ijk1 = (AngmomA + 1)*(AngmomA + 2)/2
    ijk2 = (AngmomB + 1)*(AngmomB + 2)/2
    ijk3 = (AngmomC + 1)*(AngmomC + 2)/2
    ijk4 = (AngmomD + 1)*(AngmomD + 2)/2
    nCartOrbCompP = ijk1*ijk2
    nCartOrbCompQ = ijk3*ijk4
    ijk1s = 2*AngmomA + 1
    ijk2s = 2*AngmomB + 1
    ijk3s = 2*AngmomC + 1
    ijk4s = 2*AngmomD + 1
    nOrbCompP = ijk1s*ijk2s
    nOrbCompQ = ijk3s*ijk4s
    nTUVP=(AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
    nTUVQ=(AngmomQ+1)*(AngmomQ+2)*(AngmomQ+3)/6
!    Sph1 = sphericalGTO.AND.(AngmomA.GT. 1)
!    Sph2 = sphericalGTO.AND.(AngmomB.GT. 1)
!    Sph3 = sphericalGTO.AND.(AngmomC.GT. 1)
!    Sph4 = sphericalGTO.AND.(AngmomD.GT. 1)
    Sph1 = AngmomA.GT. 1
    Sph2 = AngmomB.GT. 1
    Sph3 = AngmomC.GT. 1
    Sph4 = AngmomD.GT. 1
    SphericalTransP = Sph1.OR.Sph2
    SphericalTransQ = Sph3.OR.Sph4

    TMP1 = .TRUE.
    IF (AngmomPQ.EQ. 0) THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimQP*ntuv)
    ELSEIF (AngmomPQ.EQ. 1) THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimQP*ntuv)
    ELSEIF (AngmomPQ.EQ. 2) THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimQP*ntuv)
    ELSEIF (AngmomPQ.EQ. 3) THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimQP*ntuv)
    ELSE !AngmomPQ > 3
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimQP*ntuv)
       TMP1 = .TRUE.
       !J minimum goes from 0 to 0. For  AngmomPQ=6 J takes the values 2,1,0 
       DO j=AngmomPQ-4,0,-1 
          IF(TMP1)THEN
             TMParray2maxsize = MAX(TMParray2maxsize,nPrimQP*ntuv)
             TMP1 = .FALSE.
          ELSE
             TMParray1maxsize = MAX(TMParray1maxsize,nPrimQP*ntuv)
             TMP1 = .TRUE.
          ENDIF
       ENDDO
    ENDIF    
    IF(TMP1)THEN 
       TMParray2maxsize = MAX(TMParray2maxsize,nPrimQP*ntuvP*nCartOrbCompQ)
    ELSE 
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimQP*ntuvP*nCartOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1
    IF(TMP1)THEN 
       TMParray2maxsize = MAX(TMParray2maxsize,nContQ*nPrimP*nTUVP*nCartOrbCompQ)
    ELSE 
       TMParray1maxsize = MAX(TMParray1maxsize,nContQ*nPrimP*nTUVP*nCartOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransQ)THEN
       IF(TMP1)THEN 
          TMParray2maxsize = MAX(TMParray2maxsize,nContQ*nPrimP*nTUVP*nOrbCompQ)
       ELSE 
          TMParray1maxsize = MAX(TMParray1maxsize,nContQ*nPrimP*nTUVP*nOrbCompQ)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

    IF(TMP1)THEN
       TMParray2maxsize = MAX(TMParray2maxsize,nContQ*nPrimP*nCartOrbCompP*nOrbCompQ)
    ELSE
       TMParray1maxsize = MAX(TMParray1maxsize,nContQ*nPrimP*nCartOrbCompP*nOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(.NOT.TMP1)THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nContQ*nContP*nCartOrbCompP*nOrbCompQ)
    ELSE
       TMParray2maxsize = MAX(TMParray2maxsize,nContQ*nContP*nCartOrbCompP*nOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransP)THEN
       IF(TMP1)THEN 
          TMParray2maxsize = MAX(TMParray2maxsize,nContQ*nContP*nOrbCompP*nOrbCompQ)
       ELSE
          TMParray1maxsize = MAX(TMParray1maxsize,nContQ*nContP*nOrbCompP*nOrbCompQ)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF
  end subroutine ICI_CPU_McM_general_size
  
  subroutine reorderABCD(SCERECS,ABCD,ndim1,ndim2)
    implicit none
    integer,intent(in) :: ndim1,ndim2
    real(realk),intent(in) :: SCERECS(ndim1,ndim2)
    real(realk),intent(inout) :: ABCD(ndim2,ndim1)
    !
    integer :: I,J
    !$OMP DO COLLAPSE(2) PRIVATE(I,J)
    do I=1,ndim1
       do J=1,ndim2
          ABCD(J,I) = SCERECS(I,J)
       enddo
    enddo
    !$OMP END DO
  end subroutine ReorderABCD

  Subroutine contractBasisGenQ(RE,REC,CCC,DCC,nPrimC,nPrimD,nContC,nContD,ndim)
    implicit none
    integer,intent(in) :: nPrimC,nPrimD,ndim
    integer,intent(in) :: nContC,nContD
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: RE(nPrimC,nPrimD,ndim)
    real(realk),intent(inout) :: REC(nContC,nContD,ndim)
    !
    integer :: i,icC,icD,ipC,ipD
    real(realk) :: TMP,TMPD
    !$OMP DO COLLAPSE(3) PRIVATE(i,icC,icD,ipC,ipD,TMP,TMPD)
    do icd = 1,nContD
       do icC = 1,nContC
          do i = 1,ndim
             TMP = 0.0E0_realk          
             do ipD = 1,nPrimD
                TMPD = DCC(ipD,icD)
                do ipC = 1,nPrimC
                   TMP = TMP + RE(ipC,ipD,i)*CCC(ipC,icC)*TMPD 
                enddo
             enddo
             REC(icC,icD,i) = TMP  
          enddo
       enddo
    enddo
    !$OMP END DO
  end Subroutine contractBasisGenQ

  Subroutine contractBasisSegQ(RE,REC,nPrimC,nPrimD,ndim)
    implicit none
    integer,intent(in) :: nPrimC,nPrimD,ndim
    real(realk),intent(in) :: RE(nPrimC,nPrimD,ndim)
    real(realk),intent(inout) :: REC(ndim)
    !
    integer :: i,ipC,ipD
    real(realk) :: TMP
    !$OMP DO PRIVATE(i,ipC,ipD,TMP)
    do i = 1,ndim
       TMP = 0.0E0_realk          
       do ipD = 1,nPrimD
          do ipC = 1,nPrimC
             TMP = TMP + RE(ipC,ipD,i)
          enddo
       enddo
       REC(i) = TMP  
    enddo
    !$OMP END DO
  end Subroutine contractBasisSegQ

  subroutine SphericalTransformGenQCD(REC,RECS,SPHMATC,SPHMATD,ijk3,ijk3s,ijk4,ijk4s,ndim)
    implicit none
    integer,intent(in) :: ijk3,ijk3s,ijk4,ijk4s,ndim
    real(realk),intent(in) :: SPHMATC(ijk3,ijk3s),SPHMATD(ijk4,ijk4s)
    real(realk),intent(in) :: REC(ndim,ijk3,ijk4)
    real(realk),intent(inout) :: RECS(ndim,ijk3s,ijk4s)
    !
    integer :: c,d,cs,ds,i
    real(realk) :: TMP,TMPD
    !$OMP DO COLLAPSE(3) PRIVATE(cs,ds,i)
    do cs = 1,ijk3s
       do ds = 1,ijk4s
          do i = 1,ndim
             RECS(i,cs,ds) = 0.0E0_realk
          enddo
       enddo
    enddo
    !$OMP END DO
    !$OMP DO COLLAPSE(2) PRIVATE(c,d,cs,ds,i,TMP,TMPD)
    do cs = 1,ijk3s
       do ds = 1,ijk4s
          do d = 1,ijk4
             TMPD = SPHMATD(d,ds) 
             do c = 1,ijk3
                TMP = SPHMATC(c,cs)*TMPD
                do i = 1,ndim
                   RECS(i,cs,ds) = RECS(i,cs,ds) + REC(i,c,d)*TMP
                enddo
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO
  end subroutine SphericalTransformGenQCD

  subroutine SphericalTransformGenQC(REC,RECS,SPHMATC,ijk3,ijk3s,ijk4,ndim)
    implicit none
    integer,intent(in) :: ijk3,ijk3s,ijk4,ndim
    real(realk),intent(in) :: SPHMATC(ijk3,ijk3s)
    real(realk),intent(in) :: REC(ndim,ijk3,ijk4)
    real(realk),intent(inout) :: RECS(ndim,ijk3s,ijk4)
    !
    integer :: c,cs,ds,i
    real(realk) :: TMP
    !$OMP DO COLLAPSE(3) PRIVATE(cs,ds,i)
    do cs = 1,ijk3s
       do ds = 1,ijk4
          do i = 1,ndim
             RECS(i,cs,ds) = 0.0E0_realk
          enddo
       enddo
    enddo
    !$OMP END DO
    !$OMP DO COLLAPSE(2) PRIVATE(c,cs,ds,i,TMP)
    do cs = 1,ijk3s
       do ds = 1,ijk4
          do c = 1,ijk3
             TMP = SPHMATC(c,cs)
             do i = 1,ndim
                RECS(i,cs,ds) = RECS(i,cs,ds) + REC(i,c,ds)*TMP
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO
  end subroutine SphericalTransformGenQC

  subroutine SphericalTransformGenQD(REC,RECS,SPHMATD,ijk3,ijk4,ijk4s,ndim)
    implicit none
    integer,intent(in) :: ijk3,ijk4,ijk4s,ndim
    real(realk),intent(in) :: SPHMATD(ijk4,ijk4s)
    real(realk),intent(in) :: REC(ndim*ijk3,ijk4)
    real(realk),intent(inout) :: RECS(ndim*ijk3,ijk4s)
    !
    integer :: d,ds,i
    real(realk) :: TMP
    !$OMP DO COLLAPSE(2) PRIVATE(ds,i)
    do ds = 1,ijk4s
       do i = 1,ndim*ijk3
          RECS(i,ds) = 0.0E0_realk
       enddo
    enddo
    !$OMP END DO
    do ds = 1,ijk4s
       do d = 1,ijk4
          TMP = SPHMATD(d,ds) 
          !$OMP DO PRIVATE(i)
          do i = 1,ndim*ijk3
             RECS(i,ds) = RECS(i,ds) + REC(i,d)*TMP
          enddo
          !$OMP END DO
       enddo
    enddo
  end subroutine SphericalTransformGenQD

  SUBROUTINE contractEcoeffGenP(IN,OUT,Ecoeff,nOrbCompQ,nCartOrbCompP,nTUVP,nContQ,nPrimPassesP)
    implicit none
    Integer,intent(IN) :: nOrbCompQ,nCartOrbCompP,nTUVP,nContQ,nPrimPassesP
    Real(realk),intent(IN)  :: IN(nContQ,nPrimPassesP,nTUVP,nOrbCompQ)
    Real(realk),intent(IN)  :: Ecoeff(nPrimPassesP,nTUVP,nCartOrbCompP)
    Real(realk),intent(OUT) :: OUT(nContQ,nPrimPassesP,nCartOrbCompP,nOrbCompQ)
    !local variables
    Integer      :: ijkP,ijkQ,iP,iQ,ituvP
    Real(realk) :: TMP
    !
    IF(nTUVP.GT.1)THEN
       !$OMP DO COLLAPSE(3) PRIVATE(ijkP,ijkQ,iP,iQ,ituvP,TMP)
       DO ijkP = 1,nCartOrbCompP
          DO ijkQ = 1,nOrbCompQ
             DO iP = 1,nPrimPassesP
                TMP = Ecoeff(iP,1,ijkP)
                DO iQ = 1,nContQ
                   OUT(iQ,iP,ijkP,ijkQ) = IN(iQ,iP,1,ijkQ)*TMP
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !$OMP END DO
       !$OMP DO COLLAPSE(2) PRIVATE(ijkP,ijkQ,iP,iQ,ituvP,TMP)
       DO ijkP = 1,nCartOrbCompP
          DO ijkQ = 1,nOrbCompQ
             DO iTUVP = 2,nTUVP
                DO iP = 1,nPrimPassesP
                   TMP = Ecoeff(iP,iTUVP,ijkP)
                   DO iQ = 1,nContQ
                      OUT(iQ,iP,ijkP,ijkQ) = OUT(iQ,iP,ijkP,ijkQ)+IN(iQ,iP,iTUVP,ijkQ)*TMP
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !$OMP END DO
    ELSE
       !$OMP DO COLLAPSE(2) PRIVATE(ijkQ,iP,iQ,TMP)
       DO ijkQ = 1,nOrbCompQ
          DO iP = 1,nPrimPassesP
             TMP = Ecoeff(iP,1,1)
             DO iQ = 1,nContQ
                OUT(iQ,iP,1,ijkQ) = IN(iQ,iP,1,ijkQ)*TMP
             ENDDO
          ENDDO
       ENDDO
       !$OMP END DO
    ENDIF
  END SUBROUTINE contractEcoeffGenP
  
  Subroutine contractBasisGenP(ERECS,CERECS,ACC,BCC,nContQ,ndim,nPrimA,nPrimB,nContA,nContB)
    implicit none
    integer,intent(in) :: nPrimA,nPrimB,ndim,nContA,nContB,nContQ
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: ERECS(nContQ,nPrimA,nPrimB,ndim)
    real(realk),intent(inout) :: CERECS(nContQ,nContA,nContB,ndim)
    !
    integer :: i,icA,icB,ipA,ipB,iQ
    real(realk) :: TMP,TMPB
    !$OMP DO COLLAPSE(3) PRIVATE(i,icA,icB,ipA,ipB,iQ,TMP,TMPB)
    do icB = 1,nContB
       do icA = 1,nContA
          do i = 1,ndim
             do iQ = 1,nContQ
                CERECS(iQ,icA,icB,i) = 0.0E0_realk
             enddo
             do ipB = 1,nPrimB
                TMPB = BCC(ipB,icB)
                do ipA = 1,nPrimA
                   TMP = ACC(ipA,icA)*TMPB
                   do iQ = 1,nContQ
                      CERECS(iQ,icA,icB,i) = CERECS(iQ,icA,icB,i) + ERECS(iQ,ipA,ipB,i)*TMP
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO

  end Subroutine contractBasisGenP

  Subroutine contractBasisSegP(ERECS,CERECS,nContQ,ndim,nPrimA,nPrimB)
    implicit none
    integer,intent(in) :: nPrimA,nPrimB,ndim,nContQ
    real(realk),intent(in) :: ERECS(nContQ,nPrimA,nPrimB,ndim)
    real(realk),intent(inout) :: CERECS(nContQ,ndim)
    !
    integer :: i,ipA,ipB,iQ
    !$OMP DO PRIVATE(i,ipA,ipB,iQ)
    do i = 1,ndim
       do iQ = 1,nContQ
          CERECS(iQ,i) = 0.0E0_realk
       enddo
       do ipB = 1,nPrimB
          do ipA = 1,nPrimA
             do iQ = 1,nContQ
                CERECS(iQ,i) = CERECS(iQ,i) + ERECS(iQ,ipA,ipB,i)
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO
  end Subroutine contractBasisSegP

  subroutine SphericalTransformGenPAB(CERECS,SCERECS,SPHMATA,SPHMATB,&
       & ijk1,ijk1s,ijk2,ijk2s,ndim,nOrbCompQ)
    implicit none
    integer,intent(in) :: ijk1,ijk1s,ijk2,ijk2s,ndim,nOrbCompQ
    real(realk),intent(in) :: SPHMATA(ijk1,ijk1s),SPHMATB(ijk2,ijk2s)
    real(realk),intent(in) :: CERECS(ndim,ijk1,ijk2,nOrbCompQ)
    real(realk),intent(inout) :: SCERECS(ndim,ijk1s,ijk2s,nOrbCompQ)
    !
    integer :: a,b,as,bs,i,ijkQ
    real(realk) :: TMP,TMPB
    !$OMP DO COLLAPSE(3) PRIVATE(a,b,as,bs,i,ijkQ,TMP,TMPB)
    do ijkQ = 1,nOrbCompQ
     do as = 1,ijk1s
      do bs = 1,ijk2s
       do i = 1,ndim
        SCERECS(i,as,bs,ijkQ) = 0.0E0_realk
       enddo 
       do b = 1,ijk2
        TMPB = SPHMATB(b,bs)
        do a = 1,ijk1
         TMP = SPHMATA(a,as)*TMPB
         do i = 1,ndim
          SCERECS(i,as,bs,ijkQ) = SCERECS(i,as,bs,ijkQ) + CERECS(i,a,b,ijkQ)*TMP
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
    !$OMP END DO
  end subroutine SphericalTransformGenPAB

  subroutine SphericalTransformGenPA(CERECS,SCERECS,SPHMATA,&
       & ijk1,ijk1s,ijk2,ndim,nOrbCompQ)
    implicit none
    integer,intent(in) :: ijk1,ijk1s,ijk2,ndim,nOrbCompQ
    real(realk),intent(in) :: SPHMATA(ijk1,ijk1s)
    real(realk),intent(in) :: CERECS(ndim,ijk1,ijk2*nOrbCompQ)
    real(realk),intent(inout) :: SCERECS(ndim,ijk1s,ijk2*nOrbCompQ)
    !
    integer :: a,as,i,ijkQ
    real(realk) :: TMP
    !$OMP DO COLLAPSE(2) PRIVATE(a,as,i,ijkQ,TMP)
    do ijkQ = 1,nOrbCompQ*ijk2
       do as = 1,ijk1s
          do i = 1,ndim
             SCERECS(i,as,ijkQ) = 0.0E0_realk
          enddo
          do a = 1,ijk1
             TMP = SPHMATA(a,as)
             do i = 1,ndim
                SCERECS(i,as,ijkQ) = SCERECS(i,as,ijkQ) + CERECS(i,a,ijkQ)*TMP
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO
  end subroutine SphericalTransformGenPA

  subroutine SphericalTransformGenPB(CERECS,SCERECS,SPHMATB,&
       & ijk1,ijk2,ijk2s,ndim,nOrbCompQ)
    implicit none
    integer,intent(in) :: ijk1,ijk2,ijk2s,ndim,nOrbCompQ
    real(realk),intent(in) :: SPHMATB(ijk2,ijk2s)
    real(realk),intent(in) :: CERECS(ndim*ijk1,ijk2,nOrbCompQ)
    real(realk),intent(inout) :: SCERECS(ndim*ijk1,ijk2s,nOrbCompQ)
    !
    integer :: b,bs,i,ijkQ
    real(realk) :: TMP
    !$OMP DO COLLAPSE(2) PRIVATE(b,bs,i,ijkQ,TMP)
    do ijkQ = 1,nOrbCompQ
       do bs = 1,ijk2s
          do i = 1,ndim*ijk1
             SCERECS(i,bs,ijkQ) = 0.0E0_realk
          enddo
          do b = 1,ijk2
             TMP = SPHMATB(b,bs)
             do i = 1,ndim*ijk1
                SCERECS(i,bs,ijkQ) = SCERECS(i,bs,ijkQ) + CERECS(i,b,ijkQ)*TMP
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO
  end subroutine SphericalTransformGenPB

  subroutine SphericalTransformGen(CEREC,SCERECS,SPHMATA,SPHMATB,&
       & SPHMATC,SPHMATD,ijk1,ijk1s,ijk2,ijk2s,ijk3,ijk3s,ijk4,ijk4s,ndim)
    implicit none
    integer,intent(in) :: ijk1,ijk1s,ijk2,ijk2s,ijk3,ijk3s,ijk4,ijk4s,ndim
    real(realk),intent(in) :: SPHMATA(ijk1,ijk1s),SPHMATB(ijk2,ijk2s)
    real(realk),intent(in) :: SPHMATC(ijk3,ijk3s),SPHMATD(ijk4,ijk4s)
    real(realk),intent(in) :: CEREC(ndim,ijk1,ijk2,ijk3,ijk4)
    real(realk),intent(inout) :: SCERECS(ijk1s,ijk2s,ijk3s,ijk4s,ndim)
    !
    integer :: a,b,c,d,as,bs,cs,ds,i
    real(realk) :: TMP,TMPD,TMPCD,TMPBCD
    do i = 1,ndim
     do as = 1,ijk1s
      do bs = 1,ijk2s
       do cs = 1,ijk3s
        do ds = 1,ijk4s
         TMP = 0.0E0_realk
         do d = 1,ijk4
          TMPD = SPHMATD(d,ds)
          do c = 1,ijk3
           TMPCD = TMPD*SPHMATC(c,cs)
           do b = 1,ijk2
            TMPBCD = TMPCD*SPHMATB(b,bs)
            do a = 1,ijk1
             TMP = TMP + CEREC(i,a,b,c,d)*SPHMATA(a,as)*TMPBCD
            enddo
           enddo 
          enddo 
         enddo
         SCERECS(as,bs,cs,ds,i) = TMP
        enddo
       enddo
      enddo
     enddo
    enddo
   end subroutine SphericalTransformGen

  subroutine PreCalciChorSPHMAT(MaxAngmom)
    implicit none
    integer,intent(in) :: MaxAngmom
    !
    integer :: L,nLM,nXYZ,lupri,iprint
    iprint=0
    lupri=6
    allocate(SPH_MAT(0:MaxAngmom))
    do L=0,MaxAngmom
        nLM  = 2*L+1
        nXYZ = (L+1)*(L+2)/2
        allocate(SPH_MAT(L)%elms(nXYZ*nLM))
        call IchorSph_to_Cart_matrix(L,SPH_MAT(L)%elms,nLM,nXYZ,LUPRI,IPRINT)
     enddo
   end subroutine PreCalciChorSPHMAT

   subroutine FreeIchorSPHMAT()
    implicit none
    integer :: L
    do L=0,size(SPH_MAT)-1
        deallocate(SPH_MAT(L)%elms)
    enddo
    deallocate(SPH_MAT)
  end subroutine FreeIchorSPHMAT

  FUNCTION IchorFACULT(LUPRI,N)
    IMPLICIT NONE
    real(realk), PARAMETER :: D1=1E0_realk
    integer        :: N,I,LUPRI
    real(realk)    :: IchorFACULT
    IF (N .LT. 0) THEN
       WRITE (LUPRI,'(/,A,I10,/A)')&
            &         ' Argument less than zero in IchorFACULT:',N,&
            &         ' Program cannot continue.'
       CALL LSQUIT('Illegal argument in IchorFACULT',lupri)
    ELSE
       IchorFACULT = D1
       DO I = 1, N
          IchorFACULT = IchorFACULT*I
       ENDDO
    END IF
  END FUNCTION ICHORFACULT

  FUNCTION IchorFACUL2(LUPRI,N)
    IMPLICIT NONE
    real(realk), PARAMETER :: D1=1E0_realk
    real(realk)    :: IchorFACUL2
    integer :: N,I,LUPRI
    IF (N .LT. 0) THEN
       IchorFACUL2 = DFLOAT(N + 2)
       DO I = N + 4, 1, 2
          IchorFACUL2 = IchorFACUL2*I
       END DO
       IF (IchorFACUL2 .EQ. 0E0_realk) THEN
          WRITE (LUPRI,'(/,A,I10,/A)')&
               &            ' Double factorial undefined for ',N,&
               &            ' Program cannot continue.'
          CALL LSQUIT('Illegal argument in IchorFACUL2',lupri)
       ELSE
          IchorFACUL2 = D1/IchorFACUL2
       END IF
    ELSE IF (N.EQ. 0) THEN
       IchorFACUL2 = D1
    ELSE ! N > 0
       IchorFACUL2 = DFLOAT(N)
       DO I = N - 2, 1, -2
          IchorFACUL2 = IchorFACUL2*I
       END DO
    END IF
  END FUNCTION ICHORFACUL2

  FUNCTION IchorBINOM(LUPRI,I,J)
    real(realk), PARAMETER  :: D1=1E0_realk
    INTEGER :: I,J,LUPRI
    real(realk)    :: IchorBINOM
    IF (I .LT. J) THEN
       WRITE (LUPRI,'(/,A,2I5,/A)')&
            &         ' Second argument larger than first argument in IchorBINOM:',&
            &         I,J,' Program cannot continue.'
       CALL LSQUIT('Illegal arguments in IchorBINOM',lupri)
    ELSE
       IchorBINOM = ICHORFACULT(LUPRI,I)/(ICHORFACULT(LUPRI,I-J)*ICHORFACULT(LUPRI,J))
    END IF
  END FUNCTION ICHORBINOM

  FUNCTION IchorNCRT(I,J,K)
    IMPLICIT NONE
    INTEGER  :: I,J,K,IchorNCRT
    IchorNCRT = 1 + J + 2*K + (J + K)*(J + K - 1)/2
  END FUNCTION ICHORNCRT

  !> \brief calculate spherical to cartesian matrix
  !> \author T. Kjaergaard
  !> \date 2010
  !> \param L angular moment
  !> \param SCMAT the transformation matrix
  !> \param nLM the number of spherical components
  !> \param nXYZ the number of cartesian components
  !> \param lupri the logical unit number for the output file
  !> \param iprint the printlevel, determining how much output should be generated
  SUBROUTINE IchorSph_to_Cart_matrix(L,SCMAT,nLM,nXYZ,LUPRI,IPRINT)
    implicit none
    INTEGER,intent(in) :: L,nLM,nXYZ,LUPRI,IPRINT
    REAL(REALK),intent(inout) :: SCMAT(nXYZ,nLM)
    !
    Real(realk), parameter :: DM1 = -1.0E0_realk, D1 = 1.0E0_realk, D2 = 2.0E0_realk
    INTEGER     :: M1,MADR,MABS,V0
    REAL(realk) :: FACNRM,FAC3,FACTOR
    INTEGER     :: T,U,V,A,B,C,P,Q,R,X,Y,Z,M
    INTEGER     :: iLM,iXYZ

    IF(L .EQ. 0)THEN ! S
       SCMAT(1,1)=1E0_realk
    ELSEIF(L .EQ. 1)THEN ! P
       CALL LS_DZERO(SCMAT,9)
       SCMAT(1,1)=1E0_realk
       SCMAT(2,2)=1E0_realk
       SCMAT(3,3)=1E0_realk
    ELSEIF(L .GT. 1)THEN
       CALL LS_DZERO(SCMAT,nLM*nXYZ)
       DO M1 = 0, 2*L 
          M = M1 - L
          IF (L.EQ. 1) THEN
             IF (M .EQ. -1) MADR =  0  
             IF (M .EQ.  0) MADR =  1 
             IF (M .EQ.  1) MADR = -1 
          ELSE
             MADR = M
          END IF
          MABS = ABS(M)
          V0 = 0
          IF (M .LT. 0) V0 = 1 
          FACNRM = D1
          IF (M .NE. 0) FACNRM = SQRT(D2*IchorFACULT(LUPRI,L+MABS)*&
               &IchorFACULT(LUPRI,L-MABS))/(IchorFACULT(LUPRI,L)*(D2**MABS))
          FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
          FACNRM = FACNRM/SQRT(IchorFACUL2(LUPRI,2*L-1))
          DO T = 0, L - MABS, 2
             DO U = 0, T, 2
                DO V = V0, MABS, 2
                   !        almost 6.4.48 in the book
                   FAC3 = FACNRM*IchorBINOM(LUPRI,L,T/2)*IchorBINOM(LUPRI,L-T/2,MABS+T/2)&
                        &                    *IchorBINOM(LUPRI,T/2,U/2)*IchorBINOM(LUPRI,MABS,V)
                   DO A = 0, MIN(0,T+MABS-U-V) 
                      DO B = 0, MIN(0,U+V)
                         DO C = 0, MIN(0,L-T-MABS)
                            !           6.4.47 in the book
                            DO P = 0, - A, 2
                               DO Q = 0, - B, 2
                                  DO R = 0, - C, 2
                                     FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
                                          &   D2**(-A-B-C-P-Q-R-T)*FAC3
                                     X = T+MABS-U-V-2*A-P
                                     Y = U+V-2*B-Q
                                     Z = L-T-MABS-2*C-R
                                     iLM = 1 + L + MADR
                                     iXYZ = ichorNCRT(X,Y,Z)
                                     SCMAT(iXYZ,iLM) = SCMAT(iXYZ,iLM) + FACTOR 
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE IchorSph_to_Cart_matrix

  Subroutine contractBasisGen(ERE,CEREC,ACC,BCC,CCC,DCC,&
       & ndim,nPrimA,nPrimB,nPrimC,nPrimD,&
       & nContA,nContB,nContC,nContD)
    implicit none
    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD,ndim
    integer,intent(in) :: nContA,nContB,nContC,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: ERE(nPrimC,nPrimD,nPrimA,nPrimB,ndim)
    real(realk),intent(inout) :: CEREC(nContC,nContD,nContA,nContB,ndim)
    !
    integer :: i, icA,icB,icC,icD,ipA,ipB,ipC,ipD
    real(realk) :: TMP,TMPB,TMPAB,TMPABD
      do i = 1,ndim
       do icB = 1,nContB
        do icA = 1,nContA
         do icd = 1,nContD
          do icC = 1,nContC
           TMP = 0.0E0_realk          
           do ipB = 1,nPrimB
            TMPB = BCC(ipB,icB)
            do ipA = 1,nPrimA
             TMPAB = ACC(ipA,icA)*TMPB
             do ipD = 1,nPrimD
              TMPABD = DCC(ipD,icD)*TMPAB
              do ipC = 1,nPrimC
               TMP = TMP + ERE(ipC,ipD,ipA,ipB,i)*CCC(ipC,icC)*TMPABD 
              enddo
             enddo
            enddo
           enddo
           CEREC(icC,icD,icA,icB,i) = TMP  
          enddo
         enddo
        enddo
       enddo
      enddo
  end Subroutine contractBasisGen

SUBROUTINE contractEcoeffGen(IN,OUT,Ecoeff,&
     & nCartOrbCompQ,nCartOrbCompP,nTUVP,nPrimP,nPrimQ,nPasses)
implicit none
Integer,intent(IN) :: nCartOrbCompQ,nCartOrbCompP,nTUVP,nPrimP,nPrimQ,nPasses
Real(realk),intent(IN)  :: IN(nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ)
Real(realk),intent(IN)  :: Ecoeff(nPrimP,nPasses,nTUVP,nCartOrbCompP)
Real(realk),intent(OUT) :: OUT(nPrimQ,nPrimP,nPasses,nCartOrbCompP,nCartOrbCompQ)
!local variables
Integer      :: ijkP,ijkQ,iP,iQ,iPass,ituvP
Real(realk) :: tmp
!
DO ijkP = 1,nCartOrbCompP
 DO ijkQ = 1,nCartOrbCompQ
  DO iPass = 1,nPasses
   DO iP = 1,nPrimP
    TMP = Ecoeff(iP,iPass,1,ijkP)
    DO iQ = 1,nPrimQ
     OUT(iQ,iP,iPass,ijkP,ijkQ) = IN(iQ,iP,iPass,1,ijkQ)*TMP
    ENDDO
   ENDDO 
  ENDDO
  DO iTUVP = 2,nTUVP
   DO iPass = 1,nPasses
    DO iP = 1,nPrimP
     TMP = Ecoeff(iP,iPass,iTUVP,ijkP)
     DO iQ = 1,nPrimQ
      OUT(iQ,iP,iPass,ijkP,ijkQ) = OUT(iQ,iP,iPass,ijkP,ijkQ)+IN(iQ,iP,iPass,iTUVP,ijkQ)*TMP
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE contractEcoeffGen

SUBROUTINE DirectcontractEQgen(WTUV3,OUT,nPrim,nPrimPassP,nPrimQ,nTUV,ntuvP,ntuvQ,Ecoeffs,ijkQ,l1,l2,l3,l4,MAXJ)
implicit none
Integer,intent(in) :: nPrim,nPrimQ,ntuv,ntuvP,ntuvQ,nPrimPassP
Integer,intent(in) :: l1,l2,l3,l4,MAXJ,ijkQ
Real(realk),intent(in) :: Ecoeffs(nPrimQ,ntuvQ,ijkQ)
Real(realk),intent(in) :: WTUV3(nPrimQ,nPrimPassP,ntuv)
Real(realk),intent(inout) :: OUT(nPrimQ,nPrimPassP,ntuvP,ijkQ)
!local variables
!Real(realk)     :: Ec(nPrimQ)
!Real(realk) :: maxEcont,THRESHOLD,Etmp
!Integer     :: ntuvp2,ntuvq2
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvQP
Integer     :: iPrimQ,iprimPass,J,T,U,V,TUV,ijk,iOFF

Integer     :: startec,TUVQPindex(ntuvP*ntuvQ),ntuvp2,ntuvq2

ituvQP = 0
ituvQ = 0
DO jQ = 0,l3+l4
   DO tQ=jQ,0,-1
      DO uQ=jQ-tQ,0,-1
         vQ=jQ-tQ-uQ
         ituvQ = ituvQ+1
         ituvP = 0
         DO jP = 0,l1+l2
            DO tP=jP,0,-1
               DO uP=jP-tP,0,-1
                  vP=jP-tP-uP
                  ituvP = ituvP+1
                  iTUVQP = iTUVQP+1 
                  TUVQPindex(iTUVQP) = IchorTUVindexFuncFull(tP+tQ,uP+uQ,vP+vQ) 
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
ntuvp2=ituvp
ntuvq2=ituvq

!$OMP DO PRIVATE(ijk,ituvP,ituvQP,iPrimPass,iPrimQ,ituvQ,ioff)
DO ijk = 1, ijkQ
   DO ituvP = 1,ntuvP2
      ituvQP=TUVQPindex(ituvP)
      DO iPrimPass=1,nPrimPassP 
         DO iPrimQ=1,nPrimQ
            OUT(iPrimQ,iPrimPass,ituvP,ijk) = WTUV3(iPrimQ,iPrimPass,ituvQP)*Ecoeffs(iPrimQ,1,ijk)
         ENDDO
      ENDDO
   ENDDO
   DO ituvQ = 2,ntuvQ2
      IOFF = (ituvQ-1)*ntuvP2
      DO ituvP = 1,ntuvP2
         ituvQP=TUVQPindex(ituvP+IOFF)
         DO iPrimPass=1,nPrimPassP 
            DO iPrimQ=1,nPrimQ         
               OUT(iPrimQ,iPrimPass,ituvP,ijk) = OUT(iPrimQ,iPrimPass,ituvP,ijk) + & 
                    & WTUV3(iPrimQ,iPrimPass,ituvQP)*Ecoeffs(iPrimQ,ituvQ,ijk)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO   

!THRESHOLD = 1.0E-15_realk
!Do not use IF statement due to branches in GPU code - can use in OpenMP code
!   maxEcont = 0E0_realk
!   DO iPrimQ=1,nPrimQ
!      Ec(iPrimQ) = Ecoeffs(iPrimQ,1,ijk)
!      maxEcont = MAX(maxEcont,ABS(Ecoeffs(iPrimQ,1,ijk)))
!   ENDDO
!   IF(maxEcont .GE. THRESHOLD) THEN

!!$DO ijk = 1, ijkQ
!!$ ituvQP = 0
!!$ ituvQ = 1
!!$ ituvP = 0
!!$ DO jP = 0,l1+l2
!!$  DO tP=jP,0,-1
!!$   DO uP=jP-tP,0,-1
!!$    vP=jP-tP-uP
!!$    ituvP = ituvP+1
!!$    iTUVQP = iTUVQP+1 
!!$    DO iPrimPass=1,nPrimPassP 
!!$     DO iPrimQ=1,nPrimQ
!!$      OUT(iPrimQ,iPrimPass,ituvP,ijk) = WTUV3(iPrimQ,iPrimPass,ituvQP)*Ecoeffs(iPrimQ,1,ijk)!Ec(iPrimQ)
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$  ENDDO
!!$ ENDDO
!!$ DO jQ = 1,l3+l4
!!$  DO tQ=jQ,0,-1
!!$   DO uQ=jQ-tQ,0,-1
!!$    vQ=jQ-tQ-uQ
!!$    ituvQ = ituvQ+1
!!$    ituvP = 0
!!$    DO jP = 0,l1+l2
!!$     DO tP=jP,0,-1
!!$      DO uP=jP-tP,0,-1
!!$       vP=jP-tP-uP
!!$       ituvP = ituvP+1
!!$       iTUVQP = iTUVQP+1 
!!$       DO iPrimPass=1,nPrimPassP 
!!$        DO iPrimQ=1,nPrimQ
!!$         OUT(iPrimQ,iPrimPass,ituvP,ijk)=OUT(iPrimQ,iPrimPass,ituvP,ijk) + &
!!$              & WTUV3(iPrimQ,iPrimPass,ituvQP)*Ecoeffs(iPrimQ,ituvQ,ijk)!Ec(iPrimQ)
!!$        ENDDO
!!$       ENDDO
!!$      ENDDO
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$  ENDDO
!!$ ENDDO
!!$ENDDO
END SUBROUTINE DirectcontractEQgen

subroutine build_Rpq(nPrimQ,nPassP,nPrimP,Qcent,Pcent,Rpq,&
     & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB)
  implicit none
  integer,intent(in) :: nPrimQ,nPassP,nPrimP,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Qcent(3,nPrimQ),Pcent(3,nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Rpq(nPrimQ*nPrimP*nPassP,3)
  !
  integer :: iP,iPrimQ,iPassP,iPrimP,iAtomA,iAtomB
!$OMP DO PRIVATE(iP,iPrimQ,iPassP,iPrimP,iAtomA,iAtomB)
  DO iP=1,nPrimQ*nPrimP*nPassP
     iPrimQ = mod(IP-1,nPrimQ)+1
     iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
     iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
     rPQ(iP,1) = Pcent(1,iPrimP,iAtomA,iAtomB) - Qcent(1,iPrimQ)
     rPQ(iP,2) = Pcent(2,iPrimP,iAtomA,iAtomB) - Qcent(2,iPrimQ)
     rPQ(iP,3) = Pcent(3,iPrimP,iAtomA,iAtomB) - Qcent(3,iPrimQ)
  ENDDO
!$OMP END DO
end subroutine build_Rpq

  SUBROUTINE buildRJ000_general(nPasses,nPrimQ,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&
       & TABFJW,RJ000,JMAX,integralPrefactor,Rpq)
    IMPLICIT NONE
    INTEGER,intent(in)         :: nPrimQ,nPrimP,Jmax,nTABFJW1,nTABFJW2,nPasses
    REAL(REALK),intent(in)     :: reducedExponents(nPrimQ*nPrimP)
    REAL(REALK),intent(in)     :: integralPrefactor(nPrimQ*nPrimP)
    REAL(REALK),intent(in)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    real(realk),intent(in)     :: Rpq(nPrimQ*nPrimP*nPasses,3)     
    REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrimQ*nPrimP*nPasses)
    !local variables
    REAL(REALK)     :: D2JP36,WVAL
    REAL(REALK),PARAMETER :: HALF =0.5E0_realk,D1=1E0_realk,D2 = 2E0_realk, D4 = 4E0_realk, D100=100E0_realk
    Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk
    REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
    REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk
    Real(realk) :: WDIFF,RWVAL,REXPW,PREF,D2MALPHA
    REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
    REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
    REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
    REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
    Real(realk), parameter :: PI=3.14159265358979323846E0_realk
    REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
    REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
    REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
    Real(realk) :: W2,W3,R,pqx,pqy,pqz
    REAL(REALK), PARAMETER :: SMALL = 1E-15_realk
    Integer :: IPNT,J,iP,iPrimQP
!$OMP DO PRIVATE(iP,iPrimQP,WVAL,IPNT,WDIFF,W2,W3,R,J,REXPW,RWVAL,PREF,D2MALPHA,D2JP36)
    DO iP = 1,nPrimQ*nPrimP*nPasses
       D2JP36 = 2*JMAX + 36
       iPrimQP = mod(IP-1,nPrimQ*nPrimP)+1
       !(nPrimP,nPrimQ)*(nPrimP,nPrimQ,nPasses)        squaredDistance
       WVAL = reducedExponents(iPrimQP)*(Rpq(iP,1)*Rpq(iP,1)+Rpq(iP,2)*Rpq(iP,2)+Rpq(iP,3)*Rpq(iP,3))
       !  0 < WVAL < 0.000001
       IF (ABS(WVAL) .LT. SMALL) THEN         
          RJ000(0,iP) = D1
          DO J=1,JMAX
             RJ000(J,iP)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
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
             RJ000(J,iP) = R
          ENDDO
          !  12 < WVAL <= (2J+36) 
       ELSE IF (WVAL.LE.D2JP36) THEN
          REXPW = HALF*EXP(-WVAL)
          RWVAL = D1/WVAL
          RJ000(0,iP) = SQRPIH*SQRT(RWVAL) - &
               & REXPW*(GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3)))*RWVAL
          DO J=1,JMAX
             RJ000(J,iP) = RWVAL*((J - HALF)*RJ000(J-1,iP)-REXPW)
          ENDDO
          !  (2J+36) < WVAL 
       ELSE
          RWVAL = PID4/WVAL
          RJ000(0,iP) = SQRT(RWVAL)
          RWVAL = RWVAL*PID4I
          DO J = 1, JMAX
             RJ000(J,iP) = RWVAL*(J - HALF)*RJ000(J-1,iP)
          ENDDO
       ENDIF
       ! Scaling
       PREF = integralPrefactor(iPrimQP)
       RJ000(0,IP) = PREF*RJ000(0,IP)
       D2MALPHA = -D2*reducedExponents(iPrimQP)
       DO j=1,jmax
          PREF = PREF*D2MALPHA
          RJ000(J,IP) = PREF*RJ000(J,IP)
       ENDDO
    ENDDO
!$OMP END DO
  END SUBROUTINE buildRJ000_general

subroutine PrintRJ000(RJ000,AngmomPQ,nPrim,nPasses,lupri)
  implicit none
  integer,intent(in) :: AngmomPQ,nPrim,nPasses,lupri
  real(realk),intent(in) :: RJ000(0:AngmomPQ,nPrim,nPasses)
  !
  integer :: I,J,iPass
  WRITE(lupri,*)'Output from W000 nPrimQ*nPrimP=',nPrim
  WRITE(lupri,*)'Output from W000 nPasses      =',nPasses
  DO iPass=1,nPasses
     WRITE(lupri,*)'iPass = ', iPass
     DO I=1,nPrim
        DO J=0,AngmomPQ
           WRITE(LUPRI,'(2X,A6,I4,A1,I4,A2,ES16.8)')'RJ000(',J,',',I,')=',RJ000(J,I,iPass)
        ENDDO
     ENDDO
  ENDDO
end subroutine PrintRJ000

subroutine PrintIchorTensorRE(RE,nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ,lupri)
implicit none
integer,intent(in) :: nTUVP,nCartOrbCompQ,nPrimQ,nPrimP,nPasses,lupri
real(realk),intent(in) :: RE(nPrimQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ)
integer :: iTUV,iPass,iP,ijkQ,iQ
iTUV=0
WRITE(lupri,'(A)')'Print the RTUV tensor contracted with Ecoeff segmented: RE(nPrimP,nPasses,nTUVP,nCartOrbCompQ)'
WRITE(lupri,'(A,I7)')'nTUVP:   ',nTUVP
WRITE(lupri,'(A,I7)')'nCartOrbCompQ:',nCartOrbCompQ
WRITE(lupri,'(A,I7)')'nPrimQ:  ',nPrimQ
WRITE(lupri,'(A,I7)')'nPrimP:  ',nPrimP
WRITE(lupri,'(A,I7)')'nPasses: ',nPasses
WRITE(lupri,'(A4,1X,A6,1X,A5,1X,A5,1X,A12)')'ijkQ','iPrimQ','iPrimP','iPass','iTUV=1,nTUVP'
do iQ = 1,nPrimQ
 do iP = 1,nPrimP
  do iPass = 1,nPasses
   do ijkQ = 1,nCartOrbCompQ
    WRITE(LUPRI,'(I4,I7,I7,I6,5ES18.9/,(24X,5ES18.9))') ijkQ,iQ,iP,iPass,(RE(iQ,iP,iPass,iTUV,ijkQ),iTUV=1,nTUVP)
   enddo
  enddo
 enddo
enddo
end subroutine PrintIchorTensorRE

subroutine PrintIchorTensorREC(REC,nContQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ,lupri)
implicit none
integer,intent(in) :: nTUVP,nCartOrbCompQ,nContQ,nPrimP,nPasses,lupri
real(realk),intent(in) :: REC(nContQ,nPrimP,nPasses,nTUVP,nCartOrbCompQ)
integer :: iTUV,iPass,iP,ijkQ,iQ
iTUV=0
WRITE(lupri,'(A)')'Print REC'
WRITE(lupri,'(A,I7)')'nTUVP:   ',nTUVP
WRITE(lupri,'(A,I7)')'nCartOrbCompQ:',nCartOrbCompQ
WRITE(lupri,'(A,I7)')'nContQ:  ',nContQ
WRITE(lupri,'(A,I7)')'nPrimP:  ',nPrimP
WRITE(lupri,'(A,I7)')'nPasses: ',nPasses
WRITE(lupri,'(A4,1X,A6,1X,A5,1X,A5,1X,A12)')'ijkQ','iContQ','iPrimP','iPass','iTUV=1,nTUVP'
do iQ = 1,nContQ
 do iP = 1,nPrimP
  do iPass = 1,nPasses
   do ijkQ = 1,nCartOrbCompQ
    WRITE(LUPRI,'(I4,I7,I7,I6,5ES18.9/,(24X,5ES18.9))') ijkQ,iQ,iP,iPass,(REC(iQ,iP,iPass,iTUV,ijkQ),iTUV=1,nTUVP)
   enddo
  enddo
 enddo
enddo
end subroutine PrintIchorTensorREC

subroutine PrintIchorTensorRECS(RECS,nContQ,nPrimP,nPasses,nTUVP,nOrbCompQ,lupri)
implicit none
integer,intent(in) :: nTUVP,nOrbCompQ,nContQ,nPrimP,nPasses,lupri
real(realk),intent(in) :: RECS(nContQ,nPrimP,nPasses,nTUVP,nOrbCompQ)
integer :: iTUV,iPass,iP,ijkQ,iQ
iTUV=0
WRITE(lupri,'(A)')'Print RECS'
WRITE(lupri,'(A,I7)')'nTUVP:   ',nTUVP
WRITE(lupri,'(A,I7)')'nOrbCompQ: ',nOrbCompQ
WRITE(lupri,'(A,I7)')'nContQ:  ',nContQ
WRITE(lupri,'(A,I7)')'nPrimP:  ',nPrimP
WRITE(lupri,'(A,I7)')'nPasses: ',nPasses
WRITE(lupri,'(A7,1X,A6,1X,A5,1X,A5,1X,A12)')'nOrbCompQ','iContQ','iPrimP','iPass','iTUV=1,nTUVP'
do iQ = 1,nContQ
 do iP = 1,nPrimP
  do iPass = 1,nPasses
   do ijkQ = 1,nOrbCompQ
    WRITE(LUPRI,'(I7,I7,I7,I6,5ES18.9/,(24X,5ES18.9))') ijkQ,iQ,iP,iPass,(RECS(iQ,iP,iPass,iTUV,ijkQ),iTUV=1,nTUVP)
   enddo
  enddo
 enddo
enddo
end subroutine PrintIchorTensorRECS

subroutine PrintIchorTensorERECS(ERECS,nContQ,nPrimP,nPasses,nCartOrbCompP,nOrbCompQ,lupri)
implicit none
integer,intent(in) :: nOrbCompQ,nCartOrbCompP,nContQ,nPrimP,nPasses,lupri
real(realk),intent(in) :: ERECS(nContQ,nPrimP,nPasses,nCartOrbCompP,nOrbCompQ)
!
integer :: iPass,iP,ijkQ,ijkP,iQ
WRITE(lupri,'(A)')'Print ERECS'
WRITE(lupri,'(A,I7)')'nCartOrbCompP:',nCartOrbCompP
WRITE(lupri,'(A,I7)')'nOrbCompQ: ',nOrbCompQ
WRITE(lupri,'(A,I7)')'nContQ:  ',nContQ
WRITE(lupri,'(A,I7)')'nPrimP:  ',nPrimP
WRITE(lupri,'(A,I7)')'nPasses: ',nPasses
WRITE(lupri,'(A7,1X,A6,1X,A5,1X,A5,1X,A12)')'nOrbCompQ','iContQ','iPrimP','iPass','ijkP=1,nCartOrbCompP'
do iQ = 1,nContQ
 do iP = 1,nPrimP
  do iPass = 1,nPasses
   do ijkQ = 1,nOrbCompQ
    WRITE(LUPRI,'(I7,I7,I7,I6,5ES18.9/,(24X,5ES18.9))') ijkQ,iQ,iP,iPass,(ERECS(iQ,iP,iPass,ijkP,ijkQ),ijkP=1,nCartOrbCompP)
   enddo
  enddo
 enddo
enddo
end subroutine PrintIchorTensorERECS

subroutine PrintIchorTensorCERECS(CERECS,nContQ,nContP,nPasses,nCartOrbCompP,nOrbCompQ,lupri)
implicit none
integer,intent(in) :: nOrbCompQ,nCartOrbCompP,nContQ,nContP,nPasses,lupri
real(realk),intent(in) :: CERECS(nContQ,nContP,nPasses,nCartOrbCompP,nOrbCompQ)
!
integer :: iPass,iP,ijkQ,ijkP,iQ
WRITE(lupri,'(A)')'Print CERECS'
WRITE(lupri,'(A,I7)')'nCartOrbCompP:',nCartOrbCompP
WRITE(lupri,'(A,I7)')'nOrbCompQ:',nOrbCompQ
WRITE(lupri,'(A,I7)')'nContQ:  ',nContQ
WRITE(lupri,'(A,I7)')'nContP:  ',nContP
WRITE(lupri,'(A,I7)')'nPasses: ',nPasses
WRITE(lupri,'(A7,1X,A6,1X,A5,1X,A5,1X,A12)')'nOrbCompQ','iContQ','iContP','iPass','ijkP=1,nCartOrbCompP'
do iQ = 1,nContQ
 do iP = 1,nContP
  do iPass = 1,nPasses
   do ijkQ = 1,nOrbCompQ
    WRITE(LUPRI,'(I7,I7,I7,I6,5ES18.9/,(24X,5ES18.9))') ijkQ,iQ,iP,iPass,(CERECS(iQ,iP,iPass,ijkP,ijkQ),ijkP=1,nCartOrbCompP)
   enddo
  enddo
 enddo
enddo
end subroutine PrintIchorTensorCERECS

subroutine PrintIchorTensorSCERECS(SCERECS,nContQ,nContP,nPasses,nOrbCompP,nOrbCompQ,lupri)
implicit none
integer,intent(in) :: nOrbCompQ,nOrbCompP,nContQ,nContP,nPasses,lupri
real(realk),intent(in) :: SCERECS(nContQ,nContP,nPasses,nOrbCompP,nOrbCompQ)
!
integer :: iPass,iP,ijkQ,ijkP,iQ
WRITE(lupri,'(A)')'Print SCERECS'
WRITE(lupri,'(A,I7)')'nOrbCompP:',nOrbCompP
WRITE(lupri,'(A,I7)')'nOrbCompQ:',nOrbCompQ
WRITE(lupri,'(A,I7)')'nContQ:  ',nContQ
WRITE(lupri,'(A,I7)')'nContP:  ',nContP
WRITE(lupri,'(A,I7)')'nPasses: ',nPasses
WRITE(lupri,'(A7,1X,A6,1X,A5,1X,A5,1X,A12)')'nOrbCompQ','iContQ','iContP','iPass','ijkP=1,nOrbCompP'
do iQ = 1,nContQ
 do iP = 1,nContP
  do iPass = 1,nPasses
   do ijkQ = 1,nOrbCompQ
    WRITE(LUPRI,'(I7,I7,I7,I6,5ES18.9/,(24X,5ES18.9))') ijkQ,iQ,iP,iPass,(SCERECS(iQ,iP,iPass,ijkP,ijkQ),ijkP=1,nOrbCompP)
   enddo
  enddo
 enddo
enddo
end subroutine PrintIchorTensorSCERECS

end MODULE IchorEriCoulombintegralCPUMcMGeneralMod
