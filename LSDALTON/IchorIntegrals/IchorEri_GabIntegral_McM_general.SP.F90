!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
module SPIchorEriGabintegralCPUMcMGeneralMod
use SPIchorEriCoulombintegralCPUMcMGeneralMod
use IchorPrecisionMod
use IchorMemory
use IchorCommonMod
use SPIchorEriCoulombintegralCPUMcMGeneralEcoeffMod, only: &
     & Ichorbuild_Ecoeff_RHS,Ichorbuild_Ecoeff_LHS, printEcoeff
use SPIchorEriCoulombintegralCPUMcMGeneralWTUVMod
use IchorGaussianGeminalMod
use IchorParametersMod

CONTAINS
  subroutine SPIGI_CPU_McM_general(nPrimA,nPrimB,&
       & nPrimP,IntPrint,lupri,&
       & nContA,nContB,nContP,pexp,ACC,BCC,&
       & nOrbCompA,nOrbCompB,&
       & nCartOrbCompA,nCartOrbCompB,&
       & nCartOrbCompP,nOrbCompP,nTUVP,nTUV,&
       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,Pdistance12,PQorder,&
       & CDAB,Acenter,Bcenter,&
       & SphericalGTO,TmpArray1,TMParray1maxsize,TmpArray2,&
       & TMParray2maxsize)
    implicit none
    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: nPrimP,nPrimA,nPrimB,nTABFJW1,nTABFJW2
    integer,intent(in) :: IntPrint,lupri,nContA,nContB,nContP
    integer,intent(in) :: AngmomA,AngmomB
    integer,intent(in) :: nOrbCompA,nOrbCompB
    integer,intent(in) :: nCartOrbCompA,nCartOrbCompB
    integer,intent(in) :: nCartOrbCompP,nOrbCompP,nTUVP,nTUV
    real(reals),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    real(reals),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)
    logical,intent(in)     :: Psegmented,sphericalGTO,PQorder
    real(reals),intent(in) :: Acenter(3),Bcenter(3)
    real(reals),intent(in) :: pexp(nPrimP),pcent(3*nPrimP),PpreExpFac(nPrimP)
    real(reals),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    real(reals),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(reals),intent(inout) :: CDAB(1)    
    real(reals),intent(in) :: integralPrefactor(nPrimP*nPrimP)!integralPrefactor(nPrimQ,nPrimP)    
    real(reals),intent(in) :: reducedExponents(nPrimP*nPrimP) !reducedExponents(nPrimQ,nPrimP)
    real(reals),intent(in) :: Pdistance12(3) 
    !local variables
    integer :: AngmomP,AngmomPP
    logical :: RHS,TMP1
    logical :: Sph1,Sph2,SphericalTransP
    integer :: AngmomID,TUV,J,T,U,V,I
    integer,parameter :: nPasses=1,nPassQ=1
    integer :: IatomApass(1),IatomBpass(1)
    IatomApass(1) = 1
    IatomBpass(1) = 1
    AngmomP = AngmomA+AngmomB     
    AngmomPP  = AngmomP + AngmomP
    Sph1 = sphericalGTO.AND.(AngmomA.GT. 1)
    Sph2 = sphericalGTO.AND.(AngmomB.GT. 1)
    SphericalTransP = Sph1.OR.Sph2
#ifdef VAR_DEBUGICHOR
    IF(nTmpArray3.LT.nPrimP*nPrimP*3)call ichorquit('IchorTmp0G1',-1)
#endif
    !Rpq(nPrimP,nPrimP,nPasses,3)
    call SPbuild_Gab_Rpq(nPrimP,Pcent,TmpArray3)
    !
    !      builds RJ000(0:AngmomPP,nPrimP,nPrimP) Store in TmpArray4
    !
#ifdef VAR_DEBUGICHOR
    IF(nTmpArray4.LT.nPrimP*nPrimP*(AngmomPP+1))call ichorquit('IchorTmp1G1',-1)
#endif
    IF(.NOT.GGemOperatorCalc)THEN
       call SPbuildRJ000_general(nPasses,nPrimP,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&
            & TABFJW,TmpArray4,AngmomPP,integralPrefactor,TmpArray3)
    ELSE
       IF(GGemOperatorSpec.EQ.GGemOperator.OR.GGemOperatorSpec.EQ.GGemSqOperator)THEN
          call SPbuildGJ000_general(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TmpArray4,AngmomPP,integralPrefactor,TmpArray3)
       ELSEIF(GGemOperatorSpec.EQ.GGemGrdOperator)THEN
          call SPbuildGJ000Grad_general(nPasses,nPrimP,nPrimP,reducedExponents,&
               & TmpArray4,AngmomPP,integralPrefactor,TmpArray3)
       ELSEIF(GGemOperatorSpec.EQ.GGemCouOperator)THEN
          call SPbuildGJ000Coulomb_general(nPasses,nPrimP,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&
               & TABFJW,TmpArray4,AngmomPP,integralPrefactor,TmpArray3)
       ENDIF
    ENDIF
    IF (INTPRINT .GE. 10) call SPPrintRJ000(TmpArray4,AngmomPP,nPrimP*nPrimP,nPasses,lupri)    
    !
    !     Build WTUV(nPrimP,nPrimP,nPasses,nTUV) RJ000 = TmpArray4, Rpq = TmpArray3
    !
    TMP1 = .TRUE.
    IF (AngmomPP.EQ. 0) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3A',-1)
#endif
       call SPIchorwtuvRecurrenceJMIN0JMAX0(TmpArray4,TmpArray1,nPrimP*nPrimP)
    ELSEIF (AngmomPP.EQ. 1) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3B',-1)
#endif
       call SPIchorwtuvRecurrenceJMIN0JMAX1(TmpArray4,TmpArray1,TmpArray3,nPrimP*nPrimP)
    ELSEIF (AngmomPP.EQ. 2) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3C',-1)
#endif
       call SPIchorwtuvRecurrenceJMIN0JMAX2(TmpArray4,TmpArray1,TmpArray3,nPrimP*nPrimP)
    ELSEIF (AngmomPP.EQ. 3) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3D',-1)
#endif
       call SPIchorwtuvRecurrenceJMIN0JMAX3(TmpArray4,TmpArray1,TmpArray3,nPrimP*nPrimP)
    ELSE !AngmomPP > 3
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3a',-1)
       IF(TMParray2maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3b',-1)
#endif
       J = AngmomPP-3
       call SPIchorwtuvRecurrenceJMIN0JMAX3J(TmpArray4,TmpArray1,TmpArray3,nPrimP*nPrimP,AngmomPP,J)
       TMP1 = .TRUE.
       !J minimum goes from 0 to 0. For  AngmomPP=6 J takes the values 2,1,0 
       DO j=AngmomPP-4,0,-1 
          IF(TMP1)THEN
             call SPIchorwtuvRecurrenceCurrent(TmpArray1,TmpArray2,J,AngmomPP,nPrimP,nPrimP,nPasses,&
                  & ntuv,TmpArray3,TmpArray4)
          ELSE
             call SPIchorwtuvRecurrenceCurrent(TmpArray2,TmpArray1,J,AngmomPP,nPrimP,nPrimP,nPasses,&
                  & ntuv,TmpArray3,TmpArray4)
          ENDIF
          TMP1 = .NOT.TMP1
       ENDDO
    ENDIF    
    IF (IntPrint .GE. 25)THEN
       IF(TMP1)THEN
          call SPPrintWTUV(TmpArray1,AngmomPP,nPrimP*nPrimP,nPasses,nTUV,lupri)
       ELSE
          call SPPrintWTUV(TmpArray2,AngmomPP,nPrimP*nPrimP,nPasses,nTUV,lupri)
       ENDIF
    ENDIF

    !builds Ecoeff(nPrimP,nPasses,nTUVP,nCartOrbCompP)
    IF(TMP1)THEN !use TmpArray2 to store some tmp arrays in Ichorbuild_Ecoeff
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.9*nPrimP)call ichorquit('IchorTmp1G3QEE1R2',-1)
#endif
       call SPIchorbuild_Ecoeff_RHS(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
            & nCartOrbCompP,Aexp,Bexp,TmpArray3,Pdistance12,Ppreexpfac,intprint,lupri,TmpArray4,&
            & TmpArray2(1:3*nPrimP),TmpArray2(3*nPrimP+1:6*nPrimP),&
            & TmpArray2(6*nPrimP+1:7*nPrimP),TmpArray2(7*nPrimP+1:8*nPrimP),&
            & TmpArray2(8*nPrimP+1:9*nPrimP))
    ELSE
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.9*nPrimP)call ichorquit('IchorTmp1G3QEE1R2',-1)
#endif
       call SPIchorbuild_Ecoeff_RHS(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
            & nCartOrbCompP,Aexp,Bexp,TmpArray3,Pdistance12,Ppreexpfac,intprint,lupri,TmpArray4,&
            & TmpArray1(1:3*nPrimP),TmpArray1(3*nPrimP+1:6*nPrimP),&
            & TmpArray1(6*nPrimP+1:7*nPrimP),TmpArray1(7*nPrimP+1:8*nPrimP),&
            & TmpArray1(8*nPrimP+1:9*nPrimP))

    ENDIF

    IF (IntPrint .GE. 25)call SPprintEcoeff(TmpArray3,nTUVP,nCartOrbCompP,nPrimP,nPassQ,lupri)

    IF(TMP1)THEN !current intermediate WTUV reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nPrimP*nPrimP*ntuvP*nCartOrbCompP)call ichorquit('IchorTmp1G3P1',-1)
#endif
       !builds RE(nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP)
       call SPDirectcontractEQgen(TmpArray1,TmpArray2,nPrimP*nPrimP,nPrimP,nPrimP,nTUV,ntuvP,ntuvP,&
            & TmpArray3,nCartOrbCompP,AngmomA,AngmomB,AngmomA,AngmomB,AngmomPP)
       IF (IntPrint .GE. 25)call SPPrintIchorTensorRE(TmpArray2,nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP,lupri)
    ELSE !current intermediate WTUV reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuvP*nCartOrbCompP)call ichorquit('IchorTmp1G3Q2',-1)
#endif
       !builds RE(nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP)
       call SPDirectcontractEQgen(TmpArray2,TmpArray1,nPrimP*nPrimP,nPrimP,nPrimP,nTUV,ntuvP,ntuvP,&
            & TmpArray3,nCartOrbCompP,AngmomA,AngmomB,AngmomA,AngmomB,AngmomPP)
       IF (IntPrint .GE. 25)call SPPrintIchorTensorRE(TmpArray1,nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP,lupri)
    ENDIF
    TMP1 = .NOT.TMP1
    
    IF(TMP1)THEN !current intermediate RE(nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP) reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContP*nPrimP*nTUVP*nCartOrbCompP)call ichorquit('IchorTmp2G4',-1)
#endif
       !Build REC(nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP)
       IF(Psegmented)THEN
          call SPcontractBasisSegQ(TmpArray1,TmpArray2,nPrimA,nPrimB,nPrimP*nTUVP*nCartOrbCompP) 
       ELSE
          call SPcontractBasisGenQ(TmpArray1,TmpArray2,ACC,BCC,nPrimA,nPrimB,nContA,nContB,nPrimP*nTUVP*nCartOrbCompP) 
       ENDIF
       IF (IntPrint .GE. 25)call SPPrintIchorTensorREC(TmpArray2,nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP,lupri)
    ELSE !current intermediate RE(nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP) reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContP*nPrimP*nTUVP*nCartOrbCompP)call ichorquit('IchorTmp2G4',-1)
#endif
       !Build REC(nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP)
       IF(Psegmented)THEN
          call SPcontractBasisSegQ(TmpArray2,TmpArray1,nPrimA,nPrimB,nPrimP*nTUVP*nCartOrbCompP) 
       ELSE
          call SPcontractBasisGenQ(TmpArray2,TmpArray1,ACC,BCC,nPrimA,nPrimB,nContA,nContB,nPrimP*nTUVP*nCartOrbCompP) 
       ENDIF
       IF (IntPrint .GE. 25)call SPPrintIchorTensorREC(TmpArray1,nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP,lupri)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransP)THEN
       IF(TMP1)THEN !current intermediate REC(nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP) reside in TmpArray1
#ifdef VAR_DEBUGICHOR
          IF(TMParray2maxsize.LT.nContP*nPrimP*nTUVP*nOrbCompP)call ichorquit('IchorTmp1G5',-1)
#endif
          !Build RECS(nContP,nPrimP,nPasses,nTUVP,nOrbCompP)
          IF(Sph1.AND.Sph2)THEN
             call SPSphericalTransformGenQCD(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,SPH_MAT(AngmomB)%elms,&
                  & nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nPrimP*nTUVP)
          ELSEIF(Sph1)THEN
             call SPSphericalTransformGenQC(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContP*nPrimP*nTUVP)
          ELSE
             call SPSphericalTransformGenQD(TmpArray1,TmpArray2,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nPrimP*nTUVP)
          ENDIF
          IF (IntPrint .GE. 25) call SPPrintIchorTensorRECS(TmpArray2,nContP,nPrimP,nPasses,nTUVP,nOrbCompP,lupri)
       ELSE !current intermediate REC(nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP) reside in TmpArray2
#ifdef VAR_DEBUGICHOR
          IF(TMParray1maxsize.LT.nContP*nPrimP*nTUVP*nOrbCompP)call ichorquit('IchorTmp1G5',-1)
#endif
          !Build RECS(nContP,nPrimP,nPasses,nTUVP,nOrbCompP)
          IF(Sph1.AND.Sph2)THEN
             call SPSphericalTransformGenQCD(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,SPH_MAT(AngmomB)%elms,&
                  & nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nPrimP*nTUVP)
          ELSEIF(Sph1)THEN
             call SPSphericalTransformGenQC(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContP*nPrimP*nTUVP)
          ELSE
             call SPSphericalTransformGenQD(TmpArray2,TmpArray1,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nPrimP*nTUVP)
          ENDIF
          IF (IntPrint .GE. 25) call SPPrintIchorTensorRECS(TmpArray1,nContP,nPrimP,nPasses,nTUVP,nOrbCompP,lupri)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

    !builds Ecoeff(nPrimP,nPasses,nTUVP,nCartOrbCompP)
    IF(TMP1)THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.9*nPrimP)call ichorquit('IchorTmp1G3QEE1R2',-1)
#endif
       call SPIchorbuild_Ecoeff_LHS(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
            & nCartOrbCompP,Aexp,Bexp,TmpArray3,Pdistance12,Ppreexpfac,nPasses,&
            & 1,1,IatomApass,IatomBpass,1,intprint,lupri,TmpArray4,&
            & TmpArray2(1:3*nPrimP),TmpArray2(3*nPrimP+1:6*nPrimP),&
            & TmpArray2(6*nPrimP+1:7*nPrimP),TmpArray2(7*nPrimP+1:8*nPrimP),&
            & TmpArray2(8*nPrimP+1:9*nPrimP))
    ELSE
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.9*nPrimP)call ichorquit('IchorTmp1G3QEE1R2',-1)
#endif
       call SPIchorbuild_Ecoeff_LHS(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
            & nCartOrbCompP,Aexp,Bexp,TmpArray3,Pdistance12,Ppreexpfac,nPasses,&
            & 1,1,IatomApass,IatomBpass,1,intprint,lupri,TmpArray4,&
            & TmpArray1(1:3*nPrimP),TmpArray1(3*nPrimP+1:6*nPrimP),&
            & TmpArray1(6*nPrimP+1:7*nPrimP),TmpArray1(7*nPrimP+1:8*nPrimP),&
            & TmpArray1(8*nPrimP+1:9*nPrimP))
    ENDIF


    IF (IntPrint .GE. 25)call SPprintEcoeff(TmpArray3,nTUVP,nCartOrbCompP,nPrimP,nPasses,lupri)
    
    !builds ERECS(nContP,nPrimP,nPasses,nCartOrbCompP,nOrbCompP)
    IF(TMP1)THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContP*nPrimP*nCartOrbCompP*nOrbCompP)call ichorquit('IchorTmp2G6A',-1)
#endif
       call SPcontractEcoeffGenP(TmpArray1,TmpArray2,TmpArray3,nOrbCompP,nCartOrbCompP,nTUVP,nContP,nPrimP)

       IF (IntPrint .GE. 25)call SPPrintIchorTensorERECS(TmpArray2,nContP,nPrimP,nPasses,nCartOrbCompP,nOrbCompP,lupri)
    ELSE
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContP*nPrimP*nCartOrbCompP*nOrbCompP)call ichorquit('IchorTmp1G6B',-1)
#endif
       call SPcontractEcoeffGenP(TmpArray2,TmpArray1,TmpArray3,nOrbCompP,nCartOrbCompP,nTUVP,nContP,nPrimP)
       IF (IntPrint .GE. 25)call SPPrintIchorTensorERECS(TmpArray1,nContP,nPrimP,nPasses,nCartOrbCompP,nOrbCompP,lupri)
    ENDIF
    TMP1 = .NOT.TMP1

    !note we build CERECS(nContP,nCartOrbCompP,nOrbCompP)
    !normally we build CERECS(nContP,nContP,nPasses,nCartOrbCompP,nOrbCompP)
    !but we are only intrested in the diagonal elements iContP = iContP



    !builds CERECS(nContP,nCartOrbCompP,nOrbCompP)
    IF(.NOT.TMP1)THEN !current intermediate reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContP*nContP*nCartOrbCompP*nOrbCompP)call ichorquit('IchorTmp1G7',-1)
#endif
       IF(Psegmented)THEN
          call SPGabcontractBasisSegP(TmpArray2,TmpArray1,nContP,nCartOrbCompP*nOrbCompP,nPrimA,nPrimB) 
       ELSE
          call SPGabcontractBasisGenP(TmpArray2,TmpArray1,ACC,BCC,nContP,nCartOrbCompP*nOrbCompP,nPrimA,nPrimB,nContA,nContB) 
       ENDIF
    ELSE !current intermediate reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContP*nContP*nPasses*nCartOrbCompP*nOrbCompP)call ichorquit('IchorTmp2G7',-1)
#endif
       IF(Psegmented)THEN
          call SPGabcontractBasisSegP(TmpArray1,TmpArray2,nContP,nCartOrbCompP*nOrbCompP,nPrimA,nPrimB) 
       ELSE
          call SPGabcontractBasisGenP(TmpArray1,TmpArray2,ACC,BCC,nContP,nCartOrbCompP*nOrbCompP,nPrimA,nPrimB,nContA,nContB) 
       ENDIF
    ENDIF
    TMP1 = .NOT.TMP1

    !builds SCEREC(nContP,nOrbCompP)
    IF(SphericalTransP)THEN
       IF(TMP1)THEN !current intermediate reside in TmpArray1
#ifdef VAR_DEBUGICHOR
          IF(TMParray2maxsize.LT.nContP*nContP*nOrbCompP*nOrbCompP)call ichorquit('IchorTmp1G8',-1)
#endif
          IF(Sph1.AND.Sph2)THEN
             call SPGabSphericalTransformGenPAB(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,&
                  & SPH_MAT(AngmomB)%elms,nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,&
                  & nContP,nOrbCompP)
          ELSEIF(Sph1)THEN
             call SPGabSphericalTransformGenPA(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContP,nOrbCompP)
          ELSE
             call SPGabSphericalTransformGenPB(TmpArray1,TmpArray2,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContP,nOrbCompP)
          ENDIF
       ELSE !current intermediate reside in TmpArray2
#ifdef VAR_DEBUGICHOR
          IF(TMParray1maxsize.LT.nContP*nContP*nOrbCompP*nOrbCompP)call ichorquit('IchorTmp2G8',-1)
#endif
          !builds SCEREC(nContP,nContP,nPasses,nOrbCompP,nOrbCompP)
          IF(Sph1.AND.Sph2)THEN
             call SPGabSphericalTransformGenPAB(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,&
                  & SPH_MAT(AngmomB)%elms,nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,&
                  & nContP,nOrbCompP)
          ELSEIF(Sph1)THEN
             call SPGabSphericalTransformGenPA(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContP,nOrbCompP)
          ELSE
             call SPGabSphericalTransformGenPB(TmpArray2,TmpArray1,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContP,nOrbCompP)
          ENDIF          
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

    IF(.NOT.TMP1)THEN
       call SPextractGabElmGen(TmpArray2,CDAB,nContP,nOrbCompP)
    ELSE
       call SPextractGabElmGen(TmpArray1,CDAB,nContP,nOrbCompP)       
    ENDIF
    
  end subroutine SPIGI_CPU_McM_general

  subroutine SPIGI_CPU_McM_general_size(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,nPrimP,nContP,nPrimB,Psegmented)
    implicit none
    integer,intent(in) :: nPrimP,nContP
    integer,intent(in) :: nPrimB,AngmomA,AngmomB
    logical,intent(in) :: Psegmented
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
    AngmomQ = AngmomP
    AngmomPQ  = AngmomP + AngmomQ
    nTUV=(AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6
    ijk1 = (AngmomA + 1)*(AngmomA + 2)/2
    ijk2 = (AngmomB + 1)*(AngmomB + 2)/2
    ijk3 = ijk1
    ijk4 = ijk2
    nCartOrbCompP = ijk1*ijk2
    nCartOrbCompQ = ijk3*ijk4
    ijk1s = 2*AngmomA + 1
    ijk2s = 2*AngmomB + 1
    ijk3s = ijk1s
    ijk4s = ijk2s
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
    Sph3 = Sph1
    Sph4 = Sph2
    SphericalTransP = Sph1.OR.Sph2
    SphericalTransQ = Sph3.OR.Sph4
    !angmomPP
    TMP1 = .TRUE.
    IF (AngmomPQ.EQ. 0) THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimP*nPrimP*ntuv)
    ELSEIF (AngmomPQ.EQ. 1) THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimP*nPrimP*ntuv)
    ELSEIF (AngmomPQ.EQ. 2) THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimP*nPrimP*ntuv)
    ELSEIF (AngmomPQ.EQ. 3) THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimP*nPrimP*ntuv)
    ELSE !AngmomPQ > 3
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimP*nPrimP*ntuv)
       TMP1 = .TRUE.
       !J minimum goes from 0 to 0. For  AngmomPQ=6 J takes the values 2,1,0 
       DO j=AngmomPQ-4,0,-1 
          IF(TMP1)THEN
             TMParray2maxsize = MAX(TMParray2maxsize,nPrimP*nPrimP*ntuv)
             TMP1 = .FALSE.
          ELSE
             TMParray1maxsize = MAX(TMParray1maxsize,nPrimP*nPrimP*ntuv)
             TMP1 = .TRUE.
          ENDIF
       ENDDO
    ENDIF  
    !builds Ecoeff(nPrimP,nPasses,nTUVP,nCartOrbCompP)
    IF(TMP1)THEN !use TmpArray2 to store some tmp arrays in Ichorbuild_Ecoeff
       TMParray2maxsize = MAX(TMParray2maxsize,9*nPrimP)
    ELSE
       TMParray1maxsize = MAX(TMParray1maxsize,9*nPrimP)
    ENDIF
    !DirectcontractEQgen  
    IF(TMP1)THEN 
       TMParray2maxsize = MAX(TMParray2maxsize,nPrimP*nPrimP*ntuvP*nCartOrbCompQ)
    ELSE 
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimP*nPrimP*ntuvP*nCartOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1
    !contractBasis
    IF(TMP1)THEN 
       TMParray2maxsize = MAX(TMParray2maxsize,nContP*nPrimP*nTUVP*nCartOrbCompQ)
    ELSE 
       TMParray1maxsize = MAX(TMParray1maxsize,nContP*nPrimP*nTUVP*nCartOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1
    !SphericalTrans
    IF(SphericalTransQ)THEN
       IF(TMP1)THEN 
          TMParray2maxsize = MAX(TMParray2maxsize,nContP*nPrimP*nTUVP*nOrbCompQ)
       ELSE 
          TMParray1maxsize = MAX(TMParray1maxsize,nContP*nPrimP*nTUVP*nOrbCompQ)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

    !builds Ecoeff(nPrimP,nPasses,nTUVP,nCartOrbCompP)
    IF(TMP1)THEN
       TMParray2maxsize = MAX(TMParray2maxsize,9*nPrimP)
    ELSE
       TMParray1maxsize = MAX(TMParray1maxsize,9*nPrimP)
    ENDIF

    !contractEcoeff
    IF(TMP1)THEN
       TMParray2maxsize = MAX(TMParray2maxsize,nContP*nPrimP*nCartOrbCompP*nOrbCompQ)
    ELSE
       TMParray1maxsize = MAX(TMParray1maxsize,nContP*nPrimP*nCartOrbCompP*nOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(.NOT.TMP1)THEN
       TMParray1maxsize = MAX(TMParray1maxsize,nContP*nContP*nCartOrbCompP*nOrbCompQ)
    ELSE
       TMParray2maxsize = MAX(TMParray2maxsize,nContP*nContP*nCartOrbCompP*nOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransP)THEN
       IF(TMP1)THEN 
          TMParray2maxsize = MAX(TMParray2maxsize,nContP*nContP*nOrbCompP*nOrbCompQ)
       ELSE
          TMParray1maxsize = MAX(TMParray1maxsize,nContP*nContP*nOrbCompP*nOrbCompQ)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF
  end subroutine SPIGI_CPU_McM_general_size
  
  subroutine SPbuild_Gab_Rpq(nPrimP,Pcent,Rpq)
    implicit none
    integer,intent(in) :: nPrimP
    real(reals),intent(in) :: Pcent(3,nPrimP)
    real(reals),intent(inout) :: Rpq(nPrimP,nPrimP,3)
    !local variables
    integer :: iPrimQ,iPrimP
    !$OMP DO COLLAPSE(2) PRIVATE(iPrimQ,iPrimP)
    DO iPrimQ=1,nPrimP
       DO iPrimP=1,nPrimP
          rPQ(iPrimQ,iPrimP,1) = Pcent(1,iPrimP) - Pcent(1,iPrimQ)
          rPQ(iPrimQ,iPrimP,2) = Pcent(2,iPrimP) - Pcent(2,iPrimQ)
          rPQ(iPrimQ,iPrimP,3) = Pcent(3,iPrimP) - Pcent(3,iPrimQ)
       ENDDO
    ENDDO
    !$OMP END DO
  end subroutine SPbuild_Gab_Rpq
  
  subroutine SPextractGabElmGen(TmpArray2,CDAB,nContP,nOrbCompP)
    implicit none
    integer,intent(in) :: nContP,nOrbCompP
    real(reals),intent(in) :: TmpArray2(nContP,nOrbCompP)
    real(reals),intent(inout) :: CDAB(1)
    !
    integer :: iContP,iOrbP
    real(reals) :: TMP
    !$OMP MASTER
    TMP = ABS(TmpArray2(1,1))
    do iOrbP = 1,nOrbCompP
       do iContP = 1,nContP
          TMP = MAX(TMP,ABS(TmpArray2(iContP,iOrbP)))
       enddo
    enddo
    CDAB(1) = SQRT(TMP)
    !$OMP END MASTER
    !$OMP BARRIER 

  end subroutine SPextractGabElmGen

  subroutine SPGabcontractBasisGenP(ERECS,CERECS,ACC,BCC,nContQ,ndim,nPrimA,nPrimB,nContA,nContB)
    implicit none
    integer,intent(in) :: nPrimA,nPrimB,ndim,nContA,nContB,nContQ
    real(reals),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(reals),intent(in) :: ERECS(nContA,nContB,nPrimA,nPrimB,ndim)
    real(reals),intent(inout) :: CERECS(nContA,nContB,ndim)
    !
    integer :: i,icA,icB,ipA,ipB,iQ
    real(reals) :: TMP,TMPB,TMP2
    !$OMP DO COLLAPSE(3) PRIVATE(i,icA,icB,ipA,ipB,iQ,TMP,TMPB,TMP2)
    do i = 1,ndim
       do icB = 1,nContB
          do icA = 1,nContA
             TMP2 = 0.0E0_reals
             do ipB = 1,nPrimB
                TMPB = BCC(ipB,icB)
                do ipA = 1,nPrimA
                   TMP = ACC(ipA,icA)*TMPB
                   TMP2 = TMP2 + ERECS(icA,icB,ipA,ipB,i)*TMP                   
                enddo
             enddo
             CERECS(icA,icB,i) = TMP2 
          enddo
       enddo
    enddo
    !$OMP END DO

  end subroutine SPGabcontractBasisGenP

  subroutine SPGabcontractBasisSegP(ERECS,CERECS,nContQ,ndim,nPrimA,nPrimB)
    implicit none
    integer,intent(in) :: nPrimA,nPrimB,ndim,nContQ
    real(reals),intent(in) :: ERECS(nPrimA,nPrimB,ndim)
    real(reals),intent(inout) :: CERECS(ndim)
    !
    integer :: i,ipA,ipB,iQ
    !$OMP DO PRIVATE(i,ipA,ipB,iQ)
    do i = 1,ndim
       CERECS(i) = 0.0E0_reals
       do ipB = 1,nPrimB
          do ipA = 1,nPrimA
             CERECS(i) = CERECS(i) + ERECS(ipA,ipB,i)
          enddo
       enddo
    enddo
    !$OMP END DO
  end subroutine SPGabcontractBasisSegP

  subroutine SPGabSphericalTransformGenPAB(CERECS,SCERECS,SPHMATA,SPHMATB,&
       & ijk1,ijk1s,ijk2,ijk2s,ndim,nOrbCompQ)
    implicit none
    integer,intent(in) :: ijk1,ijk1s,ijk2,ijk2s,ndim,nOrbCompQ
    real(reals),intent(in) :: SPHMATA(ijk1,ijk1s),SPHMATB(ijk2,ijk2s)
    real(reals),intent(in) :: CERECS(ndim,ijk1,ijk2,ijk1s,ijk2s)
    real(reals),intent(inout) :: SCERECS(ndim,ijk1s,ijk2s)
    !
    integer :: a,b,as,bs,i,ijkQ
    real(reals) :: TMP,TMPB
    !$OMP DO COLLAPSE(2) PRIVATE(a,b,as,bs,i,TMP,TMPB)
     do as = 1,ijk1s
      do bs = 1,ijk2s
       do i = 1,ndim
        SCERECS(i,as,bs) = 0.0E0_reals
       enddo 
       do b = 1,ijk2
        TMPB = SPHMATB(b,bs)
        do a = 1,ijk1
         TMP = SPHMATA(a,as)*TMPB
         do i = 1,ndim
          SCERECS(i,as,bs) = SCERECS(i,as,bs) + CERECS(i,a,b,as,bs)*TMP
         enddo
        enddo
       enddo
      enddo
     enddo
    !$OMP END DO
  end subroutine SPGabSphericalTransformGenPAB

  subroutine SPGabSphericalTransformGenPA(CERECS,SCERECS,SPHMATA,&
       & ijk1,ijk1s,ijk2,ndim,nOrbCompQ)
    implicit none
    integer,intent(in) :: ijk1,ijk1s,ijk2,ndim,nOrbCompQ
    real(reals),intent(in) :: SPHMATA(ijk1,ijk1s)
    real(reals),intent(in) :: CERECS(ndim,ijk1,ijk2,ijk1s,ijk2)
    real(reals),intent(inout) :: SCERECS(ndim,ijk1s,ijk2)
    !
    integer :: a,as,i,b
    real(reals) :: TMP
    !$OMP DO COLLAPSE(2) PRIVATE(a,as,i,b,TMP)
    do b = 1,ijk2
       do as = 1,ijk1s
          do i = 1,ndim
             SCERECS(i,as,b) = 0.0E0_reals
          enddo
          do a = 1,ijk1
             TMP = SPHMATA(a,as)
             do i = 1,ndim
                SCERECS(i,as,b) = SCERECS(i,as,b) + CERECS(i,a,b,as,b)*TMP
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO
  end subroutine SPGabSphericalTransformGenPA

  subroutine SPGabSphericalTransformGenPB(CERECS,SCERECS,SPHMATB,&
       & ijk1,ijk2,ijk2s,ndim,nOrbCompQ)
    implicit none
    integer,intent(in) :: ijk1,ijk2,ijk2s,ndim,nOrbCompQ
    real(reals),intent(in) :: SPHMATB(ijk2,ijk2s)
    real(reals),intent(in) :: CERECS(ndim,ijk1,ijk2,ijk1,ijk2s)
    real(reals),intent(inout) :: SCERECS(ndim,ijk1,ijk2s)
    !
    integer :: b,bs,i,a
    real(reals) :: TMP
    !$OMP DO COLLAPSE(2) PRIVATE(b,bs,i,a,TMP)
    do bs = 1,ijk2s
       do a = 1,ijk1
          do i = 1,ndim
             SCERECS(i,a,bs) = 0.0E0_reals
          enddo
          do b = 1,ijk2
             TMP = SPHMATB(b,bs)
             do i = 1,ndim
                SCERECS(i,a,bs) = SCERECS(i,a,bs) + CERECS(i,a,b,a,bs)*TMP
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO
  end subroutine SPGabSphericalTransformGenPB

end module SPIchorEriGabintegralCPUMcMGeneralMod
