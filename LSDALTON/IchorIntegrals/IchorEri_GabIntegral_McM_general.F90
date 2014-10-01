!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriGabintegralCPUMcMGeneralMod
use IchorEriCoulombintegralCPUMcMGeneralMod
use IchorPrecisionModule
use IchorMemory
use IchorCommonModule
use IchorEriCoulombintegralCPUMcMGeneralEcoeffMod, only: &
     & Ichorbuild_Ecoeff_RHS,Ichorbuild_Ecoeff_LHS, printEcoeff
use IchorEriCoulombintegralCPUMcMGeneralWTUVMod

CONTAINS
  subroutine IGI_CPU_McM_general(nPrimA,nPrimB,&
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
    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)
    logical,intent(in)     :: Psegmented,sphericalGTO,PQorder
    real(realk),intent(in) :: Acenter(3),Bcenter(3)
    real(realk),intent(in) :: pexp(nPrimP),pcent(3*nPrimP),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(inout) :: CDAB(1)    
    real(realk),intent(in) :: integralPrefactor(nPrimP*nPrimP)!integralPrefactor(nPrimQ,nPrimP)    
    real(realk),intent(in) :: reducedExponents(nPrimP*nPrimP) !reducedExponents(nPrimQ,nPrimP)
    real(realk),intent(in) :: Pdistance12(3) 
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
    call build_Gab_Rpq(nPrimP,Pcent,TmpArray3)
    !
    !      builds RJ000(0:AngmomPP,nPrimP,nPrimP) Store in TmpArray4
    !
#ifdef VAR_DEBUGICHOR
    IF(nTmpArray4.LT.nPrimP*nPrimP*(AngmomPP+1))call ichorquit('IchorTmp1G1',-1)
#endif
    call buildRJ000_general(nPasses,nPrimP,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&
         & TABFJW,TmpArray4,AngmomPP,integralPrefactor,TmpArray3)
    IF (INTPRINT .GE. 10) CALL PrintRJ000(TmpArray4,AngmomPP,nPrimP*nPrimP,nPasses,lupri)    
    !
    !     Build WTUV(nPrimP,nPrimP,nPasses,nTUV) RJ000 = TmpArray4, Rpq = TmpArray3
    !
    TMP1 = .TRUE.
    IF (AngmomPP.EQ. 0) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3A',-1)
#endif
       call IchorwtuvRecurrenceJMIN0JMAX0(TmpArray4,TmpArray1,nPrimP*nPrimP)
    ELSEIF (AngmomPP.EQ. 1) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3B',-1)
#endif
       call IchorwtuvRecurrenceJMIN0JMAX1(TmpArray4,TmpArray1,TmpArray3,nPrimP*nPrimP)
    ELSEIF (AngmomPP.EQ. 2) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3C',-1)
#endif
       call IchorwtuvRecurrenceJMIN0JMAX2(TmpArray4,TmpArray1,TmpArray3,nPrimP*nPrimP)
    ELSEIF (AngmomPP.EQ. 3) THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3D',-1)
#endif
       call IchorwtuvRecurrenceJMIN0JMAX3(TmpArray4,TmpArray1,TmpArray3,nPrimP*nPrimP)
    ELSE !AngmomPP > 3
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3a',-1)
       IF(TMParray2maxsize.LT.nPrimP*nPrimP*ntuv)call ichorquit('IchorTmp1G3b',-1)
#endif
       J = AngmomPP-3
       call IchorwtuvRecurrenceJMIN0JMAX3J(TmpArray4,TmpArray1,TmpArray3,nPrimP*nPrimP,AngmomPP,J)
       TMP1 = .TRUE.
       !J minimum goes from 0 to 0. For  AngmomPP=6 J takes the values 2,1,0 
       DO j=AngmomPP-4,0,-1 
          IF(TMP1)THEN
             call IchorwtuvRecurrenceCurrent(TmpArray1,TmpArray2,J,AngmomPP,nPrimP,nPrimP,nPasses,&
                  & ntuv,TmpArray3,TmpArray4)
          ELSE
             call IchorwtuvRecurrenceCurrent(TmpArray2,TmpArray1,J,AngmomPP,nPrimP,nPrimP,nPasses,&
                  & ntuv,TmpArray3,TmpArray4)
          ENDIF
          TMP1 = .NOT.TMP1
       ENDDO
    ENDIF    
    IF (IntPrint .GE. 25)THEN
       IF(TMP1)THEN
          call PrintWTUV(TmpArray1,AngmomPP,nPrimP*nPrimP,nPasses,nTUV,lupri)
       ELSE
          call PrintWTUV(TmpArray2,AngmomPP,nPrimP*nPrimP,nPasses,nTUV,lupri)
       ENDIF
    ENDIF

    !builds Ecoeff(nPrimP,nPasses,nTUVP,nCartOrbCompP)
    !FOR NOW THIS IS NOT OpenMP parallized - unclear what the best method is. 
    call Ichorbuild_Ecoeff_RHS(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
         & nCartOrbCompP,Aexp,Bexp,TmpArray3,Pdistance12,Ppreexpfac,intprint,lupri,TmpArray4)

    IF (IntPrint .GE. 25)call printEcoeff(TmpArray3,nTUVP,nCartOrbCompP,nPrimP,nPassQ,lupri)

    IF(TMP1)THEN !current intermediate WTUV reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nPrimP*nPrimP*ntuvP*nCartOrbCompP)call ichorquit('IchorTmp1G3P1',-1)
#endif
       !builds RE(nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP)
       call DirectcontractEQgen(TmpArray1,TmpArray2,nPrimP*nPrimP,nPrimP,nPrimP,nTUV,ntuvP,ntuvP,&
            & TmpArray3,nCartOrbCompP,AngmomA,AngmomB,AngmomA,AngmomB,AngmomPP)
       IF (IntPrint .GE. 25)call PrintIchorTensorRE(TmpArray2,nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP,lupri)
    ELSE !current intermediate WTUV reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nPrimP*nPrimP*ntuvP*nCartOrbCompP)call ichorquit('IchorTmp1G3Q2',-1)
#endif
       !builds RE(nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP)
       call DirectcontractEQgen(TmpArray2,TmpArray1,nPrimP*nPrimP,nPrimP,nPrimP,nTUV,ntuvP,ntuvP,&
            & TmpArray3,nCartOrbCompP,AngmomA,AngmomB,AngmomA,AngmomB,AngmomPP)
       IF (IntPrint .GE. 25)call PrintIchorTensorRE(TmpArray1,nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP,lupri)
    ENDIF
    TMP1 = .NOT.TMP1
    
    IF(TMP1)THEN !current intermediate RE(nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP) reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContP*nPrimP*nTUVP*nCartOrbCompP)call ichorquit('IchorTmp2G4',-1)
#endif
       !Build REC(nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP)
       IF(Psegmented)THEN
          call contractBasisSegQ(TmpArray1,TmpArray2,nPrimA,nPrimB,nPrimP*nTUVP*nCartOrbCompP) 
       ELSE
          call contractBasisGenQ(TmpArray1,TmpArray2,ACC,BCC,nPrimA,nPrimB,nContA,nContB,nPrimP*nTUVP*nCartOrbCompP) 
       ENDIF
       IF (IntPrint .GE. 25)call PrintIchorTensorREC(TmpArray2,nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP,lupri)
    ELSE !current intermediate RE(nPrimP,nPrimP,nPasses,nTUVP,nCartOrbCompP) reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContP*nPrimP*nTUVP*nCartOrbCompP)call ichorquit('IchorTmp2G4',-1)
#endif
       !Build REC(nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP)
       IF(Psegmented)THEN
          call contractBasisSegQ(TmpArray2,TmpArray1,nPrimA,nPrimB,nPrimP*nTUVP*nCartOrbCompP) 
       ELSE
          call contractBasisGenQ(TmpArray2,TmpArray1,ACC,BCC,nPrimA,nPrimB,nContA,nContB,nPrimP*nTUVP*nCartOrbCompP) 
       ENDIF
       IF (IntPrint .GE. 25)call PrintIchorTensorREC(TmpArray1,nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP,lupri)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransP)THEN
       IF(TMP1)THEN !current intermediate REC(nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP) reside in TmpArray1
#ifdef VAR_DEBUGICHOR
          IF(TMParray2maxsize.LT.nContP*nPrimP*nTUVP*nOrbCompP)call ichorquit('IchorTmp1G5',-1)
#endif
          !Build RECS(nContP,nPrimP,nPasses,nTUVP,nOrbCompP)
          IF(Sph1.AND.Sph2)THEN
             call SphericalTransformGenQCD(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,SPH_MAT(AngmomB)%elms,&
                  & nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nPrimP*nTUVP)
          ELSEIF(Sph1)THEN
             call SphericalTransformGenQC(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContP*nPrimP*nTUVP)
          ELSE
             call SphericalTransformGenQD(TmpArray1,TmpArray2,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nPrimP*nTUVP)
          ENDIF
          IF (IntPrint .GE. 25) call PrintIchorTensorRECS(TmpArray2,nContP,nPrimP,nPasses,nTUVP,nOrbCompP,lupri)
       ELSE !current intermediate REC(nContP,nPrimP,nPasses,nTUVP,nCartOrbCompP) reside in TmpArray2
#ifdef VAR_DEBUGICHOR
          IF(TMParray1maxsize.LT.nContP*nPrimP*nTUVP*nOrbCompP)call ichorquit('IchorTmp1G5',-1)
#endif
          !Build RECS(nContP,nPrimP,nPasses,nTUVP,nOrbCompP)
          IF(Sph1.AND.Sph2)THEN
             call SphericalTransformGenQCD(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,SPH_MAT(AngmomB)%elms,&
                  & nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nPrimP*nTUVP)
          ELSEIF(Sph1)THEN
             call SphericalTransformGenQC(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContP*nPrimP*nTUVP)
          ELSE
             call SphericalTransformGenQD(TmpArray2,TmpArray1,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nPrimP*nTUVP)
          ENDIF
          IF (IntPrint .GE. 25) call PrintIchorTensorRECS(TmpArray1,nContP,nPrimP,nPasses,nTUVP,nOrbCompP,lupri)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

    !builds Ecoeff(nPrimP,nPasses,nTUVP,nCartOrbCompP)
    !currently not OpenMP parallel 
    call Ichorbuild_Ecoeff_LHS(nPrimP,nPrimA,nPrimB,AngmomP,AngmomA,AngmomB,nTUVP,&
         & nCartOrbCompP,Aexp,Bexp,TmpArray3,Pdistance12,Ppreexpfac,nPasses,&
         & 1,1,IatomApass,IatomBpass,1,intprint,lupri,TmpArray4)

    IF (IntPrint .GE. 25)call printEcoeff(TmpArray3,nTUVP,nCartOrbCompP,nPrimP,nPasses,lupri)
    
    !builds ERECS(nContP,nPrimP,nPasses,nCartOrbCompP,nOrbCompP)
    IF(TMP1)THEN
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContP*nPrimP*nCartOrbCompP*nOrbCompP)call ichorquit('IchorTmp2G6A',-1)
#endif
       call contractEcoeffGenP(TmpArray1,TmpArray2,TmpArray3,nOrbCompP,nCartOrbCompP,nTUVP,nContP,nPrimP)

       IF (IntPrint .GE. 25)call PrintIchorTensorERECS(TmpArray2,nContP,nPrimP,nPasses,nCartOrbCompP,nOrbCompP,lupri)
    ELSE
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContP*nPrimP*nCartOrbCompP*nOrbCompP)call ichorquit('IchorTmp1G6B',-1)
#endif
       call contractEcoeffGenP(TmpArray2,TmpArray1,TmpArray3,nOrbCompP,nCartOrbCompP,nTUVP,nContP,nPrimP)
       IF (IntPrint .GE. 25)call PrintIchorTensorERECS(TmpArray1,nContP,nPrimP,nPasses,nCartOrbCompP,nOrbCompP,lupri)
    ENDIF
    TMP1 = .NOT.TMP1

    !note we build CERECS(nContP,nCartOrbCompP,nOrbCompP)
    !normally we build CERECS(nContP,nContP,nPasses,nCartOrbCompP,nOrbCompP)
    !but we are only intrested in the diagonal elements iContP = iContP



    IF(.NOT.TMP1)THEN !current intermediate reside in TmpArray2
#ifdef VAR_DEBUGICHOR
       IF(TMParray1maxsize.LT.nContP*nContP*nCartOrbCompP*nOrbCompP)call ichorquit('IchorTmp1G7',-1)
#endif
       !builds CERECS(nContP,nContP,nPasses,nCartOrbCompP,nOrbCompP)
       IF(Psegmented)THEN
          call contractBasisSegP(TmpArray2,TmpArray1,nContP,nCartOrbCompP*nOrbCompP,nPrimA,nPrimB) 
       ELSE
          call contractBasisGenP(TmpArray2,TmpArray1,ACC,BCC,nContP,nCartOrbCompP*nOrbCompP,nPrimA,nPrimB,nContA,nContB) 
       ENDIF
       IF (IntPrint .GE. 25) call PrintIchorTensorCERECS(TmpArray1,nContP,nContP,nPasses,nCartOrbCompP,nOrbCompP,lupri)
    ELSE !current intermediate reside in TmpArray1
#ifdef VAR_DEBUGICHOR
       IF(TMParray2maxsize.LT.nContP*nContP*nPasses*nCartOrbCompP*nOrbCompP)call ichorquit('IchorTmp2G7',-1)
#endif
       !builds CERECS(nContP,nContP,nPasses,nCartOrbCompP,nOrbCompP)
       IF(Psegmented)THEN
          call contractBasisSegP(TmpArray1,TmpArray2,nContP,nCartOrbCompP*nOrbCompP,nPrimA,nPrimB) 
       ELSE
          call contractBasisGenP(TmpArray1,TmpArray2,ACC,BCC,nContP,nCartOrbCompP*nOrbCompP,nPrimA,nPrimB,nContA,nContB) 
       ENDIF
       IF (IntPrint .GE. 25) call PrintIchorTensorCERECS(TmpArray2,nContP,nContP,nPasses,nCartOrbCompP,nOrbCompP,lupri)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransP)THEN
       IF(TMP1)THEN !current intermediate reside in TmpArray1
#ifdef VAR_DEBUGICHOR
          IF(TMParray2maxsize.LT.nContP*nContP*nOrbCompP*nOrbCompP)call ichorquit('IchorTmp1G8',-1)
#endif
          !builds SCEREC(nContP,nContP,nPasses,nOrbCompP,nOrbCompP)
          IF(Sph1.AND.Sph2)THEN
             call SphericalTransformGenPAB(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,&
                  & SPH_MAT(AngmomB)%elms,nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,&
                  & nContP*nContP,nOrbCompP)
          ELSEIF(Sph1)THEN
             call SphericalTransformGenPA(TmpArray1,TmpArray2,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContP*nContP,nOrbCompP)
          ELSE
             call SphericalTransformGenPB(TmpArray1,TmpArray2,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nContP,nOrbCompP)
          ENDIF
          IF(IntPrint.GE.25)call PrintIchorTensorSCERECS(TmpArray2,&
               & nContP,nContP,nPasses,nOrbCompP,nOrbCompP,lupri)
       ELSE !current intermediate reside in TmpArray2
#ifdef VAR_DEBUGICHOR
          IF(TMParray1maxsize.LT.nContP*nContP*nOrbCompP*nOrbCompP)call ichorquit('IchorTmp2G8',-1)
#endif
          !builds SCEREC(nContP,nContP,nPasses,nOrbCompP,nOrbCompP)
          IF(Sph1.AND.Sph2)THEN
             call SphericalTransformGenPAB(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,&
                  & SPH_MAT(AngmomB)%elms,nCartOrbCompA,nOrbCompA,nCartOrbCompB,nOrbCompB,&
                  & nContP*nContP,nOrbCompP)
          ELSEIF(Sph1)THEN
             call SphericalTransformGenPA(TmpArray2,TmpArray1,SPH_MAT(AngmomA)%elms,&
                  & nCartOrbCompA,nOrbCompA,nOrbCompB,nContP*nContP,nOrbCompP)
          ELSE
             call SphericalTransformGenPB(TmpArray2,TmpArray1,SPH_MAT(AngmomB)%elms,&
                  & nOrbCompA,nCartOrbCompB,nOrbCompB,nContP*nContP,nOrbCompP)
          ENDIF          
          IF(IntPrint.GE.25)call PrintIchorTensorSCERECS(TmpArray1,&
               & nContP,nContP,nPasses,nOrbCompP,nOrbCompP,lupri)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

    IF(.NOT.TMP1)THEN
       call extractGabElmGen(TmpArray2,CDAB,nContP,nOrbCompP)
    ELSE
       call extractGabElmGen(TmpArray1,CDAB,nContP,nOrbCompP)       
    ENDIF
    
  end subroutine IGI_CPU_McM_general

  subroutine IGI_CPU_McM_general_size(TMParray1maxsize,&
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
    IF(TMP1)THEN 
       TMParray2maxsize = MAX(TMParray2maxsize,nPrimP*nPrimP*ntuvP*nCartOrbCompQ)
    ELSE 
       TMParray1maxsize = MAX(TMParray1maxsize,nPrimP*nPrimP*ntuvP*nCartOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1
    IF(TMP1)THEN 
       TMParray2maxsize = MAX(TMParray2maxsize,nContP*nPrimP*nTUVP*nCartOrbCompQ)
    ELSE 
       TMParray1maxsize = MAX(TMParray1maxsize,nContP*nPrimP*nTUVP*nCartOrbCompQ)
    ENDIF
    TMP1 = .NOT.TMP1

    IF(SphericalTransQ)THEN
       IF(TMP1)THEN 
          TMParray2maxsize = MAX(TMParray2maxsize,nContP*nPrimP*nTUVP*nOrbCompQ)
       ELSE 
          TMParray1maxsize = MAX(TMParray1maxsize,nContP*nPrimP*nTUVP*nOrbCompQ)
       ENDIF
       TMP1 = .NOT.TMP1
    ENDIF

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
  end subroutine IGI_CPU_McM_general_size
  
  subroutine build_Gab_Rpq(nPrimP,Pcent,Rpq)
    implicit none
    integer,intent(in) :: nPrimP
    real(realk),intent(in) :: Pcent(3,nPrimP)
    real(realk),intent(inout) :: Rpq(nPrimP,nPrimP,3)
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
  end subroutine build_Gab_Rpq
  
  subroutine extractGabElmGen(TmpArray2,CDAB,nContP,nOrbCompP)
    implicit none
    integer,intent(in) :: nContP,nOrbCompP
    real(realk),intent(in) :: TmpArray2(nContP,nContP,nOrbCompP,nOrbCompP)
    real(realk),intent(inout) :: CDAB(1)
    !
    integer :: iContP,iOrbP
    real(realk) :: TMP
    !$OMP MASTER
    TMP = ABS(TmpArray2(1,1,1,1))
    do iOrbP = 1,nOrbCompP
       do iContP = 1,nContP
          TMP = MAX(TMP,ABS(TmpArray2(iContP,iContP,iOrbP,iOrbP)))
       enddo
    enddo
    CDAB(1) = SQRT(TMP)
    !$OMP END MASTER
    !$OMP BARRIER 

  end subroutine extractGabElmGen


end MODULE IchorEriGabintegralCPUMcMGeneralMod
