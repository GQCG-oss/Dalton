MODULE IntegralInterfaceDEC
  use precision
  use TYPEDEFTYPE, only: LSSETTING, LSINTSCHEME, LSITEM, integralconfig
  use Matrix_module, only: MATRIX, MATRIXP
  use TYPEDEF, only: getNbasis, retrieve_output, gcao2ao_transform_matrixd2, &
       & retrieve_screen_output, ao2gcao_transform_matrixf
  use ls_Integral_Interface, only: ls_setDefaultFragments, ls_getIntegrals1, setaobatch
  use integraloutput_typetype, only: INTEGRALOUTPUT
  use integraloutput_type, only: initintegraloutputdims
  use lstensor_operationsmod, only: lstensor, lstensor_nullify,build_batchgab,&
       & build_batchgabk
  use screen_mod, only: DECscreenItem, nullify_decscreen, &
       & null_decscreen_and_associate_mastergab_rhs, &
       & null_decscreen_and_associate_mastergab_lhs, &
       & init_decscreen_batch
  use LSparameters
  use f12_module, only: set_ggem,stgfit
  use ao_typetype, only: aoitem, BATCHORBITALINFO
  use ao_type, only: free_aoitem, freebatchorbitalinfo, initbatchorbitalinfo, &
       & setbatchorbitalinfo
  use BUILDAOBATCH, only: BUILD_SHELLBATCH_AO,determinebatchindexandsize
  use lstiming
  use memory_handling,only: mem_alloc, mem_dealloc
  use basis_typetype, only: BASISSETINFO,RegBasParam,CABBasParam
  use IchorErimoduleHost
  PUBLIC:: II_precalc_DECScreenMat,II_getBatchOrbitalScreen,&
       & II_getBatchOrbitalScreen2,II_getBatchOrbitalScreenK,&
       & II_getBatchOrbitalScreen2K,II_GET_DECPACKED4CENTER_J_ERI,&
       & II_GET_DECPACKED4CENTER_K_ERI, II_GET_DECBATCHPACKED,&
       & II_get_eri_integralblock, II_get_eri_integralblock_inquire,&
       & II_GET_DECPACKED4CENTER_J_ERI2
  PRIVATE
CONTAINS

!> \brief Calculates and stores the screening integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
SUBROUTINE II_precalc_DECScreenMat(DECscreen,LUPRI,LUERR,SETTING,ndimA,ndimG,intspec,intThreshold)
IMPLICIT NONE
TYPE(DECscreenITEM),intent(INOUT)    :: DecScreen
TYPE(LSSETTING)      :: SETTING
INTEGER,intent(in)   :: LUPRI,LUERR,ndimA,ndimG
Character,intent(IN)          :: intSpec(5)
REAL(REALK),intent(IN)        :: intTHRESHOLD
!
TYPE(lstensor),pointer :: GAB
LOGICAL :: IntegralTransformGC,CSintsave,PSintsave
INTEGER :: I,J,nbast,natoms,Oper
integer             :: IJ,ao(4),dummy
logical,pointer     :: OLDsameBAS(:,:)
real(realk)         :: coeff(6),exponent(6),tmp
real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
integer             :: nGaussian,nG2
real(realk)         :: GGem,slater
dummy=1
IF(SETTING%SCHEME%CS_SCREEN.OR.SETTING%SCHEME%PS_SCREEN)THEN
   !set geminal
   IF (intSpec(5).NE.'C'.AND.intSpec(5).NE.'E') THEN
      nGaussian = 6
      nG2 = nGaussian*(nGaussian+1)/2
      GGem = 0E0_realk
      slater = setting%basis(1)%p%BINFO(RegBasParam)%GeminalScalingFactor
      call stgfit(slater,nGaussian,exponent,coeff)
      IJ=0
      DO I=1,nGaussian
         DO J=1,I
            IJ = IJ + 1
            coeff2(IJ) = 2E0_realk * coeff(I) * coeff(J)
            prodexponent(IJ) = exponent(I) * exponent(J)
            sumexponent(IJ) = exponent(I) + exponent(J)
         ENDDO
         coeff2(IJ) = 0.5E0_realk*coeff2(IJ)
      ENDDO
   ENDIF

   ! ***** SELECT OPERATOR TYPE *****
   IF (intSpec(5).EQ.'C') THEN
      ! Regular Coulomb operator 1/r12
      oper = CoulombOperator
   ELSE IF (intSpec(5).EQ.'E') THEN
      ! Long-Range operator erf/r12
      oper = ErfOperator
   ELSE IF (intSpec(5).EQ.'G') THEN
     ! The Gaussian geminal operator g
      oper = GGemOperator
      call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
   ELSE IF (intSpec(5).EQ.'F') THEN
      ! The Gaussian geminal divided by the Coulomb operator g/r12
      oper = GGemCouOperator
      call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
   ELSE IF (intSpec(5).EQ.'D') THEN
      ! The double commutator [[T,g],g]
      oper = GGemGrdOperator
      call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
   ELSE IF (intSpec(5).EQ.'2') THEN
      ! The Gaussian geminal operator squared g^2
      oper = GGemOperator
      call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
   ELSE
      call lsquit('Error in specification of operator in InitGaussianGeminal',-1)
   ENDIF

   
   ! ***** SELECT AO TYPES *****
   DO i=1,4
      IF (intSpec(i).EQ.'R') THEN
         !   The regular AO-basis
         ao(i) = AORegular
      ELSE IF (intSpec(i).EQ.'C') THEN
         !   The CABS AO-type basis
         ao(i) = AOdfCABS
      ELSE
         call lsquit('Error in specification of ao1 in II_precalc_DECScreenMat',-1)
      ENDIF
   ENDDO
   
   call mem_alloc(OLDsameBAS,setting%nAO,setting%nAO)
   OLDsameBAS = setting%sameBAS
   DO I=1,4
      DO J=1,4
         setting%sameBas(I,J) = setting%sameBas(I,J) .AND. ao(I).EQ.ao(J)
      ENDDO
   ENDDO
   
   SETTING%SCHEME%intTHRESHOLD=intTHRESHOLD*1.0E-5_realk

   call null_decscreen_and_associate_MasterGab_RHS(GAB,DecSCREEN)
   IntegralTransformGC=.FALSE.
   CSintsave = setting%scheme%CS_int
   PSintsave = setting%scheme%PS_int
   setting%scheme%CS_int = SETTING%SCHEME%CS_SCREEN
   setting%scheme%PS_int = SETTING%SCHEME%PS_SCREEN
   !BYPASS MPI in ls_getIntegrals because MPI is done on dec level.   
   call initIntegralOutputDims(setting%Output,dummy,dummy,1,1,1)
   CALL ls_getIntegrals1(AO(3),AO(4),AO(3),AO(4),&
        &Oper,RegularSpec,ContractedInttype,0,SETTING,LUPRI,LUERR)
   CALL retrieve_screen_Output(lupri,setting,GAB,IntegralTransformGC)   
   CALL init_DECscreen_batch(ndimA,ndimG,DecSCREEN)
   !LHS
   call null_decscreen_and_associate_MasterGab_LHS(GAB,DecSCREEN)
   !BYPASS MPI in ls_getIntegrals because MPI is done on dec level.
   call initIntegralOutputDims(setting%Output,dummy,dummy,1,1,1)
   CALL ls_getIntegrals1(AO(1),AO(2),AO(1),AO(2),&
        &Oper,RegularSpec,ContractedInttype,0,SETTING,LUPRI,LUERR)
   CALL retrieve_screen_Output(lupri,setting,GAB,IntegralTransformGC)   

   setting%scheme%CS_int = CSintsave
   setting%scheme%PS_int = PSintsave

   setting%sameBas = OLDsameBAS
   call mem_dealloc(OLDsameBAS) 
ELSE
   call nullify_decscreen(DecSCREEN)
ENDIF
END SUBROUTINE II_precalc_DECScreenMat

!> \brief build BatchOrbitalInfo
!> \author T. Kjaergaard
!> \date 2010
!> \param setting Integral evalualtion settings
!> \param AO options AORdefault,AODFdefault or AOempty
!> \param intType options contracted or primitive
!> \param nbast the number of basis functions
!> \param orbtobast for a given basis function it provides the batch index
!> \param nBatches the number of batches
!> \param lupri Default print unit
!> \param luerr Default error print unit
SUBROUTINE II_getBatchOrbitalScreen(DecScreen,Setting,nBast,nbatches1,nbatches2,&
     & batchsize1,batchsize2,batchindex1,batchindex2,batchdim1,batchdim2,&
     & INTSPEC,lupri,luerr)
implicit none
TYPE(DECscreenITEM),intent(INOUT)    :: DecScreen
TYPE(LSSETTING) :: Setting !,intent(INOUT)
Integer,intent(IN)         :: nBast,lupri,luerr
Integer,intent(IN)         :: nBatches1,nBatches2
Integer,intent(IN)         :: batchdim1(nBatches1),batchdim2(nBatches2)
Integer,intent(IN)         :: batchsize1(nBatches1),batchsize2(nBatches2)
Integer,intent(IN)         :: batchindex1(nBatches1),batchindex2(nBatches2)
Character,intent(IN)       :: intSpec(5)
!
TYPE(AOITEM)           :: AObuild3,AObuild4
Integer                :: nDim(4),AO(4),i

DO i=1,4
   IF (intSpec(i).EQ.'R') THEN
      !   The regular AO-basis
      ao(i) = AORegular
   ELSE IF (intSpec(i).EQ.'C') THEN
      !   The CABS AO-type basis
      ao(i) = AOdfCABS
   ELSE
      call lsquit('Error in specification of ao1 in II_precalc_DECScreenMat',-1)
   ENDIF
ENDDO
CALL setAObatch(AObuild3,0,1,nDim(3),AO(3),Contractedinttype,Setting%scheme,&
     & Setting%fragment(1)%p,setting%basis(1)%p,lupri,luerr)
CALL setAObatch(AObuild4,0,1,nDim(4),AO(4),Contractedinttype,Setting%scheme,&
     & Setting%fragment(2)%p,setting%basis(2)%p,lupri,luerr)
CALL II_getBatchOrbitalScreen2(DECscreen,setting,AObuild3,AObuild4,nBatches1,&
     & nBatches2,batchsize1,batchsize2,batchindex1,batchindex2,batchdim1,&
     & batchdim2,AO,lupri,luerr)
CALL free_aoitem(lupri,AObuild3)
CALL free_aoitem(lupri,AObuild4)
END SUBROUTINE II_getBatchOrbitalScreen

subroutine II_getBatchOrbitalScreen2(DECscreen,setting,AOfull1,AOfull2,&
     & nBatches1,nBatches2,batchsize1,batchsize2,batchindex1,batchindex2,&
     & batchdim1,batchdim2,AO,lupri,luerr)
implicit none
TYPE(DECscreenITEM),intent(INOUT):: DecScreen
TYPE(LSSETTING),intent(INOUT)    :: Setting
TYPE(AOITEM),intent(IN)          :: AOfull1,AOfull2
Integer,intent(IN)               :: lupri,nBatches1,nBatches2,luerr
Integer,intent(IN)               :: batchdim1(nBatches1),batchdim2(nBatches2)
Integer,intent(IN)               :: batchsize1(nBatches1),batchsize2(nBatches2)
Integer,intent(IN)               :: batchindex1(nBatches1),batchindex2(nBatches2)
Integer,intent(in)               :: AO(4)
!
TYPE(AOITEM)            :: AObatch1,AObatch2
integer :: dim1,dim2,jBatch,iBatch,I,AObatchdim1,AObatchdim2
integer :: AOT1batch2,AOT1batch1,jBatch1,iBatch1,jbat,ibat
character(len=9) :: Jlabel,Ilabel
logical :: FoundInMem,FoundOnDisk,uncont,intnrm
type(lstensor),pointer :: MasterGAB
type(lstensor),pointer :: GAB
TYPE(BASISSETINFO),pointer :: AObasis1,AObasis2
IF(ASSOCIATED(DECSCREEN%masterGabRHS))THEN
   MasterGAB => DECSCREEN%masterGabRHS
   uncont = .FALSE.
   intnrm = .FALSE.
   IF(.NOT.ASSOCIATED(DECSCREEN%batchGab))CALL LSQUIT('DECSCREEN%batchGab not associated',-1)
   AOT1batch2 = 0

   SELECT CASE(AO(3))
   CASE (AORegular)
      AObasis1 => setting%basis(3)%p%BINFO(RegBasParam)
   CASE (AOdfCABS)
      AObasis1 => setting%basis(3)%p%BINFO(CABBasParam)
   CASE DEFAULT
      print*,'case: ',AO(3)
      WRITE(luerr,*) 'case: ',AO(3)
      CALL LSQuit('Programming error: Not a case in II_getBatchOrbitalScreen2A!',lupri)
   END SELECT
   
   SELECT CASE(AO(4))
   CASE (AORegular)
      AObasis2 => setting%basis(4)%p%BINFO(RegBasParam)
   CASE (AOdfCABS)
      AObasis2 => setting%basis(4)%p%BINFO(CABBasParam)
   CASE DEFAULT
      print*,'case: ',AO(4)
      WRITE(luerr,*) 'case: ',AO(4)
      CALL LSQuit('Programming error: Not a case in II_getBatchOrbitalScreen2B!',lupri)
   END SELECT
   
   DO jBatch1=1,nBatches2
      dim2 = batchdim2(jBatch1)
      jbat = batchindex2(jBatch1)

      call BUILD_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
           & Setting%molecule(4)%p,aobasis2,AObatch2,&
           & uncont,intnrm,jBat,AObatchdim2,batchsize2(jBatch1))
      IF(AObatchdim2.NE.dim2)CALL LSQUIT(' typedef_getBatchOrbitalScreen2 mismatch A',-1)
      IF(batchsize2(jBatch1).NE.AObatch2%nbatches)CALL LSQUIT(' typedef_getBatchOrbitalScreen2 mismatch B',-1)

      
      AOT1batch1 = 0
      DO iBatch1=1,nBatches1
         dim1 = batchdim1(iBatch1)      
         ibat = batchindex1(iBatch1)

         call BUILD_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
              & Setting%molecule(3)%p,aobasis1,AObatch1,&
              & uncont,intnrm,iBat,AObatchdim1,batchsize1(iBatch1))
         IF(AObatchdim1.NE.dim1)CALL LSQUIT(' typedef_getBatchOrbitalScreen mismatch A2',-1)
         IF(batchsize1(iBatch1).NE.AObatch1%nbatches)CALL LSQUIT(' typedef_getBatchOrbitalScreen mismatch B2',-1)
         
         nullify(DECSCREEN%batchGab(iBatch1,jBatch1)%p)         
         allocate(DECSCREEN%batchGab(iBatch1,jBatch1)%p)         
         call build_BatchGab(AOfull1,AOfull2,AObatch1,AObatch2,iBatch1,jBatch1,AOT1batch1,AOT1batch2,dim1,dim2,&
              & MasterGAB,DECSCREEN%batchGab(iBatch1,jBatch1)%p)
         CALL free_aoitem(lupri,AObatch1)
         AOT1batch1 = AOT1batch1 + batchsize1(iBatch1)
      ENDDO
      CALL free_aoitem(lupri,AObatch2)
      AOT1batch2 = AOT1batch2 + batchsize2(jBatch1)
   ENDDO
ELSE
   IF(SETTING%SCHEME%CS_SCREEN.OR.SETTING%SCHEME%PS_SCREEN)THEN
      call lsquit('error in typedef_getBatchOrbitalScreen ',lupri)
   ENDIF
ENDIF

end subroutine II_getBatchOrbitalScreen2

!> \brief build BatchOrbitalInfo
!> \author T. Kjaergaard
!> \date 2010
!> \param setting Integral evalualtion settings
!> \param AO options AORdefault,AODFdefault or AOempty
!> \param intType options contracted or primitive
!> \param nbast the number of basis functions
!> \param orbtobast for a given basis function it provides the batch index
!> \param nBatches the number of batches
!> \param lupri Default print unit
!> \param luerr Default error print unit
SUBROUTINE II_getBatchOrbitalScreenK(DecScreen,Setting,nBast,nbatches1,nbatches2,&
     & batchsize1,batchsize2,batchindex1,batchindex2,batchdim1,batchdim2,INTSPEC,lupri,luerr)
implicit none
TYPE(DECscreenITEM),intent(INOUT)    :: DecScreen
TYPE(LSSETTING) :: Setting !,intent(INOUT)
Integer,intent(IN)         :: nBast,lupri,luerr
Integer,intent(IN)         :: nBatches1,nBatches2
Integer,intent(IN)         :: batchdim1(nBatches1),batchdim2(nBatches2)
Integer,intent(IN)         :: batchsize1(nBatches1),batchsize2(nBatches2)
Integer,intent(IN)         :: batchindex1(nBatches1),batchindex2(nBatches2)
Character,intent(IN)       :: intSpec(5)
!
TYPE(AOITEM)           :: AObuild1,AObuild2
Integer                :: nDim(4),AO(4),i
!
DO i=1,4
   IF (intSpec(i).EQ.'R') THEN
      !   The regular AO-basis
      ao(i) = AORegular
   ELSE IF (intSpec(i).EQ.'C') THEN
      !   The CABS AO-type basis
      ao(i) = AOdfCABS
   ELSE
      call lsquit('Error in specification of ao1 in II_precalc_DECScreenMat',-1)
   ENDIF
ENDDO

CALL setAObatch(AObuild1,0,1,nDim(1),AO(1),Contractedinttype,Setting%scheme,Setting%fragment(1)%p,&
     &          setting%basis(1)%p,lupri,luerr)
CALL II_getBatchOrbitalScreen2K_LHS(DECscreen,setting,AObuild1,AO(2),nBatches1,batchsize1,batchindex1,batchdim1,lupri)

CALL free_aoitem(lupri,AObuild1)
CALL setAObatch(AObuild1,0,1,nDim(3),AO(3),Contractedinttype,Setting%scheme,Setting%fragment(1)%p,&
     &          setting%basis(1)%p,lupri,luerr)
CALL II_getBatchOrbitalScreen2K_RHS(DECscreen,setting,AObuild1,AO(4),nBatches2,batchsize2,batchindex2,batchdim2,lupri)

CALL free_aoitem(lupri,AObuild1)

END SUBROUTINE II_getBatchOrbitalScreenK

subroutine II_getBatchOrbitalScreen2K_LHS(DECscreen,setting,AOfull1,AOSPEC,nBatches1,&
     & batchsize1,batchindex1,batchdim1,lupri)
implicit none
TYPE(DECscreenITEM),intent(INOUT)    :: DecScreen
TYPE(LSSETTING),intent(INOUT)    :: Setting
TYPE(AOITEM),intent(IN)          :: AOfull1
Integer,intent(IN)               :: lupri,nBatches1,AOSPEC
Integer,intent(IN)               :: batchdim1(nBatches1)
Integer,intent(IN)               :: batchsize1(nBatches1)
Integer,intent(IN)               :: batchindex1(nBatches1)
!
TYPE(AOITEM)            :: AObatch1
integer :: dim1,jBatch,iBatch,I,AObatchdim1
integer :: AOT1batch1,jBatch1,iBatch1,jbat,ibat
character(len=9) :: Jlabel,Ilabel
logical :: FoundInMem,FoundOnDisk,uncont,intnrm
type(lstensor),pointer :: MasterGAB
type(lstensor),pointer :: GAB
TYPE(BASISSETINFO),pointer :: AObasis
IF(ASSOCIATED(DECSCREEN%masterGabLHS))THEN
   MasterGAB => DECSCREEN%masterGabLHS
ELSE
   IF(SETTING%SCHEME%CS_SCREEN.OR.SETTING%SCHEME%PS_SCREEN)THEN
      call lsquit('error in typedef_getBatchOrbitalScreenKLHS ',lupri)
   ENDIF
ENDIF

uncont = .FALSE.
intnrm = .FALSE.
!LHS   
IF(.NOT.ASSOCIATED(DECSCREEN%batchGabKLHS))CALL LSQUIT('DECSCREEN%batchGabKLHS not assocd',-1)
AOT1batch1 = 0

SELECT CASE(AOspec)
CASE (AORegular)
   AObasis => setting%basis(1)%p%BINFO(RegBasParam)
CASE (AOdfCABS)
   AObasis => setting%basis(1)%p%BINFO(CABBasParam)
CASE DEFAULT
   print*,'case: ',AOspec
   WRITE(lupri,*) 'case: ',AOspec
   CALL LSQuit('Programming error: Not a case in II_getBatchOrbitalScreenKLHS2B!',lupri)
END SELECT

DO iBatch1=1,nBatches1
   dim1 = batchdim1(iBatch1)      
   ibat = batchindex1(iBatch1)
   call BUILD_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
        & Setting%molecule(1)%p,aobasis,AObatch1,&
        & uncont,intnrm,iBat,AObatchdim1,batchsize1(iBatch1))
   nullify(DECSCREEN%batchGabKLHS(iBatch1)%p)         
   allocate(DECSCREEN%batchGabKLHS(iBatch1)%p)         
   call build_BatchGabK(AOfull1,AObatch1,iBatch1,AOT1batch1,dim1,&
              & MasterGAB,DECSCREEN%batchGabKLHS(iBatch1)%p)
   CALL free_aoitem(lupri,AObatch1)
   AOT1batch1 = AOT1batch1 + batchsize1(iBatch1)
ENDDO

end subroutine II_getBatchOrbitalScreen2K_LHS

subroutine II_getBatchOrbitalScreen2K_RHS(DECscreen,setting,AOfull2,aospec,nBatches2,&
     & batchsize2,batchindex2,batchdim2,lupri)
implicit none
TYPE(DECscreenITEM),intent(INOUT)    :: DecScreen
TYPE(LSSETTING),intent(INOUT)    :: Setting
TYPE(AOITEM),intent(IN)          :: AOfull2
Integer,intent(IN)               :: lupri,nBatches2,aospec
Integer,intent(IN)               :: batchdim2(nBatches2)
Integer,intent(IN)               :: batchsize2(nBatches2)
Integer,intent(IN)               :: batchindex2(nBatches2)
!
TYPE(AOITEM)            :: AObatch2
integer :: dim2,jBatch,iBatch,I,AObatchdim2
integer :: AOT1batch2,jBatch1,jbat
character(len=9) :: Jlabel,Ilabel
logical :: FoundInMem,FoundOnDisk,uncont,intnrm
type(lstensor),pointer :: MasterGAB
type(lstensor),pointer :: GAB
TYPE(BASISSETINFO),pointer :: AObasis
IF(ASSOCIATED(DECSCREEN%masterGabRHS))THEN
   MasterGAB => DECSCREEN%masterGabRHS
ELSE
   IF(SETTING%SCHEME%CS_SCREEN.OR.SETTING%SCHEME%PS_SCREEN)THEN
      call lsquit('error in typedef_getBatchOrbitalScreenKRHS ',lupri)
   ENDIF
ENDIF

SELECT CASE(AOspec)
CASE (AORegular)
   AObasis => setting%basis(1)%p%BINFO(RegBasParam)
CASE (AOdfCABS)
   AObasis => setting%basis(1)%p%BINFO(CABBasParam)
CASE DEFAULT
   print*,'case: ',AOspec
   WRITE(lupri,*) 'case: ',AOspec
   CALL LSQuit('Programming error: Not a case in II_getBatchOrbitalScreenKRHS2B!',lupri)
END SELECT

uncont = .FALSE.
intnrm = .FALSE.
!RHS
IF(.NOT.ASSOCIATED(DECSCREEN%batchGabKRHS))CALL LSQUIT('DECSCREEN%batchGabKRHS not assocd',-1)
AOT1batch2 = 0
DO jBatch1=1,nBatches2
   dim2 = batchdim2(jBatch1)
   jbat = batchindex2(jBatch1)
   call BUILD_SHELLBATCH_AO(LUPRI,SETTING%SCHEME,SETTING%SCHEME%AOPRINT,&
        & Setting%molecule(1)%p,aobasis,AObatch2,&
        & uncont,intnrm,jBat,AObatchdim2,batchsize2(jBatch1))
   nullify(DECSCREEN%batchGabKRHS(jBatch1)%p)         
   allocate(DECSCREEN%batchGabKRHS(jBatch1)%p)         
   call build_BatchGabK(AOfull2,AObatch2,jBatch1,AOT1batch2,dim2,&
              & MasterGAB,DECSCREEN%batchGabKRHS(jBatch1)%p)
   CALL free_aoitem(lupri,AObatch2)
   AOT1batch2 = AOT1batch2 + batchsize2(jBatch1)
ENDDO

end subroutine II_getBatchOrbitalScreen2K_RHS


!> \brief Calculates the Batched decpacked explicit 4 center eri
!> \author T. Kjaergaard
!> \date 2014
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the output (full,full,batchC,batchD)
!> \param batchA batch index 
!> \param batchB batch index 
!> \param batchC batch index 
!> \param batchD batch index 
!> \param dim1 the dimension of batch index 
!> \param dim2 the dimension of batch index 
!> \param dim3 the dimension of batch index 
!> \param dim4 the dimension of batch index 
!> \param intSpec Specified first the four AOs and then the operator ('RRRRC' give the standard AO ERIs)
SUBROUTINE II_GET_DECBATCHPACKED(LUPRI,LUERR,SETTING,&
     & outputintegral,batchA,batchB,batchC,batchD,&
     & batchsizeA,batchSizeB,batchsizeC,batchSizeD,&
     & dim1,dim2,dim3,dim4,intSpec,intTHRESHOLD)
IMPLICIT NONE
TYPE(LSSETTING),intent(inout) :: SETTING
INTEGER,intent(in)            :: LUPRI,LUERR,dim1,dim2,dim3,dim4
INTEGER,intent(in)            :: batchA,batchB,batchC,batchD
INTEGER,intent(in)            :: batchsizeA,batchSizeB,batchsizeC,batchSizeD
REAL(REALK),target            :: outputintegral(dim1,dim2,dim3,dim4,1)
Character,intent(IN)          :: intSpec(5)
REAL(REALK),intent(IN)        :: intTHRESHOLD
!
integer               :: I,J
type(matrixp)         :: intmat(1)
logical               :: ps_screensave,saveRecalcGab,cs_screensave
integer               :: nAO,iA,iB,iC,iD
logical,pointer       :: OLDsameMOLE(:,:),OLDsameBAS(:,:)
integer               :: oper,ao(4)
real(realk)         :: coeff(6),exponent(6),tmp
real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
integer             :: iunit,k,l,IJ
integer             :: nGaussian,nG2
real(realk)         :: GGem,slater
type(lstensor),pointer :: tmpP
!call time_II_operations1()
IF (intSpec(5).NE.'C') THEN
   nGaussian = 6
   nG2 = nGaussian*(nGaussian+1)/2
   GGem = 0E0_realk
   slater = setting%basis(1)%p%BINFO(RegBasParam)%GeminalScalingFactor
   call stgfit(slater,nGaussian,exponent,coeff)
   IJ=0
   DO I=1,nGaussian
      DO J=1,I
         IJ = IJ + 1
         coeff2(IJ) = 2E0_realk * coeff(I) * coeff(J)
         prodexponent(IJ) = exponent(I) * exponent(J)
         sumexponent(IJ) = exponent(I) + exponent(J)
      ENDDO
      coeff2(IJ) = 0.5E0_realk*coeff2(IJ)
   ENDDO
ENDIF

! ***** SELECT OPERATOR TYPE *****
IF (intSpec(5).EQ.'C') THEN
   ! Regular Coulomb operator 1/r12
   oper = CoulombOperator
ELSE IF (intSpec(5).eq.'E') then
   ! Long-Range Erf operator erf(mu r_12)/r_12
   oper = ErfOperator
ELSE IF (intSpec(5).EQ.'G') THEN
   ! The Gaussian geminal operator g
   oper = GGemOperator
   call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
ELSE IF (intSpec(5).EQ.'F') THEN
   ! The Gaussian geminal divided by the Coulomb operator g/r12
   oper = GGemCouOperator
   call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
ELSE IF (intSpec(5).EQ.'D') THEN
   ! The double commutator [[T,g],g]
   oper = GGemGrdOperator
   call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
ELSE IF (intSpec(5).EQ.'2') THEN
   ! The Gaussian geminal operator squared g^2
   oper = GGemOperator
   call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
ELSE
   call lsquit('Error in specification of operator in InitGaussianGeminal',-1)
ENDIF

! ***** SELECT AO TYPES *****
DO i=1,4
  IF (intSpec(i).EQ.'R') THEN
!   The regular AO-basis
    ao(i) = AORegular
  ELSE IF (intSpec(i).EQ.'C') THEN
!   The CABS AO-type basis
    ao(i) = AOdfCABS
  ELSE
    call lsquit('Error in specification of ao1 in II_GET_DECPACKED4CENTER_J_ERI',-1)
  ENDIF
ENDDO

call mem_alloc(OLDsameBAS,setting%nAO,setting%nAO)
OLDsameBAS = setting%sameBAS
DO I=1,4
  DO J=1,4
    setting%sameBas(I,J) = setting%sameBas(I,J) .AND. ao(I).EQ.ao(J)
  ENDDO
ENDDO

!DEC wants the integrals in (nbast1,nbast2,dim3,dim4) but it is faster 
!to calculate them as (dim3,dim4,nbast1,nbast2) 
!so we calculate them as (dim3,dim4,nbast1,nbast1) but store the 
!result as (nbast1,nbast2,dim3,dim4)
SETTING%SCHEME%intTHRESHOLD=intTHRESHOLD

nAO = setting%nAO
call mem_alloc(OLDsameMOLE,nAO,nAO)
OLDsameMOLE = setting%sameMol
setting%batchindex(1)=batchA
setting%batchindex(2)=batchB
setting%batchindex(3)=batchC
setting%batchindex(4)=batchD
setting%batchSize(1)=batchSizeA
setting%batchSize(2)=batchSizeB
setting%batchSize(3)=batchSizeC
setting%batchSize(4)=batchSizeD
setting%batchdim(1)=dim1
setting%batchdim(2)=dim2
setting%batchdim(3)=dim3
setting%batchdim(4)=dim4

setting%sameMol=.FALSE.
!not 100 % safe
!DO I=1,4
!   DO J=1,4
!      setting%sameMol(I,J)=&
!           & ((setting%batchindex(I).EQ.setting%batchindex(J)).AND.&
!           & (setting%batchSize(I).EQ.setting%batchSize(J)).AND.&
!           & (setting%batchdim(I).EQ.setting%batchdim(J)))
!   ENDDO
!ENDDO
SETTING%SCHEME%intTHRESHOLD=intTHRESHOLD
nullify(setting%output%resulttensor)
call initIntegralOutputDims(setting%output,dim1,dim2,dim3,dim4,1)
setting%output%DECPACKEDK = .TRUE.
setting%output%Resultmat => outputintegral

! Set to zero
call JZERO(setting%output%ResultMat,dim1,dim2,dim3,dim4)

CALL ls_setDefaultFragments(setting)
IF(Setting%scheme%cs_screen.OR.Setting%scheme%ps_screen)THEN
   IF(.NOT.associated(SETTING%LST_GAB_RHS))THEN
      CALL LSQUIT('SETTING%LST_GAB_RHS not associated in DEC_ERI',-1)
   ENDIF
   IF(.NOT.associated(SETTING%LST_GAB_LHS))THEN
      CALL LSQUIT('SETTING%LST_GAB_LHS not associated in DEC_ERI',-1)
   ENDIF
   IF(SETTING%LST_GAB_LHS%nbast(1).NE.dim1)call lsquit('dim mismatch BGJERI1',-1)
   IF(SETTING%LST_GAB_LHS%nbast(2).NE.dim2)call lsquit('dim mismatch BGJERI2',-1)
   IF(SETTING%LST_GAB_RHS%nbast(1).NE.dim3)call lsquit('dim mismatch BGJERI3',-1)
   IF(SETTING%LST_GAB_RHS%nbast(2).NE.dim4)call lsquit('dim mismatch BGJERI4',-1)
   IF(SETTING%LST_GAB_LHS%nbatches(1).NE.batchSizeA)call lsquit('error BatchsizeA.',-1)
   IF(SETTING%LST_GAB_LHS%nbatches(2).NE.batchSizeB)call lsquit('error BatchsizeB.',-1)
   IF(SETTING%LST_GAB_RHS%nbatches(1).NE.batchSizeC)call lsquit('error BatchsizeC.',-1)
   IF(SETTING%LST_GAB_RHS%nbatches(2).NE.batchSizeD)call lsquit('error BatchsizeD.',-1)
ENDIF
CALL ls_getIntegrals1(ao(1),ao(2),ao(3),ao(4),oper,RegularSpec,ContractedInttype,0,SETTING,LUPRI,LUERR)
call mem_dealloc(setting%output%postprocess)
!back to normal

setting%output%DECPACKEDK = .FALSE.
setting%batchindex(1)=0
setting%batchindex(2)=0
setting%batchindex(3)=0
setting%batchindex(4)=0
setting%batchSize(1)=1
setting%batchSize(2)=1
setting%batchSize(3)=1
setting%batchSize(4)=1
setting%batchdim(1)=0
setting%batchdim(2)=0
setting%batchdim(3)=0
setting%batchdim(4)=0
setting%sameFrag=.TRUE.
setting%sameMol = OLDsameMOLE
call mem_dealloc(OLDsameMOLE) 
setting%sameBas = OLDsameBAS
call mem_dealloc(OLDsameBAS) 
!call time_II_operations2(JOB_II_GET_DECBATCHPACKED)
END SUBROUTINE II_GET_DECBATCHPACKED

!> \brief Calculates the decpacked explicit 4 center eri
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the output (full,full,batchC,batchD)
!> \param batchC batch index 
!> \param batchD batch index 
!> \param nbast1 full orbital dimension of ao 1
!> \param nbast2 full orbital dimension of ao 2
!> \param dim3 the dimension of batch index 
!> \param dim4 the dimension of batch index 
!> \param intSpec Specified first the four AOs and then the operator ('RRRRC' give the standard AO ERIs)
SUBROUTINE II_GET_DECPACKED4CENTER_J_ERI(LUPRI,LUERR,SETTING,&
     &outputintegral,batchC,batchD,batchsizeC,batchSizeD,nbast1,nbast2,dim3,dim4,fullRHS,intSpec,intThreshold)
IMPLICIT NONE
TYPE(LSSETTING),intent(inout) :: SETTING
INTEGER,intent(in)            :: LUPRI,LUERR,nbast1,nbast2,dim3,dim4,batchC,batchD
INTEGER,intent(in)            :: batchsizeC,batchSizeD
REAL(REALK),target            :: outputintegral(nbast1,nbast2,dim3,dim4,1)
logical,intent(in)            :: FullRhs
Character,intent(IN)          :: intSpec(5)
REAL(REALK),intent(IN)        :: intTHRESHOLD
!
integer               :: I,J
type(matrixp)         :: intmat(1)
logical               :: ps_screensave,saveRecalcGab,cs_screensave
integer               :: nAO,iA,iB,iC,iD
logical,pointer       :: OLDsameMOLE(:,:),OLDsameBAS(:,:)
integer               :: oper,ao(4)
real(realk)         :: coeff(6),exponent(6),tmp
real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
integer             :: iunit,k,l,IJ
integer             :: nGaussian,nG2
real(realk)         :: GGem,slater
type(lstensor),pointer :: tmpP
call time_II_operations1()
IF (intSpec(5).NE.'C') THEN
   nGaussian = 6
   nG2 = nGaussian*(nGaussian+1)/2
   GGem = 0E0_realk
   slater = setting%basis(1)%p%BINFO(RegBasParam)%GeminalScalingFactor
   call stgfit(slater,nGaussian,exponent,coeff)
   IJ=0
   DO I=1,nGaussian
      DO J=1,I
         IJ = IJ + 1
         coeff2(IJ) = 2E0_realk * coeff(I) * coeff(J)
         prodexponent(IJ) = exponent(I) * exponent(J)
         sumexponent(IJ) = exponent(I) + exponent(J)
      ENDDO
      coeff2(IJ) = 0.5E0_realk*coeff2(IJ)
   ENDDO
ENDIF

! ***** SELECT OPERATOR TYPE *****
IF (intSpec(5).EQ.'C') THEN
   ! Regular Coulomb operator 1/r12
   oper = CoulombOperator
ELSE IF (intSpec(5).eq.'E') then
   ! Long-Range Erf operator erf(mu r_12)/r_12
   oper = ErfOperator
ELSE IF (intSpec(5).EQ.'G') THEN
   ! The Gaussian geminal operator g
   oper = GGemOperator
   call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
ELSE IF (intSpec(5).EQ.'F') THEN
   ! The Gaussian geminal divided by the Coulomb operator g/r12
   oper = GGemCouOperator
   call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
ELSE IF (intSpec(5).EQ.'D') THEN
   ! The double commutator [[T,g],g]
   oper = GGemGrdOperator
   call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
ELSE IF (intSpec(5).EQ.'2') THEN
   ! The Gaussian geminal operator squared g^2
   oper = GGemOperator
   call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
ELSE
   call lsquit('Error in specification of operator in InitGaussianGeminal',-1)
ENDIF

! ***** SELECT AO TYPES *****
DO i=1,4
  IF (intSpec(i).EQ.'R') THEN
!   The regular AO-basis
    ao(i) = AORegular
  ELSE IF (intSpec(i).EQ.'C') THEN
!   The CABS AO-type basis
    ao(i) = AOdfCABS
  ELSE
    call lsquit('Error in specification of ao1 in II_GET_DECPACKED4CENTER_J_ERI',-1)
  ENDIF
ENDDO

call mem_alloc(OLDsameBAS,setting%nAO,setting%nAO)
OLDsameBAS = setting%sameBAS
DO I=1,4
  DO J=1,4
    setting%sameBas(I,J) = setting%sameBas(I,J) .AND. ao(I).EQ.ao(J)
  ENDDO
ENDDO

!DEC wants the integrals in (nbast1,nbast2,dim3,dim4) but it is faster 
!to calculate them as (dim3,dim4,nbast1,nbast2) 
!so we calculate them as (dim3,dim4,nbast1,nbast1) but store the 
!result as (nbast1,nbast2,dim3,dim4)
SETTING%SCHEME%intTHRESHOLD=intTHRESHOLD
IF(.NOT.FULLRHS)THEN
   nAO = setting%nAO
   call mem_alloc(OLDsameMOLE,nAO,nAO)
   OLDsameMOLE = setting%sameMol
   setting%batchindex(1)=batchC
   setting%batchindex(2)=batchD
   setting%batchSize(1)=batchSizeC
   setting%batchSize(2)=batchSizeD
   setting%batchdim(1)=dim3
   setting%batchdim(2)=dim4
   setting%sameMol(1,2)=.FALSE.
   setting%sameMol(2,1)=.FALSE.
   DO I=1,2
      DO J=3,4
         setting%sameMol(I,J)=.FALSE.
         setting%sameMol(J,I)=.FALSE.
      ENDDO
   ENDDO
ENDIF
!saveRecalcGab = setting%scheme%recalcGab
!setting%scheme%recalcGab = .TRUE.
SETTING%SCHEME%intTHRESHOLD=intTHRESHOLD
nullify(setting%output%resulttensor)
call initIntegralOutputDims(setting%output,dim3,dim4,nbast1,nbast2,1)
setting%output%DECPACKED = .TRUE.
setting%output%Resultmat => outputintegral

! Set to zero
call JZERO(setting%output%ResultMat,nbast1,nbast2,dim3,dim4)

CALL ls_setDefaultFragments(setting)
IF(Setting%scheme%cs_screen.OR.Setting%scheme%ps_screen)THEN
   IF(.NOT.associated(SETTING%LST_GAB_RHS))THEN
      CALL LSQUIT('SETTING%LST_GAB_RHS not associated in DEC_ERI',-1)
   ENDIF
   IF(.NOT.associated(SETTING%LST_GAB_LHS))THEN
      CALL LSQUIT('SETTING%LST_GAB_LHS not associated in DEC_ERI',-1)
   ENDIF
   !we switch sides as we actuall calc the integrals in (3412) order, 
   !while still storing it as (1234) aka (bachC,batchD,full,full)
   nullify(tmpP)
   tmpP => SETTING%LST_GAB_RHS
   nullify(SETTING%LST_GAB_RHS)
   SETTING%LST_GAB_RHS => SETTING%LST_GAB_LHS
   nullify(SETTING%LST_GAB_LHS)
   SETTING%LST_GAB_LHS => tmpP

   IF(SETTING%LST_GAB_RHS%nbast(1).NE.nbast1)call lsquit('dim mismatch GJERI1',-1)
   IF(SETTING%LST_GAB_RHS%nbast(2).NE.nbast2)call lsquit('dim mismatch GJERI2',-1)
   IF(SETTING%LST_GAB_LHS%nbast(1).NE.dim3)call lsquit('dim mismatch GJERI3',-1)
   IF(SETTING%LST_GAB_LHS%nbast(2).NE.dim4)call lsquit('dim mismatch GJERI4',-1)

   if(.not. fullrhs) then
     IF(SETTING%LST_GAB_LHS%nbatches(1).NE.batchSizeC)call lsquit('error BatchsizeC.',-1)
     IF(SETTING%LST_GAB_LHS%nbatches(2).NE.batchSizeD)call lsquit('error BatchsizeD.',-1)
   end if
ENDIF
CALL ls_getIntegrals1(ao(3),ao(4),ao(1),ao(2),oper,RegularSpec,ContractedInttype,0,SETTING,LUPRI,LUERR)
call mem_dealloc(setting%output%postprocess)
setting%output%DECPACKED = .FALSE.
!back to normal
IF(.NOT.FULLRHS)THEN
   setting%batchindex(1)=0
   setting%batchindex(2)=0
   setting%batchSize(1)=1
   setting%batchSize(2)=1
   setting%batchdim(1)=0
   setting%batchdim(2)=0
   setting%sameFrag=.TRUE.
   setting%sameMol = OLDsameMOLE
   call mem_dealloc(OLDsameMOLE) 
ENDIF
setting%sameBas = OLDsameBAS
call mem_dealloc(OLDsameBAS) 
!back to normal
IF(Setting%scheme%cs_screen.OR.Setting%scheme%ps_screen)THEN
   nullify(tmpP)
   tmpP => SETTING%LST_GAB_RHS
   nullify(SETTING%LST_GAB_RHS)
   SETTING%LST_GAB_RHS => SETTING%LST_GAB_LHS
   nullify(SETTING%LST_GAB_LHS)
   SETTING%LST_GAB_LHS => tmpP
ENDIF

call time_II_operations2(JOB_II_GET_DECPACKED4CENTER_J_ERI)
END SUBROUTINE II_GET_DECPACKED4CENTER_J_ERI

!> \brief Calculates the decpacked explicit 4 center eri
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the output (full,full,batchC,batchD)
!> \param batchC batch index 
!> \param batchD batch index 
!> \param nbast1 full orbital dimension of ao 1
!> \param nbast2 full orbital dimension of ao 2
!> \param dim3 the dimension of batch index 
!> \param dim4 the dimension of batch index 
!> \param intSpec Specified first the four AOs and then the operator ('RRRRC' give the standard AO ERIs)
SUBROUTINE II_GET_DECPACKED4CENTER_J_ERI2(LUPRI,LUERR,SETTING,&
     &outputintegral,batchA,batchB,batchsizeA,batchSizeB,dim1,dim2,nbast3,nbast4,fullRHS,intSpec,intThreshold)
IMPLICIT NONE
TYPE(LSSETTING),intent(inout) :: SETTING
INTEGER,intent(in)            :: LUPRI,LUERR,nbast3,nbast4,dim1,dim2,batchA,batchB
INTEGER,intent(in)            :: batchsizeA,batchSizeB
REAL(REALK),target            :: outputintegral(dim1,dim2,nbast3,nbast4,1)
logical,intent(in)            :: FullRhs
Character,intent(IN)          :: intSpec(5)
REAL(REALK),intent(IN)        :: intTHRESHOLD
!
integer               :: I,J
type(matrixp)         :: intmat(1)
logical               :: ps_screensave,saveRecalcGab,cs_screensave
integer               :: nAO,iA,iB,iC,iD
logical,pointer       :: OLDsameMOLE(:,:),OLDsameBAS(:,:)
integer               :: oper,ao(4)
real(realk)         :: coeff(6),exponent(6),tmp
real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
integer             :: iunit,k,l,IJ
integer             :: nGaussian,nG2
real(realk)         :: GGem,slater
type(lstensor),pointer :: tmpP
call time_II_operations1()
IF (intSpec(5).NE.'C') THEN
   nGaussian = 6
   nG2 = nGaussian*(nGaussian+1)/2
   GGem = 0E0_realk
   slater = setting%basis(1)%p%BINFO(RegBasParam)%GeminalScalingFactor
   call stgfit(slater,nGaussian,exponent,coeff)
   IJ=0
   DO I=1,nGaussian
      DO J=1,I
         IJ = IJ + 1
         coeff2(IJ) = 2E0_realk * coeff(I) * coeff(J)
         prodexponent(IJ) = exponent(I) * exponent(J)
         sumexponent(IJ) = exponent(I) + exponent(J)
      ENDDO
      coeff2(IJ) = 0.5E0_realk*coeff2(IJ)
   ENDDO
ENDIF

! ***** SELECT OPERATOR TYPE *****
IF (intSpec(5).EQ.'C') THEN
   ! Regular Coulomb operator 1/r12
   oper = CoulombOperator
ELSE IF (intSpec(5).eq.'E') then
   ! Long-Range Erf operator erf(mu r_12)/r_12
   oper = ErfOperator
ELSE IF (intSpec(5).EQ.'G') THEN
   ! The Gaussian geminal operator g
   oper = GGemOperator
   call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
ELSE IF (intSpec(5).EQ.'F') THEN
   ! The Gaussian geminal divided by the Coulomb operator g/r12
   oper = GGemCouOperator
   call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
ELSE IF (intSpec(5).EQ.'D') THEN
   ! The double commutator [[T,g],g]
   oper = GGemGrdOperator
   call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
ELSE IF (intSpec(5).EQ.'2') THEN
   ! The Gaussian geminal operator squared g^2
   oper = GGemOperator
   call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
ELSE
   call lsquit('Error in specification of operator in InitGaussianGeminal',-1)
ENDIF

! ***** SELECT AO TYPES *****
DO i=1,4
  IF (intSpec(i).EQ.'R') THEN
!   The regular AO-basis
    ao(i) = AORegular
  ELSE IF (intSpec(i).EQ.'C') THEN
!   The CABS AO-type basis
    ao(i) = AOdfCABS
  ELSE
    call lsquit('Error in specification of ao1 in II_GET_DECPACKED4CENTER_J_ERI',-1)
  ENDIF
ENDDO

call mem_alloc(OLDsameBAS,setting%nAO,setting%nAO)
OLDsameBAS = setting%sameBAS
DO I=1,4
  DO J=1,4
    setting%sameBas(I,J) = setting%sameBas(I,J) .AND. ao(I).EQ.ao(J)
  ENDDO
ENDDO

!DEC wants the integrals in (nbast1,nbast2,dim3,dim4) but it is faster 
!to calculate them as (dim3,dim4,nbast1,nbast2) 
!so we calculate them as (dim3,dim4,nbast1,nbast1) but store the 
!result as (nbast1,nbast2,dim3,dim4)
SETTING%SCHEME%intTHRESHOLD=intTHRESHOLD
IF(.NOT.FULLRHS)THEN
   nAO = setting%nAO
   call mem_alloc(OLDsameMOLE,nAO,nAO)
   OLDsameMOLE = setting%sameMol
   setting%batchindex(1)=batchA
   setting%batchindex(2)=batchB
   setting%batchSize(1)=batchSizeA
   setting%batchSize(2)=batchSizeB
   setting%batchdim(1)=dim1
   setting%batchdim(2)=dim2
   setting%sameMol(1,2)=.FALSE.
   setting%sameMol(2,1)=.FALSE.
   DO I=1,2
      DO J=3,4
         setting%sameMol(I,J)=.FALSE.
         setting%sameMol(J,I)=.FALSE.
      ENDDO
   ENDDO
ENDIF
!saveRecalcGab = setting%scheme%recalcGab
!setting%scheme%recalcGab = .TRUE.
SETTING%SCHEME%intTHRESHOLD=intTHRESHOLD
nullify(setting%output%resulttensor)
call initIntegralOutputDims(setting%output,dim1,dim2,nbast3,nbast4,1)
setting%output%DECPACKED2 = .TRUE.
setting%output%Resultmat => outputintegral

! Set to zero
call JZERO(setting%output%ResultMat,dim1,dim2,nbast3,nbast4)

CALL ls_setDefaultFragments(setting)
IF(Setting%scheme%cs_screen.OR.Setting%scheme%ps_screen)THEN
   IF(.NOT.associated(SETTING%LST_GAB_RHS))THEN
      CALL LSQUIT('SETTING%LST_GAB_RHS not associated in DEC_ERI2',-1)
   ENDIF
   IF(.NOT.associated(SETTING%LST_GAB_LHS))THEN
      CALL LSQUIT('SETTING%LST_GAB_LHS not associated in DEC_ERI2',-1)
   ENDIF
   IF(SETTING%LST_GAB_LHS%nbast(1).NE.dim1)call lsquit('dim mismatch GJERI23',-1)
   IF(SETTING%LST_GAB_LHS%nbast(2).NE.dim2)call lsquit('dim mismatch GJERI24',-1)
   IF(SETTING%LST_GAB_RHS%nbast(1).NE.nbast3)call lsquit('dim mismatch GJERI21',-1)
   IF(SETTING%LST_GAB_RHS%nbast(2).NE.nbast4)call lsquit('dim mismatch GJERI22',-1)
   if(.not. fullrhs) then
     IF(SETTING%LST_GAB_LHS%nbatches(1).NE.batchSizeA)call lsquit('error BatchsizeA.',-1)
     IF(SETTING%LST_GAB_LHS%nbatches(2).NE.batchSizeB)call lsquit('error BatchsizeB.',-1)
   end if
ENDIF
CALL ls_getIntegrals1(ao(1),ao(2),ao(3),ao(4),oper,RegularSpec,ContractedInttype,0,SETTING,LUPRI,LUERR)
call mem_dealloc(setting%output%postprocess)
setting%output%DECPACKED2 = .FALSE.
!back to normal
IF(.NOT.FULLRHS)THEN
   setting%batchindex(1)=0
   setting%batchindex(2)=0
   setting%batchSize(1)=1
   setting%batchSize(2)=1
   setting%batchdim(1)=0
   setting%batchdim(2)=0
   setting%sameFrag=.TRUE.
   setting%sameMol = OLDsameMOLE
   call mem_dealloc(OLDsameMOLE) 
ENDIF
setting%sameBas = OLDsameBAS
call mem_dealloc(OLDsameBAS) 
call time_II_operations2(JOB_II_GET_DECPACKED4CENTER_J_ERI)
END SUBROUTINE II_GET_DECPACKED4CENTER_J_ERI2

subroutine JZERO(ResultMat,dim1,dim2,dim3,dim4)
implicit none
integer,intent(in) :: dim1,dim2,dim3,dim4
real(realk),intent(inout) :: ResultMat(dim1,dim2,dim3,dim4)
!
integer :: i,j,k,l
!$OMP PARALLEL DO DEFAULT(none) SHARED(dim1,dim2,dim3,&
!$OMP dim4,ResultMat) PRIVATE(i,j,k,l)
do l=1,dim4
   do k=1,dim3
      do j=1,dim2
         do i=1,dim1
            ResultMat(i,j,k,l) = 0.0E0_realk
         end do
      end do
   end do
end do
!$OMP END PARALLEL DO
end subroutine JZERO

!> \brief Calculates the decpacked explicit 4 center eri
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the output (batchA,full,batchC,full)
!> \param batchA batch index 
!> \param batchC batch index 
!> \param dim1 the dimension of batch index 
!> \param nbast2 full orbital dimension of ao 2
!> \param dim3 the dimension of batch index 
!> \param nbast4 full orbital dimension of ao 4
!> \param intSpec Specified first the four AOs and then the operator ('RRRRC' give the standard AO ERIs)
SUBROUTINE II_GET_DECPACKED4CENTER_K_ERI(LUPRI,LUERR,SETTING,&
     & outputintegral,batchA,batchC,batchsizeA,batchSizeC,&
     & dim1,nbast2,dim3,nbast4,intSpec,FULLRHS,intTHRESHOLD)
  IMPLICIT NONE
  TYPE(LSSETTING),intent(inout) :: SETTING
  INTEGER,intent(in)            :: LUPRI,LUERR,nbast2,nbast4,dim1,dim3,batchA,batchC
  INTEGER,intent(in)            :: batchsizeA,batchSizeC
  REAL(REALK),target            :: outputintegral(dim1,nbast2,dim3,nbast4,1)
  Character,intent(IN)          :: intSpec(5)
  LOGICAL,intent(IN) :: FULLRHS 
  REAL(REALK),intent(IN)        :: intTHRESHOLD
  !
  integer               :: I,J
  type(matrixp)         :: intmat(1)
  logical               :: ps_screensave,saveRecalcGab,cs_screensave
  integer               :: nAO,iA,iB,iC,iD
  logical,pointer       :: OLDsameMOLE(:,:),OLDsameBAS(:,:)
  integer               :: oper,ao(4)
  integer               :: iunit,k,l,IJ
  real(realk)         :: coeff(6),exponent(6),tmp
  real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
  integer             :: nGaussian,nG2
  real(realk)         :: GGem,slater
  call time_II_operations1()
  IF(FULLRHS)THEN
     !test
     IF(dim1.NE.nbast2)call lsquit('DECK: FULLRHS but not dim1.EQ.nbast2',-1)
     IF(dim3.NE.nbast4)call lsquit('DECK: FULLRHS but not dim3.EQ.nbast4',-1)
     IF(nbast2.NE.nbast4)call lsquit('DECK: FULLRHS but not nbast2.EQ.nbast4',-1)
     WRITE(lupri,*)'Since the full 4 dim is used we call the more efficient'
     WRITE(lupri,*)'II_GET_DECPACKED4CENTER_J_ERI to calc the full 4 dim int'
     CALL II_GET_DECPACKED4CENTER_J_ERI(LUPRI,LUERR,SETTING,&
          & outputintegral,batchA,batchC,batchsizeA,batchSizeC,&
          & nbast2,nbast4,dim1,dim3,fullRHS,intSpec,intThreshold)
  ELSE
     IF (intSpec(5).NE.'C') THEN
        nGaussian = 6
        nG2 = nGaussian*(nGaussian+1)/2
        GGem = 0E0_realk
        slater = setting%basis(1)%p%BINFO(RegBasParam)%GeminalScalingFactor
        call stgfit(slater,nGaussian,exponent,coeff)
        IJ=0
        DO I=1,nGaussian
           DO J=1,I
              IJ = IJ + 1
              coeff2(IJ) = 2E0_realk * coeff(I) * coeff(J)
              prodexponent(IJ) = exponent(I) * exponent(J)
              sumexponent(IJ) = exponent(I) + exponent(J)
           ENDDO
           coeff2(IJ) = 0.5E0_realk*coeff2(IJ)
        ENDDO
     ENDIF
     
     ! ***** SELECT OPERATOR TYPE *****
     IF (intSpec(5).EQ.'C') THEN
        ! Regular Coulomb operator 1/r12
        oper = CoulombOperator
     ELSE IF (intSpec(5).eq.'E') then
        ! Long-Range Erf operator erf(mu r_12)/r_12
        oper = ErfOperator
     ELSE IF (intSpec(5).EQ.'G') THEN
        ! The Gaussian geminal operator g
        oper = GGemOperator
        call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
     ELSE IF (intSpec(5).EQ.'F') THEN
        ! The Gaussian geminal divided by the Coulomb operator g/r12
        oper = GGemCouOperator
        call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
     ELSE IF (intSpec(5).EQ.'D') THEN
        ! The double commutator [[T,g],g]
        oper = GGemGrdOperator
        call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
     ELSE IF (intSpec(5).EQ.'2') THEN
        ! The Gaussian geminal operator squared g^2
        oper = GGemOperator
        call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
     ELSE
        call lsquit('Error in specification of operator in InitGaussianGeminal',-1)
     ENDIF

     ! ***** SELECT AO TYPES *****
     DO i=1,4
        IF (intSpec(i).EQ.'R') THEN
           !   The regular AO-basis
           ao(i) = AORegular
        ELSE IF (intSpec(i).EQ.'C') THEN
           !   The CABS AO-type basis
           ao(i) = AOdfCABS
        ELSE
           call lsquit('Error in specification of ao1 in II_GET_DECPACKED4CENTER_K_ERI',-1)
        ENDIF
     ENDDO

     call mem_alloc(OLDsameBAS,setting%nAO,setting%nAO)
     OLDsameBAS = setting%sameBAS
     DO I=1,4
        DO J=1,4
           setting%sameBas(I,J) = setting%sameBas(I,J) .AND. ao(I).EQ.ao(J)
        ENDDO
     ENDDO

     !DEC wants the integrals in (dim1,nbast2,dim3,nbast4) 
     SETTING%SCHEME%intTHRESHOLD=intTHRESHOLD

     nAO = setting%nAO
     call mem_alloc(OLDsameMOLE,nAO,nAO)
     OLDsameMOLE = setting%sameMol
     setting%batchindex(1)=batchA
     setting%batchindex(3)=batchC
     setting%batchSize(1)=batchSizeA
     setting%batchSize(3)=batchSizeC
     setting%batchdim(1)=dim1
     setting%batchdim(3)=dim3

     setting%sameMol=.FALSE.
     !setting%sameMol(2,4)=.TRUE.
     !setting%sameMol(4,2)=.TRUE.

     nullify(setting%output%resulttensor)
     call initIntegralOutputDims(setting%output,dim1,nbast2,dim3,nbast4,1)
     setting%output%DECPACKEDK = .TRUE.
     setting%output%Resultmat => outputintegral

     ! Set to zero
     call JZERO(setting%output%ResultMat,dim1,nbast2,dim3,nbast4)

     CALL ls_setDefaultFragments(setting)
     IF(Setting%scheme%cs_screen.OR.Setting%scheme%ps_screen)THEN
        IF(.NOT.associated(SETTING%LST_GAB_RHS))THEN
           CALL LSQUIT('SETTING%LST_GAB_RHS not associated in DEC_ERI K',-1)
        ENDIF
        IF(.NOT.associated(SETTING%LST_GAB_LHS))THEN
           CALL LSQUIT('SETTING%LST_GAB_LHS not associated in DEC_ERI K',-1)
        ENDIF
        IF(SETTING%LST_GAB_LHS%nbatches(1).NE.batchSizeA)call lsquit('error BatchsizeA.',-1)
        IF(SETTING%LST_GAB_RHS%nbatches(1).NE.batchSizeC)call lsquit('error BatchsizeC.',-1)

        IF(SETTING%LST_GAB_LHS%nbast(1).NE.dim1)call lsquit('dim mismatch GKERI3',-1)
        IF(SETTING%LST_GAB_LHS%nbast(2).NE.nbast2)call lsquit('dim mismatch GKERI4',-1)        
        IF(SETTING%LST_GAB_RHS%nbast(1).NE.dim3)call lsquit('dim mismatch GKERI1',-1)
        IF(SETTING%LST_GAB_RHS%nbast(2).NE.nbast4)call lsquit('dim mismatch GKERI2',-1)
     ENDIF
     CALL ls_getIntegrals1(ao(1),ao(2),ao(3),ao(4),oper,RegularSpec,ContractedInttype,0,SETTING,LUPRI,LUERR)
     call mem_dealloc(setting%output%postprocess)
     setting%output%DECPACKEDK = .FALSE.
     !back to normal

     setting%batchindex(1)=0
     setting%batchindex(3)=0
     setting%batchSize(1)=1
     setting%batchSize(3)=1
     setting%batchdim(1)=0
     setting%batchdim(3)=0
     setting%sameFrag=.TRUE.

     setting%sameMol = OLDsameMOLE
     call mem_dealloc(OLDsameMOLE) 
     setting%sameBas = OLDsameBAS
     call mem_dealloc(OLDsameBAS) 
  ENDIF
  call time_II_operations2(JOB_II_GET_DECPACKED4CENTER_K_ERI)
END SUBROUTINE II_GET_DECPACKED4CENTER_K_ERI


!> \brief Calculates the decpacked explicit 4 center eri
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the output (batchA,full,batchC,full)
!> \param batchA batch index 
!> \param batchC batch index 
!> \param dim1 the dimension of batch index 
!> \param nbast2 full orbital dimension of ao 2
!> \param dim3 the dimension of batch index 
!> \param nbast4 full orbital dimension of ao 4
!> \param intSpec Specified first the four AOs and then the operator ('RRRRC' give the standard AO ERIs)
SUBROUTINE II_GET_ERI_INTEGRALBLOCK_INQUIRE(LUPRI,LUERR,SETTING,&
     & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
     & ndimAs,ndimBs,ndimCs,ndimDs,intSpec)
  IMPLICIT NONE
  INTEGER,intent(in)            :: LUPRI,LUERR,startA,startB,startC,startD
  INTEGER,intent(in)            :: ndimA,ndimB,ndimC,ndimD
  INTEGER,intent(inout)         :: ndimAs,ndimBs,ndimCs,ndimDs
  TYPE(LSSETTING),intent(inout) :: SETTING
  Character,intent(IN)          :: intSpec(5)
  !
  integer :: batchsizeA,batchsizeB,batchsizeC,batchsizeD
  integer :: batchindexA,batchindexB,batchindexC,batchindexD
  integer :: offsetA,offsetB,offsetC,offsetD,i
  DO i=1,4
     IF (intSpec(i).NE.'R')THEN
        call lsquit('AO Error in II_GET_ERI_INTEGRALBLOCK_INQUIRE',-1)
     ENDIF
  ENDDO
  call DetermineBatchIndexAndSize(lupri,setting,startA,startB,startC,startD,&
     & ndimA,ndimB,ndimC,ndimD,batchsizeA,batchsizeB,batchsizeC,batchsizeD,&
     & batchindexA,batchindexB,batchindexC,batchindexD,&
     & offsetA,offsetB,offsetC,offsetD,ndimAs,ndimBs,ndimCs,ndimDs,'R')
END SUBROUTINE II_GET_ERI_INTEGRALBLOCK_INQUIRE

SUBROUTINE II_GET_ERI_INTEGRALBLOCK(LUPRI,LUERR,SETTING,&
     & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
     & ndimAs,ndimBs,ndimCs,ndimDs,intSpec,outputintegral,&
     & ScratchArray,intThreshold,useIchor)
  IMPLICIT NONE
  INTEGER,intent(in)            :: LUPRI,LUERR,startA,startB,startC,startD
  INTEGER,intent(in)            :: ndimA,ndimB,ndimC,ndimD
  INTEGER,intent(in)            :: ndimAs,ndimBs,ndimCs,ndimDs
  TYPE(LSSETTING),intent(inout) :: SETTING
  Character,intent(IN)          :: intSpec(5)
  REAL(REALK)                   :: outputintegral(ndimA,ndimB,ndimC,ndimD)
  REAL(REALK)                   :: ScratchArray(ndimAs,ndimBs,ndimCs,ndimDs)
  real(realk),intent(in)        :: intThreshold
  logical,intent(in) :: useIchor
  !
  integer :: batchsizeA,batchsizeB,batchsizeC,batchsizeD
  integer :: batchindexA,batchindexB,batchindexC,batchindexD
  INTEGER :: ndimAs2,ndimBs2,ndimCs2,ndimDs2
  integer :: offsetA,offsetB,offsetC,offsetD,I,A,B,C,D
  real(realk) :: t1,t2
  call LSTIMER('START',t1,t2,LUPRI)
  call DetermineBatchIndexAndSize(lupri,setting,startA,startB,startC,startD,&
       & ndimA,ndimB,ndimC,ndimD,batchsizeA,batchsizeB,batchsizeC,batchsizeD,&
       & batchindexA,batchindexB,batchindexC,batchindexD,&
       & offsetA,offsetB,offsetC,offsetD,ndimAs2,ndimBs2,ndimCs2,ndimDs2,'R')
  IF(ndimAs2.NE.ndimAs)call lsquit('dim1 mismatch in II_GET_ERI_INTEGRALBLOCK')
  IF(ndimBs2.NE.ndimBs)call lsquit('dim2 mismatch in II_GET_ERI_INTEGRALBLOCK')
  IF(ndimCs2.NE.ndimCs)call lsquit('dim3 mismatch in II_GET_ERI_INTEGRALBLOCK')
  IF(ndimDs2.NE.ndimDs)call lsquit('dim4 mismatch in II_GET_ERI_INTEGRALBLOCK')

  DO i=1,4
     IF (intSpec(i).NE.'R')call lsquit('AO Error in II_GET_ERI_INTEGRALBLOCK',-1)
  ENDDO
  !calc full integral block ScratchArray(ndimAs,ndimBs,ndimCs,ndimDs)
  IF(useIchor)THEN
     call MAIN_ICHORERI_DRIVER(lupri,0,setting,ndimAs,ndimBs,ndimCs,ndimDs,&
          & ScratchArray,INTSPEC,.FALSE.,&
          & batchindexA,batchindexA+batchsizeA-1,&
          & batchindexB,batchindexB+batchsizeB-1,&
          & batchindexC,batchindexC+batchsizeC-1,&
          & batchindexD,batchindexD+batchsizeD-1,&
          & .FALSE.,ndimAs,ndimBs,ndimCs,ndimDs,&
          & .FALSE.,intThreshold)
  ELSE
     call LSTIMER('INTEGRALBLOCK INIT',t1,t2,LUPRI)
     call II_GET_DECBATCHPACKED(LUPRI,LUERR,SETTING,&
          & ScratchArray,batchindexA,batchindexB,batchindexC,batchindexD,&
          & batchsizeA,batchSizeB,batchsizeC,batchSizeD,&
          & ndimAs,ndimBs,ndimCs,ndimDs,intSpec,intThreshold)
     call LSTIMER('INTEGRALBLOCK CALC',t1,t2,LUPRI)
  ENDIF
  IF(offsetA .EQ.0.AND.offsetB .EQ.0)THEN
   DO D=1,ndimD
    DO C=1,ndimC
     DO B=1,ndimB
      DO A=1,ndimA
       outputintegral(A,B,C,D)=ScratchArray(A,B,C+offsetC,D+offsetD)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ELSEIF(offsetA .EQ.0)THEN
   DO D=1,ndimD
    DO C=1,ndimC
     DO B=1,ndimB
      DO A=1,ndimA
       outputintegral(A,B,C,D)=ScratchArray(A,B+offsetB,C+offsetC,D+offsetD)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ELSE
   DO D=1,ndimD
    DO C=1,ndimC
     DO B=1,ndimB
      DO A=1,ndimA
       outputintegral(A,B,C,D)=ScratchArray(A+offsetA,B+offsetB,C+offsetC,D+offsetD)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDIF
  call LSTIMER('INTEGRALBLOCK COPY',t1,t2,LUPRI)
END SUBROUTINE II_GET_ERI_INTEGRALBLOCK

!!$subroutine InitGaussianGeminal(intSpec,Setting,oper)
!!$  IMPLICIT NONE
!!$  Character,intent(IN)          :: intSpec(5)
!!$  TYPE(LSSETTING),intent(inout) :: SETTING
!!$  integer,intent(inout)         :: oper
!!$  !
!!$  integer               :: I,J
!!$  real(realk)         :: coeff(6),exponent(6),tmp
!!$  real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
!!$  integer             :: iunit,k,l,IJ
!!$  integer             :: nGaussian,nG2
!!$  real(realk)         :: GGem,slater
!!$     IF (intSpec(5).NE.'C') THEN
!!$        nGaussian = 6
!!$        nG2 = nGaussian*(nGaussian+1)/2
!!$        GGem = 0E0_realk
!!$        slater = setting%basis(1)%p%BINFO(RegBasParam)%GeminalScalingFactor
!!$        call stgfit(slater,nGaussian,exponent,coeff)
!!$        IJ=0
!!$        DO I=1,nGaussian
!!$           DO J=1,I
!!$              IJ = IJ + 1
!!$              coeff2(IJ) = 2E0_realk * coeff(I) * coeff(J)
!!$              prodexponent(IJ) = exponent(I) * exponent(J)
!!$              sumexponent(IJ) = exponent(I) + exponent(J)
!!$           ENDDO
!!$           coeff2(IJ) = 0.5E0_realk*coeff2(IJ)
!!$        ENDDO
!!$     ENDIF
!!$     
!!$     ! ***** SELECT OPERATOR TYPE *****
!!$     IF (intSpec(5).EQ.'C') THEN
!!$        ! Regular Coulomb operator 1/r12
!!$        oper = CoulombOperator
!!$     ELSE IF (intSpec(5).EQ.'G') THEN
!!$        ! The Gaussian geminal operator g
!!$        oper = GGemOperator
!!$        call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
!!$     ELSE IF (intSpec(5).EQ.'F') THEN
!!$        ! The Gaussian geminal divided by the Coulomb operator g/r12
!!$        oper = GGemCouOperator
!!$        call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
!!$     ELSE IF (intSpec(5).EQ.'D') THEN
!!$        ! The double commutator [[T,g],g]
!!$        oper = GGemGrdOperator
!!$        call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
!!$     ELSE IF (intSpec(5).EQ.'2') THEN
!!$        ! The Gaussian geminal operator squared g^2
!!$        oper = GGemOperator
!!$        call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
!!$     ELSE
!!$        call lsquit('Error in specification of operator in InitGaussianGeminal',-1)
!!$     ENDIF
!!$end subroutine InitGaussianGeminal

END MODULE IntegralInterfaceDEC
