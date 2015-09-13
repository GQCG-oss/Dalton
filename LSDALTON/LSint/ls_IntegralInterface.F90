!> @file
!> Contains soubroutines that bridges integral interface routines to the main Thermite driver
MODULE ls_Integral_Interface
  use precision
  use lsparameters
  use TYPEDEFTYPE, only: lssetting, LSINTSCHEME
  use integraloutput_typetype, only: INTEGRALOUTPUT
  use integral_type, only: INTEGRALINPUT
  use Matrix_module, only: matrix, matrixp
  use Matrix_Operations, only: mtype_unres_dense, matrix_type,&
       & mtype_scalapack, mat_to_full, mat_free, mat_retrieve_block,&
       & mat_init, mat_trans, mat_daxpy, mat_scal_dia, &
       & mat_setlowertriangular_zero, mat_assign, mat_print,&
       & mat_write_to_disk, mtype_pdmm, mat_add_to_fullunres
  use matrix_operations_scalapack, only: PDM_MATRIXSYNC,&
       & free_in_darray
  use matrix_util, only : matfull_get_isym,mat_get_isym,mat_same
  use ao_typetype, only: aoitem
  use ao_type, only: free_aoitem
  use basis_typetype, only: BASISINFO, BASIS_PT, BASISSETINFO, RegBasParam, &
       & AUXBasParam,CABBasParam,JKBasParam,VALBasParam,GCTBasParam,ADMBasParam
  use molecule_typetype, only: MOLECULE_PT, MOLECULEINFO
  use lstensor_typetype, only: lstensor
  use lstensor_operationsmod, only: lstensor_nullify, &
       & build_sublstensor_from_full_lstensor, lstensor_free,&
       & copy_lstensor_to_lstensor, init_pcharge_lstensor,&
       & init_gradientlstensor, init_lstensor_1dim, init_cs_lstensor,&
       & init_ps_lstensor, set_lst_maxprimgab, init_mbie_lstensor_5dim,&
       & init_lstensor_5dim, build_lstensor_from_full_3dim,&
       & build_nuclearlstensor, lstensor_buildatomicgab,&
       & set_lst_maxgabelms, set_lst_maxprimgabelms,&
       & LSTENSOR_FULL_SYMMAT_FROM_TRIANGULARMAT,&
       & add_lstensor_to_lstensor, add_sublstensor_to_full_lstensor,&
       & lstensor_zero, memdist_lstensor_setuplocalinfo,&
       & memdist_lstensor_setupfullinfo, lstensor_zero_lowertriangular,&
       & lstensor_force_symmat_to_triangularmat, lstensor_local_free,&
       & memdist_lstensor_buildfromscalapack, &
       & Build_lst_from_matarray, build_empty_sublstensor,&
       & alloc_build_empty_sublstensor
  use integraloutput_type, only: initintegraloutputdims1,nullifyintegraloutput
  use TYPEDEF, only: getNbasis, set_sameallfrag, LS_FREEMMBUF, ls_emptyibuf,&
       & ls_emptyrbuf, ls_emptynucbuf, ls_fillnucbuf, ls_initmmbuf,&
       & retrieve_screen_output, typedef_free_setting,&
       & print_reduced_screening_info, retrieve_output_slave
  use memory_handling, only: mem_alloc, mem_dealloc, &
       & mem_deallocated_mem_type_matrix
  use Integralinfo, only: init_integral_input
  use integraldriver, only: main_integral_driver, jmatclassical, gradclassical,&
       & jmatclassicalmat, gradclassicalgrad, electronnuclearclassic
  use MBIEintegraldriver, only: mbie_integral_driver
  use BUILDAOBATCH, only: build_empty_ao, build_empty_nuclear_ao,&
       & build_empty_pcharge_ao, build_ao, build_shellbatch_ao, &
       & BUILD_EMPTY_ELFIELD_AO, build_empty_single_nuclear_ao,&
       & determinenaobatches
  use lstiming, only: lstimer, print_timers
  use io, only: io_get_filename, io_get_csidentifier
  use screen_mod, only: determine_lst_in_screenlist, screen_associate,&
       & screen_add_associate_item
  use molecule_module, only: build_fragment, freeMolecularOrbitalInfo,&
       & freeDaltonFragments
  use files,only: lsclose, lsopen
  use SphCart_Matrices, only: spherical_transformation
  use Thermite_OD, only: getTotalGeoComp
#if VAR_MPI
  use lsmpi_type, only:LSGETINT,LSJENGIN,LSLINK, ls_mpibcast, lsmpi_barrier, &
       & lsmpi_reduction, get_MPI_COMM_SELF
  use lsmpi_op, only: LSTASK, LS_TASK_MANAGER, LSMPI_TASK_LIST,&
  & lsmpi_lstensor_reduction, lsmpi_probe_and_irecv_add_lstmemrealkbuf,&
  & lsmpi_isend_lstmemrealkbuf, lsmpi_blocking_recv_add_lstmemrealkbuf
  use lsmpi_mod, only: lsmpi_set_MPI_task_list, lsmpi_set_task_manager,&
       & getTaskDimension, SetTaskFragments, ls_free_tasks,&
       & lsmpi_free_MPI_task_list, lsmpi_time_harvest,&
       & set_reduced_screen_info, ls_print_tasks
  use infpar_module
#endif

  INTERFACE ls_attachDmatToSetting
     MODULE PROCEDURE ls_attachDmatToSetting_matsingle,&
                    & ls_attachDmatToSetting_matsinglep,&
                    & ls_attachDmatToSetting_matarray,&
                    & ls_attachDmatToSetting_matarrayp,&
                    & ls_attachDmatToSetting_full
  END INTERFACE
  public::  ls_getIntegrals,ls_get_exchange_mat,ls_get_exchange_mat_serial,&
       & ls_get_coulomb_mat, ls_get_coulomb_and_exchange_mat, ls_jengine,&
       & ls_multipolemoment, ls_attachDmatToSetting, ls_freeDmatFromSetting,&
       & ls_same_mats, ls_jengineclassicalgrad, ls_getnucscreenintegrals,&
       & setaobatch, ls_getscreenintegrals1, ls_setDefaultFragments,&
       & ls_getIntegrals1, ls_jengineClassicalMAT, ls_attach_gab_to_setting,&
       & ls_free_gab_from_setting, ls_subScreenAtomic, ls_subscreenfromlist,&
       & ls_LHSSameAsRHSDmatToSetting,ls_LHSSameAsRHSDmatToSetting_deactivate
  private

CONTAINS

!> \brief Does Sherical Transformation 
!> \author T. Kjaergaard
!> \date 2010
!> \param SubBlockInt the integrals to be transformed
!> \param nbast1 number of basis functions for the 1. center
!> \param nbast2 number of basis functions for the 2. center
!> \param ndMAT1 number of matrices before the transformation
!> \param ndmat2 number of matrices after the transformation
!> \param ngeoderivcomp level of derivative 
!> \param iprint level of output printing  
!> \param lupri Default print unit
SUBROUTINE ShericalMomentTransformation(SubBlockInt,nbast1,nbast2,ndMAT1,ndmat2,ngeoderivcomp,&
     &                                  iprint,lupri)
implicit none
Type(LSSETTING)     :: SETTING
Integer             :: nbast1,nbast2,ndMAT1,ndmat2,ngeoderivcomp,iprint,lupri
Real(realk),pointer :: SubBlockInt(:,:,:,:,:)
!
Real(realk),pointer :: Integrals2(:,:,:,:,:)
!
call mem_alloc(Integrals2,nbast1,nbast2,1,1,ndmat2)
CALL LS_DZERO(Integrals2,nbast1*nbast2*ndmat2)
CALL SPHERICAL_TRANSFORMATION(SubBlockInt,Integrals2,nbast1,nbast2,ndMAT1,ndmat2,ngeoderivcomp,&
     &                        iprint,lupri)
call mem_dealloc(SubBlockInt)
SubBlockInt => Integrals2
END SUBROUTINE ShericalMomentTransformation

!> \brief Adds a sub block to a full tensor
!> \author S. Reine
!> \date 2010
!> \param integrals the full tensir
!> \param SubBlockInt the subblock to be added
!> \param sameAOsLHS
!> \param sameFragmentLHS
!> \param sameAOsRHS
!> \param sameFragmentRHS
!> \param nbast1 number of basis functions for the 1. center
!> \param nbast2 number of basis functions for the 2. center
!> \param nbast3 number of basis functions for the 3. center
!> \param nbast4 number of basis functions for the 4. center
!> \param start1 the start index for the 1. dimension 
!> \param start2 the start index for the 2. dimension 
!> \param start3 the start index for the 3. dimension 
!> \param start4 the start index for the 4. dimension 
!> \param end1 the end index for the 1. dimension 
!> \param end2 the end index for the 2. dimension 
!> \param end3 the end index for the 3. dimension 
!> \param end4 the end index for the 4. dimension 
!> \param ndmat the number of matrices (the fifth dimension)
SUBROUTINE addSubBlocks(integrals,SubBlockInt,sameAOsLHS,sameFragmentLHS,&
     &                  sameAOsRHS,sameFragmentRHS,nbast1,nbast2,nbast3,nbast4,&
     &                  start1,start2,start3,start4,end1,end2,end3,end4,ndmat)
implicit none
Real(realk), pointer :: Integrals(:,:,:,:,:)
Real(realk), pointer :: SubBlockInt(:,:,:,:,:)
Integer              :: ndmat
type(matrix)         :: tmp
Integer              :: nbast1,nbast2,nbast3,nbast4
Integer              :: start1,start2,start3,start4
Integer              :: end1,end2,end3,end4
Logical              :: sameAOsLHS,sameAOsRHS
Logical              :: sameFragmentLHS,sameFragmentRHS

! Using full integral matrices - insert sub-blocks into final integrals
  CALL addSubBlockFull(integrals,SubBlockInt,sameAOsLHS,sameFragmentLHS,&
     &                 sameAOsRHS,sameFragmentRHS,nbast1,nbast2,nbast3,nbast4,&
     &                 start1,start2,start3,start4,end1,end2,end3,end4,ndmat)

END SUBROUTINE addSubBlocks

!> \brief Adds a sub gradient to a full gradient
!> \author T. Kjærgaard
!> \date 2010
!> \param GRAD the full gradient
!> \param SubGRAD the sub gradient
!> \param ndir number of directions (should be 3)
!> \param natomsFrag number of atoms in fragements (can be more than natoms)
!> \param setting Integral evalualtion settings
!> \param sameFrag Specifies if all fragments are identical or not
SUBROUTINE AddSubGradient(GRAD,ndir,natoms,SubGRAD,ndirF,natomsFrag,setting,sameFrag,&
    &                     ILHS,IRHS,lupri)
implicit none
Logical,intent(IN) :: sameFrag
Type(LSSETTING)    :: SETTING
Integer     :: ndir,natoms,ndirF,natomsFrag,lupri,ILHS,IRHS
Real(realk) :: GRAD(ndir,natoms)
Real(realk) :: SubGRAD(ndirF,natomsFrag) 
integer :: iFrag1,iFrag2,iFrag3,iFrag4,nAtoms1,nAtoms2,nAtoms3,nAtoms4

  iFrag1 = setting%FRAGMENTS%LHSblock%blocks(ILHS)%fragment1
  iFrag2 = setting%FRAGMENTS%LHSblock%blocks(ILHS)%fragment2
  iFrag3 = setting%FRAGMENTS%RHSblock%blocks(IRHS)%fragment1
  iFrag4 = setting%FRAGMENTS%RHSblock%blocks(IRHS)%fragment2
  nAtoms1 = setting%FRAGMENTS%INFO(1)%p%nAtoms(iFrag1)
  nAtoms2 = setting%FRAGMENTS%INFO(2)%p%nAtoms(iFrag2)
  nAtoms3 = setting%FRAGMENTS%INFO(3)%p%nAtoms(iFrag3)
  nAtoms4 = setting%FRAGMENTS%INFO(4)%p%nAtoms(iFrag4)
  IF(ndir.NE. 3)call lsquit('error in addSubGradient 1',-1)
  IF(ndirF.NE. 3)call lsquit('error in addSubGradient 2',-1)
! Using full integral matrices - insert sub-blocks into final integrals
  CALL addSubGradient2(GRAD,natoms,SubGrad,natomsFrag,nAtoms1,nAtoms2,nAtoms3,nAtoms4,sameFrag,&
       & setting%FRAGMENTS%INFO(1)%p%atomicindex(1:natoms,iFrag1),&
       & setting%FRAGMENTS%INFO(2)%p%atomicindex(1:natoms,iFrag2),&
       & setting%FRAGMENTS%INFO(3)%p%atomicindex(1:natoms,iFrag3),&
       & setting%FRAGMENTS%INFO(4)%p%atomicindex(1:natoms,iFrag4),lupri)

END SUBROUTINE AddSubGradient

!> \brief Adds a sub gradient to a full gradient
!> \author T. Kjærgaard
!> \date 2010
!> \param GRAD the full gradient
!> \param SubGRAD the sub gradient
!> \param ndir number of directions (should be 3)
!> \param natoms number of atoms in molecule
!> \param ndir number of directions (should be 3)
!> \param natomsFrag number of atoms in fragements (can be more than natoms)
!> \param setting Integral evalualtion settings
!> \param sameFrag Specifies if all fragments are identical or not
SUBROUTINE addSubGradient2(GRAD,natoms,SubGrad,natomsFrag,nAtoms1,nAtoms2,&
     & nAtoms3,nAtoms4,sameFrag,atomicindex1,atomicindex2,atomicindex3,atomicindex4,lupri)
implicit none
Logical,intent(IN) :: sameFrag
Integer      :: natoms,natomsFrag,lupri
integer      :: nAtoms1,nAtoms2,nAtoms3,nAtoms4
integer      :: atomicindex1(natoms),atomicindex2(natoms)
integer      :: atomicindex3(natoms),atomicindex4(natoms)
Real(realk)  :: GRAD(3,natoms)
Real(realk)  :: SubGRAD(3,natomsFrag) 
!
integer      :: iAtom1,iAtom2,iAtom3,iAtom4,iatom,ImoleculeAtom

Iatom=0
DO Iatom1=1,nAtoms1
   Iatom = Iatom+1
   ImoleculeAtom = atomicIndex1(iatom1)
   GRAD(1,imoleculeatom) = GRAD(1,imoleculeatom)+SubGrad(1,iatom)
   GRAD(2,imoleculeatom) = GRAD(2,imoleculeatom)+SubGrad(2,iatom)
   GRAD(3,imoleculeatom) = GRAD(3,imoleculeatom)+SubGrad(3,iatom)
ENDDO
IF (.NOT.sameFrag) THEN
  DO Iatom2=1,nAtoms2
     Iatom = Iatom+1
     ImoleculeAtom = atomicIndex2(iatom2)
     GRAD(1,imoleculeatom) = GRAD(1,imoleculeatom)+SubGrad(1,iatom)
     GRAD(2,imoleculeatom) = GRAD(2,imoleculeatom)+SubGrad(2,iatom)
     GRAD(3,imoleculeatom) = GRAD(3,imoleculeatom)+SubGrad(3,iatom)
  ENDDO
  DO Iatom3=1,nAtoms3
     Iatom = Iatom+1
     ImoleculeAtom = atomicIndex3(iatom3)
     GRAD(1,imoleculeatom) = GRAD(1,imoleculeatom)+SubGrad(1,iatom)
     GRAD(2,imoleculeatom) = GRAD(2,imoleculeatom)+SubGrad(2,iatom)
     GRAD(3,imoleculeatom) = GRAD(3,imoleculeatom)+SubGrad(3,iatom)
  ENDDO
  DO Iatom4=1,nAtoms4
     Iatom = Iatom+1
     ImoleculeAtom = atomicIndex4(iatom4)
     GRAD(1,imoleculeatom) = GRAD(1,imoleculeatom)+SubGrad(1,iatom)
     GRAD(2,imoleculeatom) = GRAD(2,imoleculeatom)+SubGrad(2,iatom)
     GRAD(3,imoleculeatom) = GRAD(3,imoleculeatom)+SubGrad(3,iatom)
  ENDDO
ENDIF

END SUBROUTINE addSubGradient2


!> \brief Generalized wrapper routine to get explicit integrals for given operator Oper (in case of MPI) 
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis AORegular,AOEmpty,AOdfAux,AONuclear,AOpCharge on center 1
!> \param AO2 the type of the basis AORegular,AOEmpty,AOdfAux,AONuclear,AOpCharge on center 2
!> \param AO3 the type of the basis AORegular,AOEmpty,AOdfAux,AONuclear,AOpCharge on center 3
!> \param AO4 the type of the basis AORegular,AOEmpty,AOdfAux,AONuclear,AOpCharge on center 4
!> \param Oper the label for the operator like Coulomb,Overlap,Kinetic)
!> \param Spec the label for the type of calc Regular or Gradient
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> Generalized routine to get explicit integrals for given operator Oper 
!>
!>     Coulomb       1/r_12
!>     Overlap       dirac_delta(r_1-r_2)
!>     Kinetic       -1/2 nabla
!>
!> and with given AOs according to AO1,AO2,AO3 and AO4:
!>      Name    int parameter 
!>     Regular      AORegular       Regular basis
!>     DF-Aux       AOdfAux         Auxiliary basis for density-fitting
!>     DF-CABS      AOdfCABS        Complementary Auxiliary basis for F12
!>     DF-JK        AOdfJK          Density-fitting basis set for Fock matrix for F12
!>     ADMM         AOadmm          Auxiliary Density matrix method basis set
!>     VALENCE      AOVAL           Regular Level 2  or Valence basis 
!>     Empty        AOEmpty         Empty, used for two and three-center integrals
!>
!> AO1 and AO2 belong to electron 1, and AO3 and AO4 to electron 2
!>
!> First four dimensions of Integrals should correspond to the dimensions of 
!> AO1-4, the fifth dimension is an added index allowing for example for 
!> different derivative integrals.
!>
SUBROUTINE ls_getIntegrals(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR,geoDerivOrder)
implicit none
Integer              :: LUPRI,LUERR
type(matrix)         :: tmp
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
integer,optional     :: geoDerivOrder
!
Type(LSSETTING)      :: SETTING
Integer              :: I,J
Integer              :: nbast1,nbast2,nbast3,nbast4
Integer              :: start1,start2,start3,start4,ndim
Integer              :: end1,end2,end3,end4
Logical              :: sameFragmentLHS,sameFragmentRHS
Real(realk),pointer  :: SubBlockTrans(:,:)
Real(realk),pointer  :: DblockLHS(:,:,:),DblockRHS(:,:,:)
Logical              :: FRAGMENT_SAVE,permuteLHS,permuteRHS,IntegralTransformGC
Real(realk)          :: ts,te
integer              :: ndim2(5),natoms1,natoms2,natoms3,natoms4,geoOrder
type(lstensor), pointer :: resultTensor
Logical                 :: rhsCS_created,lhsCS_created
Logical                 :: lhs_created,rhs_created,PermuteResultTensor,doscreen

!
type(lstensor),pointer  :: dmat_lhs,dmat_rhs,result_tensor
type(lstensor),pointer  :: dmat_lhs_full,dmat_rhs_full,result_tensor_full
type(lstensor),pointer  :: gabCS_rhs,gabCS_lhs
type(lstensor),pointer  :: gabCS_rhs_full,gabCS_lhs_full
character(len=7)           :: LSTYPE
!
#ifdef VAR_MPI
Integer                    :: itask,nAtoms(4),ntasks_lhs,ntasks_rhs,iAO,atomDim
Real(realk),pointer        :: Integrals(:,:,:,:,:)
Real(realk),pointer        :: SubBlockInt(:,:,:,:,:)
Real(realk),pointer        :: SubBlock(:,:)
Logical                    :: sameAOsLHS,sameAOsRHS,sameODs,sameAllMOL,sameAllFRAG
Logical                    :: lhs_aux,rhs_aux,noA,noB
Real(realk),pointer        :: timeMatrix(:,:)
integer,pointer            :: frag1(:,:),frag2(:,:)
type(ls_task_manager)      :: tasks
type(lsmpi_task_list)         :: task_list
type(lstask),pointer       :: lhs_current,rhs_current
logical                    :: sameMolSave(4,4),MBIE_SCREENINT,BOTHpartioning
Logical                    :: saveCSscreen,savePSscreen,CS_screen,PS_screen
integer                    :: ndim_full(5)
real(realk)                :: t(8),t1,t2,t3,t4
real(realk)                :: part(2)
Logical                    :: SMasterWakeSlaves
integer(kind=ls_mpik)      :: Snode,SNumnodes,SComm
IF(.NOT.setting%scheme%doMPI)THEN
 call deactivateIntegralMPI(Setting,Snode,SNumnodes,SComm,SMasterWakeSlaves)
ENDIF

CALL LS_GETTIM(t1,t2)
t = 0E0_realk
#endif
CALL LSTIMER('START',TS,TE,6)
!IF(Oper .EQ. 'Nucpot' .AND. SETTING%SCHEME%FMM)THEN
!   CALL LSTIMER('START',TS,TE,LUPRI)
!ENDIF
IF(matrix_type.EQ.mtype_scalapack.or.matrix_type.EQ.mtype_pdmm)THEN
   IF (setting%node.EQ.0) THEN
      !change into full format on master 
      call SetScalapackDmatToFull(setting,.TRUE.,.TRUE.)
   ENDIF
ENDIF

ndim2 = setting%output%ndim
IntegralTransformGC = .FALSE.
CALL ls_setDefaultFragments(setting)
LSTYPE = 'FULLINT'
IF(Spec.EQ.GeoDerivCoulombSpec)THEN
   LSTYPE = 'AB_TYPE'
ENDIF
CALL ls_create_lstensor_full(setting,LSTYPE,AO1,AO2,AO3,AO4,Oper,Spec,intType,&
     & result_tensor_full,dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created,&
     & gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created,&
     & PermuteResultTensor,doscreen,lupri,luerr,.FALSE.,.TRUE.)

geoOrder = 0
IF(present(geoDerivOrder)) geoOrder = geoDerivOrder

#ifdef VAR_MPI

!write(6,*) 'timer for ',setting%node
!CALL LSTIMER('ls_create GI',TS,TE,6)

ndim_full = setting%output%ndim
CS_screen = setting%scheme%CS_SCREEN
PS_screen = setting%scheme%PS_SCREEN

! ***************************************************************************
! *                                MPI Specific                             *
! ***************************************************************************
IF (setting%scheme%MasterWakeSlaves.AND.setting%node.EQ.infpar%master) THEN
   !if gradient test that all molecules are same otherwise quit. 
   !   Brano: Spawn here!
   call ls_mpibcast(LSGETINT,infpar%master,setting%comm)
   call lsmpi_getIntegrals_masterToSlave(AO1,AO2,AO3,AO4,Oper,Spec,intType,geoOrder,SETTING,LUPRI,LUERR)
!write(6,*) 'timer for ',setting%node
!CALL LSTIMER('ls_master',TS,TE,6)
ENDIF

sameMolSave = Setting%sameMol
Setting%sameMol = .FALSE.
Setting%sameFrag = .FALSE.

CALL LS_GETTIM(t3,t4)
t(7) = t(7) + t3 - t1 !Node CPU time  += task CPU time
t(8) = t(8) + t4 - t2 !Node wall time += task wall time
CALL LS_GETTIM(t1,t2)
!CALL LSTIMER('START',TS,TE,6)
BOTHpartioning = .FALSE.
CALL lsmpi_set_task_manager(tasks,'JENGINE',AO1,AO2,AO3,AO4,Spec,intType,setting,BOTHpartioning,lupri,luerr)
CALL lsmpi_set_MPI_task_list(task_list,tasks,setting,lupri,luerr)

part(1) = 0E0_realk
part(2) = 0E0_realk
!write(6,*) 'timer for ',setting%node
!CALL LSTIMER('mpitask',TS,TE,6)

DO itask=1,task_list%numMPItasks
  I=task_list%taskPair(1,itask)
  J=task_list%taskPair(2,itask)
  lhs_current => task_list%LHS(I)%p
  rhs_current => task_list%RHS(J)%p
! CALL LSTIMER('START',TS,TE,6)
  CALL SetTaskFragments(SETTING,lhs_current,'LHS',tasks%lhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
  CALL SetTaskFragments(SETTING,rhs_current,'RHS',tasks%rhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
  CALL getTaskDimension(tasks%orbInfo,lhs_current,AO1,AO2,'LHS',Spec,nbast1,nbast2,lupri)
  CALL getTaskDimension(tasks%orbInfo,rhs_current,AO3,AO4,'RHS',Spec,nbast3,nbast4,lupri)
  call SET_SAMEALLFRAG(sameAllFrag,setting%sameFrag,setting%nAO)
! CALL LSTIMER('gi-set-tf',TS,TE,6)
  CALL ls_create_lstensor_task(setting,result_tensor,LSTYPE,dmat_lhs,dmat_rhs,&
       & result_tensor_full,dmat_lhs_full,dmat_rhs_full,&
       & gabCS_rhs_full,gabCS_lhs_full,gabCS_rhs,gabCS_lhs,&
       & rhsCS_created,lhsCS_created,&
       & nbast1,nbast2,nbast3,nbast4,lhs_created,rhs_created,lhs_current,rhs_current,&
       & tasks%sameAOsLHS,tasks%sameAOsRHS,tasks%sameODs,sameAllFrag,CS_screen,PS_screen,doscreen)
! CALL LSTIMER('gi-mpi-clst',TS,TE,6)
  CALL LS_GETTIM(t3,t4)
  t(1) = t(1) + t3 - t1 !Node CPU time  += task CPU time
  t(2) = t(2) + t4 - t2 !Node wall time += task wall time
  CALL LS_GETTIM(t1,t2)
#else
   ! ***************************************************************************
   ! *                                  Serial                                 *
   ! ***************************************************************************
  dmat_lhs      => dmat_lhs_full
  dmat_rhs      => dmat_rhs_full
  result_tensor => result_tensor_full
  gabCS_rhs     => gabCS_rhs_full
  gabCS_lhs     => gabCS_lhs_full
#endif
  ! ***************************************************************************
  ! *       Actual calculation of integrals (Both Serial and for MPI)         *
  ! ***************************************************************************
  CALL ls_attach_lstensors_to_setting(setting,result_tensor,dmat_lhs,dmat_rhs,lhs_created,rhs_created)
  CALL ls_attach_gab_to_setting(setting,gabCS_lhs,gabCS_rhs)
  !*** CALCULATE INTEGRALS: tensor sub-block for MPI ***
  CALL ls_getIntegrals1(AO1,AO2,AO3,AO4,Oper,Spec,intType,geoOrder,SETTING,LUPRI,LUERR)
  call ls_free_lstensors_from_setting(setting,lupri)
  CALL ls_free_gab_from_setting(setting,lupri)
#ifdef VAR_MPI
  CALL LS_GETTIM(t3,t4)
  t(3) = t(3) + t3 - t1 !Node CPU time  += task CPU time
  t(4) = t(4) + t4 - t2 !Node wall time += task wall time
  CALL LS_GETTIM(t1,t2)
  part(1) = part(1) + lhs_current%task_part
  part(2) = part(2) + rhs_current%task_part
! CALL LSTIMER('gi-mpi-1',TS,TE,6)
  ! ***************************************************************************
  ! *                                  Serial                                 *
  ! ***************************************************************************
  ! ***************************************************************************
  ! *                                MPI Specific                             *
  ! ***************************************************************************
  call ls_extract_and_annihilate_lstensor_task(setting,result_tensor,LSTYPE,dmat_lhs,dmat_rhs,result_tensor_full,&
       & dmat_lhs_full,dmat_rhs_full,nbast1,nbast2,nbast3,nbast4,lhs_created,&
       & rhs_created,gabCS_rhs,gabCS_lhs,rhsCS_created,&
       & lhsCS_created,lhs_current,rhs_current,&
       & tasks%sameAOsLHS,tasks%sameAOsRHS,tasks%sameODs,SameAllFrag,CS_screen,PS_screen,doscreen)
! write(6,*) 'debug-timing:timer results for ',setting%node
! CALL LSTIMER('gi-mpi-ex_an',TS,TE,6)
ENDDO !itask

CALL lsmpi_free_MPI_task_list(task_list)
call ls_free_tasks(tasks%lhs,lupri)
call ls_free_tasks(tasks%rhs,lupri)
DO iAO=1,4
  call freeMolecularOrbitalInfo(tasks%orbInfo(iAO))
ENDDO

setting%sameMol = sameMolSave
Setting%sameFrag = sameMolSave

! ***************************************************************************
! *                                MPI Barrier                              *
! ***************************************************************************
call lsmpi_barrier(setting%comm)

IF(Setting%scheme%cs_int.or.Setting%scheme%ps_int)THEN
   call set_lst_maxgabelms(result_tensor_full)
   call set_lst_maxprimgabelms(result_tensor_full)
   call lsmpi_barrier(setting%comm)
   setting%output%screenTensor => result_tensor_full
   CALL lsmpi_lstensor_reduction(setting%output%screenTensor,infpar%master,setting%node,setting%comm)
ELSE
   setting%output%resultTensor => result_tensor_full
   !CALL LSTIMER('START',TS,TE,6)
   CALL lsmpi_lstensor_reduction(setting%output%resultTensor,infpar%master,setting%node,setting%comm)
ENDIF

!write(6,*) 'debug-timing:timer results for ',setting%node
!CALL LSTIMER('gi-mpi-red',TS,TE,6)
CALL LS_GETTIM(t3,t4)
t(5) = t(5) + t3 - t1 !Node CPU time  += task CPU time
t(6) = t(6) + t4 - t2 !Node wall time += task wall time
call lsmpi_time_harvest(setting%node,t(1),t(2),t(3),t(4),t(5),t(6),t(7),t(8),&
     & part(1),part(2),setting%numnodes,lupri,setting%scheme%intprint-1)

IF (setting%node.NE.infpar%master) THEN
   !put slave to sleep
   call ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
   if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)
   IF(Setting%scheme%cs_int.OR.Setting%scheme%ps_int)THEN
      call lstensor_free(setting%output%screenTensor)
      deallocate(setting%output%screenTensor)
      nullify(setting%output%screenTensor)
   ELSE
      call lstensor_free(setting%output%resultTensor)
      deallocate(setting%output%resultTensor)
      nullify(setting%output%resultTensor)
   ENDIF
   IF(associated(setting%output%postprocess))THEN
      call mem_dealloc(setting%output%postprocess)
   ENDIF
   call typedef_free_setting(SETTING)
!write(6,*) 'debug-timing:timer results for ',setting%node
!CALL LSTIMER('gi-mpi-free',TS,TE,6)
ELSE

!Permutational Symmetry
!In order to exploit permutational symmetry we set the lower triangular Dmatrix to zero and multiply the upper by 2. We set the lower triangular screening matrices to zero and we now copy the transpose of the upper triangular part to the lower triangular part. 
IF(PermuteResultTensor)THEN
   call lstensor_full_symMat_from_triangularMat(setting%output%resultTensor)
ENDIF

Call freeDaltonFragments(SETTING)
setting%output%ndim = ndim_full
!write(6,*) 'debug-timing:timer results for ',setting%node
!CALL LSTIMER('gi-mpi-permute',TS,TE,6)
#else
IF(Setting%scheme%cs_int.OR.Setting%scheme%ps_int)THEN
   call lstensor_free(result_tensor_full)
   deallocate(result_tensor_full)
   nullify(result_tensor_full)
ENDIF
#endif

IF(Oper .EQ. NucpotOperator .AND. SETTING%SCHEME%FMM)THEN
   IF(Spec .EQ. RegularSpec)THEN
      CALL LSTIMER('ea-non',TS,TE,LUPRI)
      CALL ls_attach_lstensors_to_setting(setting,result_tensor_full,dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
      CALL ls_electronNuclearClassic(AO1,AO2,AO3,AO4,Oper,intType,SETTING,LUPRI,LUERR)
      call ls_free_lstensors_from_setting(setting,lupri)
      CALL LSTIMER('ea-cls  ',TS,TE,LUPRI)
   ELSE
      ! Nuclear-electron attraction gradient contribution together with electron-electron repulsion
   ENDIF
ENDIF

CALL ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)

#ifdef VAR_MPI
IF(.NOT.setting%scheme%doMPI)THEN
   call ReactivateIntegralMPI(Setting,Snode,SNumnodes,SComm,SMasterWakeSlaves)
ENDIF

ENDIF
#endif
END SUBROUTINE ls_getIntegrals

#ifdef VAR_MPI
subroutine DeactivateIntegralMPI(Setting,Savenode,SaveNumnodes,SaveComm,&
     & SaveMasterWakeSlaves)
  implicit none
  Type(LSSETTING),intent(inout)       :: SETTING
  Logical,intent(inout)               :: SaveMasterWakeSlaves
  integer(kind=ls_mpik),intent(inout) :: Savenode,SaveNumnodes,SaveComm   
  !this means that this subroutine (ls_getIntegrals) have been
  !call from a slave and it wants to calculate the full contribution
  SaveMasterWakeSlaves = setting%scheme%MasterWakeSlaves
  Savenode = setting%node
  SaveNumnodes = setting%numnodes
  SaveComm = setting%comm

  setting%scheme%MasterWakeSlaves = .FALSE.
  setting%node = infpar%master
  setting%numnodes = 1_ls_mpik
  !MPI_COMM_SELF is the local comm which only contains the rank itself
  call GET_MPI_COMM_SELF(setting%comm) 
END subroutine DeactivateIntegralMPI

subroutine ReactivateIntegralMPI(Setting,Savenode,SaveNumnodes,SaveComm,&
     & SaveMasterWakeSlaves)
  implicit none
  Type(LSSETTING),intent(inout)    :: SETTING
  Logical,intent(in)               :: SaveMasterWakeSlaves
  integer(kind=ls_mpik),intent(in) :: Savenode,SaveNumnodes,SaveComm   
  !Revert Back
   setting%scheme%MasterWakeSlaves = SaveMasterWakeSlaves
   setting%node = Savenode
   setting%numnodes = SaveNumnodes
   setting%comm = SaveComm
END subroutine ReactivateIntegralMPI
#endif

!> \brief Generalized routine to get explicit integrals for given operator Oper 
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!>
!> Generalized routine to get explicit integrals for given operator Oper 
!>
!>     'Coulomb'       1/r_12
!>     'Overlap'       dirac_delta(r_1-r_2)
!>     'Kinetic'       -1/2 nabla
!>
!> and with given AOs according to AO1,AO2,AO3 and AO4:
!>
!>     'Regular'       Regular basis
!>     'DF-Aux'        Auxiliary basis for density-fitting
!>     AOEmpty         Empty, used for two and three-center integrals
!>
!> AO1 and AO2 belong to electron 1, and AO3 and AO4 to electron 2
!>
!> First four dimensions of Integrals should correspond to the dimensions of 
!> AO1-4, the fifth dimension is an added index allowing for example for 
!> different derivative integrals.
!>
SUBROUTINE ls_getIntegrals1(AO1,AO2,AO3,AO4,Oper,Spec,intType,geoOrder,SETTING,LUPRI,LUERR)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,Spec,intType,geoOrder
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,nMAT,nAtoms,iMol
Logical              :: Classical,buildGab,doscreen,primScreen,batchOnly
Real(realk),pointer  :: derivativeIntegrals(:,:,:,:,:)
integer :: ndim2(5),nAtoms1,nAtoms2,nAtoms3,nAtoms4,I,J,ndim25,CMorder
logical :: sameAllFrag,MBIE_SCREENINT,dograd,FullAlphaCD

FullAlphaCD = setting%output%FullAlphaCD
setting%output%FullAlphaCD = .FALSE.
ndim2 = setting%output%ndim 
IF(Setting%scheme%cs_int.OR.Setting%scheme%ps_int)THEN
   MBIE_SCREENINT = Setting%scheme%MBIE_SCREEN.AND.&
        & (Oper.EQ.CoulombOperator.OR.Oper.EQ.ErfcOperator.OR.Oper.EQ.CAMOperator)
   CALL ls_getScreenIntegrals(AO1,AO2,Oper,Setting%scheme%cs_int,&
        & Setting%scheme%ps_int,MBIE_SCREENINT,SETTING,LUPRI,LUERR)
   RETURN
ENDIF

buildGab = .FALSE.
CALL init_integral_input(INT_INPUT,SETTING)

INT_INPUT%operator = Oper

IF(INT_INPUT%DO_PROP)CALL SET_PROPINFO1(Oper,INT_INPUT)
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%NonClassical_SCREEN = .FALSE.
IF (INT_INPUT%operator .EQ. MulmomOperator) THEN
   Int_input%hermiteEcoeff = .TRUE.
ENDIF

! Specifies if integrals are limited to one ore more AO-batches
batchOnly = (setting%batchindex(1).NE. 0).OR.(setting%batchindex(2).NE. 0) &
     &      .OR.(setting%batchindex(3).NE. 0).OR.(setting%batchindex(4).NE. 0)

! Screen only when calculationg three- and four-center integrals
! Currently turned off when calculating integrals only for a specific batch
doscreen = ((((Oper.EQ.CoulombOperator).OR.(Oper.EQ.NucleiOperator)).OR.((Oper.EQ.NucpotOperator).OR.(Oper.EQ.ErfcOperator)))&
     & .OR.(Oper.EQ.CAMOperator)).OR.(Oper.EQ.GGemOperator).OR.(Oper.EQ.GGemCouOperator).OR.(Oper.EQ.GGemGrdOperator) &
     & .OR. (Oper.EQ.GGemQuaOperator) &
     & .AND.(SETTING%SCHEME%CS_SCREEN.OR.SETTING%SCHEME%PS_SCREEN)

IF(INT_INPUT%DO_PROP)doscreen=.FALSE.

IF (doscreen) THEN
  Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
  Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
  IF(batchonly)THEN
     IF(associated(setting%LST_GAB_LHS).AND.associated(setting%LST_GAB_RHS))THEN
        IF(associated(setting%LST_GAB_LHS%maxprimgab).AND.associated(setting%LST_GAB_RHS%maxprimgab))THEN
           !the correct gab matrices have been set so we can activate primitve screening. 
        ELSE
           Int_input%PS_SCREEN = .FALSE. 
!           print*,'turn of Int_input%PS_SCREEN maxgab'
        ENDIF
     ELSE
        Int_input%CS_SCREEN = .FALSE. 
        Int_input%PS_SCREEN = .FALSE. 
!        print*,'turn of Int_input%PS_SCREEN and Int_input%CS_SCREEN gab'
     ENDIF
  ENDIF
  IF(Oper .EQ. CoulombOperator)Int_input%MBIE_SCREEN = SETTING%SCHEME%MBIE_SCREEN
  IF(Oper .EQ. ErfcOperator)Int_input%MBIE_SCREEN = SETTING%SCHEME%MBIE_SCREEN
  IF(Oper .EQ. CAMOperator)Int_input%MBIE_SCREEN = SETTING%SCHEME%MBIE_SCREEN
  CALL AttachScreeningMatricesToInput(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
ENDIF

IF(Oper .EQ. NucpotOperator) THEN

  INT_INPUT%AddToIntegral = .TRUE.
  
  Classical = SETTING%SCHEME%FMM .AND. (.NOT. SETTING%SCHEME%MM_NO_ONE)
  IF (Classical) THEN
    Int_input%NonClassical_SCREEN = SETTING%SCHEME%FMM
    Int_input%DO_FMM              = SETTING%SCHEME%FMM
    INT_INPUT%MM_TLMAX            = SETTING%SCHEME%MM_TLMAX
    INT_INPUT%MM_LMAX             = SETTING%SCHEME%MM_LMAX
    INT_INPUT%MM_SCREENTHR        = SETTING%SCHEME%MM_SCREEN*SETTING%SCHEME%intTHRESHOLD
  ENDIF
ENDIF

CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

IF (Oper.EQ.CarmomOperator) THEN
   INT_INPUT%nCartesianMomentComp = ndim2(5)   
   INT_INPUT%CMorder = setting%scheme%CMorder
   INT_INPUT%CMimat = setting%scheme%CMimat
   IF(INT_INPUT%CMimat.NE.0)THEN
      INT_INPUT%nCartesianMomentComp = (INT_INPUT%CMorder+1)*&
           & (INT_INPUT%CMorder+2)*(INT_INPUT%CMorder+3)/6
   ENDIF
ENDIF

IF(INT_INPUT%DO_PROP)CALL SET_PROPINFO2(setting%scheme%propoper,INT_INPUT)

call set_input_from_spec(INT_INPUT,SPEC,AO1,AO2,AO3,AO4,Oper,lupri,dograd,.TRUE.,geoOrder)
setting%output%dograd = dograd

IF (Spec.EQ.EcontribSpec) THEN
  DO I=1,4
     IF(ndim2(I).NE.1)THEN
       print*,'ndim2:',ndim2(1),ndim2(2),ndim2(3),ndim2(4)
       call lsquit('Error dim mismatch in Econtrib in ls_getIntegrals1',lupri)
     ENDIF
  ENDDO
ENDIF

CALL ls_setDensityDimensions(INT_INPUT,SETTING,lupri)
setting%output%FullAlphaCD = FullAlphaCD

call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)

IF(FullAlphaCD)CAll mem_dealloc(setting%output%postprocess)

CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
IF(doscreen) CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)
setting%output%ndim = ndim2

END SUBROUTINE ls_getIntegrals1

!> \brief Adds a sub block to a full tensor
!> \author S. Reine
!> \date 2010
!> \param integrals the full tensir
!> \param SubBlock the subblock to be added
!> \param sameAOsLHS if both LHS AOs are the same
!> \param sameFragmentLHS if both LHS fragments are the same
!> \param sameAOsRHS if both RHS AOs are the same
!> \param sameFragmentRHS if both RHS fragments are the same
!> \param nbast1 number of basis functions for the 1. center
!> \param nbast2 number of basis functions for the 2. center
!> \param nbast3 number of basis functions for the 3. center
!> \param nbast4 number of basis functions for the 4. center
!> \param start1 the start index for the 1. dimension 
!> \param start2 the start index for the 2. dimension 
!> \param start3 the start index for the 3. dimension 
!> \param start4 the start index for the 4. dimension 
!> \param end1 the end index for the 1. dimension 
!> \param end2 the end index for the 2. dimension 
!> \param end3 the end index for the 3. dimension 
!> \param end4 the end index for the 4. dimension 
!> \param ndmat the number of matrices (the fifth dimension)
SUBROUTINE addSubBlockFull(integrals,SubBlock,sameAOsLHS,sameFragmentLHS,&
     &                  sameAOsRHS,sameFragmentRHS,nbast1,nbast2,nbast3,nbast4,&
     &                  start1,start2,start3,start4,end1,end2,end3,end4,nMat)
implicit none
integer             :: nbast1,nbast2,nbast3,nbast4,nMat
integer             :: start1,start2,start3,start4
integer             :: end1,end2,end3,end4
Real(realk),pointer :: integrals(:,:,:,:,:)
Real(realk),pointer :: SubBlock(:,:,:,:,:)
Real(realk),pointer :: SubBlockPerm(:,:,:,:,:)
Logical             :: sameAOsLHS,sameAOsRHS,sameFragmentLHS,sameFragmentRHS
!
Integer :: I1,I2,I3,I4,s2,s4,orb1,orb2,orb3,orb4,iMat
Logical :: permuteLHS,permuteRHS

permuteLHS = sameAOsLHS.AND..NOT.sameFragmentLHS
permuteRHS = sameAOsRHS.AND..NOT.sameFragmentRHS
! Regular block
!integrals(start1:end1,start2:end2,start3:end3,start4:end4,1:nMat) = &
!     &  integrals(start1:end1,start2:end2,start3:end3,start4:end4,1:nMat) &
!     &+ SubBlock(:,:,:,:,1:nMat)
DO iMat=1,nMat
 DO I4=1,nbast4
  orb4 = start4 + I4 - 1
  DO I3=1,nbast3
   orb3 = start3 + I3 -1
   DO I2=1,nbast2
    orb2 = start2 + I2 - 1
    DO I1=1,nbast1
     orb1 = start1 + I1 - 1
     integrals(orb1,orb2,orb3,orb4,iMat) = integrals(orb1,orb2,orb3,orb4,iMat)&
          & + SubBlock(I1,I2,I3,I4,iMat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

! Integral blocks using permutational symmetry
IF (permuteLHS) THEN
  call mem_alloc(SubBlockPerm,nbast2,nbast1,nbast3,nbast4,nMat)
  DO iMat=1,nMat
    DO I3=1,nbast3
      orb3 = start3 + I3 -1
      s4 = 1
!      IF (sameFragmentRHS) s4=I3
      DO I4=s4,nbast4
        orb4 = start4 + I4 - 1                     
        DO I1=1,nbast1
          orb1 = start1 + I1 - 1
          s2 = 1
!          IF (sameFragmentLHS) s2=I1 !this will never happen
          DO I2=s2,nbast2
            orb2 = start2 + I2 - 1
            integrals(orb2,orb1,orb3,orb4,iMat) = integrals(orb2,orb1,orb3,orb4,iMat) + SubBlock(I1,I2,I3,I4,iMat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  call mem_dealloc(SubBlockPerm)
ENDIF
IF (permuteRHS) THEN
  call mem_alloc(SubBlockPerm,nbast2,nbast1,nbast3,nbast4,nMat)
  DO iMat=1,nMat
    DO I3=1,nbast3
      orb3 = start3 + I3 - 1
      s4 = 1
!      IF (sameFragmentRHS) s4=I3
      DO I4=s4,nbast4
        orb4 = start4 + I4 - 1
        DO I1=1,nbast1
          orb1 = start1 + I1 - 1
          s2 = 1
!          IF (sameFragmentLHS) s2=I1
          DO I2=s2,nbast2
            orb2 = start2 + I2 - 1
            integrals(orb1,orb2,orb4,orb3,iMat) = integrals(orb1,orb2,orb4,orb3,iMat) + SubBlock(I1,I2,I3,I4,iMat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  call mem_dealloc(SubBlockPerm)
ENDIF
IF (permuteLHS.AND.permuteRHS) THEN
  call mem_alloc(SubBlockPerm,nbast2,nbast1,nbast3,nbast4,nMat)
  DO iMat=1,nMat
    DO I3=1,nbast3
      orb3 = start3 + I3 - 1
      s4 = 1
!      IF (sameFragmentRHS) s4=I3
      DO I4=s4,nbast4
        orb4 = start4 + I4 - 1 
        DO I1=1,nbast1
          orb1 = start1 + I1 - 1 
          s2 = 1
!          IF (sameFragmentLHS) s2=I1
          DO I2=s2,nbast2
            orb2 = start2 + I2 - 1
            integrals(orb2,orb1,orb4,orb3,iMat) = integrals(orb2,orb1,orb4,orb3,iMat) + SubBlock(I1,I2,I3,I4,iMat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  call mem_dealloc(SubBlockPerm)
ENDIF

END SUBROUTINE addSubBlockFull

!> \brief Calculate the classical part (in an FMM sense) of the electron nuclear attraction contribution to the fock matrix
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param ndmat the size of the fifth dimension - number of matrices
!> \param ngeoderivcomp the order of derivative 
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_electronNuclearClassic(AO1,AO2,AO3,AO4,Oper,intType,SETTING,LUPRI,LUERR)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,nMAT
Logical              :: Classical,ODscreen

CALL init_integral_input(INT_INPUT,SETTING)

INT_INPUT%operator = Oper
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%NonClassical_SCREEN = .FALSE.
  
Classical = SETTING%SCHEME%FMM .AND. (.NOT. SETTING%SCHEME%MM_NO_ONE)
IF (Classical) THEN
   Int_input%NonClassical_SCREEN = SETTING%SCHEME%FMM
   Int_input%DO_FMM              = SETTING%SCHEME%FMM
   INT_INPUT%MM_TLMAX            = SETTING%SCHEME%MM_TLMAX
   INT_INPUT%MM_LMAX             = SETTING%SCHEME%MM_LMAX
   INT_INPUT%MM_SCREENTHR        = SETTING%SCHEME%MM_SCREEN*SETTING%SCHEME%intTHRESHOLD
ENDIF

CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

ODscreen=.FALSE.
CALL ls_setIntegralInputDensities(INT_INPUT,SETTING,lupri,.true.,.true.,.false.,ODscreen)

CALL initIntegralOutputDims1(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,setting%Output%ndim(5))
IF (INT_INPUT%DO_FMM) THEN
   call electronNuclearClassic(setting%Output,setting%Output%ndim(1),setting%Output%ndim(2),INT_INPUT)
ENDIF
call ls_freeIntegralInputDensities(INT_INPUT)
!Call ls_freeDmatFromSetting(setting)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)

END SUBROUTINE ls_electronNuclearClassic

!> \brief Calculate the exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_get_exchange_mat(AO1,AO3,AO2,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
#ifdef VAR_MPI
integer                    :: ndim_full(5),itask,jtask,I,J
integer                    :: nbast1,nbast2,nbast3,nbast4
logical                    :: CS_SCREEN,PS_SCREEN,BOTHpartioning
#endif
type(lstensor),pointer :: dmat_lhs,dmat_rhs,Kmat
type(lstensor),pointer :: dmat_lhs_full,dmat_rhs_full,Kmat_full
type(lstensor),pointer :: gabCS_rhs,gabCS_lhs
type(lstensor),pointer :: gabCS_rhs_full,gabCS_lhs_full
integer :: iAO,natoms,numnodes
logical :: rhsCS_created,doscreen,UseMPI
logical :: lhs_created,rhs_created,lhsCS_created,PermuteResultTensor
#ifdef VAR_MPI
Logical                    :: SMasterWakeSlaves
integer(kind=ls_mpik)      :: Snode,SNumnodes,SComm
IF(.NOT.setting%scheme%doMPI)THEN
 call deactivateIntegralMPI(Setting,Snode,SNumnodes,SComm,SMasterWakeSlaves)
ENDIF
IF(matrix_type.EQ.mtype_pdmm)THEN
   IF (setting%node.EQ.infpar%master) THEN
      !change into full format on master so that it can be bcast
      call SetScalapackDmatToFull(setting,.TRUE.,.TRUE.)
   ENDIF
ENDIF
#endif
#ifdef VAR_SCALAPACK
IF(matrix_type.EQ.mtype_scalapack)THEN
   IF (setting%node.EQ.infpar%master) THEN
      !change into full format on master so that it can be bcast
      call SetScalapackDmatToFull(setting,.TRUE.,.TRUE.)
   ENDIF
ENDIF
#endif
CALL ls_setDefaultFragments(setting)
UseMPI = .FALSE. !this is done inside ThermiteDriver for Link!
CALL ls_create_lstensor_full(setting,'AC_TYPE',AO1,AO3,AO2,AO4,Oper,Spec,&
     & intType,Kmat_full,dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created,&
     & gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created,&
     & PermuteResultTensor,doscreen,lupri,luerr,.TRUE.,useMPI)
#ifdef VAR_MPI
! *******************************************************************************
! *  MPI Specific: Wake slaves and call this subroutine                         *
! *******************************************************************************
IF (setting%node.EQ.infpar%master) THEN
   natoms = MAX(setting%molecule(1)%p%nAtoms,setting%molecule(2)%p%nAtoms,&
        &setting%molecule(3)%p%nAtoms,setting%molecule(4)%p%nAtoms)
   IF(setting%scheme%MasterWakeSlaves.AND.natoms.GT.1.AND.setting%scheme%LINK)THEN
      call ls_mpibcast(LSLINK,infpar%master,setting%comm)
      call lsmpi_link_masterToSlave(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
   ELSE
      !do not wake slaves for atomic calculation
      numnodes = setting%numnodes
      setting%numnodes = 1
   ENDIF
ENDIF
#endif
dmat_lhs => dmat_lhs_full
dmat_rhs => dmat_rhs_full
Kmat => Kmat_full
gabCS_rhs => gabCS_rhs_full
gabCS_lhs => gabCS_lhs_full
CALL ls_attach_lstensors_to_setting(setting,Kmat,dmat_lhs,dmat_rhs,lhs_created,rhs_created)
CALL ls_attach_gab_to_setting(setting,gabCS_lhs,gabCS_rhs)   
Call ls_get_exchange_mat1(AO1,AO3,AO2,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
call ls_free_lstensors_from_setting(setting,lupri)
CALL ls_free_gab_from_setting(setting,lupri)
#ifdef VAR_MPI
! *******************************************************************************
! *  MPI Specific: Reduction off the Exchange Matrix                            *
! *******************************************************************************
IF (setting%node.EQ.infpar%master) THEN
   IF(natoms.GT.1.AND.setting%scheme%LINK)then
      !only do reduction if the master woke up the slaves
      call lsmpi_barrier(setting%comm)
      setting%output%resultTensor => kmat_full
      CALL lsmpi_lstensor_reduction(setting%output%resultTensor,infpar%master,setting%node,setting%comm)
   ELSE
      setting%numnodes = numnodes
   ENDIF
ELSE
   !the slaves always does reduction
   call lsmpi_barrier(setting%comm)
   setting%output%resultTensor => kmat_full
   CALL lsmpi_lstensor_reduction(setting%output%resultTensor,infpar%master,setting%node,setting%comm)
ENDIF
#endif

IF (setting%node.NE.0) THEN
#ifdef VAR_MPI
   ! *******************************************************************************
   ! *  MPI Slaves Specific: Free mem and put slaves to sleep                      *
   ! *******************************************************************************
   !free mem and put slave to sleep
   call ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
   if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)
   call lstensor_free(setting%output%resultTensor)
   deallocate(setting%output%resultTensor)
   nullify(setting%output%resultTensor)
   IF(associated(setting%output%postprocess))THEN
      call mem_dealloc(setting%output%postprocess)
   ENDIF
   call typedef_free_setting(SETTING)
#endif
ELSE
   CALL ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
   if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)
ENDIF

#ifdef VAR_MPI
IF(.NOT.setting%scheme%doMPI)THEN
 call ReactivateIntegralMPI(Setting,Snode,SNumnodes,SComm,SMasterWakeSlaves)
ENDIF
#endif
END SUBROUTINE ls_get_exchange_mat

!!$!> \brief Calculate the exchange matrix
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!!$!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!!$!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!!$!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!!$!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!!$!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!!$!> \param intType the label for primitive or contracted calc
!!$!> \param SETTING Integral evalualtion settings
!!$!> \param LUPRI logical unit number of the Default output file
!!$!> \param LUERR logical unit number of the Default error file
!!$SUBROUTINE ls_get_exchange_mat(AO1,AO3,AO2,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
!!$implicit none
!!$integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
!!$Type(LSSETTING)      :: SETTING
!!$Integer              :: LUPRI,LUERR
!!$!
!!$#ifdef VAR_MPI
!!$integer                    :: ndim_full(5),itask,jtask,I,J
!!$integer                    :: nbast1,nbast2,nbast3,nbast4
!!$logical                    :: CS_SCREEN,PS_SCREEN,BOTHpartioning,SameAllFrag
!!$type(ls_task_manager)      :: tasks
!!$type(lsmpi_task_list)      :: task_list
!!$type(lstask),pointer       :: lhs_current,rhs_current
!!$#endif
!!$type(lstensor),pointer :: dmat_lhs,dmat_rhs,Kmat
!!$type(lstensor),pointer :: dmat_lhs_full,dmat_rhs_full,Kmat_full
!!$type(lstensor),pointer :: gabCS_rhs,gabCS_lhs
!!$type(lstensor),pointer :: gabCS_rhs_full,gabCS_lhs_full
!!$integer :: iAO
!!$logical :: rhsCS_created,doscreen
!!$logical :: lhs_created,rhs_created,lhsCS_created,PermuteResultTensor
!!$logical :: sameMolSave(4,4)
!!$
!!$#ifdef VAR_MPI
!!$!HACK Turn off permutational symmetries in case of MPI
!!$!ToDo Make MPI work with permutational symmetries
!!$sameMolSave = setting%sameMol
!!$setting%sameMol = .FALSE.
!!$DO iAO=1,4
!!$  setting%sameMol(iAO,iAO) = .TRUE.
!!$ENDDO
!!$#endif
!!$
!!$#ifdef VAR_SCALAPACK
!!$   IF(matrix_type.EQ.mtype_scalapack.or.matrix_type.EQ.mtype_pdmm)THEN
!!$      IF (setting%node.EQ.infpar%master) THEN
!!$         !change into full format on master 
!!$         call SetScalapackDmatToFull(setting,.TRUE.,.TRUE.)
!!$      ENDIF
!!$   ENDIF
!!$#endif
!!$
!!$CALL ls_setDefaultFragments(setting)
!!$CALL ls_create_lstensor_full(setting,'AC_TYPE',AO1,AO3,AO2,AO4,Oper,Spec,&
!!$     & intType,Kmat_full,dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created,&
!!$     & gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created,&
!!$     & PermuteResultTensor,doscreen,lupri,luerr,.TRUE.)
!!$
!!$#ifdef VAR_MPI
!!$ndim_full = setting%output%ndim
!!$CS_screen = setting%scheme%CS_SCREEN
!!$PS_screen = setting%scheme%PS_SCREEN
!!$! ***************************************************************************
!!$! *                                MPI Specific                             *
!!$! ***************************************************************************
!!$IF (setting%node.EQ.infpar%master) THEN
!!$   call ls_mpibcast(LSLINK,infpar%master,setting%comm)
!!$   call lsmpi_link_masterToSlave(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
!!$ENDIF
!!$
!!$!!!!!Symmetry for MPI is handled by the CS-screening matrices
!!$!!!!sameMolSave = Setting%sameMol
!!$!!!!Setting%sameMol = .FALSE.
!!$!!!!Setting%sameFrag = .FALSE.
!!$
!!$BOTHpartioning=.FALSE.
!!$CALL lsmpi_set_task_manager(tasks,'LinK',AO1,AO2,AO3,AO4,Spec,intType,setting,BOTHpartioning,lupri,luerr)
!!$CALL lsmpi_set_MPI_task_list(task_list,tasks,setting,lupri,luerr)
!!$
!!$DO itask=1,task_list%numMPItasks
!!$  I=task_list%taskPair(1,itask)
!!$  J=task_list%taskPair(2,itask)
!!$  lhs_current => task_list%LHS(I)%p
!!$  rhs_current => task_list%RHS(J)%p
!!$  CALL SetTaskFragments(SETTING,lhs_current,'LHS',tasks%lhs_aux,.FALSE.,.FALSE.,.TRUE.,lupri)
!!$  CALL SetTaskFragments(SETTING,rhs_current,'RHS',tasks%rhs_aux,.FALSE.,.FALSE.,.TRUE.,lupri)
!!$  CALL getTaskDimension(tasks%orbInfo,lhs_current,AO1,AO3,'LHS',Spec,nbast1,nbast3,lupri)
!!$  CALL getTaskDimension(tasks%orbInfo,rhs_current,AO2,AO4,'RHS',Spec,nbast2,nbast4,lupri)
!!$  call SET_SAMEALLFRAG(sameAllFrag,setting%sameFrag,setting%nAO)
!!$  CALL ls_create_lstensor_task(setting,kmat,'AC_TYPE',dmat_lhs,dmat_rhs,&
!!$       & kmat_full,dmat_lhs_full,dmat_rhs_full,&
!!$       & gabCS_rhs_full,gabCS_lhs_full,gabCS_rhs,gabCS_lhs,&
!!$       & rhsCS_created,lhsCS_created,&
!!$       & nbast1,nbast2,nbast3,nbast4,lhs_created,rhs_created,lhs_current,rhs_current,&
!!$       & tasks%sameAOsLHS,tasks%sameAOsRHS,tasks%sameODs,sameAllFrag,CS_screen,PS_SCREEN,doscreen)
!!$#else
!!$  ! ***************************************************************************
!!$  ! *                                  Serial                                 *
!!$  ! ***************************************************************************
!!$  dmat_lhs => dmat_lhs_full
!!$  dmat_rhs => dmat_rhs_full
!!$  Kmat => Kmat_full
!!$  gabCS_rhs => gabCS_rhs_full
!!$  gabCS_lhs => gabCS_lhs_full
!!$#endif
!!$  ! ***************************************************************************
!!$  ! *      Both for  Serial and MPI                                           *
!!$  ! ***************************************************************************
!!$  CALL ls_attach_lstensors_to_setting(setting,Kmat,dmat_lhs,dmat_rhs,lhs_created,rhs_created)
!!$  CALL ls_attach_gab_to_setting(setting,gabCS_lhs,gabCS_rhs)   
!!$  Call ls_get_exchange_mat1(AO1,AO3,AO2,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
!!$  call ls_free_lstensors_from_setting(setting,lupri)
!!$  CALL ls_free_gab_from_setting(setting,lupri)
!!$#ifdef VAR_MPI
!!$  call ls_extract_and_annihilate_lstensor_task(setting,kmat,'AC_TYPE',dmat_lhs,dmat_rhs,kmat_full,&
!!$         & dmat_lhs_full,dmat_rhs_full,nbast1,nbast2,nbast3,nbast4,lhs_created,&
!!$         & rhs_created,gabCS_rhs,gabCS_lhs,rhsCS_created,&
!!$         & lhsCS_created,lhs_current,rhs_current,&
!!$         & tasks%sameAOsLHS,tasks%sameAOsRHS,tasks%sameODs,SameAllFrag,CS_screen,PS_SCREEN,doscreen)
!!$ENDDO !task
!!$CALL lsmpi_free_MPI_task_list(task_list)
!!$call ls_free_tasks(tasks%lhs,lupri)
!!$call ls_free_tasks(tasks%rhs,lupri)
!!$DO iAO=1,4
!!$  call freeMolecularOrbitalInfo(tasks%orbInfo(iAO))
!!$ENDDO
!!$setting%sameMol = sameMolSave !Restores permutational symmetries
!!$call lsmpi_barrier(setting%comm)
!!$setting%output%resultTensor => kmat_full
!!$CALL lsmpi_lstensor_reduction(setting%output%resultTensor,infpar%master,setting%node,setting%comm)
!!$
!!$IF (setting%node.NE.infpar%master) THEN
!!$   !put slave to sleep
!!$   call ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
!!$   if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)
!!$   call lstensor_free(setting%output%resultTensor)
!!$   deallocate(setting%output%resultTensor)
!!$   nullify(setting%output%resultTensor)
!!$   IF(associated(setting%output%postprocess))THEN
!!$      call mem_dealloc(setting%output%postprocess)
!!$   ENDIF
!!$   call typedef_free_setting(SETTING)
!!$   RETURN
!!$ENDIF
!!$
!!$IF(PermuteResultTensor)THEN
!!$   call lstensor_full_symMat_from_triangularMat(setting%output%resultTensor)
!!$ENDIF
!!$
!!$Call freeDaltonFragments(SETTING)
!!$setting%output%ndim = ndim_full
!!$#endif
!!$
!!$CALL ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
!!$if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)
!!$
!!$END SUBROUTINE ls_get_exchange_mat

!debug version DO NOT USE 
subroutine ls_get_exchange_mat_serial(AO1,AO3,AO2,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
  implicit none
  integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
  Type(LSSETTING)      :: SETTING
  Integer              :: LUPRI,LUERR
  type(lstensor),pointer :: dmat_lhs,dmat_rhs,Kmat
  type(lstensor),pointer :: dmat_lhs_full,dmat_rhs_full,Kmat_full
  type(lstensor),pointer :: gabCS_rhs,gabCS_lhs
  type(lstensor),pointer :: gabCS_rhs_full,gabCS_lhs_full
  integer :: iAO
  logical :: rhsCS_created,doscreen
  logical :: lhs_created,rhs_created,lhsCS_created,PermuteResultTensor
  logical :: sameMolSave(4,4)
  CALL ls_setDefaultFragments(setting)
  CALL ls_create_lstensor_full(setting,'AC_TYPE',AO1,AO3,AO2,AO4,Oper,Spec,&
       & intType,Kmat_full,dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created,&
       & gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created,&
       & PermuteResultTensor,doscreen,lupri,luerr,.TRUE.,.TRUE.)
  dmat_lhs => dmat_lhs_full
  dmat_rhs => dmat_rhs_full
  Kmat => Kmat_full
  gabCS_rhs => gabCS_rhs_full
  gabCS_lhs => gabCS_lhs_full
  CALL ls_attach_lstensors_to_setting(setting,Kmat,dmat_lhs,dmat_rhs,lhs_created,rhs_created)
  CALL ls_attach_gab_to_setting(setting,gabCS_lhs,gabCS_rhs)   
  Call ls_get_exchange_mat1(AO1,AO3,AO2,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
  call ls_free_lstensors_from_setting(setting,lupri)
  CALL ls_free_gab_from_setting(setting,lupri)   
  CALL ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
  if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)
end subroutine ls_get_exchange_mat_serial

!> \brief Calculate the exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_get_exchange_mat1(AO1,AO3,AO2,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,idmat,idmat2,lupdmat
logical              :: dograd

CALL init_integral_input(INT_INPUT,SETTING)
CALL set_symmetry(INT_INPUT%iPQxyz,setting,lupri)
INT_INPUT%DO_FOCK = .TRUE.
INT_INPUT%exchangeFactor = SETTING%SCHEME%exchangeFactor
INT_INPUT%DO_EXCHANGE = (INT_INPUT%exchangeFactor.GT. 0.0E0_realk)&
     & .OR. (SETTING%SCHEME%CAMbeta .GT. 0.0E0_realk)

INT_INPUT%DO_LINK = SETTING%SCHEME%LINK
INT_INPUT%DO_DALINK = SETTING%SCHEME%DALINK
IF(.NOT.INT_INPUT%DO_EXCHANGE)CALL LSQUIT('ERROR',-1)

IF(INT_INPUT%DO_LINK)THEN 
   !The Input Matrix (normally density matrix) is symmetric
   !or it is split into symmetric and antisymmetric part
   !the anti symmetric part is still treated as symmetric when
   !construction the Fock matrix contribution but 
   !it is followed by anti symmetrization in Post Processing 
   INT_INPUT%DRHS_SYM=.TRUE.
   DO idmat = 1,setting%nDmatRHS
      idmat2 = idmat
      IF(matrix_type.EQ.mtype_unres_dense)then
         IF(MOD(idmat,2).EQ.0)THEN !2,4,6
            idmat2=idmat/2
         ELSE !1,3,5
            idmat2=idmat/2+1
         ENDIF
      ENDIF
      IF(setting%DsymRHS(idmat).EQ.1)THEN
         !Symmetric D => Symmetric K
         setting%output%postprocess(idmat2) = SymmetricPostprocess
      ELSEIF(setting%DsymRHS(idmat).EQ.2)THEN
         !AntiSymmetric D => AntiSymmetric K
         setting%output%postprocess(idmat2) = AntiSymmetricPostprocess
      ELSEIF(setting%DsymRHS(idmat).EQ.4)THEN
         !zero matrix 
         setting%output%postprocess(idmat2) = 0
      ELSE
         print*,'the code can handle nonsym densities but'
         print*,'from a perfomance perspective it is better'
         print*,'to divide matrix up into a Sym and antiSym Part'
         print*,'This should be done per default unless bypassed'
         print*,'by specifying that D is symmetric. '
         print*,'This Error statement can occur if Symmetry threshold'
         print*,'too high and the Symmetric Dmat is judged to be nonsymmetric.'
         print*,'The Matrix Judged to be non symmetric:'
         call lsquit('Exchange Called with nonsym Dmat',-1)
      ENDIF
      IF(Spec.EQ.MagDerivSpec)THEN
         setting%output%postprocess(idmat2) = 0         
      ENDIF
   ENDDO
   IF(AO1.NE.AO2)then
      !non square Exchange matrix
      INT_INPUT%DRHS_SYM=.FALSE.
      setting%output%postprocess=0
   ENDIF
   IF (setting%LHSdmat) THEN
      DO idmat = 1,setting%nDmatLHS
         idmat2 = idmat
         IF(matrix_type.EQ.mtype_unres_dense)then
            IF(MOD(idmat,2).EQ.0)THEN !2,4,6
               idmat2=idmat/2
            ELSE !1,3,5
               idmat2=idmat/2+1
            ENDIF
         ENDIF
         IF(setting%DsymLHS(idmat).EQ.1.OR.setting%DsymLHS(idmat).EQ.2)THEN
            INT_INPUT%DLHS_SYM=.TRUE.
         ENDIF
         IF(Spec.EQ.MagDerivSpec)THEN
            setting%output%postprocess(idmat2) = 0
            INT_INPUT%DLHS_SYM=.FALSE.
         ENDIF
      ENDDO
   ENDIF
   IF(Spec.EQ.GradientSpec)THEN
      setting%output%postprocess = 0
   ENDIF
   IF(Spec.EQ.EcontribSpec)THEN
      setting%output%postprocess = 0
   ENDIF
   IF(Spec.EQ.magderivEcontribSpec)THEN
      setting%output%postprocess = 0
   ENDIF
ENDIF
Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
Int_input%MBIE_SCREEN = SETTING%SCHEME%MBIE_SCREEN
CALL AttachScreeningMatricesToInput(INT_INPUT,AO1,AO3,AO2,AO4,Oper,SETTING,LUPRI,LUERR)          

IF(Oper .EQ. NucpotOperator) INT_INPUT%AddToIntegral = .TRUE.
INT_INPUT%operator = Oper
CALL SetInputAO(INT_INPUT,AO1,AO3,AO2,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
call set_input_from_spec(INT_INPUT,SPEC,AO1,AO3,AO2,AO4,Oper,lupri,dograd,.FALSE.)
CALL ls_setDensityDimensions(int_input,setting,lupri)
IF(.NOT.(INT_INPUT%sameODs.AND.INT_INPUT%sameLHSaos.AND.INT_INPUT%sameRHSaos))THEN
   INT_INPUT%DRHS_SYM=.FALSE.
   setting%output%postprocess = 0
ENDIF
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)

END SUBROUTINE ls_get_exchange_mat1

subroutine set_input_from_spec(INT_INPUT,SPEC,AO1,AO2,AO3,AO4,Oper,lupri,dograd,getint,geoDerivOrder)
implicit none
integer              :: AO1,AO2,AO3,AO4,lupri,Oper
integer,optional     :: geoDerivOrder
TYPE(INTEGRALINPUT)  :: INT_INPUT
integer,intent(in)   :: spec
logical,intent(out)  :: dograd
logical,intent(IN)   :: getint
!
integer :: I,geoOrder
logical :: emptyP,emptyQ,singleP,singleQ

geoOrder = 1
IF (present(geoDerivOrder)) geoOrder = geoDerivOrder
dograd  = .false.
emptyP  = (AO1.EQ.AOEmpty).AND.(AO2.EQ.AOEmpty)
emptyQ  = (AO3.EQ.AOEmpty).AND.(AO4.EQ.AOEmpty)
singleP = ((AO1.EQ.AOEmpty).OR.(AO2.EQ.AOEmpty)).AND..NOT.emptyP
singleQ = ((AO3.EQ.AOEmpty).OR.(AO4.EQ.AOEmpty)).AND..NOT.emptyQ

SELECT CASE(Spec)
CASE(RegularSpec)
   !Default undifferentiated case      
CASE(pChargeSpec)
  !Do nothing
CASE(GradientSpec)
   IF(( (AO1.NE.AOEmpty).AND.(AO2.EQ.AOEmpty) ).AND.((AO3.NE.AOEmpty).AND.(AO4.EQ.AOEmpty)))THEN
      Int_input%AC_center=.TRUE.
   ENDIF
   dograd = .true.
   !Molecular gradient
   INT_INPUT%DO_GRADIENT = .TRUE.

   IF (emptyP) CALL LSQUIT('Error in set_input_from_spec. GradientSpec and AO1=AO2="Empty"',-1)
   INT_INPUT%geoderOrderP  = 1
   INT_INPUT%geoderOrderQ  = 0
   INT_INPUT%sameODs       = .FALSE.
   INT_INPUT%NGEODERIVCOMP = getTotalGeoComp(1,.TRUE.,.FALSE.,singleP,singleQ,emptyP,emptyQ)
   INT_INPUT%AddToIntegral = .TRUE.
   INT_INPUT%GEODERIVORDER = max(INT_INPUT%geoderOrderP,INT_INPUT%geoderOrderQ)

   IF (Oper.EQ.KineticOperator) INT_INPUT%sameODs = .FALSE.
CASE(GeoDerivLHSSpec)
   IF (emptyP) CALL LSQUIT('Error in set_input_from_spec. GeoDerivLHSSpec and AO1=AO2="Empty"',-1)
   INT_INPUT%geoderOrderP  = 1
   INT_INPUT%geoderOrderQ  = 0
   INT_INPUT%sameODs       = .FALSE.
   INT_INPUT%NGEODERIVCOMP = getTotalGeoComp(1,.TRUE.,.FALSE.,singleP,singleQ,emptyP,emptyQ)
   INT_INPUT%AddToIntegral = .TRUE.
   INT_INPUT%GEODERIVORDER = max(INT_INPUT%geoderOrderP,INT_INPUT%geoderOrderQ)
CASE(GeoDerivRHSSpec)
   IF (emptyQ) CALL LSQUIT('Error in set_input_from_spec. GeoDerivRHSSpec and AO3=AO4="Empty"',-1)
   INT_INPUT%geoderOrderP  = 0
   INT_INPUT%geoderOrderQ  = 1
   INT_INPUT%sameODs       = .FALSE.
   INT_INPUT%NGEODERIVCOMP = getTotalGeoComp(1,.FALSE.,.TRUE.,singleP,singleQ,emptyP,emptyQ)
   INT_INPUT%AddToIntegral = .TRUE.
   INT_INPUT%GEODERIVORDER = max(INT_INPUT%geoderOrderP,INT_INPUT%geoderOrderQ)
CASE(GeoDerivSpec)
   IF (emptyP.AND.emptyQ) CALL LSQUIT('Error in set_input_from_spec. GeoDerivSpec and AO1=AO2=AO3=AO4="Empty"',-1)
   INT_INPUT%geoderOrderP  = geoOrder
   INT_INPUT%geoderOrderQ  = geoOrder
!  INT_INPUT%sameODs       = .FALSE.
   INT_INPUT%NGEODERIVCOMP = getTotalGeoComp(geoOrder,.TRUE.,.TRUE.,singleP,singleQ,emptyP,emptyQ)
   INT_INPUT%AddToIntegral = .TRUE.
   INT_INPUT%GEODERIVORDER = max(INT_INPUT%geoderOrderP,INT_INPUT%geoderOrderQ)
CASE(GeoDerivCoulombSpec)
   IF (emptyP) CALL LSQUIT('Error in set_input_from_spec. GeoDerivCoulombSpec and AO1=AO2="Empty"',-1)
   INT_INPUT%geoderOrderP  = 1
   INT_INPUT%geoderOrderQ  = 0
   INT_INPUT%sameODs       = .FALSE.
   INT_INPUT%NGEODERIVCOMP = getTotalGeoComp(1,.TRUE.,.FALSE.,singleP,singleQ,emptyP,emptyQ)
   INT_INPUT%AddToIntegral = .TRUE.
   INT_INPUT%GEODERIVORDER = max(INT_INPUT%geoderOrderP,INT_INPUT%geoderOrderQ)
   INT_INPUT%DO_FOCK = .TRUE.
   INT_INPUT%DO_Coulomb  = .TRUE.
CASE(MagDerivSpec) 
   IF(Oper.EQ.OverlapOperator)THEN
      INT_INPUT%LinComCarmomType=1
      !magnetic derivative integrals is a special case of integrals
      !which can be written as a linear combination of carmom integrals
      INT_INPUT%operator = CarmomOperator
      INT_INPUT%nCartesianMomentComp = 4
      INT_INPUT%CMorder = 1
      INT_INPUT%sameLHSAOs  = .FALSE.
   ELSE
      ! we do the derivative on the LHS 
      INT_INPUT%magderOrderP  = 1
      INT_INPUT%magderOrderQ  = 0
      INT_INPUT%NMAGDERIVCOMPP = 3
      INT_INPUT%NMAGDERIVCOMPQ = 1
!      INT_INPUT%sphericalEcoeff = .TRUE.!.FALSE.!
      INT_INPUT%HermiteEcoeff = .FALSE.
      INT_INPUT%MAGDERIVORDER = 1
   ENDIF
   INT_INPUT%doMagScreen=.TRUE.
   INT_INPUT%sameODs  = .FALSE.
   INT_INPUT%AddToIntegral = .TRUE.
CASE(MagGradSpec) 
   IF( ( (AO1.NE.AOEmpty).AND.(AO2.EQ.AOEmpty) ).AND.((AO3.NE.AOEmpty).AND.(AO4.EQ.AOEmpty)))THEN
      Int_input%AC_center=.TRUE.
   ENDIF
   dograd = .true.
   IF(Oper.EQ.OverlapOperator)THEN
      INT_INPUT%LinComCarmomType=1
      !magnetic derivative integrals is a special case of integrals
      !which can be written as a linear combination of carmom integrals
      INT_INPUT%operator = CarmomOperator
      INT_INPUT%nCartesianMomentComp = 4
      INT_INPUT%CMorder = 1
      INT_INPUT%sameLHSAOs  = .FALSE.
   ELSE
      ! we do the derivative on the LHS 
      INT_INPUT%magderOrderP  = 1
      INT_INPUT%magderOrderQ  = 0
      INT_INPUT%NMAGDERIVCOMPP = 3
      INT_INPUT%NMAGDERIVCOMPQ = 1
      INT_INPUT%HermiteEcoeff = .FALSE.
      INT_INPUT%MAGDERIVORDER = 1
   ENDIF

!  Special cases depending on which routine this is called from
   IF (getint) THEN  !ls_getIntegrals
     INT_INPUT%sameLHSaos  = .FALSE.
     INT_INPUT%sameRHSaos  = .FALSE.
   ELSE              !ls_get_coulomb_mat,ls_get_exchange_mat1
     INT_INPUT%doMagScreen=.TRUE.
   ENDIF

   INT_INPUT%sameODs  = .FALSE.
   INT_INPUT%AddToIntegral = .TRUE.
   INT_INPUT%DO_GRADIENT = .TRUE.
CASE(MagDerivLSpec) 
   dograd = .false.
   IF(Oper.EQ.OverlapOperator)THEN
      INT_INPUT%LinComCarmomType=2
      !magnetic derivative integrals is a special case of integrals
      !which can be written as a linear combination of carmom integrals
      INT_INPUT%operator = CarmomOperator
      INT_INPUT%nCartesianMomentComp = 4
      INT_INPUT%CMorder = 1
      INT_INPUT%sameLHSAOs  = .FALSE.
      INT_INPUT%sameODs  = .FALSE.
      INT_INPUT%AddToIntegral = .TRUE.
   ELSE
      ! we do the derivative on the LHS 
      INT_INPUT%magderOrderP  = 1
      INT_INPUT%magderOrderQ  = 0
      INT_INPUT%NMAGDERIVCOMPP = 3
      INT_INPUT%NMAGDERIVCOMPQ = 1
      INT_INPUT%HermiteEcoeff = .FALSE.
      INT_INPUT%MAGDERIVORDER = 1
      INT_INPUT%doMagScreen=.TRUE.
      INT_INPUT%sameODs  = .FALSE.
      INT_INPUT%AddToIntegral = .TRUE.
   ENDIF
CASE(MagDerivRSpec) 
   dograd = .false.
   IF(Oper.EQ.OverlapOperator)THEN
      INT_INPUT%LinComCarmomType=3
      !magnetic derivative integrals is a special case of integrals
      !which can be written as a linear combination of carmom integrals
      INT_INPUT%operator = CarmomOperator
      INT_INPUT%nCartesianMomentComp = 4
      INT_INPUT%CMorder = 1
      INT_INPUT%sameLHSAOs  = .FALSE.
      INT_INPUT%sameODs  = .FALSE.
      INT_INPUT%AddToIntegral = .TRUE.
   ELSE
      ! we do the derivative on the RHS 
      INT_INPUT%magderOrderP  = 0
      INT_INPUT%magderOrderQ  = 1
      INT_INPUT%NMAGDERIVCOMPP = 1
      INT_INPUT%NMAGDERIVCOMPQ = 3
      INT_INPUT%HermiteEcoeff = .FALSE.
      INT_INPUT%MAGDERIVORDER = 1
      INT_INPUT%doMagScreen=.TRUE.
      INT_INPUT%sameODs  = .FALSE.
      INT_INPUT%sameRHSaos  = .FALSE.
      INT_INPUT%AddToIntegral = .TRUE.
!     Special cases depending on which routine this is called from
      IF (getint) THEN ! From ls_getIntegrals
        INT_INPUT%do_passes = .FALSE. 
      ELSE             ! From ls_get_Coulomb_mat, ls_get_exchange_mat
        INT_INPUT%sameLHSaos  = .FALSE.
        INT_INPUT%do_passes = .FALSE. 
      ENDIF
   ENDIF
CASE (EcontribSpec)
   INT_INPUT%fullcontraction = .TRUE.
CASE (magderivEcontribSpec)
   INT_INPUT%fullcontraction = .TRUE.
   ! we do the derivative on the LHS 
   INT_INPUT%magderOrderP  = 1
   INT_INPUT%magderOrderQ  = 0
   INT_INPUT%NMAGDERIVCOMPP = 3
   INT_INPUT%NMAGDERIVCOMPQ = 1
   INT_INPUT%HermiteEcoeff = .FALSE.
   INT_INPUT%MAGDERIVORDER = 1
   INT_INPUT%doMagScreen=.TRUE.
   INT_INPUT%sameODs  = .FALSE. !false since only derivative on LHS 
   INT_INPUT%AddToIntegral = .TRUE.
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Error: Wrong case in set_input_from_spec =',Spec
   CALL LSQUIT('Wrong case in set_input_from_spec',lupri)
END SELECT

IF (INT_INPUT%geoderOrderP.EQ.1) THEN
  IF (singleP) THEN
    INT_INPUT%NGEODERIVCOMPP = 3
  ELSE
    INT_INPUT%NGEODERIVCOMPP = 6
  ENDIF
ENDIF

IF (INT_INPUT%geoderOrderQ.EQ.1) THEN
  IF (singleQ) THEN
    INT_INPUT%NGEODERIVCOMPQ = 3
  ELSE
    INT_INPUT%NGEODERIVCOMPQ = 6
  ENDIF
ENDIF

end subroutine set_input_from_spec

!> \brief Calculate the coulomb matrix
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_get_coulomb_mat(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,intType,spec
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
Logical             :: FRAGMENT_SAVE,dograd,ODscreen
integer :: ndim2(5)

CALL ls_setDefaultFragments(setting)
ndim2 = setting%output%ndim 

!CALL ls_setDfull(setting)

FRAGMENT_SAVE = SETTING%SCHEME%FRAGMENT
SETTING%SCHEME%FRAGMENT = .FALSE.

CALL init_integral_input(INT_INPUT,SETTING)

INT_INPUT%DO_FOCK = .TRUE.
INT_INPUT%DO_Coulomb  = .TRUE.
INT_INPUT%DO_DACOULOMB = SETTING%SCHEME%DACOULOMB

Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
Int_input%MBIE_SCREEN = SETTING%SCHEME%MBIE_SCREEN
CALL AttachScreeningMatricesToInput(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)           
setting%output%ndim = ndim2
INT_INPUT%operator = Oper

CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

call set_input_from_spec(INT_INPUT,SPEC,AO1,AO2,AO3,AO4,Oper,lupri,dograd,.FALSE.)

nullify(setting%output%resultTensor)
allocate(setting%output%resultTensor)
IF(Spec.EQ.EcontribSpec.OR.Spec.EQ.magderivEcontribSpec)THEN
   call init_lstensor_1dim(setting%output%resultTensor,ndim2(5),lupri)
ELSE
   call init_lstensor_5dim(setting%output%resultTensor,Int_Input%AO(1)%p,Int_Input%AO(2)%p,&
        & Int_Input%AO(3)%p,Int_Input%AO(4)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),&
     & .TRUE.,.TRUE.,.FALSE.,.FALSE.,Int_input%OD_SCREEN,.FALSE.,lupri)
ENDIF
ODscreen=.TRUE.
CALL ls_setintegralinputdensities(INT_INPUT,SETTING,lupri,.TRUE.,.TRUE.,.FALSE.,ODscreen)
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
call ls_freeIntegralInputDensities(INT_INPUT)
!Call ls_freeDmatFromSetting(setting)
CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)

SETTING%SCHEME%FRAGMENT = FRAGMENT_SAVE
!CALL ls_freeDfull(setting)
setting%output%ndim = ndim2

END SUBROUTINE ls_get_coulomb_mat

!> \brief Calculate the coulomb and exchange matrix
!> \author T.Kjaergaard and S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_get_coulomb_and_exchange_mat(AO1,AO2,AO3,AO4,Oper,intType,SETTING,LUPRI,LUERR)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: ndmat,nAObuilds,I
Logical             :: FRAGMENT_SAVE,ODscreen
integer :: ndim2(5)

CALL ls_setDefaultFragments(setting)
ndim2 = setting%output%ndim 

FRAGMENT_SAVE = SETTING%SCHEME%FRAGMENT
SETTING%SCHEME%FRAGMENT = .FALSE.
CALL init_integral_input(INT_INPUT,SETTING)
INT_INPUT%DO_FOCK = .TRUE.
INT_INPUT%DO_Coulomb  = .TRUE.
INT_INPUT%exchangeFactor = SETTING%SCHEME%exchangeFactor
INT_INPUT%DO_EXCHANGE = (INT_INPUT%exchangeFactor.GT. 0.0E0_realk)&
     & .OR. (SETTING%SCHEME%CAMbeta .GT. 0.0E0_realk)

Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
Int_input%MBIE_SCREEN = SETTING%SCHEME%MBIE_SCREEN
CALL AttachScreeningMatricesToInput(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)           
setting%output%ndim = ndim2
INT_INPUT%operator = Oper
CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

CALL initIntegralOutputDims1(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,SETTING%ndmatRHS)
nullify(setting%output%resultTensor)
allocate(setting%output%resultTensor)
! OD_SCREENING FALSE BECAUSE K(b,d) = (ab|cd)*D(ac) it would test an ODscreen between b and d which is not in same OD
ODscreen = .FALSE.
call init_lstensor_5dim(setting%output%resultTensor,Int_Input%AO(1)%p,Int_Input%AO(2)%p,&
     & Int_Input%AO(3)%p,Int_Input%AO(4)%p,ndim2(1),ndim2(2),ndim2(3),ndim2(4),ndim2(5),&
     & .TRUE.,.TRUE.,.FALSE.,.FALSE.,ODscreen,ODscreen,lupri)
CALL ls_setintegralinputdensities(INT_INPUT,SETTING,lupri,.TRUE.,.TRUE.,.FALSE.,ODscreen)
INT_INPUT%exchangeFactor = -0.5E0_realk*SETTING%SCHEME%exchangeFactor
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
call ls_freeIntegralInputDensities(INT_INPUT)
!Call ls_freeDmatFromSetting(setting)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)

SETTING%SCHEME%FRAGMENT = FRAGMENT_SAVE
setting%output%ndim = ndim2

END SUBROUTINE ls_get_coulomb_and_exchange_mat

!> \brief Calculate the coulomb matrix using the jengine method
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!>
!> and with given AOs according to AO1,AO2,AO3 and AO4:
!>      Name    int parameter 
!>     Regular      AORegular       Regular basis
!>     DF-Aux       AOdfAux         Auxiliary basis for density-fitting
!>     DF-CABS'     AOdfCABS        Complementary Auxiliary basis for F12
!>     DF-JK'       AOdfJK          Density-fitting basis set for Fock matrix for F12
!>     ADMM         AOadmm          Auxiliary Density matrix method basis set
!>     VALENCE      AOVAL           Regular Level 2  or Valence basis 
!>     Empty        AOEmpty         Empty, used for two and three-center integrals
SUBROUTINE ls_jengine(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
Integer              :: LUPRI,LUERR
Type(LSSETTING)      :: SETTING
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
!
TYPE(INTEGRALINPUT)  :: INT_INPUT

Integer              :: nbast1,nbast2,nbast3,nbast4
Integer              :: ndim
Integer              :: I1,I2,I3,I4,I,J,nbast
Integer              :: idmatRHS
real(realk),parameter :: ONE = 1E0_realk
Logical             :: permuteLHS,permuteRHS,FRAGMENT_SAVE
Logical             :: sameFragmentLHS,sameFragmentRHS
logical             :: l1,l2,lhs_created,rhs_created
integer             :: n1,n2,n3,n4,n5,iMol
Real(realk)         :: ts,te
type(lstensor),pointer :: dmat_lhs,dmat_rhs,jmat
type(lstensor),pointer     :: dmat_lhs_full,dmat_rhs_full,jmat_full
type(lstensor),pointer :: gabCS_rhs,gabCS_lhs
type(lstensor),pointer :: gabCS_rhs_full,gabCS_lhs_full
!
logical :: IntegralTransformGC
Logical :: rhsCS_created,lhsCS_created,rhsPS_created,lhsPS_created
Logical                    :: PermuteResultTensor,doscreen
!
#ifdef VAR_MPI
Integer                    :: itask,nAtoms(4),ntasks_lhs,ntasks_rhs,iAO,atomDim
Real(realk),pointer        :: Integrals(:,:,:,:,:)
Real(realk),pointer        :: SubBlockInt(:,:,:,:,:)
Real(realk),pointer        :: SubBlock(:,:)
Logical                    :: sameAOsLHS,sameAOsRHS,sameODs,sameAllMOL,sameAllFRAG
Logical                    :: lhs_aux,rhs_aux,noA,noB
Real(realk),pointer        :: timeMatrix(:,:)
integer,pointer            :: frag1(:,:),frag2(:,:)
type(ls_task_manager)      :: tasks
type(lsmpi_task_list)      :: task_list
type(lstask),pointer       :: lhs_current,rhs_current
logical                    :: sameMolSave(4,4),BOTHpartioning
Logical                    :: saveCSscreen,savePSscreen,CS_screen,PS_screen
integer                    :: ndim_full(5),iatom,jatom,ilsao,iatom2,jatom2,node
real(realk)                :: t(8),t1,t2,t3,t4
real(realk)                :: part(2)
Logical                    :: SMasterWakeSlaves
integer(kind=ls_mpik)      :: Snode,SNumnodes,SComm
IF(.NOT.setting%scheme%doMPI)THEN
   call deactivateIntegralMPI(Setting,Snode,SNumnodes,SComm,SMasterWakeSlaves)
ENDIF
#endif
!type(matrix)               :: tmp
!CALL LSTIMER('START',TS,TE,6)
IF(setting%scheme%memdist.AND.matrix_type.EQ.mtype_scalapack)THEN
   !Tempory fix for gradients and stuff 
   CALL ls_jengine_memdist(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
ELSE
#ifdef VAR_MPI
   IF(matrix_type.EQ.mtype_pdmm)THEN
      IF (setting%node.EQ.infpar%master) THEN
         !change into full format on master 
         call SetScalapackDmatToFull(setting,.TRUE.,.TRUE.)
      ENDIF
   ENDIF
#endif
#ifdef VAR_SCALAPACK
   IF(matrix_type.EQ.mtype_scalapack)THEN
      IF (setting%node.EQ.infpar%master) THEN
         !change into full format on master 
         call SetScalapackDmatToFull(setting,.TRUE.,.TRUE.)
      ENDIF
   ENDIF
#endif

#ifdef VAR_MPI
CALL LS_GETTIM(t1,t2)
t = 0E0_realk
#endif

IntegralTransformGC = .FALSE.
CALL ls_setDefaultFragments(setting)
CALL ls_create_lstensor_full(setting,'AB_TYPE',AO1,AO2,AO3,AO4,Oper,Spec,intType,&
     & jmat_full,dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created,&
     & gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created,&
     & PermuteResultTensor,doscreen,lupri,luerr,.FALSE.,.TRUE.)

#ifdef VAR_MPI

ndim_full = setting%output%ndim 
CS_screen = setting%scheme%CS_SCREEN
PS_screen = setting%scheme%PS_SCREEN
! ***************************************************************************
! *                                MPI Specific                             *
! ***************************************************************************
IF (setting%scheme%MasterWakeSlaves.AND.setting%node.EQ.infpar%master) THEN
   call ls_mpibcast(LSJENGIN,infpar%master,setting%comm)
   call lsmpi_jengine_masterToSlave(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
ENDIF

!Symmetry for MPI is handled by the CS-screening matrices
sameMolSave = Setting%sameMol 
Setting%sameMol = .FALSE.
Setting%sameFrag = .FALSE.

CALL LS_GETTIM(t3,t4)
t(7) = t(7) + t3 - t1 !Node CPU time  += task CPU time
t(8) = t(8) + t4 - t2 !Node wall time += task wall time
CALL LS_GETTIM(t1,t2)

BOTHpartioning=.FALSE.
CALL lsmpi_set_task_manager(tasks,'JENGINE',AO1,AO2,AO3,AO4,Spec,intType,setting,BOTHpartioning,lupri,luerr)
CALL lsmpi_set_MPI_task_list(task_list,tasks,setting,lupri,luerr)

part(1) = 0E0_realk
part(2) = 0E0_realk
DO itask=1,task_list%numMPItasks
  I=task_list%taskPair(1,itask)
  J=task_list%taskPair(2,itask)
  lhs_current => task_list%LHS(I)%p
  rhs_current => task_list%RHS(J)%p
  CALL SetTaskFragments(SETTING,lhs_current,'LHS',tasks%lhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
  CALL SetTaskFragments(SETTING,rhs_current,'RHS',tasks%rhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
  CALL getTaskDimension(tasks%orbInfo,lhs_current,AO1,AO2,'LHS',Spec,nbast1,nbast2,lupri)
  CALL getTaskDimension(tasks%orbInfo,rhs_current,AO3,AO4,'RHS',Spec,nbast3,nbast4,lupri)
  call SET_SAMEALLFRAG(sameAllFrag,setting%sameFrag,setting%nAO)
  CALL ls_create_lstensor_task(setting,jmat,'AB_TYPE',dmat_lhs,dmat_rhs,&
       & jmat_full,dmat_lhs_full,dmat_rhs_full,&
       & gabCS_rhs_full,gabCS_lhs_full,gabCS_rhs,gabCS_lhs,&
       & rhsCS_created,lhsCS_created,&
       & nbast1,nbast2,nbast3,nbast4,lhs_created,rhs_created,lhs_current,rhs_current,&
       & tasks%sameAOsLHS,tasks%sameAOsRHS,tasks%sameODs,sameAllFrag,CS_screen,PS_SCREEN,doscreen)
  CALL LS_GETTIM(t3,t4)
  t(1) = t(1) + t3 - t1 !Node CPU time  += task CPU time
  t(2) = t(2) + t4 - t2 !Node wall time += task wall time
  CALL LS_GETTIM(t1,t2)
#else
  ! ***************************************************************************
  ! *                                  Serial                                 *
  ! ***************************************************************************
  ! Set pointers to the full lstensors
  dmat_lhs => dmat_lhs_full
  dmat_rhs => dmat_rhs_full
  jmat => jmat_full
  gabCS_rhs => gabCS_rhs_full
  gabCS_lhs => gabCS_lhs_full
#endif
  ! ***************************************************************************
  ! *      Both for  Serial and MPI                                           *
  ! ***************************************************************************
  CALL ls_attach_lstensors_to_setting(setting,jmat,dmat_lhs,dmat_rhs,lhs_created,rhs_created)
  CALL ls_attach_gab_to_setting(setting,gabCS_lhs,gabCS_rhs)
  !*** CALCULATE INTEGRALS: matrix sub-block for MPI ***
  CALL ls_jengine1(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
  call ls_free_lstensors_from_setting(setting,lupri)
  CALL ls_free_gab_from_setting(setting,lupri)
#ifdef VAR_MPI
  CALL LS_GETTIM(t3,t4)
  t(3) = t(3) + t3 - t1 !Node CPU time  += task CPU time
  t(4) = t(4) + t4 - t2 !Node wall time += task wall time
! task_list%taskCPU(itask)  = t3 - t1
! task_list%taskWall(itask) = t4 - t2
  CALL LS_GETTIM(t1,t2)
  part(1) = part(1) + lhs_current%task_part
  part(2) = part(2) + rhs_current%task_part
  ! ***************************************************************************
  ! *                                MPI Specific                              *
  ! ***************************************************************************
  call ls_extract_and_annihilate_lstensor_task(setting,jmat,'AB_TYPE',dmat_lhs,dmat_rhs,jmat_full,&
       & dmat_lhs_full,dmat_rhs_full,nbast1,nbast2,nbast3,nbast4,lhs_created,&
       & rhs_created,gabCS_rhs,gabCS_lhs,rhsCS_created,&
       & lhsCS_created,lhs_current,rhs_current,&
       & tasks%sameAOsLHS,tasks%sameAOsRHS,tasks%sameODs,SameAllFrag,CS_screen,PS_SCREEN,doscreen)
ENDDO !task
CALL lsmpi_free_MPI_task_list(task_list)
call ls_free_tasks(tasks%lhs,lupri)
call ls_free_tasks(tasks%rhs,lupri)
DO iAO=1,4
  call freeMolecularOrbitalInfo(tasks%orbInfo(iAO))
ENDDO
Setting%sameMol = sameMolSave !HACK END
Setting%sameFrag = sameMolSave
call lsmpi_barrier(setting%comm)
CALL LS_GETTIM(t1,t2)

setting%output%resultTensor => jmat_full
CALL lsmpi_lstensor_reduction(setting%output%resultTensor,infpar%master,setting%node,setting%comm)

CALL LS_GETTIM(t3,t4)
t(5) = t(5) + t3 - t1 !Node CPU time  += task CPU time
t(6) = t(6) + t4 - t2 !Node wall time += task wall time
CALL LS_GETTIM(t1,t2)
call lsmpi_time_harvest(setting%node,t(1),t(2),t(3),t(4),t(5),t(6),t(7),t(8),&
     & part(1),part(2),setting%numnodes,lupri,setting%scheme%intprint)

IF (setting%node.NE.infpar%master) THEN
   !put slave to sleep
   call ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
   if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)
   call lstensor_free(setting%output%resultTensor)
   deallocate(setting%output%resultTensor)
   nullify(setting%output%resultTensor)
   IF(associated(setting%output%postprocess))THEN
      call mem_dealloc(setting%output%postprocess)
   ENDIF
   call typedef_free_setting(SETTING)
ELSE

CALL LS_GETTIM(t3,t4)
t(5) = t3 - t1 !Node CPU time  += task CPU time
t(6) = t4 - t2 !Node wall time += task wall time
CALL LS_GETTIM(t1,t2)
!Permutational Symmetry
!In order to exploit permutational symmetry we set the lower triangular Dmatrix to zero and multiply the upper by 2. We set the lower triangular screening matrices to zero and we now copy the transpose of the upper triangular part to the lower triangular part. 

IF(PermuteResultTensor)THEN
   call lstensor_full_symMat_from_triangularMat(setting%output%resultTensor)
ENDIF

Call freeDaltonFragments(SETTING)
setting%output%ndim = ndim_full
#endif

CALL ls_free_lstensors(dmat_lhs_full,dmat_rhs_full,lhs_created,rhs_created)
if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)

#ifdef VAR_MPI
CALL LS_GETTIM(t3,t4)
t(5) = t3 - t1 !Node CPU time  += task CPU time
t(6) = t4 - t2 !Node wall time += task wall time
ENDIF
#endif
ENDIF !memdist_jengine
#ifdef VAR_MPI
IF(.NOT.setting%scheme%doMPI)THEN
   call ReactivateIntegralMPI(Setting,Snode,SNumnodes,SComm,SMasterWakeSlaves)
ENDIF
#endif
END SUBROUTINE ls_jengine

!> \brief Calculate the coulomb matrix using the jengine method
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!>
!> and with given AOs according to AO1,AO2,AO3 and AO4:
!>      Name    int parameter 
!>     Regular      AORegular       Regular basis
!>     DF-Aux       AOdfAux         Auxiliary basis for density-fitting
!>     DF-CABS'     AOdfCABS        Complementary Auxiliary basis for F12
!>     DF-JK'       AOdfJK          Density-fitting basis set for Fock matrix for F12
!>     ADMM         AOadmm          Auxiliary Density matrix method basis set
!>     VALENCE      AOVAL           Regular Level 2  or Valence basis 
!>     Empty        AOEmpty         Empty, used for two and three-center integrals
SUBROUTINE ls_jengine_memdist(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
Integer              :: LUPRI,LUERR
Type(LSSETTING)      :: SETTING
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(INTEGRALINPUT)  :: INT_INPUT2

Integer              :: nbast1,nbast2,nbast3,nbast4
Integer              :: ndim
Integer              :: I1,I2,I3,I4,I,J,nbast
Integer              :: idmatRHS
real(realk),parameter :: ONE = 1E0_realk
Logical             :: permuteLHS,permuteRHS,FRAGMENT_SAVE
Logical             :: sameFragmentLHS,sameFragmentRHS
logical             :: l1,l2,lhs_created,rhs_created
integer             :: n1,n2,n3,n4,n5,iMol
Real(realk)         :: ts,te
type(lstensor),pointer :: dmat_lhs,dmat_rhs,jmat,jmat_other
type(lstensor),pointer :: gabCS_rhs,gabCS_lhs,gabCS_lhs_other
type(lstensor),pointer :: gabCS_rhs_full,gabCS_lhs_full
!
logical :: IntegralTransformGC
Logical :: rhsCS_created,lhsCS_created,rhsPS_created,lhsPS_created
Logical                    :: PermuteResultTensor,doscreen
!
#ifdef VAR_MPI
Integer                    :: itask,nAtoms(4),ntasks_lhs,ntasks_rhs,iAO,atomDim
Logical                    :: sameAOsLHS,sameAOsRHS,sameODs,sameAllMOL,sameAllFRAG
Logical                    :: lhs_aux,rhs_aux,noA,noB,LHSpartioning,RHSpartioning
logical                    :: sameMolSave(4,4),ForceRHSsymDMAT,ForceLHSsymDMAT,LOCALINFOSETUP
Logical                    :: saveCSscreen,savePSscreen,CS_screen,PS_screen,BOTHpartioning
integer                    :: ndim_full(5),iatom,jatom,ilsao,iatom2,jatom2,node
type(ls_task_manager)      :: tasks
type(lsmpi_task_list)      :: task_list
type(lstask),pointer       :: lhs_current,rhs_current
real(realk)                :: t(8),t1,t2,t3,t4,part(2)
!THIS SHOULD BE DELETED WHEN SIMEN HAVE CLEANEDUP
integer :: natom1,natom2,Dummyatomlist1(1),Dummyatomlist2(1),idmat
integer,pointer :: atoms1(:),atoms2(:)
type(matrix) :: DmatrixRHS,DmatrixLHS
integer(kind=long) :: nsizeFULL,nsizeLOCAL,nbuffer34
logical,pointer :: MessageRecieved(:)
real(realk),pointer :: buffer34(:)
integer :: nbuffer34T
CALL LS_GETTIM(t1,t2)
t = 0E0_realk

!we have 3 possibilities 
!1. J_AB = (AB|CD) D_CD        Partitioning on the LHS and RHS (all memdist matrices and isend/irecv)
!2. G_ALPHA = (ALPHA|CD) D_CD  Only partitioning on the RHS (memdist DmatRHS)
!3. J_AB = (AB|ALPHA) C_ALPHA  Only partitioning on the LHS (memdist Jmat and DmatLHS)
! so we set up the RHSpartitioning and LHSpartitioning
!print*,'AO1,AO2,AO3,AO4',AO1,AO2,AO3,AO4
! THIS IS THE WAY IT SHOULD BE 
!!$IF((AO1.EQ.AORdefault).AND.(AO2.EQ.AORdefault))THEN
!!$   LHSpartioning = .TRUE. 
!!$ELSEIF((AO1.EQ.AODFdefault).AND.(AO2.EQ.AOempty))THEN
!!$   LHSpartioning = .FALSE.
!!$ENDIF
!!$IF((AO3.EQ.AORdefault).AND.(AO4.EQ.AORdefault))THEN
!!$   RHSpartioning = .TRUE.
!!$ELSEIF((AO3.EQ.AODFdefault).AND.(AO4.EQ.AOempty))THEN
!!$   RHSpartioning = .FALSE.
!!$ENDIF
!FOR NOW WE FOLLOW the lsmpi.f90 file
BOTHpartioning = .FALSE.
IF((AO3.EQ.AORdefault).AND.(AO4.EQ.AORdefault))THEN
!   IF((AO1.EQ.AORdefault).AND.(AO2.EQ.AORdefault))THEN
!      LHSpartioning = .TRUE.
!      RHSpartioning = .TRUE.
!      BOTHpartioning = .TRUE.
!   ELSE
      LHSpartioning = .FALSE.
      RHSpartioning = .TRUE.
!   ENDIF
ELSEIF((AO1.EQ.AORdefault).AND.(AO2.EQ.AORdefault))THEN
   LHSpartioning = .TRUE.
   RHSpartioning = .FALSE.
ELSE
   call lsquit('need testing',-1)
ENDIF
IF(setting%LHSdmat.AND..NOT.LHSpartioning)then
   call SetScalapackDmatToFull(setting,.TRUE.,.FALSE.)
ENDIF
IF(setting%RHSdmat.AND..NOT.RHSpartioning)then
   call SetScalapackDmatToFull(setting,.FALSE.,.TRUE.)
ENDIF
IntegralTransformGC = .FALSE.
CALL ls_setDefaultFragments(setting)
!setup full lstensors if requires and setup global info for memory distributed lstensors 
CALL ls_memdist_lstensor_SetupFullinfo(setting,'AB_TYPE',AO1,AO2,AO3,AO4,Oper,Spec,intType,&
     & jmat,jmat_other,dmat_lhs,dmat_rhs,lhs_created,rhs_created,&
     & gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created,&
     & PermuteResultTensor,doscreen,ForceRHSsymDMAT,ForceLHSsymDMAT,&
     & LHSpartioning,RHSpartioning,BOTHpartioning,lupri,luerr)
IF(lhs_created.AND.BOTHpartioning)call lsquit('not tested',-1)
ndim_full = setting%output%ndim 
CS_screen = setting%scheme%CS_SCREEN
PS_screen = setting%scheme%PS_SCREEN
! ***************************************************************************
! *                                MPI Specific                             *
! ***************************************************************************
IF (setting%node.EQ.infpar%master) THEN
   IF(RHSpartioning.AND.rhs_created)THEN
      IF(ForceRHSsymDMAT)THEN
         !build DmatrixRHS = D + D^T
         call mat_init(DmatrixRHS,setting%DmatRHS(1)%p%ncol,setting%DmatRHS(1)%p%nrow)
         call mat_trans(setting%DmatRHS(1)%p,DmatrixRHS)
         call mat_daxpy(1.E0_realk,setting%DmatRHS(1)%p,DmatrixRHS)
         call mat_scal_dia(0.5E0_realk,DmatrixRHS)
         call mat_setlowertriangular_zero(DmatrixRHS)
      ELSE
         call mat_init(DmatrixRHS,setting%DmatRHS(1)%p%nrow,setting%DmatRHS(1)%p%ncol)
         call mat_assign(DmatrixRHS,setting%DmatRHS(1)%p)
      ENDIF
!      print*,'SCALAPACK PRINT AFTER possible(setlowertriangular_zero)',ForceRHSsymDMAT
!      call mat_print(DmatrixRHS,1,DmatrixRHS%nrow,1,DmatrixRHS%ncol,6)
   ENDIF
   IF(LHSpartioning.AND.lhs_created)THEN
      IF(ForceLHSsymDMAT)THEN
         !build DmatrixLHS = D + D^T
         call mat_init(DmatrixLHS,setting%DmatLHS(1)%p%ncol,setting%DmatLHS(1)%p%nrow)
         call mat_trans(setting%DmatLHS(1)%p,DmatrixLHS)
         call mat_daxpy(1.E0_realk,setting%DmatLHS(1)%p,DmatrixLHS)
         call mat_scal_dia(0.5E0_realk,DmatrixLHS)
         call mat_setlowertriangular_zero(DmatrixLHS)
      ELSE
         call mat_init(DmatrixLHS,setting%DmatLHS(1)%p%nrow,setting%DmatLHS(1)%p%ncol)
         call mat_assign(DmatrixLHS,setting%DmatLHS(1)%p)
      ENDIF
!      print*,'SCALAPACK PRINT AFTER possible(setlowertriangular_zero)',ForceLHSsymDMAT
!      call mat_print(DmatrixLHS,1,DmatrixLHS%nrow,1,DmatrixLHS%ncol,6)
   ENDIF
   call ls_mpibcast(LSJENGIN,infpar%master,setting%comm)
   call lsmpi_jengine_masterToSlave(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
ENDIF

!Symmetry for MPI is handled by the CS-screening matrices
sameMolSave = Setting%sameMol 
Setting%sameMol = .FALSE.
Setting%sameFrag = .FALSE.

CALL LS_GETTIM(t3,t4)
t(7) = t(7) + t3 - t1 !Node CPU time  += task CPU time
t(8) = t(8) + t4 - t2 !Node wall time += task wall time
CALL LS_GETTIM(t1,t2)
CALL lsmpi_set_task_manager(tasks,'JENGINE',AO1,AO2,AO3,AO4,Spec,intType,setting,BOTHpartioning,lupri,luerr)
!CALL lsmpi_set_MPI_task_list(task_list,tasks,setting,lupri,luerr)
part(1) = 0E0_realk
part(2) = 0E0_realk

!=====================================================================================================
!
! Setup Local Dmat and/or Jmat 
!
!=====================================================================================================

IF(RHSpartioning)THEN   
   !setting mynum info
   ! CLEARLY THIS IS A HACK AND CRAPPY CODE HOPEFULLY SIMEN WILL CLEAN THIS UP WHEN 
   ! THE INFO ABOUT WHICH NODE HAVE WHICH ATOMS IS A LITTLE MORE EASILY AVALIBLE 
   !STEP 1 : Determine which nodes own which part of the lstensor 
   rhs_current => tasks%rhs%first
   IF(tasks%rhs%ntasks.GT.setting%numnodes)call lsquit('MPI MEMDIST RHS error',-1)
   DO J=1,tasks%rhs%ntasks
      ! I always set the mynum in the global info of the lstensor
      call set_single_atomspointers(rhs_current,natom1,natom2,atoms1,atoms2)
      do IATOM2=1,natom1
         IATOM = atoms1(IATOM2)
         do JATOM2=1,natom2
            JATOM = atoms2(JATOM2)
            ILSAO = dmat_rhs%G_INDEX(IATOM,JATOM)
            IF(ILSAO.NE.0) dmat_rhs%G_LSAO(ILSAO)%mynum = J-1
         enddo
      enddo
      IF(rhs_current%full_mol_row)call mem_dealloc(atoms1)
      IF(rhs_current%full_mol_col)call mem_dealloc(atoms2)
      rhs_current => rhs_current%next
   ENDDO !J
   IF(rhs_created)Dmat_rhs%LowerDiagZero = ForceRHSsymDMAT
   !STEP 2 : setup node local info for the memory distributed lstensor
   rhs_current => tasks%rhs%first
   LOCALINFOSETUP = .FALSE.
   DO J=1,tasks%rhs%ntasks
      IF(J-1.EQ.setting%node)THEN
         LOCALINFOSETUP = .TRUE.
         CALL SetTaskFragments(SETTING,rhs_current,'RHS',tasks%rhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
         CALL getTaskDimension(tasks%orbInfo,rhs_current,AO3,AO4,'RHS',Spec,nbast3,nbast4,lupri)
         ! I always set the mynum in the global info of the lstensor
         call set_single_atomspointers(rhs_current,natom1,natom2,atoms1,atoms2)
         !this task belongs to my so I allocate the block - Special if no partitioning
         call ls_memdist_lstensor_SetupLocalinfoDRHS(setting,AO1,AO2,AO3,AO4,&
              & intType,dmat_rhs,rhs_created,nAtom1,nAtom2,atoms1,atoms2,RHSpartioning,lupri,luerr)
         !setup Gab for the RHS
         IF(rhsCS_created)THEN
            nullify(gabCS_rhs)
            allocate(gabCS_rhs)
            call build_sublstensor_from_full_lstensor(gabCS_rhs,gabCS_rhs_full,natom1,natom2,1,1,&
                 & atoms1,atoms2,Dummyatomlist1,Dummyatomlist2,nbast3,nbast4,1,1,.TRUE.)
         ENDIF
         IF(rhs_current%full_mol_row)call mem_dealloc(atoms1)
         IF(rhs_current%full_mol_col)call mem_dealloc(atoms2)
      ENDIF
      rhs_current => rhs_current%next
   ENDDO !J
   IF(.NOT.LOCALINFOSETUP)THEN
      print*,'WARNING: no RHS tasks for mynum',setting%node
      IF(rhs_created)THEN
         call build_empty_sublstensor(dmat_rhs)
      ENDIF
      IF(rhsCS_created)THEN
         nullify(gabCS_rhs)
         allocate(gabCS_rhs)
         call alloc_build_empty_sublstensor(gabCS_rhs)
      ENDIF
   ENDIF
   if(dmat_rhs%ndim5.GT.1) call lsquit('error',-1)
   !STEP 3 : build memory distributed lstensor from scalapack formattet Density matrix
   IF(rhs_created)THEN
      IF (setting%node.EQ.infpar%master) THEN
         call memdist_lstensor_BuildFromScalapack(dmat_rhs,setting%comm,&
              & setting%node,setting%numnodes,DmatrixRHS)
      ELSE
         call memdist_lstensor_BuildFromScalapack(dmat_rhs,setting%comm,&
              & setting%node,setting%numnodes)
      ENDIF
   ENDIF
ELSE
   IF(rhsCS_created)THEN
      gabCS_rhs => gabCS_rhs_full
   ENDIF
ENDIF

IF(LHSpartioning)THEN
   !setting mynum info
   ! CLEARLY THIS IS A HACK AND CRAPPY CODE HOPEFULLY SIMEN WILL CLEAN THIS UP WHEN 
   ! THE INFO ABOUT WHICH NODE HAVE WHICH ATOMS IS A LITTLE MORE EASILY AVALIBLE 
   !STEP 1 : Determine which nodes own which part of the lstensor 
   lhs_current => tasks%lhs%first
   IF(tasks%lhs%ntasks.GT.setting%numnodes)call lsquit('MPI MEMDIST LHS error',-1)
   DO I=1,tasks%lhs%ntasks
      ! I always set the mynum in the global info of the lstensor
      call set_single_atomspointers(lhs_current,natom1,natom2,atoms1,atoms2)
      do IATOM2=1,natom1
         IATOM = atoms1(IATOM2)
         do JATOM2=1,natom2
            JATOM = atoms2(JATOM2)
            ILSAO = jmat%G_INDEX(IATOM,JATOM)
            IF(ILSAO.NE.0) jmat%G_LSAO(ILSAO)%mynum = I-1
         enddo
      enddo
      IF(lhs_current%full_mol_row)call mem_dealloc(atoms1)
      IF(lhs_current%full_mol_col)call mem_dealloc(atoms2)
      lhs_current => lhs_current%next
   ENDDO !I
   IF(lhs_created)Dmat_lhs%LowerDiagZero = ForceLHSsymDMAT
   !STEP 2 : setup node local info for the memory distributed lstensor
   lhs_current => tasks%lhs%first
   LOCALINFOSETUP = .FALSE.
   DO I=1,tasks%lhs%ntasks
      IF(I-1.EQ.setting%node)THEN
         LOCALINFOSETUP = .TRUE.
         CALL SetTaskFragments(SETTING,lhs_current,'LHS',tasks%lhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
         CALL getTaskDimension(tasks%orbInfo,lhs_current,AO3,AO4,'LHS',Spec,nbast3,nbast4,lupri)
         ! I always set the mynum in the global info of the lstensor
         call set_single_atomspointers(lhs_current,natom1,natom2,atoms1,atoms2)
         !LHS DMAT
         call ls_memdist_lstensor_SetupLocalinfoDLHS(setting,AO1,AO2,AO3,AO4,&
              & intType,dmat_lhs,lhs_created,nAtom1,nAtom2,atoms1,atoms2,LHSpartioning,lupri,luerr)
         !The Coulomb Matrix
         call ls_memdist_lstensor_SetupLocalinfoRESLHS(setting,AO1,AO2,AO3,AO4,&
              & intType,jmat,nAtom1,nAtom2,atoms1,atoms2,1,2,LHSpartioning,lupri,luerr)
         !setup Gab for the LHS
         IF(lhsCS_created)THEN
            nullify(gabCS_lhs)
            allocate(gabCS_lhs)
            call build_sublstensor_from_full_lstensor(gabCS_lhs,gabCS_lhs_full,natom1,natom2,1,1,&
                 & atoms1,atoms2,Dummyatomlist1,Dummyatomlist2,nbast1,nbast2,1,1,.TRUE.)
         ENDIF
         IF(lhs_current%full_mol_row)call mem_dealloc(atoms1)
         IF(lhs_current%full_mol_col)call mem_dealloc(atoms2)
      ENDIF
      lhs_current => lhs_current%next
   ENDDO !I
   IF(.NOT.LOCALINFOSETUP)THEN
      print*,'WARNING: no LHS tasks for mynum',setting%node
      call build_empty_sublstensor(jmat)
      IF(lhs_created)THEN
         call build_empty_sublstensor(dmat_lhs)
      ENDIF
      IF(lhsCS_created)THEN
         nullify(gabCS_lhs)
         allocate(gabCS_lhs)
         call alloc_build_empty_sublstensor(gabCS_lhs)
      ENDIF
   ENDIF
   if(jmat%ndim5.GT.1) call lsquit('error',-1)
   !STEP 3 : build memory distributed lstensor from scalapack formattet Density matrix
   IF(lhs_created)THEN
      IF (setting%node.EQ.0) THEN
         call memdist_lstensor_BuildFromScalapack(dmat_lhs,setting%comm,&
              & setting%node,setting%numnodes,DmatrixLHS)
      ELSE
         call memdist_lstensor_BuildFromScalapack(dmat_lhs,setting%comm,&
              & setting%node,setting%numnodes)
      ENDIF
   ENDIF
ELSE
   IF(lhsCS_created)THEN
      gabCS_lhs => gabCS_lhs_full
   ENDIF
ENDIF
!=====================================================================================================
!
! Actual Calculation
!
!=====================================================================================================
IF(Bothpartioning)THEN
   call mem_alloc(MessageRecieved,setting%numnodes)
   DO I=1,setting%numnodes
      MessageRecieved(I)=.FALSE.
   ENDDO
   MessageRecieved(setting%node+1) = .TRUE.
ENDIF

!IF Bothpartioning then ALL nodes will start building LHS belonging to master and send to master 
!that is probably not a good Idea 
!We should probably change this to start from mynum and go around to mynum-1 or something
itask = -1
lhs_current => tasks%lhs%first
DO I=1,tasks%lhs%ntasks
   IF((.NOT.LHSpartioning).OR.I-1.EQ.setting%node)THEN      
      rhs_current => tasks%rhs%first
      DO J=1,tasks%rhs%ntasks
         IF((.NOT.RHSpartioning).OR.J-1.EQ.setting%node)THEN      
            CALL SetTaskFragments(SETTING,lhs_current,'LHS',tasks%lhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
            CALL SetTaskFragments(SETTING,rhs_current,'RHS',tasks%rhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
            CALL getTaskDimension(tasks%orbInfo,lhs_current,AO1,AO2,'LHS',Spec,nbast1,nbast2,lupri)
            CALL getTaskDimension(tasks%orbInfo,rhs_current,AO3,AO4,'RHS',Spec,nbast3,nbast4,lupri)
            call SET_SAMEALLFRAG(sameAllFrag,setting%sameFrag,setting%nAO)
            CALL LS_GETTIM(t3,t4)
            t(1) = t(1) + t3 - t1 !Node CPU time  += task CPU time
            t(2) = t(2) + t4 - t2 !Node wall time += task wall time
            CALL LS_GETTIM(t1,t2)
            CALL ls_attach_lstensors_to_setting(setting,jmat,dmat_lhs,dmat_rhs,lhs_created,rhs_created)
            CALL ls_attach_gab_to_setting(setting,gabCS_lhs,gabCS_rhs)
            !*** CALCULATE INTEGRALS: matrix sub-block for MPI ***
!            print*,'I setting%node',setting%node,'I,J',I,J,'tasks%lhs%ntasks,tasks%rhs%ntasks',tasks%lhs%ntasks,tasks%rhs%ntasks
            CALL ls_jengine1(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
            call ls_free_lstensors_from_setting(setting,lupri)
            CALL ls_free_gab_from_setting(setting,lupri)            
            CALL LS_GETTIM(t3,t4)
            t(3) = t(3) + t3 - t1 !Node CPU time  += task CPU time
            t(4) = t(4) + t4 - t2 !Node wall time += task wall time
            CALL LS_GETTIM(t1,t2)
            part(1) = part(1) + lhs_current%task_part
            part(2) = part(2) + rhs_current%task_part
         ENDIF
         rhs_current => rhs_current%next
      ENDDO !J
   ELSEIF(Bothpartioning.AND.I-1.NE.setting%node)THEN
      !node setting%node calculates a contribution that do not belong to itself
      !but needs to be sendt to another node
      rhs_current => tasks%rhs%first
      DO J=1,tasks%rhs%ntasks
         IF((.NOT.RHSpartioning).OR.J-1.EQ.setting%node)THEN      
            CALL SetTaskFragments(SETTING,lhs_current,'LHS',tasks%lhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
            CALL SetTaskFragments(SETTING,rhs_current,'RHS',tasks%rhs_aux,.FALSE.,.FALSE.,.FALSE.,lupri)
            CALL getTaskDimension(tasks%orbInfo,lhs_current,AO1,AO2,'LHS',Spec,nbast1,nbast2,lupri)
            CALL getTaskDimension(tasks%orbInfo,rhs_current,AO3,AO4,'RHS',Spec,nbast3,nbast4,lupri)
            call SET_SAMEALLFRAG(sameAllFrag,setting%sameFrag,setting%nAO)
            !still need to build tmp result_tensor and gabCS_lhs from gabCS_lhs_full
            !if we have LHS and RHS partitioning, so that node I build a part of the 
            !result tensor belonging to node J. 
            call set_single_atomspointers(lhs_current,natom1,natom2,atoms1,atoms2)
            call ls_memdist_lstensor_SetupLocalinfoRESLHS(setting,AO1,AO2,AO3,AO4,&
                 & intType,jmat_other,nAtom1,nAtom2,atoms1,atoms2,1,2,LHSpartioning,lupri,luerr)
            !setup Gab for the LHS
            IF(lhsCS_created)THEN
               nullify(gabCS_lhs_other)
               allocate(gabCS_lhs_other)
               call build_sublstensor_from_full_lstensor(gabCS_lhs_other,gabCS_lhs_full,natom1,natom2,1,1,&
                    & atoms1,atoms2,Dummyatomlist1,Dummyatomlist2,nbast1,nbast2,1,1,.TRUE.)
            ENDIF
            IF(lhs_current%full_mol_row)call mem_dealloc(atoms1)
            IF(lhs_current%full_mol_col)call mem_dealloc(atoms2)

            CALL LS_GETTIM(t3,t4)
            t(1) = t(1) + t3 - t1 !Node CPU time  += task CPU time
            t(2) = t(2) + t4 - t2 !Node wall time += task wall time
            CALL LS_GETTIM(t1,t2)
            CALL ls_attach_lstensors_to_setting(setting,jmat_other,dmat_lhs,dmat_rhs,lhs_created,rhs_created)
            CALL ls_attach_gab_to_setting(setting,gabCS_lhs_other,gabCS_rhs)
            !*** CALCULATE INTEGRALS: matrix sub-block for MPI ***
!            print*,'I setting%node',setting%node,'I,J',I,J,'tasks%lhs%ntasks,tasks%rhs%ntasks',tasks%lhs%ntasks,tasks%rhs%ntasks
            CALL ls_jengine1(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
            call ls_free_lstensors_from_setting(setting,lupri)
            CALL ls_free_gab_from_setting(setting,lupri)
            
            CALL LS_GETTIM(t3,t4)
            t(3) = t(3) + t3 - t1 !Node CPU time  += task CPU time
            t(4) = t(4) + t4 - t2 !Node wall time += task wall time
            CALL LS_GETTIM(t1,t2)
            part(1) = part(1) + lhs_current%task_part
            part(2) = part(2) + rhs_current%task_part

            !send to node I-1
            call lsmpi_isend_lstmemrealkbuf(jmat_other%lstmem_index,I-1,setting%node,setting%comm)!,request(I))
            ! check for incoming messages and if any messages add to jmat realk buffer 
            call lsmpi_probe_and_irecv_add_lstmemrealkbuf(jmat%lstmem_index,setting%node,&
                 & setting%comm,MessageRecieved,setting%numnodes)
            !free LOCAL INFO from lstensor
            call lstensor_local_free(jmat_other)
            IF(lhsCS_created)THEN
               call lstensor_free(gabCS_lhs_other)
               deallocate(gabCS_lhs_other)
               nullify(gabCS_lhs_other)
            ENDIF         
         ENDIF
         rhs_current => rhs_current%next
      ENDDO !J
      !do nothing for now but otherwise build sublstensor an Isend to other
   ENDIF
   lhs_current => lhs_current%next
ENDDO !I
IF(Bothpartioning)THEN
   !if the node have not yet recieved a message from the other nodes 
   !it will have to be a blocking operation so we have to wait for it 
   call lsmpi_blocking_recv_add_lstmemrealkbuf(jmat%lstmem_index,setting%node,setting%comm,MessageRecieved,setting%numnodes)
   call mem_dealloc(MessageRecieved)
ENDIF

call ls_free_lstensors(dmat_lhs,dmat_rhs,lhs_created,rhs_created)
if (doscreen)then
   IF(RHSpartioning)THEN
      !local gab created in addition to the full which is deallocated later
      Call ls_free_screeninglstensors(gabCS_rhs,gabCS_lhs,rhsCS_created,.FALSE.)
   ENDIF
   IF(LHSpartioning)THEN
      !local gab created in addition to the full which is deallocated later
      Call ls_free_screeninglstensors(gabCS_rhs,gabCS_lhs,.FALSE.,lhsCS_created)
   ENDIF
endif
call ls_free_tasks(tasks%lhs,lupri)
call ls_free_tasks(tasks%rhs,lupri)
DO iAO=1,4
  call freeMolecularOrbitalInfo(tasks%orbInfo(iAO))
ENDDO
Setting%sameMol = sameMolSave !HACK END
Setting%sameFrag = sameMolSave

call lsmpi_barrier(setting%comm)
CALL LS_GETTIM(t1,t2)

setting%output%resultTensor => jmat
IF(LHSpartioning)THEN
   !each node hold their own piece of the lstensor so we can now build the full fock matrix from this
   !for now we do nothing 
!   print*,'no reduction'
ELSE
   !a simple reduction sould give us the full jmat
   CALL lsmpi_lstensor_reduction(setting%output%resultTensor,infpar%master,setting%node,setting%comm)
ENDIF
CALL LS_GETTIM(t3,t4)
t(5) = t(5) + t3 - t1 !Node CPU time  += task CPU time
t(6) = t(6) + t4 - t2 !Node wall time += task wall time
CALL LS_GETTIM(t1,t2)
call lsmpi_time_harvest(setting%node,t(1),t(2),t(3),t(4),t(5),t(6),t(7),t(8),&
     & part(1),part(2),setting%numnodes,lupri,setting%scheme%intprint)

IF (setting%node.NE.infpar%master) THEN
   !put slave to sleep
   if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)

   IF(LHSpartioning)then
      IF(RHSpartioning.AND.rhs_created)THEN
         !Directly scalapack code because the slave is alive 
         !and normal mat_free wakes it up.      
         call PDM_MATRIXSYNC(DmatrixRHS)
         CALL FREE_IN_DARRAY(DmatrixRHS)
      ENDIF
      IF(lhs_created)THEN
         !Directly scalapack code because the slave is alive 
         !and normal mat_free wakes it up.          
         call PDM_MATRIXSYNC(DmatrixLHS)
         CALL FREE_IN_DARRAY(DmatrixLHS)
      ENDIF
      IF(PermuteResultTensor)THEN
         setting%output%resultTensor%PermuteResultTensor=.TRUE.
      ENDIF
      call retrieve_output_slave(lupri,setting)
   ENDIF
   call lstensor_free(setting%output%resultTensor)
   deallocate(setting%output%resultTensor)
   nullify(setting%output%resultTensor)
   IF(BOTHpartioning)THEN      
      call lstensor_free(jmat_other)
      deallocate(jmat_other)
      nullify(jmat_other)
   ENDIF
   IF(associated(setting%output%postprocess))THEN
      call mem_dealloc(setting%output%postprocess)
   ENDIF
   call typedef_free_setting(SETTING)
ELSE
   IF(RHSpartioning.AND.rhs_created)THEN
      IF(LHSpartioning)then
         call PDM_MATRIXSYNC(DmatrixRHS)
         !Directly scalapack code because the slave is alive 
         !and normal mat_free wakes it up.          
#ifdef VAR_SCALAPACK
         nsizeFULL = DmatrixRHS%nrow*DmatrixRHS%ncol
         nsizeLOCAL= DmatrixRHS%localnrow*DmatrixRHS%localncol
#endif
         call mem_deallocated_mem_type_matrix(nsizeLOCAL,nsizeFULL)
         CALL FREE_IN_DARRAY(DmatrixRHS)
      else
         call mat_free(DmatrixRHS)
      endif
   ENDIF
   IF(LHSpartioning.AND.lhs_created)THEN
      !Directly scalapack code because the slave is alive 
      !and normal mat_free wakes it up.          
      call PDM_MATRIXSYNC(DmatrixLHS)
#ifdef VAR_SCALAPACK
      nsizeFULL = DmatrixLHS%nrow*DmatrixLHS%ncol
      nsizeLOCAL= DmatrixLHS%localnrow*DmatrixLHS%localncol
#endif
      call mem_deallocated_mem_type_matrix(nsizeLOCAL,nsizeFULL)
      CALL FREE_IN_DARRAY(DmatrixLHS)
   ENDIF
   CALL LS_GETTIM(t3,t4)
   t(5) = t3 - t1 !Node CPU time  += task CPU time
   t(6) = t4 - t2 !Node wall time += task wall time
   CALL LS_GETTIM(t1,t2)
   IF(BOTHpartioning)THEN
      call lstensor_free(jmat_other)
      deallocate(jmat_other)
      nullify(jmat_other)
   ENDIF
   !Permutational Symmetry is done using matrix operations in retrieve_output
   IF(PermuteResultTensor)THEN
      IF(LHSpartioning)THEN
         setting%output%resultTensor%PermuteResultTensor=.TRUE.
         setting%output%postprocess = SymFromTriangularPostprocess
      ELSE
         call lstensor_full_symMat_from_triangularMat(setting%output%resultTensor)
      ENDIF
   ENDIF
   IF(LHSpartioning)THEN   
      setting%output%memdistResultTensor = .TRUE.
   ENDIF

   Call freeDaltonFragments(SETTING)
   setting%output%ndim = ndim_full
   
   if (doscreen) Call ls_free_screeninglstensors(gabCS_rhs_full,gabCS_lhs_full,rhsCS_created,lhsCS_created)
   
   CALL LS_GETTIM(t3,t4)
   t(5) = t3 - t1 !Node CPU time  += task CPU time
   t(6) = t4 - t2 !Node wall time += task wall time
ENDIF
#else
  call lsquit('ls_jengine_memdist but no MPI',-1) 
#endif
END SUBROUTINE ls_jengine_memdist

!> \brief retrive a full block from a type matrix
!> \author S. Reine
!> \date 2010
!> \param Dmat the matrix type
!> \param Dblock the output block 
!> \param ndmat the number of matrices
!> \param n1 the number of functions in the small block dim 1
!> \param n2 the number of functions in the small block dim 2
!> \param s1 the start index of the small block dim 1
!> \param s2 the start index of the small block dim 2
!> \param permute if the permutational symmetry is used 
SUBROUTINE ls_mat_retrive_block(Dmat,Dblock,ndmat,n1,n2,s1,s2,permute)
implicit none
Integer       :: ndmat,n1,n2,s1,s2
type(matrixp) :: Dmat(ndmat)
Real(realk)   :: Dblock(n1,n2,ndmat)
Logical       :: permute
!
Integer :: idmat
DO idmat = 1, ndmat
  call mat_retrieve_block(Dmat(idmat)%p,Dblock(:,:,iDmat),n1,n2,s1,s2)
  IF (permute) THEN
    call ls_add_sym_block(Dmat(idmat)%p,Dblock(:,:,iDmat),n1,n2,s1,s2)
  ENDIF
ENDDO
END SUBROUTINE ls_mat_retrive_block

!> \brief add symmetry block from full block 
!> \author S. Reine
!> \date 2010
!> \param Dmat the matrix type
!> \param Dblock the output block 
!> \param n1 the number of functions in the small block dim 1
!> \param n2 the number of functions in the small block dim 2
!> \param s1 the start index of the small block dim 1
!> \param s2 the start index of the small block dim 2
SUBROUTINE ls_add_sym_block(Dmat,Dblock,n1,n2,s1,s2)
implicit none
TYPE(Matrix) :: Dmat
Integer      :: n1,n2,s1,s2
Real(realk)  :: Dblock(n1,n2)
!
Real(realk),pointer :: Dtrans(:,:)
Integer                 :: i1,i2

call mem_alloc(Dtrans,n2,n1)
call mat_retrieve_block(Dmat,Dtrans,n2,n1,s2,s1)
DO i2=1,n2
  DO i1=1,n1
     Dblock(i1,i2) = Dblock(i1,i2) + Dtrans(i2,i1)
  ENDDO
ENDDO
call mem_dealloc(Dtrans)
END SUBROUTINE ls_add_sym_block

!> \brief calculate the coulomb matrix using the jengine method
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_jengine1(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)    :: INT_INPUT
TYPE(AOITEM),target    :: AObuild(4)
Integer                :: nAObuilds,idmat,nAtoms,iMol,nmat
integer                :: nAtoms1,nAtoms2,nAtoms3,nAtoms4,ndmat
Logical                :: SameAllFrag,dograd
type(lstensor),pointer :: result_tensor_save
real(realk)            :: ts,te

CALL init_integral_input(INT_INPUT,SETTING)

Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
! OBS: For the calculation of screeningmatrices, the initIntegralOutputDims is called and thus
!      allocates a new resultTensor. Saved here. 
! ToDo: Find a better solution!!
result_tensor_save => setting%output%resultTensor
NULLIFY(setting%output%resultTensor)
CALL AttachScreeningMatricesToInput(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
IF (associated(setting%output%resultTensor)) DEALLOCATE(setting%output%resultTensor)
setting%output%resultTensor => result_tensor_save
Int_input%DO_FMM = SETTING%SCHEME%FMM .AND. (Oper.EQ.CoulombOperator)
Int_input%NonClassical_SCREEN = Int_input%DO_FMM
Int_input%OE_SCREEN = SETTING%SCHEME%OE_SCREEN .AND. (Oper.EQ.OverlapOperator)

INT_INPUT%operator = Oper

CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

INT_INPUT%sameODs      = .FALSE.
INT_INPUT%DO_FOCK      = .TRUE.
INT_INPUT%DO_COULOMB   = .TRUE.
INT_INPUT%DO_JENGINE   = .TRUE.
INT_INPUT%DO_DAJENGINE = SETTING%SCHEME%DAJENGINE

call set_input_from_spec(INT_INPUT,SPEC,AO1,AO2,AO3,AO4,Oper,lupri,dograd,.FALSE.)

IF (Spec.EQ.MagDerivSpec) THEN
!ToDo: Fix this - will not work
  dograd = .false.
ENDIF

setting%output%dograd = dograd

CALL ls_setDensityDimensions(int_input,setting,lupri)
call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
CALL free_screening_matrices(INT_INPUT,SETTING,LUPRI,LUERR)

END SUBROUTINE ls_jengine1

!> \brief calculate the classical part of the coulomb matrix using the jengine method (in a FMM sense)
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular ' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_jengineClassical(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR,ndmatRHS
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,idmat,nAtoms,iMol

CALL init_integral_input(INT_INPUT,SETTING)


Int_input%DO_FMM = SETTING%SCHEME%FMM .AND. (Oper.EQ.CoulombOperator)


IF(Int_input%DO_FMM)THEN
   Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
   INT_INPUT%DO_GRADIENT = Spec.EQ.GradientSpec
   Int_input%NonClassical_SCREEN = Int_input%DO_FMM
   Int_input%OE_SCREEN = SETTING%SCHEME%OE_SCREEN .AND. (Oper.EQ.OverlapOperator)
   INT_INPUT%operator = Oper
   CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
   INT_INPUT%sameODs      = .FALSE.
   INT_INPUT%DO_FOCK      = .TRUE.
   INT_INPUT%DO_COULOMB   = .TRUE.
   INT_INPUT%DO_JENGINE   = .TRUE.
   
   IF (INT_INPUT%DO_GRADIENT) THEN
     IF(INT_INPUT%DO_MMGRD) THEN
        setting%Output%dograd = .true.
        nAtoms =  setting%molecule(1)%p%nAtoms
        DO iMol=2,4
         IF(setting%molecule(iMol)%p%nAtoms.NE.nAtoms) THEN 
          CALL LSQUIT('Error in ls_jengine1. nAtoms inconsistency!',lupri)
         ENDIF
        ENDDO
        !Molecular gradient
        CALL ls_setDensityDimensions(INT_INPUT,SETTING,lupri)
        CALL GradClassical(setting%Output,nAtoms,INT_INPUT)
     ENDIF
   ELSE
     CALL ls_setDensityDimensions(INT_INPUT,SETTING,lupri)
     CALL JmatClassical(setting%Output,setting%Output%ndim(1),setting%Output%ndim(2),INT_INPUT)
   ENDIF
   CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
ENDIF

END SUBROUTINE ls_jengineClassical

!> \brief calculate the classical part of the coulomb matrix using the jengine method (in a FMM sense)
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular ' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_jengineClassicalGRAD(gradient,AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR,natoms)
implicit none
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec,natoms
Real(realk)          :: gradient(3*natoms)
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR,ndmatRHS
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,idmat,iMol
CALL init_integral_input(INT_INPUT,SETTING)
Int_input%DO_FMM = SETTING%SCHEME%FMM .AND. (Oper.EQ.CoulombOperator)
IF(Int_input%DO_FMM)THEN
   Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
   INT_INPUT%DO_GRADIENT = Spec.EQ.GradientSpec
   Int_input%NonClassical_SCREEN = Int_input%DO_FMM
   Int_input%OE_SCREEN = SETTING%SCHEME%OE_SCREEN.AND.(Oper.EQ.OverlapOperator)
   INT_INPUT%operator = Oper
   CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,&
        & SETTING,LUPRI,LUERR,.TRUE.)
   INT_INPUT%sameODs      = .FALSE.
   INT_INPUT%DO_FOCK      = .TRUE.
   INT_INPUT%DO_COULOMB   = .TRUE.
   INT_INPUT%DO_JENGINE   = .TRUE.   
   IF(INT_INPUT%DO_MMGRD) THEN
      setting%Output%dograd = .true.
      nAtoms =  setting%molecule(1)%p%nAtoms
      DO iMol=2,4
         IF(setting%molecule(iMol)%p%nAtoms.NE.nAtoms) THEN 
            CALL LSQUIT('Error in ls_jengine1. nAtoms inconsistency!',lupri)
         ENDIF
      ENDDO
      !Molecular gradient
      CALL ls_setDensityDimensions(INT_INPUT,SETTING,lupri)
      CALL GradClassicalGRAD(gradient,nAtoms,INT_INPUT)
   ENDIF
   CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
ENDIF
END SUBROUTINE ls_jengineClassicalGRAD

!> \brief calculate the classical part of the coulomb matrix using the jengine method (in a FMM sense)
!> \author S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param Spec the label for the type of calc 'Regular ' or 'Gradient'
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_jengineClassicalMat(MAT,AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
implicit none
type(matrix)         :: MAT
integer              :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR,ndmatRHS
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds,idmat,nAtoms,iMol
CALL init_integral_input(INT_INPUT,SETTING)
Int_input%DO_FMM = SETTING%SCHEME%FMM .AND. (Oper.EQ.CoulombOperator)
IF(Int_input%DO_FMM)THEN
   Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
   INT_INPUT%DO_GRADIENT = Spec.EQ.GradientSpec
   Int_input%NonClassical_SCREEN = Int_input%DO_FMM
   Int_input%OE_SCREEN = SETTING%SCHEME%OE_SCREEN.AND.(Oper.EQ.OverlapOperator)
   INT_INPUT%operator = Oper
   CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,&
        & SETTING,LUPRI,LUERR,.TRUE.)
   INT_INPUT%sameODs      = .FALSE.
   INT_INPUT%DO_FOCK      = .TRUE.
   INT_INPUT%DO_COULOMB   = .TRUE.
   INT_INPUT%DO_JENGINE   = .TRUE.   
   CALL ls_setDensityDimensions(INT_INPUT,SETTING,lupri)
   CALL JmatClassicalMAT(MAT,MAT%nrow,MAT%ncol,INT_INPUT)
   CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
ENDIF
END SUBROUTINE ls_jengineClassicalMat

!> \brief wrapper to calculate the screening matrices required for quadratic or better scaling
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param Input the integral input specification
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE AttachScreeningMatricesToInput(Input,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
implicit none
TYPE(INTEGRALINPUT)  :: Input
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!logical              :: family
integer              :: AO1,AO2,AO3,AO4,Oper

INPUT%sameODs    = (AO1.EQ.AO3) .AND. (AO2.EQ.AO4).AND.setting%samefrag(1,3).AND.setting%samefrag(2,4)
!clean 
nullify(input%LST_GAB_RHS)
nullify(input%LST_GAB_LHS)
! Contracted screening matrix/matrices
IF(Input%CS_SCREEN.OR.Input%PS_SCREEN.OR.Input%MBIE_SCREEN) THEN
  IF(SETTING%SCHEME%INTPRINT .GT. 10) WRITE(LUPRI,'(2X,A)')&
          &'call GET_GAB_MATRIX - the  screening matrix (THE BOOK 9.12.26)'
  IF (Oper.EQ.NucpotOperator.AND.AO1.EQ.AORdefault) THEN
!???????? SIMEN ???????? Alternatives to Coulomb here?
    CALL AttachScreenMatrixToInput('LHS',Input,AO1,AO2,CoulombOperator,SETTING,LUPRI,LUERR)
    CALL AttachScreenMatrixToInput('RHS',Input,AO3,AO4,NucleiOperator,SETTING,LUPRI,LUERR)
  ELSEIF (Oper.EQ.NucpotOperator.AND.AO1.EQ.AOdfCABS) THEN
    CALL AttachScreenMatrixToInput('LHS',Input,AO1,AO2,CoulombOperator,SETTING,LUPRI,LUERR)
    CALL AttachScreenMatrixToInput('RHS',Input,AO3,AO4,NucleiOperator,SETTING,LUPRI,LUERR)
  ELSEIF (Oper.EQ.NucpotOperator.AND.AO3.EQ.AORdefault) THEN
    CALL AttachScreenMatrixToInput('LHS',Input,AO1,AO2,NucleiOperator,SETTING,LUPRI,LUERR)
    CALL AttachScreenMatrixToInput('RHS',Input,AO3,AO4,CoulombOperator,SETTING,LUPRI,LUERR)
  ELSE
    CALL AttachScreenMatrixToInput('LHS',Input,AO1,AO2,Oper,SETTING,LUPRI,LUERR)
    CALL AttachScreenMatrixToInput('RHS',Input,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
  ENDIF
ENDIF

END SUBROUTINE AttachScreeningMatricesToInput

!> \brief free the calculated screening matrices
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param Input the integral input specification
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE free_screening_matrices(Input,SETTING,LUPRI,LUERR)
implicit none
TYPE(INTEGRALINPUT)  :: Input
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
IF(Input%CS_SCREEN.OR.Input%PS_SCREEN.OR.Input%MBIE_SCREEN) THEN
   IF(.NOT.Input%GAB_LHSusePointer)THEN
      call lstensor_free(input%LST_GAB_LHS)      
      deallocate(input%LST_GAB_LHS)
   ENDIF
   IF(.NOT.Input%GAB_RHSusePointer)THEN
      call lstensor_free(input%LST_GAB_RHS)
      deallocate(input%LST_GAB_RHS)
   ENDIF
   nullify(input%LST_GAB_RHS)
   nullify(input%LST_GAB_LHS)
ENDIF
END SUBROUTINE free_screening_matrices

!> \brief calculate the screening matrix
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE ls_getScreenIntegrals(AO1,AO2,Oper,CS_INT,PS_INT,MBIE_INT,SETTING,LUPRI,LUERR)
implicit none
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
integer              :: AO1,AO2,Oper!,intType
Logical              :: CS_INT,PS_INT,MBIE_INT
!
Integer             :: n1,n2,iAO,THR
Type(BASIS_PT)      :: BasisBackup(4)
LOGICAL             :: LHS,MBIE_SCREENINT
!
LHS=.TRUE.
SETTING%FRAGMENT(3)%p => SETTING%FRAGMENT(1)%p
SETTING%FRAGMENT(4)%p => SETTING%FRAGMENT(2)%p
DO iAO=1,4
  BasisBackup(iAO)%p => Setting%basis(iAO)%p
ENDDO
Setting%basis(3)%p => Setting%basis(1)%p
Setting%basis(4)%p => Setting%basis(2)%p
MBIE_SCREENINT = MBIE_INT.AND.(Oper.EQ.CoulombOperator.OR.Oper.EQ.ErfcOperator.OR.Oper.EQ.CAMOperator)
CALL ls_getScreenIntegrals1(AO1,AO2,Oper,CS_INT,PS_INT,MBIE_SCREENINT,SETTING,LUPRI,LUERR,LHS)
DO iAO=1,4
  Setting%basis(iAO)%p => BasisBackup(iAO)%p
ENDDO

END SUBROUTINE ls_getScreenIntegrals

!> \brief calculate the screening matrix
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param IntegralSide The side of the screening matrix (LHS or RHS) 
!> \param Input the integral input specification
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1/3
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2/4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
SUBROUTINE AttachScreenMatrixToInput(integralSide,Input,&
     & AO1,AO2,Oper,SETTING,LUPRI,LUERR)
implicit none
TYPE(INTEGRALINPUT),target  :: Input
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
integer        :: AO1,AO2,Oper
Character*(*) :: integralSide
!
integer             :: intType
type(lstensor),pointer :: GAB
Integer             :: iprint,n1,n2,i1,i2,s1,s2
Logical             :: nuclei,LHS,RHS
integer(kind=short) :: maxGabElm,maxprimGabElm,maxgabelm2
Character(80)       :: Filename
Character(53)       :: identifier
Character(22)       :: label1,label2
Real(realk),pointer :: screenMat(:,:)
integer             :: f1,f2,iAO,AOA,AOB,THR,ilst,I,J,THR2,ndim2(5)
Type(BASIS_PT)      :: BasisBackup(4)
logical             :: IntegralTransformGC,PrimScreenIntType,RHSGABusePointer
logical             :: attach_to_input,saveGABtoMem
logical             :: FoundInMem,recalcGab,MBIE_SCREENINT
real(realk)         :: tstart,tend
integer,pointer     :: postprocess(:)
integer(kind=short),pointer :: MAT2(:,:)
ndim2 = setting%output%ndim
! Multipole Based Integral Estimate Screening   
MBIE_SCREENINT = setting%scheme%MBIE_SCREEN.AND.(Oper.EQ.CoulombOperator.OR.Oper.EQ.ErfcOperator.OR.Oper.EQ.CAMOperator)
LHS = integralSide.EQ.'LHS'
!PrimScreenIntType=intType.EQ.'Primitive'
intType=Contractedinttype
attach_to_input=.true.
! see if screening matrix have been attached to setting
! if that is the case, we point the inputGab to the settingGab
IF(LHS)THEN
   IF(associated(setting%LST_GAB_LHS))THEN
      input%LST_GAB_LHS => setting%LST_GAB_LHS
      IF(Input%PS_SCREEN)Input%PS_MAXELM_LHS = setting%PS_MAXELM_LHS 
      IF(Input%CS_SCREEN)Input%CS_MAXELM_LHS = setting%CS_MAXELM_LHS 
      Input%GAB_LHSusePointer = .TRUE. 
      attach_to_input=.false.
   ENDIF
ELSE
   IF(associated(setting%LST_GAB_RHS))THEN
      input%LST_GAB_RHS => setting%LST_GAB_RHS
      IF(Input%PS_SCREEN)Input%PS_MAXELM_RHS = setting%PS_MAXELM_RHS 
      IF(Input%CS_SCREEN)Input%CS_MAXELM_RHS = setting%CS_MAXELM_RHS 
      Input%GAB_RHSusePointer = .TRUE. 
      attach_to_input=.false.
   ENDIF
ENDIF

IF(attach_to_input)THEN
   ! We attach screening directly to input
   saveGABtoMem = SETTING%SCHEME%saveGABtoMem
   recalcGab = SETTING%SCHEME%recalcGab
   IF(.NOT.saveGABtoMem) ReCalcGab=.TRUE. 
   IF(AO1.EQ.AOpCharge) ReCalcGab = .TRUE.
   IntegralTransformGC = .FALSE.
   iprint = SETTING%SCHEME%intPrint

   IF(LHS)THEN
      AOA=1
      AOB=2
   ELSE
      AOA=3
      AOB=4
   ENDIF
   IF((SETTING%batchindex(AOA).NE. 0).AND.SETTING%molID(AOA).NE. 0)THEN
      WRITE(LUPRI,'(1X,A)') 'Input error in AttachScreenMatrixToInput: only made for either batchindex use og molID use '
      CALL LSQUIT('Input error in AttachScreenMatrixToInput',lupri)      
   ENDIF

   s1 = SETTING%molID(AOA) !not 0, when the molecule pointers are different molecules
   s2 = SETTING%molID(AOB)
   IF(SETTING%batchindex(AOA).NE. 0)THEN
      s1=SETTING%batchindex(AOA) !not 0, only when a single AObatch is computed
   ENDIF
   IF(SETTING%batchindex(AOB).NE. 0)THEN
      s2=SETTING%batchindex(AOB)
   ENDIF
   IF(ASSOCIATED(Setting%FRAGMENT(AOA)%p)) THEN
      n1 = getNbasis(AO1,intType,SETTING%FRAGMENT(AOA)%p,LUPRI)
   ELSE
      n1 = 1
   ENDIF
   IF(ASSOCIATED(Setting%FRAGMENT(AOB)%p)) THEN
      n2 = getNbasis(AO2,intType,SETTING%FRAGMENT(AOB)%p,LUPRI)
   ELSE
      n2=1
   ENDIF
   IF(SETTING%batchindex(AOA).NE. 0)THEN
      n1=SETTING%batchdim(AOA) 
      IF (intType.EQ.Primitiveinttype)CALL LSQUIT('PrimScreening and batchindex not implemented',lupri)
   ENDIF
   IF(SETTING%batchindex(AOB).NE. 0)THEN
      n2=SETTING%batchdim(AOB) 
      IF (intType.EQ.Primitiveinttype)CALL LSQUIT('PrimScreening and batchindex not implemented',lupri)
   ENDIF

   IF (iprint.gt. 5) THEN
      WRITE(LUPRI,'(1X,A)') 'Output from AttachScreenMatrixToInput'
      WRITE(LUPRI,'(3X,A,I3)')   'AO1:          ',AO1
      WRITE(LUPRI,'(3X,A,I3)')   'AO2:          ',AO2
      WRITE(LUPRI,'(3X,A,I3)')   'intType:      ',intType
      WRITE(LUPRI,'(3X,2A)')   'integralSide: ',integralSide
      WRITE(LUPRI,'(3X,A,I3)') 'n1:           ',n1
      WRITE(LUPRI,'(3X,A,I3)') 'n2:           ',n2
      WRITE(LUPRI,'(3X,A,I3)') 's1:           ',s1
      WRITE(LUPRI,'(3X,A,I3)') 's2:           ',s2
      WRITE(LUPRI,'(3X,A,I3)') 'Oper          ',Oper
   ENDIF

   !determine if we can use a pointer to the LHS GAB for the RHS GAB 
   !we only do this for the right hand side, and if sameOD and if LHS is associated 
   call determine_RHSGABusePointer(RHSGABusePointer,LHS,Input,intType,Oper)

   IF(RHSGABusePointer)THEN
      input%LST_GAB_RHS => input%LST_GAB_LHS
      Input%GAB_RHSusePointer = .TRUE. 
      IF(Input%PS_SCREEN) Input%PS_MAXELM_RHS = Input%PS_MAXELM_LHS
      IF(Input%CS_SCREEN) Input%CS_MAXELM_RHS = Input%CS_MAXELM_LHS
   ELSE
      Input%GAB_RHSusePointer = .FALSE. 
      THR = ABS(NINT(LOG10(SETTING%SCHEME%intTHRESHOLD)))
      !Find Filename
      CALL io_get_CSidentifier(identifier,THR,Setting%FRAGMENT(AOA)%p,Setting%FRAGMENT(AOB)%p,&
     &                         Input%CS_SCREEN,Input%PS_SCREEN)
      CALL io_get_filename(Filename,identifier,AO1,AO2,AO1,AO2,s1,s2,Oper,intType,SETTING%SCHEME%FRAGMENT,LUPRI,LUERR)
      call determine_lst_in_screenlist(Filename,FoundInMem,SETTING%IO)

      IF(.NOT.FoundInMem)THEN
         DO I=1,10
            THR2=THR+I
            CALL io_get_CSidentifier(identifier,THR2,Setting%FRAGMENT(AOA)%p,Setting%FRAGMENT(AOB)%p,&
     &                               Input%CS_SCREEN,Input%PS_SCREEN)
            CALL io_get_filename(Filename,identifier,AO1,AO2,AO1,AO2,s1,s2,Oper,intType,SETTING%SCHEME%FRAGMENT,LUPRI,LUERR)
            call determine_lst_in_screenlist(Filename,FoundInMem,SETTING%IO)
            IF(FoundInMem)EXIT
         ENDDO

         IF(.NOT.FoundInMem)THEN
            CALL io_get_CSidentifier(identifier,THR,Setting%FRAGMENT(AOA)%p,Setting%FRAGMENT(AOB)%p,&
     &                               Input%CS_SCREEN,Input%PS_SCREEN)
            CALL io_get_filename(Filename,identifier,AO1,AO2,AO1,AO2,s1,s2,Oper,intType,SETTING%SCHEME%FRAGMENT,LUPRI,LUERR)
            call determine_lst_in_screenlist(Filename,FoundInMem,SETTING%IO)
         ENDIF
      ENDIF
      if(FoundInMem.AND..NOT.ReCalcGab)then 
         !IF FoundInMem We have already calculated this so set pointers
         IF(LHS)THEN
            call screen_associate(input%LST_GAB_LHS,Filename,FoundInMem)
            iF(Input%PS_SCREEN)Input%PS_MAXELM_LHS = input%LST_GAB_LHS%maxprimGabElm
            iF(Input%CS_SCREEN)Input%CS_MAXELM_LHS = input%LST_GAB_LHS%maxGabElm
            Input%GAB_LHSusePointer = .TRUE.
         ELSE
            call screen_associate(input%LST_GAB_RHS,Filename,FoundInMem)
            iF(Input%PS_SCREEN)Input%PS_MAXELM_RHS = input%LST_GAB_RHS%maxprimGabElm
            iF(Input%CS_SCREEN)Input%CS_MAXELM_RHS = input%LST_GAB_RHS%maxGabElm
            Input%GAB_RHSusePointer = .TRUE.
         ENDIF
      else 
         !Calc Screening lstensor
         IF(saveGABtoMem.AND..NOT.ReCalcGab)then
            call screen_add_associate_item(GAB,filename)
         ENDIF
         IF(ReCalcGab)then
            IF(LHS)THEN
               IF(.NOT.associated(input%LST_GAB_LHS))THEN
                  nullify(input%LST_GAB_LHS)
                  allocate(input%LST_GAB_LHS)
                  call lstensor_nullify(input%LST_GAB_LHS)
               ELSE
                  call lsquit('A input%LST_GAB_LHS ass',-1)
               ENDIF
               Input%GAB_LHSusePointer = .FALSE.
               GAB => input%LST_GAB_LHS
            ELSE      
               IF(.NOT.associated(input%LST_GAB_RHS))THEN
                  nullify(input%LST_GAB_RHS)
                  allocate(input%LST_GAB_RHS)
                  call lstensor_nullify(input%LST_GAB_RHS)
               ELSE
                  call lsquit('A input%LST_GAB_RHS ass',-1)
               ENDIF
               Input%GAB_RHSusePointer = .FALSE.
               GAB => input%LST_GAB_RHS
            ENDIF
         ENDIF
         DO iAO=1,4
            BasisBackup(iAO)%p => Setting%basis(iAO)%p
         ENDDO
         IF (LHS) THEN
            Setting%basis(3)%p => Setting%basis(1)%p
            Setting%basis(4)%p => Setting%basis(2)%p
         ELSE
            Setting%basis(1)%p => Setting%basis(3)%p
            Setting%basis(2)%p => Setting%basis(4)%p
         ENDIF
         IF ((AO1.EQ.AOEmpty).AND.(AO2.EQ.AOEmpty)) THEN
            call mem_alloc(screenMat,n1,n2)
            screenMat(1,1) = 1.0E0_realk
            call build_Nuclearlstensor(ScreenMat,GAB,1)
            call mem_dealloc(screenMat)
         ELSEIF (Oper.EQ.NucleiOperator) THEN
            call mem_alloc(screenMat,n1,n2)
            CALL ls_getNucScreenIntegrals(GAB,ScreenMat,Filename,AO1,AO2,Oper,&
                 & SETTING,LUPRI,LUERR,LHS)
            call mem_dealloc(screenMat)
         ELSE
            call mem_alloc(postprocess,size(setting%Output%postprocess))
            postprocess = setting%Output%postprocess
            setting%Output%postprocess = 0
            call initIntegralOutputDims1(setting%Output,n1,n2,1,1,1)
            IF(.NOT.ReCalcGab)CALL LSTIMER('START ',tstart,tend,lupri)
            CALL ls_getScreenIntegrals1(AO1,AO2,Oper,Input%CS_SCREEN,&
                 & Input%PS_SCREEN,MBIE_SCREENINT,SETTING,LUPRI,LUERR,LHS)
            IF(.NOT.ReCalcGab)CALL LSTIMER(Filename(1:15),tstart,tend,lupri)
            CALL retrieve_screen_output(lupri,setting,GAB,IntegralTransformGC)
            call mem_alloc(setting%Output%postprocess,size(postprocess))
            setting%Output%postprocess = postprocess
            call mem_dealloc(postprocess)
         ENDIF
         IF(Input%CS_SCREEN) maxGABElm = GAB%maxgabelm
         IF(Input%PS_SCREEN) maxprimGABElm = GAB%maxprimgabelm
!         call determine_maxelm(maxGABElm,GAB) needs to be done in 
!         build_Nuclearlstensor,ls_getNucScreenIntegrals,ls_getScreenIntegrals1
         DO iAO=1,4
            Setting%basis(iAO)%p => BasisBackup(iAO)%p
         ENDDO
         !Place it somewhere
         IF(ReCalcGab)THEN
            IF(LHS)THEN
               IF(Input%PS_SCREEN)Input%PS_MAXELM_LHS = maxprimGABElm 
               IF(Input%CS_SCREEN)Input%CS_MAXELM_LHS = maxGABElm
            ELSE      
               IF(Input%PS_SCREEN)Input%PS_MAXELM_RHS = maxprimGABElm
               IF(Input%CS_SCREEN)Input%CS_MAXELM_RHS = maxGABElm
            ENDIF
         ENDIF
         IF(saveGABtoMem.AND..NOT.ReCalcGab)then
            IF(LHS)THEN
               call screen_associate(input%LST_GAB_LHS,Filename,FoundInMem)
               IF(Input%PS_SCREEN)Input%PS_MAXELM_LHS = input%LST_GAB_LHS%maxprimGabElm
               IF(Input%CS_SCREEN)Input%CS_MAXELM_LHS = input%LST_GAB_LHS%maxGabElm
               Input%GAB_LHSusePointer = .TRUE.
            ELSE      
               call screen_associate(input%LST_GAB_RHS,Filename,FoundInMem)
               IF(Input%PS_SCREEN)Input%PS_MAXELM_RHS = input%LST_GAB_RHS%maxprimGabElm
               IF(Input%CS_SCREEN)Input%CS_MAXELM_RHS = input%LST_GAB_RHS%maxGabElm
               Input%GAB_RHSusePointer = .TRUE.
            ENDIF
         endif
      ENDIF
   ENDIF
ENDIF
setting%output%ndim = ndim2

END SUBROUTINE AttachScreenMatrixToInput

subroutine determine_RHSGABusePointer(RHSGABusePointer,LHS,input,intType,Oper)
implicit none
LOGICAL,intent(out) :: RHSGABusePointer
LOGICAL,intent(in)  :: LHS
TYPE(INTEGRALINPUT),intent(in)  :: Input
integer        :: Oper,intType

RHSGABusePointer = .TRUE.
IF(LHS)RHSGABusePointer = .FALSE.
IF(.NOT.Input%sameODs)RHSGABusePointer = .FALSE.
IF (intType.EQ.Primitiveinttype) THEN
   IF(.NOT.ASSOCIATED(input%LST_GAB_LHS))RHSGABusePointer = .FALSE.
ELSE
   IF(.NOT.ASSOCIATED(input%LST_GAB_LHS))RHSGABusePointer = .FALSE.
ENDIF
IF(Oper.EQ.NucleiOperator)RHSGABusePointer = .FALSE.
END subroutine DETERMINE_RHSGABUSEPOINTER

subroutine Determine_GAB_filename(SETTING,AO1,AO2,s1,s2,Oper,intType,LUPRI,LUERR,LHSGAB,CS_SCREEN,PS_SCREEN,Filename)
  implicit none
  Type(LSSETTING)     :: SETTING
  Integer             :: LUPRI,LUERR,AO1,AO2,Oper,intType,s1,s2
  LOGICAL             :: LHSGAB,CS_SCREEN,PS_SCREEN
  Character(80)       :: Filename
  !
  Character(53)       :: identifier
  Character(22)       :: label1,label2
  integer :: THR,AOA,AOB
  IF(LHSGAB)THEN
     AOA=1
     AOB=2
  ELSE
     AOA=3
     AOB=4
  ENDIF
  !Find Filename
  THR = ABS(NINT(LOG10(SETTING%SCHEME%intTHRESHOLD)))
  IF (ASSOCIATED(Setting%FRAGMENT(AOA)%p)) THEN
     label1 = Setting%FRAGMENT(AOA)%p%label
  ELSE
     label1 = 'Empty_________________'
  ENDIF
  IF (ASSOCIATED(Setting%FRAGMENT(AOB)%p)) THEN
     label2 = Setting%FRAGMENT(AOB)%p%label
  ELSE
     label2 = 'Empty_________________'
  ENDIF
  !$OMP CRITICAL (ifortwrite)
  IF(THR.GT.9)THEN 
     write(identifier,'(A2,I2,A1,A22,A1,A22,A1,L1,L1)')'CS',THR,'_',label1,'_',label2,'_',CS_SCREEN,PS_SCREEN
  ELSE
     write(identifier,'(A2,I1,A1,A22,A1,A22,A1,L1,L1)')'CS0',THR,'_',label1,'_',label2,'_',CS_SCREEN,PS_SCREEN
  ENDIF
  !$OMP END CRITICAL (ifortwrite)
  CALL io_get_filename(Filename,identifier,AO1,AO2,AO1,AO2,s1,s2,Oper,intType,SETTING%SCHEME%FRAGMENT,LUPRI,LUERR)
END subroutine DETERMINE_GAB_FILENAME

!> \brief calculate the screening integrals
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param Filename the filename used to store the screening matrix
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1/3
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2/4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param LHSGAB if this is a LHS screening matrix
SUBROUTINE ls_getScreenIntegrals1(AO1,AO2,Oper,CS_SCREEN,PS_SCREEN,MBIE_SCREEN,&
     & SETTING,LUPRI,LUERR,LHSGAB)
implicit none
integer              :: AO1,AO2,Oper!,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
LOGICAL              :: LHSGAB,CS_SCREEN,PS_SCREEN,MBIE_SCREEN
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
!TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
Real(realk)          :: tstart,tend
Integer              :: iAtom,iFragment,iAO,iprint
integer :: ndim2(5)
integer(kind=short) :: CS_THRLOG
logical :: MBIE_INT
ndim2 = setting%output%ndim 
!CALL ls_setDefaultFragments(setting)
CALL init_integral_input(INT_INPUT,SETTING)

INT_INPUT%operator = Oper
INT_INPUT%CS_int=.TRUE.
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
INT_INPUT%CS_SCREEN = .FALSE.
INT_INPUT%PS_SCREEN = .FALSE.
IPRINT = SETTING%SCHEME%INTPRINT
MBIE_INT = MBIE_SCREEN
IF (MBIE_INT.AND.((AO1.EQ.AOEmpty).OR.(AO2.EQ.AOEmpty))) THEN
   WRITE(lupri,*)'WARNING: MBIE not working for densfit basis (Aux Empty|Aux Empty) but maybe it should'
   MBIE_INT = .FALSE.   
ENDIF
IF (Oper.EQ.NucleiOperator) THEN
   CALL LSQUIT('error in getScreenIntegrals',lupri)
ELSEIF ((AO1.EQ.AOEmpty).AND.(AO2.EQ.AOEmpty)) THEN
   CALL LSQUIT('result is screenMat(1,1) = 1.0E0_realk but this should have been done a different place',lupri)
ELSE
 IF(setting%Output%RealGabMatrix)Then
    IF(CS_SCREEN)THEN   
       CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,Contractedinttype,AObuild,nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
       CALL initIntegralOutputDims1(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,1)
       Int_Input%CS_int = .TRUE.
       Int_Input%PS_int = .FALSE.
       nullify(setting%output%resultTensor)
       allocate(setting%output%resultTensor)
       call init_lstensor_5dim(setting%output%resultTensor,Int_Input%AO(1)%p,Int_Input%AO(2)%p,&
            & Int_Input%AO(3)%p,Int_Input%AO(4)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,&
            & 1,1,.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,lupri)
!       call init_cs_lstensor(setting%output%ScreenTensor,Int_Input%AO(1)%p,&
!           &Int_Input%AO(2)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),INT_INPUT%OD_SCREEN,lupri)
      call MAIN_INTEGRAL_DRIVER(LUPRI,IPRINT,INT_INPUT,setting%OUTPUT)
!      call set_lst_maxgabelms(setting%output%ScreenTensor)
      CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
   ENDIF
 ELSE
   nullify(setting%output%ScreenTensor)
   allocate(setting%output%ScreenTensor)
   call lstensor_nullify(setting%output%ScreenTensor)
   IF(CS_SCREEN)THEN   
      CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,Contractedinttype,AObuild,nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
      CALL initIntegralOutputDims1(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,1)
      Int_Input%CS_int = .TRUE.
      Int_Input%PS_int = .FALSE.
      call init_cs_lstensor(setting%output%ScreenTensor,Int_Input%AO(1)%p,&
           &Int_Input%AO(2)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),INT_INPUT%OD_SCREEN,lupri)
      call MAIN_INTEGRAL_DRIVER(LUPRI,IPRINT,INT_INPUT,setting%OUTPUT)
      call set_lst_maxgabelms(setting%output%ScreenTensor)
      CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
   ENDIF
   IF(PS_SCREEN)THEN
      CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,Primitiveinttype,AObuild,&
           & nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
      CALL initIntegralOutputDims1(setting%Output,INT_INPUT%AOdim(1),&
           & INT_INPUT%AOdim(2),1,1,1)
      Int_Input%CS_int = .FALSE.
      Int_Input%PS_int = .TRUE.
      call init_ps_lstensor(setting%output%ScreenTensor,Int_Input%AO(1)%p,&
           &Int_Input%AO(2)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),INT_INPUT%OD_SCREEN,lupri)
      call MAIN_INTEGRAL_DRIVER(LUPRI,IPRINT,INT_INPUT,setting%OUTPUT)
      call set_lst_maxprimgab(setting%output%ScreenTensor)
      call set_lst_maxprimgabelms(setting%output%ScreenTensor)
      CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
   ENDIF
   IF(MBIE_INT)THEN
      CALL SetInputAO(INT_INPUT,AO1,AO2,AOEmpty,AOEmpty,Contractedinttype,AObuild,&
           & nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
      CALL initIntegralOutputDims1(setting%Output,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),1,1,1)
      call init_MBIE_lstensor_5dim(setting%output%ScreenTensor,Int_Input%AO(1)%p,&
           &Int_Input%AO(2)%p,.FALSE.,lupri)
      Int_Input%CS_int = .FALSE.
      Int_Input%PS_int = .FALSE.
      INT_INPUT%MBIE_INT=.TRUE.
      call MBIE_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,setting%OUTPUT)
      INT_INPUT%MBIE_INT=.FALSE.
      CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
   ENDIF
 ENDIF
ENDIF
setting%output%ndim = ndim2
CS_THRLOG = Int_Input%CS_Thrlog
!write(lupri,*)'THE COMBINED GAB MATRIX'
!call lstensor_print(setting%output%ScreenTensor,lupri)
!call cleanup_gabmatrix(setting%output%ScreenTensor,CS_THRLOG,CS_SCREEN,PS_SCREEN,lupri)
!write(lupri,*)'THE CLEANED UP GAB MATRIX'
!call lstensor_print(setting%output%ScreenTensor,lupri)
END SUBROUTINE ls_getScreenIntegrals1

subroutine ls_init_screenlstensor(fullGAB,AO1,AO2,intType,setting,lupri,luerr)
implicit none
type(lstensor)       :: fullgab
integer              :: AO1,AO2,intType
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
LOGICAL              :: LHSGAB
LHSGAB=.TRUE.
CALL init_integral_input(INT_INPUT,SETTING)
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
IF(intType.EQ.Primitiveinttype)THEN
   Int_Input%PS_int = .TRUE.
   call init_ps_lstensor(fullGAB,Int_Input%AO(1)%p,&
        &Int_Input%AO(2)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),INT_INPUT%OD_SCREEN,lupri)
ELSE
   call init_cs_lstensor(fullGAB,Int_Input%AO(1)%p,&
        &Int_Input%AO(2)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),INT_INPUT%OD_SCREEN,lupri)
ENDIF
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
end subroutine ls_init_screenlstensor

!> \brief calculate the nuclear screening integrals
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param Integrals the output matrix of nuclear screening integrals
!> \param Filename the filename used to store the screening matrix
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1/3
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2/4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param OUTPUT_INPUT the integral input 
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param LHSGAB if this is a LHS screening matrix
SUBROUTINE ls_getNucScreenIntegrals(GAB,Integrals,Filename,AO1,AO2,Oper,SETTING,LUPRI,LUERR,LHSGAB)
implicit none
type(lstensor)      :: GAB
Real(realk)          :: Integrals(:,:)
Character*(*)        :: Filename
integer              :: AO1,AO2,Oper
Type(LSSETTING)      :: SETTING
Integer              :: LUPRI,LUERR
LOGICAL              :: LHSGAB
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(INTEGRALOUTPUT) :: INT_OUTPUT
Integer              :: AOA,AOB
Real(realk)          :: tstart,tend
Integer              :: iAtom,iFragment,iAO,natoms
TYPE(MOLECULE_PT),pointer   :: FRAGMENTS(:)
!CALL ls_setDefaultFragments(setting)
NULLIFY(FRAGMENTS)
ALLOCATE(FRAGMENTS(4))
IF(LHSGAB)THEN
   AOA=1
   AOB=2
ELSE
   AOA=3
   AOB=4
ENDIF
CALL LSTIMER('START ',tstart,tend,lupri)
DO iAO=1,setting%nAO
   FRAGMENTS(iAO)%p => SETTING%FRAGMENT(iAO)%p
ENDDO
IF (LHSGAB) THEN
   FRAGMENTS(3)%p => FRAGMENTS(1)%p
   FRAGMENTS(4)%p => FRAGMENTS(2)%p
   IF (AO1.EQ.AOEmpty) THEN
      iFragment = 2
   ELSEIF (AO2.EQ.AOEmpty) THEN
      iFragment = 1
   ELSE
      CALL LSQUIT('Error in ls_getNucScreenIntegrals, FRAGMENT,LHSGAB',lupri)
   ENDIF
ELSE
   FRAGMENTS(1)%p => FRAGMENTS(3)%p
   FRAGMENTS(2)%p => FRAGMENTS(4)%p
   IF (AO1.EQ.AOEmpty) THEN
      iFragment = 4
   ELSEIF (AO2.EQ.AOEmpty) THEN
      iFragment = 3
   ELSE
      CALL LSQUIT('Error in ls_getNucScreenIntegrals, FRAGMENT,RHSGAB',lupri)
   ENDIF
ENDIF
DO iAtom=1,FRAGMENTS(iFragment)%p%nAtoms
   Integrals(iAtom,1) = ABS(FRAGMENTS(iFragment)%p%ATOM(iAtom)%Charge)
ENDDO
nAtoms = FRAGMENTS(iFragment)%p%nAtoms
call build_Nuclearlstensor(Integrals,GAB,nAtoms)

CALL LSTIMER(Filename(1:15),tstart,tend,lupri)
DEALLOCATE(FRAGMENTS)

END SUBROUTINE ls_getNucScreenIntegrals

!> \brief Sets up the input AOs
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param INT_INPUT the integral input 
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param intType the label for primitive or contracted calc
!> \param AObuild the list of AOITEMS to be build
!> \param nAObuilds the number of AOITEMS build
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param LHSGAB if this is a LHS screening matrix
SUBROUTINE SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,LHSGAB)
implicit none
integer              :: AO1,AO2,AO3,AO4,intType
Type(LSSETTING)      :: SETTING
TYPE(INTEGRALINPUT)  :: INT_INPUT
Integer              :: LUPRI,LUERR,nAObuilds
TYPE(AOITEM),target  :: AObuild(4)
TYPE(MOLECULE_PT),pointer  :: FRAGMENTS(:)
LOGICAL              :: LHSGAB !THIS ONLY HAS AN EFFECT WHEN USING FRAGMENTS 
!
Integer                    :: IAO,JAO,indAO
integer                    :: AOstring(4)
TYPE(BASIS_PT)             :: AObasis(4)
Logical                    :: uniqueAO,emptyAO,intnrm,sameFrag(4,4)
Integer                    :: ndim(4),indexUnique(4),AObatchdim,batchindex(4),batchsize(4)
IF (setting%nAO.NE. 4) CALL LSQUIT('Error in SetInputAO. nAO .NE. 4',lupri)
NULLIFY(FRAGMENTS)
ALLOCATE(FRAGMENTS(4))

!Used for testing if the AOs are identical (depends also in the AO-strings)
sameFrag = setting%sameFrag
batchindex = SETTING%batchindex
batchsize = SETTING%batchsize

AOstring(1) = AO1
AOstring(2) = AO2
AOstring(3) = AO3
AOstring(4) = AO4

DO iAO=1,setting%nAO
  FRAGMENTS(iAO)%p => SETTING%FRAGMENT(iAO)%p
ENDDO
DO iAO=1,setting%nAO
   AObasis(iAO)%p => Setting%BASIS(iAO)%p
ENDDO
INT_INPUT%sameLHSaos = (AO1.EQ.AO2) .AND. (.NOT. AO1.EQ.AOEmpty).AND.samefrag(1,2)
INT_INPUT%sameRHSaos = (AO3.EQ.AO4) .AND. (.NOT. AO3.EQ.AOEmpty).AND.samefrag(3,4)
INT_INPUT%sameODs    = (AO1.EQ.AO3) .AND. (AO2.EQ.AO4).AND.samefrag(1,3).AND.samefrag(2,4)
!Specical settings for Cauchy-Schwarz screening integrals
IF(INT_INPUT%CS_int)THEN
   IF(LHSGAB)THEN
      FRAGMENTS(3)%p => FRAGMENTS(1)%p
      FRAGMENTS(4)%p => FRAGMENTS(2)%p
      AObasis(3)%p => Setting%BASIS(1)%p
      AObasis(4)%p => Setting%BASIS(2)%p
      AOstring(3) = AO1
      AOstring(4) = AO2
      sameFrag(3,4) = sameFrag(1,2)
      sameFrag(4,3) = sameFrag(2,1)
      INT_INPUT%sameRHSaos = INT_INPUT%sameLHSaos
      batchindex(3) = SETTING%batchindex(1) 
      batchindex(4) = SETTING%batchindex(2) 
      IF(batchindex(3).NE.0.OR.batchindex(4).NE.0)THEN
         batchsize(1) = SETTING%batchsize(1) 
         batchsize(2) = SETTING%batchsize(2) 
         batchsize(3) = SETTING%batchsize(1) 
         batchsize(4) = SETTING%batchsize(2) 
      ENDIF
   ELSE
      FRAGMENTS(1)%p => FRAGMENTS(3)%p
      FRAGMENTS(2)%p => FRAGMENTS(4)%p
      AObasis(1)%p => Setting%BASIS(3)%p
      AObasis(2)%p => Setting%BASIS(4)%p
      AOstring(1) = AO3
      AOstring(2) = AO4
      sameFrag(1,2) = sameFrag(3,4)
      sameFrag(2,1) = sameFrag(4,3)
      INT_INPUT%sameLHSaos = INT_INPUT%sameRHSaos
      batchindex(1) = SETTING%batchindex(3) 
      batchindex(2) = SETTING%batchindex(4) 
      IF(batchindex(1).NE.0.OR.batchindex(2).NE.0)THEN
         batchsize(1) = SETTING%batchsize(3) 
         batchsize(2) = SETTING%batchsize(4) 
         batchsize(3) = SETTING%batchsize(3) 
         batchsize(4) = SETTING%batchsize(4) 
      ENDIF
   ENDIF
   sameFrag(1,3) = .TRUE.
   sameFrag(2,4) = .TRUE.
   sameFrag(3,1) = .TRUE.
   sameFrag(4,2) = .TRUE.
   IF(sameFrag(1,2))then
      sameFrag(1,4) = .TRUE.
      sameFrag(2,3) = .TRUE.
      sameFrag(3,2) = .TRUE.
      sameFrag(4,1) = .TRUE.
   ELSE
      sameFrag(1,4) = .FALSE.
      sameFrag(2,3) = .FALSE.
      sameFrag(3,2) = .FALSE.
      sameFrag(4,1) = .FALSE.
   ENDIF
   INT_INPUT%sameODs = .TRUE.
ENDIF

DO iAO=1,setting%nAO
  if(batchindex(iAO).NE. 0)THEN
     FRAGMENTS(iAO)%p => SETTING%MOLECULE(iAO)%p
  endif
ENDDO
!ENDIF
   
IF (intType.EQ.Primitiveinttype) THEN
  intnrm = .TRUE.
ELSEIF (intType.EQ.Contractedinttype) THEN
  intnrm = .FALSE.
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Wrong case in SetInputAO, intType =',intType
  CALL LSQUIT('Error - wrong intType in SetInputAO',lupri)
ENDIF

!Simen Should be set another place!
IF (AOstring(3).NE.AORdefault.OR.AOstring(4).NE.AORdefault) INT_INPUT%CoulombFactor = 1.0E0_realk

nAObuilds = 0
DO IAO=1,4
   uniqueAO = .true.
   DO JAO=1,IAO-1
      IF ((AOstring(IAO).EQ.AOstring(JAO)).AND.sameFrag(iAO,jAO)) THEN
         uniqueAO = .false.
         indAO=indexUnique(JAO)
        EXIT
     ENDIF
  ENDDO
  IF (SETTING%SCHEME%FRAGMENT) uniqueAO = .true.
  IF (uniqueAO) THEN
    nAObuilds = nAObuilds+1
    indAO   = nAObuilds
    indexUnique(IAO) = indAO
    CALL SetAObatch(AObuild(indAO),batchindex(iAO),batchsize(iAO),ndim(indAO),AOstring(iAO),intType,&
     &              SETTING%Scheme,FRAGMENTS(iAO)%p,AObasis(iAO)%p,LUPRI,LUERR)
    IF (AOstring(IAO).EQ.AOpCharge) THEN
      INT_INPUT%sameLHSaos = .FALSE.
      INT_INPUT%sameRHSaos = .FALSE.
    ENDIF
  ENDIF
  Int_Input%AO(IAO)%p => AObuild(indAO)
  Int_Input%AOdim(IAO) = ndim(indAO)
ENDDO

DEALLOCATE(FRAGMENTS)

END SUBROUTINE SetInputAO

!> \brief Sets up the AO-batch
!> \author S. Reine
!> \date 18-03-2010
!> \param AObatch The AO-batch
!> \param AO Specifying what basis set to use: 'Regular', 'DF-Aux', ...
!> \param intType Specifying contracted or primitive basis
!> \param AOtype Specifying basis-function type: 'Hermite', 'Cartesian', 'Default'
!> \param Scheme Specifies integral scheme
!> \param Molecule The molecule
!> \param Basis The basis set
!> \param LUPRI Default output unit
!> \param LUERR Deafult error unit
SUBROUTINE SetAObatch(AObatch,batchindex,batchsize,nDim,AO,intType,Scheme,Molecule,Basis,LUPRI,LUERR)
implicit none
TYPE(AOITEM),intent(INOUT)        :: AObatch
Integer,intent(IN)                :: batchindex
Integer,intent(IN)                :: batchsize
Integer,intent(OUT)               :: nDim
integer,intent(IN)                :: AO,intType
Type(LSINTSCHEME),intent(IN)      :: Scheme
Type(MOLECULEINFO),pointer        :: Molecule
Type(BASISINFO),pointer           :: Basis
Integer,intent(IN)                :: LUPRI,LUERR
!
TYPE(BASISSETINFO),pointer :: AObasis,AObasis2
Logical :: uncont,intnrm,emptyAO,AddBasis2
integer :: AObatchdim,iATOM,nAObatches,batchindex2
Character(len=8)     :: AOstring
uncont = scheme%uncont
IF (intType.EQ.Primitiveinttype) THEN
  intnrm = .TRUE.
ELSEIF (intType.EQ.Contractedinttype) THEN
  intnrm = .FALSE.
ELSE
  WRITE(LUPRI,'(1X,2A)') 'Wrong case in SetAObatch, intType =',intType
  CALL LSQUIT('Error - wrong intType in SetAObatch',lupri)
ENDIF

AddBasis2 = .FALSE.
emptyAO = .false.
SELECT CASE(AO)
CASE (AORegular)
   AObasis => Basis%BINFO(RegBasParam)  !Regular Basis
CASE (AOdfAux)
   AObasis => Basis%BINFO(AUXBasParam)  !AUXILIARY Basis
CASE (AOdfCABS)
   AObasis => Basis%BINFO(RegBasParam)  !Regular Basis
   AObasis2 => Basis%BINFO(CABBasParam)  !CABS Basis
   AddBasis2 = .TRUE.
CASE (AOdfCABO)
   AObasis => Basis%BINFO(CABBasParam)  !CABS Basis only 
CASE (AOdfJK)
   AObasis => Basis%BINFO(JKBasParam)   !JK Basis
CASE (AOVAL)
   AObasis => Basis%BINFO(VALBasParam)  !VALENCE Basis
CASE (AOadmm)
   AObasis => Basis%BINFO(ADMBasParam)  !ADMM Basis
CASE (AOEmpty)
   emptyAO = .true.
   CALL BUILD_EMPTY_AO(AObatch,LUPRI)
   nDim = 1
CASE (AONuclear)
   emptyAO = .true.
   CALL BUILD_EMPTY_NUCLEAR_AO(AObatch,Molecule,LUPRI)
   nDim = 1
CASE (AONuclearSpec)
   !specific nuclei 
   emptyAO = .true.
   IATOM = scheme%AONuclearSpecID 
   CALL BUILD_EMPTY_SINGLE_NUCLEAR_AO(AObatch,Molecule,LUPRI,IATOM)
   nDim = 1
CASE (AOpCharge)
   emptyAO = .true.
   CALL BUILD_EMPTY_PCHARGE_AO(AObatch,Molecule,LUPRI)
   nDim = Molecule%nAtoms
CASE (AOelField)
   emptyAO = .true.
   CALL BUILD_EMPTY_ELFIELD_AO(AObatch,Molecule,LUPRI)
   nDim = 3
CASE DEFAULT
   print*,'Programming error: Not a case in SetAObatch: case: ',AO
   WRITE(lupri,*) 'Programming error: Not a case in SetAObatch: case: ',AO
   call param_AO_Stringfromparam(AOstring,AO)
   WRITE(luerr,*) 'Programming error: Not a case in SetAObatch: case: ',AOstring
   CALL LSQuit('Programming error: Not a case in SetAObatch!',lupri)
END SELECT
IF (.not.emptyAO) THEN
   IF(batchindex.EQ. 0)THEN
      CALL BUILD_AO(LUPRI,SCHEME,SCHEME%AOPRINT,&
           &              Molecule,AObasis,AObatch,&
           &              uncont,intnrm,.FALSE.)
      IF(AddBasis2)THEN
         CALL BUILD_AO(LUPRI,SCHEME,SCHEME%AOPRINT,&
              &              Molecule,AObasis2,AObatch,&
              &              uncont,intnrm,AddBasis2)
      ENDIF
      nDim = getNbasis(AO,intType,Molecule,LUPRI)
   ELSE
      IF(AddBasis2)THEN
         Call determinenAObatches(nAObatches,LUPRI,SCHEME,&
              & SCHEME%AOPRINT,molecule,AObasis,uncont,intnrm)
         !print*,'batchindex',batchindex,'nAObatches',nAObatches
         IF(batchindex.LE.nAObatches)THEN
            CALL BUILD_SHELLBATCH_AO(LUPRI,SCHEME,&
                 & SCHEME%AOPRINT,molecule,AObasis,AObatch,&
                 & uncont,intnrm,batchindex,AObatchdim,batchsize)
         ELSE
            batchindex2=batchindex-nAObatches
            !print*,'batchindex2',batchindex2,'nAObatches',nAObatches
            CALL BUILD_SHELLBATCH_AO(LUPRI,SCHEME,&
                 & SCHEME%AOPRINT,molecule,AObasis2,AObatch,&
                 & uncont,intnrm,batchindex2,AObatchdim,batchsize)
         ENDIF
      ELSE
         CALL BUILD_SHELLBATCH_AO(LUPRI,SCHEME,&
              & SCHEME%AOPRINT,molecule,AObasis,AObatch,&
              & uncont,intnrm,batchindex,AObatchdim,batchsize)
      ENDIF
      nDim = AObatchdim
   ENDIF
ENDIF

END SUBROUTINE SetAObatch

!> \brief frees the input AOs
!> \author T. Kjaergaard and S. Reine
!> \date 2010
!> \param AObuild the list of AOITEMS to be build
!> \param nAObuilds the number of AOITEMS build
!> \param LUPRI logical unit number of the Default output file
SUBROUTINE FreeInputAO(AObuild,nAObuilds,LUPRI)
implicit none
Integer              :: LUPRI,nAObuilds
TYPE(AOITEM),target  :: AObuild(4)
!
integer :: iAO

DO iAO=1,nAObuilds
  CALL free_aoitem(lupri,AObuild(iAO))
ENDDO
 
END SUBROUTINE FreeInputAO

!> \brief calculates the multipole moments required for FMM
!> \author T. Kjaergaard (modified by A. Krapp)
!> \date 2010
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param SETTING Integral evalualtion settings
!> \param nbast the number of basis functions 
!> \param nbastaux the number of auxillary basis functions 
!> \param ndim1 the size of the 1. dimension 
!> \param ndim2 the size of the 2. dimension 
!> \param ndim3 the size of the 3. dimension 
!> \param ndim4 the size of the 4. dimension 
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper the label for the operator like ('Coulomb','Overlap','Kinetic')
!> \param intType the label for primitive or contracted calc
!> \param COULOMB if this is a coulomb calc (or overlap)
SUBROUTINE ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,nbastaux,ndim1,ndim2,ndim3,ndim4,AO1,AO2,AO3,AO4,Spec,intType,COULOMB)
IMPLICIT NONE
! input/output variables
INTEGER               :: LUPRI,LUERR,ndim1,ndim2,ndim3,ndim4,nbast,nbastaux
TYPE(LSSETTING)       :: SETTING
integer               :: AO1,AO2,AO3,AO4,intType,spec
LOGICAL               :: COULOMB
! local variables
Real(realk),pointer   :: integrals(:,:,:,:,:)
REAL(REALK)           :: TS,TE
LOGICAL               :: LHSDENSFIT,RHSDENSFIT,FOURcenter
LOGICAL               :: fragment2
TYPE(AOITEM),target   :: AObuild(4)
Integer               :: nderiv,nmat,nsphmat,I
Integer               :: nAObuilds,LU_DENS,LUINTM2
Integer               :: MMunique_ID1,MMunique_ID2,IDUMMY,N,T,J,IDER,start1,start2
Integer               :: nrowLHS,ncolLHS,nrowRHS,ncolRHS
!
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
!
!initialize the Integral Output structure.
call nullifyIntegralOutput(INT_OUTPUT)
IDUMMY=1
LHSDENSFIT = .FALSE.
RHSDENSFIT = .FALSE.
FOURcenter = .FALSE.
IF(AO1 .EQ.Aodfdefault) LHSDENSFIT = .TRUE.
IF(AO3 .EQ.Aodfdefault) RHSDENSFIT = .TRUE.
IF(.NOT.LHSDENSFIT .AND. .NOT.RHSDENSFIT)THEN
   IF(AO1 .EQ. AORdefault .AND. AO2 .EQ. AORdefault)THEN
      IF(AO3 .EQ. AORdefault .AND. AO4 .EQ. AORdefault) FOURcenter = .TRUE.
   ENDIF
ENDIF

IF(Spec.EQ.GradientSpec) THEN
  INT_OUTPUT%DOGRAD = .TRUE.
  SETTING%OUTPUT%DOGRAD = .true.
  SETTING%SCHEME%CREATED_MMFILES=.false.
ELSE
  INT_OUTPUT%DOGRAD = .false.
  SETTING%OUTPUT%DOGRAD = .false.
ENDIF

IF(SETTING%SCHEME%NO_MMFILES)THEN
   WRITE(LUPRI,*)'You have chosen not to print the MM_DATA and the'
   WRITE(LUPRI,*)'MM_CNTS files which means that you take these files from the'
   WRITE(LUPRI,*)'traditionel dalton and not the new Integral-interface'
   IF(Spec.EQ.GradientSpec)THEN
      CALL LSQUIT('For a gradient run you must print the MM_DATA files',LUPRI)
   ENDIF
ELSE

   ! build info file MM_CNT0
   CALL BUILD_MM_CNT0(LUPRI,LHSDENSFIT, RHSDENSFIT)

   IF(.NOT.SETTING%SCHEME%CREATED_MMFILES)THEN
      CALL LSTIMER('START ',TS,TE,LUPRI)
      ! get info if buffer should be used
      INT_OUTPUT%USEBUFMM = SETTING%SCHEME%USEBUFMM
      IF(INT_OUTPUT%USEBUFMM) THEN
         ! IDER SHOULD BE SET TO 2 OR HIGHER FOR GRADIENT AND HIGHER DERIVATIVES 
         ! (higher than IDER=2 not yet supported)
         IDER = 1
         IF(INT_OUTPUT%DOGRAD) IDER = 2
         ! initialize buffer, allocate memory
         call LS_INITMMBUF(INT_OUTPUT,IDER)
      END IF

      ! open files for multipole moments : MM_DATA, MM_DATR, MM_DADA, MM_DADR
      CALL OPEN_MM_DATA(SETTING,INT_OUTPUT,LUPRI)

      MMunique_ID1 = 0
      MMunique_ID2 = 0
      start1 = 0
      start2 = 0

      ! calculate moments
      IF(FOURCENTER)THEN
       CALL MM_calculation(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
      ELSEIF(RHSDENSFIT)THEN
       CALL MM_calculation(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
       CALL MM_calculation(AO3,AO4,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
      ELSEIF(LHSDENSFIT)THEN
       CALL MM_calculation(AO3,AO4,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
       CALL MM_calculation(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
      ENDIF

      ! empty buffer if used, write stop signals, write nuclear info to MM data files
      CALL WRITE_FINAL_MM_DATA(INT_OUTPUT,SETTING)

      ! close MM data files
      CALL CLOSE_MM_DATA(SETTING,INT_OUTPUT)

      ! timer
      IF(FOURCENTER)THEN
       CALL LSTIMER('MULMOM',TS,TE,LUPRI)
      ELSEIF(RHSDENSFIT)THEN
       CALL LSTIMER('MULMOM-RHSDENS',TS,TE,LUPRI)
      ELSEIF(LHSDENSFIT)THEN
       CALL LSTIMER('MULMOM-LHSDENS',TS,TE,LUPRI)
      ENDIF

      ! build MM_CNTS info file
      CALL BUILD_MM_CNTS(LUPRI,SETTING%SCHEME%MM_LMAX,NBAST,SETTING%SCHEME%MMunique_ID1,&
         &SETTING%MOLECULE(1)%p%NATOMS,IDUMMY,nbastaux,LHSDENSFIT,RHSDENSFIT,INT_OUTPUT%DOGRAD)

      IF(SETTING%SCHEME%CFG_LSDALTON .AND. COULOMB .AND. .NOT. INT_OUTPUT%DOGRAD) THEN
       SETTING%SCHEME%CREATED_MMFILES = .TRUE.
      ELSE
       !WE HAVE TO OVERWRITE EVERY TIME BECAUSE "FCK3" would overwrite the file
      ENDIF

   ELSE
      ! build MM_CNTS info file
      CALL BUILD_MM_CNTS(LUPRI,SETTING%SCHEME%MM_LMAX,NBAST,SETTING%SCHEME%MMunique_ID1,&
       &SETTING%MOLECULE(1)%p%NATOMS,IDUMMY,nbastaux,LHSDENSFIT,RHSDENSFIT,INT_OUTPUT%DOGRAD)
   ENDIF

    ! write density to file
   CALL LS_WRITE_MM_DENS(SETTING,ndim1,ndim2,ndim3,ndim4,nbast,nbastaux,FOURCENTER,RHSDENSFIT,LHSDENSFIT,LUPRI)

ENDIF

END SUBROUTINE ls_multipolemoment
      

!> \brief opens the file to write the multipoles, which is then read by the FMM driver.
!> \author T. Kjaergaard
!> \date 2010
!> \param LUINTM logical unit number of the MM file
!> \param FILENAME the name of the MM file
!> \param LUPRI logical unit number of the Default output file
SUBROUTINE OPENMMFILE(LUINTM,FILENAME,LUPRI)
Character*(*)   :: FILENAME
INTEGER         :: LUINTM!logic unit number
INTEGER         :: LUPRI,IDUMMY
LOGICAL         :: fileexist  

INQUIRE(file=FILENAME,EXIST=fileexist) 
IF(fileexist)THEN
   CALL lsOPEN(LUINTM,FILENAME,'UNKNOWN','UNFORMATTED')
   call lsclose(LUINTM,'DELETE')
   INQUIRE(file=FILENAME,EXIST=fileexist) 
   IF(fileexist)THEN
      WRITE(LUPRI,*)'ERROR',FILENAME ,' not deleted in ls_getmultipolemom'
      WRITE(*,*)'ERROR',FILENAME ,' not deleted in ls_getmultipolemom'
      CALL LSQUIT('ERROR'//FILENAME //' not deleted in ls_getmultipolemom',lupri)
   ELSE
!     WRITE(LUPRI,*)FILENAME ,' exist so we overwrite'
      CALL lsOPEN(LUINTM,FILENAME,'UNKNOWN','UNFORMATTED')
   ENDIF
ELSE
   WRITE(LUPRI,*)FILENAME ,' do not exist so we make a one'
   CALL lsOPEN(LUINTM,FILENAME,'UNKNOWN','UNFORMATTED')
ENDIF

END SUBROUTINE OPENMMFILE

!> \brief wrapper to write density to file, which is then read by the FMM driver.
!> \author T. Kjaergaard (modified by A. Krapp)
!> \date 2010
SUBROUTINE LS_WRITE_MM_DENS(SETTING,ndim1,ndim2,ndim3,ndim4,nbast,nbastaux,FOURCENTER,RHSDENSFIT,LHSDENSFIT,LUPRI)
IMPLICIT NONE
INTEGER               :: LUPRI,ndim1,ndim2,ndim3,ndim4,nbast,nbastaux
TYPE(LSSETTING)       :: SETTING
Real(realk),pointer   :: Dfull(:,:,:),RHSDENS(:,:,:),DmatRHS(:,:,:),DmatLHS(:,:,:)
LOGICAL               :: LHSDENSFIT,RHSDENSFIT,FOURCENTER
!
Integer               :: nrowLHS,ncolLHS,nrowRHS,ncolRHS

IF(SETTING%LHSdmat)THEN
   nrowLHS = SETTING%DmatLHS(1)%p%nrow
   ncolLHS = SETTING%DmatLHS(1)%p%ncol
   call mem_alloc(Dfull,nrowLHS,nrowLHS,1)
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      CALL DCOPY(nrowLHS*ncolLHS,SETTING%DmatLHS(1)%p%elms,1,Dfull(:,:,1),1)
      CALL DAXPY(nrowLHS*ncolLHS,1E0_realk,SETTING%DmatLHS(1)%p%elmsb,1,Dfull(:,:,1),1)
      CALL DSCAL(nrowLHS*ncolLHS,0.5E0_realk,Dfull(:,:,1),1)
   ELSE
      call mat_to_full(SETTING%DmatLHS(1)%p,1E0_realk,Dfull(:,:,1))
   ENDIF
ELSEIF(SETTING%LHSdfull)THEN
   nrowLHS = ndim1
   ncolLHS = ndim2
   call mem_alloc(Dfull,nrowLHS,nrowLHS,1)
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      call lsquit('this use of ls_multipole have not been tested',-1)
      CALL DCOPY(nrowLHS*ncolLHS,SETTING%DfullLHS,1,Dfull(:,:,1),1)
      CALL DSCAL(nrowLHS*ncolLHS,0.5E0_realk,Dfull(:,:,1),1)
   ELSE
      CALL DCOPY(nrowLHS*ncolLHS,SETTING%DfullLHS,1,Dfull(:,:,1),1)
   ENDIF
ELSE
      call lsquit('LHS density matrix not attached to setting in ls_multipolemoment',-1)
ENDIF
DmatLHS => Dfull

IF(SETTING%RHSdmat)THEN
   nrowRHS = SETTING%DmatRHS(1)%p%nrow
   ncolRHS = SETTING%DmatRHS(1)%p%ncol
   call mem_alloc(RHSDENS,nrowRHS,ncolRHS,1)
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      CALL DCOPY(nrowRHS*ncolRHS,SETTING%DmatRHS(1)%p%elms,1,RHSDENS(:,:,1),1)
      CALL DAXPY(nrowRHS*ncolRHS,1E0_realk,SETTING%DmatRHS(1)%p%elmsb,1,RHSDENS(:,:,1),1)
      CALL DSCAL(nrowRHS*ncolRHS,0.5E0_realk,RHSDENS(:,:,1),1)
   ELSE
      call mat_to_full(SETTING%DmatRHS(1)%p,1E0_realk,RHSDENS(:,:,1))
   ENDIF
ELSEIF(SETTING%RHSdfull)THEN
   nrowRHS = ndim3
   ncolRHS = ndim4
   call mem_alloc(RHSDENS,nrowRHS,nrowRHS,1)
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      call lsquit('this use of ls_multipole have not been tested',-1)
      CALL DCOPY(nrowRHS*ncolRHS,SETTING%DfullRHS,1,RHSDENS,1)
      CALL DSCAL(nrowRHS*ncolRHS,0.5E0_realk,RHSDENS,1)
   ELSE
      call dcopy(nrowRHS*ncolRHS,SETTING%DfullRHS,1,RHSDENS,1)
   ENDIF
ELSE
      call lsquit('RHS density matrix not attached to setting in ls_multipolemoment',-1)
ENDIF
DmatRHS => RHSDENS

CALL LS_WRITE_MM_DENS1(NBAST,NBASTAUX,DmatLHS(:,:,1),NROWLHS,NCOLLHS,DmatRHS(:,:,1),&
     &NROWRHS,NCOLRHS,FOURCENTER,RHSDENSFIT,LHSDENSFIT,LUPRI)

IF(SETTING%LHSdmat .OR. SETTING%LHSdfull) call mem_dealloc(DmatLHS)
IF(SETTING%RHSdmat .OR. SETTING%RHSdfull) call mem_dealloc(DmatRHS)

END SUBROUTINE  LS_WRITE_MM_DENS

!> \brief write density to file, which is then read by the FMM driver.
!> \author T. Kjaergaard (modified by A. Krapp)
!> \date 2010
SUBROUTINE LS_WRITE_MM_DENS1(NBAST,NBASTAUX,DFULL,DIM1,DIM2,RHSDENS,DIM21,DIM22,FOURCENTER,RHSDENSFIT,LHSDENSFIT,LUPRI)
IMPLICIT NONE
!
INTEGER                :: NBAST,NBASTAUX,LUPRI,DIM1,DIM2,DIM21,DIM22
Real(realk)            :: Dfull(DIM1,DIM2),RHSDENS(DIM21,DIM22)
LOGICAL                :: LHSDENSFIT,RHSDENSFIT,FOURCENTER
!
INTEGER                :: LU_DENS, N, T, J, I
REAL(realk), pointer   :: WORK(:)
real(realk), parameter :: TWO=2E0_realk

 ! open file
LU_DENS = -1
IF(FOURCENTER)THEN
   CALL OPENMMFILE(LU_DENS,'MM_DENS',LUPRI)
ELSEIF(RHSDENSFIT) THEN
   CALL OPENMMFILE(LU_DENS,'MM_DENR',LUPRI)
ELSEIF(LHSDENSFIT) THEN
   CALL OPENMMFILE(LU_DENS,'MM_DENL',LUPRI)
ELSE
   CALL LSQUIT('ERROR IN LS_WRITE_MM_DENS',LUPRI)
ENDIF

 ! write to file 
N = (NBAST*(NBAST-1))/2+NBAST
call mem_alloc(WORK,N)
T = 0
DO I = 1, NBAST
   DO J = 1, I
      WORK(T+J) = TWO*DFULL(I,J)+TWO*DFULL(J,I)
   END DO
   T = T + I
   WORK(I*(I-1)/2+I)=TWO*DFULL(I,I)
END DO
REWIND(LU_DENS)
WRITE(LU_DENS) NBAST
CALL ls_write(LU_DENS,WORK,N)
call mem_dealloc(WORK)

 ! close file 
call lsclose(LU_DENS,'KEEP')

 ! fitting density
IF(RHSDENSFIT)THEN
   CALL OPENMMFILE(LU_DENS,'MM_FRDN',LUPRI)
   REWIND(LU_DENS)
   WRITE(LU_DENS) NBASTAUX
   CALL ls_write(LU_DENS,RHSDENS,NBASTAUX)
   call lsclose(LU_DENS,'KEEP')
ELSEIF(LHSDENSFIT)THEN
   CALL OPENMMFILE(LU_DENS,'MM_FLDN',LUPRI)
   REWIND(LU_DENS)
   WRITE(LU_DENS) NBASTAUX
   CALL ls_write(LU_DENS,RHSDENS,NBASTAUX)
   call lsclose(LU_DENS,'KEEP')
ENDIF

END SUBROUTINE LS_WRITE_MM_DENS1

!> \brief wrapper for the calculation of the multipoles moments used by FMM
!> \author T. Kjaergaard
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param start1 the start index for the 1 dimension
!> \param start2 the start index for the 2 dimension
!> \param MMunique_ID1 a unique identifier used in the file,required by FMM code
!> \param MMunique_ID2 a unique identifier used in the file,required by FMM code
!> \param integral_output the output specifications of the integral eval
SUBROUTINE MM_calculation(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
integer               :: AO1,AO2,intType
Integer               :: start1,start2,MMunique_ID1,MMunique_ID2
Integer               :: s1,s2,I,J,inode
logical               :: samefragment,Ldummy

Call MM_kernel(AO1,AO2,intType,SETTING,LUPRI,LUERR,start1,start2,&
  &            MMunique_ID1,MMunique_ID2,INT_OUTPUT)
SETTING%SCHEME%MMunique_ID1 = MMunique_ID1

end subroutine MM_calculation

!> \brief calculate the multipoles and write them to a file, which is then read by the FMM driver.
!> \author T. Kjaergaard
!> \date 2010
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param intType the label for primitive or contracted calc
!> \param SETTING Integral evalualtion settings
!> \param LUPRI logical unit number of the Default output file
!> \param LUERR logical unit number of the Default error file
!> \param sA the start index for the 1 dimension
!> \param sB the start index for the 2 dimension
!> \param MMunique_ID1 a unique identifier used in the file,required by FMM code
!> \param MMunique_ID2 a unique identifier used in the file,required by FMM code
!> \param integral_output the output specifications of the integral eval
subroutine MM_kernel(AO1,AO2,intType,SETTING,LUPRI,LUERR,sA,sB,MMunique_ID1,MMunique_ID2,INT_OUTPUT)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR
TYPE(LSSETTING)       :: SETTING
integer               :: AO1,AO2,intType
Real(realk),pointer   :: integrals(:,:,:,:,:)
INTEGER               :: nderiv,nmat,mmorder
TYPE(AOITEM),target   :: AObuild(4)
Integer               :: nAObuilds
Integer               :: MMunique_ID1,MMunique_ID2,sA,sB
TYPE(INTEGRALINPUT)   :: INT_INPUT
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
Logical               :: singleP,emptyP

CALL init_integral_input(INT_INPUT,SETTING)
INT_INPUT%LU_MMDATA = SETTING%SCHEME%LU_LUINTM
IF (INT_OUTPUT%USEBUFMM) INT_INPUT%LU_MMDATR = SETTING%SCHEME%LU_LUINTR
INT_INPUT%operator = MulmomOperator
INT_INPUT%DO_MULMOM  = .TRUE.
INT_INPUT%OD_MOM     = .TRUE.
INT_INPUT%OD_SCREEN = SETTING%SCHEME%OD_SCREEN        
IF(INT_OUTPUT%DOGRAD) INT_INPUT%DO_GRADIENT = .TRUE.
CALL SetInputAO(INT_INPUT,AO1,AO2,AOEmpty,AOEmpty,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
! We screen the multipole moments because if the 
! integral (ab|cd)D_{cd} is not calculated then there is no need to 
! build the multipole moments.
! very generally then we do not know if we have the expansion
! M_{ab}*M_{cd}*D_{cd}
! or 
! M_{ab}*M_{alpha}*C_{alpha}
! so when screening we set GAB_RHS = 1 and build M_{ab} if
! G_{ab} > modifiedThreshold = Threshold*1E-2_realk
! so we set Gab_{RHS} = 1 inside SET_OVERLAP inside MAININTDRIVER
IF(sA .EQ. sB)THEN
   INT_INPUT%sameLHSaos = (AO1.EQ.AO2)
ELSE
   INT_INPUT%sameLHSaos = .FALSE.
ENDIF

Int_Output%ndim = 1

INT_INPUT%MM_NOSCREEN = SETTING%SCHEME%MM_NOSCREEN
MMORDER               = SETTING%SCHEME%MM_LMAX
INT_INPUT%MMORDER     = MMORDER
INT_INPUT%nMultipoleMomentComp = (MMORDER+1)*(MMORDER+2)*(MMORDER+3)/6

INT_INPUT%GEODERORDERP   = 0
INT_INPUT%GEODERORDERQ   = 0
INT_INPUT%NgeoDERIVcompP = 1
INT_INPUT%NgeoDERIVcompQ = 1
INT_INPUT%GEODERIVORDER  = 0

INT_INPUT%sphericalEcoeff = .FALSE.

IF(INT_OUTPUT%DOGRAD) THEN
   INT_INPUT%LU_MMDADA = SETTING%SCHEME%LU_LUINDM
   IF(INT_OUTPUT%USEBUFMM) INT_INPUT%LU_MMDADR = SETTING%SCHEME%LU_LUINDR

   INT_INPUT%GEODERORDERP  = 1
   INT_INPUT%GEODERIVORDER = 1
   singleP                 = ((AO1.EQ.AOEmpty).OR.(AO2.EQ.AOEmpty)).AND.&
     &                        .NOT.((AO1.EQ.AOEmpty).AND.(AO2.EQ.AOEmpty))
   emptyP = (AO1.EQ.AOEmpty).AND.(AO2.EQ.AOEmpty)
   INT_INPUT%NGEODERIVCOMP = getTotalGeoComp(1,.TRUE.,.FALSE.,singleP,.FALSE.,emptyP,.TRUE.)

   IF ((AO1.EQ.AOEmpty).AND.(AO2.EQ.AOEmpty)) THEN
     CALL LSQUIT('Error in MM_Kernel DO_GRADIENT and AO1=AO2="Empty"',-1)
   ELSEIF ((AO1.EQ.AOEmpty).OR.(AO2.EQ.AOEmpty)) THEN
     INT_INPUT%NgeoDERIVcompP = 3
   ELSE
     INT_INPUT%NgeoDERIVcompP = 6
   ENDIF
   SETTING%OUTPUT%DOGRAD = .TRUE.
ENDIF
!
INT_INPUT%MMunique_ID1 = MMunique_ID1
INT_INPUT%MMunique_ID2 = MMunique_ID2
INT_INPUT%MMstartA     = sA
INT_INPUT%MMstartB     = sB

call MAIN_INTEGRAL_DRIVER(LUPRI,SETTING%SCHEME%INTPRINT,INT_INPUT,INT_OUTPUT)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)

MMunique_ID1 = INT_INPUT%MMunique_ID1
MMunique_ID2 = INT_INPUT%MMunique_ID2

end subroutine MM_kernel

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_matarrayp(Dmat,ndmat,setting,side,AOindex1,AOindex2,unres_J,lupri)
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrixp),intent(IN)        :: Dmat(ndmat)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
logical :: unres_J
!
integer                :: idmat,n1,n2,ndmat2
real(realk), pointer :: Dfull(:,:,:)

IF(matrix_type .EQ. mtype_unres_dense)THEN !FIXME BYKOV This should be for all unres types! 
   n1=Dmat(1)%p%nrow
   n2=Dmat(1)%p%ncol
   IF(unres_J)THEN  ! Coulomb type handling: We can add the alpha and beta
      !The alpha and beta part of the coulomb matrix is the same.
      ndmat2 = ndmat
      call mem_alloc(Dfull,n1,n2,ndmat)
      call ls_dzero(Dfull,n1*n2*ndmat)
      DO Idmat =1,ndmat
         call mat_add_to_fullunres(Dmat(idmat)%p,0.5E0_realk,Dfull(:,:,idmat),1) !alpha part
         call mat_add_to_fullunres(Dmat(idmat)%p,0.5E0_realk,Dfull(:,:,idmat),2) !beta part
      ENDDO
   ELSE ! Exchange type handling: We treat alpha and beta seperate
      ndmat2 = 2*ndmat
      call mem_alloc(Dfull,n1,n2,ndmat2)
      call ls_dzero(Dfull,n1*n2*ndmat2)
      DO Idmat =1,ndmat
         call mat_add_to_fullunres(Dmat(idmat)%p,1.0E0_realk,Dfull(:,:,2*idmat-1),1) !alpha part
         call mat_add_to_fullunres(Dmat(idmat)%p,1.0E0_realk,Dfull(:,:,2*idmat),2)   !beta part
      ENDDO
   ENDIF
   call ls_attachDmatToSetting2full(Dfull,n1,n2,ndmat2,setting,side,AOindex1,AOindex2,lupri)
   call mem_dealloc(Dfull)
ELSE
   call ls_attachDmatToSetting2matp(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
ENDIF
END SUBROUTINE ls_attachDmatToSetting_matarrayp

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_matarray(Dmat,ndmat,setting,side,AOindex1,AOindex2,Unres_J,lupri)
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrix),intent(IN),target  :: Dmat(ndmat)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
logical :: unres_J
!
TYPE(matrixp)  :: Dmatp(ndmat)
integer :: idmat,n1,n2,ndmat2
real(realk), pointer :: Dfull(:,:,:)
IF(matrix_type .EQ. mtype_unres_dense)THEN !FIXME BYKOV This should be for all unres types! 
   n1=Dmat(1)%nrow
   n2=Dmat(1)%ncol
   IF(unres_J)THEN  ! Coulomb type handling: We can add the alpha and beta
      !The alpha and beta part of the coulomb matrix is the same.
      ndmat2 = ndmat
      call mem_alloc(Dfull,n1,n2,ndmat)
      call ls_dzero(Dfull,n1*n2*ndmat)
      DO Idmat =1,ndmat
         call mat_add_to_fullunres(Dmat(idmat),0.5E0_realk,Dfull(:,:,idmat),1) !alpha part
         call mat_add_to_fullunres(Dmat(idmat),0.5E0_realk,Dfull(:,:,idmat),2) !beta part
      ENDDO
   ELSE ! Exchange type handling: We treat alpha and beta seperate
      ndmat2 = 2*ndmat
      call mem_alloc(Dfull,n1,n2,ndmat2)
      call ls_dzero(Dfull,n1*n2*ndmat2)
      DO Idmat =1,ndmat
         call mat_add_to_fullunres(Dmat(idmat),1.0E0_realk,Dfull(:,:,2*idmat-1),1) !alpha part
         call mat_add_to_fullunres(Dmat(idmat),1.0E0_realk,Dfull(:,:,2*idmat),2)   !beta part
      ENDDO
   ENDIF
   call ls_attachDmatToSetting2full(Dfull,n1,n2,ndmat2,setting,side,AOindex1,AOindex2,lupri)
   call mem_dealloc(Dfull)
ELSE
   DO idmat=1,ndmat
      Dmatp(idmat)%p => Dmat(idmat)
   ENDDO
   call ls_attachDmatToSetting2matp(Dmatp,ndmat,setting,side,AOindex1,AOindex2,lupri)
ENDIF
END SUBROUTINE ls_attachDmatToSetting_matarray

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_matsinglep(Dmat,ndmat,setting,side,AOindex1,AOindex2,unres_J,lupri)
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrixp),intent(IN),target :: Dmat
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
logical :: unres_J
!
integer                :: idmat,n1,n2,ndmat2
TYPE(matrixp)          :: Dmatp(1)
real(realk) :: thresh
real(realk), pointer :: Dfull(:,:,:)
IF(matrix_type .EQ. mtype_unres_dense)THEN !FIXME BYKOV This should be for all unres types!
   n1=Dmat%p%nrow
   n2=Dmat%p%ncol
   IF(unres_J)THEN  ! Coulomb type handling: We can add the alpha and beta
      !The alpha and beta part of the coulomb matrix is the same.
      ndmat2 = ndmat
      call mem_alloc(Dfull,n1,n2,ndmat)
      call ls_dzero(Dfull,n1*n2*ndmat)
      call mat_add_to_fullunres(Dmat%p,0.5E0_realk,Dfull(:,:,1),1) !alpha part
      call mat_add_to_fullunres(Dmat%p,0.5E0_realk,Dfull(:,:,1),2) !beta part
   ELSE ! Exchange type handling: We treat alpha and beta seperate
      ndmat2 = 2*ndmat
      call mem_alloc(Dfull,n1,n2,ndmat2)
      call ls_dzero(Dfull,n1*n2*ndmat2)
      call mat_add_to_fullunres(Dmat%p,1.0E0_realk,Dfull(:,:,1),1) !alpha part
      call mat_add_to_fullunres(Dmat%p,1.0E0_realk,Dfull(:,:,2),2)   !beta part
   ENDIF
   call ls_attachDmatToSetting2full(Dfull,n1,n2,ndmat2,setting,side,AOindex1,AOindex2,lupri)
   call mem_dealloc(Dfull)
ELSE
   Dmatp(1)%p => Dmat%p
   call ls_attachDmatToSetting2matp(Dmatp,ndmat,setting,side,AOindex1,AOindex2,lupri)
ENDIF
END SUBROUTINE ls_attachDmatToSetting_matsinglep

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_matsingle(Dmat,ndmat,setting,side,AOindex1,AOindex2,unres_J,lupri)
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrix),intent(IN),target  :: Dmat
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
logical :: unres_J
!
TYPE(matrixp)          :: Dmatp(1)
integer                :: idmat,n1,n2,ndmat2
real(realk), pointer :: Dfull(:,:,:)
IF(matrix_type .EQ. mtype_unres_dense)THEN !FIXME BYKOV This should be for all unres types!
   n1=Dmat%nrow
   n2=Dmat%ncol
   IF(unres_J)THEN  ! Coulomb type handling: We can add the alpha and beta
      !The alpha and beta part of the coulomb matrix is the same.
      ndmat2 = ndmat
      call mem_alloc(Dfull,n1,n2,ndmat)
      call ls_dzero(Dfull,n1*n2*ndmat)
      call mat_add_to_fullunres(Dmat,0.5E0_realk,Dfull(:,:,1),1) !alpha part
      call mat_add_to_fullunres(Dmat,0.5E0_realk,Dfull(:,:,1),2) !beta part
   ELSE ! Exchange type handling: We treat alpha and beta seperate
      ndmat2 = 2*ndmat
      call mem_alloc(Dfull,n1,n2,ndmat2)
      call ls_dzero(Dfull,n1*n2*ndmat2)
      call mat_add_to_fullunres(Dmat,1.0E0_realk,Dfull(:,:,1),1) !alpha part
      call mat_add_to_fullunres(Dmat,1.0E0_realk,Dfull(:,:,2),2)   !beta part
   ENDIF
   call ls_attachDmatToSetting2full(Dfull,n1,n2,ndmat2,setting,side,AOindex1,AOindex2,lupri)
   call mem_dealloc(Dfull)
ELSE
   Dmatp(1)%p => Dmat
   call ls_attachDmatToSetting2matp(Dmatp,ndmat,setting,side,AOindex1,AOindex2,lupri)
ENDIF
END SUBROUTINE ls_attachDmatToSetting_matsingle

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting_full(Dmat,dim1,dim2,ndmat,setting,side,AOindex1,AOindex2,lupri)
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2,dim1,dim2
real(realk),intent(IN)          :: Dmat(:,:,:)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
#if VAR_SCALAPACK
IF(setting%scheme%memdist.AND.matrix_type.EQ.mtype_scalapack)then
   call lsquit('programming error use matrix type in attachDmatTosetting',-1)
ENDIF
#endif
call ls_attachDmatToSetting2full(Dmat,dim1,dim2,ndmat,setting,side,AOindex1,AOindex2,lupri)

END SUBROUTINE ls_attachDmatToSetting_full

!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting2matp(Dmat,ndmat,setting,side,AOindex1,AOindex2,lupri)
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2
TYPE(matrixp),intent(IN),target :: Dmat(ndmat)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
integer                :: idmat
SELECT CASE(side)
CASE('LHS')
  IF (setting%LHSdmat) CALL LSQUIT('Error in ls_attachDmatToSetting. LHS',lupri)
  setting%nDmatLHS = ndmat
  setting%LHSdmat  = .TRUE.
  call mem_alloc(setting%DmatLHS,ndmat)
  setting%LHSdmatAlloc=.FALSE.
  do idmat=1,ndmat
     setting%DmatLHS(idmat)%p => Dmat(idmat)%p
  enddo
  call mem_alloc(setting%DsymLHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymLHS(idmat) = mat_get_isym(Dmat(idmat)%p)
  ENDDO
  setting%LHSdmatAOindex1 = AOindex1
  setting%LHSdmatAOindex2 = AOindex2
CASE('RHS')
  IF (setting%RHSdmat) CALL LSQUIT('Error in ls_attachDmatToSetting. RHS',lupri)
  setting%nDmatRHS = ndmat
  setting%RHSdmat  = .TRUE.
  call mem_alloc(setting%DmatRHS,ndmat)
  setting%RHSdmatAlloc=.FALSE.
  do idmat=1,ndmat
     setting%DmatRHS(idmat)%p => Dmat(idmat)%p
  enddo
  call mem_alloc(setting%DsymRHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymRHS(idmat) = mat_get_isym(Dmat(idmat)%p)
  ENDDO
  setting%RHSdmatAOindex1 = AOindex1
  setting%RHSdmatAOindex2 = AOindex2
CASE DEFAULT
  WRITE(LUPRI,'(1X,2A)') 'Error in ls_attachDmatToSetting. Side =',side
  CALL LSQUIT('Error in ls_attachDmatToSetting. Wrong side',lupri)
END SELECT

END SUBROUTINE ls_attachDmatToSetting2matp


!> \brief Attaches a (number of) density matrix(ces) to setting
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param Dmat The density matrices to be assigned
!> \param ndmat The number of density matrices to be assigned
!> \param setting Integral settings
!> \param side Selects which 'side' (LHS or RHS) to contract density matrices
!> \param lupri Default print unit
SUBROUTINE ls_attachDmatToSetting2full(Dmat,dim1,dim2,ndmat,setting,side,AOindex1,AOindex2,lupri)
implicit none
Integer,intent(IN)              :: ndmat,lupri,AOindex1,AOindex2,dim1,dim2
real(realk),intent(IN),target   :: Dmat(:,:,:)
TYPE(LSSETTING),intent(INOUT)   :: setting
character*(*)                   :: side
!
integer                :: idmat

SELECT CASE(side)
CASE('LHS')
  IF (setting%LHSdfull) CALL LSQUIT('Error in ls_attachDmatToSetting. LHS',lupri)
  setting%nDmatLHS = ndmat
  setting%LHSdmat  = .FALSE.
  setting%LHSdfull  = .TRUE.
#if VAR_SCALAPACK
  setting%LHSdAlloc=.TRUE.
  call mem_alloc(setting%DfullLHS,dim1,dim2,ndmat)
  call dcopy(dim1*dim2*ndmat,Dmat,1,setting%DfullLHS,1)
#else
  IF(matrix_type .EQ. mtype_unres_dense.OR.matrix_type .EQ. mtype_pdmm)THEN
     setting%LHSdAlloc=.TRUE.
     call mem_alloc(setting%DfullLHS,dim1,dim2,ndmat)
     call dcopy(dim1*dim2*ndmat,Dmat,1,setting%DfullLHS,1)
  ELSE
     setting%LHSdAlloc=.FALSE.
     setting%DfullLHS => Dmat
  ENDIF
#endif
  setting%LHSdmatAOindex1 = AOindex1
  setting%LHSdmatAOindex2 = AOindex2
  call mem_alloc(setting%DsymLHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymLHS(idmat) = matfull_get_isym(Dmat(:,:,idmat),dim1,dim2)
  ENDDO
CASE('RHS')
  IF (setting%RHSdfull) CALL LSQUIT('Error in ls_attachDmatToSetting. RHS',lupri)
  setting%nDmatRHS = ndmat
  setting%RHSdmat  = .FALSE.
  setting%RHSdfull  = .TRUE.
#if VAR_SCALAPACK
  setting%RHSdAlloc=.TRUE.
  call mem_alloc(setting%DfullRHS,dim1,dim2,ndmat)
  call dcopy(dim1*dim2*ndmat,Dmat,1,setting%DfullRHS,1)
#else
  IF(matrix_type .EQ. mtype_unres_dense.OR.matrix_type .EQ. mtype_pdmm)THEN
     setting%RHSdAlloc=.TRUE.
     call mem_alloc(setting%DfullRHS,dim1,dim2,ndmat)
     call dcopy(dim1*dim2*ndmat,Dmat,1,setting%DfullRHS,1)
  ELSE
     setting%RHSdAlloc=.FALSE.
     setting%DfullRHS => Dmat
  ENDIF
#endif
  setting%RHSdmatAOindex1 = AOindex1
  setting%RHSdmatAOindex2 = AOindex2
  call mem_alloc(setting%DsymRHS,ndmat)
  DO idmat = 1,ndmat
    setting%DsymRHS(idmat) = matfull_get_isym(Dmat(:,:,idmat),dim1,dim2)
  ENDDO
CASE DEFAULT
  WRITE(LUPRI,'(1X,2A)') 'Error in ls_attachDmatToSetting. Side =',side
  CALL LSQUIT('Error in ls_attachDmatToSetting. Wrong side',lupri)
END SELECT

END SUBROUTINE ls_attachDmatToSetting2full

!> \brief Deassosiate density matrices from setting (if they have been associated)
!> \author S.Reine and P.Merlot
!> \date 2010-03-03
!> \param setting Integral settings
SUBROUTINE ls_freeDmatFromSetting(setting)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
!
Integer :: I

IF (setting%LHSdmat) THEN
  IF(setting%LHSdmatAlloc)THEN
     DO I=1,setting%nDmatLHS
        call mat_free(setting%DmatLHS(I)%p)
        deallocate(setting%DmatLHS(I)%p)
        nullify(setting%DmatLHS(I)%p)
     ENDDO
     setting%LHSdmatAlloc=.FALSE.
  ENDIF
  setting%LHSdmat = .FALSE.
  setting%nDmatLHS = 1
  IF (.NOT.ASSOCIATED(setting%DsymLHS)) THEN
     CALL LSQUIT('Error in call to ls_freeDmatFromSetting. DsymLHS',-1)
  ELSE
     call mem_dealloc(setting%DsymLHS)
  ENDIF
  if(associated(setting%DmatLHS))then
     call mem_dealloc(setting%DmatLHS)
  endif
ENDIF
!
IF (setting%RHSdmat) THEN
  IF(setting%RHSdmatAlloc)THEN
     DO I=1,setting%nDmatRHS
        call mat_free(setting%DmatRHS(I)%p)
        deallocate(setting%DmatRHS(I)%p)
        nullify(setting%DmatRHS(I)%p)
     ENDDO
     setting%RHSdmatAlloc=.FALSE.
  ENDIF
  setting%RHSdmat = .FALSE.
  setting%nDmatRHS = 1
  IF (.NOT.ASSOCIATED(setting%DsymRHS)) THEN
     CALL LSQUIT('Error in call to ls_freeDmatFromSetting. DsymRHS',-1)
  ELSE
     call mem_dealloc(setting%DsymRHS)
  ENDIF
  if(associated(setting%DmatRHS))then
     call mem_dealloc(setting%DmatRHS)
  endif
ENDIF
!

if(setting%LHSdalloc)then
   call mem_dealloc(Setting%DfullLHS)
   setting%LHSdalloc=.FALSE.   
endif
IF (setting%LHSdfull) THEN
   setting%LHSdfull = .FALSE.
   setting%nDmatLHS = 1
   nullify(setting%DfullLHS)
ENDIF
IF (ASSOCIATED(setting%DsymLHS)) call mem_dealloc(setting%DsymLHS)
!
if(setting%RHSdalloc)then
   call mem_dealloc(Setting%DfullRHS)
   setting%RHSdalloc=.FALSE.
endif
IF (setting%RHSdfull) THEN
   setting%RHSdfull = .FALSE.
   setting%nDmatRHS = 1
   nullify(setting%DfullRHS)
ENDIF
IF (ASSOCIATED(setting%DsymRHS)) call mem_dealloc(setting%DsymRHS)

setting%LHSdmatAlloc=.FALSE.
setting%RHSdmatAlloc=.FALSE.
setting%LHSdalloc=.FALSE.   
setting%RHSdalloc=.FALSE.   
setting%LHSdmat=.FALSE.   
setting%RHSdmat=.FALSE.   
setting%LHSdfull = .FALSE.
setting%RHSdfull = .FALSE.
END SUBROUTINE ls_freeDmatFromSetting

!> \brief Determines whether the matrices are identical (one-by-one)
!> \author S. Reine
!> \date 2010-04-19
!> \param D First set of matrices
!> \param P Second set of matrices
!> \param nd The number of D matrices
!> \param np The number of P matrices
!> \param ls_same_mats True if D and P matrices are identical
FUNCTION ls_same_mats(D,P,nd,np)
implicit none
Integer,intent(IN)       :: nd,np
TYPE(matrixP),intent(IN) :: D(nd),P(np)
Logical                  :: ls_same_mats
!
logical :: same
integer :: id
real(realk),parameter :: thresh = 1E-12_realk

IF (nd.EQ.np) THEN
  same = .TRUE.
  DO id = 1,nd
    IF (.NOT.mat_same(D(id)%p,P(id)%p,thresh)) same = .FALSE.
  ENDDO
ELSE
  same = .FALSE.
ENDIF

ls_same_mats = same

END FUNCTION ls_same_mats

!> \brief Determines the symmetry of a matrix
!> \author S. Reine
!> \date 2010-07-01
!> \param D Matrix
!> \param ls_mat_sym 1 if symmetric, 2 if anti-symmetric and 3 if no symmetry, 4 zero matrix
FUNCTION ls_mat_sym(D)
implicit none
TYPE(matrix),intent(IN) :: D
Integer                 :: ls_mat_sym

ls_mat_sym = mat_get_isym(D)

END FUNCTION ls_mat_sym

SUBROUTINE set_symmetry(iPQxyz,setting,lupri)
implicit none
integer :: iPQXYZ ,lupri
type(lssetting) :: setting
!
logical :: sameMOL
integer :: A,B,IX,IY,IZ
real(realk) :: DX,DY,DZ,XA,YA,ZA

sameMOL=.TRUE.
DO B=1,setting%nAO
   DO A=1,setting%nAO
      IF(.NOT.setting%sameMOL(A,B))THEN
         sameMOL=.FALSE.
         EXIT
      ENDIF
   ENDDO
ENDDO

IF(sameMOL)THEN
   DX = 0E0_realk
   DY = 0E0_realk
   DZ = 0E0_realk
   XA = setting%molecule(1)%p%ATOM(1)%CENTER(1)
   YA = setting%molecule(1)%p%ATOM(1)%CENTER(2)
   ZA = setting%molecule(1)%p%ATOM(1)%CENTER(3)
   DO B=1,setting%molecule(1)%p%nAtoms
      DX = DX + ABS(XA-setting%molecule(1)%p%ATOM(B)%CENTER(1))
      DY = DY + ABS(YA-setting%molecule(1)%p%ATOM(B)%CENTER(2))
      DZ = DZ + ABS(ZA-setting%molecule(1)%p%ATOM(B)%CENTER(3))
   ENDDO
   IX=1
   IY=1
   IZ=1
   IF(ABS(DX).GT. 1.0E-14_realk)IX=0
   IF(ABS(DY).GT. 1.0E-14_realk)IY=0
   IF(ABS(DZ).GT. 1.0E-14_realk)IZ=0
   iPQxyz = IX + 2*IY + 4*IZ
   !iPQxyz  Meaning
   ! 0      no symmetry   
   IF(setting%scheme%intprint.GT.0)THEN
      IF(setting%node.EQ.0)THEN
         IF(iPQxyz.EQ. 0)WRITE(lupri,*)'SET_SYMMETRY: No Symmetry' 
         IF(iPQxyz.EQ. 1)WRITE(lupri,*)'SET_SYMMETRY: YZ planar molecule' 
         IF(iPQxyz.EQ. 2)WRITE(lupri,*)'SET_SYMMETRY: XZ planar molecule' 
         IF(iPQxyz.EQ. 3)WRITE(lupri,*)'SET_SYMMETRY: Linear molecule along Z axis' 
         IF(iPQxyz.EQ. 4)WRITE(lupri,*)'SET_SYMMETRY: XY planar molecule' 
         IF(iPQxyz.EQ. 5)WRITE(lupri,*)'SET_SYMMETRY: Linear molecule along Y axis' 
         IF(iPQxyz.EQ. 6)WRITE(lupri,*)'SET_SYMMETRY: Linear molecule along X axis' 
         IF(iPQxyz.EQ. 7)WRITE(lupri,*)'SET_SYMMETRY: A single X,Y,Z coordinate' 
      ELSE
         IF(iPQxyz.EQ. 0)WRITE(6,*)'SET_SYMMETRY: No Symmetry' 
         IF(iPQxyz.EQ. 1)WRITE(6,*)'SET_SYMMETRY: YZ planar molecule' 
         IF(iPQxyz.EQ. 2)WRITE(6,*)'SET_SYMMETRY: XZ planar molecule' 
         IF(iPQxyz.EQ. 3)WRITE(6,*)'SET_SYMMETRY: Linear molecule along Z axis' 
         IF(iPQxyz.EQ. 4)WRITE(6,*)'SET_SYMMETRY: XY planar molecule' 
         IF(iPQxyz.EQ. 5)WRITE(6,*)'SET_SYMMETRY: Linear molecule along Y axis' 
         IF(iPQxyz.EQ. 6)WRITE(6,*)'SET_SYMMETRY: Linear molecule along X axis' 
         IF(iPQxyz.EQ. 7)WRITE(6,*)'SET_SYMMETRY: A single X,Y,Z coordinate' 
      ENDIF
   ENDIF
ELSE
! for now assume no symmetry
   iPQxyz = 0
ENDIF

END SUBROUTINE SET_SYMMETRY

!SUBROUTINE ls_init_IntCalcCounter()
!implicit none
!
!call th_init_IntCalcCounter()
!
!END SUBROUTINE ls_init_IntCalcCounter

!SUBROUTINE ls_get_IntCalcCounter(nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO)
!implicit none
!integer :: nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO
!
!call th_get_IntCalcCounter(nCalcIntO,nCalcIntZeroO,nCalcIntZeroContribO)
!
!END SUBROUTINE ls_get_IntCalcCounter

!> \brief Sets the four fragments for point to the four molecules
!> \author S. Reine
!> \date 2010-08-23
!> \param setting The ls-settings
SUBROUTINE ls_setDefaultFragments(setting)
implicit none
TYPE(LSSETTING),intent(inout) :: setting
Integer :: iao
DO iao=1,4
  setting%fragment(iao)%p => setting%molecule(iao)%p
ENDDO
setting%sameFrag = setting%sameMol
END SUBROUTINE ls_setDefaultFragments

!> \brief Open and write to the MM info file on fitting, read in by the FMM routines
!> \author T. Kjaeergaard (modififed by A. Krapp)
!> \date 2010
SUBROUTINE BUILD_MM_CNT0(LUPRI,LHSDENSFIT,RHSDENSFIT)
IMPLICIT NONE
INTEGER LUINTM2, LUPRI
LOGICAL LHSDENSFIT, RHSDENSFIT
!
LUINTM2 = -1
CALL OPENMMFILE(LUINTM2,'MM_CNT0',LUPRI)
WRITE (LUINTM2) LHSDENSFIT, RHSDENSFIT
call lsclose(LUINTM2,'KEEP')
!
END SUBROUTINE BUILD_MM_CNT0

!> \brief opene the MM data file(s), read in by the FMM routines
!> \author T. Kjaeergaard (modififed by A. Krapp)
!> \date 2010
SUBROUTINE OPEN_MM_DATA(SETTING,INT_OUTPUT,LUPRI)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
INTEGER               :: LUPRI
!
IF(INT_OUTPUT%DOGRAD) THEN
   SETTING%SCHEME%LU_LUINDM = -1
   CALL OPENMMFILE(SETTING%SCHEME%LU_LUINDM,'MM_DADA',LUPRI)
   IF (INT_OUTPUT%USEBUFMM)THEN
      SETTING%SCHEME%LU_LUINDR = -1
      CALL OPENMMFILE(SETTING%SCHEME%LU_LUINDR,'MM_DADR',LUPRI)
   ENDIF
ELSE
   SETTING%SCHEME%LU_LUINTM = -1
   CALL OPENMMFILE(SETTING%SCHEME%LU_LUINTM,'MM_DATA',LUPRI)
   IF (INT_OUTPUT%USEBUFMM)THEN
      SETTING%SCHEME%LU_LUINTR = -1
      CALL OPENMMFILE(SETTING%SCHEME%LU_LUINTR,'MM_DATR',LUPRI)
   ENDIF
ENDIF
!
END SUBROUTINE OPEN_MM_DATA

!> \brief write the final part of the MM data file(s), read in by the FMM routines
!> \author T. Kjaeergaard (modififed by A. Krapp)
!> \date 2010
SUBROUTINE WRITE_FINAL_MM_DATA(INT_OUTPUT,SETTING)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
INTEGER               :: I, LU1,LU2
!
IF(INT_OUTPUT%DOGRAD) THEN
   LU1 = SETTING%SCHEME%LU_LUINDM
   LU2 = SETTING%SCHEME%LU_LUINDR
ELSE
   LU1 = SETTING%SCHEME%LU_LUINTM
   LU2 = SETTING%SCHEME%LU_LUINTR
ENDIF
!
IF(INT_OUTPUT%USEBUFMM) THEN
   ! we empty the buffer and write the stop signal
   IF (INT_OUTPUT%IBUFI .GT. 1) THEN
      CALL LS_EMPTYIBUF(INT_OUTPUT,INT_OUTPUT%IBUF,LU1)
      CALL LS_EMPTYRBUF(INT_OUTPUT,INT_OUTPUT%RBUF,LU2)
   END IF
   WRITE(LU1) -1
   WRITE(LU2) -1
   ! we buffer also the saving of the nuclear information
   INT_OUTPUT%IBUFN = 1
   DO I = 1, SETTING%MOLECULE(1)%p%nATOMS
      IF (INT_OUTPUT%IBUFN .GE. INT_OUTPUT%MMBUFLEN-1) THEN
         CALL LS_EMPTYNUCBUF(INT_OUTPUT,INT_OUTPUT%NBUF,LU2)
      END IF
      CALL LS_FILLNUCBUF(INT_OUTPUT,SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE,  &
          &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1),SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2),&
          &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3))
   END DO
   CALL LS_EMPTYNUCBUF(INT_OUTPUT,INT_OUTPUT%NBUF,LU2)
   WRITE(LU2) -2
ELSE
   WRITE(LU1) -1,0,0,0,0,0,0,0,0E0_realk,0E0_realk,0E0_realk,0E0_realk,0E0_realk,0,0,0,0E0_realk
   DO I = 1,SETTING%MOLECULE(1)%p%nATOMS
      WRITE(LU1) SETTING%MOLECULE(1)%p%ATOM(I)%CHARGE,&
            &SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(:)
   ENDDO
ENDIF
!
END SUBROUTINE WRITE_FINAL_MM_DATA

!> \brief close the MM data files, read in by the FMM routines
!> \author T. Kjaeergaard (modififed by A. Krapp)
!> \date 2010
SUBROUTINE CLOSE_MM_DATA(SETTING,INT_OUTPUT)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALOUTPUT)  :: INT_OUTPUT
!
IF(INT_OUTPUT%DOGRAD) THEN
   call lsclose(SETTING%SCHEME%LU_LUINDM,'KEEP')
   IF(INT_OUTPUT%USEBUFMM)  THEN
      call lsclose(SETTING%SCHEME%LU_LUINDR,'KEEP')
      call LS_FREEMMBUF(INT_OUTPUT)
   ENDIF
ELSE
   call lsclose(SETTING%SCHEME%LU_LUINTM,'KEEP')
   IF(INT_OUTPUT%USEBUFMM)  THEN
      call lsclose(SETTING%SCHEME%LU_LUINTR,'KEEP')
      call LS_FREEMMBUF(INT_OUTPUT)
   ENDIF
ENDIF
!
END SUBROUTINE CLOSE_MM_DATA

!> \brief Open and write to the MM information file, read in by the FMM routines
!> \author T. Kjaeergaard (modififed by A. Krapp)
!> \date 2010
SUBROUTINE BUILD_MM_CNTS(LU2,LMAX,NBAST,ID,NATOMS,I,NBASTAUX,LHSFIT,RHSFIT,DO_GRAD)
IMPLICIT NONE
INTEGER LU1,LU2,LMAX,NBAST,NATOMS,NBASTAUX,I,ID
LOGICAL LHSFIT,RHSFIT,DO_GRAD
!
LU1 = -1
IF (DO_GRAD) THEN
   CALL OPENMMFILE(LU1,'MM_CNTD',LU2)
ELSE
   CALL OPENMMFILE(LU1,'MM_CNTS',LU2)
ENDIF
WRITE (LU1) LMAX, NBAST, ID,NATOMS,I,NBASTAUX,LHSFIT,RHSFIT
CALL LSCLOSE(LU1,'KEEP')
!
END SUBROUTINE BUILD_MM_CNTS

subroutine SetScalapackDmatToFull(setting,LHS,RHS)
  implicit none
  Type(LSSETTING)      :: SETTING
  logical              :: LHS,RHS
  !
  integer :: n1,n2,i,ndmat
  !LHS
  IF(LHS)THEN
     IF (setting%LHSdmat)THEN
        n1=setting%DmatLHS(1)%p%nrow
        n2=setting%DmatLHS(1)%p%ncol
        ndmat = setting%nDmatLHS
        setting%LHSdAlloc=.TRUE.
        call mem_alloc(setting%DfullLHS,n1,n2,ndmat)
        DO i=1,ndmat
           call mat_to_full(setting%DmatLHS(i)%p, 1E0_realk,setting%DfullLHS(:,:,i))
        ENDDO
        IF(setting%LHSdmatAlloc)THEN
           DO I=1,setting%nDmatLHS
              call mat_free(setting%DmatLHS(I)%p)
              deallocate(setting%DmatLHS(I)%p)
              nullify(setting%DmatLHS(I)%p)
           ENDDO
           setting%LHSdmatAlloc=.FALSE.
        ENDIF
        if(associated(setting%DmatLHS))then
           call mem_dealloc(setting%DmatLHS)
        endif
        SETTING%LHSdfull = .TRUE.
        setting%LHSdmat = .FALSE.
     ENDIF
  ENDIF
  !RHS
  IF(RHS)THEN
     IF (setting%RHSdmat)THEN
        n1=setting%DmatRHS(1)%p%nrow
        n2=setting%DmatRHS(1)%p%ncol
        ndmat = setting%nDmatRHS
        setting%RHSdAlloc=.TRUE.
        call mem_alloc(setting%DfullRHS,n1,n2,ndmat)
        DO i=1,ndmat
           call mat_to_full(setting%DmatRHS(i)%p, 1E0_realk,setting%DfullRHS(:,:,i))
        ENDDO
        IF(setting%RHSdmatAlloc)THEN
           DO I=1,setting%nDmatRHS
              call mat_free(setting%DmatRHS(I)%p)
              deallocate(setting%DmatRHS(I)%p)
              nullify(setting%DmatRHS(I)%p)
           ENDDO
           setting%RHSdmatAlloc=.FALSE.
        ENDIF
        if(associated(setting%DmatRHS))then
           call mem_dealloc(setting%DmatRHS)
        endif
        SETTING%RHSdfull = .TRUE.
        setting%RHSdmat = .FALSE.
     ENDIF
  ENDIF
end subroutine SetScalapackDmatToFull

SUBROUTINE ls_LHSSameAsRHSDmatToSetting(setting)
implicit none
TYPE(LSSETTING)   :: SETTING
#ifdef VAR_MPI
!due to fragmentation this is not the case. 
setting%LHSSameAsRHSDmat  = .FALSE.
#else
setting%LHSSameAsRHSDmat  = .TRUE.
#endif
END SUBROUTINE LS_LHSSAMEASRHSDMATTOSETTING

SUBROUTINE ls_LHSSameAsRHSDmatToSetting_deactivate(setting)
implicit none
TYPE(LSSETTING)   :: SETTING
setting%LHSSameAsRHSDmat  = .FALSE.
END SUBROUTINE LS_LHSSAMEASRHSDMATTOSETTING_DEACTIVATE

!> \brief set the densities in the integralinput structure correctly
!> \author S. Reine and T. Kjaeergaard
!> \date 2010
!> \param intinput the integral input to be set 
!> \param SETTING Integral evalualtion settings containing the densities
SUBROUTINE ls_setIntegralInputDensities(intinput,setting,lupri,LHS,RHS,LinK,ODscreenInput)
implicit none
integer,intent(in) :: lupri
TYPE(INTEGRALINPUT) :: INTINPUT
TYPE(LSSETTING)   :: SETTING
logical :: LHS,RHS,LinK,ODscreenInput
!
INTEGER :: LHS1,LHS2,RHS1,RHS2,nbastLHS1,nbastLHS2,nbastRHS1,nbastRHS2
INTEGER :: nDmatLHS,nDmatRHS,idmat
LOGICAL :: useAO1,useAO2,TMPalloc,permuteLHS,permuteRHS,ODscreen
Integer :: subdim(4),substart(4),I,n1,n2,n3,n4,nbast(4)
real(realk),pointer :: DmatTMP(:,:,:),DmatTMP2(:,:,:)
real(realk) :: Dthr
Dthr = setting%scheme%THRESHOLD*1.0E-4_realk
IF ((INTINPUT%operator.EQ.KineticOperator) .AND.setting%RHSdfull)&
     &CALL LSQUIT('Error in ls_setintegralinputdensities. Kinetic and RHS density',-1)
!At this point there are several options for both LHS and RHS matrices
IF(RHS)THEN
   IF (setting%RHSdfull.OR.setting%RHSdmat)THEN
      nDmatRHS = setting%nDmatRHS
      INTINPUT%NDMAT_RHS = ndmatRHS
      RHS1 = setting%RHSdmatAOindex1
      RHS2 = setting%RHSdmatAOindex2
      nbastRHS1 = INTINPUT%AOdim(RHS1)
      nbastRHS2 = INTINPUT%AOdim(RHS2)
      useAO1 = .NOT.INTINPUT%AO(RHS1)%p%empty
      useAO2 = .NOT.INTINPUT%AO(RHS2)%p%empty
      ALLOCATE(INTINPUT%LST_DRHS)

      IF(nDmatRHS.GT.size(setting%output%postprocess))THEN
         CAll mem_dealloc(setting%output%postprocess)
         CAll mem_alloc(setting%output%postprocess,nDmatRHS)
         setting%output%postprocess = 0
      ENDIF
   ENDIF
   Tmpalloc = .FALSE.
   IF (setting%RHSdfull) THEN
      !We have attached a full Dmat(:,:,:) in ls_attachDmatToSetting (ls_attachDmatToSetting2full)
      !setting%DfullRHS => Dmat 
      !This means that Setting%DfullRHS is currently pointing to this Density matrix, but is not allocated
      
      !IF setting%RHSdalloc is true, it means that the settings have been copied or something else
      !but it should only affect deallocation, not the features of this subroutine
      
      !we need the full matrix so we build INTINPUT%lst_DRHS directly 
      DmatTMP => Setting%DfullRHS
      !      CALL Build_Dmatlstensor_from_full_3dim(INTINPUT%lst_dRHS,DmatTMP,INTINPUT%AO(RHS1)%p,&
      !           &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,Dthr,lupri)
      IF(INTINPUT%uselst_DRHS)Call lsquit('Error in ls_setIntegralInputDensities RHS lstensor already used ',-1)
      CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dRHS,DmatTMP,INTINPUT%AO(RHS1)%p,&
           &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,lupri)
      INTINPUT%uselst_DRHS = .TRUE.
      IF(Tmpalloc)call mem_dealloc(DmatTMP)         
   ELSEIF(setting%RHSdmat)THEN
      !This is the usual option. We have attached a matrix 
      !pointer Dmat(ndmat) in ls_attachDmatToSetting (ls_attachDmatToSetting2mat(p/single))
      !setting%DmatRHS(idmat)%p => Dmat(idmat)%p
      !This means that Setting%DmatRHS(idmat)%p is currently pointing to these Density matrices, but is not allocated
      IF((INTINPUT%DO_EXCHANGE .AND. INTINPUT%DO_DALINK).AND. matrix_type .EQ. mtype_unres_dense)THEN
         !the unrestricted case is taked care of outside by using the RHSdfull option
         call LSquit(' exchange error unres in ls_setIntegralInputDensities',-1)
      ELSE
         IF((matrix_type .NE. mtype_unres_dense).AND.&
              & (nbastRHS1.EQ.setting%DmatRHS(1)%p%nrow.AND.nbastRHS2.EQ.setting%DmatRHS(1)%p%ncol))THEN
            !coulomb then screen AO
            !FIXME the force triangluar sym matrix thing expects a full non OD screened lstensor
            ODscreen = .FALSE.
            !(ODscreenInput.AND.INTINPUT%sameRHSaos).AND.(((RHS1.EQ.3).AND.(RHS2.EQ.4)).AND.(useAO1.AND.useAO2))
            call Build_lst_from_matarray(INTINPUT%lst_DRHS,setting%DmatRHS,INTINPUT%AO(RHS1)%p,&
                 & INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,ODscreen,lupri)
         ELSE
            call mem_alloc(DmatTMP,nbastRHS1,nbastRHS2,ndmatRHS)
            Tmpalloc = .TRUE.
            IF(matrix_type .EQ. mtype_unres_dense)THEN
               n3 = setting%DmatRHS(1)%p%nrow
               n4 = setting%DmatRHS(1)%p%ncol
               IF(n3.NE.nbastRHS1)CALL LSQUIT(' n3.NE.nbastRHS1 ',lupri)
               IF(n4.NE.nbastRHS2)CALL LSQUIT(' n4.NE.nbastRHS2 ',lupri)
               DO I =1,setting%ndmatRHS
                  CALL DCOPY(n3*n4,setting%DmatRHS(I)%p%elms,1,DmatTMP(:,:,I),1)
                  CALL DAXPY(n3*n4,1E0_realk,setting%DmatRHS(I)%p%elmsb,1,DmatTMP(:,:,I),1)
                  CALL DSCAL(n3*n4,0.5E0_realk,DmatTMP(:,:,I),1)
               ENDDO
            ELSE
               CALL ls_mat_retrive_block(setting%DmatRHS,DmatTMP,ndmatRHS,nbastRHS1,nbastRHS2,1,1,.FALSE.)
            ENDIF
            !         CALL Build_Dmatlstensor_from_full_3dim(INTINPUT%lst_dRHS,DmatTMP,INTINPUT%AO(RHS1)%p,&
            !              &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,Dthr,lupri)
            IF(INTINPUT%uselst_DRHS)Call lsquit('Error in ls_setIntegralInputDensities RHS lstensor already used ',-1)
            CALL Build_lstensor_from_full_3dim(INTINPUT%lst_dRHS,DmatTMP,INTINPUT%AO(RHS1)%p,&
                 &INTINPUT%AO(RHS2)%p,nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,lupri)
         ENDIF
         INTINPUT%uselst_DRHS = .TRUE.
         IF(Tmpalloc)call mem_dealloc(DmatTMP)         
      ENDIF
   ENDIF
ENDIF
!We need to set up the lstensor ( INTINPUT%lst_DLHS )
IF(LHS)THEN
   IF(setting%LHSSameAsRHSDmat)THEN
      IF(INTINPUT%uselst_DRHS)THEN
         nullify(INTINPUT%lst_dLHS)
         INTINPUT%lst_dLHS => INTINPUT%lst_dRHS
         INTINPUT%NDMAT_LHS = INTINPUT%NDMAT_RHS
         INTINPUT%uselst_DLHS = .FALSE.
      ELSE
         call lsquit('LHSSameAsRHSDmat but no RHS set',-1)
      ENDIF
   ELSE
      IF (setting%LHSdfull.OR.setting%LHSdmat)THEN
         ndmatLHS = setting%nDmatLHS
         INTINPUT%NDMAT_LHS = ndmatLHS
         LHS1 = setting%LHSdmatAOindex1
         LHS2 = setting%LHSdmatAOindex2
         nbastLHS1 = INTINPUT%AOdim(LHS1)
         nbastLHS2 = INTINPUT%AOdim(LHS2)
         useAO1 = .NOT.INTINPUT%AO(LHS1)%p%empty
         useAO2 = .NOT.INTINPUT%AO(LHS2)%p%empty
         ALLOCATE(INTINPUT%LST_DLHS)
      ENDIF
      Tmpalloc = .FALSE.
      
      IF (setting%LHSdfull) THEN
         !We have attached a full Dmat(:,:,:) in ls_attachDmatToSetting (ls_attachDmatToSetting2full)
         !setting%DfullLHS => Dmat 
         !This means that Setting%DfullLHS is currently pointing to this Density matrix, but is not allocated
         
         !IF setting%LHSdalloc is true, it means that the settings have been copied or something else
         !but it should only affect deallocation, not the features of this subroutine
         
         !we need the full matrix so we build INTINPUT%lst_DLHS directly 
         DmatTMP => Setting%DfullLHS
         IF(INTINPUT%uselst_DLHS)Call lsquit('Error in ls_setIntegralInputDensities LHS lstensor already used ',-1)
         !      CALL Build_Dmatlstensor_from_full_3dim(INTINPUT%lst_DLHS,DmatTMP,INTINPUT%AO(LHS1)%p,&
         !           &INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,Dthr,lupri)
         CALL Build_lstensor_from_full_3dim(INTINPUT%lst_DLHS,DmatTMP,INTINPUT%AO(LHS1)%p,&
              &INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,lupri)
         INTINPUT%uselst_DLHS = .TRUE.
         IF(Tmpalloc)call mem_dealloc(DmatTMP)
      ELSEIF(setting%LHSdmat)THEN
         IF(matrix_type .EQ. mtype_unres_dense)THEN
            CALL LSQUIT('Error in ls_setIntegralInputDensities. Unrestricted RHSdmat not yet implemented',-1)
         ENDIF
         !This is the usual option. We have attached a matrix 
         !pointer Dmat(ndmat) in ls_attachDmatToSetting (ls_attachDmatToSetting2mat(p/single))
         !setting%DmatLHS(idmat)%p => Dmat(idmat)%p
         !This means that Setting%DmatLHS(idmat)%p is currently pointing to these Density matrices, but is not allocated
         IF((INTINPUT%DO_EXCHANGE .AND. INTINPUT%DO_DALINK).AND. matrix_type .EQ. mtype_unres_dense)THEN
            call LSquit(' exchange error unres in ls_setIntegralInputDensities',-1)
         ELSE
            IF(matrix_type .EQ. mtype_unres_dense)THEN
               call lsquit('matrix_type .EQ. mtype_unres_dense ERROR ',-1)
            ENDIF
            IF(nbastLHS1.EQ.setting%DmatLHS(1)%p%nrow.AND.nbastLHS2.EQ.setting%DmatLHS(1)%p%ncol)THEN
            !FIXME the force triangluar sym matrix thing expects a full non OD screened lstensor
               ODscreen = .FALSE.
               !(ODscreenInput.AND.INTINPUT%sameLHSaos).AND.(((LHS1.EQ.1).AND.(LHS2.EQ.2)).AND.(useAO1.AND.useAO2))
               call Build_lst_from_matarray(INTINPUT%lst_DLHS,setting%DmatLHS,INTINPUT%AO(LHS1)%p,&
                    & INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,ODscreen,lupri)
            ELSE
               call mem_alloc(DmatTMP,nbastLHS1,nbastLHS2,ndmatLHS)
               Tmpalloc = .TRUE.
               CALL ls_mat_retrive_block(setting%DmatLHS,DmatTMP,ndmatLHS,nbastLHS1,nbastLHS2,1,1,.FALSE.)
               IF(INTINPUT%uselst_DLHS)Call lsquit('Error in ls_setIntegralInputDensities LHS lstensor already used ',-1)
               CALL Build_lstensor_from_full_3dim(INTINPUT%lst_DLHS,DmatTMP,INTINPUT%AO(LHS1)%p,&
                    &INTINPUT%AO(LHS2)%p,nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,lupri)
            ENDIF
            INTINPUT%uselst_DLHS = .TRUE.
            IF(Tmpalloc)call mem_dealloc(DmatTMP)
         ENDIF
      ENDIF
   ENDIF
ENDIF
END SUBROUTINE ls_setIntegralInputDensities

SUBROUTINE retrieve_fullblock(BigMat,n1,n2,SmallMat,smalln1,smalln2,s1,s2,ndmat,permute)
implicit none
integer,intent(in) :: n1,n2,ndmat,smalln1,smalln2,s1,s2
real(realk) :: BigMat(n1,n2,ndmat),SmallMat(smalln1,smalln2,ndmat)
logical,intent(in) :: permute
!
integer  :: IDMAT,A,B,start1,start2 

start1=s1-1
start2=s2-1
DO IDMAT=1,ndmat
   DO B=1,smalln2
      DO A=1,smalln1
         SmallMat(A,B,IDMAT) = BigMat(start1+A,start2+B,IDMAT)
      ENDDO
   ENDDO
ENDDO
IF(permute)then
   DO IDMAT=1,ndmat
      DO B=1,smalln2
         DO A=1,smalln1
            SmallMat(A,B,IDMAT) = SmallMat(A,B,IDMAT) + BigMat(start2+B,start1+A,IDMAT)
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE RETRIEVE_FULLBLOCK

!> \brief free the densities in the integralinput structure
!> \author S. Reine and T. Kjaeergaard
!> \date 2010
!> \param intinput the integral input
SUBROUTINE ls_freeIntegralInputDensities(intinput)
implicit none
TYPE(INTEGRALINPUT) :: INTINPUT
!
IF (INTINPUT%uselst_DLHS) THEN
   CALL lstensor_free(INTINPUT%lst_DLHS)
   DEALLOCATE(INTINPUT%LST_DLHS)
   INTINPUT%uselst_DLHS = .FALSE.
ENDIF
IF (INTINPUT%uselst_DRHS)THEN
   CALL lstensor_free(INTINPUT%lst_dRHS)
   DEALLOCATE(INTINPUT%LST_DRHS)
   INTINPUT%uselst_DRHS = .FALSE.
ENDIF

END SUBROUTINE ls_freeIntegralInputDensities

!> \brief Set up lstensor density and result matrices (before parallell loop)
!> \author Simen Reine
!> \date 2011-10-25
!> \param SETTING Integral evalualtion settings
!> \param AO1 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 1
!> \param AO2 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 2
!> \param AO3 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 3
!> \param AO4 the type of the basis ('Regular','DF-Aux',AOEmpty) on center 4
!> \param Oper specifies the operator
!> \param Spec integral specifier ('Regular','Gradient')
!> \param intType the label for primitive or contracted calc
!> \param result_tensor the tensor to store the results
!> \param dmat_lhs the tensor for "lhs" density-matrices (lhs AOs (1,2) for Coulomb, 
!>        AOs (1,3) for exchange, etc.)
!> \param dmat_rhs the tensor for "rhs" density-matrices (rhs AOs (3,4) for Coulomb, 
!>        AOs (2,4) for exchange, etc.)
!> \param lhs_created specifies whether the lhs dmat was created or not
!> \param rhs_created specifies whether the rhs dmat was created or not
!> \param CS_rhs_full the rhs CAUCHY-SCHWARZ (CS) screening tensor
!> \param CS_lhs_full the lhs CAUCHY-SCHWARZ (CS) screening tensor
!> \param PS_rhs_full the rhs Primitive CAUCHY-SCHWARZ (PS) screening tensor
!> \param PS_lhs_full the lhs Primitive CAUCHY-SCHWARZ (PS) screening tensor
!> \param rhsCS_created specifies whether the rhs CS tensor was created or not
!> \param lhsCS_created specifies whether the lhs CS tensor was created or not
!> \param rhsPS_created specifies whether the rhs PS tensor was created or not
!> \param lhsPS_created specifies whether the lhs PS tensor was created or not
!> \param LUPRI logical unit number of the default output file
!> \param LUERR logical unit number of the default error file
SUBROUTINE ls_create_lstensor_full(setting,lstype,AO1,AO2,AO3,AO4,Oper,Spec,intType,&
     &                             result_tensor,dmat_lhs,dmat_rhs,lhs_created,rhs_created,&
     &                             CS_rhs_full,CS_lhs_full,rhsCS_created,lhsCS_created,&
     &                             PermuteResultTensor,doscreen,lupri,luerr,exchangeOrder,&
     &                             UseMPI)
implicit none
Type(LSSETTING),intent(inout) :: SETTING
Character*(*)                 :: lstype
integer                       :: AO1,AO2,AO3,AO4,Oper,intType,Spec
TYPE(LSTENSOR),pointer        :: result_tensor,dmat_lhs,dmat_rhs
TYPE(LSTENSOR),pointer        :: CS_rhs_full,CS_lhs_full
INTEGER,intent(IN)            :: lupri,luerr
Logical,intent(OUT)           :: lhs_created,rhs_created
Logical,intent(OUT)           :: rhsCS_created,lhsCS_created
Logical,intent(OUT)           :: PermuteResultTensor,doscreen
Logical,intent(IN)            :: ExchangeOrder,UseMPI
!
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4),AObuild2(4)
Integer              :: nAObuilds,nAObuilds2,nAtoms1,nAtoms2,nAtoms3,nAtoms4
integer              :: ndim2(5),nmat
Logical              :: SameAllFRAG,useAO(4),fiveDim,ODscreen,batchOnly,FourCenterERI,ODscreen13
Logical              :: MBIE_SCREENINT,LHSGAB,LinK
type(matrix)         :: tmp
integer              :: AO_order(4)

ODscreen13 = .FALSE.
LinK = lstype.EQ.'AC_TYPE'

NULLIFY(result_tensor)
NULLIFY(dmat_lhs)
NULLIFY(dmat_rhs)

ndim2 = setting%output%ndim
nmat = setting%output%ndim(5)
CALL init_integral_input(INT_INPUT,SETTING)

Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN

INT_INPUT%operator = Oper

IF(INT_INPUT%DO_PROP)CALL SET_PROPINFO1(Oper,INT_INPUT)

! Set up AOs for the full system
CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

!call lstimer('cr-ini',ts,te,6)
! Set up the input densities (if any)
ODscreen = .FALSE.
IF (lstype.EQ.'AB_TYPE')ODscreen = .TRUE.
CALL ls_setIntegralInputDensities(INT_INPUT,setting,lupri,.TRUE.,.TRUE.,LinK,ODscreen)

lhs_created = INT_INPUT%uselst_DLHS
rhs_created = INT_INPUT%uselst_DRHS
dmat_lhs => INT_INPUT%lst_DLHS
dmat_rhs => INT_INPUT%lst_DRHS

!call lstimer('cr-id',ts,te,6)
#ifdef VAR_MPI
IF(UseMPI)THEN
   !Permutational Symmetry
   !In order to exploit permutational symmetry we set the lower triangular density matrix to 0
   !and multiply the upper triangular matrix by 2. This must be done in combination with the screening matrix
   !and a extraction of lower triangular from the upper triangluar result matrix
   !In the case of sameODs the two matrices will point to the same memory location, so we only do it once
   !
   ! FixMe Currently turned off for four-center eri integrals. Needs proper implementation in lstensor_full_symMat_from_triangularMat1
   !  
   FourCenterERI = (AO1.EQ.AORdefault).AND.(AO2.EQ.AORdefault).AND.(AO3.EQ.AORdefault).AND.(AO4.EQ.AORdefault).AND..NOT.rhs_created
   !
   IF (FourCenterERI.OR.INT_INPUT%DO_PROP.OR.Spec.EQ.MagDerivSpec&
        & .OR.Spec.EQ.MagDerivLSpec.OR.Spec.EQ.MagDerivRSpec.OR.Spec.EQ.MagGradSpec) THEN
      PermuteResultTensor = .FALSE.
   ELSE
      PermuteResultTensor = INT_INPUT%sameLHSaos
      IF(rhs_created.AND.INT_INPUT%sameRHSaos)THEN
         call lstensor_force_symMat_to_triangularMat(INT_INPUT%lst_DRHS)
      ENDIF
      IF(lhs_created.AND.INT_INPUT%sameLHSaos)THEN
         call lstensor_force_symMat_to_triangularMat(INT_INPUT%lst_DLHS)
      ENDIF
   ENDIF
ELSE
   PermuteResultTensor = .FALSE.
ENDIF
!call lstimer('cr-symm',ts,te,6)
#endif

! Specifies if integrals are limited to one ore more AO-batches
batchOnly = (setting%batchindex(1).NE. 0).OR.(setting%batchindex(2).NE. 0) &
     &      .OR.(setting%batchindex(3).NE. 0).OR.(setting%batchindex(4).NE. 0)

! Screen only when calculationg three- and four-center integrals
! Currently turned off when calculating integrals only for a specific batch
doscreen = ((Oper.EQ.CoulombOperator).OR.(Oper.EQ.NucpotOperator)) &
     & .AND.(SETTING%SCHEME%CS_SCREEN.OR.SETTING%SCHEME%PS_SCREEN)
doscreen = doscreen.OR.(((Oper.EQ.GGemCouOperator).OR.(Oper.EQ.GGemOperator).OR.(Oper.EQ.GGemGrdOperator).OR.&
     &(Oper.EQ.ErfcOperator).OR.(Oper.EQ.CAMOperator).OR.(Oper.EQ.GGemQuaOperator)).AND.SETTING%SCHEME%MBIE_SCREEN)

IF(INT_INPUT%DO_PROP)doscreen=.FALSE.
IF(SETTING%SCHEME%CS_INT.OR.SETTING%SCHEME%PS_INT)doscreen=.FALSE.

IF (doscreen) THEN
   Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
   Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
   Int_input%MBIE_SCREEN = .FALSE.
   IF(Oper .EQ. CoulombOperator.OR.Oper .EQ. ErfcOperator.OR.Oper.EQ.CAMOperator)&
        & Int_input%MBIE_SCREEN = SETTING%SCHEME%MBIE_SCREEN
!  Set up the screening tensors (if any)
!call lstimer('cr-scr0',ts,te,6)
   CALL AttachScreeningMatricesToInput(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
!call lstimer('cr-scr1',ts,te,6)
   setting%output%ndim = ndim2
   IF(.NOT.associated(INT_INPUT%LST_GAB_RHS))THEN
      CALL LSQUIT('RHS GAB ERROR in ls_create_lstensor_full',-1)
   ENDIF
!  If pointers are used in AttachScreeningMatricesToInput make a copy here before symmetrizing below
   IF (INT_Input%GAB_RHSusePointer) THEN
     NULLIFY(CS_rhs_full)
     ALLOCATE(CS_rhs_full)
     call copy_lstensor_to_lstensor(INT_INPUT%LST_GAB_RHS,CS_rhs_full)
     rhsCS_created               = .TRUE.
     INT_Input%GAB_RHSusePointer = .FALSE.
   ELSE
     CS_rhs_full => INT_INPUT%LST_GAB_RHS
     rhsCS_created = .TRUE.
   ENDIF
   IF (INT_Input%GAB_LHSusePointer) THEN
     NULLIFY(CS_lhs_full)
     ALLOCATE(CS_lhs_full)
     call copy_lstensor_to_lstensor(INT_INPUT%LST_GAB_LHS,CS_lhs_full)
     lhsCS_created               = .TRUE.
     INT_Input%GAB_LHSusePointer = .FALSE.
   ELSE
     CS_lhs_full => INT_INPUT%LST_GAB_LHS
     lhsCS_created = .TRUE.
   ENDIF
#ifdef VAR_MPI
   !Permutational Symmetry
   !In order to exploit permutational symmetry we set the lower triangular screening matrices
   !to zero. In the case of sameODs the two matrices will point to the same memory location, 
   !so we only do it once
   IF (PermuteResultTensor) THEN
     IF(INT_INPUT%sameRHSaos)THEN
        call lstensor_zero_lowertriangular(CS_rhs_full)
     ENDIF
     IF((INT_INPUT%sameLHSaos).AND..NOT.INT_INPUT%sameODs)THEN
        call lstensor_zero_lowertriangular(CS_lhs_full)
     ENDIF
   ENDIF
!  call lstimer('cr-scrm',ts,te,6)
#endif
ELSE
   nullify(CS_lhs_full)
   nullify(CS_rhs_full)
   lhsCS_created = .FALSE.
   rhsCS_created = .FALSE.
ENDIF

allocate(result_tensor)
call LSTENSOR_nullify(result_tensor)

fiveDim  = .FALSE.
IF (lstype.EQ.'AB_TYPE') THEN
  useAO(1) = .TRUE.
  useAO(2) = .TRUE.
  useAO(3) = .FALSE.
  useAO(4) = .FALSE.
ELSE IF (lstype.EQ.'CD_TYPE') THEN
  useAO(1) = .FALSE.
  useAO(2) = .FALSE.
  useAO(3) = .TRUE.
  useAO(4) = .TRUE.
ELSE IF (lstype.EQ.'AC_TYPE') THEN
  useAO(1) = .TRUE.
  useAO(2) = .TRUE.
  useAO(3) = .FALSE.
  useAO(4) = .FALSE.
ELSE IF (lstype.EQ.'FULLINT') THEN
  useAO(1) = .TRUE.
  useAO(2) = .TRUE.
  useAO(3) = .TRUE.
  useAO(4) = .TRUE.
ELSE
  write(lupri,*) 'Error in ls_create_lstensor_full, no corresponding case for lstype =',lstype
  CALL LSQUIT('Error in ls_create_lstensor_full, no corresponding case for lstype',-1)
ENDIF
ODscreen = SETTING%SCHEME%OD_SCREEN

!call lstimer('cr-B',ts,te,6)
! Set up result tensor based on Spec and AOs
SELECT CASE(Spec)
CASE(RegularSpec)
   fiveDim = .TRUE.
  !Default undifferentiated case

!WARNING THIS SHOULD BE SET BETTER SOMEWHERE ELSE
   !<a|c>
   IF (ODscreen) THEN
     IF( ((AO1.EQ.AORdefault).AND.(AO2.EQ.AOEmpty)).AND.((AO3.EQ.AORdefault).AND.(AO4.EQ.AOEmpty)))THEN
        IF((Oper.EQ.OverlapOperator).OR.(Oper.EQ.KineticOperator))THEN
           ODscreen13 = .TRUE.
        ENDIF
     ENDIF
     !<alpha|cd>
     IF( ((AO1.EQ.AODFdefault).AND.(AO2.EQ.AOEmpty)).AND.((AO3.EQ.AORdefault).AND.(AO4.EQ.AORdefault)))THEN
        IF(Oper.EQ.OverlapOperator)THEN
           ODscreen13 = .TRUE.
        ENDIF
     ENDIF
   ELSE
     ODscreen13 = .FALSE.
   ENDIF
CASE(pChargeSPec)
   call init_pcharge_lstensor(result_tensor,ndim2(1),nmat,lupri)
  !Default undifferentiated case
   PermuteResultTensor = .FALSE.
CASE(GradientSpec)
   nAtoms1 =  setting%fragment(1)%p%nAtoms
   nAtoms2 =  setting%fragment(2)%p%nAtoms
   nAtoms3 =  setting%fragment(3)%p%nAtoms
   nAtoms4 =  setting%fragment(4)%p%nAtoms
   call SET_SAMEALLFRAG(sameAllFrag,setting%sameFrag,setting%nAO)
   call init_gradientlstensor(result_tensor,nAtoms1,nAtoms2,nAtoms3,nAtoms4,sameAllFrag,&
        & nmat,3,.FALSE.,lupri)
  PermuteResultTensor = .FALSE.
CASE(GeoDerivLHSSpec)
   ODscreen = .FALSE.
   fiveDim = .TRUE.
   IF (Oper.EQ.KineticOperator) INT_INPUT%sameODs = .FALSE.
CASE(GeoDerivRHSSpec)
   ODscreen = .FALSE.
   fiveDim = .TRUE.
   IF (Oper.EQ.KineticOperator) INT_INPUT%sameODs = .FALSE.
CASE(GeoDerivSpec)
   ODscreen = .FALSE.
   fiveDim = .TRUE.
   IF (Oper.EQ.KineticOperator) INT_INPUT%sameODs = .FALSE.
CASE(GeoDerivCoulombSpec)
   ODscreen = .FALSE.
   fiveDim = .TRUE.
   INT_INPUT%sameODs = .FALSE.
   useAO(1) = .TRUE.
   useAO(2) = .TRUE.
   useAO(3) = .FALSE.
   useAO(4) = .FALSE.
CASE(MagGradSpec)
   call init_gradientlstensor(result_tensor,1,1,1,1,.TRUE.,nmat,1,.TRUE.,lupri)
   IF(Oper.EQ.OverlapOperator)THEN
      INT_INPUT%sameLHSAOs  = .FALSE.
   ENDIF
   INT_INPUT%sameODs  = .FALSE.
CASE(MagDerivSpec)
   fiveDim = .TRUE.
   INT_INPUT%sameODs  = .FALSE.
CASE(MagDerivLSpec)
   fiveDim = .TRUE.
   INT_INPUT%sameODs  = .FALSE.
CASE(MagDerivRSpec)
   fiveDim = .TRUE.
   INT_INPUT%sameLHSAOs  = .FALSE.
   INT_INPUT%sameODs  = .FALSE.
CASE (EcontribSpec)
   call init_lstensor_1dim(result_tensor,ndim2(5),lupri)
   PermuteResultTensor = .FALSE.
CASE (magderivEcontribSpec)
   call init_lstensor_1dim(result_tensor,ndim2(5),lupri)
   PermuteResultTensor = .FALSE.
CASE DEFAULT
  CALL LSQUIT('Error in ls_create_lstensor_full. Wrong Spec case',lupri)
END SELECT

LHSGAB=.TRUE.
IF(.NOT.ExchangeOrder)THEN
   AO_order(1) = 1
   AO_order(2) = 2
   AO_order(3) = 3
   AO_order(4) = 4
ELSE
   AO_order(1) = 1
   AO_order(2) = 3
   AO_order(3) = 2
   AO_order(4) = 4
   ODscreen=.FALSE.
   ODscreen13=.FALSE.
ENDIF

MBIE_SCREENINT = Setting%scheme%MBIE_SCREEN.AND.(Oper.EQ.CoulombOperator.OR.Oper.EQ.ErfcOperator.OR.Oper.EQ.CAMOperator)
IF((AO1.EQ.AOEmpty).OR.(AO2.EQ.AOEmpty))MBIE_SCREENINT=.false.
IF (setting%scheme%cs_int.OR.setting%scheme%ps_int)THEN
   CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
   IF(setting%scheme%cs_int)THEN   
      CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,Contractedinttype,AObuild2,nAObuilds2,SETTING,LUPRI,LUERR,LHSGAB)
      call init_cs_lstensor(result_tensor,Int_Input%AO(1)%p,&
           &Int_Input%AO(2)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),INT_INPUT%OD_SCREEN,lupri)
      CALL FreeInputAO(AObuild2,nAObuilds2,LUPRI)
   ENDIF
   IF(setting%scheme%ps_int)THEN
      CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,Primitiveinttype,AObuild2,nAObuilds2,SETTING,LUPRI,LUERR,LHSGAB)
      call init_ps_lstensor(result_tensor,Int_Input%AO(1)%p,&
           &Int_Input%AO(2)%p,INT_INPUT%AOdim(1),INT_INPUT%AOdim(2),INT_INPUT%OD_SCREEN,lupri)
      call set_lst_maxprimgab(result_tensor)
      CALL FreeInputAO(AObuild2,nAObuilds2,LUPRI)
   ENDIF
   IF(MBIE_SCREENINT)THEN
      CALL SetInputAO(INT_INPUT,AO1,AO2,AOEmpty,AOEmpty,Contractedinttype,AObuild2,nAObuilds2,SETTING,LUPRI,LUERR,LHSGAB)
      call init_MBIE_lstensor_5dim(result_tensor,Int_Input%AO(1)%p,Int_Input%AO(2)%p,.FALSE.,lupri)
      CALL FreeInputAO(AObuild2,nAObuilds2,LUPRI)
   ENDIF
   CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,&
        & SETTING,LUPRI,LUERR,.TRUE.)
ELSE
   IF (fiveDim) THEN
      call init_lstensor_5dim(result_tensor,Int_Input%AO(AO_order(1))%p,Int_Input%AO(AO_order(2))%p,&
           & Int_Input%AO(AO_order(3))%p,Int_Input%AO(AO_order(4))%p,ndim2(1),ndim2(2),ndim2(3),&
           & ndim2(4),ndim2(5),useAO(1),useAO(2),useAO(3),useAO(4),ODscreen,ODscreen13,lupri)
      IF(useAO(1))THEN
         IF(Int_Input%AO(AO_order(1))%p%nbast.NE.ndim2(1))call lsquit('dim mismatch lstensor',-1)
      ENDIF
      IF(useAO(2))THEN
         IF(Int_Input%AO(AO_order(2))%p%nbast.NE.ndim2(2))call lsquit('dim mismatch lstensor',-1)
      ENDIF
      IF(useAO(3))THEN
         IF(Int_Input%AO(AO_order(3))%p%nbast.NE.ndim2(3))call lsquit('dim mismatch lstensor',-1)
      ENDIF
      IF(useAO(4))THEN
         IF(Int_Input%AO(AO_order(4))%p%nbast.NE.ndim2(4))call lsquit('dim mismatch lstensor',-1)
      ENDIF
   ENDIF
ENDIF
!call lstimer('cr-5-dim',ts,te,6)

#ifdef VAR_MPI
IF(UseMPI)THEN
   CALL set_reduced_screen_info(setting%redCS,INT_INPUT%AO,dmat_lhs,dmat_rhs,lhs_created,rhs_created,&
        &                      CS_rhs_full,CS_lhs_full,INT_INPUT%CS_THRLOG,doscreen,lupri)
   !call lstimer('cr-red',ts,te,6)
   call print_reduced_screening_info(setting%redCS,setting%scheme%intprint,lupri)
   !call lstimer('cr-pred',ts,te,6)
ENDIF
#endif

! Do not free the input densities, this should be done in ls_free_lstensor_full
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
setting%output%ndim = ndim2
END SUBROUTINE ls_create_lstensor_full

subroutine ls_memdist_lstensor_SetupFullinfo(setting,lstype,AO1,AO2,AO3,AO4,&
     & Oper,Spec,intType,&
     & result_tensor,result_tensor_other,dmat_lhs,dmat_rhs,lhs_created,rhs_created,&
     & CS_rhs_full,CS_lhs_full,rhsCS_created,lhsCS_created,&
     & PermuteResultTensor,doscreen,ForceRHSsymDMAT,ForceLHSsymDMAT,&
     & LHSpartioning,RHSpartioning,BOTHpartioning,lupri,luerr)
implicit none
Type(LSSETTING),intent(inout) :: SETTING
Character*(*)                 :: lstype
integer                       :: AO1,AO2,AO3,AO4,Oper,intType,Spec
TYPE(LSTENSOR),pointer        :: result_tensor,result_tensor_other,dmat_lhs,dmat_rhs
TYPE(LSTENSOR),pointer        :: CS_rhs_full,CS_lhs_full
INTEGER,intent(IN)            :: lupri,luerr
Logical,intent(INOUT)         :: lhs_created,rhs_created,ForceRHSsymDMAT,ForceLHSsymDMAT
Logical,intent(INOUT)         :: rhsCS_created,lhsCS_created
Logical,intent(INOUT)         :: PermuteResultTensor,doscreen
Logical,intent(IN)            :: LHSpartioning,RHSpartioning,BOTHpartioning
!
#ifdef VAR_MPI
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4),AObuild2(4)
Integer              :: nAObuilds,nAObuilds2,nAtoms1,nAtoms2,nAtoms3,nAtoms4
integer              :: ndim2(5),nmat
Logical              :: SameAllFRAG,useAO(4),fiveDim,ODscreen,batchOnly,FourCenterERI,ODscreen13
Logical              :: MBIE_SCREENINT,LHSGAB
type(matrix)         :: tmp
integer              :: AO_order(4)
INTEGER :: LHS1,LHS2,RHS1,RHS2,nbastLHS1,nbastLHS2,nbastRHS1,nbastRHS2
INTEGER :: nDmatLHS,nDmatRHS,idmat
LOGICAL :: useAO1,useAO2,LinK
ODscreen13 = .FALSE.
LinK = lstype.EQ.'AC_TYPE'


NULLIFY(result_tensor)
NULLIFY(result_tensor_other)

ndim2 = setting%output%ndim
nmat = setting%output%ndim(5)
CALL init_integral_input(INT_INPUT,SETTING)
Int_input%OD_SCREEN = SETTING%SCHEME%OD_SCREEN
Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
INT_INPUT%operator = Oper
IF(INT_INPUT%DO_PROP)CALL SET_PROPINFO1(Oper,INT_INPUT)

! Set up AOs for the full system
CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)

!call lstimer('cr-ini',ts,te,6)
!=====================================================
!   The Density Matrices    
!=====================================================
! LHS Density Matrix
IF (setting%LHSdfull.OR.setting%LHSdmat)THEN   
   nullify(dmat_lhs)
   IF(LHSpartioning)THEN
      !use memory distributed Density matrices as we do a LHS partitioning
      nDmatLHS = setting%nDmatLHS
      INT_INPUT%NDMAT_LHS = ndmatLHS
      LHS1 = setting%LHSdmatAOindex1
      LHS2 = setting%LHSdmatAOindex2
      nbastLHS1 = INT_INPUT%AOdim(LHS1)
      nbastLHS2 = INT_INPUT%AOdim(LHS2)
      useAO1 = .NOT.INT_INPUT%AO(LHS1)%p%empty
      useAO2 = .NOT.INT_INPUT%AO(LHS2)%p%empty
      ODscreen = .FALSE. 
      ALLOCATE(dmat_lhs)
      call LSTENSOR_nullify(dmat_lhs)
      call memdist_lstensor_SetupFullinfo(dmat_lhs,&
           &INT_INPUT%AO(LHS1)%p,INT_INPUT%AO(LHS2)%p,&
           & nbastLHS1,nbastLHS2,nDmatLHS,useAO1,useAO2,ODscreen,lupri)
      lhs_created = .TRUE.
   ELSE
      !FULL Densities 
      CALL ls_setIntegralInputDensities(INT_INPUT,setting,lupri,.TRUE.,.FALSE.,LinK,ODscreen)
      lhs_created = INT_INPUT%uselst_DLHS
      dmat_lhs => INT_INPUT%lst_DLHS
   ENDIF
ELSE
   lhs_created = .FALSE.
ENDIF

IF (setting%RHSdfull.OR.setting%RHSdmat)THEN
   nullify(dmat_rhs)
   IF(RHSpartioning)THEN
      !use memory distributed Density matrices as we do a RHS partitioning
      nDmatRHS = setting%nDmatRHS
      INT_INPUT%NDMAT_RHS = ndmatRHS
      RHS1 = setting%RHSdmatAOindex1
      RHS2 = setting%RHSdmatAOindex2
      nbastRHS1 = INT_INPUT%AOdim(RHS1)
      nbastRHS2 = INT_INPUT%AOdim(RHS2)
      useAO1 = .NOT.INT_INPUT%AO(RHS1)%p%empty
      useAO2 = .NOT.INT_INPUT%AO(RHS2)%p%empty
      ALLOCATE(dmat_rhs)
      call LSTENSOR_nullify(dmat_rhs)
      call memdist_lstensor_SetupFullinfo(dmat_rhs,&
           &INT_INPUT%AO(RHS1)%p,INT_INPUT%AO(RHS2)%p,&
           & nbastRHS1,nbastRHS2,nDmatRHS,useAO1,useAO2,ODscreen,lupri)
      rhs_created = .TRUE.
   ELSE
      !FULL Densities 
      CALL ls_setIntegralInputDensities(INT_INPUT,setting,lupri,.FALSE.,.TRUE.,LinK,ODscreen)
      rhs_created = INT_INPUT%uselst_DRHS
      dmat_rhs => INT_INPUT%lst_DRHS
   ENDIF
ELSE
   rhs_created = .FALSE.
ENDIF

!call lstimer('cr-id',ts,te,6)

!Permutational Symmetry
!In order to exploit permutational symmetry we set the lower triangular density matrix to 0
!and multiply the upper triangular matrix by 2. This must be done in combination with the screening matrix
!and a extraction of lower triangular from the upper triangluar result matrix
!In the case of sameODs the two matrices will point to the same memory location, so we only do it once
!
! FixMe Currently turned off for four-center eri integrals. Needs proper implementation in lstensor_full_symMat_from_triangularMat1
!  
FourCenterERI = (AO1.EQ.AORdefault).AND.(AO2.EQ.AORdefault).AND.(AO3.EQ.AORdefault).AND.(AO4.EQ.AORdefault).AND..NOT.rhs_created
!
!WARNING THE DENSITIES IS NOT YET SETUP UP - ONLY full molecular info
! this must be done later - and in a MPI distributed fashion. 
ForceRHSsymDMAT = .FALSE.
ForceLHSsymDMAT = .FALSE.
IF (FourCenterERI.OR.INT_INPUT%DO_PROP) THEN
  PermuteResultTensor = .FALSE.
ELSE
  PermuteResultTensor = INT_INPUT%sameLHSaos
  IF(rhs_created.AND.INT_INPUT%sameRHSaos)THEN
     IF(RHSpartioning)THEN
        !done with matrix operations outside
        ForceRHSsymDMAT = .TRUE.
     ELSE     
        call lstensor_force_symMat_to_triangularMat(INT_INPUT%lst_DRHS)
     ENDIF
  ENDIF
  IF(lhs_created.AND.INT_INPUT%sameLHSaos)THEN
     IF(LHSpartioning)THEN
        ForceLHSsymDMAT = .TRUE.
     ELSE
        call lstensor_force_symMat_to_triangularMat(INT_INPUT%lst_DLHS)
     ENDIF
  ENDIF
ENDIF
!call lstimer('cr-symm',ts,te,6)

! Specifies if integrals are limited to one ore more AO-batches
batchOnly = (setting%batchindex(1).NE. 0).OR.(setting%batchindex(2).NE. 0) &
     &      .OR.(setting%batchindex(3).NE. 0).OR.(setting%batchindex(4).NE. 0)

! Screen only when calculationg three- and four-center integrals
! Currently turned off when calculating integrals only for a specific batch
doscreen = ((Oper.EQ.CoulombOperator).OR.(Oper.EQ.NucpotOperator)) &
     & .AND.(SETTING%SCHEME%CS_SCREEN.OR.SETTING%SCHEME%PS_SCREEN)
doscreen = doscreen.OR.(((Oper.EQ.GGemCouOperator).OR.(Oper.EQ.GGemOperator).OR.(Oper.EQ.GGemGrdOperator)&
     &.OR.(Oper.EQ.ErfcOperator).OR.(Oper.EQ.GGemQuaOperator).OR.(Oper.EQ.CAMOperator))&
     & .AND.SETTING%SCHEME%MBIE_SCREEN)

IF(INT_INPUT%DO_PROP)doscreen=.FALSE.
IF(SETTING%SCHEME%CS_INT.OR.SETTING%SCHEME%PS_INT)doscreen=.FALSE.

IF (doscreen) THEN
   Int_input%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
   Int_input%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
   Int_input%MBIE_SCREEN = .FALSE.
   IF(Oper .EQ. CoulombOperator.OR.Oper .EQ. ErfcOperator.OR.Oper.EQ.CAMOperator)&
        & Int_input%MBIE_SCREEN = SETTING%SCHEME%MBIE_SCREEN
!  Set up the screening tensors (if any)
   CALL AttachScreeningMatricesToInput(INT_INPUT,AO1,AO2,AO3,AO4,Oper,SETTING,LUPRI,LUERR)
   setting%output%ndim = ndim2
   IF(.NOT.associated(INT_INPUT%LST_GAB_RHS))THEN
      CALL LSQUIT('RHS GAB ERROR in ls_create_lstensor_full',-1)
   ENDIF
!  If pointers are used in AttachScreeningMatricesToInput make a copy here before symmetrizing below
   IF (INT_Input%GAB_RHSusePointer) THEN
     NULLIFY(CS_rhs_full)
     ALLOCATE(CS_rhs_full)
     call LSTENSOR_nullify(CS_rhs_full)     
     call copy_lstensor_to_lstensor(INT_INPUT%LST_GAB_RHS,CS_rhs_full)
     rhsCS_created               = .TRUE.
     INT_Input%GAB_RHSusePointer = .FALSE.
   ELSE
     CS_rhs_full => INT_INPUT%LST_GAB_RHS
     rhsCS_created = .TRUE.
   ENDIF
   IF (INT_Input%GAB_LHSusePointer) THEN
     NULLIFY(CS_lhs_full)
     ALLOCATE(CS_lhs_full)
     call LSTENSOR_nullify(CS_lhs_full)     
     call copy_lstensor_to_lstensor(INT_INPUT%LST_GAB_LHS,CS_lhs_full)
     lhsCS_created               = .TRUE.
     INT_Input%GAB_LHSusePointer = .FALSE.
   ELSE
     CS_lhs_full => INT_INPUT%LST_GAB_LHS
     lhsCS_created = .TRUE.
   ENDIF
   !Permutational Symmetry
   !In order to exploit permutational symmetry we set the lower triangular screening matrices
   !to zero. In the case of sameODs the two matrices will point to the same memory location, 
   !so we only do it once
   IF (PermuteResultTensor) THEN
     IF(INT_INPUT%sameRHSaos)THEN
        call lstensor_zero_lowertriangular(CS_rhs_full)
     ENDIF
     IF((INT_INPUT%sameLHSaos).AND..NOT.INT_INPUT%sameODs)THEN
        call lstensor_zero_lowertriangular(CS_lhs_full)
     ENDIF
   ENDIF
!  call lstimer('cr-scrm',ts,te,6)
ELSE
   nullify(CS_lhs_full)
   nullify(CS_rhs_full)
   lhsCS_created = .FALSE.
   rhsCS_created = .FALSE.
ENDIF

allocate(result_tensor)
call LSTENSOR_nullify(result_tensor)
IF(BOTHpartioning)THEN
   allocate(result_tensor_other)
   call LSTENSOR_nullify(result_tensor_other)
ENDIF
fiveDim  = .FALSE.
IF (lstype.EQ.'AB_TYPE') THEN
  useAO(1) = .TRUE.
  useAO(2) = .TRUE.
  useAO(3) = .FALSE.
  useAO(4) = .FALSE.
ELSE IF (lstype.EQ.'CD_TYPE') THEN
  useAO(1) = .FALSE.
  useAO(2) = .FALSE.
  useAO(3) = .TRUE.
  useAO(4) = .TRUE.
ELSE IF (lstype.EQ.'AC_TYPE') THEN
  useAO(1) = .TRUE.
  useAO(2) = .TRUE.
  useAO(3) = .FALSE.
  useAO(4) = .FALSE.
ELSE IF (lstype.EQ.'FULLINT') THEN
  useAO(1) = .TRUE.
  useAO(2) = .TRUE.
  useAO(3) = .TRUE.
  useAO(4) = .TRUE.
ELSE
  write(lupri,*) 'Error in ls_create_lstensor_full, no corresponding case for lstype =',lstype
  CALL LSQUIT('Error in ls_create_lstensor_full, no corresponding case for lstype',-1)
ENDIF
ODscreen = SETTING%SCHEME%OD_SCREEN

!call lstimer('cr-B',ts,te,6)
! Set up result tensor based on Spec and AOs
SELECT CASE(Spec)
CASE(RegularSpec)
   fiveDim = .TRUE.
  !Default undifferentiated case

!WARNING THIS SHOULD BE SET BETTER SOMEWHERE ELSE
   !<a|c>
   IF( ((AO1.EQ.AORdefault).AND.(AO2.EQ.AOEmpty)).AND.((AO3.EQ.AORdefault).AND.(AO4.EQ.AOEmpty)))THEN
      IF((Oper.EQ.OverlapOperator).OR.(Oper.EQ.KineticOperator))THEN
         ODscreen13 = .TRUE.
      ENDIF
   ENDIF
   !<alpha|cd>
   IF( ((AO1.EQ.AODFdefault).AND.(AO2.EQ.AOEmpty)).AND.((AO3.EQ.AORdefault).AND.(AO4.EQ.AORdefault)))THEN
      IF(Oper.EQ.OverlapOperator)THEN
         ODscreen13 = .TRUE.
      ENDIF
   ENDIF
CASE(pChargeSPec)
   call init_pcharge_lstensor(result_tensor,ndim2(1),nmat,lupri)
  !Default undifferentiated case
   PermuteResultTensor = .FALSE.
CASE(GradientSpec)
   nAtoms1 =  setting%fragment(1)%p%nAtoms
   nAtoms2 =  setting%fragment(2)%p%nAtoms
   nAtoms3 =  setting%fragment(3)%p%nAtoms
   nAtoms4 =  setting%fragment(4)%p%nAtoms
   call SET_SAMEALLFRAG(sameAllFrag,setting%sameFrag,setting%nAO)
   call init_gradientlstensor(result_tensor,nAtoms1,nAtoms2,nAtoms3,nAtoms4,sameAllFrag,&
        & nmat,3,.FALSE.,lupri)
  PermuteResultTensor = .FALSE.
CASE(GeoDerivLHSSpec)
   ODscreen = .FALSE.
   fiveDim = .TRUE.
   IF (Oper.EQ.KineticOperator) INT_INPUT%sameODs = .FALSE.
CASE(MagGradSpec)
   call init_gradientlstensor(result_tensor,1,1,1,1,.TRUE.,nmat,1,.TRUE.,lupri)
   IF(Oper.EQ.OverlapOperator)THEN
      INT_INPUT%sameLHSAOs  = .FALSE.
   ENDIF
   INT_INPUT%sameODs  = .FALSE.
CASE(MagDerivSpec)
   fiveDim = .TRUE.
   INT_INPUT%sameODs  = .FALSE.
CASE(MagDerivLSpec)
   fiveDim = .TRUE.
   INT_INPUT%sameODs  = .FALSE.
CASE(MagDerivRSpec)
   fiveDim = .TRUE.
   INT_INPUT%sameLHSAOs  = .FALSE.
   INT_INPUT%sameODs  = .FALSE.
CASE (EcontribSpec)
   call init_lstensor_1dim(result_tensor,ndim2(5),lupri)
   PermuteResultTensor = .FALSE.
CASE (magderivEcontribSpec)
   call init_lstensor_1dim(result_tensor,ndim2(5),lupri)
   PermuteResultTensor = .FALSE.
CASE DEFAULT
  CALL LSQUIT('ls_memdist_lstensor_SetupFullinfo',lupri)
END SELECT

AO_order(1) = 1
AO_order(2) = 2
AO_order(3) = 3
AO_order(4) = 4
IF (fiveDim) THEN      
   IF(LHSpartioning)THEN
      IF(useAO(3))CALL LSQUIT('useAO3: memdist LHS result require 2dim resulttensor',-1)
      IF(useAO(4))CALL LSQUIT('useAO4: memdist LHS result require 2dim resulttensor',-1)
      !use memory distributed matrix for result as we do a LHS partitioning
      call memdist_lstensor_SetupFullinfo(result_tensor,&
           & Int_Input%AO(AO_order(1))%p,Int_Input%AO(AO_order(2))%p,&
           & ndim2(1),ndim2(2),ndim2(5),useAO(1),useAO(2),ODscreen,lupri)
      IF(BOTHpartioning)THEN
         call memdist_lstensor_SetupFullinfo(result_tensor_other,&
              & Int_Input%AO(AO_order(1))%p,Int_Input%AO(AO_order(2))%p,&
              & ndim2(1),ndim2(2),ndim2(5),useAO(1),useAO(2),ODscreen,lupri)
      ENDIF
   ELSE
      call init_lstensor_5dim(result_tensor,Int_Input%AO(AO_order(1))%p,Int_Input%AO(AO_order(2))%p,&
           & Int_Input%AO(AO_order(3))%p,Int_Input%AO(AO_order(4))%p,ndim2(1),ndim2(2),ndim2(3),&
           & ndim2(4),ndim2(5),useAO(1),useAO(2),useAO(3),useAO(4),ODscreen,ODscreen13,lupri)
   ENDIF
ENDIF
CALL set_reduced_screen_info(setting%redCS,INT_INPUT%AO,dmat_lhs,dmat_rhs,&
     & .FALSE.,.FALSE.,CS_rhs_full,CS_lhs_full,INT_INPUT%CS_THRLOG,doscreen,lupri)
call print_reduced_screening_info(setting%redCS,setting%scheme%intprint,lupri)

! Do not free the input densities, this should be done in ls_free_lstensor_full
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
setting%output%ndim = ndim2
#endif
end subroutine ls_memdist_lstensor_SetupFullinfo

subroutine ls_memdist_lstensor_SetupLocalinfoDLHS(setting,AO1,AO2,AO3,AO4,&
     & intType,dmat_lhs,lhs_created,nAtom1,nAtom2,atoms1,&
     & atoms2,LHSpartioning,lupri,luerr)
implicit none
Type(LSSETTING),intent(inout) :: SETTING
integer                       :: AO1,AO2,AO3,AO4,intType
logical                       :: lhs_created
TYPE(LSTENSOR),pointer        :: dmat_lhs
INTEGER,intent(IN)            :: lupri,luerr
integer,intent(in) :: natom1,natom2
integer,intent(in) :: atoms1(natom1),atoms2(natom2)
logical :: LHSpartioning
!
#ifdef VAR_MPI
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
INTEGER :: LHS1,LHS2,nbastLHS1,nbastLHS2,nDmatLHS,idmat
LOGICAL :: useAO1,useAO2

CALL init_integral_input(INT_INPUT,SETTING)
! Set up AOs for the SUB system! Each node does this for itself once. 
CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
!=====================================================
!   The Density Matrices    
!=====================================================
IF (lhs_created.AND.LHSpartioning)THEN   
   LHS1 = setting%LHSdmatAOindex1
   LHS2 = setting%LHSdmatAOindex2
   nbastLHS1 = INT_INPUT%AOdim(LHS1)
   nbastLHS2 = INT_INPUT%AOdim(LHS2)
   useAO1 = .NOT.INT_INPUT%AO(LHS1)%p%empty
   useAO2 = .NOT.INT_INPUT%AO(LHS2)%p%empty
   call memdist_lstensor_SetupLocalinfo(dmat_lhs,&
        &INT_INPUT%AO(LHS1)%p,INT_INPUT%AO(LHS2)%p,&
        & nbastLHS1,nbastLHS2,useAO1,useAO2,&
        & nAtom1,nAtom2,atoms1,atoms2,lupri)
ENDIF
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
#endif
end subroutine ls_memdist_lstensor_SetupLocalinfoDLHS

subroutine ls_memdist_lstensor_SetupLocalinfoDRHS(setting,AO1,AO2,AO3,AO4,&
     & intType,dmat_rhs,rhs_created,nAtom1,nAtom2,atoms1,atoms2,RHSpartioning,lupri,luerr)
implicit none
Type(LSSETTING),intent(inout) :: SETTING
integer,intent(in)            :: AO1,AO2,AO3,AO4,intType
logical,intent(in)            :: rhs_created
TYPE(LSTENSOR),pointer        :: dmat_rhs
INTEGER,intent(IN)            :: lupri,luerr
integer,intent(in) :: natom1,natom2,atoms1(natom1),atoms2(natom2)
logical :: RHSpartioning
!
#ifdef VAR_MPI
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
INTEGER :: RHS1,RHS2,nbastRHS1,nbastRHS2,nDmatRHS,idmat
LOGICAL :: useAO1,useAO2

CALL init_integral_input(INT_INPUT,SETTING)
! Set up AOs for the SUB system! Each node does this for itself once. 
CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
!=====================================================
!   The Density Matrices    
!=====================================================
IF (rhs_created.AND.RHSpartioning)THEN   
   RHS1 = setting%RHSdmatAOindex1
   RHS2 = setting%RHSdmatAOindex2
   nbastRHS1 = INT_INPUT%AOdim(RHS1)
   nbastRHS2 = INT_INPUT%AOdim(RHS2)
   useAO1 = .NOT.INT_INPUT%AO(RHS1)%p%empty
   useAO2 = .NOT.INT_INPUT%AO(RHS2)%p%empty
   call memdist_lstensor_SetupLocalinfo(dmat_rhs,&
        &INT_INPUT%AO(RHS1)%p,INT_INPUT%AO(RHS2)%p,&
        & nbastRHS1,nbastRHS2,useAO1,useAO2,&
        & nAtom1,nAtom2,atoms1,atoms2,lupri)
ENDIF
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
#endif
end subroutine ls_memdist_lstensor_SetupLocalinfoDRHS

subroutine ls_memdist_lstensor_SetupLocalinfoRESLHS(setting,AO1,AO2,AO3,AO4,&
     & intType,result_tensor,nAtom1,nAtom2,atoms1,atoms2,AOindex1,AOindex2,&
     & LHSpartioning,lupri,luerr)
implicit none
Type(LSSETTING),intent(inout) :: SETTING
integer,intent(IN)            :: AO1,AO2,AO3,AO4,intType,AOindex1,AOindex2
TYPE(LSTENSOR),pointer        :: result_tensor
INTEGER,intent(IN)            :: lupri,luerr
integer,intent(IN) :: natom1,natom2
integer,intent(IN) :: atoms1(natom1),atoms2(natom2)
logical,intent(IN) :: LHSpartioning
!
#ifdef VAR_MPI
TYPE(INTEGRALINPUT)  :: INT_INPUT
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
INTEGER :: nbast1,nbast2
LOGICAL :: useAO1,useAO2

CALL init_integral_input(INT_INPUT,SETTING)
! Set up AOs for the SUB system! Each node does this for itself once. 
CALL SetInputAO(INT_INPUT,AO1,AO2,AO3,AO4,intType,AObuild,nAObuilds,SETTING,LUPRI,LUERR,.TRUE.)
!=====================================================
!   The Density Matrices    
!=====================================================
IF (LHSpartioning)THEN   
   nbast1 = INT_INPUT%AOdim(AOindex1)
   nbast2 = INT_INPUT%AOdim(AOindex2)
   useAO1 = .NOT.INT_INPUT%AO(AOindex1)%p%empty
   useAO2 = .NOT.INT_INPUT%AO(AOindex2)%p%empty
   call memdist_lstensor_SetupLocalinfo(result_tensor,&
        &INT_INPUT%AO(AOindex1)%p,INT_INPUT%AO(AOindex2)%p,&
        & nbast1,nbast2,useAO1,useAO2,nAtom1,nAtom2,atoms1,atoms2,lupri)
ENDIF
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
#endif
end subroutine ls_memdist_lstensor_SetupLocalinfoRESLHS

subroutine ls_BuildAtomicGab(atomGab,nAtomA,nAtomB,AO1,AO2,&
     & Oper,SETTING,LUPRI,LUERR)
implicit none
integer :: nAtomA,nAtomB,AO1,AO2,Oper,LUPRI,LUERR
logical :: LHSGAB
integer(kind=short) :: atomGab(nAtomA,nAtomB)
type(lssetting) :: SETTING
!
TYPE(AOITEM),target  :: AObuild(4)
Integer              :: nAObuilds
TYPE(INTEGRALINPUT)  :: INT_INPUT
CALL init_integral_input(INT_INPUT,SETTING)
CALL SetInputAO(INT_INPUT,AO1,AO2,AO1,AO2,Contractedinttype,AObuild,nAObuilds,&
     & SETTING,LUPRI,LUERR,.TRUE.)
INT_INPUT%CS_SCREEN = SETTING%SCHEME%CS_SCREEN
INT_INPUT%PS_SCREEN = SETTING%SCHEME%PS_SCREEN
nullify(setting%lst_gab_LHS)
nullify(setting%lst_gab_RHS)
CALL AttachScreeningMatricesToInput(INT_INPUT,AO1,AO2,AO1,AO2,Oper,&
     & SETTING,LUPRI,LUERR)
CALL lstensor_BuildAtomicGab(AtomGab,nAtomA,nAtomB,INT_INPUT%LST_GAB_LHS,&
     & Int_input%AO(1)%p,Int_input%AO(2)%p,lupri)
CALL FreeInputAO(AObuild,nAObuilds,LUPRI)
call free_screening_matrices(Int_Input,SETTING,LUPRI,LUERR)
end subroutine ls_BuildAtomicGab

!> \brief Free lstensor density and result matrices
!> \author Simen Reine
!> \date 2011-10-25
!> \param result_tensor the tensor to store the results
!> \param dmat_lhs the tensor for "lhs" density-matrices
!> \param dmat_rhs the tensor for "rhs" density-matrices
!> \param lhs_created specifies whether the lhs dmat was created or not
!> \param rhs_created specifies whether the rhs dmat was created or not
SUBROUTINE ls_free_lstensors(dmat_lhs,dmat_rhs,lhs_created,rhs_created)
implicit none
TYPE(LSTENSOR),pointer     :: dmat_lhs,dmat_rhs
Logical,intent(IN)        :: lhs_created,rhs_created

IF (lhs_created) THEN
  CALL lstensor_free(dmat_lhs)
  deallocate(dmat_lhs)
  nullify(dmat_lhs)
ENDIF
IF (rhs_created) THEN
   CALL lstensor_free(dmat_rhs)
   deallocate(dmat_rhs)
   nullify(dmat_rhs)
ENDIF
END SUBROUTINE ls_free_lstensors

subroutine ls_free_screeninglstensors(CS_rhs,CS_lhs,rhsCS_created,lhsCS_created)
implicit none
TYPE(LSTENSOR),pointer   :: CS_rhs,CS_lhs
Logical,intent(IN)       :: rhsCS_created,lhsCS_created
IF (rhsCS_created) THEN
  call lstensor_free(CS_rhs)      
  deallocate(CS_rhs)
  nullify(CS_rhs)
ENDIF
IF (lhsCS_created) THEN
  call lstensor_free(CS_lhs)      
  deallocate(CS_lhs)
  nullify(CS_lhs)
ENDIF
end subroutine ls_free_screeninglstensors

!> \brief Attach lstensor density and result matrices to setting
!> \author Simen Reine
!> \date 2011-10-25
!> \param SETTING Integral evalualtion settings
!> \param result_tensor the tensor to store the results
!> \param dmat_lhs the tensor for "lhs" density-matrices
!> \param dmat_rhs the tensor for "rhs" density-matrices
!> \param lhs_created specifies whether the lhs dmat was created or not
!> \param rhs_created specifies whether the rhs dmat was created or not
SUBROUTINE ls_attach_lstensors_to_setting(setting,result_tensor,dmat_lhs,dmat_rhs,&
     &                                    lhs_created,rhs_created)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
TYPE(LSTENSOR),pointer        :: result_tensor,dmat_lhs,dmat_rhs
Logical,intent(IN)            :: lhs_created,rhs_created
!
Integer :: ndim1,ndim2,ndim3,ndim4,ndim5

! Set up proper dimensions
ndim1 = result_tensor%nbast(1)
ndim2 = result_tensor%nbast(2)
ndim3 = result_tensor%nbast(3)
ndim4 = result_tensor%nbast(4)
ndim5 = result_tensor%ndim5
call initIntegralOutputDims1(setting%Output,ndim1,ndim2,ndim3,ndim4,ndim5)

! Attach pointer
setting%lstensor_attached = .TRUE.
setting%output%resultTensor => result_tensor
setting%lst_dLHS     => dmat_lhs
setting%lst_dRHS     => dmat_rhs

END SUBROUTINE ls_attach_lstensors_to_setting

!> \brief Free lstensor density and result matrices from setting
!> \author Simen Reine
!> \date 2011-10-25
!> \param SETTING Integral evalualtion settings
!> \param lupri default print unit
SUBROUTINE ls_free_lstensors_from_setting(setting,lupri)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
INTEGER,INTENT(IN)            :: lupri

IF (.NOT.setting%lstensor_attached) &
     & CALL LSQUIT('Error in ls_free_lstensors_from_setting. Lstensor not attached',lupri)

setting%lstensor_attached = .FALSE.
setting%output%ndim = 0
!NULLIFY(setting%output%resultTensor)
NULLIFY(setting%lst_dLHS)
NULLIFY(setting%lst_dRHS)

END SUBROUTINE ls_free_lstensors_from_setting

subroutine ls_attach_gab_to_setting(setting,gabCS_lhs,gabCS_rhs)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
TYPE(LSTENSOR),pointer       :: gabCS_rhs,gabCS_lhs

IF(setting%scheme%CS_SCREEN)THEN
   setting%LST_GAB_RHS => gabCS_rhs
   setting%LST_GAB_LHS => gabCS_lhs
ENDIF
end subroutine ls_attach_gab_to_setting

subroutine ls_free_gab_from_setting(setting,lupri)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
INTEGER,INTENT(IN)            :: lupri
nullify(setting%LST_GAB_RHS)
nullify(setting%LST_GAB_LHS)
end subroutine ls_free_gab_from_setting

#ifdef VAR_MPI
SUBROUTINE ls_create_lstensor_task(setting,result_tensor,lstype,dmat_lhs,dmat_rhs,&
     & result_tensor_full,dmat_lhs_full,dmat_rhs_full,&
     & gabCS_rhs_full,gabCS_lhs_full,gabCS_rhs,gabCS_lhs,&
     & rhsCS_created,lhsCS_created,&
     & nbast1,nbast2,nbast3,nbast4,lhs_created,rhs_created,lhs,rhs,sameAOsLHS,sameAOsRHS,&
     & sameODs,sameAllFrag,CS_SCREEN,PS_SCREEN,doscreen)
implicit none
integer,intent(in)      :: nbast1,nbast2,nbast3,nbast4
character(len=7)        :: lstype
TYPE(LSSETTING)         :: setting
TYPE(LSTENSOR),pointer  :: result_tensor,dmat_lhs,dmat_rhs,result_tensor_full,dmat_lhs_full,dmat_rhs_full
TYPE(LSTASK),intent(IN) :: lhs,rhs
TYPE(LSTENSOR),pointer  :: gabCS_rhs_full,gabCS_lhs_full,gabCS_rhs,gabCS_lhs
Logical,intent(IN)      :: sameAOsLHS,sameAOsRHS,sameODs,lhs_created,rhs_created,sameAllFrag
Logical,intent(IN)      :: rhsCS_created,lhsCS_created,cs_screen,PS_SCREEN,doscreen
!
integer :: iao1,iao2
logical :: LinK
!integer                 :: n1,n2,n3,n4


nullify(dmat_lhs)
nullify(dmat_rhs)
nullify(gabCS_lhs)
nullify(gabCS_rhs)

LinK = .FALSE.

!result matrix
IF (lstype.EQ.'AB_TYPE') THEN
   call ls_create_lstensor_task_ab(result_tensor,result_tensor_full,lhs,rhs,nbast1,nbast2,sameAllFrag)
ELSE IF (lstype.EQ.'AC_TYPE') THEN
   call ls_create_lstensor_task_ab(result_tensor,result_tensor_full,lhs,rhs,nbast1,nbast3,sameAllFrag)
   LinK = .TRUE.
ELSE IF (lstype.EQ.'FULLINT') THEN
   call ls_create_lstensor_task_full(result_tensor,result_tensor_full,lhs,rhs,nbast1,nbast2,nbast3,nbast4,&
        &                              sameAllFrag)
ELSE
   write(*,*) 'Error in ls_create_lstensor_task. Wrong lstype: ',lstype
   call lsquit('Error in ls_create_lstensor_task. Wrong lstype',-1)
ENDIF
call lstensor_zero(result_tensor)
!density matrices
IF(lhs_created)THEN
   IF (LinK) THEN
     iao1 = 1
     iao2 = 2
   ELSE
     iao1 = setting%LHSdmatAOindex1
     iao2 = setting%LHSdmatAOindex2
   ENDIF
   call ls_create_lstensor_task1(dmat_lhs,dmat_lhs_full,lhs,rhs,nbast1,nbast2,nbast3,nbast4,iao1,iao2,sameAllFrag,.FALSE.)
ENDIF
IF(rhs_created)THEN
   IF (linK) THEN
     iao1 = 3
     iao2 = 4
   ELSE
     iao1 = setting%RHSdmatAOindex1
     iao2 = setting%RHSdmatAOindex2
   ENDIF
   call ls_create_lstensor_task1(dmat_rhs,dmat_rhs_full,lhs,rhs,nbast1,nbast2,nbast3,nbast4,iao1,iao2,sameAllFrag,.FALSE.)
ENDIF
!screening matrices
IF((CS_SCREEN.OR.PS_SCREEN).AND.doscreen)THEN
   call ls_create_lstensor_task1(gabCS_lhs,gabCS_lhs_full,lhs,rhs,nbast1,nbast2,nbast3,nbast4,1,2,sameAllFrag,LinK)
   call ls_create_lstensor_task1(gabCS_rhs,gabCS_rhs_full,lhs,rhs,nbast1,nbast2,nbast3,nbast4,3,4,sameAllFrag,LinK)
ENDIF
END SUBROUTINE ls_create_lstensor_task

SUBROUTINE ls_create_lstensor_task1(result_tensor,result_full,lhs,rhs,nbast1,nbast2,nbast3,nbast4,iao1,iao2,sameAllFrag,type13)
implicit none
integer,intent(in)      :: nbast1,nbast2,nbast3,nbast4,iao1,iao2
TYPE(LSTENSOR),pointer  :: result_tensor,result_full
TYPE(LSTASK),intent(IN) :: lhs,rhs
logical :: sameAllFrag,type13
!
integer :: Dummyatomlist1(1),Dummyatomlist2(1),i,natom1,natom2
integer,pointer :: rowatoms(:),colatoms(:)
integer :: natom(4),natom_full(4),nbast(4),n1,n2
logical :: full(4),full1,full2
type integer_pt 
integer,pointer :: p(:)
end type integer_pt
type(integer_pt) :: atom_list(4)
integer :: AO(4)

IF (type13) THEN
  AO(1) = 1
  AO(2) = 3
  AO(3) = 2
  AO(4) = 4
ELSE
  AO(1) = 1
  AO(2) = 2
  AO(3) = 3
  AO(4) = 4
ENDIF

full(AO(1)) = lhs%full_mol_row
full(AO(2)) = lhs%full_mol_col
full(AO(3)) = rhs%full_mol_row
full(AO(4)) = rhs%full_mol_col

nbast(AO(1)) = nbast1
nbast(AO(2)) = nbast2
nbast(AO(3)) = nbast3
nbast(AO(4)) = nbast4

natom(AO(1)) = lhs%nrow
natom(AO(2)) = lhs%ncol
natom(AO(3)) = rhs%nrow
natom(AO(4)) = rhs%ncol

natom_full(AO(1)) = lhs%nrow_full
natom_full(AO(2)) = lhs%ncol_full
natom_full(AO(3)) = rhs%nrow_full
natom_full(AO(4)) = rhs%ncol_full

atom_list(AO(1))%p => lhs%row_atoms
atom_list(AO(2))%p => lhs%col_atoms
atom_list(AO(3))%p => rhs%row_atoms
atom_list(AO(4))%p => rhs%col_atoms

natom1 = natom(iao1)
natom2 = natom(iao2)

n1 = nbast(iao1)
n2 = nbast(iao2)

full1 = full(iao1)
full2 = full(iao2)

nullify(result_tensor)
allocate(result_tensor)
call LSTENSOR_nullify(result_tensor)
IF(full1)THEN
   natom1 = natom_full(iao1)
   call mem_alloc(rowatoms,natom1)
   DO I=1,natom1
      rowatoms(I) = I
   ENDDO
ELSE
   rowatoms => atom_list(iao1)%p
ENDIF
IF(full2)THEN
   natom2 = natom_full(iao2)
   call mem_alloc(colatoms,natom2)
   DO I=1,natom2
      colatoms(I) = I
   ENDDO
ELSE
   colatoms => atom_list(iao2)%p
ENDIF
Dummyatomlist1(1) = 1
Dummyatomlist2(1) = 1
call build_sublstensor_from_full_lstensor(result_tensor,result_full,natom1,natom2,1,1,&
     & rowatoms,colatoms,Dummyatomlist1,Dummyatomlist2,n1,n2,1,1,sameAllFrag)
IF(full1)call mem_dealloc(rowatoms)
IF(full2)call mem_dealloc(colatoms)

END SUBROUTINE ls_create_lstensor_task1

SUBROUTINE ls_create_lstensor_task_ab(jmat,jmat_full,lhs,rhs,nbast1,nbast2,sameAllFrag)
implicit none
integer,intent(in)      :: nbast1,nbast2
TYPE(LSTENSOR),pointer  :: jmat,jmat_full
TYPE(LSTASK),intent(IN) :: lhs,rhs
logical :: sameAllFrag
!
integer,target :: Dummyatomlist1(1),Dummyatomlist2(1)
integer,pointer :: atoms1(:),atoms2(:),atoms3(:),atoms4(:)
integer :: natom1,natom2,natom3,natom4 ,i

nullify(jmat)
allocate(jmat)
call LSTENSOR_nullify(jmat)
call set_atomspointers(lhs,rhs,natom1,natom2,natom3,natom4,atoms1,atoms2,atoms3,atoms4,&
     &                 Dummyatomlist1,Dummyatomlist2,jmat_full%gradienttensor,jmat_full%gradienttensor)
call build_sublstensor_from_full_lstensor(jmat,jmat_full,nAtom1,nAtom2,nAtom3,nAtom4,&
     & atoms1,atoms2,atoms3,atoms4,nbast1,nbast2,1,1,sameAllFrag)
IF(lhs%full_mol_row)call mem_dealloc(atoms1)
IF(lhs%full_mol_col)call mem_dealloc(atoms2)
IF(jmat_full%gradienttensor)THEN
   IF(rhs%full_mol_row)call mem_dealloc(atoms3)
   IF(rhs%full_mol_col)call mem_dealloc(atoms4)
ENDIF
END SUBROUTINE ls_create_lstensor_task_ab

SUBROUTINE ls_create_lstensor_task_full(result_tensor,result_full,lhs,rhs,nbast1,nbast2,nbast3,nbast4,sameAllFrag)
implicit none
integer,intent(in)      :: nbast1,nbast2,nbast3,nbast4
TYPE(LSTENSOR),pointer  :: result_tensor,result_full
TYPE(LSTASK),intent(IN) :: lhs,rhs
logical :: sameAllFrag
!
integer,target :: Dummyatomlist1(1),Dummyatomlist2(1)
integer,pointer :: atoms1(:),atoms2(:),atoms3(:),atoms4(:)
integer :: natom1,natom2,natom3,natom4 ,i

nullify(result_tensor)
allocate(result_tensor)
call LSTENSOR_nullify(result_tensor)
call set_atomspointers(lhs,rhs,natom1,natom2,natom3,natom4,atoms1,atoms2,atoms3,atoms4,&
     &                 Dummyatomlist1,Dummyatomlist2,.TRUE.,result_full%gradienttensor)
call build_sublstensor_from_full_lstensor(result_tensor,result_full,nAtom1,nAtom2,nAtom3,nAtom4,&
     & atoms1,atoms2,atoms3,atoms4,nbast1,nbast2,nbast3,nbast4,sameAllFrag)
IF(lhs%full_mol_row)call mem_dealloc(atoms1)
IF(lhs%full_mol_col)call mem_dealloc(atoms2)
IF(rhs%full_mol_row)call mem_dealloc(atoms3)
IF(rhs%full_mol_col)call mem_dealloc(atoms4)
END SUBROUTINE ls_create_lstensor_task_full

subroutine ls_extract_and_annihilate_lstensor_task(setting,result_tensor,lstype,dmat_lhs,dmat_rhs,&
     & result_full,dmat_lhs_full,dmat_rhs_full,nbast1,nbast2,nbast3,nbast4,&
     & lhs_created,rhs_created,gabCS_rhs,gabCS_lhs,rhsCS_created,&
     & lhsCS_created,lhs,rhs,sameAOsLHS,sameAOsRHS,sameODs,SameAllFrag,CS_SCREEN,PS_SCREEN,doscreen)
implicit none
TYPE(LSSETTING),intent(INOUT) :: setting
Character(len=7),intent(IN) :: lstype
integer,intent(in)          :: nbast1,nbast2,nbast3,nbast4
TYPE(LSTENSOR),pointer      :: result_tensor,dmat_lhs,dmat_rhs,result_full,dmat_lhs_full,dmat_rhs_full
TYPE(LSTASK),intent(IN)     :: lhs,rhs
TYPE(LSTENSOR),pointer      :: gabCS_rhs,gabCS_lhs
Logical,intent(IN)          :: sameAOsLHS,sameAOsRHS,sameODs,lhs_created,rhs_created,SameAllFrag
Logical,intent(IN)          :: rhsCS_created,lhsCS_created,CS_SCREEN,PS_SCREEN,doscreen
!
integer,target :: Dummyatomlist1(1),Dummyatomlist2(1)
integer,pointer :: atoms1(:),atoms2(:),atoms3(:),atoms4(:)
integer :: natom1,natom2,natom3,natom4,i,n1,n2,n3,n4
logical :: full,grad

grad = result_full%gradienttensor
IF (lstype.EQ.'AB_TYPE') THEN
  full = grad
ELSE IF (lstype.EQ.'AC_TYPE') THEN
  full = grad
ELSE IF (lstype.EQ.'FULLINT') THEN
  full = .TRUE.
ELSE
  write(*,*) 'Error in ls_extract_and_annihilate_lstensor_task. Not a valid lstype: ',lstype
  call lsquit('Error in ls_extract_and_annihilate_lstensor_task. Not a valid lstype',-1)
ENDIF

IF (full) THEN
  n1 = nbast1 
  n2 = nbast2
  n3 = nbast3 
  n4 = nbast4 
ELSE
  IF (lstype.EQ.'AC_TYPE') THEN
    n1 = nbast1 
    n2 = nbast3
    n3 = 1 
    n4 = 1 
  ELSE
    n1 = nbast1 
    n2 = nbast2
    n3 = 1 
    n4 = 1 
  ENDIF
ENDIF
call set_atomspointers(lhs,rhs,natom1,natom2,natom3,natom4,atoms1,atoms2,atoms3,atoms4,&
     &                    Dummyatomlist1,Dummyatomlist2,full,grad)

IF(Setting%scheme%cs_int.OR.Setting%scheme%ps_int)then
   call add_sublstensor_to_full_lstensor(setting%output%ScreenTensor,result_full,nAtom1,nAtom2,nAtom3,nAtom4,&
        & atoms1,atoms2,atoms3,atoms4,n1,n2,n3,n4,sameAllFrag)
   call lstensor_free(setting%output%ScreenTensor)
   deallocate(setting%output%ScreenTensor)
   nullify(setting%output%ScreenTensor)
ELSE
   call add_sublstensor_to_full_lstensor(result_tensor,result_full,nAtom1,nAtom2,nAtom3,nAtom4,&
        & atoms1,atoms2,atoms3,atoms4,n1,n2,n3,n4,sameAllFrag)
ENDIF
call lstensor_free(result_tensor)
deallocate(result_tensor)
nullify(result_tensor)

IF(lhs%full_mol_row)call mem_dealloc(atoms1)
IF(lhs%full_mol_col)call mem_dealloc(atoms2)
IF(full)THEN
   IF(rhs%full_mol_row)call mem_dealloc(atoms3)
   IF(rhs%full_mol_col)call mem_dealloc(atoms4)
ENDIF


call ls_free_lstensors(dmat_lhs,dmat_rhs,lhs_created,rhs_created)
IF((CS_SCREEN.OR.PS_SCREEN).AND.doscreen)THEN
   Call ls_free_screeninglstensors(gabCS_rhs,gabCS_lhs,.TRUE.,.TRUE.)
ENDIF

end subroutine ls_extract_and_annihilate_lstensor_task

subroutine ls_extract_and_annihilate_lstensor_task2(jmat,dmat_lhs,dmat_rhs,&
     & jmat_full,dmat_lhs_full,dmat_rhs_full,nbast1,nbast2,nbast3,nbast4,&
     & lhs_created,rhs_created,lhs,rhs,sameAOsLHS,sameAOsRHS,sameODs,SameAllFrag)
implicit none
integer,intent(in)      :: nbast1,nbast2,nbast3,nbast4
TYPE(LSTENSOR),pointer  :: jmat,jmat_full
TYPE(LSTENSOR),pointer  :: dmat_lhs,dmat_rhs,dmat_lhs_full,dmat_rhs_full
TYPE(LSTASK),intent(IN) :: lhs,rhs
Logical,intent(IN)      :: sameAOsLHS,sameAOsRHS,sameODs,lhs_created,rhs_created,SameAllFrag
!
integer,target :: Dummyatomlist1(1),Dummyatomlist2(1)
integer,pointer :: atoms1(:),atoms2(:),atoms3(:),atoms4(:)
integer :: natom1,natom2,natom3,natom4 ,i

IF(lhs%full_mol_row.AND.lhs%full_mol_col)THEN
   call add_lstensor_to_lstensor(jmat,jmat_full)
ELSE
   call set_atomspointers(lhs,rhs,natom1,natom2,natom3,natom4,atoms1,atoms2,atoms3,atoms4,&
     &                 Dummyatomlist1,Dummyatomlist2,.FALSE.,jmat_full%gradienttensor)

   call add_sublstensor_to_full_lstensor(jmat,jmat_full,nAtom1,nAtom2,1,1,&
        & atoms1,atoms2,atoms3,atoms4,nbast1,nbast2,1,1,sameAllFrag)
   IF(lhs%full_mol_row)call mem_dealloc(atoms1)
   IF(lhs%full_mol_col)call mem_dealloc(atoms2)
ENDIF
call lstensor_free(jmat)
deallocate(jmat)
call ls_free_lstensors(dmat_lhs,dmat_rhs,lhs_created,rhs_created)

end subroutine ls_extract_and_annihilate_lstensor_task2

subroutine set_atomspointers(lhs,rhs,natom1,natom2,natom3,natom4,atoms1,atoms2,atoms3,atoms4,&
     &                       Dummyatomlist1,Dummyatomlist2,full,grad)
implicit none
TYPE(LSTASK),intent(IN) :: lhs,rhs
integer,target :: Dummyatomlist1(1),Dummyatomlist2(1)
integer,pointer :: atoms1(:),atoms2(:),atoms3(:),atoms4(:)
integer :: natom1,natom2,natom3,natom4,i
logical :: full,grad

IF(lhs%full_mol_row)THEN
   call mem_alloc(atoms1,lhs%nrow)
   DO I=1,lhs%nrow
      atoms1(I) = I
   ENDDO
ELSE
   atoms1 => lhs%row_atoms
ENDIF
IF(lhs%full_mol_col)THEN
   call mem_alloc(atoms2,lhs%ncol)
   DO I=1,lhs%ncol
      atoms2(I) = I
   ENDDO
ELSE
   atoms2 => lhs%col_atoms
ENDIF
nAtom1=lhs%nrow
nAtom2=lhs%ncol

IF(full)THEN
   nAtom3=rhs%nrow
   nAtom4=rhs%ncol
   IF(rhs%full_mol_row)THEN
      IF (grad) nAtom3 = rhs%nrow_full
      call mem_alloc(atoms3,nAtom3)
      DO I=1,nAtom3
         atoms3(I) = I
      ENDDO
   ELSE
      atoms3 => rhs%row_atoms
   ENDIF
   IF(rhs%full_mol_col)THEN
      IF (grad) nAtom4 = rhs%ncol_full
      call mem_alloc(atoms4,nAtom4)
      DO I=1,nAtom4
         atoms4(I) = I
      ENDDO
   ELSE
      atoms4 => rhs%col_atoms
   ENDIF
ELSE
   nAtom3=1
   nAtom4=1
   Dummyatomlist1(1) = 1
   Dummyatomlist2(1) = 1
   atoms3 => Dummyatomlist1
   atoms4 => Dummyatomlist2
ENDIF
end subroutine set_atomspointers

subroutine set_single_atomspointers(lhs,natom1,natom2,atoms1,atoms2)
implicit none
TYPE(LSTASK),intent(IN) :: lhs
integer,pointer :: atoms1(:),atoms2(:)
integer :: natom1,natom2,i

IF(lhs%full_mol_row)THEN
   call mem_alloc(atoms1,lhs%nrow)
   DO I=1,lhs%nrow
      atoms1(I) = I
   ENDDO
ELSE
   atoms1 => lhs%row_atoms
ENDIF
IF(lhs%full_mol_col)THEN
   call mem_alloc(atoms2,lhs%ncol)
   DO I=1,lhs%ncol
      atoms2(I) = I
   ENDDO
ELSE
   atoms2 => lhs%col_atoms
ENDIF
nAtom1=lhs%nrow
nAtom2=lhs%ncol
end subroutine set_single_atomspointers
#endif

!> \brief set the densities in the integralinput structure correctly
!> \author S. Reine and T. Kjaeergaard
!> \date 2010
!> \param intinput the integral input to be set 
!> \param SETTING Integral evalualtion settings containing the densities
SUBROUTINE ls_setDensityDimensions(intinput,setting,lupri)
implicit none
integer,intent(in) :: lupri
TYPE(INTEGRALINPUT) :: INTINPUT
TYPE(LSSETTING)   :: SETTING
!
INTEGER :: LHS1,LHS2,RHS1,RHS2,nbastLHS1,nbastLHS2,nbastRHS1,nbastRHS2
INTEGER :: nDmatLHS,nDmatRHS,idmat
LOGICAL :: useAO1,useAO2,TMPalloc,permuteLHS,permuteRHS
Integer :: subdim(4),substart(4),I,n1,n2,n3,n4,nbast(4)
real(realk),pointer :: DmatTMP(:,:,:),DmatTMP2(:,:,:)

IF ((INTINPUT%operator.EQ.KineticOperator) .AND.setting%RHSdfull)&
     &CALL LSQUIT('Error in ls_setDensityDimensions. Kinetic and RHS density',-1)

IF (setting%LHSdfull.OR.setting%LHSdmat)THEN
   INTINPUT%NDMAT_LHS = setting%nDmatLHS
   INTINPUT%uselst_DLHS = .TRUE.
ENDIF

IF (setting%RHSdfull.OR.setting%RHSdmat)THEN
   INTINPUT%NDMAT_RHS = setting%nDmatRHS
   INTINPUT%uselst_DRHS = .TRUE.
ENDIF

END SUBROUTINE ls_setDensityDimensions

SUBROUTINE SET_PROPINFO1(Oper,INPUT)
implicit none
!type(LSINTSCHEME) :: scheme
integer              :: Oper
TYPE(INTEGRALINPUT)  :: INPUT

!DEFAULT IS
!INPUT%PROP_ORIGIN(1) = D0
!INPUT%PROP_ORIGIN(2) = D0
!INPUT%PROP_ORIGIN(3) = D0
!INPUT%PropDerivEcoeff = .TRUE.
!INPUT%PropKineticEcoeff = .FALSE.
!INPUT%PropMomentEcoeff =.TRUE.
!INPUT%PropRequireBoys = -1
!Fix me: determine which of these can use INPUT%doMagScreen true.
INPUT%PropMaxD = 1
INPUT%PropMaxM = 1
INPUT%sphericalEcoeff = .FALSE.
INPUT%operator = Oper
INPUT%do_passes = .TRUE.

IF(Oper .EQ. OVERLAPOperator)THEN
   INPUT%PROPTYPE = 1
   INPUT%PropAnti = .FALSE.
ELSEIF(Oper .EQ. DIPLENOperator)THEN
   INPUT%PROPTYPE = 2
   INPUT%PropAnti = .FALSE.
ELSEIF(Oper .EQ. DIPVELOperator)THEN
   INPUT%PROPTYPE = 3
   INPUT%PropAnti = .TRUE.
ELSEIF(Oper .EQ. THETAOperator)THEN
   INPUT%PROPTYPE = 7
   INPUT%PropAnti = .FALSE.
   INPUT%PropMaxD = 0
   INPUT%PropMaxM = 2
   INPUT%PropDerivEcoeff = .FALSE.
   INPUT%PropMomentEcoeff =.TRUE.   
   INPUT%do_passes = .FALSE. !FIXME: Test that this is needed
ELSEIF(Oper .EQ. PSOOperator)THEN
   INPUT%operator = NucpotOperator
   INPUT%PROPTYPE = 10
   INPUT%PropAnti = .FALSE.
   INPUT%PropMaxD = 1
   INPUT%PropRequireBoys = 4
   INPUT%addtoIntegral = .FALSE.
   INPUT%PropDerivEcoeff = .TRUE.
   INPUT%PropMomentEcoeff =.TRUE.
   INPUT%do_passes = .FALSE. !FIXME: Test that this is needed
ELSEIF(Oper .EQ. ANGLONOperator)THEN
   INPUT%PROPTYPE = 17
   INPUT%PropAnti = .FALSE.
ELSEIF(Oper .EQ. ANGMOMOperator)THEN
   INPUT%PROPTYPE = 18
   INPUT%PropAnti = .TRUE.
ELSEIF(Oper .EQ. LONMOM1Operator)THEN
   INPUT%PROPTYPE = 19
   INPUT%PropAnti = .FALSE.
   INPUT%PropMaxD = 2
ELSEIF(Oper .EQ. NSTNOLOperator)THEN
   INPUT%operator = NucpotOperator
   INPUT%PROPTYPE = 26
   INPUT%PropAnti = .FALSE.
   INPUT%PropMaxD = 1
   INPUT%PropMaxM = 1
   INPUT%PropRequireBoys = 7
   INPUT%addtoIntegral = .FALSE.
   INPUT%PropDerivEcoeff = .FALSE.
   INPUT%PropMomentEcoeff =.TRUE.
   INPUT%do_passes = .FALSE. !FIXME: Test that this is needed
ELSEIF(Oper .EQ. NSTLONOperator)THEN
   INPUT%operator = NucpotOperator
   INPUT%PROPTYPE = 27
   INPUT%PropAnti = .FALSE.
   INPUT%PropMaxD = 1
   INPUT%PropMaxM = 1
   INPUT%PropRequireBoys = 7
   INPUT%addtoIntegral = .FALSE.
   INPUT%PropDerivEcoeff = .TRUE.!false?
   INPUT%PropMomentEcoeff =.TRUE.
   INPUT%do_passes = .FALSE. !FIXME: Test that this is needed
ELSEIF(Oper .EQ. DCM1Operator)THEN
   INPUT%PROPTYPE = 42
   INPUT%PropAnti = .TRUE.
   INPUT%PropMaxD = 1
   INPUT%PropMaxM = 2
   INPUT%addtoIntegral = .FALSE.
   INPUT%PropDerivEcoeff = .TRUE. !false
   INPUT%PropMomentEcoeff =.TRUE.
   INPUT%do_passes = .FALSE. !FIXME: Test that this is needed
ELSEIF(Oper .EQ. DCM2Operator)THEN
   INPUT%PROPTYPE = 43
   INPUT%PropAnti = .FALSE.
   INPUT%PropMaxD = 1
   INPUT%PropMaxM = 3
   INPUT%addtoIntegral = .FALSE.
   INPUT%PropDerivEcoeff = .TRUE.!false
   INPUT%PropMomentEcoeff =.TRUE.
   INPUT%do_passes = .FALSE. !FIXME: Test that this is needed
ELSEIF(Oper .EQ. ROTSTROperator)THEN
   INPUT%PROPTYPE = 55
   INPUT%PropAnti = .TRUE.
   INPUT%PropMaxD = 1
   INPUT%PropMaxM = 1
   INPUT%PropDerivEcoeff = .TRUE.
   INPUT%PropMomentEcoeff =.TRUE.   
   INPUT%do_passes = .FALSE. !FIXME: Test that this is needed
ELSEIF(Oper .EQ. LONMOM2Operator)THEN
   INPUT%operator = NucpotOperator
   INPUT%PROPTYPE = 190
   INPUT%PropAnti = .FALSE.
   INPUT%PropMaxD = 2
   INPUT%PropRequireBoys = 2
   INPUT%addtoIntegral = .TRUE. 
ELSEIF(Oper .EQ. ELPOTOperator)THEN
   INPUT%operator = NucpotOperator
   INPUT%PROPTYPE = 64
   INPUT%PropAnti = .FALSE.
   INPUT%PropMaxD = 0
   INPUT%PropMaxM = 0
   INPUT%PropMomentEcoeff =.FALSE.
   INPUT%PropRequireBoys = 0
   INPUT%addtoIntegral = .TRUE.
   INPUT%do_passes = .FALSE. !FIXME: Test that this is needed
ENDIF

END SUBROUTINE SET_PROPINFO1

SUBROUTINE SET_PROPINFO2(PropOper,INPUT)
implicit none
integer     :: PropOper
TYPE(INTEGRALINPUT)  :: INPUT

IF(PropOper .EQ. ANGLONOperator)THEN
   INPUT%sameLHSaos = .FALSE.
ELSEIF(PropOper .EQ. PSOOperator)THEN
   INPUT%sameLHSaos = .FALSE.
ELSEIF(PropOper .EQ. NSTNOLOperator)THEN
   INPUT%sameLHSaos = .FALSE.
ELSEIF(PropOper .EQ. LONMOM1Operator)THEN
   INPUT%sameLHSaos = .FALSE.
ELSEIF(PropOper .EQ. LONMOM2Operator)THEN
   INPUT%sameLHSaos = .FALSE.
ELSEIF(PropOper .EQ. NSTLONOperator)THEN
   INPUT%sameLHSaos = .FALSE.
ENDIF

END SUBROUTINE SET_PROPINFO2

#ifdef VAR_MPI
!> \brief Adds a sub gradient to a full gradient
!> \author T. Kjærgaard
!> \date 2010
!> \param GRAD the full gradient
!> \param SubGRAD the sub gradient
!> \param ndir number of directions (should be 3)
!> \param natoms number of atoms in molecule
!> \param ndirFrag number of directions for the fragments (should be 3)
!> \param natomsFrag number of atoms in fragements (can be more than natoms)
!> \param lhs The lhs task
!> \param rhs The rhs task
!> \param sameFrag Specifies if all fragments are identical or not
!> \param lupri The print unit number
!> Two possibilities:
!>   1) all molecules/fragments are identical (for the four AO's in the four-index integrals)
!>   2) the molecules are different (sameFrag is false)
!> In case 1) it is straightforward to add the calculated gradient elements into the
!> full molecular gradient. In case four, we need to add the gradient contributions for the
!> four different fragments separately; with the nAtoms1+nAtoms2+nAtoms3+nAtoms4 
!> elements ordered consecutively for the four fragments
SUBROUTINE AddSubGradientTensor(GRAD,ndir,natoms,SubGRAD,nDirFrag,natomsFrag,lhs,rhs,&
    &                           sameFrag,lupri)
implicit none
Logical,intent(IN)        :: sameFrag
TYPE(LSTASK),intent(IN)   :: lhs,rhs
Integer,intent(IN)        :: ndir,natoms,nDirFrag,natomsFrag,lupri
Real(realk),intent(INOUT) :: GRAD(ndir,natoms)
Real(realk),intent(IN)    :: SubGRAD(nDirFrag,natomsFrag) 
!
TYPE IP
Integer,pointer :: p(:)
END TYPE IP
TYPE(IP)    :: atoms(4)
Integer     :: iatom,ifull,ifrag,idir,nA(4),iAO,nAO
Logical     :: full(4)

IF(ndir.NE. 3)call lsquit('error in addSubGradientTensor 1',lupri)
IF(nDirFrag.NE. 3)call lsquit('error in addSubGradientTensor 2',lupri)

nA(1) = lhs%nrow
nA(2) = lhs%ncol
nA(3) = rhs%nrow
nA(4) = rhs%ncol
full(1) = lhs%full_mol_row
full(2) = lhs%full_mol_col
full(3) = rhs%full_mol_row
full(4) = rhs%full_mol_col
atoms(1)%p => lhs%row_atoms
atoms(2)%p => lhs%col_atoms
atoms(3)%p => rhs%row_atoms
atoms(4)%p => rhs%col_atoms

nAO = 1
IF (.NOT.sameFrag) nAO=4

ifrag = 0
! This is a dummy loop if all fragemtns are identical (sameFrag is true)
DO iAO=1,nAO
  !Loop over all atoms for which the gradient has been calculated
  DO iatom=1,nA(iAO)
    ifrag = ifrag + 1

    ifull = iatom
    ! If the fragment is not the full molecule we need the proper atomic index
    IF (.NOT.full(iAO)) ifull = atoms(iAO)%p(iatom)
    DO idir=1,3
      GRAD(idir,ifull) = GRAD(idir,ifull) + SubGRAD(iDir,ifrag) 
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE AddSubGradientTensor
#endif

SUBROUTINE ls_subScreenFromList(GABsub,GABfull,atoms1,atoms2,natom1,natom2,nbast1,nbast2,sameFrag)
implicit none
integer,intent(in)      :: natom1,natom2,nbast1,nbast2
TYPE(LSTENSOR),pointer  :: GABsub,GABfull
logical, intent(IN)     :: sameFrag
integer,intent(IN)      :: atoms1(natom1),atoms2(natom2)
!
integer :: Dummyatomlist1(1),Dummyatomlist2(1)

Dummyatomlist1(1) = 1
Dummyatomlist2(1) = 1
call build_sublstensor_from_full_lstensor(GABsub,GABfull,natom1,natom2,1,1,&
     & atoms1,atoms2,Dummyatomlist1,Dummyatomlist2,nbast1,nbast2,1,1,sameFrag)

END SUBROUTINE ls_subScreenFromList

SUBROUTINE ls_subScreenAtomic(GABsub,GABfull,iatom1,iatom2,nbast1,nbast2,sameFrag)
implicit none
integer,intent(in)      :: iatom1,iatom2,nbast1,nbast2
TYPE(LSTENSOR),pointer  :: GABsub,GABfull
logical, intent(IN)     :: sameFrag
!
integer :: Dummyatomlist1(1),Dummyatomlist2(1),atoms1(1),atoms2(1)

Dummyatomlist1(1) = 1
Dummyatomlist2(1) = 1
atoms1(1) = iatom1
atoms2(1) = iatom2
!nullify(GABsub)
!allocate(GABsub)
call LSTENSOR_nullify(GABsub)
call build_sublstensor_from_full_lstensor(GABsub,GABfull,1,1,1,1,&
     & atoms1,atoms2,Dummyatomlist1,Dummyatomlist2,nbast1,nbast2,1,1,sameFrag)

END SUBROUTINE ls_subScreenAtomic

END MODULE ls_Integral_Interface

#ifdef VAR_MPI
SUBROUTINE lsmpi_getIntegrals_masterToSlave(AO1,AO2,AO3,AO4,Oper,Spec,intType,geoOrder,SETTING,LUPRI,LUERR)
use lsmpi_op, only: mpicopy_setting
use lsparameters
use lsmpi_type, only: ls_mpiFinalizeBuffer, ls_mpiInitBuffer, ls_mpi_buffer, &
     & LSMPIBROADCAST
use infpar_module
use typedeftype, only: lssetting
implicit none
Integer,intent(in)            :: LUPRI,LUERR
integer                       :: AO1,AO2,AO3,AO4,Oper,intType,Spec,geoOrder
Type(LSSETTING),intent(inout) :: SETTING
!
Integer :: Input(10)

call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,setting%comm)
Input(1) = AO1
Input(2) = AO2
Input(3) = AO3
Input(4) = AO4
Input(5) = Oper
Input(6) = intType
Input(7) = spec
Input(8) = geoOrder
Input(9) = LUPRI
Input(10)= LUERR
CALL ls_mpi_buffer(Input,10,infpar%master)
CALL mpicopy_setting(setting,setting%comm,.FALSE.)
call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,setting%comm)

END SUBROUTINE lsmpi_getIntegrals_masterToSlave

SUBROUTINE lsmpi_getIntegrals_Slave(comm)
use lsmpi_op, only: mpicopy_setting
use lsmpi_type, only: ls_mpiFinalizeBuffer, ls_mpiInitBuffer, ls_mpi_buffer, &
     & LSMPIBROADCAST
use lsparameters
use infpar_module
use typedeftype, only: lssetting
use ls_Integral_Interface, only: ls_getIntegrals
implicit none
Integer         :: LUPRI,LUERR
integer(kind=ls_mpik) :: comm
integer:: AO1,AO2,AO3,AO4,Oper,intType,Spec,geoOrder
Type(LSSETTING) :: SETTING
!
Integer :: Input(10)
call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)

CALL ls_mpi_buffer(Input,10,infpar%master)
AO1 = Input(1)
AO2 = Input(2)
AO3 = Input(3)
AO4 = Input(4)
Oper = Input(5)
IntType = Input(6)
Spec = Input(7)
geoOrder = Input(8)
LUPRI = Input(9)
LUERR = Input(10)
CALL mpicopy_setting(setting,comm,.FALSE.)
call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)

call ls_getIntegrals(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR,geoOrder)

END SUBROUTINE lsmpi_getIntegrals_Slave

SUBROUTINE lsmpi_jengine_masterToSlave(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
use precision
use lsmpi_op, only: mpicopy_setting
use lsmpi_type, only: ls_mpiFinalizeBuffer, ls_mpiInitBuffer, ls_mpi_buffer, &
     & LSMPIBROADCAST
use lsparameters
use infpar_module
use typedeftype, only: lssetting
!use lstiming, only: lstimer
implicit none
Integer,intent(in)            :: LUPRI,LUERR
integer                 :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING),intent(inout) :: SETTING
!
real(realk) :: ts,te

!call lstimer('START',ts,te,lupri)

call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,setting%comm)
CALL ls_mpi_buffer(AO1,infpar%master)
CALL ls_mpi_buffer(AO2,infpar%master)
CALL ls_mpi_buffer(AO3,infpar%master)
CALL ls_mpi_buffer(AO4,infpar%master)
CALL ls_mpi_buffer(Oper,infpar%master)
CALL ls_mpi_buffer(Spec,infpar%master)
CALL ls_mpi_buffer(intType,infpar%master)
CALL ls_mpi_buffer(LUPRI,infpar%master)
CALL ls_mpi_buffer(LUERR,infpar%master)
!call lstimer('master1',ts,te,lupri)
CALL mpicopy_setting(setting,setting%comm,.FALSE.)
!call lstimer('master2',ts,te,lupri)
call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,setting%comm)
!call lstimer('master3',ts,te,lupri)

END SUBROUTINE lsmpi_jengine_masterToSlave

SUBROUTINE lsmpi_jengine_Slave(comm)
use lsmpi_op, only: mpicopy_setting
use lsmpi_type, only: ls_mpiFinalizeBuffer, ls_mpiInitBuffer, ls_mpi_buffer, &
     & LSMPIBROADCAST
use lsparameters
use infpar_module
use typedeftype, only: lssetting
use ls_Integral_Interface, only: ls_jengine
implicit none
Integer         :: LUPRI,LUERR
integer(kind=ls_mpik) :: comm
integer         :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING) :: SETTING
call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
CALL ls_mpi_buffer(AO1,infpar%master)
CALL ls_mpi_buffer(AO2,infpar%master)
CALL ls_mpi_buffer(AO3,infpar%master)
CALL ls_mpi_buffer(AO4,infpar%master)
CALL ls_mpi_buffer(Oper,infpar%master)
CALL ls_mpi_buffer(Spec,infpar%master)
CALL ls_mpi_buffer(intType,infpar%master)
CALL ls_mpi_buffer(LUPRI,infpar%master)
CALL ls_mpi_buffer(LUERR,infpar%master)
CALL mpicopy_setting(setting,comm,.FALSE.)
call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
call ls_jengine(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)

END SUBROUTINE lsmpi_jengine_Slave

SUBROUTINE lsmpi_LinK_masterToSlave(AO1,AO2,AO3,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)
use lsmpi_op, only: mpicopy_setting
use lsmpi_type, only: ls_mpiFinalizeBuffer, ls_mpiInitBuffer, ls_mpi_buffer, &
     & LSMPIBROADCAST
use lsparameters
use infpar_module
use typedeftype, only: lssetting
use ls_Integral_Interface, only: ls_get_exchange_mat
implicit none
Integer,intent(in)            :: LUPRI,LUERR
integer                       :: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING),intent(inout) :: SETTING
!
Integer :: Input(9)

call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,setting%comm)
Input(1) = AO1
Input(2) = AO2
Input(3) = AO3
Input(4) = AO4
Input(5) = Oper
Input(6) = intType
Input(7) = spec
Input(8) = LUPRI
Input(9) = LUERR
CALL ls_mpi_buffer(Input,9,infpar%master)
CALL mpicopy_setting(setting,setting%comm,.FALSE.)
call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,setting%comm)

END SUBROUTINE lsmpi_LinK_masterToSlave

SUBROUTINE lsmpi_LinK_Slave(comm)
use lsmpi_op, only: mpicopy_setting
use lsmpi_type, only: ls_mpiFinalizeBuffer, ls_mpiInitBuffer, ls_mpi_buffer, &
     & LSMPIBROADCAST
use infpar_module
use lsparameters
use typedeftype, only: lssetting
use ls_Integral_Interface, only: ls_get_exchange_mat
implicit none
Integer         :: LUPRI,LUERR
integer(kind=ls_mpik) :: comm
integer:: AO1,AO2,AO3,AO4,Oper,intType,Spec
Type(LSSETTING) :: SETTING
!
Integer :: Input(9)
call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)

CALL ls_mpi_buffer(Input,9,infpar%master)
AO1 = Input(1)
AO2 = Input(2)
AO3 = Input(3)
AO4 = Input(4)
Oper = Input(5)
IntType = Input(6)
Spec = Input(7)
LUPRI = Input(8)  ! KK: I believe this should be removed
LUERR = Input(9)  ! KK: I believe this should be removed
! KK QUICK FIX: It makes no sense to set lupri and luerr for the slave
! equal to lupri and luerr for the master!!!
! Thomas/Simen: Please fix this properly.
if(infpar%mynum /= infpar%master) then
   LUPRI=0
   LUERR=0
end if
CALL mpicopy_setting(setting,comm,.FALSE.)
call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)

call ls_get_exchange_mat(AO1,AO3,AO2,AO4,Oper,Spec,intType,SETTING,LUPRI,LUERR)

END SUBROUTINE lsmpi_LinK_Slave

#endif

