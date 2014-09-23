MODULE ProfileIchorMod
  use precision
  use TYPEDEFTYPE, only: LSSETTING, LSINTSCHEME, LSITEM, integralconfig,&
       & BASISSETLIBRARYITEM
  use basis_type, only: free_basissetinfo
  use basis_typetype,only: BASISSETINFO,BASISINFO,RegBasParam,nBasisBasParam
  use BuildBasisSet, only: Build_BASIS
  use Matrix_module, only: MATRIX, MATRIXP
  use Integralparameters
  use molecule_typetype, only: MOLECULEINFO, MOLECULE_PT, ATOMITEM
  use molecule_type, only: build_pointmolecule, DETERMINE_MAXCOOR, &
       & free_moleculeinfo, build_atomicmolecule, print_mol
  use integral_type, only: INTEGRALINPUT
  use TYPEDEF, only: getNbasis, retrieve_output, gcao2ao_transform_matrixd2, &
       & retrieve_screen_output, ao2gcao_transform_matrixf, &
       & gcao2ao_transform_fulld, ao2gcao_transform_fullf, &
       & ao2gcao_half_transform_matrix
  use molecule_module, only: DETERMINE_NBAST2
  use matrix_operations, only: mat_dotproduct, matrix_type, mtype_unres_dense,&
       & mat_daxpy, mat_init, mat_free, mat_write_to_disk, mat_print, mat_zero,&
       & mat_scal, mat_mul, mat_assign, mat_trans
  use matrix_util, only: mat_get_isym, util_get_symm_part,util_get_antisymm_part, matfull_get_isym
  use memory_handling, only: mem_alloc, mem_dealloc, debug_mem_stats
  use IntegralInterfaceMOD, only: II_get_4center_eri
  use IchorErimoduleHost!,only: MAIN_ICHORERI_DRIVER
  use lsmatrix_operations_dense
  use LSmatrix_type
  use configurationType
  use lstiming
  public::  profile_Ichor
  private
CONTAINS
!config%prof%Ichor
SUBROUTINE profile_Ichor(LUPRI,LUERR,SETTING,config)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
type(configItem)    :: config
#ifdef VAR_ICHOR
!
real(realk),pointer   :: integralsII(:,:,:,:),integralsIchor(:,:,:,:)
integer :: dim1,dim2,dim3,dim4,A,B,C,D,iprint,nbast(4),ibasiselm(4)
integer :: iBasis1,ibasis2,ibasis3,ibasis4,icharge,nbasis,nPass,ipass,itest,perm
logical :: dirac,doprint,debug
real(realk)         :: TIMSTR,TIMEND,normII,normIchorNoscreen,normIchorScreen
TYPE(BASISSETLIBRARYITEM) :: LIBRARY(nBasisBasParam)
CHARACTER(len=9)     :: BASISLABEL
TYPE(BASISINFO),pointer :: unittestBASIS(:)
TYPE(BASISINFO),pointer :: originalBASIS
CHARACTER(len=80)    :: BASISSETNAME
logical      :: spherical,savedospherical,SameMOL,COMPARE,ForcePrint,sameBAS(4,4)
logical :: MoTrans,NoSymmetry

Character    :: intSpec(5)
NoSymmetry = .FALSE. !activate permutational symmetry
MoTrans=.FALSE.
intSpec(1) = 'R'
intSpec(2) = 'R'
intSpec(3) = 'R'
intSpec(4) = 'R'
intSpec(5) = 'C'
ForcePrint = .TRUE.
WRITE(lupri,*)'IchorProfile Routine'
!call debug_mem_stats(LUPRI)
do A=1,80
   BASISSETNAME(A:A) = ' '
enddo
originalBasis => SETTING%BASIS(1)%p
dirac = .FALSE.
doprint = .FALSE.
allocate(UNITTESTBASIS(4))
BASISLABEL='REGULAR  '
spherical = .TRUE.!.FALSE.
IPRINT = 0
   
IF(config%prof%IchorProfInputBasis)THEN
   sameBAS = SETTING%sameBAS
   do A = 1,4       
      BASISSETNAME(1:20) = config%prof%IchorProfInputBasisString(A)
      WRITE(lupri,*)'Using Input Basis:',BASISSETNAME(1:20)
      CALL Build_basis(LUPRI,IPRINT,&
           &SETTING%MOLECULE(A)%p,UNITTESTBASIS(A)%BINFO(RegBasParam),LIBRARY,&
           &BASISLABEL,.FALSE.,.FALSE.,doprint,spherical,RegBasParam,BASISSETNAME)
      SETTING%BASIS(A)%p => UNITTESTBASIS(A)
      call determine_nbast2(SETTING%MOLECULE(A)%p,SETTING%BASIS(A)%p%BINFO(RegBasParam),spherical,.FALSE.,nbast(A))
   enddo
   do A = 1,4       
    do B = 1,4       
     SETTING%sameBAS(A,B)=(config%prof%IchorProfInputBasisString(A).EQ.config%prof%IchorProfInputBasisString(B))&
          &.AND.(nbast(A).EQ.nbast(B))
    enddo
   enddo
   dim1 = nbast(1); dim2 = nbast(2); dim3 = nbast(3); dim4 = nbast(4)
ELSE
   do A = 1,4       
      call determine_nbast2(SETTING%MOLECULE(A)%p,SETTING%BASIS(A)%p%BINFO(RegBasParam),spherical,.FALSE.,nbast(A))
   enddo
   dim1 = nbast(1); dim2 = nbast(2); dim3 = nbast(3); dim4 = nbast(4)
ENDIF
COMPARE=.FALSE.
IF(config%prof%IchorProfdoLink)THEN
   !generate random non symmetric DensityMatrix
!   QQQQQQQQQQQQQQQQQQQQQQQQQ
ELSE
   COMPARE = config%prof%IchorProfDoThermite.AND.config%prof%IchorProfDoIchor
ENDIF

WRITE(lupri,*)'Dims:',nbast(1:4)
IF(config%prof%IchorProfDoThermite)THEN
   IF(config%prof%IchorProfdoLink)THEN
!QQQQQQQQQQQQQQQQQQQQ
   ELSE
      COMPARE = .TRUE.
      WRITE(lupri,*)'Performing Thermite Profiling'
      call mem_alloc(integralsII,dim1,dim2,dim3,dim4)
      savedospherical = setting%scheme%dospherical
      setting%scheme%dospherical = spherical
         Setting%sameMol = .FALSE.
         Setting%sameFrag = .FALSE.
      !   setting%scheme%OD_SCREEN = .FALSE.
      !   setting%scheme%CS_SCREEN = .FALSE.
      !   setting%scheme%PS_SCREEN = .FALSE.
      CALL LSTIMER('START',TIMSTR,TIMEND,lupri)
      call II_get_4center_eri(LUPRI,LUERR,SETTING,integralsII,dim1,dim2,dim3,dim4,dirac)
      CALL LSTIMER('Thermite4Center',TIMSTR,TIMEND,lupri,ForcePrint)
      call determine_norm(IntegralsII,normII,dim1,dim2,dim3,dim4)
      WRITE(lupri,*)'Norm of Thermite:',normII
      setting%scheme%dospherical = savedospherical
      !   setting%scheme%OD_SCREEN = .TRUE.
      !   setting%scheme%CS_SCREEN = .TRUE.
      !   setting%scheme%PS_SCREEN = .TRUE.
         Setting%sameMol = .TRUE.
         Setting%sameFrag = .TRUE.
      IF(.NOT.COMPARE)THEN
         call mem_dealloc(integralsII)
      ENDIF
   ENDIF
ENDIF

IF(config%prof%IchorProfDoIchor)THEN
   IF(config%prof%IchorProfdoLink)THEN

   ELSE
      WRITE(lupri,*)'Performing Ichor Profiling'
      call mem_alloc(integralsIchor,dim1,dim2,dim3,dim4)
      integralsIchor = 0.0E0_realk
      CALL LSTIMER('START',TIMSTR,TIMEND,lupri)
      SameMOL = .TRUE.
      !   setting%scheme%OD_SCREEN = .FALSE.
      !   setting%scheme%CS_SCREEN = .FALSE.
      !   setting%scheme%PS_SCREEN = .FALSE.
      call SCREEN_ICHORERI_DRIVER(LUPRI,iprint,setting,INTSPEC,SameMOL)
      CALL LSTIMER('IchorScreen',TIMSTR,TIMEND,lupri,ForcePrint)
      CALL LSTIMER('START',TIMSTR,TIMEND,lupri)
      call MAIN_ICHORERI_DRIVER(LUPRI,iprint,setting,dim1,dim2,dim3,dim4,&
           & integralsIchor,intspec,.TRUE.,1,1,1,1,1,1,1,1,MoTrans,&
           & dim1,dim2,dim3,dim4,NoSymmetry)
      call FREE_SCREEN_ICHORERI
      !   setting%scheme%OD_SCREEN = .TRUE.
      !   setting%scheme%CS_SCREEN = .TRUE.
      !   setting%scheme%PS_SCREEN = .TRUE.
      CALL LSTIMER('Ichor4Center',TIMSTR,TIMEND,lupri,ForcePrint)
      call determine_norm(integralsIchor,normIchorScreen,dim1,dim2,dim3,dim4)
      WRITE(lupri,*)'Norm of Ichor4Center:',normIchorScreen
      IF(config%prof%IchorProfDoThermite.AND.(ABS(normIchorScreen-normII).LT.1.0E-12_realk))THEN
         WRITE(lupri,*)'Norm the same screen'
         print*,'Norm the same screen'
      ELSE
         print*,'normIchorScreen',normIchorScreen
         IF(config%prof%IchorProfDoThermite)THEN
            print*,'normII   ',normII
            print*,'ABS(normIchor-normII)',ABS(normIchorScreen-normII)
         ENDIF
      ENDIF
      IF(.NOT.COMPARE)THEN
         call mem_dealloc(integralsIchor)
      ENDIF
   ENDIF
ENDIF
IF(COMPARE)THEN
 DO D=1,dim4
   DO C=1,dim3
      DO B=1,dim2
         DO A=1,dim1
            WRITE(lupri,'(A,I4,A,I4,A,I4,A,I4,A,2F16.8)')'INT(',A,',',B,',',C,',',D,')',integralsIchor(a,b,c,d),integralsII(a,b,c,d)
         ENDDO
      ENDDO
   ENDDO
 ENDDO
 write(lupri,*)'integralsIchor:'
 call ls_output(integralsIchor,1,dim1*dim2,1,dim3*dim4,dim1*dim2,dim3*dim4,1,lupri)
 write(lupri,*)'integralsThermie:'
 call ls_output(integralsII,1,dim1*dim2,1,dim3*dim4,dim1*dim2,dim3*dim4,1,lupri)
 call mem_dealloc(integralsIchor)
 call mem_dealloc(integralsII)
ENDIF

IF(config%prof%IchorProfInputBasis)THEN
   call free_basissetinfo(UNITTESTBASIS(1)%BINFO(RegBasParam))
   call free_basissetinfo(UNITTESTBASIS(2)%BINFO(RegBasParam))
   call free_basissetinfo(UNITTESTBASIS(3)%BINFO(RegBasParam))
   call free_basissetinfo(UNITTESTBASIS(4)%BINFO(RegBasParam))
endif
deallocate(UNITTESTBASIS)
SETTING%BASIS(1)%p => originalBasis
SETTING%BASIS(2)%p => originalBasis
SETTING%BASIS(3)%p => originalBasis
SETTING%BASIS(4)%p => originalBasis
IF(config%prof%IchorProfInputBasis)THEN
   SETTING%sameBAS = sameBAS
ENDIF
WRITE(lupri,*)'Done IchorUnitTest'
!call debug_mem_stats(LUPRI)
#else
call lsquit('Profile_Ichor requires VAR_ICHOR',-1)
#endif
END SUBROUTINE Profile_Ichor

#ifdef VAR_ICHOR
subroutine determine_norm(Integrals,norm,dim1,dim2,dim3,dim4)
implicit none
integer :: dim1,dim2,dim3,dim4
real(realk) :: Integrals(dim1,dim2,dim3,dim4)
real(realk) :: norm
!
integer A,B,C,D
norm = 0.0E0_realk
DO D=1,dim4
   DO C=1,dim3
      DO B=1,dim2
         DO A=1,dim1
            norm = norm + Integrals(A,B,C,D)
         ENDDO
      ENDDO
   ENDDO
ENDDO
end subroutine determine_norm
#endif

End MODULE ProfileIchorMod

