MODULE HODItest_module
  use precision
  use TYPEDEFTYPE, only: LSSETTING
  use LSTIMING
  use memory_handling, only: mem_alloc, mem_dealloc
  use IntegralInterfaceMOD
  use matrix_module
  use matrix_operations, only: mat_dotproduct, matrix_type, mtype_unres_dense,&
       & mat_daxpy, mat_init, mat_free, mat_write_to_disk, mat_print, mat_zero,&
       & mat_scal, mat_mul, mat_assign, mat_trans, mat_copy, mat_add, mat_trAB,&
       & mat_sqnorm2,mat_tr,mat_max_elm

  public :: debugTestHODI
  private
CONTAINS

SUBROUTINE debugTestHODI(lupri,luerr,setting,D,nbast,nAtoms)
implicit none
Integer,intent(IN)            :: lupri,luerr,nbast,nAtoms
TYPE(LSSETTING),intent(INOUT) :: setting
TYPE(Matrix),intent(IN)       :: D
!
Integer,parameter        :: ncontract = 2
Integer                  :: i
TYPE(Matrix)             :: ContractMat1(ncontract)
Real(realk),parameter    :: D1=1.E0_realk
Real(realk),parameter    :: D2=2.E0_realk

DO i=1,ncontract
  call mat_init(ContractMat1(i),nbast,nbast)
  call mat_assign(ContractMat1(i),D)
  call mat_scal(D1*i,ContractMat1(i))
ENDDO

!******** First order *********
!Coulomb-type matrix - first-order geometrical derivative
call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,3,4,nbast,nAtoms,1,.FALSE.,D2,&
     &'first-order geometrical-derivative Coulomb matrix')

!Exchange-type matrix - first-order geometrical derivative
call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,2,4,nbast,nAtoms,1,.FALSE.,D1,&
     &'first-order geometrical-derivative exchange matrix')

!******** Second order *********
!Coulomb-type matrix - seconde-order geometrical derivative
call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,3,4,nbast,nAtoms,2,.FALSE.,D2,&
     &'second-order geometrical-derivative Coulomb matrix')

!!!Exchange-type matrix - seconde-order geometrical derivative
!!call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,2,4,nbast,nAtoms,2,.FALSE.,D1,&
!!     &'second-order geometrical-derivative exchange matrix')

!******** Third order *********
!Coulomb-type matrix - third-order geometrical derivative
call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,3,4,nbast,nAtoms,3,.FALSE.,D2,&
     &'third-order geometrical-derivative Coulomb matrix')

!!!Exchange-type matrix - third-order geometrical derivative
!!call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,2,4,nbast,nAtoms,3,.FALSE.,D1,&
!!     &'third-order geometrical-derivative exchange matrix')


!******** Fourth order *********
!!!Coulomb-type matrix - forth-order geometrical derivative
!!call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,3,4,nbast,nAtoms,4,.FALSE.,D2,&
!!     &'fourth-order geometrical-derivative Coulomb matrix')

!!!Exchange-type matrix - third-order geometrical derivative
!!call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,2,4,nbast,nAtoms,4,.FALSE.,D1,&
!!     &'fourth-order geometrical-derivative exchange matrix')


DO i=1,ncontract
  call mat_free(ContractMat1(i))
ENDDO

END SUBROUTINE debugTestHODI

SUBROUTINE debugTestHODIcontract1(lupri,luerr,setting,ContractMat1,ncontract,contract1,contract2,&
     &                            nbast,nAtoms,geoderiv,add,factor,txt)
implicit none
Integer,intent(IN)            :: lupri,luerr,ncontract,contract1,contract2,nbast,nAtoms,geoderiv
TYPE(LSSETTING),intent(INOUT) :: setting
TYPE(Matrix)                  :: ContractMat1(ncontract)
Logical,intent(IN)            :: add
real(realk),intent(IN)        :: factor
Character*(*)                 :: txt
!
TYPE(Matrix),allocatable :: ResultMat(:)
Integer                  :: nGeoderivComp,nc,i,j,k
real(realk),pointer      :: expectation(:)
real(realk),pointer      :: RMS(:)

IF (geoderiv.EQ.1) THEN
  nGeoderivComp = 3*nAtoms
ELSE IF (geoderiv.EQ.2) THEN
  nGeoderivComp = 9*nAtoms*nAtoms
ELSE IF (geoderiv.EQ.3) THEN
  nGeoderivComp = (3*nAtoms+1)*(3*nAtoms+2)*(3*nAtoms+3)/6
ELSE
  write(lupri,'(X,A,I3)') 'Error in debugTestHODIcontract1 - unknown case geoderiv=',geoderiv
ENDIF

nc=ncontract
IF (add) nc=1

ALLOCATE(ResultMat(nGeoderivComp*nc))

call mem_alloc(expectation,nGeoderivComp*nc)
call mem_alloc(RMS,nc)

DO i=1,nGeoderivComp*nc
  call mat_init(ResultMat(i),nbast,nbast)
  call mat_zero(ResultMat(i))
ENDDO

call II_get_hodi_eri_contract1(LUPRI,LUERR,SETTING,ResultMat,ContractMat1,ncontract,contract1,contract2,&
     &                         nbast,nbast,nbast,nbast,nGeoderivComp,geoderiv,add)

RMS(:) = 0E0_realk
i=0
write(lupri,'(1X,A,A)') 'HODI contract1 expectation values: ',txt
DO j=1,nc
  write(lupri,'(3X,A,I3)') 'Contraction number ',j
  DO k=1,nGeoderivComp
    i=i+1
    expectation(i) = factor*mat_trAB(ResultMat(i),ContractMat1(j))
    RMS(j) = RMS(j) + expectation(i)*expectation(i)
    write(lupri,'(5X,I3,F19.8)') k,expectation(i)
    call mat_free(ResultMat(i))
  ENDDO
  write(lupri,'(1X,A,A,X,F19.8)') 'RMS value for ',txt,sqrt(RMS(j)/nGeoderivComp)
ENDDO


call mem_dealloc(expectation)
DEALLOCATE(ResultMat)

END SUBROUTINE debugTestHODIcontract1

END MODULE HODItest_module
