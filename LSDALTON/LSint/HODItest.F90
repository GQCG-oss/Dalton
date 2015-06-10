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
!Coulomb-type matrix - third-order geometrical derivative, single contraction
call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,3,4,nbast,nAtoms,3,.FALSE.,D2,&
     &'third-order geometrical-derivative Coulomb matrix')

!Exchange-type matrix - third-order geometrical derivative
call debugTestHodiContract1(LUPRI,LUERR,SETTING,ContractMat1,ncontract,2,4,nbast,nAtoms,3,.FALSE.,D1,&
     &'third-order geometrical-derivative exchange matrix')

!Coulomb-type expectation - third-order geometrical derivative, double contraction
call debugTestHodiContract2(LUPRI,LUERR,SETTING,ContractMat1,ContractMat1,ncontract,1,2,3,4,nbast,nAtoms,3,.FALSE.,D2,&
     &'third-order geometrical-derivative Coulomb expectation')

!Exchange-type expectation - third-order geometrical derivative, double contraction
call debugTestHodiContract2(LUPRI,LUERR,SETTING,ContractMat1,ContractMat1,ncontract,1,3,2,4,nbast,nAtoms,3,.FALSE.,D1,&
     &'third-order geometrical-derivative exchange expectation')

!Overlap-type matrix - third-order geometrical derivative
call debugTestHodiOne(LUPRI,LUERR,SETTING,'overlap',ContractMat1,ncontract,nbast,nAtoms,3,.FALSE.,D1,&
     &'third-order overlap derivative matrix')

!Nuclear-electron repulsion-type matrix - third-order geometrical derivative
call debugTestHodiOne(LUPRI,LUERR,SETTING,'nucel',ContractMat1,ncontract,nbast,nAtoms,3,.FALSE.,D1,&
     &'third-order nuclear-electron repulsion derivative matrix')

!kinetic-integral type matrix - third-order geometrical derivative
call debugTestHodiOne(LUPRI,LUERR,SETTING,'kinetic',ContractMat1,ncontract,nbast,nAtoms,3,.FALSE.,D1,&
     &'third-order kinetic-integral derivative matrix')

!Overlap-type matrix type expectation - third-order geometrical derivative
call debugTestHODIoneContract(LUPRI,LUERR,SETTING,'overlap',ContractMat1,ncontract,nbast,nAtoms,3,.FALSE.,D1,&
     &'third-order overlap derivative expectation')

!Nuclear-electron repulsion type expectation - third-order geometrical derivative
call debugTestHODIoneContract(LUPRI,LUERR,SETTING,'nucel',ContractMat1,ncontract,nbast,nAtoms,3,.FALSE.,D1,&
     &'third-order nuclear-electron repulsion derivative expectation')

!kinetic-integral type expectation - third-order geometrical derivative
call debugTestHODIoneContract(LUPRI,LUERR,SETTING,'kinetic',ContractMat1,ncontract,nbast,nAtoms,3,.FALSE.,D1,&
     &'third-order kinetic-integral derivative expectation')

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
  nGeoderivComp = (3*nAtoms)*(3*nAtoms+1)*(3*nAtoms+2)/6
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
call mem_dealloc(RMS)
DEALLOCATE(ResultMat)

END SUBROUTINE debugTestHODIcontract1

SUBROUTINE debugTestHODIcontract2(lupri,luerr,setting,ContractMat1,ContractMat2,ncontract,contract1,contract2,&
     &                            contract3,contract4,nbast,nAtoms,geoderiv,add,factor,txt)
implicit none
Integer,intent(IN)            :: lupri,luerr,ncontract,contract1,contract2,contract3,contract4,nbast,nAtoms,geoderiv
TYPE(LSSETTING),intent(INOUT) :: setting
TYPE(Matrix)                  :: ContractMat1(ncontract),ContractMat2(ncontract)
Logical,intent(IN)            :: add
real(realk),intent(IN)        :: factor
Character*(*)                 :: txt
!
Integer                  :: nGeoderivComp,nc,i,j,k
real(realk),pointer      :: expectation(:)
real(realk),pointer      :: RMS(:)

IF (geoderiv.EQ.1) THEN
  nGeoderivComp = 3*nAtoms
ELSE IF (geoderiv.EQ.2) THEN
  nGeoderivComp = 9*nAtoms*nAtoms
ELSE IF (geoderiv.EQ.3) THEN
  nGeoderivComp = (3*nAtoms)*(3*nAtoms+1)*(3*nAtoms+2)/6
ELSE
  write(lupri,'(X,A,I3)') 'Error in debugTestHODIcontract2 - unknown case geoderiv=',geoderiv
ENDIF

nc=ncontract
IF (add) nc=1

ALLOCATE(expectation(nGeoderivComp*nc))
call ls_dzero(expectation,nGeoderivComp*nc)

call mem_alloc(RMS,nc)

call II_get_hodi_eri_contract2(LUPRI,LUERR,SETTING,expectation,ContractMat1,ContractMat2,ncontract,contract1,contract2,&
     &                         contract3,contract4,nbast,nbast,nbast,nbast,nGeoderivComp,geoderiv,add)

call ls_dzero(RMS,nc)
i=0
write(lupri,'(1X,A,A)') 'HODI contract1 expectation values: ',txt
DO j=1,nc
  write(lupri,'(3X,A,I3)') 'Contraction number ',j
  DO k=1,nGeoderivComp
    i=i+1
    RMS(j) = RMS(j) + factor*factor*expectation(i)*expectation(i)
    write(lupri,'(5X,I3,F19.8)') k,factor*expectation(i)
  ENDDO
  write(lupri,'(1X,A,A,X,F19.8)') 'RMS value for ',txt,sqrt(RMS(j)/nGeoderivComp)
ENDDO


call mem_dealloc(RMS)

END SUBROUTINE debugTestHODIcontract2

SUBROUTINE debugTestHODIone(lupri,luerr,setting,oneElType,ContractMat,ncontract,nbast,nAtoms,geoderiv,add,factor,txt)
implicit none
Integer,intent(IN)            :: lupri,luerr,nbast,nAtoms,geoderiv,ncontract
TYPE(Matrix),intent(IN)       :: ContractMat(ncontract)
TYPE(LSSETTING),intent(INOUT) :: setting
Logical,intent(IN)            :: add
real(realk),intent(IN)        :: factor
Character*(*)                 :: txt,oneElType
!
TYPE(Matrix),pointer :: oneElMat(:)
real(realk),pointer  :: expectation(:),RMS(:)
Integer :: n,m,mn,nGeoderivComp,nc


IF (geoderiv.EQ.1) THEN
  nGeoderivComp = 3*nAtoms
ELSE IF (geoderiv.EQ.2) THEN
  nGeoderivComp = 9*nAtoms*nAtoms
ELSE IF (geoderiv.EQ.3) THEN
  nGeoderivComp = (3*nAtoms)*(3*nAtoms+1)*(3*nAtoms+2)/6
ELSE
  write(lupri,'(X,A,I3)') 'Error in debugTestHODIone - unknown case geoderiv=',geoderiv
ENDIF

nc=ncontract
IF (add) nc=1
call mem_alloc(expectation,nc*nGeoderivComp)
call mem_alloc(RMS,nc)
call ls_dzero(expectation,nc*nGeoderivComp)
call ls_dzero(RMS,nc)

allocate(oneElMat(nGeoderivComp))
DO n=1,nGeoderivComp
  CALL mat_init(oneElMat(n),nbast,nbast)
  CALL mat_zero(oneElMat(n))
ENDDO

call II_get_hodi_1el_mat(LUPRI,LUERR,SETTING,oneElMat,oneElType,nbast,nbast,nGeoderivComp,geoderiv)

mn=0
DO m=1,ncontract
  DO n=1,nGeoderivComp
    mn=mn+1
    IF (add) mn=n
    expectation(mn)= expectation(mn) + factor*mat_trab(ContractMat(m),oneElMat(n))
  ENDDO
ENDDO

mn=0
DO m=1,nc
  DO n=1,nGeoderivComp
    mn=mn+1
    write(lupri,'(5X,I3,F19.8)') mn,expectation(mn)
    RMS(m)= RMS(m) + expectation(mn)*expectation(mn)
  ENDDO
  write(lupri,'(1X,A,A,X,F19.8)') 'RMS value for ',txt,sqrt(RMS(m)/nGeoderivComp)
ENDDO

call mem_dealloc(expectation)
call mem_dealloc(RMS)
DO n=1,nGeoderivComp
  CALL mat_free(oneElMat(n))
ENDDO


!CALL II_get_hodi_1el_contract(LUPRI,LUERR,SETTING,expectation,ContractMat,ncontract,oneElType,dim1,dim2,dim5,geoderiv,add)

END SUBROUTINE debugTestHODIone

SUBROUTINE debugTestHODIoneContract(lupri,luerr,setting,oneElType,ContractMat,ncontract,nbast,nAtoms,geoderiv,add,factor,txt)
implicit none
Integer,intent(IN)            :: lupri,luerr,nbast,nAtoms,geoderiv,ncontract
TYPE(Matrix),intent(IN)       :: ContractMat(ncontract)
TYPE(LSSETTING),intent(INOUT) :: setting
Logical,intent(IN)            :: add
real(realk),intent(IN)        :: factor
Character*(*)                 :: txt,oneElType
!
real(realk),pointer  :: expectation(:),RMS(:)
Integer :: n,m,mn,nGeoderivComp,nc

IF (geoderiv.EQ.1) THEN
  nGeoderivComp = 3*nAtoms
ELSE IF (geoderiv.EQ.2) THEN
  nGeoderivComp = 9*nAtoms*nAtoms
ELSE IF (geoderiv.EQ.3) THEN
  nGeoderivComp = (3*nAtoms)*(3*nAtoms+1)*(3*nAtoms+2)/6
ELSE
  write(lupri,'(X,A,I3)') 'Error in debugTestHODIone - unknown case geoderiv=',geoderiv
ENDIF

nc=ncontract
IF (add) nc=1

call mem_alloc(expectation,nc*nGeoderivComp)
call mem_alloc(RMS,nc)
call ls_dzero(expectation,nc*nGeoderivComp)
call ls_dzero(RMS,nc)

CALL II_get_hodi_1el_contract(LUPRI,LUERR,SETTING,expectation,ContractMat,ncontract,oneElType,&
     &                        nbast,nbast,nGeoderivComp,geoderiv,add)
 
mn=0
DO m=1,nc
  DO n=1,nGeoderivComp
    mn=mn+1
    write(lupri,'(5X,I3,F19.8)') mn,expectation(mn)
    RMS(m)= RMS(m) + expectation(mn)*expectation(mn)
  ENDDO
  write(lupri,'(1X,A,A,X,F19.8)') 'RMS value for ',txt,sqrt(RMS(m)/nGeoderivComp)
ENDDO

call mem_dealloc(expectation)
call mem_dealloc(RMS)



END SUBROUTINE debugTestHODIoneContract

END MODULE HODItest_module
