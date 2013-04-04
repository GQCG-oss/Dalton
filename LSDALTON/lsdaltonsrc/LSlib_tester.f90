!#define LSLIB_RESTART
PROGRAM lslib_test
use files
use lsmpi_type, only: lsmpi_finalize
#ifdef LSLIB_RESTART
  use configuration
  use TYPEDEF
  use TYPEDEFTYPE
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
  use integralinterfaceMod
  use memory_handling
#endif
implicit none
Integer             :: nbast,natoms,nelectrons,lupri,luerr,i,j,k,l,n,x,iGrad,ij
#ifdef LSLIB_RESTART
  type(matrix) :: D
  logical :: dens_exsist, OnMaster=.true., gcbasis
  integer :: restart_lun
#else
Integer,parameter   :: realk = 8
#endif
Real(realk),pointer :: Smat(:,:),Dmat(:,:,:),TempMat(:,:,:),TempGrad(:,:,:),DFD(:,:,:),h1(:,:),Fmat(:,:,:)
Real(realk)         :: tmp1,tmp2,EXC(2),constant
Real(realk),pointer :: eri(:,:,:,:,:)
Integer,external    :: LSlib_get_nbasis

call lsinit_all()

LUPRI=-1
LUERR=-1
CALL LSOPEN(LUPRI,'LSDALTON.OUT','NEW','FORMATTED')
CALL LSOPEN(LUERR,'LSDALTON.ERR','UNKNOWN','FORMATTED')

CALL LSlib_get_dimensions(nbast,natoms,nelectrons,lupri,luerr)
write(lupri,'(A,I8,A,I8,A,I8)') 'Starting lslib_test with nbast =',nbast,', natoms =', natoms,&
     &                          ' and nelectrons =', nelectrons


!* Allocations
nullify(h1)
nullify(TempMat)
nullify(Smat)
nullify(Dmat)
nullify(Fmat)
nullify(TempGrad)
nullify(DFD)
nullify(eri)
allocate(h1(nbast,nbast))
allocate(TempMat(nbast,nbast,2))
allocate(Smat(nbast,nbast))
allocate(Dmat(nbast,nbast,2))
allocate(Fmat(nbast,nbast,2))
allocate(DFD(nbast,nbast,2))
allocate(TempGrad(3,nAtoms,2))
allocate(eri(nbast,nbast,nbast,nbast,1))

#ifdef LSLIB_RESTART
   ! Get inital density for file (requires a full calculation using lsdalton.x prior to
   ! the lslib_tester.x excecution)
   INQUIRE(file='dens.restart',EXIST=dens_exsist) 
   if (dens_exsist) then
      call mat_init(D,nbast,nbast)
      restart_lun = -1  !initialization
      call lsopen(restart_lun,'dens.restart','OLD','UNFORMATTED')
      rewind restart_lun
      call mat_read_from_disk(restart_lun,D,OnMaster)
      call mat_read_info_from_disk(restart_lun,gcbasis)
      call lsclose(restart_lun,'KEEP')
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '*** RESTART FROM DENSITY ON DISK - READ FROM dens.restart  ***'
      WRITE(LUPRI,*)
      WRITE(*,*)
      WRITE(*,*) '*** RESTART FROM DENSITY ON DISK - READ FROM dens.restart  ***'
      WRITE(*,*)
      if (gcbasis) then
         WRITE(LUPRI,*) 'Your dens.restart was constructed using the grand-canonical (GC) basis,'
         WRITE(LUPRI,*) 'while LSlib_build_exchange uses the standard basis. '
         WRITE(LUPRI,*) 'The GC basis is default when running TRILEVEL or ATOMS.' 
         WRITE(LUPRI,*) 'Contruct a new dens.restart with '
         WRITE(LUPRI,*) '.START'
         WRITE(LUPRI,*) 'H1DIAG'
         call lsquit('Calculation in standard basis, dens.restart in GC basis!',lupri)
      endif
      CALL DCOPY(nbast*nbast,D%elms,1,Dmat(1,1,1),1)
      CALL DCOPY(nbast*nbast,D%elms,1,Dmat(1,1,2),1)
      call mat_free(D)
   else
      call lsquit('Calculation using LSlib_tester requires a dens.restart!',lupri)
   endif
#endif
!*****************************************************************************
!******                             OVERLAP MATRIX
!*****************************************************************************

CALL LSlib_get_overlap(Smat,nbast,lupri,luerr)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
!write(*,*) 'Overlap matrix'
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+Smat(i,j)*Smat(i,j)
    tmp2=tmp2+Smat(i,j)*ij
  ENDDO
  !write(*,'(15F12.5)') (Smat(i,j),i=1,nbast)
ENDDO
write(lupri,'(A80,2F18.10)') 'AO overlap matrix: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast
#ifndef LSLIB_RESTART
! Default is to use AO overlap matrix as density matrix
CALL DCOPY(nbast*nbast,Smat,1,Dmat(1,1,1),1)
CALL DCOPY(nbast*nbast,Smat,1,Dmat(1,1,2),1)

#endif
CALL DSCAL(nbast*nbast,2.0_realk,Dmat(1,1,2),1)
!*****************************************************************************
!******                             1-electron MATRIX
!*****************************************************************************

CALL LSlib_get_h1(h1,nbast,lupri,luerr)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+h1(i,j)*h1(i,j)
    tmp2=tmp2+h1(i,j)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'One-electron AO-matrix matrix: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

!*****************************************************************************
!******                             Coulomb MATRIX
!*****************************************************************************

! Get the Coulomb matrix
CALL LSlib_get_Coulomb(Fmat,Dmat,nbast,2,lupri,luerr)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+Fmat(i,j,1)*Fmat(i,j,1)
    tmp2=tmp2+Fmat(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb AO-matrix matrix number 1: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+Fmat(i,j,2)*Fmat(i,j,2)
    tmp2=tmp2+Fmat(i,j,2)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb AO-matrix matrix number 2: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

!*****************************************************************************
!******                             Exchange MATRIX
!*****************************************************************************

! Get the Exchange matrix
CALL LSlib_get_Exchange(TempMat,Dmat,nbast,2,.TRUE.,lupri,luerr,.FALSE.)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+TempMat(i,j,1)*TempMat(i,j,1)
    tmp2=tmp2+TempMat(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Exchange AO-matrix matrix number 1: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+TempMat(i,j,2)*TempMat(i,j,2)
    tmp2=tmp2+TempMat(i,j,2)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Exchange AO-matrix matrix number 2: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

call daxpy(2*nbast*nbast,1.0_realk,TempMat,1,Fmat,1)

!*****************************************************************************
!******                             XC MATRIX
!*****************************************************************************

! Get the XC matrix (not a valid AO density matrix so do not test #of electrons)
CALL LSlib_get_XC(TempMat,Dmat,EXC,nbast,2,lupri,luerr,.FALSE.)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+TempMat(i,j,1)*TempMat(i,j,1)
    tmp2=tmp2+TempMat(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'XC AO-matrix matrix number 1: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+TempMat(i,j,2)*TempMat(i,j,2)
    tmp2=tmp2+TempMat(i,j,2)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'XC AO-matrix matrix number 2: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

call daxpy(2*nbast*nbast,1.0_realk,TempMat,1,Fmat,1)
call daxpy(nbast*nbast,1.0_realk,h1,1,Fmat,1)
call daxpy(nbast*nbast,1.0_realk,h1,1,Fmat(1,1,2),1)

!*****************************************************************************
!******                             Fock MATRIX
!*****************************************************************************

! Get the Fock matrix (not a valid AO density matrix so do not test # of electrons)
!ToDo Make II_get_Fock work with ADMM
!CALL LSlib_get_Fock(TempMat,Dmat,nbast,2,.TRUE.,lupri,luerr,.FALSE.)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+Fmat(i,j,1)*Fmat(i,j,1)
    tmp2=tmp2+Fmat(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Fock/KS AO-matrix matrix number 1: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+Fmat(i,j,2)*Fmat(i,j,2)
    tmp2=tmp2+Fmat(i,j,2)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Fock/KS AO-matrix matrix number 2: RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

! -2DFD
CALL DGEMM('N','N',nbast,nbast,nbast,-2.0_realk,Dmat(1,1,1),nbast,Fmat(1,1,1),nbast,&
     &     0.0_realk,TempMat,nbast)
CALL DGEMM('N','N',nbast,nbast,nbast,1.0_realk,TempMat,nbast,Dmat(1,1,1),nbast,&
     &     0.0_realk,DFD(1,1,1),nbast)
CALL DGEMM('N','N',nbast,nbast,nbast,-2.0_realk,Dmat(1,1,2),nbast,Fmat(1,1,1),nbast,&
     &     0.0_realk,TempMat,nbast)
CALL DGEMM('N','N',nbast,nbast,nbast,1.0_realk,TempMat,nbast,Dmat(1,1,2),nbast,&
     &     0.0_realk,DFD(1,1,2),nbast)


!*****************************************************************************
!******                             4-center eri integrals (Mulliken notation)
!*****************************************************************************
CALL LSlib_get_4center_eri(eri,nbast,.FALSE.,lupri,luerr)

CALL LS_DZERO(TempMat,nbast*nbast)
DO l=1,nbast
  DO k=1,nbast
   DO j=1,nbast
     DO i=1,nbast
       TempMat(i,j,1) = TempMat(i,j,1) + 2.0_realk*eri(i,j,k,l,1)*Dmat(k,l,1)
     ENDDO
   ENDDO
 ENDDO
ENDDO
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+TempMat(i,j,1)*TempMat(i,j,1)
    tmp2=tmp2+TempMat(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb AO-matrix matrix from eri (Mulliken): RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

CALL LS_DZERO(TempMat,nbast*nbast)
DO l=1,nbast
  DO k=1,nbast
   DO j=1,nbast
     DO i=1,nbast
       TempMat(i,k,1) = TempMat(i,k,1) - eri(i,j,k,l,1)*Dmat(j,l,1)
     ENDDO
   ENDDO
 ENDDO
ENDDO
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+TempMat(i,j,1)*TempMat(i,j,1)
    tmp2=tmp2+TempMat(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Exchange AO-matrix matrix from eri (Mulliken): RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

!*****************************************************************************
!******                             4-center eri integrals (Dirac notation)
!*****************************************************************************
CALL LSlib_get_4center_eri(eri,nbast,.TRUE.,lupri,luerr)

CALL LS_DZERO(TempMat,nbast*nbast)
DO l=1,nbast
  DO k=1,nbast
   DO j=1,nbast
     DO i=1,nbast
       TempMat(i,j,1) = TempMat(i,j,1) + 2.0_realk*eri(i,k,j,l,1)*Dmat(k,l,1)
     ENDDO
   ENDDO
 ENDDO
ENDDO
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+TempMat(i,j,1)*TempMat(i,j,1)
    tmp2=tmp2+TempMat(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb AO-matrix matrix from eri (Dirac): RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

CALL LS_DZERO(TempMat,nbast*nbast)
DO l=1,nbast
  DO k=1,nbast
   DO j=1,nbast
     DO i=1,nbast
       TempMat(i,k,1) = TempMat(i,k,1) - eri(i,k,j,l,1)*Dmat(j,l,1)
     ENDDO
   ENDDO
 ENDDO
ENDDO
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,nbast
  DO i=1,nbast
    ij=ij+1
    tmp1=tmp1+TempMat(i,j,1)*TempMat(i,j,1)
    tmp2=tmp2+TempMat(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Exchange AO-matrix matrix from eri (Dirac): RMS and index-weighted sum',&
     &                     sqrt(tmp1/nbast/nbast),tmp2/nbast/nbast

!For all gradient contributions except the XC gradient we need to include the factor two 
!in the density matrix for closed shell systems
!CALL DSCAL(nbast*nbast,2.0_realk,Dmat(1,1,1),1)
!CALL DSCAL(nbast*nbast,2.0_realk,Dmat(1,1,2),1)
!*****************************************************************************
!******                             nn gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_nn_gradient(TempGrad(1,1,1),natoms,lupri,luerr)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'nn gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


!*****************************************************************************
!******                             kinetic gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_kinetic_gradient(TempGrad,Dmat,nbast,natoms,lupri,luerr)
CALL DSCAL(natoms*3,2.0_realk,TempGrad,1)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'kinetic gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


!*****************************************************************************
!******                             ne gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_ne_gradient(TempGrad,Dmat,nbast,natoms,lupri,luerr)
CALL DSCAL(natoms*3,2.0_realk,TempGrad,1)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'ne gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


!*****************************************************************************
!******                             1-electron gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_oneElectron_gradient(TempGrad,Dmat,nbast,natoms,lupri,luerr)
CALL DSCAL(natoms*3,2.0_realk,TempGrad,1)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') '1-electron gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


!*****************************************************************************
!******                             Coulomb gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_J_gradient(TempGrad,Dmat,Dmat,nbast,2,2,natoms,lupri,luerr)
CALL DSCAL(natoms*3,4.0_realk,TempGrad,1)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3



!*****************************************************************************
!******                             Exchange gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_K_gradient(TempGrad,Dmat,Dmat,nbast,2,2,natoms,lupri,luerr)
CALL DSCAL(natoms*3,4.0_realk,TempGrad,1)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Exchange gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3



!*****************************************************************************
!******                             2-electron gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_twoElectron_gradient(TempGrad,Dmat,Dmat,nbast,2,2,natoms,lupri,luerr,.FALSE.)
CALL DSCAL(natoms*3,4.0_realk,TempGrad,1)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') '2-electron gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


!For all gradient contributions except the XC gradient we need to include the factor two 
!in the density matrix for closed shell systems - rescale to get factor 1
!CALL DSCAL(nbast*nbast,0.5_realk,Dmat(1,1,1),1)
!CALL DSCAL(nbast*nbast,0.5_realk,Dmat(1,1,2),1)
!*****************************************************************************
!******                             Exchange-correlation gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_xc_gradient(TempGrad,Dmat,nbast,natoms,lupri,luerr,.FALSE.)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'XC gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3



!*****************************************************************************
!******                             Reorthonormalization gradient
!*****************************************************************************

call ls_dzero(TempGrad,natoms*3)
CALL LSlib_get_reorthoNormalization(TempGrad,DFD,nbast,1,natoms,lupri,luerr)

tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Reorthonormalization gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


!*****************************************************************************
!******                             Coulomb gradient from differentiated 4-center ERI (Mulliken)
!*****************************************************************************

deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,nbast,nbast,3*nAtoms))

CALL LSlib_get_4center_eri_geoderiv(eri,nbast,1,3*nAtoms,.FALSE.,lupri,luerr)

call ls_dzero(TempGrad,natoms*3)

iGrad = 0
DO n=1,nAtoms
  DO x=1,3
   iGrad = iGrad+1
   DO l=1,nbast
    DO k=1,nbast
     DO j=1,nbast
      DO i=1,nbast
       TempGrad(x,n,1) = TempGrad(x,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iGrad)*Dmat(k,l,1)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
!write(*,*) 'Coulomb gradient components',j,(TempGrad(i,j,1),i=1,3)
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb gradient from diff-eri (Mulliken): RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


#if 0
!*****************************************************************************
!******                             Coulomb gradient from differentiated 4-center ERI (Dirac)
!*****************************************************************************

deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,nbast,nbast,3*nAtoms))

CALL LSlib_get_4center_eri_geoderiv(eri,nbast,1,3*nAtoms,.TRUE.,lupri,luerr)

call ls_dzero(TempGrad,natoms*3)

iGrad = 0
DO n=1,nAtoms
  DO x=1,3
   iGrad = iGrad+1
   DO l=1,nbast
    DO k=1,nbast
     DO j=1,nbast
      DO i=1,nbast
       TempGrad(x,n,1) = TempGrad(x,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,k,j,l,iGrad)*Dmat(k,l,1)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
ij=0
DO j=1,natoms
  DO i=1,3
    ij=ij+1
    tmp1=tmp1+TempGrad(i,j,1)*TempGrad(i,j,1)
    tmp2=tmp2+TempGrad(i,j,1)*ij
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb gradient from diff-eri (Dirac): RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3
#endif

deallocate(eri)
deallocate(TempGrad)
deallocate(TempMat)
deallocate(h1)
deallocate(Fmat)
deallocate(Dmat)
deallocate(Smat)
deallocate(DFD)

call lsfree_all()
call lsmpi_finalize(lupri,.FALSE.)

write(lupri,'(A)') ''
write(lupri,'(A)') '*** LSlib tester completed ***'
write(lupri,'(A)') ''

CALL LSCLOSE(LUPRI,'KEEP')
CALL LSCLOSE(LUERR,'KEEP')

END PROGRAM lslib_test
