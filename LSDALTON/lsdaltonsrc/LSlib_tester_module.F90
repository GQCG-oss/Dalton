MODULE LSlib_tester_mod
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
  use lslib_state, only: LSlibconfig

CONTAINS

SUBROUTINE LSlib_test_driver(OnMaster,lupri,luerr,meminfo_slaves)
  implicit none
   logical, intent(in) :: OnMaster
   integer, intent(inout) :: lupri, luerr
   logical, intent(out) :: meminfo_slaves
  Integer             :: nbast,natoms,nelectrons,i,j,k,l,n,m,o,x,y,z,iGrad,iHess,iCubic,ij,nDerivPacked
#ifdef LSLIB_RESTART
  type(matrix) :: D
  logical :: dens_exsist, DiskOnMaster=.true., gcbasis
  integer :: restart_lun
#else
  Integer,parameter   :: realk = 8
#endif
  Real(realk),pointer :: Smat(:,:),Dmat(:,:,:),TempMat(:,:,:),TempGrad(:,:,:),DFD(:,:,:),h1(:,:),Fmat(:,:,:)
  Real(realk),pointer :: TempHess(:,:,:,:,:)
  Real(realk),pointer :: TempCubic(:,:,:,:)
  Real(realk)         :: tmp1,tmp2,EXC(2),constant
  Real(realk),pointer :: eri(:,:,:,:,:)
  Integer,external    :: LSlib_get_nbasis
  logical :: diff,i1,j1,k1,l1

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
      call mat_read_from_disk(restart_lun,D,DiskOnMaster)
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
CALL LSlib_get_Fock(Fmat,Dmat,nbast,2,.TRUE.,lupri,luerr,.FALSE.,h1)

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
call ls_dzero(eri,nbast*nbast*nbast*nbast*3*nAtoms)

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
!write(*,*) 'Coulomb gradient contribution',j,(TempGrad(i,j,1),i=1,3)
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

#if 0
!*****************************************************************************
!******                             Coulomb Hessian from differentiated 4-center ERI (Mulliken)
!*****************************************************************************

deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,nbast,nbast,9*nAtoms*nAtoms))
call ls_dzero(eri,nbast*nbast*nbast*nbast*9*nAtoms*nAtoms)

CALL LSlib_get_4center_eri_geoderiv(eri,nbast,2,9*nAtoms*nAtoms,.FALSE.,lupri,luerr)

nullify(TempHess)
allocate(TempHess(3,nAtoms,3,nAtoms,1))
call ls_dzero(TempHess,natoms*3*nAtoms*3)

iHess = 0
DO n=1,nAtoms
  DO x=1,3
    DO m=1,nAtoms
      DO y=1,3
        iHess = iHess+1
        DO l=1,nbast
         DO k=1,nbast
          DO j=1,nbast
           DO i=1,nbast
            TempHess(y,m,x,n,1) = TempHess(y,m,x,n,1) + &
     &         2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iHess)*Dmat(k,l,1)
!           TempHess(y,m,x,n,1) = TempHess(y,m,x,n,1) + 0.25_realk*eri(i,j,k,l,iHess)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
iHess = 0
DO n=1,nAtoms
  DO x=1,3
    DO m=1,nAtoms
      DO y=1,3
        iHess = iHess+1
        tmp1=tmp1+TempHess(y,m,x,n,1)*TempHess(y,m,x,n,1)
        tmp2=tmp2+TempHess(y,m,x,n,1)*iHess
 write(*,'(A,I4,F21.9)') 'Hessian components:',iHess,TempHess(y,m,x,n,1)
      ENDDO
    ENDDO
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb Hessian from diff-eri (Mulliken): RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/natoms/9),tmp2/natoms/natoms/9
deallocate(TempHess)

#endif
#if 0
!*****************************************************************************
!******                             Coulomb third derivative from differentiated 4-center ERI (Mulliken)
!*****************************************************************************

deallocate(eri)
nullify(eri)
nDerivPacked = 3*nAtoms*(3*nAtoms+1)*(3*nAtoms+2)/6
allocate(eri(nbast,nbast,nbast,nbast,nDerivPacked))
call ls_dzero(eri,nbast*nbast*nbast*nbast*nDerivPacked)

CALL LSlib_get_4center_eri_geoderiv(eri,nbast,3,nDerivPacked,.FALSE.,lupri,luerr)

nullify(TempCubic)
allocate(TempCubic(3*nAtoms,3*nAtoms,3*nAtoms,1))
call ls_dzero(TempCubic,natoms*3*nAtoms*3*nAtoms*3)

iCubic = 0
DO n=1,3*nAtoms
  DO m=n,3*nAtoms
    DO o=m,3*nAtoms
      iCubic = iCubic+1
      DO l=1,nbast
       DO k=1,nbast
        DO j=1,nbast
         DO i=1,nbast
          TempCubic(o,m,n,1) = TempCubic(o,m,n,1) + &
     &       2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iCubic)*Dmat(k,l,1)
!         TempCubic(o,m,n,1) = TempCubic(o,m,n,1) + 0.25_realk*eri(i,j,k,l,iCubic)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      TempCubic(o,n,m,1) = TempCubic(o,m,n,1)
      TempCubic(m,o,n,1) = TempCubic(o,m,n,1)
      TempCubic(m,n,o,1) = TempCubic(o,m,n,1)
      TempCubic(n,m,o,1) = TempCubic(o,m,n,1)
      TempCubic(n,o,m,1) = TempCubic(o,m,n,1)
    ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
iCubic = 0
DO n=1,3*nAtoms
  DO m=1,3*nAtoms
    DO o=1,3*nAtoms
      iCubic = iCubic+1
      tmp1=tmp1+TempCubic(o,m,n,1)*TempCubic(o,m,n,1)
      tmp2=tmp2+TempCubic(o,m,n,1)*iCubic
      write(*,'(A,4I4,F21.9)') 'Cubic force components:',iCubic,o,m,n,TempCubic(o,m,n,1)
    ENDDO
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Coulomb Cubic force from diff-eri (Mulliken): RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/natoms/natoms/27),tmp2/natoms/natoms/natoms/27
deallocate(TempCubic)
#endif

!*****************************************************************************
!******                             First derivative nuclear-electron attraction integrals
!*****************************************************************************

deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,1,1,3*nAtoms))
call ls_dzero(eri,nbast*nbast*1*1*3*nAtoms)

CALL LSlib_get_1el_geoderiv(eri,'nucel',nbast,nAtoms,1,3*nAtoms,lupri,luerr)

call ls_dzero(TempGrad,natoms*3)

iGrad = 0
DO n=1,nAtoms
  DO x=1,3
   iGrad = iGrad+1
   DO l=1,1
    DO k=1,1
     DO j=1,nbast
      DO i=1,nbast
       TempGrad(x,n,1) = TempGrad(x,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iGrad)
!      TempGrad(x,n,1) = TempGrad(x,n,1) + eri(i,j,k,l,iGrad)
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
!write(*,*) 'debug:TempGrad',i,j,TempGrad(i,j,1)
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Nuclear-electron attraction gradient: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


!*****************************************************************************
!******                             Second derivative nuclear-electron attraction integrals
!*****************************************************************************

nullify(TempHess)
allocate(TempHess(3,nAtoms,3,nAtoms,1))
deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,1,1,9*nAtoms*nAtoms))
call ls_dzero(eri,nbast*nbast*1*1*9*nAtoms*nAtoms)

CALL LSlib_get_1el_geoderiv(eri,'nucel',nbast,nAtoms,2,9*nAtoms*nAtoms,lupri,luerr)

call ls_dzero(TempHess,natoms*3*nAtoms*3)

iHess = 0
DO n=1,nAtoms
  DO x=1,3
    DO m=1,nAtoms
      DO y=1,3
        iHess = iHess+1
        DO l=1,1
         DO k=1,1
          DO j=1,nbast
           DO i=1,nbast
            TempHess(y,m,x,n,1) = TempHess(y,m,x,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iHess)
!           TempHess(y,m,x,n,1) = TempHess(y,m,x,n,1) + eri(i,j,k,l,iHess)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
iHess = 0
DO n=1,nAtoms
  DO x=1,3
    DO m=1,nAtoms
      DO y=1,3
        iHess = iHess+1
        tmp1=tmp1+TempHess(y,m,x,n,1)*TempHess(y,m,x,n,1)
        tmp2=tmp2+TempHess(y,m,x,n,1)*iHess
!write(*,*) 'debug:TempHess',y,m,x,n,TempHess(y,m,x,n,1)
      ENDDO
    ENDDO
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Nuclear-electron attraction Hessian: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/natoms/9),tmp2/natoms/natoms/9

#if 0
!*****************************************************************************
!******                             Third derivative nuclear-electron attraction integrals 
!*****************************************************************************

deallocate(eri)
nullify(eri)
nDerivPacked = 3*nAtoms*(3*nAtoms+1)*(3*nAtoms+2)/6
allocate(eri(nbast,nbast,1,1,nDerivPacked))
call ls_dzero(eri,nbast*nbast*1*1*nDerivPacked)

CALL LSlib_get_1el_geoderiv(eri,'nucel',nbast,nAtoms,3,nDerivPacked,lupri,luerr)

nullify(TempCubic)
allocate(TempCubic(3*nAtoms,3*nAtoms,3*nAtoms,1))
call ls_dzero(TempCubic,natoms*3*nAtoms*3*nAtoms*3)

iCubic = 0
DO n=1,3*nAtoms
  DO m=n,3*nAtoms
    DO o=m,3*nAtoms
      iCubic = iCubic+1
      DO l=1,1
        DO k=1,1
          DO j=1,nbast
            DO i=1,nbast
            TempCubic(o,m,n,1) = TempCubic(o,m,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iCubic)*Dmat(k,l,1)
!           TempCubic(o,m,n,1) = TempCubic(o,m,n,1) + eri(i,j,k,l,iCubic)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      TempCubic(o,n,m,1) = TempCubic(o,m,n,1)
      TempCubic(m,o,n,1) = TempCubic(o,m,n,1)
      TempCubic(m,n,o,1) = TempCubic(o,m,n,1)
      TempCubic(n,m,o,1) = TempCubic(o,m,n,1)
      TempCubic(n,o,m,1) = TempCubic(o,m,n,1)
    ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
iCubic = 0
DO n=1,3*nAtoms
  DO m=1,3*nAtoms
    DO o=1,3*nAtoms
        iCubic = iCubic+1
        tmp1=tmp1+TempCubic(o,m,n,1)*TempCubic(o,m,n,1)
        tmp2=tmp2+TempCubic(o,m,n,1)*iCubic
!write(*,*) 'debug:TempCubic',o,m,n,TempCubic(o,m,n,1)
    ENDDO
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Cubic nuclear-electron attraction derivative integrals: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/natoms/natoms/27),tmp2/natoms/natoms/natoms/27
deallocate(TempCubic)
#endif

!*****************************************************************************
!******                             First derivative overlap integrals
!*****************************************************************************

deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,1,1,3*nAtoms))
call ls_dzero(eri,nbast*nbast*1*1*3*nAtoms)

CALL LSlib_get_1el_geoderiv(eri,'overlap',nbast,nAtoms,1,3*nAtoms,lupri,luerr)

call ls_dzero(TempGrad,natoms*3)

iGrad = 0
DO n=1,nAtoms
  DO x=1,3
   iGrad = iGrad+1
   DO l=1,1
    DO k=1,1
     DO j=1,nbast
      DO i=1,nbast
!      TempGrad(x,n,1) = TempGrad(x,n,1) + DFD(i,j,1)*eri(i,j,k,l,iGrad)
       TempGrad(x,n,1) = TempGrad(x,n,1) + Dmat(i,j,1)*eri(i,j,k,l,iGrad)
       TempGrad(x,n,1) = TempGrad(x,n,1) + eri(i,j,k,l,iGrad)
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
write(lupri,'(A80,2F18.10)') 'Reorthonormalization gradient using 1el_diff: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3

!*****************************************************************************
!******                             Second derivative overlap integrals
!*****************************************************************************

deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,1,1,9*nAtoms*nAtoms))
call ls_dzero(eri,nbast*nbast*1*1*9*nAtoms*nAtoms)

CALL LSlib_get_1el_geoderiv(eri,'overlap',nbast,nAtoms,2,9*nAtoms*nAtoms,lupri,luerr)

call ls_dzero(TempHess,natoms*3*nAtoms*3)

iHess = 0
DO n=1,nAtoms
  DO x=1,3
    DO m=1,nAtoms
      DO y=1,3
        iHess = iHess+1
        DO l=1,1
         DO k=1,1
          DO j=1,nbast
           DO i=1,nbast
            TempHess(y,m,x,n,1) = TempHess(y,m,x,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iHess)
!           TempHess(y,m,x,n,1) = TempHess(y,m,x,n,1) + eri(i,j,k,l,iHess)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
iHess = 0
DO n=1,nAtoms
  DO x=1,3
    DO m=1,nAtoms
      DO y=1,3
        iHess = iHess+1
        tmp1=tmp1+TempHess(y,m,x,n,1)*TempHess(y,m,x,n,1)
        tmp2=tmp2+TempHess(y,m,x,n,1)*iHess
      ENDDO
    ENDDO
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Overlap Hessian: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/natoms/9),tmp2/natoms/natoms/9


#if 0
!*****************************************************************************
!******                             Third derivative overlap integrals 
!*****************************************************************************

deallocate(eri)
nullify(eri)
nDerivPacked = 3*nAtoms*(3*nAtoms+1)*(3*nAtoms+2)/6
allocate(eri(nbast,nbast,1,1,nDerivPacked))
call ls_dzero(eri,nbast*nbast*1*1*nDerivPacked)

CALL LSlib_get_1el_geoderiv(eri,'overlap',nbast,nAtoms,3,nDerivPacked,lupri,luerr)

nullify(TempCubic)
allocate(TempCubic(3*nAtoms,3*nAtoms,3*nAtoms,1))
call ls_dzero(TempCubic,natoms*3*nAtoms*3*nAtoms*3)

iCubic = 0
DO n=1,3*nAtoms
  DO m=n,3*nAtoms
    DO o=m,3*nAtoms
      iCubic = iCubic+1
      DO l=1,1
       DO k=1,1
        DO j=1,nbast
         DO i=1,nbast
          TempCubic(o,m,n,1) = TempCubic(o,m,n,1) + &
   &         2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iCubic)*Dmat(k,l,1)
!         TempCubic(o,m,n,1) = TempCubic(o,m,n,1) + eri(i,j,k,l,iCubic)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      TempCubic(o,n,m,1) = TempCubic(o,m,n,1)
      TempCubic(m,o,n,1) = TempCubic(o,m,n,1)
      TempCubic(m,n,o,1) = TempCubic(o,m,n,1)
      TempCubic(n,m,o,1) = TempCubic(o,m,n,1)
      TempCubic(n,o,m,1) = TempCubic(o,m,n,1)
    ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
iCubic = 0
DO n=1,3*nAtoms
  DO m=1,3*nAtoms
    DO o=1,3*nAtoms
      iCubic = iCubic+1
      tmp1=tmp1+TempCubic(o,m,n,1)*TempCubic(o,m,n,1)
      tmp2=tmp2+TempCubic(o,m,n,1)*iCubic
    ENDDO
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Cubic overlap derivative integrals: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/natoms/natoms/27),tmp2/natoms/natoms/natoms/27
deallocate(TempCubic)
#endif

!*****************************************************************************
!******                             First derivative kinetic energy integrals
!*****************************************************************************

deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,1,1,3*nAtoms))
call ls_dzero(eri,nbast*nbast*1*1*3*nAtoms)

CALL LSlib_get_1el_geoderiv(eri,'kinetic',nbast,nAtoms,1,3*nAtoms,lupri,luerr)

call ls_dzero(TempGrad,natoms*3)

iGrad = 0
DO n=1,nAtoms
  DO x=1,3
   iGrad = iGrad+1
   DO l=1,1
    DO k=1,1
     DO j=1,nbast
      DO i=1,nbast
       TempGrad(x,n,1) = TempGrad(x,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iGrad)
!      TempGrad(x,n,1) = TempGrad(x,n,1) + eri(i,j,k,l,iGrad)
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
write(lupri,'(A80,2F18.10)') 'Kinetic gradient using 1el_diff: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/3),tmp2/natoms/3


!*****************************************************************************
!******                             Second derivative kinetic energy integrals
!*****************************************************************************

deallocate(eri)
nullify(eri)
allocate(eri(nbast,nbast,1,1,9*nAtoms*nAtoms))
call ls_dzero(eri,nbast*nbast*1*1*9*nAtoms*nAtoms)

CALL LSlib_get_1el_geoderiv(eri,'kinetic',nbast,nAtoms,2,9*nAtoms*nAtoms,lupri,luerr)

call ls_dzero(TempHess,9*nAtoms*nAtoms)

iHess = 0
DO n=1,nAtoms
  DO x=1,3
    DO m=1,nAtoms
      DO y=1,3
       iHess = iHess+1
       DO l=1,1
        DO k=1,1
         DO j=1,nbast
          DO i=1,nbast
    !      TempHess(y,m,x,n,1) = TempHess(y,m,x,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iHess)
           TempHess(y,m,x,n,1) = TempHess(y,m,x,n,1) + eri(i,j,k,l,iHess)
          ENDDO
         ENDDO
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
    DO m=1,nAtoms
      DO y=1,3
        ij=ij+1
        tmp1=tmp1+TempHess(y,m,i,j,1)*TempHess(y,m,i,j,1)
        tmp2=tmp2+TempHess(y,m,i,j,1)*ij
      ENDDO
    ENDDO
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Kinetic Hessian using 1el_diff: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/natoms/9),tmp2/natoms/natoms/9


#if 0
!*****************************************************************************
!******                             Third derivative kinetic energy integrals 
!*****************************************************************************

deallocate(eri)
nullify(eri)
nDerivPacked = 3*nAtoms*(3*nAtoms+1)*(3*nAtoms+2)/6
allocate(eri(nbast,nbast,1,1,nDerivPacked))
call ls_dzero(eri,nbast*nbast*1*1*nDerivPacked)

CALL LSlib_get_1el_geoderiv(eri,'kinetic',nbast,nAtoms,3,nDerivPacked,lupri,luerr)

nullify(TempCubic)
allocate(TempCubic(3*nAtoms,3*nAtoms,3*nAtoms,1))
call ls_dzero(TempCubic,natoms*3*nAtoms*3*nAtoms*3)

iCubic = 0
DO n=1,3*nAtoms
  DO m=n,3*nAtoms
    DO o=m,3*nAtoms
      iCubic = iCubic+1
      DO l=1,1
       DO k=1,1
        DO j=1,nbast
         DO i=1,nbast
          TempCubic(o,m,n,1) = TempCubic(o,m,n,1) + 2.0_realk*Dmat(i,j,1)*eri(i,j,k,l,iCubic)*Dmat(k,l,1)
!         TempCubic(o,m,n,1) = TempCubic(o,m,n,1) + eri(i,j,k,l,iCubic)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      TempCubic(o,n,m,1) = TempCubic(o,m,n,1)
      TempCubic(m,o,n,1) = TempCubic(o,m,n,1)
      TempCubic(m,n,o,1) = TempCubic(o,m,n,1)
      TempCubic(n,m,o,1) = TempCubic(o,m,n,1)
      TempCubic(n,o,m,1) = TempCubic(o,m,n,1)
    ENDDO
  ENDDO
ENDDO
   
tmp1 = 0.0_realk
tmp2 = 0.0_realk
iCubic = 0
DO n=1,3*nAtoms
  DO m=1,3*nAtoms
    DO o=1,3*nAtoms
      iCubic = iCubic+1
      tmp1=tmp1+TempCubic(o,m,n,1)*TempCubic(o,m,n,1)
      tmp2=tmp2+TempCubic(o,m,n,1)*iCubic
    ENDDO
  ENDDO
ENDDO
write(lupri,'(A80,2F18.10)') 'Cubic kinetic energy derivative integrals: RMS and index-weighted sum',&
     &                     sqrt(tmp1/natoms/natoms/natoms/27),tmp2/natoms/natoms/natoms/27
deallocate(TempCubic)
#endif

!*****************************************************************************
!******                             First derivative kinetic energy integrals
!*****************************************************************************

deallocate(eri)
deallocate(TempHess)
deallocate(TempGrad)
deallocate(TempMat)
deallocate(h1)
deallocate(Fmat)
deallocate(Dmat)
deallocate(Smat)
deallocate(DFD)

meminfo_slaves = LSlibconfig%mpi_mem_monitor

write(lupri,'(A)') ''
write(lupri,'(A)') '*** LSlib tester completed ***'
write(lupri,'(A)') ''

END SUBROUTINE LSLIB_TEST_DRIVER

END MODULE LSlib_tester_mod
