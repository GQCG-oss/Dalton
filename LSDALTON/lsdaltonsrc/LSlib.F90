!> @file
!> Contains general library routines for integral evaluation

!> \brief Used to initialize LSlib functionality- must be called prior to other LSlib calls
!> \author S. Reine
!> \date 2014-01-22
!> \param lupri The print unit
!> \param luerr The error unit
SUBROUTINE LSlib_init(OnMaster,lupri,luerr)
use precision
use lslib_state
use configuration, only: config_set_default_config, config_read_input
use daltoninfo, only: ls_init
implicit none
logical,intent(OUT)       :: OnMaster
integer,intent(IN)        :: lupri,luerr
!
logical :: dodft

call lsinit_all(OnMaster,lupri,luerr,tstart,tend)

IF (OnMaster) THEN
  call config_set_default_config(LSlibconfig)
  call config_read_input(LSlibconfig,lupri,luerr)
  doDFT = LSlibconfig%opt%calctype.EQ.LSlibconfig%opt%dftcalc
  call ls_init(ls,lupri,luerr,nbasis,LSlibconfig%integral,dodft,.false.,.false.)
  state_set = .TRUE.
ENDIF

END SUBROUTINE LSlib_init

SUBROUTINE LSlib_free(OnMaster,lupri,luerr,meminfo_slaves)
use precision
use lslib_state
use configuration, only: config_shutdown, config_free
use daltoninfo, only: ls_free
implicit none
logical,intent(IN) :: OnMaster
integer,intent(IN) :: lupri,luerr
Logical,optional   :: meminfo_slaves
logical :: memslave

IF (OnMaster) THEN
  call ls_free(ls)
  call config_shutdown(LSlibconfig)
  call config_free(LSlibconfig)
  state_set = .FALSE.
ENDIF
memslave = .FALSE.
IF (present(meminfo_slaves)) memslave = meminfo_slaves
call lsfree_all(OnMaster,lupri,luerr,tstart,tend,memslave)

END SUBROUTINE LSlib_free

#ifdef VAR_MPI
SUBROUTINE LSlib_set_external_comm(external_comm)
use infpar_module, only : call_mpi_init
use lsmpi_param,    only : MPI_COMM_LSDALTON
use precision,     only : ls_mpik
implicit none
integer(kind=ls_mpik),intent(IN) :: external_comm

MPI_COMM_LSDALTON = external_comm
call_mpi_init = .FALSE.

END SUBROUTINE LSlib_set_external_comm

SUBROUTINE LSlib_exit_slave()
use infpar_module,      only : infpar
use lsmpi_type,         only : ls_mpibcast
use LSparameters, only : LSMPIQUIT
use lslib_state
implicit none

CALL ls_mpibcast(LSMPIQUIT,infpar%master,ls%setting%comm)

END SUBROUTINE LSlib_exit_slave
#endif

!> \brief Returns the number of basis functions
!> \author S. Reine
!> \date 2013-01-20
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_dimensions(nbast,natoms,nelectrons,lupri,luerr)
  use precision
  use TYPEDEF  
  use TYPEDEFTYPE  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
  use integralinterfaceMod
  use memory_handling
  use lslib_state
IMPLICIT NONE
Integer,intent(OUT) :: nbast,natoms,nelectrons
Integer,intent(IN)  :: lupri,luerr
!
IF (.NOT.state_set) CALL LSQUIT('LSlib_get_dimensions error: LSlib not initialized',lupri)
!
nbast      = nbasis
natoms     = ls%input%MOLECULE%nAtoms
nelectrons = ls%input%MOLECULE%nelectrons

END SUBROUTINE LSlib_get_dimensions

!> \brief Calculates the overlap AO-matrix
!> \author S. Reine
!> \date 2013-01-20
!> \param S The overlap AO-matrix
!> \param nbast The number of oribtals
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_overlap(S,nbast,lupri,luerr)
  use precision
  use Matrix_module
  use Matrix_Operations
  use integralinterfaceMod
  use memory_handling
  use lslib_state
IMPLICIT NONE
Integer,intent(in)      :: nbast,lupri,luerr
Real(realk),intent(out) :: S(nbast,nbast)
!
TYPE(MATRIX),target :: h
Integer             :: mtype_save, mynum, nodtot, ierr

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_overlap error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_overlap. Basis-function mismatch',lupri)

!this will most likely fail for MPI ask Thomas (or call call mat_select_type)
mtype_save = matrix_type
matrix_type = mtype_dense
CALL mat_init(h,nbast,nbast)
CALL mat_zero(h)

CALL II_get_overlap(lupri,luerr,ls%setting,h)

call dcopy(nbast*nbast,h%elms,1,S,1)

CALL mat_free(h)
matrix_type = mtype_save

END SUBROUTINE LSlib_get_overlap

!> \brief Calculates the Fock/Kohn-Sham AO-matrix
!> \author S. Reine
!> \date 2013-01-22
!> \param F The Fock/KS AO-matrix
!> \param D The AO density-matrix
!> \param nbast The number of oribtals
!> \param ndmat The number of density matrices
!> \param dsym True if all density matrices are symmetric
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The default error-unit. If -1 on input it opens default LSDALTON.ERR
!> \param testElectrons Test # of electons for XC quadrature
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_Fock(F,D,nbast,ndmat,dsym,lupri,luerr,testElectrons,h1)
  use precision
  use Matrix_module
  use Matrix_Operations
  use integralinterfaceMod
  use memory_handling
  use lslib_state
  use dal_interface, only: di_get_fock_lsdalton
IMPLICIT NONE
Integer,intent(in)      :: nbast,ndmat,lupri,luerr
Real(realk),intent(out) :: F(nbast,nbast,ndmat)
Real(realk),intent(in)  :: D(nbast,nbast,ndmat)
Real(realk),optional    :: h1(nbast,nbast)
Logical,intent(IN)      :: Dsym
Logical,optional        :: testElectrons
!
TYPE(MATRIX),target  :: Fmat(ndmat),Dmat(ndmat)
Integer              :: mtype_save, idmat
Real(realk)          :: Etotal(ndmat)
TYPE(MATRIX)         :: h1mat

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_Fock error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_Fock. Basis-function mismatch',lupri)

!this will most likely fail for ScaLapack
mtype_save = matrix_type
matrix_type = mtype_dense
DO idmat=1,ndmat
  CALL mat_init(Fmat(idmat),nbast,nbast)
  call mat_zero(Fmat(idmat))
  CALL mat_init(Dmat(idmat),nbast,nbast)
  call dcopy(nbast*nbast,D(1,1,idmat),1,Dmat(idmat)%elms,1)
ENDDO

CALL mat_init(h1mat,nbast,nbast)
IF (PRESENT(h1)) THEN
  call dcopy(nbast*nbast,h1,1,h1mat%elms,1)
ELSE
  call II_get_h1(lupri,luerr,ls%setting,h1mat) 
ENDIF
IF (PRESENT(testElectrons)) ls%setting%scheme%dft%testNelectrons = testElectrons

call di_get_fock_LSDALTON(Dmat,h1mat,Fmat,ndmat,Etotal,lupri,luerr,ls)
CALL mat_free(h1mat)

DO idmat=1,ndmat
  call dcopy(nbast*nbast,Fmat(idmat)%elms,1,F(1,1,idmat),1)
  CALL mat_free(Fmat(idmat))
  CALL mat_free(Dmat(idmat))
ENDDO

matrix_type = mtype_save

END SUBROUTINE LSlib_get_Fock

!> \brief Calculates the Coulomb AO-matrix
!> \author S. Reine
!> \date 2013-01-22
!> \param J The Coulomb AO-matrix
!> \param D The AO density-matrix
!> \param nbast The number of oribtals
!> \param ndmat The number of density matrices
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_Coulomb(J,D,nbast,ndmat,lupri,luerr)
  use precision
  use Matrix_module
  use Matrix_Operations
  use integralinterfaceMod
  use memory_handling
  use lslib_state
IMPLICIT NONE
Integer,intent(in)      :: nbast,ndmat,lupri,luerr
Real(realk),intent(out) :: J(nbast,nbast,ndmat)
Real(realk),intent(in)  :: D(nbast,nbast,ndmat)
!
TYPE(MATRIX),target :: Jmat(ndmat),Dmat(ndmat)
Integer             :: mtype_save, idmat

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_Coulomb error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_Coulomb. Basis-function mismatch',lupri)

!this will most likely fail for ScaLapack
mtype_save = matrix_type
matrix_type = mtype_dense
DO idmat=1,ndmat
  CALL mat_init(Jmat(idmat),nbast,nbast)
  CALL mat_init(Dmat(idmat),nbast,nbast)
  call dcopy(nbast*nbast,D(1,1,idmat),1,Dmat(idmat)%elms,1)
ENDDO

CALL II_get_coulomb_mat(lupri,luerr,ls%setting,Dmat,Jmat,ndmat)

DO idmat=1,ndmat
  call dcopy(nbast*nbast,Jmat(idmat)%elms,1,J(1,1,idmat),1)
  CALL mat_free(Jmat(idmat))
  CALL mat_free(Dmat(idmat))
ENDDO

matrix_type = mtype_save

END SUBROUTINE LSlib_get_Coulomb

!> \brief Calculates the Exchange AO-matrix
!> \author S. Reine
!> \date 2013-01-22
!> \param K The Exchange AO-matrix
!> \param D The AO density-matrix
!> \param nbast The number of oribtals
!> \param ndmat The number of density matrices
!> \param dsym True if all density matrices are symmetric
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_Exchange(K,D,nbast,ndmat,dsym,lupri,luerr,testElectrons)
  use precision
!  use TYPEDEF  
!  use TYPEDEFTYPE  
  use Matrix_module
  use Matrix_Operations
!  use ls_Integral_Interface
!  use daltonInfo
  use integralinterfaceMod
  use memory_handling
  use lslib_state
IMPLICIT NONE
Integer,intent(in)      :: nbast,ndmat,lupri,luerr
Real(realk),intent(out) :: K(nbast,nbast,ndmat)
Real(realk),intent(in)  :: D(nbast,nbast,ndmat)
Logical,intent(IN)      :: Dsym
Logical,optional        :: testElectrons
!
TYPE(MATRIX),target  :: Kmat(ndmat),Dmat(ndmat),dXCmat
Integer              :: mtype_save, idmat
real(realk)          :: EdXC,Econt(5)

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_Exchange error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_Exchange. Basis-function mismatch',lupri)

!this will most likely fail for ScaLapack
mtype_save = matrix_type
matrix_type = mtype_dense
DO idmat=1,ndmat
  CALL mat_init(Kmat(idmat),nbast,nbast)
  call mat_zero(Kmat(idmat))
  CALL mat_init(Dmat(idmat),nbast,nbast)
  call dcopy(nbast*nbast,D(1,1,idmat),1,Dmat(idmat)%elms,1)
ENDDO

IF (ls%input%dalton%ADMM_EXCHANGE) THEN
  IF (present(testElectrons)) ls%setting%scheme%dft%testNelectrons = testElectrons
  call mat_init(dXCmat,nbast,nbast)
  DO idmat=1,ndmat
    CALL II_get_admm_exchange_mat(LUPRI,LUERR,ls%SETTING,ls%optlevel,Dmat(idmat),Kmat(idmat),dXCmat,1,EdXC,Dsym,Econt,.FALSE.)
    call mat_daxpy(1.E0_realk,dXCmat,Kmat(idmat))
  ENDDO
  CALL mat_free(dXCmat)
ELSE
  CALL II_get_exchange_mat(lupri,luerr,ls%setting,Dmat,ndmat,Dsym,Kmat)
ENDIF

DO idmat=1,ndmat
  call dcopy(nbast*nbast,Kmat(idmat)%elms,1,K(1,1,idmat),1)
  CALL mat_free(Kmat(idmat))
  CALL mat_free(Dmat(idmat))
ENDDO

matrix_type = mtype_save

END SUBROUTINE LSlib_get_Exchange

!> \brief Calculates the exchange-correlation AO-matrix
!> \author S. Reine
!> \date 2013-01-23
!> \param XC The AO XC-matrix
!> \param D The AO density-matrix
!> \param EXC The XC energy
!> \param testElectrons The number of oribtals
!> \param nbast The number of oribtals
!> \param ndmat The number of density matrices
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_XC(XC,D,EXC,nbast,ndmat,lupri,luerr,testElectrons)
  use precision
  use Matrix_module
  use Matrix_Operations
  use memory_handling
  use II_XC_interfaceModule
  use lslib_state
IMPLICIT NONE
Integer,intent(in)      :: nbast,ndmat,lupri,luerr
Real(realk),intent(out) :: XC(nbast,nbast,ndmat),EXC(ndmat)
Real(realk),intent(in)  :: D(nbast,nbast,ndmat)
logical,optional        :: testElectrons
!
TYPE(MATRIX),target  :: XCmat(ndmat),Dmat(ndmat)
Integer              :: mtype_save, idmat

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_XC error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_XC. Basis-function mismatch',lupri)

!this will most likely fail for ScaLapack
mtype_save = matrix_type
matrix_type = mtype_dense
DO idmat=1,ndmat
  CALL mat_init(XCmat(idmat),nbast,nbast)
  call mat_zero(XCmat(idmat))
  CALL mat_init(Dmat(idmat),nbast,nbast)
  call dcopy(nbast*nbast,D(1,1,idmat),1,Dmat(idmat)%elms,1)
ENDDO

IF (present(testElectrons)) ls%setting%scheme%dft%testNelectrons = testElectrons
IF (ls%setting%do_dft) THEN
  CALL II_get_xc_Fock_mat(lupri,luerr,ls%setting,nbast,Dmat,XCmat,EXC,ndmat)
ENDIF

DO idmat=1,ndmat
  call dcopy(nbast*nbast,XCmat(idmat)%elms,1,XC(1,1,idmat),1)
  CALL mat_free(XCmat(idmat))
  CALL mat_free(Dmat(idmat))
ENDDO

matrix_type = mtype_save

END SUBROUTINE LSlib_get_XC


!> \brief Calculates the one-electron AO-matrix
!> \author S. Reine
!> \date 2010-02-26
!> \param h1 The one-electron AO-matrix
!> \param nbast The number of oribtals
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_h1(h1,nbast,lupri,luerr)
  use precision
  use TYPEDEF  
  use TYPEDEFTYPE  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
  use integralinterfaceMod
  use memory_handling
  use lslib_state
IMPLICIT NONE
Integer,intent(in)      :: nbast,lupri,luerr
Real(realk),intent(out) :: h1(nbast,nbast)
!
TYPE(MATRIX),target :: h
Integer             :: mtype_save

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_h1 error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_h1. Basis-function mismatch',lupri)

!this will most likely fail for MPI ask Thomas (or call call mat_select_type)
mtype_save = matrix_type
matrix_type = mtype_dense
CALL mat_init(h,nbast,nbast)

CALL II_get_h1(lupri,luerr,ls%setting,h)

call dcopy(nbast*nbast,h%elms,1,h1,1)

CALL mat_free(h)
matrix_type = mtype_save

END SUBROUTINE LSlib_get_h1

!> \brief Calculates the 4-center 2-electron AO electron-repulsion integrals (ERI)
!> \author S. Reine
!> \date 2013-01-20
!> \param eri The 4-center 2-electron AO integrals
!> \param nbast The number of oribtals
!> \param dirac Specifies Dirac or Mulliken notation
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_4center_eri(eri,nbast,dirac,lupri,luerr)
  use precision
  use TYPEDEF  
  use TYPEDEFTYPE  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
  use integralinterfaceMod
  use memory_handling
  use lslib_state
IMPLICIT NONE
Integer,intent(in)      :: nbast,lupri,luerr
Real(realk),intent(out) :: eri(nbast,nbast,nbast,nbast,1)
Logical,intent(IN)      :: dirac
character :: intspec(5)
intspec(1) = 'R'
intspec(2) = 'R'
intspec(3) = 'R'
intspec(4) = 'R'
intspec(5) = 'C'

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_4center_eri error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_4center_eri. Basis-function mismatch',lupri)

CALL II_get_4center_eri(lupri,luerr,ls%setting,eri,&
     & nbast,nbast,nbast,nbast,intspec,dirac)

END SUBROUTINE LSlib_get_4center_eri

!> \brief Calculates the differentiated 4-center 2-electron AO electron-repulsion integrals (ERI)
!> \author S. Reine
!> \date 2013-01-20
!> \param eri The 4-center 2-electron AO integrals
!> \param nbast The number of oribtals
!> \param geoOrder The the geometry differential order (0 regular integrals, 1 first order integral derivatives, etc.)
!> \param nGeoComp The number of derivative components (1 for 0th order, 3*nAtoms for 1st order, (3*nAtoms)**2 for 2nd, etc.)
!> \param dirac Specifies Dirac or Mulliken notation
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_4center_eri_geoderiv(eri,nbast,geoOrder,nGeoComp,dirac,lupri,luerr)
  use precision
  use TYPEDEF  
  use TYPEDEFTYPE  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
  use integralinterfaceMod
  use memory_handling
  use lslib_state
IMPLICIT NONE
Integer,intent(in)      :: nbast,lupri,luerr
Real(realk),intent(out) :: eri(nbast,nbast,nbast,nbast,nGeoComp)
Integer,intent(IN)      :: geoOrder,nGeoComp
Logical,intent(IN)      :: dirac

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_4center_eri_geoderiv error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_4center_eri_geoderiv. Basis-function mismatch',lupri)

CALL II_get_hodi_eri_4center(lupri,luerr,ls%setting,eri,nbast,nbast,nbast,nbast,nGeoComp,dirac=dirac,geoderiv=geoOrder)

END SUBROUTINE LSlib_get_4center_eri_geoderiv

!> \brief Calculates the differentiated 1-electron AO integrals
!> \author S. Reine
!> \date 2013-01-20
!> \param 1el The 2-center 1-electron AO integrals (for nucel the nucei are implictly summed)
!> \param nbast The number of oribtals
!> \param geoOrder The the geometry differential order (0 regular integrals, 1 first order integral derivatives, etc.)
!> \param nGeoComp The number of derivative components (1 for 0th order, 3*nAtoms for 1st order, (3*nAtoms)**2 for 2nd, etc.)
!> \param dirac Specifies Dirac or Mulliken notation
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_1el_geoderiv(oneEl,oneElType,nbast,nAtoms,geoOrder,nGeoComp,lupri,luerr)
  use precision
  use TYPEDEF  
  use TYPEDEFTYPE  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
  use integralinterfaceMod
  use memory_handling
  use lslib_state
IMPLICIT NONE
Integer,intent(in)      :: nbast,nAtoms,lupri,luerr
Real(realk),intent(out) :: oneEl(nbast,nbast,1,1,nGeoComp)
Integer,intent(IN)      :: geoOrder,nGeoComp
Character(*),intent(IN) :: oneElType

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_1el_geoderiv error: LSlib not initialized',lupri)
IF (nbasis.NE.nbast) CALL lsQUIT('Error in LSlib_get_1el_geoderiv. Basis-function mismatch',lupri)

CALL II_get_hodi_1el(lupri,luerr,ls%setting,oneEl,oneElType,nbast,nbast,nGeoComp,geoderiv=geoOrder)

END SUBROUTINE LSlib_get_1el_geoderiv

!> \brief Calculates the gradient of the nuclear potential
!> \author S. Reine
!> \date 2010-02-26
!> \param nucGrad The nuclear potential gradient
!> \param nAtoms The number of atoms
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_nn_gradient(nucGrad,nAtoms,lupri,luerr)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use integralinterfaceMod
use memory_handling
use lslib_state
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr
Real(realk),intent(OUT) :: nucGrad(3,nAtoms)

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_nn_gradient error: LSlib not initialized',lupri)

call II_get_nn_gradient(nucGrad,ls%setting,lupri,luerr)

END SUBROUTINE LSlib_get_nn_gradient

!> \brief Calculates the one-electron gradient contribution
!> \author S. Reine
!> \date 2010-03-22
!> \param oneGrad Nuclear-electronic-attraction gradient
!> \param D The density matrix
!> \param nbast The number of orbitals/basis functions
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_oneElectron_gradient(oneGrad,D,nbast,nAtoms,lupri,luerr)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use matrix_module
use matrix_operations
use integralinterfaceMod
use memory_handling
use lslib_state
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast
Integer,parameter       :: ndmat = 1
Real(realk),intent(OUT) :: oneGrad(3,nAtoms)
Real(realk),intent(IN)  :: D(nbast,nbast,ndmat)
!
type(matrixp)       :: Dmat(ndmat)
type(matrix),target :: Dtarget(ndmat)
integer             :: idmat

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_oneElectron_gradient error: LSlib not initialized',lupri)
IF (ndmat.lt. 1) CALL lsQUIT('ndmat<1 in LSlib_get_oneElectron_gradient',lupri)
IF (nAtoms.lt.ls%setting%molecule(1)%p%nAtoms) CALL lsQUIT('nAtoms inconsistency in LSlib_get_oneElectron_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_oneElectron_gradient',lupri)

!Copy density-matrices from D to Dmat (using Dtarget)
DO idmat=1,ndmat
  Dmat(idmat)%p => Dtarget(idmat)
  call mat_init(Dmat(idmat)%p,nbast,nbast)
  call mat_set_from_full(D(:,:,idmat),1E0_realk,Dmat(idmat)%p)
ENDDO

!Calculate the nuclear-electron attraction gradient
CALL II_get_oneElectron_gradient(oneGrad,Dmat,ndmat,ls%setting,lupri,luerr)

!Free density-matrices Dmat
DO idmat=1,ndmat
  call mat_free(Dmat(idmat)%p)
ENDDO

END SUBROUTINE LSlib_get_oneElectron_gradient

!> \brief Calculates the nuclear-electronic-attraction gradient
!> \author S. Reine
!> \date 2010-03-22
!> \param neGrad Nuclear-electronic-attraction gradient
!> \param D The density matrix
!> \param nbast The number of orbitals/basis functions
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_ne_gradient(neGrad,D,nbast,nAtoms,lupri,luerr)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use matrix_module
use matrix_operations
use integralinterfaceMod
use memory_handling
use lslib_state
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast
Integer,parameter       :: ndmat = 1
Real(realk),intent(OUT) :: neGrad(3,nAtoms)
Real(realk),intent(IN)  :: D(nbast,nbast,ndmat)
!
type(matrixp)       :: Dmat(ndmat)
type(matrix),target :: Dtarget(ndmat)
integer             :: idmat

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_ne_gradient error: LSlib not initialized',lupri)
IF (ndmat.lt. 1) CALL lsQUIT('ndmat<1 in LSlib_get_ne_gradient',lupri)
IF (nbast.lt. 1) CALL lsQUIT('nbast<1 in LSlib_get_ne_gradient',lupri)
IF (nAtoms.lt. 1) CALL lsQUIT('nAtoms<1 in LSlib_get_ne_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_ne_gradient',lupri)

!Copy density-matrices from D to Dmat (using Dtarget)
DO idmat=1,ndmat
  Dmat(idmat)%p => Dtarget(idmat)
  call mat_init(Dmat(idmat)%p,nbast,nbast)
  call mat_set_from_full(D(:,:,idmat),1E0_realk,Dmat(idmat)%p)
ENDDO

!Calculate the nuclear-electron attraction gradient
CALL II_get_ne_gradient(neGrad,Dmat,ndmat,ls%setting,lupri,luerr)

!Free density-matrices Dmat
DO idmat=1,ndmat
  call mat_free(Dmat(idmat)%p)
ENDDO

END SUBROUTINE LSlib_get_ne_gradient

!> \brief Calculates the electron-electron repulsion gradient
!> \author S. Reine
!> \date 2010-04-06
!> \param eeGrad The exchange gradient contribution
!> \param DLHS The density matrix of the first electron (or left-hand-side)
!> \param DRHS The density matrix of the second electron (or right-hand-side)
!> \param nbast The number of orbitals/basis functions
!> \param ndlhs The number of density LHS matrices
!> \param ndrhs The number of density RHS matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_twoElectron_gradient(eeGrad,DLHS,DRHS,nbast,ndlhs,ndrhs,nAtoms,lupri,luerr,testElectrons)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use matrix_module
use matrix_operations
use integralinterfaceMod
use memory_handling
use io
use lslib_state
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndlhs,ndrhs
Real(realk),intent(OUT) :: eeGrad(3,nAtoms)
Real(realk),intent(IN)  :: DLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)  :: DRHS(nbast,nbast,ndrhs)
Logical,optional        :: testElectrons
!
type(matrixp)       :: DmatLHS(ndlhs),DmatRHS(ndrhs)
type(matrix),target :: DtargetLHS(ndlhs),DtargetRHS(ndrhs)
integer             :: idmat
logical             :: coeff
Character(80)       :: Filename

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_twoElectron_gradient error: LSlib not initialized',lupri)
IF (ndlhs.lt. 1) CALL lsQUIT('ndlhs<1 in LSlib_get_twoElectron_gradient',lupri)
IF (ndrhs.lt. 1) CALL lsQUIT('ndrhs<1 in LSlib_get_twoElectron_gradient',lupri)
IF (nbast.lt. 1) CALL lsQUIT('nbast<1 in LSlib_get_twoElectron_gradient',lupri)
IF (nAtoms.lt. 1) CALL lsQUIT('nAtoms<1 in LSlib_get_twoElectron_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_twoElectron_gradient',lupri)

!Copy density-matrices from DLHS to DmatLHS (using DtargetLHS)
DO idmat=1,ndlhs
  DmatLHS(idmat)%p => DtargetLHS(idmat)
  call mat_init(DmatLHS(idmat)%p,nbast,nbast)
  call mat_set_from_full(DLHS(:,:,idmat),1E0_realk,DmatLHS(idmat)%p)
ENDDO
!Copy density-matrices from DRHS to DmatRHS (using DtargetRHS)
DO idmat=1,ndrhs
  DmatRHS(idmat)%p => DtargetRHS(idmat)
  call mat_init(DmatRHS(idmat)%p,nbast,nbast)
  call mat_set_from_full(DRHS(:,:,idmat),1E0_realk,DmatRHS(idmat)%p)
ENDDO

IF (present(testElectrons)) ls%setting%scheme%dft%testNelectrons = testElectrons

CALL II_get_twoElectron_gradient(eeGrad,nAtoms,DmatLHS,DmatRHS,ndlhs,ndrhs,ls%setting,lupri,luerr)

!Free density-matrices DmatLHS and DmatRHS
DO idmat=1,ndlhs
  call mat_free(DmatLHS(idmat)%p)
ENDDO
DO idmat=1,ndrhs
  call mat_free(DmatRHS(idmat)%p)
ENDDO

END SUBROUTINE LSlib_get_twoElectron_gradient

!> \brief Calculates the exchange-correlation contribution to the molecular gradient
!> \author T. Kjaergaard
!> \date 2010-04-26
!> \param exGrad The exchange-correlation gradient contribution
!> \param DMAT The density matrix
!> \param nbast The number of orbitals/basis functions
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_xc_gradient(exGrad,Dmat,nbast,nAtoms,lupri,luerr,testElectrons)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use matrix_module
use matrix_operations
use integralinterfaceMod
use memory_handling
use II_XC_interfaceModule
use lslib_state
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast
Real(realk),intent(OUT) :: exGrad(3,nAtoms)
Real(realk),intent(IN)  :: DMAT(nbast,nbast)
Logical,optional        :: testElectrons
!
type(matrix) :: D
logical      :: doDFT,testNel_save

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_xc_gradient error: LSlib not initialized',lupri)
IF (nbast.lt. 1) CALL lsQUIT('nbast<1 in LSlib_get_xc_gradient',lupri)
IF (nAtoms.lt. 1) CALL lsQUIT('nAtoms<1 in LSlib_get_xc_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_xc_gradient',lupri)

call mat_init(D,nbast,nbast)
call mat_set_from_full(DMAT,1E0_realk,D)

IF (present(testElectrons)) THEN
  testNel_save = ls%setting%scheme%dft%testNelectrons
  ls%setting%scheme%dft%testNelectrons = testElectrons
ENDIF

!Calculate the exchange-correlation gradient
doDFT = LSlibconfig%opt%calctype.EQ.LSlibconfig%opt%dftcalc
IF (doDFT) THEN
  CALL II_get_xc_geoderiv_molgrad(lupri,luerr,ls%setting,nbast,D,exGrad,natoms)
ENDIF

IF (present(testElectrons)) THEN
  ls%setting%scheme%dft%testNelectrons = testNel_save
ENDIF

END SUBROUTINE LSlib_get_xc_gradient

!> \brief Calculates the Coulomb contribution to the gradient
!> \author S. Reine
!> \date 2010-04-06
!> \param coulombGrad The Coulomb gradient contribution
!> \param DLHS The density matrix of the first electron (or left-hand-side)
!> \param DRHS The density matrix of the second electron (or right-hand-side)
!> \param nbast The number of orbitals/basis functions
!> \param ndlhs The number of density LHS matrices
!> \param ndrhs The number of density RHS matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_J_gradient(coulombGrad,DLHS,DRHS,nbast,ndlhs,ndrhs,nAtoms,lupri,luerr)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use matrix_module
use matrix_operations
use integralinterfaceMod
use memory_handling
use io
use lslib_state
implicit none
Integer,intent(IN)        :: nAtoms,lupri,luerr,nbast,ndlhs,ndrhs
Real(realk),intent(INOUT) :: coulombGrad(3,nAtoms)
Real(realk),intent(IN)    :: DLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)    :: DRHS(nbast,nbast,ndrhs)
!
type(matrixp)       :: DmatLHS(ndlhs),DmatRHS(ndrhs)
type(matrix),target :: DtargetLHS(ndlhs),DtargetRHS(ndrhs)
integer             :: idmat
Character(80)       :: Filename
Logical             :: coeff

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_J_gradient error: LSlib not initialized',lupri)
IF (ndlhs.lt. 1) CALL lsQUIT('ndlhs<1 in LSlib_get_J_gradient',lupri)
IF (ndrhs.lt. 1) CALL lsQUIT('ndrhs<1 in LSlib_get_J_gradient',lupri)
IF (nbast.lt. 1) CALL lsQUIT('nbast<1 in LSlib_get_J_gradient',lupri)
IF (nAtoms.lt. 1) CALL lsQUIT('nAtoms<1 in LSlib_get_J_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_J_gradient',lupri)

!Copy density-matrices from DLHS to DmatLHS (using DtargetLHS)
DO idmat=1,ndlhs
  call mat_init(DtargetLHS(idmat),nbast,nbast)
  call mat_set_from_full(DLHS(:,:,idmat),1E0_realk,DtargetLHS(idmat))
  DmatLHS(idmat)%p => DtargetLHS(idmat)
ENDDO
!Copy density-matrices from DRHS to DmatRHS (using DtargetRHS)
DO idmat=1,ndrhs
  call mat_init(DtargetRHS(idmat),nbast,nbast)
  call mat_set_from_full(DRHS(:,:,idmat),1E0_realk,DtargetRHS(idmat))
  DmatRHS(idmat)%p => DtargetRHS(idmat)
ENDDO

!Calculate the nuclear-electron attraction gradient
CALL II_get_J_gradient(coulombGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,ls%setting,lupri,luerr)

!Free density-matrices DmatLHS and DmatRHS
DO idmat=1,ndlhs
  call mat_free(DtargetLHS(idmat))
ENDDO
DO idmat=1,ndrhs
  call mat_free(DtargetRHS(idmat))
ENDDO

END SUBROUTINE LSlib_get_J_gradient

!> \brief Calculates the exchange contribution to the gradient
!> \author S. Reine
!> \date 2010-04-06
!> \param exchangeGrad The exchange gradient contribution
!> \param DLHS The density matrix of the first electron (or left-hand-side)
!> \param DRHS The density matrix of the second electron (or right-hand-side)
!> \param nbast The number of orbitals/basis functions
!> \param ndlhs The number of density LHS matrices
!> \param ndrhs The number of density RHS matrices
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_K_gradient(exchangeGrad,DLHS,DRHS,nbast,ndlhs,ndrhs,nAtoms,lupri,luerr)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use matrix_module
use matrix_operations
use integralinterfaceMod
use memory_handling
use lslib_state
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndlhs,ndrhs
Real(realk),intent(OUT) :: exchangeGrad(3,nAtoms)
Real(realk),intent(IN)  :: DLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)  :: DRHS(nbast,nbast,ndrhs)
!
type(matrixp)       :: DmatLHS(ndlhs),DmatRHS(ndrhs)
type(matrix),target :: DtargetLHS(ndlhs),DtargetRHS(ndrhs)
integer             :: idmat

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_K_gradient error: LSlib not initialized',lupri)
IF (nbast.lt. 1) CALL lsQUIT('nbast<1 in LSlib_get_K_gradient',lupri)
IF (nAtoms.lt. 1) CALL lsQUIT('nAtoms<1 in LSlib_get_K_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_K_gradient',lupri)

!Copy density-matrices from DLHS to DmatLHS (using DtargetLHS)
DO idmat=1,ndlhs
  DmatLHS(idmat)%p => DtargetLHS(idmat)
  call mat_init(DmatLHS(idmat)%p,nbast,nbast)
  call mat_set_from_full(DLHS(:,:,idmat),1E0_realk,DmatLHS(idmat)%p)
ENDDO
!Copy density-matrices from DRHS to DmatRHS (using DtargetRHS)
DO idmat=1,ndrhs
  DmatRHS(idmat)%p => DtargetRHS(idmat)
  call mat_init(DmatRHS(idmat)%p,nbast,nbast)
  call mat_set_from_full(DRHS(:,:,idmat),1E0_realk,DmatRHS(idmat)%p)
ENDDO

!Calculate the exchange contribution gradient
CALL II_get_K_gradient(exchangeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,ls%setting,lupri,luerr)

!Free density-matrices DmatLHS and DmatRHS
DO idmat=1,ndlhs
  call mat_free(DmatLHS(idmat)%p)
ENDDO
DO idmat=1,ndrhs
  call mat_free(DmatRHS(idmat)%p)
ENDDO

END SUBROUTINE LSlib_get_K_gradient

!> \brief Calculates the kinetic energy gradient
!> \author S. Reine
!> \date 2010-03-22
!> \param kinGrad Kinetic energy gradient
!> \param D The density matrix
!> \param nbast The number of orbitals/basis functions
!> \param nAtoms The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error output unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_kinetic_gradient(kinGrad,D,nbast,nAtoms,lupri,luerr)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use matrix_module
use matrix_operations
use integralinterfaceMod
use memory_handling
use lslib_state
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast
Real(realk),intent(OUT) :: kinGrad(3,nAtoms)
Integer,parameter       :: ndmat = 1
Real(realk),intent(IN)  :: D(nbast,nbast,ndmat)
!
type(matrixp)       :: Dmat(ndmat)
type(matrix),target :: Dtarget(ndmat)
integer             :: idmat

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_kinetic_gradient error: LSlib not initialized',lupri)
IF (ndmat.lt. 1) CALL lsQUIT('ndmat<1 in LSlib_get_kinetic_gradient',lupri)
IF (nbast.lt. 1) CALL lsQUIT('nbast<1 in LSlib_get_kinetic_gradient',lupri)
IF (nAtoms.lt. 1) CALL lsQUIT('nAtoms<1 in LSlib_get_kinetic_gradient',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_kinetic_gradient',lupri)

!Copy density-matrices from D to Dmat (using Dtarget)
DO idmat=1,ndmat
  Dmat(idmat)%p => Dtarget(idmat)
  call mat_init(Dmat(idmat)%p,nbast,nbast)
  call mat_set_from_full(D(:,:,idmat),1E0_realk,Dmat(idmat)%p)
ENDDO

!Calculate the kinetic energy gradient
CALL II_get_kinetic_gradient(kinGrad,Dmat,ndmat,ls%setting,lupri,luerr)

!Free density-matrices Dmat
DO idmat=1,ndmat
  call mat_free(Dmat(idmat)%p)
ENDDO

END SUBROUTINE LSlib_get_kinetic_gradient

!> \brief Calculates the reorthonormalization gradient term
!> \author S. Reine
!> \date 2010-03-22
!> \param reOrtho The reorthonormalization gradient term
!> \param DFD The matrix product DFD = - D F D
!> \param nbast The number of orbital (basis functions)
!> \param ndmat The number of density matrices
!> \param nAtom The number of atoms
!> \param lupri Default output unit
!> \param luerr Default error print unit
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_get_reorthoNormalization(reOrtho,DFD,nbast,ndmat,nAtoms,lupri,luerr)
use TYPEDEF
use TYPEDEFTYPE
use daltonInfo
use matrix_operations
use matrix_module
use integralinterfaceMod
use memory_handling
use lslib_state
implicit none
Integer,intent(IN)      :: nAtoms,lupri,luerr,nbast,ndmat
Real(realk),intent(OUT) :: reOrtho(3,nAtoms)
Real(realk),intent(IN)  :: DFD(nbast,nbast,ndmat)
!
type(matrixp)       :: DFDmat(ndmat)
type(matrix),target :: DFDtarget(ndmat)
integer             :: idmat

IF (.NOT.state_set) CALL LSQUIT('LSlib_get_reorthoNormalization error: LSlib not initialized',lupri)
IF (ndmat.lt. 1) CALL lsQUIT('ndmat<1 in LSlib_get_reorthoNormalization',lupri)
IF (nbast.lt. 1) CALL lsQUIT('nbast<1 in LSlib_get_reorthoNormalization',lupri)
IF (nAtoms.lt. 1) CALL lsQUIT('nAtoms<1 in LSlib_get_reorthoNormalization',lupri)
IF (nbasis.ne.nbast) CALL lsQUIT('nbasis .NE. nbast in LSlib_get_reorthoNormalization',lupri)

!Copy density-matrices from D to Dmat (using Dtarget)
DO idmat=1,ndmat
  call mat_init(DFDtarget(idmat),nbast,nbast)
  call mat_set_from_full(DFD(:,:,idmat),1E0_realk,DFDtarget(idmat))
  DFDmat(idmat)%p => DFDtarget(idmat)
ENDDO

!Calculate the nuclear-electron attraction gradient
CALL II_get_reorthoNormalization(reOrtho,DFDmat,ndmat,ls%setting,lupri,luerr)

!Free density-matrices Dmat
DO idmat=1,ndmat
  call mat_free(DFDtarget(idmat))
ENDDO

END SUBROUTINE LSlib_get_reorthoNormalization


!> \brief Calculates the exchange matrix and nothing else. Used for performance testing
!> \author T. Kjaergaard
!> \date 2010-03-17
! this routine can be call directly from lsdalton_wrapper
SUBROUTINE LSlib_build_exchange()
  use configurationType, only: configitem
  use configuration, only: config_set_default_config, config_read_input, config_shutdown, config_free, set_final_config_and_print
  use TYPEDEF  
  use TYPEDEFTYPE  
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
  use integralinterfaceMod
  use memory_handling
  use files
  use lstiming
IMPLICIT NONE
Integer      :: lupri,luerr
!
integer      :: nbast,restart_lun
TYPE(MATRIX) :: K,D
type(lsitem) :: ls
logical :: Dsym,dens_exsist,gcbasis,doDFT,OnMaster
real(realk) :: t1,t2
type(configItem)    :: config
OnMaster=.TRUE.
LUPRI=-1
LUERR=-1
CALL LSOPEN(LUPRI,'DALTON.OUT','NEW','FORMATTED')
CALL LSOPEN(LUERR,'DALTON.ERR','UNKNOWN','FORMATTED')

call LSTIMER('START',t1,t2,LUPRI)

call config_set_default_config(config)
config%opt%cfg_start_guess = 'ATOMS'
call config_read_input(config,lupri,luerr)
ls%input%dalton = config%integral
doDFT = config%opt%calctype.EQ.config%opt%dftcalc
call ls_init(ls,lupri,luerr,nbast,config%integral,doDFT,.false.,.true.)
call set_final_config_and_print(lupri,config,ls,nbast)

!get inital density 
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
else
   call lsquit('Calculation using LSlib_build_exchange requires a dens.restart!',lupri)
endif

Dsym = .TRUE. !symmetric Density matrix
call mat_init(K,nbast,nbast)
call mat_zero(K)
call II_get_exchange_mat(LUPRI,LUERR,ls%SETTING,D,1,Dsym,K)
WRITE(lupri,*)'The Exchange energy:',mat_dotproduct(D,K)
CALL mat_free(D)
CALL mat_free(K)
call ls_free(ls)
call stats_mem(lupri)
call LSTIMER('LSlib_build_exc',t1,t2,LUPRI)
call config_shutdown(config)
call config_free(config)

END SUBROUTINE LSlib_build_exchange

!> \brief Debug the LSlib routines in this file
!> \author T. Kjaergaard
!> \date 2012-08
!> \param Overlap The one-electron AO-matrix
!> \param nbast The number of oribtals
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE LSlib_debug(lupri,luerr,setting,nbast)
  use precision
  use typedeftype
!  use typedef
  use matrix_util
  use matrix_module
  use matrix_operations
  use memory_handling
  use integralinterfaceMod
IMPLICIT NONE
Integer,intent(in)      :: lupri,luerr,nbast
type(lssetting)   :: setting
!
TYPE(MATRIX)      :: OverlapMat,S
type(matrix)      :: DIPLENmat,DIPLEN(3)
real(realk) :: THR
integer :: nAtoms,I
real(realk),pointer :: Coord(:,:),Charge(:)
THR = 1E-15_realk
nAtoms = setting%molecule(1)%p%nAtoms
call mem_alloc(Coord,3,nAtoms)
call mem_alloc(Charge,nAtoms)
DO I = 1,nAtoms
   Coord(1,I) = setting%molecule(1)%p%ATOM(I)%CENTER(1)
   Coord(2,I) = setting%molecule(1)%p%ATOM(I)%CENTER(2)
   Coord(3,I) = setting%molecule(1)%p%ATOM(I)%CENTER(3)
   Charge(I) = setting%molecule(1)%p%ATOM(I)%CHARGE
ENDDO

call mat_init(S,nbast,nbast)
CALL II_get_overlap(lupri,luerr,setting,S)
call mat_init(OverlapMat,nbast,nbast)
call lsdalton_get_Overlap(Overlapmat%elms,nbast,nAtoms,Coord,Charge,'ccD',lupri)
call VerifyMatrices(S,OverlapMat,'LSlib_debug: lsdalton_get_overlap',THR,lupri)
call mat_free(S)
call mat_free(OverlapMat)

call mat_init(DIPLEN(1),nbast,nbast)
call mat_init(DIPLEN(2),nbast,nbast)
call mat_init(DIPLEN(3),nbast,nbast)
call II_get_prop(LUPRI,LUERR,SETTING,DIPLEN,3,'DIPLEN ')

call mat_init(DIPLENMat,nbast,nbast)
call lsdalton_get_XDIPLEN(DIPLENmat%elms,nbast,nAtoms,Coord,Charge,'ccD',lupri)
call VerifyMatrices(DIPLEN(1),DIPLENmat,'LSlib_debug: lsdalton_get_XDIPLEN',THR,lupri)
call lsdalton_get_YDIPLEN(DIPLENmat%elms,nbast,nAtoms,Coord,Charge,'ccD',lupri)
call VerifyMatrices(DIPLEN(2),DIPLENmat,'LSlib_debug: lsdalton_get_YDIPLEN',THR,lupri)
call lsdalton_get_ZDIPLEN(DIPLENmat%elms,nbast,nAtoms,Coord,Charge,'ccD',lupri)
call VerifyMatrices(DIPLEN(3),DIPLENmat,'LSlib_debug: lsdalton_get_ZDIPLEN',THR,lupri)

call mat_free(DIPLENMat)
call mat_free(DIPLEN(1))
call mat_free(DIPLEN(2))
call mat_free(DIPLEN(3))
call mem_dealloc(Coord)
call mem_dealloc(Charge)
WRITE(lupri,'(A)')'LSlib Test Successful'

END SUBROUTINE LSlib_debug


!> \brief Calculates the AO-Overlap matrix
!> \author T. Kjaergaard
!> \date 2012-08
!> \param Overlap The one-electron AO-matrix
!> \param nbast The number of oribtals
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE lsdalton_get_Overlap(Overlapmat,nbast,nAtoms,Coord,Charge,BasisString,lupri)
  use precision
  use Matrix_module
  use Matrix_Operations
  use typedeftype
  use typedef
  use daltonInfo
  use integralinterfaceMod
IMPLICIT NONE
character(len=3) :: BasisString
Integer,intent(in)      :: nbast,lupri,nAtoms
Real(realk),intent(inout) :: OverlapMat(nbast,nbast)
Real(realk),intent(in) :: Charge(nAtoms),Coord(3,nAtoms)
!
TYPE(DALTONINPUT) :: input
TYPE(MATRIX)      :: OverlapMat1
Integer           :: mtype_save
type(lssetting)   :: setting
call build_setting_from_scratch(input,setting,nbast,nAtoms,Coord,Charge,BasisString,lupri)
!CALL PRINT_MOLECULE_AND_BASIS(LUPRI,input%MOLECULE,input%BASIS%REGULAR)

mtype_save = matrix_type
!matrix_type = mtype_dense
CALL mat_select_type(mtype_dense,lupri)

CALL mat_init(OverlapMat1,nbast,nbast)
CALL II_get_Overlap(lupri,lupri,setting,OverlapMat1)
call dcopy(nbast*nbast,OverlapMat1%elms,1,OverlapMat,1)
CALL mat_free(OverlapMat1)

CALL dalton_finalize(input,lupri,lupri)
call typedef_free_setting(SETTING)
CALL mat_select_type(mtype_save,lupri)
!matrix_type = mtype_save
END SUBROUTINE Lsdalton_get_Overlap

!> \brief Calculates the AO-Overlap matrix
!> \author T. Kjaergaard
!> \date 2012-08
!> \param Overlap The one-electron AO-matrix
!> \param nbast The number of oribtals
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE lsdalton_get_XDIPLEN(XDIPLENmat,nbast,nAtoms,Coord,Charge,BasisString,lupri)
  use precision
  use Matrix_module
  use Matrix_Operations
  use typedef
  use typedeftype
  use daltonInfo
  use integralinterfaceMod
IMPLICIT NONE
character(len=3) :: BasisString
Integer,intent(in)      :: nbast,lupri,nAtoms
Real(realk),intent(inout) :: XDIPLENMat(nbast,nbast)
Real(realk),intent(in) :: Charge(nAtoms),Coord(3,nAtoms)
!
TYPE(DALTONINPUT) :: input
TYPE(MATRIX)      :: Mat1(3)
Integer           :: mtype_save
type(lssetting)   :: setting
call build_setting_from_scratch(input,setting,nbast,nAtoms,Coord,Charge,BasisString,lupri)
!CALL PRINT_MOLECULE_AND_BASIS(LUPRI,input%MOLECULE,input%BASIS%REGULAR)

mtype_save = matrix_type
!matrix_type = mtype_dense
CALL mat_select_type(mtype_dense,lupri)
CALL mat_init(Mat1(1),nbast,nbast)
CALL mat_init(Mat1(2),nbast,nbast)
CALL mat_init(Mat1(3),nbast,nbast)

call II_get_prop(LUPRI,LUPRI,SETTING,mat1,3,'DIPLEN ')

call dcopy(nbast*nbast,Mat1(1)%elms,1,XdiplenMat,1)
CALL mat_free(Mat1(1))
CALL mat_free(Mat1(2))
CALL mat_free(Mat1(3))

CALL dalton_finalize(input,lupri,lupri)
call typedef_free_setting(SETTING)
!matrix_type = mtype_save
CALL mat_select_type(mtype_save,lupri)
END SUBROUTINE Lsdalton_get_XDIPLEN

!> \brief Calculates the AO-Overlap matrix
!> \author T. Kjaergaard
!> \date 2012-08
!> \param Overlap The one-electron AO-matrix
!> \param nbast The number of oribtals
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE lsdalton_get_YDIPLEN(YDIPLENmat,nbast,nAtoms,Coord,Charge,BasisString,lupri)
  use precision
  use Matrix_module
  use Matrix_Operations
  use daltonInfo
  use typedef
  use typedeftype
  use integralinterfaceMod
IMPLICIT NONE
character(len=3) :: BasisString
Integer,intent(in)      :: nbast,lupri,nAtoms
Real(realk),intent(inout) :: YDIPLENMat(nbast,nbast)
Real(realk),intent(in) :: Charge(nAtoms),Coord(3,nAtoms)
!
TYPE(DALTONINPUT) :: input
TYPE(MATRIX)      :: Mat1(3)
Integer           :: mtype_save
type(lssetting)   :: setting
call build_setting_from_scratch(input,setting,nbast,nAtoms,Coord,Charge,BasisString,lupri)
!CALL PRINT_MOLECULE_AND_BASIS(LUPRI,input%MOLECULE,input%BASIS%REGULAR)

mtype_save = matrix_type
!matrix_type = mtype_dense
CALL mat_select_type(mtype_dense,lupri)
CALL mat_init(Mat1(1),nbast,nbast)
CALL mat_init(Mat1(2),nbast,nbast)
CALL mat_init(Mat1(3),nbast,nbast)

call II_get_prop(LUPRI,LUPRI,SETTING,mat1,3,'DIPLEN ')

call dcopy(nbast*nbast,Mat1(2)%elms,1,YdiplenMat,1)
CALL mat_free(Mat1(1))
CALL mat_free(Mat1(2))
CALL mat_free(Mat1(3))

CALL dalton_finalize(input,lupri,lupri)
call typedef_free_setting(SETTING)
!matrix_type = mtype_save
CALL mat_select_type(mtype_save,lupri)
END SUBROUTINE Lsdalton_get_YDIPLEN

!> \brief Calculates the AO-Overlap matrix
!> \author T. Kjaergaard
!> \date 2012-08
!> \param Overlap The one-electron AO-matrix
!> \param nbast The number of oribtals
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE lsdalton_get_ZDIPLEN(ZDIPLENmat,nbast,nAtoms,Coord,Charge,BasisString,lupri)
  use precision
  use Matrix_module
  use Matrix_Operations
  use typedef
  use typedeftype
  use daltonInfo
  use integralinterfaceMod
IMPLICIT NONE
character(len=3) :: BasisString
Integer,intent(in)      :: nbast,lupri,nAtoms
Real(realk),intent(inout) :: ZDIPLENMat(nbast,nbast)
Real(realk),intent(in) :: Charge(nAtoms),Coord(3,nAtoms)
!
TYPE(DALTONINPUT) :: input
TYPE(MATRIX)      :: Mat1(3)
Integer           :: mtype_save
type(lssetting)   :: setting
call build_setting_from_scratch(input,setting,nbast,nAtoms,Coord,Charge,BasisString,lupri)
!CALL PRINT_MOLECULE_AND_BASIS(LUPRI,input%MOLECULE,input%BASIS%REGULAR)

mtype_save = matrix_type
!matrix_type = mtype_denseC
CALL mat_select_type(mtype_dense,lupri)
CALL mat_init(Mat1(1),nbast,nbast)
CALL mat_init(Mat1(2),nbast,nbast)
CALL mat_init(Mat1(3),nbast,nbast)

call II_get_prop(LUPRI,LUPRI,SETTING,mat1,3,'DIPLEN ')

call dcopy(nbast*nbast,Mat1(3)%elms,1,ZdiplenMat,1)
CALL mat_free(Mat1(1))
CALL mat_free(Mat1(2))
CALL mat_free(Mat1(3))

CALL dalton_finalize(input,lupri,lupri)
call typedef_free_setting(SETTING)
!matrix_type = mtype_save
CALL mat_select_type(mtype_save,lupri)
END SUBROUTINE Lsdalton_get_ZDIPLEN

subroutine build_setting_from_scratch(input,setting,nbast,nAtoms,Coord,Charge,&
     & BasisString,lupri)
  use precision
  use configurationType, only: configitem
  use configuration, only: config_set_default_config, config_read_input, config_shutdown, config_free
  use Matrix_module
  use Matrix_Operations
  use ls_Integral_Interface
  use daltonInfo
  use molecule_type
  use BUILDBASISSET
  use typedef
  use typedeftype
  use integral_type
  use memory_handling
  use io
#ifdef VAR_MPI
  use lsmpi_type
#endif
  implicit none
  TYPE(DALTONINPUT)    :: input
  character(len=3) :: BasisString
  Integer,intent(in)      :: nbast,lupri,nAtoms
  type(lssetting) :: setting
  Real(realk),intent(in) :: Charge(nAtoms),Coord(3,nAtoms)
  !
  TYPE(integralconfig)   :: integral
  character(len=80) :: BasisString2
  TYPE(BASISSETLIBRARYITEM) :: LIBRARY(nBasisBasParam)  
  integer,pointer     :: UNIQUECHARGES(:)
  Integer             :: nbas,I,J,nUCharge
  integer(kind=ls_mpik)::ierr,mynum,nodtot
  call integral_set_default_config(integral)
  input%DALTON = integral
  do I=1,80
     BasisString2(I:I) = ' '
  enddo
  NULLIFY(input%AUXMOLECULE)
  NULLIFY(input%MOLECULE)
  NULLIFY(input%BASIS)
  ALLOCATE(input%AUXMOLECULE)
  ALLOCATE(input%MOLECULE)
  ALLOCATE(input%BASIS)
  input%AUXMOLECULE%nAtoms=0
  input%AUXMOLECULE%nSubSystems=0
  CALL io_init(input%IO)
  call build_Molecule_From_coordList(input%MOLECULE,coord,natoms,charge,lupri)
  SELECT CASE(BasisString)
  CASE('ccd'); BasisString2(1:7) = 'cc-pVDZ'
  CASE('ccD'); BasisString2(1:7) = 'cc-pVDZ'
  CASE('CCD'); BasisString2(1:7) = 'cc-pVDZ'
  CASE('cct'); BasisString2(1:7) = 'cc-pVTZ'
  CASE('ccT'); BasisString2(1:7) = 'cc-pVTZ'
  CASE('CCT'); BasisString2(1:7) = 'cc-pVTZ'
  CASE('ccq'); BasisString2(1:7) = 'cc-pVQZ'
  CASE('ccQ'); BasisString2(1:7) = 'cc-pVQZ'
  CASE('CCQ'); BasisString2(1:7) = 'cc-pVQZ'
  CASE('cc5'); BasisString2(1:7) = 'cc-pV5Z'
  CASE('CC5'); BasisString2(1:7) = 'cc-pV5Z'
  CASE('cc6'); BasisString2(1:7) = 'cc-pV6Z'
  CASE('CC6'); BasisString2(1:7) = 'cc-pV6Z'
  CASE DEFAULT
   print*,'STRING GIVEN:',BasisString
   print*,'unknown basis only ccd to cc6 implemented'
   call lsquit('unknown basis',-1)
  END SELECT
  nUCharge=0
  DO J=1,50
     ATOMloop: DO I=1,nAtoms
        IF(NINT(CHARGE(I)).EQ.J)THEN
           nUCharge=nUCharge+1
           EXIT ATOMloop
        ENDIF
     ENDDO ATOMloop
  ENDDO
  call mem_alloc(UNIQUECHARGES,nUCharge)
  nUCharge=0
  DO J=1,50
     ATOMloop2: DO I=1,nAtoms
        IF(NINT(CHARGE(I)).EQ.J)THEN
           nUCharge=nUCharge+1
           UNIQUECHARGES(nUCharge) = J
           EXIT ATOMloop2
        ENDIF
     ENDDO ATOMloop2
  ENDDO

  call nullifyMainBasis(input%BASIS)
  LIBRARY(RegBasParam)%BASISSETNAME = BasisString2
  LIBRARY(RegBasParam)%nbasissets = 1
  LIBRARY(RegBasParam)%nCharges(1) = nUCharge
  LIBRARY(RegBasParam)%Charges(1,1:nUCharge) = UNIQUECHARGES(1:nUCharge)
  LIBRARY(RegBasParam)%pointcharges = .FALSE.
  LIBRARY(RegBasParam)%phantom = .FALSE.
  LIBRARY(RegBasParam)%DunningsBasis = .TRUE.
  LIBRARY(RegBasParam)%GeminalScalingFactor = 1.0E0_realk
  call mem_dealloc(UNIQUECHARGES)  
  call Build_BASIS(LUPRI,0,input%MOLECULE,input%BASIS%BINFO(RegBasParam),LIBRARY,&
       &'REGULAR  ',.FALSE.,.FALSE.,.FALSE.,.TRUE.,RegBasParam)
  input%BASIS%WBASIS(RegBasParam)=.TRUE.

  CALL typedef_init_setting(setting)
  CALL typedef_set_default_setting(setting,input)
  setting%SCHEME%DoSpherical = .TRUE.
#ifdef VAR_MPI
  call get_rank_for_comm(MPI_COMM_LSDALTON,mynum)
  CALL get_size_for_comm(MPI_COMM_LSDALTON,nodtot)
  setting%numNodes = nodtot 
  setting%node = mynum
  setting%comm = MPI_COMM_LSDALTON
#endif

end subroutine build_setting_from_scratch
