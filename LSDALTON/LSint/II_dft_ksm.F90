!> Module containing soubroutine for calculation of the exchange-correlation contribution to KS-matrix
MODULE IIDFTKSM
use precision
use TYPEDEFTYPE, only: LSSETTING
use IIDFTINT, only: II_DFTINT, TEST_NELECTRONS
use IIDFTD, only: II_DFT_DISP
use dft_type
use dft_memory_handling
use IIDFTKSMWORK
use IIABSVALINT
!WARNING you must not add memory_handling, all memory goes through 
!grid_memory_handling  module so as to determine the memory used in this module.
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_mod
  use LSparameters,only: LSMPI_IIDFTKSM,IIDFTGEO,IIDFTLIN,IIDFTQRS,&
       & IIDFTMAG,IIDFTMAL,IIDFTGKS,IIDFTGLR,LSMPI_IIDFTKSME
#endif
  use xcfun_host,only: xcfun_type_gga,xcfun_type_lda,xcfun_type_metagga,&
       & xcfun_host_set_order, USEXCFUN
CONTAINS
!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Main kohn-sham matrix driver - calls the nummerical integrater with the 
!> name of a worker routine (II_DFT_KSMGGA,II_DFT_KSMLDA,II_DFT_KSMGGAUNRES,
!> II_DFT_KSMLDAUNRES) which does the work for each grid point
!>
SUBROUTINE II_DFT_KSM(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,ENERGY,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY(ndmat)
!> Unrestricted
LOGICAL :: UNRES
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE
LOGICAL          :: DOGGA,DOLND,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP,IDMAT,NGEODRV,IDMAT1
REAL(REALK)      :: SUM,NELE,ERROR,DFTELS,valMPI(2*ndmat)
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk

call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
NGEODRV=0
DOLND=.FALSE.
USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
!MPI Specific
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(LSMPI_IIDFTKSM,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif
call mem_dft_alloc(DFTDATA%ENERGY,ndmat)
DFTDATA%ENERGY = 0E0_realk
IF(DOGGA) THEN
   IF(USEXCFUN)call xcfun_host_set_order(1,UNRES,xcfun_type_gga)
   DFTDATA%nWorkNactBastNblen = 5
   DFTDATA%nWorkNactBast = 1
   DFTDATA%nWorkNactBastNactBast = 1
   IF(UNRES)THEN !UNRESTRICTED
      IF(DFTDATA%LB94.OR.DFTDATA%CS00)THEN
         CALL LSQUIT('LB94 and CS00 does currently not work with unrestricted',-1)
      ENDIF
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
           & II_DFT_KSMGGAUNRES,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ELSE
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
           & II_DFT_KSMGGA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ENDIF
#ifdef VAR_MPI
!=================MPI Specific================================
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%FKSM,NBAST,NBAST,NDMAT,infpar%master,setting%comm)   
      DO IDMAT=1,ndmat
         valMPI(IDMAT)=DFTDATA%ENERGY(IDMAT) 
      ENDDO
      DO IDMAT=1,ndmat
         valMPI(IDMAT+NDMAT)=ELECTRONS(IDMAT)
      ENDDO
      CALL lsmpi_reduction(valMPI,2*ndmat,infpar%master,setting%comm)
      DO IDMAT=1,ndmat
         DFTDATA%ENERGY(IDMAT) = valMPI(IDMAT)
      ENDDO
      DO IDMAT=1,ndmat
         ELECTRONS(IDMAT) = valMPI(IDMAT+NDMAT)
      ENDDO
   ENDIF
!=============================================================
#endif
   DO IDMAT = 1,NDMAT
      DO I = 1, NBAST
         DO J = 1, I-1
            AVERAG = 0.5*(DFTDATA%FKSM(J,I,IDMAT) + DFTDATA%FKSM(I,J,IDMAT))
            DFTDATA%FKSM(J,I,IDMAT) = AVERAG
            DFTDATA%FKSM(I,J,IDMAT) = AVERAG
         END DO
      END DO
   END DO
ELSE
   IF(DFTDATA%LB94.OR.DFTDATA%CS00)THEN
      CALL LSQUIT('LB94 and CS00 does not work with a pure LDA functional',-1)
   ENDIF
   IF(USEXCFUN)call xcfun_host_set_order(1,UNRES,xcfun_type_lda)
   DFTDATA%nWorkNactBastNblen = 2
   DFTDATA%nWorkNactBast = 1
   DFTDATA%nWorkNactBastNactBast = 1
   IF(UNRES)THEN !UNRESTRICTED
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
           & II_DFT_KSMLDAUNRES,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ELSE !DEFAULT (RESTRICTED)
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
           & II_DFT_KSMLDA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   END IF
#ifdef VAR_MPI
!=================MPI Specific================================
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%FKSM,NBAST,NBAST,NDMAT,infpar%master,setting%comm)
      DO IDMAT=1,ndmat
         valMPI(IDMAT)=DFTDATA%ENERGY(IDMAT) 
      ENDDO
      DO IDMAT=1,ndmat
         valMPI(IDMAT+NDMAT)=ELECTRONS(IDMAT)
      ENDDO
      CALL lsmpi_reduction(valMPI,2*ndmat,infpar%master,setting%comm)
      DO IDMAT=1,ndmat
         DFTDATA%ENERGY(IDMAT) = valMPI(IDMAT)
      ENDDO
      DO IDMAT=1,ndmat
         ELECTRONS(IDMAT) = valMPI(IDMAT+NDMAT)
      ENDDO
   ENDIF
!=============================================================
#endif
ENDIF

#ifdef VAR_MPI
IF(setting%node.EQ.infpar%master) THEN
   !master test for the correct number of electrons
   IF (SETTING%SCHEME%DFT%testNelectrons) THEN
     CALL TEST_NELECTRONS(ELECTRONS,SETTING%MOLECULE(1)%p%NELECTRONS,&
        &                 ndmat,SETTING%SCHEME%DFT%DFTELS,unres,LUPRI)
   ENDIF
#endif
   NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
   IF(IPRINT.GE. 0)THEN
      IF(UNRES)THEN
         DO IDMAT=1,NDMAT/2
            IDMAT1 = 1 + (IDMAT-1)*2
            WRITE(LUPRI,'(A,2F20.14,A,E9.2)')&
                 &     '  KS electrons/energy:', ELECTRONS(IDMAT1), DFTDATA%ENERGY(IDMAT1),&
                 &     ' rel.err:', (ELECTRONS(IDMAT1)-NELE)/(NELE)
         ENDDO
      ELSE
         DO IDMAT=1,NDMAT
            WRITE(LUPRI,'(A,2F20.14,A,E9.2)')&
                 &     '  KS electrons/energy:', ELECTRONS(IDMAT), DFTDATA%ENERGY(IDMAT),&
                 &     ' rel.err:', (ELECTRONS(IDMAT)-NELE)/(NELE)
         ENDDO
      ENDIF
   ENDIF
   
   ENERGY = DFTDATA%ENERGY
   call mem_dft_dealloc(DFTDATA%ENERGY)
#ifdef VAR_MPI
ENDIF
#endif

END SUBROUTINE II_DFT_KSM

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Main kohn-sham matrix driver - calls the nummerical integrater with the 
!> name of a worker routine (II_DFT_KSMGGA,II_DFT_KSMLDA,II_DFT_KSMGGAUNRES,
!> II_DFT_KSMLDAUNRES) which does the work for each grid point
!>
SUBROUTINE II_DFT_ABSVAL_OVERLAP(SETTING,LUPRI,IPRINT,nbast,nMO1,nMO2,CMAT1,CMAT2,ABSVALOVERLAP,SameCmat)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of MOcoef (nocc or nvirt or nbast)
INTEGER,intent(in) :: NMO1,NMO2
!> MO coef matrix 1
REAL(REALK),intent(in) :: CMAT1(NBAST,NMO1)
!> MO coef matrix 2
REAL(REALK),intent(in) :: CMAT2(NBAST,NMO2)
!> absolute valued Overlap matrix
logical,intent(in) :: SameCmat
!> absolute valued Overlap matrix
REAL(REALK),intent(inout) :: ABSVALOVERLAP(NMO1,NMO2)
!
LOGICAL          :: USE_MPI
REAL(REALK)      :: DFTHRI
USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
!MPI Specific
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(LSMPI_IIDFTABSVALOVERLAP,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_ABSVAL_masterToSlave(SETTING,LUPRI,IPRINT,nbast,nMO1,nMO2,&
        & CMAT1,CMAT2,ABSVALOVERLAP,SameCmat,setting%comm)
ENDIF
#endif

DFTHRI = 1.0E-16_realk

call II_ABSVALINT(LUPRI,IPRINT,SETTING,CMAT1,CMAT2,NBAST,NMO1,NMO2,&
     & ABSVALOVERLAP,USE_MPI,DFTHRI,SameCmat)

#ifdef VAR_MPI
!=================MPI Specific================================
IF(USE_MPI)THEN
   call lsmpi_barrier(setting%comm)    
   CALL lsmpi_reduction(ABSVALOVERLAP,NMO1,NMO2,infpar%master,setting%comm)
ENDIF
!=============================================================
#endif
END SUBROUTINE II_DFT_ABSVAL_OVERLAP

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Main kohn-sham matrix driver - calls the nummerical integrater with the 
!> name of a worker routine (II_DFT_KSMGGA,II_DFT_KSMLDA,II_DFT_KSMGGAUNRES,
!> II_DFT_KSMLDAUNRES) which does the work for each grid point
!>
SUBROUTINE II_DFT_KSME(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,ENERGY,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> kohn-sham energy
REAL(REALK),intent(inout) :: ENERGY(ndmat)
!> Unrestricted
LOGICAL :: UNRES
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE
LOGICAL          :: DOGGA,DOLND,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP,IDMAT,NGEODRV,IDMAT1
REAL(REALK)      :: SUM,NELE,ERROR,DFTELS,valMPI(2*ndmat)
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk

call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
NGEODRV=0
DOLND=.FALSE.

USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
!MPI Specific
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(LSMPI_IIDFTKSME,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif
call mem_dft_alloc(DFTDATA%ENERGY,ndmat)
DFTDATA%ENERGY = 0E0_realk
IF(DOGGA) THEN
   IF(USEXCFUN)call xcfun_host_set_order(0,UNRES,xcfun_type_gga)
   DFTDATA%nWorkNactBastNblen = 0
   DFTDATA%nWorkNactBast = 0
   DFTDATA%nWorkNactBastNactBast = 0
   IF(UNRES)THEN !UNRESTRICTED
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
           & II_DFT_KSMEGGAUNRES,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ELSE
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
           & II_DFT_KSMEGGA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ENDIF
#ifdef VAR_MPI
!=================MPI Specific================================
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      DO IDMAT=1,ndmat
         valMPI(IDMAT)=DFTDATA%ENERGY(IDMAT) 
      ENDDO
      DO IDMAT=1,ndmat
         valMPI(IDMAT+NDMAT)=ELECTRONS(IDMAT)
      ENDDO
      CALL lsmpi_reduction(valMPI,2*ndmat,infpar%master,setting%comm)
      DO IDMAT=1,ndmat
         DFTDATA%ENERGY(IDMAT) = valMPI(IDMAT)
      ENDDO
      DO IDMAT=1,ndmat
         ELECTRONS(IDMAT) = valMPI(IDMAT+NDMAT)
      ENDDO
   ENDIF
!=============================================================
#endif
ELSE 
   IF(USEXCFUN)call xcfun_host_set_order(0,UNRES,xcfun_type_gga)
   DFTDATA%nWorkNactBastNblen = 0
   DFTDATA%nWorkNactBast = 0
   DFTDATA%nWorkNactBastNactBast = 0
   IF(UNRES)THEN !UNRESTRICTED
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
           & II_DFT_KSMELDAUNRES,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ELSE !DEFAULT (RESTRICTED)
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
           & II_DFT_KSMELDA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   END IF
#ifdef VAR_MPI
!=================MPI Specific================================
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      DO IDMAT=1,ndmat
         valMPI(IDMAT)=DFTDATA%ENERGY(IDMAT) 
      ENDDO
      DO IDMAT=1,ndmat
         valMPI(IDMAT+NDMAT)=ELECTRONS(IDMAT)
      ENDDO
      CALL lsmpi_reduction(valMPI,2*ndmat,infpar%master,setting%comm)
      DO IDMAT=1,ndmat
         DFTDATA%ENERGY(IDMAT) = valMPI(IDMAT)
      ENDDO
      DO IDMAT=1,ndmat
         ELECTRONS(IDMAT) = valMPI(IDMAT+NDMAT)
      ENDDO
   ENDIF
!=============================================================
#endif
ENDIF

#ifdef VAR_MPI
IF(setting%node.EQ.infpar%master) THEN
   !master test for the correct number of electrons
   CALL TEST_NELECTRONS(ELECTRONS,SETTING%MOLECULE(1)%p%NELECTRONS,&
        &ndmat,SETTING%SCHEME%DFT%DFTELS,UNRES,LUPRI)
#endif
   NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
   IF(IPRINT.GE. 0)THEN
      IF(UNRES)THEN
         DO IDMAT=1,NDMAT/2
            IDMAT1 = 1 + (IDMAT-1)*2
            WRITE(LUPRI,'(A,2F20.14,A,E9.2)')&
                 &     '  KS electrons/energy:', ELECTRONS(IDMAT1), DFTDATA%ENERGY(IDMAT1),&
                 &     ' rel.err:', (ELECTRONS(IDMAT1)-NELE)/(NELE)
         ENDDO
      ELSE
         DO IDMAT=1,NDMAT
            WRITE(LUPRI,'(A,2F20.14,A,E9.2)')&
                 &     '  KS electrons/energy:', ELECTRONS(IDMAT), DFTDATA%ENERGY(IDMAT),&
                 &     ' rel.err:', (ELECTRONS(IDMAT)-NELE)/(NELE)
         ENDDO
      ENDIF
   ENDIF
   
   ENERGY = DFTDATA%ENERGY
   call mem_dft_dealloc(DFTDATA%ENERGY)

#ifdef VAR_MPI
ENDIF !(setting%node.EQ.infpar%master) THEN
#endif

END SUBROUTINE II_DFT_KSME

!> \brief main kohn-sham matrix driver
!> \author T. Kjaergaard
!> \date 2008
!>
!> Main molecular gradient driver - calls the nummerical integrater with the 
!> name of a worker routine (II_geoderiv_molgrad_worker) which 
!> does work for each grid point
!>
SUBROUTINE II_dft_geoderiv_molgrad(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDmat
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> Unrestricted
LOGICAL :: UNRES
!
REAL(REALK)      :: ELE
INTEGER          :: I,J,KMAX,MXPRIM
INTEGER          :: NDERIV
LOGICAL          :: DOGGA,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP
REAL(REALK)      :: SUM,NELE
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk
USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(IIDFTGEO,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif

DFTDATA%nWorkNactBastNblen = 0
DFTDATA%nWorkNactBast = 0
DFTDATA%nWorkNactBastNactBast = 0

call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
IF(DOGGA) THEN
   IF(USEXCFUN)call xcfun_host_set_order(1,UNRES,xcfun_type_gga)
   CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,1,.FALSE.,&
        & II_geoderiv_molgrad_worker_gga,DFTDATA,UNRES,ELECTRONS,USE_MPI)
ELSE
   IF(USEXCFUN)call xcfun_host_set_order(1,UNRES,xcfun_type_lda)
   CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,1,.FALSE.,&
        & II_geoderiv_molgrad_worker_lda,DFTDATA,UNRES,ELECTRONS,USE_MPI)
ENDIF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%GRAD,3,SETTING%MOLECULE(1)%p%NATOMS,infpar%master,setting%comm)
   ENDIF
#else
NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
IF(IPRINT.GT. 1) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
     &     'KS electrons:', ELECTRONS(1),' rel.err:', (ELECTRONS(1)-NELE)/(NELE)
#endif

#ifdef VAR_MPI
IF (setting%node.EQ.infpar%master) THEN
#endif
   ! add eventually empirical dispersion correction
   NDERIV=1
   CALL II_DFT_DISP(SETTING,DFTDATA%GRAD,3,SETTING%MOLECULE(1)%p%NATOMS,NDERIV,LUPRI)
#ifdef VAR_MPI
ENDIF
#endif

END SUBROUTINE II_DFT_GEODERIV_MOLGRAD

!> \brief Main linear response driver
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE ii_dft_linrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> Unrestricted
LOGICAL :: UNRES
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL          :: DOGGA,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP,IDMAT,IBMAT
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk

USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(IIDFTLIN,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif

call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
IF(DOGGA) THEN
   IF(USEXCFUN)call xcfun_host_set_order(2,UNRES,xcfun_type_gga)
   IF(UNRES)THEN !UNRESTRICTED
      DFTDATA%nWorkNactBastNblen = 0
      DFTDATA%nWorkNactBast = 0
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_LINRSPGGAUNRES,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ELSE
      DFTDATA%nWorkNactBastNblen = 5
      DFTDATA%nWorkNactBast = 1
      DFTDATA%nWorkNactBastNactBast = 1
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_LINRSPGGA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ENDIF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%FKSM,NBAST,NBAST,DFTDATA%NBMAT,infpar%master,setting%comm)
   ENDIF
#endif
   DO IBMAT = 1,DFTDATA%NBMAT
      DO I = 1, NBAST
         DO J = 1, I-1
            AVERAG = 0.5*(DFTDATA%FKSM(J,I,IBMAT) + DFTDATA%FKSM(I,J,IBMAT))
            DFTDATA%FKSM(J,I,IBMAT) = AVERAG
            DFTDATA%FKSM(I,J,IBMAT) = AVERAG
         END DO
      END DO
   END DO
ELSE 
   IF(USEXCFUN)call xcfun_host_set_order(2,UNRES,xcfun_type_lda)
   IF(UNRES)THEN !UNRESTRICTED
      DFTDATA%nWorkNactBastNblen = 0
      DFTDATA%nWorkNactBast = 0
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_LINRSPLDAUNRES,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ELSE !DEFAULT (RESTRICTED)
      DFTDATA%nWorkNactBastNblen = 2
      DFTDATA%nWorkNactBast = 1
      DFTDATA%nWorkNactBastNactBast = 1
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_LINRSPLDA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   END IF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%FKSM,NBAST,NBAST,DFTDATA%NBMAT,infpar%master,setting%comm)
   ENDIF
#endif
ENDIF
#ifndef VAR_MPI
NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
IF(IPRINT.GT. 1) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
     &     'KS electrons:', ELECTRONS(1),&
     &     ' rel.err:', (ELECTRONS(1)-NELE)/(NELE)
#endif
END SUBROUTINE II_DFT_LINRSP

!> \brief Main quadratic response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ii_dft_quadrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> Unrestricted
LOGICAL :: UNRES
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL          :: DOGGA,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP,IDMAT,IBMAT
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk

USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(IIDFTQRS,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif

call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
IF(DOGGA) THEN
   IF(UNRES)THEN !UNRESTRICTED
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
      !             & II_DFT_QUADRSPGGAUNRES,DFTDATA)
      CALL LSQUIT('II_DFT_QUADRSPGGAUNRES NOT DONE YET',lupri)
   ELSE
      DFTDATA%nWorkNactBastNblen = 5
      DFTDATA%nWorkNactBast = 1
      DFTDATA%nWorkNactBastNactBast = 1
      IF(USEXCFUN)call xcfun_host_set_order(3,UNRES,xcfun_type_gga)
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_QUADRSPGGA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
      !        CALL LSQUIT('II_DFT_QUADRSPGGA NOT DONE YET',lupri)
   ENDIF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%FKSM,NBAST,NBAST,NDMAT,infpar%master,setting%comm)
   ENDIF
#endif
   DO IDMAT = 1,DFTDATA%NDMAT
      DO I = 1, NBAST
         DO J = 1, I-1
            AVERAG = 0.5*(DFTDATA%FKSM(J,I,IDMAT) + DFTDATA%FKSM(I,J,IDMAT))
            DFTDATA%FKSM(J,I,IDMAT) = AVERAG
            DFTDATA%FKSM(I,J,IDMAT) = AVERAG
         END DO
      END DO
   END DO
ELSE 
   IF(UNRES)THEN !UNRESTRICTED
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
      !             & II_DFT_QUADRSPLDAUNRES,DFTDATA)
      CALL LSQUIT('II_DFT_QUADRSPLDAUNRES NOT DONE YET',lupri)
   ELSE !DEFAULT (RESTRICTED)
      DFTDATA%nWorkNactBastNblen = 2
      DFTDATA%nWorkNactBast = 1
      DFTDATA%nWorkNactBastNactBast = 1
      IF(USEXCFUN)call xcfun_host_set_order(3,UNRES,xcfun_type_lda)
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,.FALSE.,&
           & II_DFT_QUADRSPLDA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   END IF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%FKSM,NBAST,NBAST,NDMAT,infpar%master,setting%comm)
   ENDIF
#endif
ENDIF

#ifndef VAR_MPI
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GT. 1) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
       &     'KS electrons:', ELECTRONS(1),&
       &     ' rel.err:', (ELECTRONS(1)-NELE)/(NELE)
#endif

END SUBROUTINE II_DFT_QUADRSP

!> \brief main magnetic derivative kohn-sham matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_magderiv_kohnsham_mat(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> Unrestricted
LOGICAL :: UNRES
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL          :: DOGGA,DOLONDON,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP,IDMAT
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk

USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(IIDFTMAG,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif

call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
DOLONDON = .TRUE.
IF(DOGGA) THEN
   IF(USEXCFUN)call xcfun_host_set_order(1,UNRES,xcfun_type_gga)
   IF(UNRES)THEN !UNRESTRICTED
      CALL LSQUIT('II_DFT_magderiv_kohnshamGGAUNRES not implemented',lupri)
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
      !             & II_DFT_magderiv_kohnshamGGAUNRES,DFTDATA)
   ELSE
      DFTDATA%nWorkNactBastNblen = 20
      DFTDATA%nWorkNactBast = 1
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
           & II_DFT_magderiv_kohnshamGGA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ENDIF
ELSE 
   IF(USEXCFUN)call xcfun_host_set_order(1,UNRES,xcfun_type_lda)
   IF(UNRES)THEN !UNRESTRICTED
      CALL LSQUIT('II_DFT_magderiv_kohnshamLDAUNRES not implemented',lupri)
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
      !             & II_DFT_magderiv_kohnshamLDAUNRES,DFTDATA)
   ELSE !DEFAULT (RESTRICTED)
      DFTDATA%nWorkNactBastNblen = 0
      DFTDATA%nWorkNactBast = 0
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
           & II_DFT_magderiv_kohnshamLDA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   END IF
ENDIF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%FKSM,NBAST,NBAST,NDMAT*3,infpar%master,setting%comm)
   ENDIF
#endif
!ANTISYMMETRIZE
DO IDMAT = 1,NDMAT*3
   DO I = 1, NBAST
      DO J = 1, I-1
         AVERAG = 0.5*(DFTDATA%FKSM(J,I,IDMAT) - DFTDATA%FKSM(I,J,IDMAT))
         DFTDATA%FKSM(J,I,IDMAT) = AVERAG
         DFTDATA%FKSM(I,J,IDMAT) = -AVERAG
      END DO
      DFTDATA%FKSM(I,I,IDMAT) = 0E0_realk  
   END DO
END DO
END SUBROUTINE II_DFT_MAGDERIV_KOHNSHAM_MAT

!> \brief Main driver for the calculation of magnetic derivative of the linear response
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ii_dft_magderiv_linrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> Unrestricted
LOGICAL :: UNRES
!
INTEGER          :: I,J,KMAX,natoms,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL          :: DOGGA,DOLONDON,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP,IDMAT,IBMAT,nbmat
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk

USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(IIDFTMAL,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif

NBMAT = DFTDATA%NBMAT
call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
DOLONDON = .TRUE.
IF(DOGGA) THEN
   IF(USEXCFUN)call xcfun_host_set_order(2,UNRES,xcfun_type_gga)
   IF(UNRES)THEN !UNRESTRICTED
      CALL LSQUIT('II_DFT_MAGDERIV_LINRSPGGAUNRES not implemented',lupri)
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
      !             & II_DFT_MAGDERIV_LINRSPGGAUNRES,DFTDATA)
   ELSE
      DFTDATA%nWorkNactBastNblen = 0
      DFTDATA%nWorkNactBast = 0
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
           & II_DFT_MAGDERIV_LINRSPGGA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ENDIF
ELSE 
   IF(USEXCFUN)call xcfun_host_set_order(2,UNRES,xcfun_type_lda)
   IF(UNRES)THEN !UNRESTRICTED
      CALL LSQUIT('II_DFT_MAGDERIV_LINRSPLDAUNRES not implemented',lupri)
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
      !             & II_DFT_MAGDERIV_LINRSPLDAUNRES,DFTDATA)
   ELSE !DEFAULT (RESTRICTED)
      DFTDATA%nWorkNactBastNblen = 0
      DFTDATA%nWorkNactBast = 0
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,0,DOLONDON,&
           & II_DFT_MAGDERIV_LINRSPLDA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   END IF
ENDIF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%FKSM,NBAST,NBAST,NBMAT*3,infpar%master,setting%comm)
      CALL lsmpi_reduction(DFTDATA%FKSMS,NBAST,NBAST,NBMAT*3,infpar%master,setting%comm)
   ENDIF
#endif
  !ANTISYMMETRIZE FKSM
  DO IBMAT = 1,NBMAT*3
     DO I = 1, NBAST
        DO J = 1, I-1
           AVERAG = 0.5*(DFTDATA%FKSM(J,I,IBMAT) - DFTDATA%FKSM(I,J,IBMAT))
           DFTDATA%FKSM(J,I,IBMAT) =  AVERAG
           DFTDATA%FKSM(I,J,IBMAT) = -AVERAG
        END DO
        DFTDATA%FKSM(I,I,IBMAT) = 0E0_realk  
     END DO
  END DO
  !SYMMETRIZE FKSMS
  DO IBMAT = 1,NBMAT*3
     DO I = 1, NBAST
        DO J = 1, I-1
           AVERAG = 0.5*(DFTDATA%FKSMS(J,I,IBMAT) + DFTDATA%FKSMS(I,J,IBMAT))
           DFTDATA%FKSMS(J,I,IBMAT) =  AVERAG
           DFTDATA%FKSMS(I,J,IBMAT) =  AVERAG
        END DO
     END DO
  END DO
#ifndef VAR_MPI
  NELE = REAL(SETTING%MOLECULE(1)%p%NELECTRONS)
  IF(IPRINT.GT. 1) WRITE(LUPRI,'(A,F20.14,A,E9.2)')&
       &     'KS electrons:', ELECTRONS(1),&
       &     ' rel.err:', (ELECTRONS(1)-NELE)/(NELE)
#endif

END SUBROUTINE II_DFT_MAGDERIV_LINRSP

!> \brief main geometrical derivative of the kohn-sham matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_kohnsham_mat(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> Unrestricted
LOGICAL :: UNRES
!
INTEGER          :: I,J,KMAX,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL          :: DOGGA,DOLONDON,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP,IDMAT,NGEODERIV
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk

USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(IIDFTGKS,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif

call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
DOLONDON = .FALSE.
NGEODERIV = 1
IF(DOGGA) THEN
   IF(USEXCFUN)call xcfun_host_set_order(2,UNRES,xcfun_type_gga)
   IF(UNRES)THEN !UNRESTRICTED
      CALL LSQUIT('II_DFT_geoderiv_kohnshamGGAUNRES not implemented',lupri)
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
      !             & II_DFT_geoderiv_kohnshamGGAUNRES,DFTDATA)
   ELSE
      !        CALL LSQUIT('II_DFT_geoderiv_kohnshamGGA not implemented',lupri)
      DFTDATA%nWorkNactBastNblen = 14
      DFTDATA%nWorkNactBast = 1
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
           & II_DFT_geoderiv_kohnshamGGA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ENDIF
ELSE 
   IF(USEXCFUN)call xcfun_host_set_order(2,UNRES,xcfun_type_lda)
   IF(UNRES)THEN !UNRESTRICTED
      CALL LSQUIT('II_DFT_geoderiv_kohnshamLDAUNRES not implemented',lupri)
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
      !             & II_DFT_geoderiv_kohnshamLDAUNRES,DFTDATA)
   ELSE !DEFAULT (RESTRICTED)
      DFTDATA%nWorkNactBastNblen = 0
      DFTDATA%nWorkNactBast = 0
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
           & II_DFT_geoderiv_kohnshamLDA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   END IF
ENDIF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%GRAD,3,SETTING%MOLECULE(1)%p%NATOMS,infpar%master,setting%comm)
   ENDIF
#endif

END SUBROUTINE II_DFT_GEODERIV_KOHNSHAM_MAT

!> \brief main geometrical derivative linear response driver
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFT_geoderiv_linrspgrad(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
  IMPLICIT NONE
!> contains info about the molecule,basis and dft grid parameters
TYPE(LSSETTING),intent(inout)  :: SETTING
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of density matrices
INTEGER,intent(in) :: NDMAT
!> Density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> contains all info required which is not directly related to the integration
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> Unrestricted
LOGICAL :: UNRES
!
INTEGER          :: I,J,KMAX,MXPRIM
REAL(REALK)      :: AVERAG,ELE,SUM,NELE
LOGICAL          :: DOGGA,DOLONDON,USE_MPI,DOMETA
INTEGER          :: GRDONE,NHTYP,IDMAT,NGEODERIV
REAL(REALK)      :: ELECTRONS(ndmat)
ELECTRONS = 0E0_realk

USE_MPI = SETTING%SCHEME%doMPI
IF(SETTING%MOLECULE(1)%p%NATOMS.EQ.1)USE_MPI=.FALSE.
#ifdef VAR_MPI
IF (setting%node.EQ.infpar%master) THEN
   IF(USE_MPI)call ls_mpibcast(IIDFTGLR,infpar%master,setting%comm)
   IF(USE_MPI)call lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,IPRINT,nbast,&
        & ndmat,DMAT,DFTDATA,UNRES,setting%comm)
ENDIF
#endif

call DFT_DOGGA_DOMETA(DOGGA,DOMETA)
DOLONDON = .FALSE.
NGEODERIV = 1
IF(DOGGA) THEN
   IF(USEXCFUN)call xcfun_host_set_order(3,UNRES,xcfun_type_gga)
   IF(UNRES)THEN !UNRESTRICTED
      CALL LSQUIT('II_DFT_geoderiv_linrspGGAUNRES not implemented',lupri)
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
      !             & II_DFT_geoderiv_linrspGGAUNRES,DFTDATA)
   ELSE
      !        CALL LSQUIT('II_DFT_geoderiv_linrspGGA not implemented',lupri)
      DFTDATA%nWorkNactBastNblen = 5
      DFTDATA%nWorkNactBast = 1
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
           & II_DFT_geoderiv_linrspGGA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   ENDIF
ELSE 
   IF(USEXCFUN)call xcfun_host_set_order(3,UNRES,xcfun_type_lda)
   IF(UNRES)THEN !UNRESTRICTED
      CALL LSQUIT('II_DFT_geoderiv_linrspLDAUNRES not implemented',lupri)
      !        CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
      !             & II_DFT_geoderiv_linrspLDAUNRES,DFTDATA)
   ELSE !DEFAULT (RESTRICTED)
      DFTDATA%nWorkNactBastNblen = 0
      DFTDATA%nWorkNactBast = 0
      DFTDATA%nWorkNactBastNactBast = 0
      CALL II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODERIV,DOLONDON,&
           & II_DFT_geoderiv_linrspLDA,DFTDATA,UNRES,ELECTRONS,USE_MPI)
   END IF
ENDIF
#ifdef VAR_MPI
   IF(USE_MPI)THEN
      call lsmpi_barrier(setting%comm)    
      CALL lsmpi_reduction(DFTDATA%GRAD,3,SETTING%MOLECULE(1)%p%NATOMS,infpar%master,setting%comm)
   ENDIF
#endif

END SUBROUTINE II_DFT_GEODERIV_LINRSPGRAD

END MODULE IIDFTKSM

#ifdef VAR_MPI
SUBROUTINE lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
use lsmpi_mod
use infpar_module
use typedef
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
TYPE(DFTDATATYPE) :: DFTDATA
!INTEGERS
call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
CALL ls_mpi_buffer(LUPRI,infpar%master)
CALL ls_mpi_buffer(IPRINT,infpar%master)
CALL ls_mpi_buffer(NBAST,infpar%master)
CALL ls_mpi_buffer(NDMAT,infpar%master)
CALL ls_mpi_buffer(UNRES,infpar%master)
END SUBROUTINE LSMPI_BUFFERDFTESSENTIALS1

SUBROUTINE lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
use lsmpi_mod
use infpar_module
use typedef
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
TYPE(LSSETTING)  :: SETTING
REAL(REALK)      :: DMAT(NBAST,NBAST,NDMAT)
TYPE(DFTDATATYPE) :: DFTDATA
!REAL(REALK), pointer
CALL ls_mpi_buffer(DMAT,NBAST,NBAST,NDMAT,infpar%master)
!SETTING
CALL mpicopy_setting(setting,comm,.FALSE.)
!DFTDATA
call mpicopy_dftdata(dftdata,setting%node)
call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
END SUBROUTINE LSMPI_BUFFERDFTESSENTIALS2

!===========================================================
! Generic MPI routine for master 
! As It always broadcast the same stuff
!==========================================================
subroutine lsmpi_XCgeneric_masterToSlave(SETTING,LUPRI,&
     & IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES,comm)
use lsmpi_mod
use infpar_module
use typedef
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK)      :: DMAT(NBAST,NBAST,NDMAT)
TYPE(DFTDATATYPE) :: DFTDATA
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
end subroutine lsmpi_XCgeneric_masterToSlave

subroutine lsmpi_ABSVAL_masterToSlave(SETTING,LUPRI,IPRINT,nbast,nMO1,nMO2,&
        & CMAT1,CMAT2,ABSVALOVERLAP,SameCmat,comm)
use lsmpi_mod
use infpar_module
use typedef
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,nMO1,nMO2
integer(kind=ls_mpik) :: comm
TYPE(LSSETTING)  :: SETTING
REAL(REALK)      :: CMAT1(NBAST,NMO1)
REAL(REALK)      :: CMAT2(NBAST,NMO2),ABSVALOVERLAP(NMO1,NMO2)
logical          :: SameCmat

call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
CALL ls_mpi_buffer(LUPRI,infpar%master)
CALL ls_mpi_buffer(IPRINT,infpar%master)
CALL ls_mpi_buffer(NBAST,infpar%master)
CALL ls_mpi_buffer(NMO1,infpar%master)
CALL ls_mpi_buffer(NMO2,infpar%master)
CALL ls_mpi_buffer(SameCmat,infpar%master)
CALL ls_mpi_buffer(CMAT1,NBAST,NMO1,infpar%master)
IF(.NOT.SameCmat)THEN
   CALL ls_mpi_buffer(CMAT2,NBAST,NMO2,infpar%master)
ENDIF
CALL ls_mpi_buffer(ABSVALOVERLAP,NMO1,NMO2,infpar%master)
!SETTING
CALL mpicopy_setting(setting,comm,.FALSE.)
call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
end subroutine lsmpi_ABSVAL_masterToSlave

!===========================================================
! MPI II_DFT_KSM Slave
!==========================================================
subroutine lsmpi_II_DFT_KSM_Slave(comm)
use lsmpi_mod
use infpar_module
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
REAL(REALK),pointer :: ENERGY(:)
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
call mem_dft_alloc(ENERGY,NDMAT)
ENERGY=0E0_realk
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call II_DFT_KSM(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,ENERGY,UNRES)
call mem_dft_dealloc(ENERGY)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_II_DFT_KSM_Slave

!LSMPI_IIDFTABSVALOVERLAP
subroutine lsmpi_II_DFT_ABSVALOVERLAP_Slave(comm)
  use lsmpi_mod
  use infpar_module
  use typedef
  use IIDFTKSM
  implicit none
  INTEGER :: LUPRI,IPRINT,Nbast,NDMAT,NMO1,NMO2
  integer(kind=ls_mpik) :: comm
  TYPE(LSSETTING)  :: SETTING
  logical :: SameCmat
  REAL(REALK),pointer :: CMAT1(:,:),CMAT2(:,:),ABSVALOVERLAP(:,:)
  call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
  CALL ls_mpi_buffer(LUPRI,infpar%master)
  CALL ls_mpi_buffer(IPRINT,infpar%master)
  CALL ls_mpi_buffer(NBAST,infpar%master)
  CALL ls_mpi_buffer(NMO1,infpar%master)
  CALL ls_mpi_buffer(NMO2,infpar%master)
  CALL ls_mpi_buffer(SameCmat,infpar%master)
  call mem_dft_alloc(CMAT1,NBAST,NMO1)
  CALL ls_mpi_buffer(CMAT1,NBAST,NMO1,infpar%master)
  IF(.NOT.SameCmat)THEN
     call mem_dft_alloc(CMAT2,NBAST,NMO2)
     CALL ls_mpi_buffer(CMAT2,NBAST,NMO2,infpar%master)
  ENDIF
  call mem_dft_alloc(ABSVALOVERLAP,NMO1,NMO2)
  CALL ls_mpi_buffer(ABSVALOVERLAP,NMO1,NMO2,infpar%master)
  CALL mpicopy_setting(setting,comm,.FALSE.)
  call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
  IF(SameCmat)THEN
     call II_DFT_ABSVAL_OVERLAP(SETTING,LUPRI,IPRINT,nbast,nMO1,nMO1,&
          & CMAT1,CMAT1,ABSVALOVERLAP,SameCmat)
  ELSE
     call II_DFT_ABSVAL_OVERLAP(SETTING,LUPRI,IPRINT,nbast,nMO1,nMO2,&
          & CMAT1,CMAT2,ABSVALOVERLAP,SameCmat)
  ENDIF
  call mem_dft_dealloc(CMAT1)
  IF(.NOT.SameCmat)THEN
     call mem_dft_dealloc(CMAT2)
  ENDIF
  call mem_dft_dealloc(ABSVALOVERLAP)
  call typedef_free_setting(SETTING)
end subroutine lsmpi_II_DFT_ABSVALOVERLAP_Slave

subroutine lsmpi_II_DFT_KSME_Slave(comm)
use lsmpi_mod
use infpar_module
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
REAL(REALK),pointer :: ENERGY(:)
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
call mem_dft_alloc(ENERGY,NDMAT)
ENERGY=0E0_realk
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call II_DFT_KSME(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,ENERGY,UNRES)
call mem_dft_dealloc(ENERGY)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_II_DFT_KSME_Slave
!===========================================================
! MPI II_dft_geoderiv_molgrad Slave
!==========================================================
subroutine lsmpi_geoderiv_molgrad_Slave(comm)
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
TYPE(LSSETTING)  :: SETTING
LOGICAL :: UNRES
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call II_dft_geoderiv_molgrad(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_geoderiv_molgrad_Slave
!===========================================================
! MPI ii_dft_linrsp Slave
!==========================================================
subroutine lsmpi_linrsp_Slave(comm)
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call ii_dft_linrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_linrsp_Slave
!===========================================================
! MPI ii_dft_quadrsp Slave
!==========================================================
subroutine lsmpi_quadrsp_Slave(comm)
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call ii_dft_quadrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_quadrsp_Slave
!===========================================================
! MPI II_DFT_magderiv_kohnsham_mat Slave
!==========================================================
subroutine lsmpi_magderiv_Slave(comm)
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call II_DFT_magderiv_kohnsham_mat(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_magderiv_Slave
!===========================================================
! MPI II_dft_magderiv_linrsp Slave
!==========================================================
subroutine lsmpi_magderiv_linrsp_Slave(comm)
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call ii_dft_magderiv_linrsp(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_magderiv_linrsp_Slave
!===========================================================
! MPI II_DFT_geoderiv_kohnsham_mat Slave
!==========================================================
subroutine lsmpi_geoderiv_F_Slave(comm)
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call II_DFT_geoderiv_kohnsham_mat(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_geoderiv_F_Slave
!===========================================================
! MPI II_DFT_geoderiv_linrspgrad Slave
!==========================================================
subroutine lsmpi_geoderiv_G_Slave(comm)
use typedef
use IIDFTKSM
implicit none
INTEGER :: LUPRI,IPRINT,Nbast,NDMAT
integer(kind=ls_mpik) :: comm
LOGICAL :: UNRES
TYPE(LSSETTING)  :: SETTING
REAL(REALK),pointer :: DMAT(:,:,:)
TYPE(DFTDATATYPE) :: DFTDATA
CALL lsmpi_bufferDFTessentials1(LUPRI,IPRINT,NBAST,NDMAT,UNRES,comm)
call mem_dft_alloc(DMAT,NBAST,NBAST,NDMAT)
CALL lsmpi_bufferDFTessentials2(SETTING,LUPRI,IPRINT,NBAST,NDMAT,DMAT,DFTdata,comm)
!CALL FUNCTION
call II_DFT_geoderiv_linrspgrad(SETTING,LUPRI,IPRINT,nbast,ndmat,DMAT,DFTDATA,UNRES)
call mem_dft_dealloc(DMAT)
call typedef_free_setting(SETTING)
call free_dftdata(dftdata)

end subroutine lsmpi_geoderiv_G_Slave
#endif
