module II_XC_interfaceModule
  use precision
  use TYPEDEFTYPE, only: lssetting
  use TYPEDEF, only: GCAO2AO_transform_fullD,AO2GCAO_transform_matrixF
  use Matrix_module, only: matrix
  use Matrix_Operations, only: mat_to_full, mat_init, mat_daxpy, &
       & mat_set_from_full, mat_assign, &
       & mat_free, mat_scal, mtype_unres_dense, matrix_type
  use dft_memory_handling
!  use memory_handling
  use lstiming
  use IIDFTKSM
  use IIDFTINT, only: II_DFTDISP
  use DFT_type
  private
  public :: II_get_xc_Fock_mat,&
       & II_get_AbsoluteValue_overlap, II_get_xc_energy,&
       & II_get_xc_geoderiv_molgrad, II_get_xc_linrsp,&
       & II_get_xc_quadrsp, II_get_xc_magderiv_kohnsham_mat,&
       & II_get_xc_magderiv_linrsp, II_get_xc_geoderiv_FxDgrad,&
       & II_get_xc_geoderiv_GxDgrad
!> @file
!> Interface subroutines for exchange-correlation contributions

!> \brief Calculates the xc contribution to the Kohn-Sham matrix
!> \author T. Kjaergaard
!> \date 2008
INTERFACE II_get_xc_Fock_mat
   MODULE PROCEDURE II_get_xc_Fock_mat_array,II_get_xc_Fock_mat_single
END INTERFACE

CONTAINS

SUBROUTINE II_get_xc_Fock_mat_single(LUPRI,LUERR,SETTING,nbast,D,F,EDFT,ndmat)
IMPLICIT NONE
INTEGER,intent(in)    :: LUPRI
INTEGER,intent(in)    :: LUERR
TYPE(LSSETTING)       :: SETTING
INTEGER,intent(in)    :: nbast
INTEGER,intent(in)    :: ndmat
TYPE(MATRIX),intent(in) :: D
TYPE(MATRIX),intent(inout) :: F
REAL(REALK),intent(inout)  :: EDFT(ndmat)
!
TYPE(MATRIX)  :: Farray(ndmat)
!NOT OPTIMAL USE OF MEMORY OR ANYTHING
call mat_init(Farray(1),F%nrow,F%ncol)
call mat_assign(Farray(1),F)
call II_get_xc_Fock_mat_array(LUPRI,LUERR,SETTING,nbast,(/D/),Farray,EDFT,ndmat)
call mat_assign(F,Farray(1))
call mat_free(Farray(1))
end SUBROUTINE II_get_xc_Fock_mat_single

SUBROUTINE II_get_xc_Fock_mat_array(LUPRI,LUERR,SETTING,nbast,D,F,EDFT,ndmat)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in)    :: LUPRI
!> the logical unit number for the error file
INTEGER,intent(in)    :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER,intent(in)    :: nbast
!> number of Densitymatrices
INTEGER,intent(in)    :: ndmat
!> The density matrix
TYPE(MATRIX),intent(in) :: D(ndmat)
!> The Kohn-Sham matrix
TYPE(MATRIX),intent(inout) :: F(ndmat)
!> The xc contribution to the energy
REAL(REALK),intent(inout)           :: EDFT(ndmat)
!
TYPE(MATRIX)          :: temp
INTEGER               :: i,j,ndmat2,idmat
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:),DmatAO(:,:,:),EDFT2(:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2,CPUTIME,WALLTIME
REAL(REALK)   :: DUMMY(1,1)
LOGICAL               :: UNRES
call time_II_operations1
IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
CALL LS_GETTIM(CPU1,WALL1)
call initDFTdatatype(DFTDATA)
DFTDATA%LB94=SETTING%SCHEME%DFT%LB94
DFTDATA%CS00=SETTING%SCHEME%DFT%CS00
DFTDATA%CS00shift=SETTING%SCHEME%DFT%CS00shift
DFTDATA%CS00eHOMO=SETTING%SCHEME%DFT%CS00eHOMO
DFTDATA%CS00ZND1=SETTING%SCHEME%DFT%CS00ZND1
DFTDATA%CS00ZND2=SETTING%SCHEME%DFT%CS00ZND2
DFTDATA%HFexchangeFac=SETTING%SCHEME%DFT%HFexchangeFac
IF(DFTDATA%CS00)THEN
   IF(ABS(DFTDATA%CS00shift).LT.1.0E-12_realk)THEN
      IF(ABS(DFTDATA%CS00eHOMO).LT.1.0E-12_realk)THEN
         call lsquit('The CS00 keyword only work with a method which calculates the HOMO energy',-1)
      ENDIF
   ENDIF
ENDIF
DFTDATA%nbast = nbast
ndmat2=ndmat
IF(UNRES)ndmat2=2*ndmat
call mem_dft_alloc(EDFT2,ndmat2)
DFTDATA%ndmat = ndmat2
call mem_dft_alloc(Dmat,nbast,nbast,ndmat2)
DFTDATA%nfmat = ndmat2
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,ndmat2)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*ndmat2)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(UNRES)THEN
   do idmat = 1,ndmat
      CALL DCOPY(nbast*nbast,D(idmat)%elms,1,Dmat(:,:,1+(idmat-1)*2),1)
      CALL DCOPY(nbast*nbast,D(idmat)%elmsb,1,Dmat(:,:,2+(idmat-1)*2),1)
   enddo
ELSE !CLOSED_SHELL
   do idmat = 1,ndmat
      call mat_to_full(D(idmat),1E0_realk,Dmat(:,:,idmat))
      CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,idmat),1)
   enddo
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat2,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_KSM(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,EDFT2,UNRES)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

IF(SETTING%SCHEME%DFT%DODISP) THEN 
    ! add empirical dispersion correction \Andreas Krapp
   CALL II_DFTDISP(SETTING,DUMMY,1,1,0,LUPRI,1)
   do idmat = 1,ndmat2
      EDFT2(idmat) = EDFT2(idmat) + SETTING%EDISP
   enddo
ENDIF
call mem_dft_dealloc(DmatAO)

IF(UNRES)THEN
   do idmat = 1,ndmat
      EDFT(idmat) = EDFT2(1+(idmat-1)*2)
   enddo
ELSE
   EDFT = EDFT2
ENDIF
call mem_dft_dealloc(EDFT2)

!WARNING: 
! For a closed shell molecule calculated using an unrestricted and a 
! closed shell aproach  would fullfill
! F(from closed shell) = F_alpha(from onres) + F_beta(from onres) 
! but for some reason this is not what they want in SCF-loop so we 
! multiply the unrestriced result with 2

call mat_init(temp,nbast,nbast)
IF(UNRES)THEN
   do idmat = 1,ndmat
      CALL DCOPY(nbast*nbast,DFTDATA%FKSM(:,:,1+(idmat-1)*2),1,temp%elms,1)
      CALL DCOPY(nbast*nbast,DFTDATA%FKSM(:,:,2+(idmat-1)*2),1,temp%elmsb,1)
      IF(setting%IntegralTransformGC)THEN
         call AO2GCAO_transform_matrixF(temp,setting,lupri)
      ENDIF
      CALL mat_DAXPY(1.0E0_realk,temp,F(idmat))
   enddo
ELSE !CLOSED_SHELL
   do idmat = 1,ndmat
      CALL mat_set_from_full(DFTDATA%FKSM(:,:,idmat),1E0_realk,temp,'XCmat')
      IF(setting%IntegralTransformGC)THEN
         call AO2GCAO_transform_matrixF(temp,setting,lupri)
      ENDIF
      CALL mat_DAXPY(0.5E0_realk,temp,F(idmat))
   enddo
ENDIF
call mat_free(temp)
call mem_dft_dealloc(DFTDATA%FKSM)

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_Fock_mat is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_Fock_mat is  ',WALLTIME,LUPRI)
call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_Fock_mat)

END SUBROUTINE II_get_xc_Fock_mat_array

SUBROUTINE II_get_AbsoluteValue_overlap(LUPRI,LUERR,SETTING,nbast,CMO,S)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in)    :: LUPRI
!> the logical unit number for the error file
INTEGER,intent(in)    :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER,intent(in)    :: nbast
!> The density matrix
TYPE(MATRIX),intent(in) :: CMO
!> The Absolute Valued overlap  matrix
TYPE(MATRIX),intent(inout) :: S
#if MOD_UNRELEASED
!
REAL(REALK),pointer   :: Cmat(:,:),ABSVALOVERLAP(:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2,CPUTIME,WALLTIME
LOGICAL               :: UNRES
!call time_II_operations1
UNRES=.FALSE.
IF(matrix_type .EQ. mtype_unres_dense)UNRES=.TRUE.
call init_dftmemvar
CALL LS_GETTIM(CPU1,WALL1)
call mem_dft_alloc(Cmat,nbast,nbast)
call mem_dft_alloc(ABSVALOVERLAP,nbast,nbast)
CALL LS_DZERO(ABSVALOVERLAP,nbast*nbast)

IF(UNRES)THEN
   call lsquit('not implemeted',-1)
ELSE !CLOSED_SHELL
   call mat_to_full(CMO,1E0_realk,Cmat)
ENDIF

IF(setting%IntegralTransformGC)THEN
   call lsquit('IntegralTransformGC must be false in II_get_AbsoluteValue_overlap',-1)
ENDIF

!chose ABSVAL grid
SETTING%scheme%DFT%igrid = Grid_ABSVAL

!default for fine
SETTING%scheme%DFT%GridObject(Grid_ABSVAL)%RADINT = 2.15443E-17_realk
SETTING%scheme%DFT%GridObject(Grid_ABSVAL)%ANGINT = 47

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_ABSVAL_OVERLAP(SETTING,LUPRI,1,nbast,CMAT,ABSVALOVERLAP)
CALL LSTIMER('ABSVAL-Overlap',TS,TE,LUPRI)

!revert to default grid
SETTING%scheme%DFT%igrid = Grid_Default

call mem_dft_dealloc(Cmat)
CALL mat_set_from_full(ABSVALOVERLAP,1E0_realk,S,'ABSVAL')
call mem_dft_dealloc(ABSVALOVERLAP)

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>> CPU  Time used II_get_AbsoluteValue_overlap',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>> WALL Time used II_get_AbsoluteValue_overlap',WALLTIME,LUPRI)
call stats_dft_mem(lupri)
!call time_II_operations2(JOB_II_get_xc_Fock_mat)
#endif
END SUBROUTINE II_get_AbsoluteValue_overlap

!> \brief Calculates the xc contribution to the Kohn-Sham energy
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_get_xc_energy(LUPRI,LUERR,SETTING,nbast,D,EDFT,ndmat)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D(ndmat)
!> The xc contribution to the energy
REAL(REALK)           :: EDFT(ndmat)
!> number of Densitymatrices
INTEGER               :: ndmat
!
TYPE(MATRIX)          :: temp
INTEGER               :: i,j,ndmat2,idmat
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:),DmatAO(:,:,:),EDFT2(:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2,CPUTIME,WALLTIME
REAL(REALK)   :: DUMMY(1,1),RHOTHR,DFTHRI
LOGICAL               :: UNRES

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
CALL LS_GETTIM(CPU1,WALL1)
call initDFTdatatype(DFTDATA)
DFTDATA%nbast = nbast
ndmat2=ndmat
IF(UNRES)ndmat2=2*ndmat
call mem_dft_alloc(EDFT2,ndmat2)
DFTDATA%ndmat = ndmat2
call mem_dft_alloc(Dmat,nbast,nbast,ndmat2)
!reduced accuracy !!
!SETTING%SCHEME%DFT%RHOTHR = SETTING%SCHEME%DFT%RHOTHR*100.0E0_realk
!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(UNRES)THEN
   do idmat = 1,ndmat
      CALL DCOPY(nbast*nbast,D(idmat)%elms,1,Dmat(:,:,1+(idmat-1)*2),1)
      CALL DCOPY(nbast*nbast,D(idmat)%elmsb,1,Dmat(:,:,2+(idmat-1)*2),1)
   enddo
ELSE !CLOSED_SHELL
   do idmat = 1,ndmat
      call mat_to_full(D(idmat),1E0_realk,Dmat(:,:,idmat))
      CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,idmat),1)
   enddo
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat2,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

!save parameters
DFTHRI = SETTING%scheme%DFT%DFTHRI 
RHOTHR = SETTING%scheme%DFT%RHOTHR 
!chose new loosend parameters
SETTING%scheme%DFT%DFTHRI = DFTHRI*5000
SETTING%scheme%DFT%RHOTHR = RHOTHR*5000

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_KSME(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,EDFT2,UNRES)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

!revert to default grid
SETTING%scheme%DFT%DFTHRI = DFTHRI
SETTING%scheme%DFT%RHOTHR = RHOTHR

IF(SETTING%SCHEME%DFT%DODISP) THEN 
    ! add empirical dispersion correction \Andreas Krapp
   CALL II_DFTDISP(SETTING,DUMMY,1,1,0,LUPRI,1)
   do idmat = 1,ndmat2
      EDFT2(idmat) = EDFT2(idmat) + SETTING%EDISP
   enddo
ENDIF
call mem_dft_dealloc(DmatAO)

IF(UNRES)THEN
   do idmat = 1,ndmat
      EDFT(idmat) = EDFT2(1+(idmat-1)*2)
   enddo
ELSE
   EDFT = EDFT2
ENDIF
call mem_dft_dealloc(EDFT2)

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_energy is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_energy is  ',WALLTIME,LUPRI)
call stats_dft_mem(lupri)

END SUBROUTINE II_get_xc_energy

!> \brief Calculates the xc contribution to the molecular gradient
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_get_xc_geoderiv_molgrad(LUPRI,LUERR,SETTING,nbast,D,grad,natoms)
use Integralparameters
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc contribution to the molecular grad
REAL(REALK)           :: GRAD(3,natoms)
!> Number of atoms
INTEGER               :: natoms
!
INTEGER               :: i,j
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,norb,norbitals
LOGICAL               :: UNRES
call time_II_operations1

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
ndmat = 1
call initDFTdatatype(DFTDATA)
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
call mem_dft_alloc(Dmat,nbast,nbast,ndmat)
call mem_dft_alloc(DFTDATA%grad,3,natoms)
DFTDATA%grad = 0E0_realk

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
DFTDATA%natoms = natoms

call mem_dft_alloc(DFTDATA%orb2atom,nbast)
nOrbitals = 0
DO I=1,nAtoms
   select case(AORdefault)
   case(AORegular)
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbREG
   case(AOdfAux)
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbAUX
   case(AOdfCABS)
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbCABS
   case(AOdfJK)
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbJK
   case(AOVAL)
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbVAL
   case default
     CALL LSQUIT('Non-valid AORdefault in II_get_xc_geoderiv_molgrad',-1)
   end select
   DO J = nOrbitals+1,nOrbitals+nOrb
      DFTDATA%orb2atom(J)=I
   ENDDO
   nOrbitals = nOrbitals + nOrb
ENDDO
IF(nOrbitals .NE. nbast)&
     &CALL LSQUIT('mismatch in orbital dimension in II_get_xc_geoderiv_molgrad',-1)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('not working - not implemented',-1)
   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1E0_realk,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,1),1)
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,D%nrow,D%ncol,ndmat)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_dft_geoderiv_molgrad(setting,LUPRI,1,nbast,ndmat,DmatAO,DFTDATA,UNRES)
CALL LSTIMER('geoderiv_molgrd',TS,TE,LUPRI)
GRAD = DFTDATA%grad
CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_geoderiv_molgrad is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_geoderiv_molgrad is  ',WALLTIME,LUPRI)
call mem_dft_dealloc(DFTDATA%orb2atom)
call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%grad)
call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_geoderiv_molgrad)

END SUBROUTINE II_get_xc_geoderiv_molgrad

!> \brief Calculates the xc contribution to the linear response
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_get_xc_linrsp(LUPRI,LUERR,SETTING,nbast,b,D,G,nbmat)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> Number of b mat
INTEGER,intent(in)    :: nbmat
!> number of basisfunctions
INTEGER,intent(in)    :: nbast
!> The b matrix G(b)
TYPE(MATRIX),intent(in) :: b(nbmat)
!> The density matrix
TYPE(MATRIX),intent(in) :: D
!> The xc cont to the linear response
TYPE(MATRIX),intent(inout) :: G(nbmat)
!
INTEGER               :: i,j,ndmat
TYPE(MATRIX)          :: temp
TYPE(DFTDATATYPE)  :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ibmat
LOGICAL               :: UNRES
call time_II_operations1

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
CALL LS_GETTIM(CPU1,WALL1)
call initDFTdatatype(DFTDATA)
DFTDATA%nbast = nbast
ndmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
DFTDATA%nbmat = nbmat
call mem_dft_alloc(Dmat,nbast,nbast,ndmat)
DFTDATA%nfmat = nbmat
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,nbmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*nbmat)

!WARNING: the densitymatrix for an unres calc fullfill Tr(DS)= N   
!while the closed shell calculation fullfill Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
!   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)   
!   call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,2*nbmat)
!   CALL DCOPY(b%nrow*b%ncol,b%elms,1,DFTDATA%bmat(:,:,1),1)
!   CALL DCOPY(b%nrow*b%ncol,b%elmsb,1,DFTDATA%bmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1E0_realk,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,1),1)

   call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   IF(setting%IntegralTransformGC)THEN
      !we use DmatAO as a temporary array
      call mem_dft_alloc(DmatAO,nbast,nbast,nbmat)
      DO IBMAT=1,nbmat
         call mat_to_full(b(IBMAT),1E0_realk,DmatAO(:,:,IBMAT))
      ENDDO
      CALL GCAO2AO_transform_fullD(DmatAO,DFTDATA%BMAT,nbast,nbmat,setting,lupri)
      call mem_dft_dealloc(DmatAO)
   ELSE
      DO IBMAT=1,nbmat
         call mat_to_full(b(IBMAT),1E0_realk,DFTDATA%bmat(:,:,IBMAT))
      ENDDO
   ENDIF
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,D%nrow,D%ncol,ndmat)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_LINRSP(SETTING,LUPRI,1,nbast,ndmat,DmatAO,DFTDATA,UNRES)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,G%elms,1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,G%elmsb,1)
ELSE !CLOSED_SHELL
   DO IBMAT=1,nbmat
      CALL mat_set_from_full(DFTDATA%FKSM(:,:,IBMAT),1E0_realk,G(IBMAT),'XCmat')
      IF(setting%IntegralTransformGC)THEN
         call AO2GCAO_transform_matrixF(G(IBMAT),setting,lupri)
      ENDIF
      CALL mat_scal(0.5E0_realk,G(IBMAT))
   ENDDO
ENDIF

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_linrsp is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_linrsp is  ',WALLTIME,LUPRI)

call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%FKSM)
call mem_dft_dealloc(DFTDATA%BMAT)
call stats_dft_mem(lupri)

call time_II_operations2(JOB_II_get_xc_linrsp)
END SUBROUTINE II_get_xc_linrsp

!> \brief Calculates the xc contribution to the quadratic response
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_quadrsp(LUPRI,LUERR,SETTING,nbast,b,c,D,T)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The b matrix T(b,c)
TYPE(MATRIX)          :: b
!> The b matrix T(b,c)
TYPE(MATRIX)          :: c
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc cont to the quadratic response
TYPE(MATRIX)          :: T
!
INTEGER               :: i,j,ndmat,nbmat
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
LOGICAL               :: UNRES
call time_II_operations1

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
call initDFTdatatype(DFTDATA)
CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
ndmat = 1
nbmat = 2
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
DFTDATA%nbmat = nbmat
call mem_dft_alloc(Dmat,nbast,nbast,ndmat)
DFTDATA%nfmat = ndmat
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,ndmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*ndmat)

!WARNING: the densitymatrix for an unres calc fullfill Tr(DS)= N   
!while the closed shell calculation fullfill Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
!   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)   
!   call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,2*nbmat)
!   CALL DCOPY(b%nrow*b%ncol,b%elms,1,DFTDATA%bmat(:,:,1),1)
!   CALL DCOPY(b%nrow*b%ncol,b%elmsb,1,DFTDATA%bmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1E0_realk,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,1),1)
   call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   IF(setting%IntegralTransformGC)THEN
      !we use DmatAO as temporary array
      call mem_dft_alloc(DmatAO,nbast,nbast,nbmat)
      call mat_to_full(b,1E0_realk,DmatAO(:,:,1))
      call mat_to_full(c,1E0_realk,DmatAO(:,:,2))
      CALL GCAO2AO_transform_fullD(DmatAO,DFTDATA%BMAT,nbast,nbmat,setting,lupri)
      call mem_dft_dealloc(DmatAO)
   ELSE
      call mat_to_full(b,1E0_realk,DFTDATA%bmat(:,:,1))
      call mat_to_full(c,1E0_realk,DFTDATA%bmat(:,:,2))
   ENDIF
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,D%nrow,D%ncol,ndmat)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_QUADRSP(SETTING,LUPRI,1,nbast,ndmat,DmatAO,DFTDATA,UNRES)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,T%elms,1)
   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,T%elmsb,1)
   IF(setting%IntegralTransformGC)THEN
      call AO2GCAO_transform_matrixF(T,setting,lupri)
   ENDIF
ELSE !CLOSED_SHELL
   CALL mat_set_from_full(DFTDATA%FKSM(:,:,1),1E0_realk,T,'XCmat')
   IF(setting%IntegralTransformGC)THEN
      call AO2GCAO_transform_matrixF(T,setting,lupri)
   ENDIF
   CALL mat_scal(0.5E0_realk,T)
ENDIF

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_quadrsp is ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_quadrsp is ',WALLTIME,LUPRI)

call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%FKSM)
call mem_dft_dealloc(DFTDATA%BMAT)
call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_quadrsp)

END SUBROUTINE II_get_xc_quadrsp

!> \brief Calculates the xc contribution to the magnetic derivative kohn-sham matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,SETTING,nbast,D,F)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc cont to the magnetic deriv of F
TYPE(MATRIX)          :: F(3) !x,y and z components
INTEGER               :: i,j,ndmat
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
LOGICAL               :: UNRES
call time_II_operations1

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
call initDFTdatatype(DFTDATA)
CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
ndmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
call mem_dft_alloc(Dmat,nbast,nbast,ndmat)
DFTDATA%nfmat = 3*ndmat
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,3*ndmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*3*ndmat)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('II_get_xc_magderiv_kohnsham_mat not implemented for unrestricted yet',lupri)
!   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
!   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1E0_realk,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,1),1)
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,D%nrow,D%ncol,ndmat)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_dft_magderiv_kohnsham_mat(SETTING,LUPRI,1,nbast,ndmat,DmatAO,DFTDATA,UNRES)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

!WARNING: 
! For a closed shell molecule calculated using an unrestricted and a 
! closed shell aproach  would fullfill
! F(from closed shell) = F_alpha(from onres) + F_beta(from onres) 
! but for some reason this is not what they want in SCF-loop so we 
! multiply the unrestriced result with 2

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('II_get_xc_magderiv_kohnsham_mat not implemented for unrestricted yet',lupri)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,temp%elms,1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,temp%elmsb,1)
!   CALL mat_DAXPY(1.0E0_realk,temp,F)
ELSE !CLOSED_SHELL
   do I=1,3*ndmat
      CALL mat_set_from_full(DFTDATA%FKSM(:,:,I),1E0_realk,F(I),'XCmat')
      IF(setting%IntegralTransformGC)THEN
         call AO2GCAO_transform_matrixF(F(I),setting,lupri)
      ENDIF
      CALL mat_scal(0.5E0_realk,F(I))
   enddo
ENDIF

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_magderiv_kohnsham_mat is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_magderiv_kohnsham_mat is  ',WALLTIME,LUPRI)

call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%FKSM)
call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_magderiv_kohnsham_mat)

END SUBROUTINE II_get_xc_magderiv_kohnsham_mat

!> \brief Calculates the xc contribution to the mag derivative of the linear rsp
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_magderiv_linrsp(LUPRI,LUERR,SETTING,nbast,b,D,G,nbmat)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> the B mat, perturb density
TYPE(MATRIX)          :: b(nbmat)
!> The density matrix
TYPE(MATRIX)          :: D
!> The xc cont to the magnetic deriv of G(b)
TYPE(MATRIX)          :: G(3*nbmat)
!> number of B matrices
INTEGER               :: nbmat
INTEGER               :: i,j,ndmat
TYPE(MATRIX)          :: temp
TYPE(DFTDATATYPE)  :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ibmat!,ISYM
LOGICAL               :: UNRES!,ALLSYM
!INTEGER,external :: matfull_get_isym
call time_II_operations1

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
call initDFTdatatype(DFTDATA)
WRITE(lupri,*)'STARTING II_get_xc_magderiv_linrsp'
CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
ndmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
DFTDATA%nbmat = nbmat
call mem_dft_alloc(Dmat,nbast,nbast,ndmat)
DFTDATA%nfmat = 3*nbmat
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,3*nbmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*3*nbmat)
DFTDATA%dosympart = .TRUE. !should be a test
!NONE OF THE BMATS ARE SYMMETRIC SO THIS WILL ALWAYS BE TRUE
!ALLSYM = .TRUE.
!DO IBMAT=1,nBmat
!   !1 = symmetric, 2 = anti-symmetric, 3 = no symmetry, 4 = zero matrix
!   ISYM = matfull_get_isym(Dmat(:,:,IBMAT),nbast,nbast,1.0E-10_realk)   
!   IF(ISYM.NE.1)ALLSYM = .FALSE.
!ENDDO
!IF(ALLSYM)THEN
!   !all the Bmatrices are symmetric
!   !which means that the will be no symmetric 
!   !contribution to the output matrix
!   DFTDATA%dosympart=.FALSE.
!ENDIF

IF(DFTDATA%dosympart)THEN
   call mem_dft_alloc(DFTDATA%FKSMS,nbast,nbast,3*nbmat)
   CALL LS_DZERO(DFTDATA%FKSMS,nbast*nbast*3*nbmat)
ENDIF
!WARNING: the densitymatrix for an unres calc fullfill Tr(DS)= N   
!while the closed shell calculation fullfill Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
!   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)   
!   call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,2*nbmat)
!   CALL DCOPY(b%nrow*b%ncol,b%elms,1,DFTDATA%bmat(:,:,1),1)
!   CALL DCOPY(b%nrow*b%ncol,b%elmsb,1,DFTDATA%bmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1E0_realk,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,1),1)
   call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   IF(setting%IntegralTransformGC)THEN
      !we use DmatAO as a temporary array
      call mem_dft_alloc(DmatAO,nbast,nbast,nbmat)
      DO IBMAT=1,nbmat
         call mat_to_full(b(IBMAT),1E0_realk,DmatAO(:,:,IBMAT))
      ENDDO
      CALL GCAO2AO_transform_fullD(DmatAO,DFTDATA%BMAT,nbast,nbmat,setting,lupri)
      call mem_dft_dealloc(DmatAO)
   ELSE
      DO IBMAT=1,nbmat
         call mat_to_full(b(IBMAT),1E0_realk,DFTDATA%bmat(:,:,IBMAT))
      ENDDO
   ENDIF
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,D%nrow,D%ncol,ndmat)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_MAGDERIV_LINRSP(SETTING,LUPRI,1,nbast,ndmat,DmatAO,DFTDATA,UNRES)
CALL LSTIMER('xc-Fock',TS,TE,LUPRI)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,1),1,G%elms,1)
!   CALL DCOPY(D%nrow*D%ncol,DFTDATA%FKSM(:,:,2),1,G%elmsb,1)
ELSE !CLOSED_SHELL
   DO IBMAT=1,nbmat*3
      CALL mat_set_from_full(DFTDATA%FKSM(:,:,IBMAT),1E0_realk,G(IBMAT),'XCmat')
      CALL mat_scal(1E0_realk,G(IBMAT))
   ENDDO
   IF(DFTDATA%dosympart)THEN
      call mat_init(temp,nbast,nbast)
      DO IBMAT=1,nbmat*3
         CALL mat_set_from_full(DFTDATA%FKSMS(:,:,IBMAT),1E0_realk,temp,'XCmat')
         CALL mat_daxpy(0.5E0_realk,temp,G(IBMAT))
      ENDDO
      call mat_free(temp)
   ENDIF
   IF(setting%IntegralTransformGC)THEN
      DO IBMAT=1,nbmat*3
         call AO2GCAO_transform_matrixF(G(IBMAT),setting,lupri)
      ENDDO
   ENDIF
ENDIF

CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_linrsp is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_linrsp is  ',WALLTIME,LUPRI)

call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%FKSM)
IF(DFTDATA%dosympart)call mem_dft_dealloc(DFTDATA%FKSMS)
call mem_dft_dealloc(DFTDATA%BMAT)
call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_magderiv_linrsp)

END SUBROUTINE II_get_xc_magderiv_linrsp

!> \brief Calculates the xc contribution to the geo derivative of the KohnSham matrix
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_geoderiv_FxDgrad(LUPRI,LUERR,SETTING,nbast,D,b,grad,natoms)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> the B mat, perturb density
TYPE(MATRIX)          :: b
!> The xc cont to the geomderiv of F 
REAL(REALK)           :: GRAD(3,natoms)
!> number of atoms
INTEGER               :: natoms
!
INTEGER               :: i,j
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,norb,norbitals,nbmat,ibmat
LOGICAL               :: UNRES
call time_II_operations1

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
call initDFTdatatype(DFTDATA)

ndmat = 1
nbmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
DFTDATA%nbmat = nbmat
DFTDATA%nfmat = 0
call mem_dft_alloc(Dmat,nbast,nbast,ndmat)
call mem_dft_alloc(DFTDATA%grad,3,natoms)
DFTDATA%grad = 0E0_realk

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
DFTDATA%natoms = natoms

call mem_dft_alloc(DFTDATA%orb2atom,nbast)
nOrbitals = 0
DO I=1,nAtoms
   nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbREG
   DO J = nOrbitals+1,nOrbitals+nOrb
      DFTDATA%orb2atom(J)=I
   ENDDO
   nOrbitals = nOrbitals + nOrb
ENDDO
IF(nOrbitals .NE. nbast)&
     &CALL LSQUIT('mismatch in orbital dimension in II_get_xc_geoderiv_FxDgrad',-1)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('not working - not implemented',-1)
   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1E0_realk,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,1),1)
   call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   IF(setting%IntegralTransformGC)THEN
      !we use DmatAO as a temporary array
      call mem_dft_alloc(DmatAO,nbast,nbast,nbmat)
      call mat_to_full(b,1E0_realk,DmatAO(:,:,1))
      CALL GCAO2AO_transform_fullD(DmatAO,DFTDATA%BMAT,nbast,nbmat,setting,lupri)
      call mem_dft_dealloc(DmatAO)
   ELSE
      call mat_to_full(b,1E0_realk,DFTDATA%bmat(:,:,1))
   ENDIF
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,D%nrow,D%ncol,ndmat)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_geoderiv_kohnsham_mat(setting,LUPRI,1,nbast,1,DmatAO,DFTDATA,UNRES)
CALL LSTIMER('geoderiv_FxDgrd',TS,TE,LUPRI)
GRAD = DFTDATA%grad
CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_geoderiv_FxDgrad is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_geoderiv_FxDgrad is  ',WALLTIME,LUPRI)

call mem_dft_dealloc(DFTDATA%orb2atom)
call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%grad)
call mem_dft_dealloc(DFTDATA%BMAT)
call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_geoderiv_FxDgrad)

END SUBROUTINE II_get_xc_geoderiv_FxDgrad

!> \brief Calculates the xc contribution to the geo derivative of the linear response
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_get_xc_geoderiv_GxDgrad(LUPRI,LUERR,SETTING,nbast,D,a,b,grad,natoms)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the logical unit number for the error file
INTEGER                   :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER               :: nbast
!> The density matrix
TYPE(MATRIX)          :: D
!> the A mat, perturb density A*G(B)
TYPE(MATRIX)          :: a
!> the B mat, perturb density A*G(B)
TYPE(MATRIX)          :: b
!> The xc cont to the geomderiv of G 
REAL(REALK)           :: GRAD(3,natoms)
!> number of atoms
INTEGER               :: natoms
!
INTEGER               :: i,j
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: Dmat(:,:,:)
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,norb,norbitals,nbmat,ibmat
LOGICAL               :: UNRES
call time_II_operations1
IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
call initDFTdatatype(DFTDATA)

ndmat = 1
nbmat = 2
DFTDATA%nbmat = nbmat
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat
call mem_dft_alloc(Dmat,nbast,nbast,ndmat)
call mem_dft_alloc(DFTDATA%grad,3,natoms)
DFTDATA%grad = 0E0_realk

CALL LS_GETTIM(CPU1,WALL1)
DFTDATA%nbast = nbast
DFTDATA%natoms = natoms

call mem_dft_alloc(DFTDATA%orb2atom,nbast)
nOrbitals = 0
DO I=1,nAtoms
   nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbREG
   DO J = nOrbitals+1,nOrbitals+nOrb
      DFTDATA%orb2atom(J)=I
   ENDDO
   nOrbitals = nOrbitals + nOrb
ENDDO
IF(nOrbitals .NE. nbast)&
     &CALL LSQUIT('mismatch in orbital dimension in II_get_xc_geoderiv_FxDgrad',-1)

!WARNING: the densitymatrix for an unrestriced calculation fullfill 
! Tr(DS)= N   
!while the closed shell calculation fullfill
! Tr(DS)=N/2
!we therefore have to multiply the closed shell density with 2 in
!order to get a physical electron density $ \int \rho(r) dr = N $

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('not working - not implemented',-1)
   CALL DCOPY(D%nrow*D%ncol,D%elms,1,Dmat(:,:,1),1)
   CALL DCOPY(D%nrow*D%ncol,D%elmsb,1,Dmat(:,:,2),1)
ELSE !CLOSED_SHELL
   call mat_to_full(D,1E0_realk,Dmat(:,:,1))
   CALL DSCAL(nbast*nbast,2E0_realk,Dmat(:,:,1),1)
   call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
   IF(setting%IntegralTransformGC)THEN
      !we use DmatAO as temporary array
      call mem_dft_alloc(DmatAO,nbast,nbast,nbmat)
      call mat_to_full(a,1E0_realk,DmatAO(:,:,1))
      call mat_to_full(b,1E0_realk,DmatAO(:,:,2))
      CALL GCAO2AO_transform_fullD(DmatAO,DFTDATA%BMAT,nbast,nbmat,setting,lupri)
      call mem_dft_dealloc(DmatAO)
   ELSE
      call mat_to_full(a,1E0_realk,DFTDATA%bmat(:,:,1))
      call mat_to_full(b,1E0_realk,DFTDATA%bmat(:,:,2))
   ENDIF
ENDIF
IF(setting%IntegralTransformGC)THEN
   call mem_dft_alloc(DmatAO,D%nrow,D%ncol,ndmat)
   CALL GCAO2AO_transform_fullD(Dmat,DmatAO,nbast,ndmat,setting,lupri)
   call mem_dft_dealloc(Dmat)
ELSE
   DmatAO => Dmat
ENDIF

CALL LSTIMER('START',TS,TE,LUPRI)
CALL II_DFT_geoderiv_linrspgrad(setting,LUPRI,1,nbast,1,DmatAO,DFTDATA,UNRES)
CALL LSTIMER('geoderiv_FxDgrd',TS,TE,LUPRI)
GRAD = DFTDATA%grad
CALL LS_GETTIM(CPU2,WALL2)
CPUTIME = CPU2-CPU1
WALLTIME = WALL2-WALL1
CALL ls_TIMTXT('>>>  CPU  Time used II_get_xc_geoderiv_FxDgrad is  ',CPUTIME,LUPRI)
CALL ls_TIMTXT('>>>  WALL Time used II_get_xc_geoderiv_FxDgrad is  ',WALLTIME,LUPRI)

call mem_dft_dealloc(DFTDATA%orb2atom)
call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%grad)
call mem_dft_dealloc(DFTDATA%BMAT)
call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_geoderiv_GxDgrad)

END SUBROUTINE II_get_xc_geoderiv_GxDgrad

end module II_XC_interfaceModule
