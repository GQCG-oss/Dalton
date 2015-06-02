module II_XC_interfaceModule
  use precision
  use TYPEDEFTYPE, only: lssetting
  use TYPEDEF, only: GCAO2AO_transform_fullD,AO2GCAO_transform_matrixF
  use Matrix_module, only: matrix
  use Matrix_Operations, only: mat_to_full, mat_init, mat_daxpy, &
       & mat_set_from_full, mat_assign, &
       & mat_free, mat_scal, mtype_unres_dense, matrix_type,&
       & mat_to_full3D
  use dft_memory_handling
!  use memory_handling
  use lstiming
  use IIDFTKSM
  use IIDFTD, only: II_DFT_DISP
  use DFT_type
  use GCtransMod
  private
  public :: II_get_xc_Fock_mat, II_get_xc_Fock_mat_full,&
       & II_get_AbsoluteValue_overlap, II_get_AbsoluteValue_overlapSame,&
       & II_get_xc_geoderiv_molgrad, II_get_xc_linrsp,&
       & II_get_xc_quadrsp, II_get_xc_magderiv_kohnsham_mat,&
       & II_get_xc_magderiv_linrsp, II_get_xc_geoderiv_FxDgrad,&
       & II_get_xc_geoderiv_GxDgrad, II_get_xc_energy
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
REAL(REALK),pointer   :: DmatAO(:,:,:),EDFT2(:)
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
CALL LSTIMER('START',TS,TE,LUPRI)
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
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
DFTDATA%nfmat = ndmat2
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,ndmat2)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*ndmat2)

call II_XC_TransformDmatToAOFull(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)

CALL II_DFT_KSM(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,EDFT2,UNRES)

call mem_dft_dealloc(DmatAO)

IF(SETTING%SCHEME%DFT%DODISP) THEN 
    ! add empirical dispersion correction \Andreas Krapp
   CALL II_DFT_DISP(SETTING,DUMMY,1,1,0,LUPRI)

   do idmat = 1,ndmat2
      EDFT2(idmat) = EDFT2(idmat) + SETTING%EDISP
   enddo
ENDIF

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

CALL LSTIMER('II_get_xc_Fock_mat',TS,TE,LUPRI)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_Fock_mat)

END SUBROUTINE II_get_xc_Fock_mat_array

SUBROUTINE II_get_xc_Fock_mat_full(LUPRI,LUERR,SETTING,nbast,DmatAO,FmatAO)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in)    :: LUPRI
!> the logical unit number for the error file
INTEGER,intent(in)    :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER,intent(in)    :: nbast
!!> number of Densitymatrices
!INTEGER,intent(in)    :: ndmat
!> The density matrix
real(realk),intent(in) :: DmatAO(nbast,nbast)
!> The Kohn-Sham matrix
real(realk),intent(inout) :: FmatAO(nbast,nbast)
!
INTEGER               :: i,j,ndmat2,idmat
TYPE(DFTDATATYPE)     :: DFTDATA
REAL(REALK),pointer   :: EDFT2(:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2,CPUTIME,WALLTIME,DUMMY(1,1)
LOGICAL               :: UNRES
call time_II_operations1
ndmat2 = 1
UNRES=.FALSE.
IF(SETTING%SCHEME%DFT%DODISP) THEN 
   call lsquit('DODISP not implemented II_get_xc_Fock_mat_full',-1)
ENDIF
IF(setting%IntegralTransformGC)THEN
   call lsquit('IntegralTransformGC not implemented for II_get_xc_Fock_mat_full',-1)
ENDIF
call init_dftmemvar
CALL LSTIMER('START',TS,TE,LUPRI)
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
call mem_dft_alloc(EDFT2,ndmat2)
DFTDATA%ndmat = ndmat2
DFTDATA%nfmat = ndmat2
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,ndmat2)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*ndmat2)
CALL II_DFT_KSM(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,EDFT2,UNRES)
call mem_dft_dealloc(EDFT2)
call DAXPY(nbast*nbast*ndmat2,0.5E0_realk,DFTDATA%FKSM,1,FmatAO,1)
call mem_dft_dealloc(DFTDATA%FKSM)
CALL LSTIMER('II_get_xc_Fock_mat_full',TS,TE,LUPRI)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_Fock_mat)
END SUBROUTINE II_get_xc_Fock_mat_full

SUBROUTINE II_get_AbsoluteValue_overlap(LUPRI,LUERR,SETTING,nbast,nmo1,nmo2,Cmat1,Cmat2,S)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in)    :: LUPRI
!> the logical unit number for the error file
INTEGER,intent(in)    :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER,intent(in)    :: nbast
!> number of occupied basisfunctions
INTEGER,intent(in)    :: nmo1,nmo2
!> The MO coef matrix for 1. 
real(realk),intent(in) :: Cmat1(nbast,nmo1)
!> The MO coef matrix for 2.
real(realk),intent(in) :: Cmat2(nbast,nmo2)
!> The Absolute Valued overlap  matrix
real(realk),intent(inout) :: S(nmo1,nmo2)
!
REAL(REALK)           :: TS,TE
REAL(REALK),pointer   :: CMAT1AO(:,:),CMAT2AO(:,:)
LOGICAL               :: UNRES,SameCmat
!call time_II_operations1
UNRES=.FALSE.
IF(matrix_type .EQ. mtype_unres_dense)UNRES=.TRUE.
call init_dftmemvar
CALL LSTIMER('START',TS,TE,LUPRI)
CALL LS_DZERO(S,nmo1*nmo2)

!chose ABSVAL grid
SETTING%scheme%DFT%igrid = Grid_ABSVAL

!default for fine
SETTING%scheme%DFT%GridObject(Grid_ABSVAL)%RADINT = 2.15443E-17_realk
SETTING%scheme%DFT%GridObject(Grid_ABSVAL)%ANGINT = 47
SameCmat = .FALSE. 

IF(setting%IntegralTransformGC)THEN
!   call lsquit('IntegralTransformGC must be false in II_get_AbsoluteValue_overlap',-1)
   call mem_dft_alloc(Cmat1AO,nbast,nmo1)
!   Cmat1AO = Cmat1
   CALL DCOPY(nbast*nMO1,Cmat1,1,Cmat1AO,1)
   !Cmat1AO(nbast,nmo1) = CAO2GCAO(nbast,nbast)*Cmat1AO(nbast,nmo1)
   call GCAO2AO_half_transform_matrixFull(Cmat1AO,nbast,nmo1,setting,lupri,1)
   call mem_dft_alloc(Cmat2AO,nbast,nmo2)
!   Cmat2AO = Cmat2
   CALL DCOPY(nbast*nMO2,Cmat2,1,Cmat2AO,1)
   !Cmat1AO(nbast,nmo1) = CAO2GCAO(nbast,nbast)*Cmat1AO(nbast,nmo1)
   call GCAO2AO_half_transform_matrixFull(Cmat2AO,nbast,nmo2,setting,lupri,1)
   CALL II_DFT_ABSVAL_OVERLAP(SETTING,LUPRI,1,nbast,nmo1,nmo2,CMAT1AO,CMAT2AO,S,SameCmat)
   call mem_dft_dealloc(Cmat1AO)
   call mem_dft_dealloc(Cmat2AO)
ELSE
   CALL II_DFT_ABSVAL_OVERLAP(SETTING,LUPRI,1,nbast,nmo1,nmo2,CMAT1,CMAT2,S,SameCmat)
ENDIF


!revert to default grid
SETTING%scheme%DFT%igrid = Grid_Default
CALL LSTIMER('II_get_AbsoluteValue        ',TS,TE,LUPRI)
call stats_dft_mem(lupri)
!call time_II_operations2(JOB_II_get_xc_Fock_mat)
END SUBROUTINE II_get_AbsoluteValue_overlap

SUBROUTINE II_get_AbsoluteValue_overlapSame(LUPRI,LUERR,SETTING,nbast,nmo,Cmat,S)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in)    :: LUPRI
!> the logical unit number for the error file
INTEGER,intent(in)    :: LUERR
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!> number of basisfunctions
INTEGER,intent(in)    :: nbast
!> number of occupied basisfunctions
INTEGER,intent(in)    :: nmo
!> The MO coef matrix for 1. 
real(realk),intent(in) :: Cmat(nbast,nmo)
!> The Absolute Valued overlap  matrix
real(realk),intent(inout) :: S(nmo,nmo)
!
REAL(REALK)           :: TS,TE
REAL(REALK),pointer   :: CMAT1AO(:,:)
LOGICAL               :: UNRES,SameCmat
!call time_II_operations1
UNRES=.FALSE.
IF(matrix_type .EQ. mtype_unres_dense)UNRES=.TRUE.
call init_dftmemvar
CALL LSTIMER('START',TS,TE,LUPRI)
CALL LS_DZERO(S,nmo*nmo)

!chose ABSVAL grid
SETTING%scheme%DFT%igrid = Grid_ABSVAL

!default for fine
SETTING%scheme%DFT%GridObject(Grid_ABSVAL)%RADINT = 2.15443E-17_realk
SETTING%scheme%DFT%GridObject(Grid_ABSVAL)%ANGINT = 47
SameCmat = .TRUE. 

IF(setting%IntegralTransformGC)THEN
!   call lsquit('IntegralTransformGC must be false in II_get_AbsoluteValue_overlapsame',-1)
   call mem_dft_alloc(Cmat1AO,nbast,nmo)
!   Cmat1AO = Cmat
   CALL DCOPY(nbast*nMO,Cmat,1,Cmat1AO,1)
   !Cmat1AO(nbast,nmo1) = CAO2GCAO(nbast,nbast)*Cmat1AO(nbast,nmo1)
   call GCAO2AO_half_transform_matrixFull(Cmat1AO,nbast,nmo,setting,lupri,1)
   CALL II_DFT_ABSVAL_OVERLAP(SETTING,LUPRI,1,nbast,nmo,nmo,CMAT1AO,CMAT1AO,S,SameCmat)
   call mem_dft_dealloc(Cmat1AO)
ELSE
   CALL II_DFT_ABSVAL_OVERLAP(SETTING,LUPRI,1,nbast,nmo,nmo,CMAT,CMAT,S,SameCmat)
ENDIF

!revert to default grid
SETTING%scheme%DFT%igrid = Grid_Default
CALL LSTIMER('II_get_AbsoluteValueSame    ',TS,TE,LUPRI)
call stats_dft_mem(lupri)
!call time_II_operations2(JOB_II_get_xc_Fock_mat)
END SUBROUTINE II_get_AbsoluteValue_overlapSame

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
REAL(REALK),pointer   :: DmatAO(:,:,:),EDFT2(:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2,CPUTIME,WALLTIME
REAL(REALK)   :: DUMMY(1,1),RHOTHR,DFTHRI
LOGICAL               :: UNRES

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
CALL LSTIMER('START',TS,TE,LUPRI)
call initDFTdatatype(DFTDATA)
DFTDATA%nbast = nbast
ndmat2=ndmat
IF(UNRES)ndmat2=2*ndmat
call mem_dft_alloc(EDFT2,ndmat2)
DFTDATA%ndmat = ndmat2
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
!reduced accuracy !!
!SETTING%SCHEME%DFT%RHOTHR = SETTING%SCHEME%DFT%RHOTHR*100.0E0_realk

call II_XC_TransformDmatToAOFull(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)

!save parameters
DFTHRI = SETTING%scheme%DFT%DFTHRI 
RHOTHR = SETTING%scheme%DFT%RHOTHR 
!chose new loosend parameters
SETTING%scheme%DFT%DFTHRI = DFTHRI*5000
SETTING%scheme%DFT%RHOTHR = RHOTHR*5000

CALL II_DFT_KSME(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,EDFT2,UNRES)

call mem_dft_dealloc(DmatAO)

!revert to default grid
SETTING%scheme%DFT%DFTHRI = DFTHRI
SETTING%scheme%DFT%RHOTHR = RHOTHR

IF(SETTING%SCHEME%DFT%DODISP) THEN 
   ! add empirical dispersion correction
   CALL II_DFT_DISP(SETTING,DUMMY,1,1,0,LUPRI)
   do idmat = 1,ndmat2
      EDFT2(idmat) = EDFT2(idmat) + SETTING%EDISP
   enddo
ENDIF

IF(UNRES)THEN
   do idmat = 1,ndmat
      EDFT(idmat) = EDFT2(1+(idmat-1)*2)
   enddo
ELSE
   EDFT = EDFT2
ENDIF
call mem_dft_dealloc(EDFT2)

CALL LSTIMER('II_get_xc_energy',TS,TE,LUPRI)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)

END SUBROUTINE II_get_xc_energy

!> \brief Calculates the xc contribution to the molecular gradient
!> \author T. Kjaergaard
!> \date 2008
SUBROUTINE II_get_xc_geoderiv_molgrad(LUPRI,LUERR,SETTING,nbast,D,grad,natoms)
use LSparameters
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
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,norb,norbitals,ndmat2
LOGICAL               :: UNRES
call time_II_operations1

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
ndmat = 1
ndmat2 = 1
call initDFTdatatype(DFTDATA)
IF(matrix_type .EQ. mtype_unres_dense)ndmat2=2
DFTDATA%ndmat = ndmat2
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
call mem_dft_alloc(DFTDATA%grad,3,natoms)
CALL LS_DZERO(DFTDATA%grad,3*natoms)

CALL LSTIMER('START',TS,TE,LUPRI)
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
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbCABS &
          & + setting%molecule(1)%p%ATOM(I)%nContOrbREG
   case(AOdfCABO)!CABS only 
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbCABS
   case(AOdfJK)
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbJK
   case(AOVAL)
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbVAL
   case(AOadmm)
     nOrb = setting%molecule(1)%p%ATOM(I)%nContOrbADMM
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

call II_XC_TransformDmatToAOFull_single(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)

CALL II_dft_geoderiv_molgrad(setting,LUPRI,1,nbast,ndmat,DmatAO,DFTDATA,UNRES)
CALL DCOPY(3*natoms,DFTDATA%grad,1,GRAD,1)
CALL LSTIMER('II_get_xc_geoderiv_molgrad',TS,TE,LUPRI)
call mem_dft_dealloc(DFTDATA%orb2atom)
call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%grad)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)
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
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE
REAL(REALK)           :: CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ibmat,ndmat2
LOGICAL               :: UNRES
call time_II_operations1

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF
call init_dftmemvar
CALL LSTIMER('START',TS,TE,LUPRI)
call initDFTdatatype(DFTDATA)
DFTDATA%nbast = nbast
ndmat = 1
ndmat2 = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat2=2
DFTDATA%ndmat = ndmat2
DFTDATA%nbmat = nbmat
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
DFTDATA%nfmat = nbmat
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,nbmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*nbmat)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
ELSE
   call II_XC_TransformDmatToAOFull_single(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)
ENDIF

call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
call II_XC_TransformDmatToAOFull(b,nbmat,nbast,DFTDATA%BMAT,nbmat,.FALSE.,setting,lupri)

CALL II_DFT_LINRSP(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,UNRES)

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

CALL LSTIMER('II_get_xc_linrsp',TS,TE,LUPRI)
call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%FKSM)
call mem_dft_dealloc(DFTDATA%BMAT)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)

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
INTEGER               :: i,j,ndmat,ndmat2,nbmat
TYPE(DFTDATATYPE)     :: DFTDATA
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
CALL LSTIMER('START',TS,TE,LUPRI)
DFTDATA%nbast = nbast
ndmat = 1
ndmat2 = 1
nbmat = 2
IF(matrix_type .EQ. mtype_unres_dense)ndmat=2
DFTDATA%ndmat = ndmat2
DFTDATA%nbmat = nbmat
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat)
DFTDATA%nfmat = ndmat
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,ndmat)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*ndmat)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
ELSE
   call II_XC_TransformDmatToAOFull_single(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)
ENDIF

call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
call II_XC_TransformDmatToAOFull_bc(b,c,nbast,DFTDATA%BMAT,nbmat,.FALSE.,setting,lupri)

CALL II_DFT_QUADRSP(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,UNRES)

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

CALL LSTIMER('II_get_xc_quadrsp',TS,TE,LUPRI)

call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%FKSM)
call mem_dft_dealloc(DFTDATA%BMAT)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)
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
INTEGER               :: i,j,ndmat,ndmat2
TYPE(DFTDATATYPE)     :: DFTDATA
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
CALL LSTIMER('START',TS,TE,LUPRI)
DFTDATA%nbast = nbast
ndmat = 1
ndmat2 = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat2=2
DFTDATA%ndmat = ndmat2
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
DFTDATA%nfmat = 3*ndmat2
call mem_dft_alloc(DFTDATA%FKSM,nbast,nbast,3*ndmat2)
CALL LS_DZERO(DFTDATA%FKSM,nbast*nbast*3*ndmat2)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('II_get_xc_magderiv_kohnsham_mat not implemented for unrestricted yet',lupri)
ELSE
   call II_XC_TransformDmatToAOFull_single(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)
ENDIF

CALL II_dft_magderiv_kohnsham_mat(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,UNRES)

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
   do I=1,3*ndmat2
      CALL mat_set_from_full(DFTDATA%FKSM(:,:,I),1E0_realk,F(I),'XCmat')
      IF(setting%IntegralTransformGC)THEN
         call AO2GCAO_transform_matrixF(F(I),setting,lupri)
      ENDIF
      CALL mat_scal(0.5E0_realk,F(I))
   enddo
ENDIF

CALL LSTIMER('II_get_xc_magderiv_kohnsham',TS,TE,LUPRI)

call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%FKSM)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)
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
INTEGER               :: i,j,ndmat,ndmat2
TYPE(MATRIX)          :: temp
TYPE(DFTDATATYPE)  :: DFTDATA
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
CALL LSTIMER('START',TS,TE,LUPRI)
DFTDATA%nbast = nbast
ndmat = 1
ndmat2 = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat2=2
DFTDATA%ndmat = ndmat2
DFTDATA%nbmat = nbmat
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
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

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('NOT IMPLEMENTED YET',-1)
ELSE
   call II_XC_TransformDmatToAOFull_single(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)
ENDIF

call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
call II_XC_TransformDmatToAOFull(b,nbmat,nbast,DFTDATA%BMAT,nbmat,.FALSE.,setting,lupri)

CALL II_DFT_MAGDERIV_LINRSP(SETTING,LUPRI,1,nbast,ndmat2,DmatAO,DFTDATA,UNRES)

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

CALL LSTIMER('II_get_xc_magderiv_linrsp',TS,TE,LUPRI)

call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%FKSM)
IF(DFTDATA%dosympart)call mem_dft_dealloc(DFTDATA%FKSMS)
call mem_dft_dealloc(DFTDATA%BMAT)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)
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
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,ndmat2,norb,norbitals,nbmat,ibmat
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
ndmat2 = 1
nbmat = 1
IF(matrix_type .EQ. mtype_unres_dense)ndmat2=2
DFTDATA%ndmat = ndmat2
DFTDATA%nbmat = nbmat
DFTDATA%nfmat = 0
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
call mem_dft_alloc(DFTDATA%grad,3,natoms)
CALL LS_DZERO(DFTDATA%grad,3*natoms)

CALL LSTIMER('START',TS,TE,LUPRI)
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

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('not working - not implemented',-1)
ELSE
   call II_XC_TransformDmatToAOFull_single(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)
ENDIF
call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
call II_XC_TransformDmatToAOFull_single(B,nbmat,nbast,DFTDATA%BMAT,nbmat,.FALSE.,setting,lupri)

CALL II_DFT_geoderiv_kohnsham_mat(setting,LUPRI,1,nbast,1,DmatAO,DFTDATA,UNRES)
CALL DCOPY(3*natoms,DFTDATA%grad,1,GRAD,1)
CALL LSTIMER('II_get_xc_geoderiv_FxDgrad',TS,TE,LUPRI)

call mem_dft_dealloc(DFTDATA%orb2atom)
call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%grad)
call mem_dft_dealloc(DFTDATA%BMAT)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)
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
REAL(REALK),pointer   :: DmatAO(:,:,:)
REAL(REALK)           :: TS,TE,CPU1,CPU2,WALL1,WALL2
REAL(REALK)           :: CPUTIME,WALLTIME
INTEGER               :: ndmat,ndmat2,norb,norbitals,nbmat,ibmat
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
ndmat2 = 1
nbmat = 2
DFTDATA%nbmat = nbmat
IF(matrix_type .EQ. mtype_unres_dense)ndmat2=2
DFTDATA%ndmat = ndmat2
call mem_dft_alloc(DmatAO,nbast,nbast,ndmat2)
call mem_dft_alloc(DFTDATA%grad,3,natoms)
CALL LS_DZERO(DFTDATA%grad,3*natoms)

CALL LSTIMER('START',TS,TE,LUPRI)
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

IF(matrix_type .EQ. mtype_unres_dense)THEN
   CALL LSQUIT('not working - not implemented',-1)
ELSE
   call II_XC_TransformDmatToAOFull_single(D,ndmat,nbast,DmatAO,ndmat2,.TRUE.,setting,lupri)
ENDIF
call mem_dft_alloc(DFTDATA%BMAT,nbast,nbast,nbmat)
call II_XC_TransformDmatToAOFull_bc(A,B,nbast,DFTDATA%BMAT,nbmat,.FALSE.,setting,lupri)

CALL II_DFT_geoderiv_linrspgrad(setting,LUPRI,1,nbast,1,DmatAO,DFTDATA,UNRES)
CALL DCOPY(3*natoms,DFTDATA%grad,1,GRAD,1)
CALL LSTIMER('II_get_xc_geoderiv_GxDgrad',TS,TE,LUPRI)

call mem_dft_dealloc(DFTDATA%orb2atom)
call mem_dft_dealloc(DmatAO)
call mem_dft_dealloc(DFTDATA%grad)
call mem_dft_dealloc(DFTDATA%BMAT)
IF (setting%scheme%intprint.GT.0.OR.PrintDFTmem) call stats_dft_mem(lupri)
call time_II_operations2(JOB_II_get_xc_geoderiv_GxDgrad)

END SUBROUTINE II_get_xc_geoderiv_GxDgrad

subroutine II_XC_TransformDmatToAOFull(D,ndmat,nbast,Dmat,ndmat2,Dens2scal,setting,lupri)
implicit none
!> logical unit number for printing
INTEGER,intent(in)    :: lupri
!> number of basis functions
INTEGER,intent(in)    :: nbast
!> number of Densitymatrices
INTEGER,intent(in)    :: ndmat,ndmat2
!> The density matrix
TYPE(MATRIX),intent(in) :: D(ndmat)
!> The full density matrix
real(realk) :: Dmat(nbast,nbast,ndmat2)
logical,intent(in)    :: Dens2scal
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!
LOGICAL               :: UNRES
TYPE(MATRIX) :: DAO
integer :: idmat
real(realk) :: Factor

IF(Dens2scal)THEN
   Factor = 2.0E0_realk
ELSE
   Factor = 1.0E0_realk
ENDIF

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF

IF(setting%IntegralTransformGC)THEN
   call mat_init(DAO,nbast,nbast)
   do idmat=1,ndmat
      call GCAO2AO_transform_matrixD2(D(idmat),DAO,setting,lupri)
      IF(UNRES)THEN
         call mat_to_full3D(DAO, 1.0E0_realk, Dmat,nbast,nbast,ndmat2,1+(idmat-1)*2,2+(idmat-1)*2)
      ELSE !CLOSED_SHELL
         call mat_to_full3D(DAO, Factor, Dmat,nbast,nbast,ndmat2,idmat,idmat)
      endif
   enddo
   call mat_free(DAO)
ELSE
   IF(UNRES)THEN
      do idmat = 1,ndmat
         call mat_to_full3D(D(idmat), 1.0E0_realk, Dmat,nbast,nbast,ndmat2,1+(idmat-1)*2,2+(idmat-1)*2)
      enddo
   ELSE !CLOSED_SHELL
      do idmat = 1,ndmat
         call mat_to_full3D(D(idmat), Factor, Dmat,nbast,nbast,ndmat2,idmat,idmat)
      enddo
   ENDIF
ENDIF
end subroutine II_XC_TransformDmatToAOFull

subroutine II_XC_TransformDmatToAOFull_bc(b,c,nbast,Dmat,ndmat2,Dens2scal,setting,lupri)
implicit none
!> logical unit number for printing
INTEGER,intent(in)    :: lupri
!> number of basis functions
INTEGER,intent(in)    :: nbast
!> number of Densitymatrices
INTEGER,intent(in)    :: ndmat2
!> The Bmat matrix
TYPE(MATRIX),intent(in) :: B
!> The Cmat matrix
TYPE(MATRIX),intent(in) :: C
!> The full density matrix
real(realk) :: Dmat(nbast,nbast,ndmat2)
logical,intent(in)    :: dens2scal
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!
LOGICAL               :: UNRES
TYPE(MATRIX) :: DAO
integer :: idmat
real(realk) :: Factor

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
   print*,'ndmat2',ndmat2
   IF(ndmat2.NE.4)call lsquit('Error II_XC_TransformDmatToAOFull_bc unres',-1)
ELSE
   UNRES=.FALSE.
   print*,'ndmat2',ndmat2
   IF(ndmat2.NE.2)call lsquit('Error II_XC_TransformDmatToAOFull_bc',-1)
ENDIF

IF(Dens2scal)THEN
   Factor = 2.0E0_realk
ELSE
   Factor = 1.0E0_realk
ENDIF

IF(setting%IntegralTransformGC)THEN
   call mat_init(DAO,nbast,nbast)
   !B vec
   call GCAO2AO_transform_matrixD2(B,DAO,setting,lupri)
   IF(UNRES)THEN
      call mat_to_full3D(DAO, 1.0E0_realk, Dmat,nbast,nbast,ndmat2,1,2)
   ELSE !CLOSED_SHELL
      call mat_to_full3D(DAO, Factor, Dmat,nbast,nbast,ndmat2,1,1)
   endif
   !C vec
   call GCAO2AO_transform_matrixD2(C,DAO,setting,lupri)
   IF(UNRES)THEN
      call mat_to_full3D(DAO, 1.0E0_realk, Dmat,nbast,nbast,ndmat2,3,4)
   ELSE !CLOSED_SHELL
      call mat_to_full3D(DAO, Factor, Dmat,nbast,nbast,ndmat2,2,2)
   endif
   call mat_free(DAO)
ELSE
   IF(UNRES)THEN
      call mat_to_full3D(B, 1.0E0_realk, Dmat,nbast,nbast,ndmat2,1,2)
      call mat_to_full3D(C, 1.0E0_realk, Dmat,nbast,nbast,ndmat2,3,4)
   ELSE !CLOSED_SHELL
      call mat_to_full3D(B, Factor, Dmat,nbast,nbast,ndmat2,1,1)
      call mat_to_full3D(C, Factor, Dmat,nbast,nbast,ndmat2,2,2)
   ENDIF
ENDIF
end subroutine II_XC_TransformDmatToAOFull_bc


subroutine II_XC_TransformDmatToAOFull_single(D,ndmat,nbast,Dmat,ndmat2,Dens2scal,setting,lupri)
implicit none
!> logical unit number for printing
INTEGER,intent(in)    :: lupri
!> number of basis functions
INTEGER,intent(in)    :: nbast
!> number of Densitymatrices
INTEGER,intent(in)    :: ndmat,ndmat2
!> The density matrix
TYPE(MATRIX),intent(in) :: D
!> The full density matrix
real(realk) :: Dmat(nbast,nbast,ndmat2)
logical,intent(in)    :: dens2scal
!> info about molecule,basis and dft parameters
TYPE(LSSETTING)       :: SETTING
!
LOGICAL               :: UNRES
TYPE(MATRIX) :: DAO
integer :: idmat
real(realk) :: Factor
IF(ndmat.NE.1)call lsquit('error in II_XC_TransformDmatToAOFull_single',-1)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   UNRES=.TRUE.
ELSE
   UNRES=.FALSE.
ENDIF

IF(Dens2scal)THEN
   Factor = 2.0E0_realk
ELSE
   Factor = 1.0E0_realk
ENDIF

IF(setting%IntegralTransformGC)THEN
   call mat_init(DAO,nbast,nbast)
   call GCAO2AO_transform_matrixD2(D,DAO,setting,lupri)
   IF(UNRES)THEN
      call mat_to_full3D(DAO, 1.0E0_realk, Dmat,nbast,nbast,ndmat2,1,2)
   ELSE !CLOSED_SHELL
      call mat_to_full3D(DAO, Factor, Dmat,nbast,nbast,ndmat2,1,1)
   endif
   call mat_free(DAO)
ELSE
   IF(UNRES)THEN
      call mat_to_full3D(D, 1.0E0_realk,Dmat,nbast,nbast,ndmat2,1,2)
   ELSE !CLOSED_SHELL
      call mat_to_full3D(D, Factor, Dmat,nbast,nbast,ndmat2,1,1)
   ENDIF
ENDIF
end subroutine II_XC_TransformDmatToAOFull_single

end module II_XC_interfaceModule
