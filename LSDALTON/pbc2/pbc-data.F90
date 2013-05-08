module pbc_data
  use precision
  implicit none

  !private
  public

  ! Unit number for output file writing
!  INTEGER, SAVE :: LUPRI = 2

  type latticegeom
     logical :: latvec_ok           ! indicates whether there are meaningful data in the lattice vector variables
     logical :: dim_is_active(3)    ! indicates whether a lattice dimension is active (i.e. the system is to replicated along this dimension)
     real(realk) :: latvec(3,3)     ! latvec(k,:) is the k:th lattice vector
     real(realk) :: reclatvec(3,3)  ! reciprocal lattice vectors
     real(realk) :: invlatvec(3,3)  ! the matrix inverse of latvec(:,:)
     real(realk) :: cell_volume     ! volume of the unit cell

     integer :: num_kpoints   ! number of k points to sample in BZ
     integer :: num_k1, num_k2, num_k3
     !real(realk) :: realspc_thres   ! defines when to neglect Fourier comp.

     !integer :: nfdebug     ! manual NF size for debugging

     integer :: qfict_update_cnt ! counts the number of updates of Qfict
     real(realk) :: qfict(4)           ! fictitious charges
     real(realk) :: qfict_pos_std(3,4) ! position (standard coord.)
  end type latticegeom

  type(latticegeom), save :: lat_data

  type occ_scheme_typ
     logical insulator_occ
     logical force_symmetric
  end type occ_scheme_typ

  type(occ_scheme_typ), save :: occ_scheme

  type scf_scheme_typ
     logical :: use_tc_dmat
     logical :: use_Cmax2_tol
     integer :: nproject
     real(realk) :: lindep_tol
     real(realk) :: Cmax2_tol
  end type scf_scheme_typ

  type(scf_scheme_typ), save :: scf_scheme

  type lindep_data_t
     logical :: do_projection
     logical :: use_nproject
     logical :: be_silent
     integer(realk) :: num_occ
     integer(realk) :: nproject
     real(realk) :: Seig_tol
     real(realk) :: Cmax2_tol
  end type lindep_data_t

  !============================================
  ! data structures to handle the k-point grid

  integer, parameter :: Max_bassiz = 1000
  integer, parameter :: Max_kpoints = 1500

  type BZpoint_t
     ! This is a data structure to represent a single k-point
     logical :: self_dual, is_gamma
     integer :: ix_orig
     integer :: n(3)
     integer :: ninv(3)
     real(realk) :: weight
     real(realk) :: lambda(3)
  end type BZpoint_t

  type BZgrid_t
     ! This is a data structure to represent a sampling grid in
     ! the first Brillouin zone
     logical :: use_invsym
     integer :: Nk_dim1, Nk_dim2, Nk_dim3, Nk
     integer :: Nk_nosym
     real(realk) :: reclvec(3,3)
     type(BZpoint_t) :: kpnt(Max_kpoints)
  end type BZgrid_t

  type splitBZgrid_t
     ! This data structure represents coarse-grained grids that
     ! when joined form the full fine-grained grid
     integer :: Nsplit(3)
     type(BZgrid_t), pointer :: subBZ(:,:,:)
     type(BZgrid_t), pointer :: fullBZ
  end type splitBZgrid_t

!contains

end module pbc_data
