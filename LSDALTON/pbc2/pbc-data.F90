MODULE pbc_data
  USE precision
  IMPLICIT NONE

  !PRIVATE
  PUBLIC

  TYPE latticegeom
	  !> indicates whether there are meaningful data in the lattice vectors
     LOGICAL :: latvec_ok           
	  !> indicates whether a lattice dimension is active (i.e. the system is to 
	  !> replicated along this dimension)
     LOGICAL :: dim_is_active(3)    
	  !> latvec(k,:) is the k:th lattice vector
     REAL(realk) :: latvec(3,3)     
	  !> reciprocal lattice vectors
     REAL(realk) :: reclatvec(3,3)  
	  !> the matrix inverse of latvec(:,:)
     REAL(realk) :: invlatvec(3,3)  
	  !> volume of the unit cell
     REAL(realk) :: cell_volume     
	  !> number of k points to sample in BZ
     INTEGER :: num_kpoints   
	  !> Number of k points in each spatial direction
     INTEGER :: num_k1, num_k2, num_k3
	  !> counts the number of updates of Qfict
     INTEGER :: qfict_update_cnt 
	  !> fictitious charges
     REAL(realk) :: qfict(4)           
	  !> position (standard coord.)
     REAL(realk) :: qfict_pos_std(3,4) 
  END TYPE latticegeom

  TYPE(latticegeom), SAVE :: lat_data

  TYPE occ_scheme_typ
     LOGICAL insulator_occ
     LOGICAL force_symmetric
  END TYPE occ_scheme_typ

  TYPE(occ_scheme_typ), SAVE :: occ_scheme

  TYPE scf_scheme_typ
     LOGICAL :: use_tc_dmat
     LOGICAL :: use_Cmax2_tol
     INTEGER :: nproject
     REAL(realk) :: lindep_tol
     REAL(realk) :: Cmax2_tol
  END TYPE scf_scheme_typ

  TYPE(scf_scheme_typ), SAVE :: scf_scheme

  TYPE lindep_data_t
     LOGICAL :: do_projection
     LOGICAL :: use_nproject
     LOGICAL :: be_silent
     INTEGER(realk) :: num_occ
     INTEGER(realk) :: nproject
     REAL(realk) :: Seig_tol
     REAL(realk) :: Cmax2_tol
  END TYPE lindep_data_t

  !============================================
  ! data structures to handle the k-point grid
  !============================================

  INTEGER, PARAMETER :: Max_bassiz = 1000
  INTEGER, PARAMETER :: Max_kpoints = 1500

  TYPE BZpoint_t
     !> This is a data structure to represent a single k-point.
     LOGICAL :: self_dual, is_gamma
     INTEGER :: ix_orig
     INTEGER :: n(3)
     INTEGER :: ninv(3)
     REAL(realk) :: weight
     REAL(realk) :: lambda(3)
  END TYPE BZpoint_t

  TYPE BZgrid_t
     !> This is a data structure to represent a sampling grid in
     !> the first Brillouin zone.
     LOGICAL :: use_invsym
     INTEGER :: Nk_dim1, Nk_dim2, Nk_dim3, Nk
     INTEGER :: Nk_nosym
     REAL(realk) :: reclvec(3,3)
     TYPE(BZpoint_t) :: kpnt(Max_kpoints)
  END TYPE BZgrid_t

  TYPE splitBZgrid_t
     !> This data structure represents coarse-grained grids that
     !> when joined form the full fine-grained grid.
     INTEGER :: Nsplit(3)
     TYPE(BZgrid_t), POINTER :: subBZ(:,:,:)
     TYPE(BZgrid_t), POINTER :: fullBZ
  END TYPE splitBZgrid_t

END MODULE pbc_data
