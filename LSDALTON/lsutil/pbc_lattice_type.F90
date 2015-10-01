MODULE lattice_type
use precision
use matrix_module, only: matrix
use molecule_typetype, only: MOLECULEINFO
!use dft_type
Implicit none
INTEGER, public :: max_layer
INTEGER, public :: n_neighbour
INTEGER, public :: nearfield

integer, parameter :: MaxPBCOpTypes = 4

type, public ::rspace_lat_info_t
  logical :: is_active(1:3)
  real(realk) :: avec(1:3,1:3)      ! lattice vectors as *column* vectors
end type

type,public :: rspcopdata_t
  logical :: should_be_computed
  logical :: has_been_computed
  character(len=40) :: filename
  character(len=100)  :: basename
end type

type,public ::lvec_data_t
  real(realk) :: lat_coord(1:3)   ! lattice coordinates of a cell
  real(realk) :: std_coord(1:3)      ! standard coordinates
  REAL(realk),pointer :: d_vec(:),d_mat(:,:), elms(:)
  LOGICAL      :: is_redundant,g2_computed,f1_computed,J_computed,Kx_computed
  LOGICAL      :: ovl_computed,Vz_computed,dm_computed
  TYPE(Moleculeinfo) :: molecule    
  TYPE(matrix) :: oper(MaxPBCOpTypes)         !Not sure 
  !type(rspcopdata_t) :: opdat(MaxPBCOpTypes+1)
  ! other data?
  !type(matrix), pointer :: density(:)
  real(realk),pointer :: C(:,:)
  INTEGER(short)      :: maxGab
end type

type, public :: lvec_list_t
  type(rspace_lat_info_t) :: ldef   ! definition of real-space lattice
  integer :: num_entries,nf_entries
  integer :: fdim(3)
  type(lvec_data_t), pointer :: lvec(:)
  type(lvec_data_t), pointer :: nflvec(:)
  type(rspcopdata_t) :: opdat(MaxPBCOpTypes+5)
  logical            :: COMP_PBC
  logical            :: compare_elmnts
  logical            :: read_file
  logical            :: store_mats
  logical            :: setup_pbclatt
  logical            :: testcase
  character(len=30) ::  debugdensfile
  integer            :: max_layer
  integer            :: fc1,fc2,fc3 !how many layers included in fockmatrix
  integer            :: oneop1,oneop2,oneop3!layers included in 1-operators
  integer            :: col1,col2,col3!how many layers included in J-matrix
  integer            :: Kx1,Kx2,Kx3!how many layers included in K-matrix
  integer            :: lmax,Tlmax
  integer            :: num_its,num_store
  Real(realk)        :: error,intthr
  integer            :: nneighbour,nf,ndmat
  INTEGER(short)     :: realthr
  integer      :: nk1,nk2,nk3
  integer      :: num_its_densmat
  character(len=10)  :: wannier_direct
  character(len=100) :: basename
end type

type, public :: lattice_cell_info_t
  type(latticecell_atom_pos_t),pointer :: atom(:)
  type(matrix),pointer :: getmultipole(:)
  LOGICAL              :: is_defined = .false. 
end type

type, public :: latticecell_atom_pos_t
  real(realk) :: center(3)
end type

TYPE(lvec_list_t), save :: pbc_control

public :: pbc_control,MaxPBCOpTypes

private 
contains
subroutine lattice_type_dummy()
implicit none
end subroutine lattice_type_dummy

END MODULE lattice_type
