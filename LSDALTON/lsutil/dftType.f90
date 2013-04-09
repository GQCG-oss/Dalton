!> dft type module
!> \author T.Kjaergaard
!> \date 2010-02-21
MODULE dft_typetype
use precision
integer,save      :: GRIDITERATIONS
! choose reasonably large. Exceeding this limit means that boxes are too large.
Integer, parameter :: MXBLLEN=128
Integer,save :: BoxMemRequirement
! Structure for different dft grid (example in II_get_admm_exchange_mat)
TYPE dft_grid
  logical                :: isset
  CHARACTER(len=80)      :: dftfunc
  character(21)          :: filename
  integer                :: GRIDITERATIONS
  integer                :: nbuflen
  integer                :: maxNactBAST
  TYPE(dft_grid),pointer :: next
END TYPE dft_grid

!> Keywords and input to gridgeneration and exchange-correlation calculation
TYPE DFTparam
INTEGER           :: RADIALGRID !(1 = GC2, 2 = LMG, 3 = TURBO)
INTEGER           :: PARTITIONING !integer in following list
!(1=SSF, 2=Becke, 3=Becke-original, 4=block, 5=blockssf, 6=cartesian)
LOGICAL           :: ZdependenMaxAng !(default 0)

INTEGER           :: GRDONE !IF GRID HAS BEEN CREATED 
REAL(REALK)       :: RADINT
INTEGER           :: ANGMIN
INTEGER           :: ANGINT
INTEGER           :: HRDNES ! hardness of the partition function in the becke schemes
!
REAL(REALK)       :: DFTELS !       DEFAULT: 1.0E-3_realk
REAL(REALK)       :: DFTHR0 !                1.0E-9_realk not used
REAL(REALK)       :: DFTHRI !                2.0E-12_realk
REAL(REALK)       :: DFTHRL !                1.0E-10_realk
REAL(REALK)       :: RHOTHR !                2.0E-15_realk
LOGICAL           :: NOPRUN !                .FALSE.
LOGICAL           :: DFTASC !                .FALSE.
LOGICAL           :: DFTPOT !                .FALSE.
LOGICAL           :: DODISP !                .FALSE. ; empirical dispersion correction following Grimme
REAL(REALK)       :: DFTIPT !                1.0E-20_realk      
REAL(REALK)       :: DFTBR1 !                1.0E-20_realk
REAL(REALK)       :: DFTBR2 !                1.0E-20_realk
LOGICAL           :: DFTADD !                .TRUE.
LOGICAL           :: DISPDONE !              .FALSE.
INTEGER           :: TURBO
INTEGER           :: NBUFLEN                 !0
INTEGER           :: maxNactBAST             !0
LOGICAL           :: NEWGRID                 !.FALSE.
Logical           :: testNelectrons          !.TRUE.
Logical           :: LB94      !van Leeuwen-Baerends correction
Logical           :: CS00      !Casida-Salahub asymptotic correction
REAL(REALK)       :: CS00shift !Casida-Salahub shift - if 0 use Zhan-Nichols-Dixon shift
REAL(REALK)       :: CS00eHOMO !energy of HOMO to use with Zhan-Nichols-Dixon shift
REAL(REALK)       :: CS00ZND1   !Zhan-Nichols-Dixon shift parameter
REAL(REALK)       :: CS00ZND2   !Zhan-Nichols-Dixon shift parameter
REAL(REALK)       :: HFexchangeFac
type(dft_grid)    :: L2GRID     !Grid parameters for level 2/ADMM grid
type(dft_grid)    :: L3GRID     !Grid parameters for level 3/regular grid
CHARACTER(len=80) :: dftfunc                 !""
LOGICAL           :: XCFUN                   !.FALSE.
END type DFTPARAM

!> USED IN II_DFTINT TO save data like Fockmatrix and gradients
!> calculated at each gridpoint 
TYPE DFTDATATYPE
INTEGER             :: nbast
INTEGER             :: nfmat !# Result matrices
INTEGER             :: ndmat !# Density matrices
INTEGER             :: nbmat !# response vectors
INTEGER             :: nWorkNactBastNblen
INTEGER             :: nWorkNactBast 
INTEGER             :: nWorkNactBastNactBast 
real(realk),pointer :: energy(:)     !ndmat
real(realk),pointer :: electrons(:)  !ndmat
real(realk),pointer :: BMAT(:,:,:)!nbast,nbast,nbmat
real(realk),pointer :: FKSM(:,:,:)!nbast,nbast,ndmat
real(realk),pointer :: FKSMS(:,:,:)!nbast,nbast,ndmat
logical             :: dosympart
INTEGER             :: natoms
INTEGER,pointer     :: orb2atom(:)
real(realk),pointer :: grad(:,:)
Logical             :: LB94      !van Leeuwen-Baerends correction
Logical             :: CS00      !Casida-Salahub asymptotic correction
REAL(REALK)         :: CS00shift !Casida-Salahub shift - if 0 use Zhan-Nichols-Dixon shift
REAL(REALK)         :: CS00eHOMO !energy of HOMO to use with Zhan-Nichols-Dixon shift
REAL(REALK)         :: CS00ZND1   !Zhan-Nichols-Dixon shift parameter
REAL(REALK)         :: CS00ZND2   !Zhan-Nichols-Dixon shift parameter
REAL(REALK)         :: HFexchangeFac
END TYPE DFTDATATYPE

CONTAINS
SUBROUTINE store_dft_grid(grid,filename,dft)
TYPE(dft_grid) :: grid
TYPE(dftparam) :: dft
character(21)  :: filename
grid%filename       = filename
grid%GRIDITERATIONS = GRIDITERATIONS
grid%dftfunc        = dft%dftfunc
grid%nbuflen        = dft%nbuflen
grid%maxNactBAST    = dft%maxNactBAST
END SUBROUTINE store_dft_grid

SUBROUTINE get_dft_grid(grid,filename,dft)
TYPE(dft_grid) :: grid
TYPE(dftparam) :: dft
character(21)  :: filename

IF (filename.EQ.grid%filename) THEN
  GRIDITERATIONS  = grid%GRIDITERATIONS
  dft%dftfunc     = grid%dftfunc       
  dft%nbuflen     = grid%nbuflen       
  dft%maxNactBAST = grid%maxNactBAST   
ELSE
  call lsquit('Filename error in get_dft_grid',-1)
ENDIF
END SUBROUTINE get_dft_grid

END MODULE dft_typetype
