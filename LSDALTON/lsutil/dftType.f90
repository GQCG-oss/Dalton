!> dft type module
!> \author T.Kjaergaard
!> \date 2010-02-21
MODULE dft_typetype
use precision

! choose reasonably large. Exceeding this limit means that boxes are too large.
Integer, parameter :: MXBLLEN=128  
!DFT%NBUFLEN=1024 or maxNBUFLEN inside GenerateGrid is maximum number of 
!gridpoints in a box, these are read in and processed in chunks of MXBLLEN.

Integer,save :: BoxMemRequirement
!MIN(MXBLLEN,GridBox%npoints)*NactBast should be less then BoxMemRequirement
!otherwise the memory requirements are too big and the boxes are subdevided

!NactBast = number of active basis functions
!integers to determine the different grids in use
Integer,parameter :: Grid_Default = 1
Integer,parameter :: Grid_ADMML2 = 2
Integer,parameter :: Grid_ABSVAL = 3
!MPI node specific and grid specific 
integer,save      :: dft_GRIDITERATIONS(3)
INTEGER,save      :: dft_maxNactBAST(3)   

!integers to determine the different functionals in use
Integer,parameter :: dftfunc_Default = 1
Integer,parameter :: dftfunc_ADMML2 = 2

! Structure for different dft grid 
TYPE GridItem
INTEGER           :: RADIALGRID   !(1 = GC2, 2 = LMG, 3 = TURBO)
INTEGER           :: PARTITIONING !integer in following list
!(1=SSF, 2=Becke, 3=Becke-original, 4=block, 5=blockssf, 6=cartesian)
LOGICAL           :: ZdependenMaxAng !(default 0)
INTEGER           :: GRIDDONE !IF GRID HAS BEEN CREATED 
REAL(REALK)       :: RADINT
INTEGER           :: ANGMIN
INTEGER           :: ANGINT
INTEGER           :: HRDNES ! hardness of the partition function in the becke schemes
!INTEGER           :: IPRUNE
LOGICAL           :: NOPRUN                  !.TRUE.
INTEGER           :: TURBO
INTEGER           :: NBUFLEN                 !0
LOGICAL           :: NEWGRID                 !.FALSE.
integer           :: Id !1,2 or 3 corresponds to (Grid_Default,Grid_ADMML2,..)
integer           :: NBAST
END TYPE GridItem

!> Keywords and input to gridgeneration and exchange-correlation calculation
TYPE DFTparam
INTEGER           :: iGrid !which gridObject to use
TYPE(Griditem)    :: GridObject(3) !3 different grids (Grid_Default,Grid_ADMML2,..)
INTEGER           :: iDFTtype
CHARACTER(len=80) :: dftfuncObject(2)    !2 different functionals (dftfunc_Default,dftfunc_ADMML2,..)
!
INTEGER           :: RADIALGRID !(1 = GC2, 2 = LMG, 3 = TURBO)
INTEGER           :: PARTITIONING !integer in following list
!(1=SSF, 2=Becke, 3=Becke-original, 4=block, 5=blockssf, 6=cartesian)
LOGICAL           :: ZdependenMaxAng !(default 0)

INTEGER           :: GRIDDONE !IF GRID HAS BEEN CREATED 
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
LOGICAL           :: NEWGRID                 !.FALSE.
Logical           :: testNelectrons          !.TRUE.
Logical           :: LB94      !van Leeuwen-Baerends correction
Logical           :: CS00      !Casida-Salahub asymptotic correction
REAL(REALK)       :: CS00shift !Casida-Salahub shift - if 0 use Zhan-Nichols-Dixon shift
REAL(REALK)       :: CS00eHOMO !energy of HOMO to use with Zhan-Nichols-Dixon shift
REAL(REALK)       :: CS00ZND1   !Zhan-Nichols-Dixon shift parameter
REAL(REALK)       :: CS00ZND2   !Zhan-Nichols-Dixon shift parameter
REAL(REALK)       :: HFexchangeFac
!type(dft_grid)    :: L2GRID     !Grid parameters for level 2/ADMM grid
!type(dft_grid)    :: L3GRID     !Grid parameters for level 3/regular grid
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
subroutine init_gridObject(dft,gridObject)
TYPE(GridItem) :: gridObject(:)
TYPE(dftparam) :: dft
!
integer :: iGrid

do iGrid=1,size(gridObject)
   GridObject(iGrid)%RADIALGRID = dft%RADIALGRID
   GridObject(iGrid)%PARTITIONING = dft%PARTITIONING
   GridObject(iGrid)%ZdependenMaxAng = dft%ZdependenMaxAng
   GridObject(iGrid)%griddone = 0
   GridObject(iGrid)%radint = dft%radint
   GridObject(iGrid)%angmin = dft%angmin
   GridObject(iGrid)%angint = dft%angint
   GridObject(iGrid)%HRDNES = dft%HRDNES
   GridObject(iGrid)%NOPRUN = dft%NOPRUN
   GridObject(iGrid)%TURBO = dft%TURBO
   GridObject(iGrid)%nbuflen = dft%nbuflen
   GridObject(iGrid)%NewGrid = dft%NewGrid
   GridObject(iGrid)%Id = iGrid
   GridObject(iGrid)%NBAST = 0
   !module parameters
   DFT_GRIDITERATIONS(iGrid) = 0
   DFT_MaxNactBast(iGrid) = 0
enddo
end subroutine init_gridObject

subroutine init_dftfunc(dft)
TYPE(dftparam) :: dft
!
integer :: iDFT
do iDFT=1,size(dft%dftfuncObject)
   DFT%DFTfuncObject(iDFT) = dft%dftfunc
enddo
end subroutine init_dftfunc

!!$SUBROUTINE store_dft_grid(grid,filename,dft)
!!$TYPE(dft_grid) :: grid
!!$TYPE(dftparam) :: dft
!!$character(21)  :: filename
!!$grid%filename       = filename
!!$grid%GRIDITERATIONS = GRIDITERATIONS
!!$grid%dftfunc        = dft%dftfunc
!!$grid%nbuflen        = dft%nbuflen
!!$grid%maxNactBAST    = dft%maxNactBAST
!!$END SUBROUTINE store_dft_grid
!!$
!!$SUBROUTINE get_dft_grid(grid,filename,dft)
!!$TYPE(dft_grid) :: grid
!!$TYPE(dftparam) :: dft
!!$character(21)  :: filename
!!$
!!$IF (filename.EQ.grid%filename) THEN
!!$  GRIDITERATIONS  = grid%GRIDITERATIONS
!!$  dft%dftfunc     = grid%dftfunc       
!!$  dft%nbuflen     = grid%nbuflen       
!!$  dft%maxNactBAST = grid%maxNactBAST   
!!$ELSE
!!$  call lsquit('Filename error in get_dft_grid',-1)
!!$ENDIF
!!$END SUBROUTINE get_dft_grid

END MODULE dft_typetype