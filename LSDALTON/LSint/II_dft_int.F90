!> @file
!> Module containing main exchange-correlation integral driver, and routines to evaluate AOs and electron-densities
!> \brief 
!> \author T. Kjaergaard
!> \date 2009 
MODULE IIDFTINT
use gridgenerationmodule
use dft_memory_handling
use Integralparameters
!use memory_handling
!WARNING you must not add memory_handling, all memory goes through 
!grid_memory_handling  module so as to determine the memory used in this module.
use precision
use TYPEDEF
use dft_type
use dft_typetype
use IIDFTKSMWORK
use LS_UTIL,only: DGEMM_TS
private
public :: II_DFTINT, TEST_NELECTRONS, II_DFTDISP, II_DFTsetFunc, II_DFTaddFunc
! Notes 
! OLD NOTATION !   NEW NOTATION
!==========================================================
! NUCO         !   SHELLNPRIM   !nprimitives for the shell
! NHKT         !   SHELLANGMOM  !angmom for the shell
! NCENT        !   SHELL2ATOM   !Atom for the shell
! KMAX         !   MAXNSHELL    !maximum number of shells
! NHTYP        !   maxAngmom    !maximum angmom+1
! JSTRT        !   PRIEXPSTART  !index to start in PRIEXP for shell
!========================================================== 
CONTAINS
!> \brief wrapper exchange-correlation integral routine that build basinf.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_DFTINT(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND,&
     & CB, DFTDATA,UNRES,ELECTRONS,USE_MPI)
use BUILDAOBATCH
IMPLICIT NONE
INTEGER,intent(in)     :: LUPRI,IPRINT,NBAST,NDMAT,NGEODRV
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
EXTERNAL CB !NAME OF SUBROUTINE TO CALL 
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
TYPE(LSSETTING) :: SETTING
REAL(REALK)     :: ELECTRONS(NDMAT)
LOGICAL         :: DOLND !do london
LOGICAL         :: UNRES !unrestricted calc
LOGICAL         :: USE_MPI !use MPI ?
!
TYPE(BASINF)  :: BAS
LOGICAL       :: DOGGA,DOMETA
INTEGER       :: NDER,NTYPSO,IGEODRV,NSOB,igrid,GRIDDONE
INTEGER       :: maxNactbast,GRIDITERATIONS
REAL(REALK),parameter :: D0=0E0_realk

igrid = SETTING%scheme%DFT%igrid
GRIDDONE = SETTING%scheme%DFT%GridObject(igrid)%GRIDDONE
IF(GRIDDONE.EQ.0)THEN
   maxNactbast = 0
   GRIDITERATIONS = 0
ELSE
   maxNactbast = dft_maxNactbast(igrid)
   GRIDITERATIONS = dft_GRIDITERATIONS(igrid)
ENDIF
IF (SETTING%SCHEME%CONTANG) CALL LSQUIT('Error in II_DFTINT. ContAng option not implemented!',-1)

ELECTRONS=d0
CALL BUILD_BASINF(LUPRI,IPRINT,BAS,SETTING,GRIDDONE,.FALSE.)

CALL DFT_DOGGA_DOMETA(DOGGA,DOMETA)
IGEODRV = NGEODRV
IF (DOGGA) IGEODRV = IGEODRV + 1
CALL II_SETUPSOS(IGEODRV,DOLND,NBAST,NDER,NTYPSO,NSOB)
CALL II_DFTINT1(LUPRI,IPRINT,DMAT,NBAST,SETTING%SCHEME%DFT,BAS,&
     & NDMAT,NGEODRV,DOLND,CB,DFTDATA,ELECTRONS,NDER,NTYPSO,NSOB,DOGGA,&
     & SETTING%MOLECULE(1)%p%nELECTRONS,SETTING%SCHEME%noOMP,&
     & UNRES,USE_MPI,setting%numnodes,setting%node,&
     & SETTING%scheme%DFT%GridObject(igrid),maxNactbast,&
     & GRIDITERATIONS)
CALL FREE_BASINF(BAS)

dft_maxNactbast(igrid) = maxNactbast
dft_GRIDITERATIONS(igrid) = GRIDITERATIONS

END SUBROUTINE II_DFTINT

!> \brief main exchange-correlation integral routine
!> \author T. Kjaergaard
!> \date 2010
!>
!>  The routine first generates a grid which means that it
!>  generates a number of gridpoints/coordinates with an associated weight
!>  then follows a loop over all gridpoints - the actual implementation is 
!>  a loop where in each iteration a batch of gridpoints+weights is 
!>  read from file. In addition to the coordinates and weights there is
!>  also read which shells contribute to this batch of gridpoints. 
!>  The shells are transformed into orbitalindexes (BLOCKS) which state
!>  which orbitals contribute to the gridpoints (which orbitals are active). 
!>  the number of these active orbitals are determined (NactBas) 
!>  INXACT is then created which for a given active index give the orbitalindex
!>  An active Dmat is then created from the full Dmat.
!>  Next the (active) gaussian atomic orbitals are created (GAO) for each gridpoint 
!>  The GAO have dimensions (ngridpoints,NactBAS,NTYPSO)
!>  NTYPSO is 1 for LDA and 4 for GGA (because we both need the atomic orbitals and 
!>  the first geometrical derivatives). In general  
!>
!>  GAO: evaluated orbitals for a batch of grid points.
!>  GAO(:,:,1) contains orbital values.
!>  GAO(:,:,2:4) contains first geom. derivatives. - if requested.
!>  GAO(:,:,5:10) contains second derivatives - if requested.
!>  GAO(:,:,11:20) contains third derivatives - if requested.
!>  After requested geometric derivatives, london related derivatives are placed.
!>  
!>  From the GAO the density (RHO) is calculated for each grid point and the 
!>  gradient of the density (GRAD) is calculated. 
!>  The number of electrons is calculated (an estimate of the error of the grid)
!>  and a generic function call CB is called. CB is an input argument, so when 
!>  you call II_DFTINT1 you can call it with for instance II_DFT_KSMGGA, and 
!>  II_DFTINT1 therefore calls II_DFT_KSMGGA. II_DFT_KSMGGA and all other 
!>  'worker' routines are for now in the II_dft_ksm.f90 file. 
SUBROUTINE II_DFTINT1(LUPRI,IPRINT,DMAT,NBAST,DFT,BAS,&
     & NDMAT,NGEODRV,DOLND,CB,DFTDATA,ELECTRONS,&
     & NDER,NTYPSO,NSOB,DOGGA,NELECTRONS,noOMP,UNRES,&
     & USE_MPI,numnodes,node,GridObject,maxNactbast,&
     & GRIDITERATIONS) 
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)  :: IPRINT
!> the density matrix
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> number of basis functions
INTEGER,intent(in)  :: NBAST
!> DFT parameters
type(DFTparam),intent(INOUT) :: DFT
!> Basis-set information
TYPE(BASINF),intent(INOUT)  :: BAS
!> number of density matrices
INTEGER,intent(in)  :: NDMAT
!> the order of geometrical derivative
INTEGER,intent(in)  :: NGEODRV
!> do london derivative GAOs?
LOGICAL,intent(in)  :: DOLND 
!> Unrestricted 
LOGICAL,intent(in)  :: UNRES 
!> the NAME OF external SUBROUTINE TO CALL 
EXTERNAL CB 
!> contains the data that must be given to the CB routine
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> number of electrons from nummerical integration
REAL(REALK),intent(inout) :: ELECTRONS(NDMAT)
!> level of geometrical derivatives
INTEGER,intent(in)  :: NDER
!> number of GAOs and geomderiv derivativ GAOs
INTEGER,intent(in)  :: NTYPSO
!> place in GAO to place LONDON contrib
INTEGER,intent(in)  :: NSOB
!> is it a GGA calc or a LDA calc
LOGICAl,intent(in)  :: DOGGA
!> the number of electrons we should have
INTEGER,intent(in)  :: NELECTRONS
!> shoould OpenMP be deactivated
LOGICAL,intent(in)  :: noOMP
!> should we use MPI
LOGICAL,intent(in) :: USE_MPI
!> MPI info
INTEGER(kind=ls_mpik),intent(in) :: numnodes,node
type(Griditem) :: GridObject
integer :: maxNactBAST,GRIDITERATIONS
!
integer :: NPOINTS,NLEN,NSHELLBLOCKS,NROW,NCOL
INTEGER :: iprune,L_prev,L_curr,IPT,spSIZE,L,NCURLEN,I,J
INTEGER,pointer     :: SPINDEX(:)
INTEGER   :: KCKTA,KHKTA,SHELLANGMOMA,IT,IBUF,IBUF_PREV,IDUM,ILEN,K,NactBAS,NRED,NRED2,IDMAT
INTEGER   :: spsize2,IJ
REAL(REALK) :: TELECTRONS(ndmat),ERROR,TS,TE,DMAX,GAOMAX
REAL(REALK),pointer :: SPHMAT(:)
LOGICAL     :: CHECKELS,SETIT,LDUM,PRINTTIM
integer,pointer :: LVALUE(:,:),MVALUE(:,:),NVALUE(:,:)
integer :: GC2,mmm

TELECTRONS(1:ndmat) = 0E0_realk
IT=0
IF(GridObject%GRIDDONE .EQ. 0)THEN
   GRIDITERATIONS=0
   GridObject%NBUFLEN=0
   maxNactBAST = 0
   SETIT=.TRUE.
   PRINTTIM=.TRUE.
ELSE
   SETIT=.FALSE.
   PRINTTIM=.FALSE.
ENDIF
!pruning: per default on
IPRUNE = 1
IF (GridObject%NOPRUN) IPRUNE = 0

CALL LSTIMER('START',TS,TE,LUPRI)
GridObject%NBUFLEN= 1024
BoxMemRequirement = 160000
CALL GenerateGrid(NBAST,GridObject%radint,GridObject%angmin,GridObject%angint,&
     & GridObject%HRDNES,iprune,BAS%natoms,BAS%X,BAS%Y,BAS%Z,BAS%Charge,GridObject%GRIDDONE,&
     & BAS%SHELL2ATOM,BAS%SHELLANGMOM,BAS%SHELLNPRIM,BAS%MAXANGMOM,&
     & BAS%MAXNSHELL,BAS%MXPRIM,BAS%PRIEXP,BAS%PRIEXPSTART,BAS%RSHEL,IT,GridObject%TURBO,&
     & GridObject%nbuflen,GridObject%RADIALGRID,GridObject%ZdependenMaxAng,&
     & GridObject%PARTITIONING,BAS%nstart,MaxNactBast,LUPRI,&
     & IPRINT,USE_MPI,numnodes,node,GridObject%Id,GridObject%numnodes)
IF(GridObject%numnodes.NE.numnodes)THEN
   print*,'GridObject%numnodes',GridObject%numnodes
   print*,'numnodes',numnodes
   call lsquit('dim mismatch in Number of nodes used to construct/calc grid',-1)
ENDIF
IF(PRINTTIM)THEN
#ifdef VAR_MPI
   IF (infpar%mynum.EQ.infpar%master) THEN
#endif
      CALL LSTIMER('gridgeneration2',TS,TE,LUPRI)
#ifdef VAR_MPI
   ENDIF
#endif
ENDIF
IF(SETIT)THEN 
   !set global variable concerning how many time we should read from disk
   GRIDITERATIONS=IT
ELSE
   !set local variable from global variable
   IT=GRIDITERATIONS
ENDIF

! Setting spherical transformation quantities
IF(BAS%MAXANGMOM .GE. 3)THEN
   spSIZE=0
   DO L=2,BAS%MAXANGMOM-1
      spSIZE=spSIZE+(2*L+1)*(L+1)*(L+2)/2
   ENDDO
   call mem_dft_alloc(SPHMAT,spSIZE)
   CALL LS_DZERO(SPHMAT,spSIZE)
   call mem_dft_alloc(SPINDEX,BAS%MAXANGMOM)
   spsize2 = BAS%maxangmom
   CALL Build_PRECALCULATED_SPHMAT(LUPRI,BAS%MAXANGMOM-1,spSIZE,SPHMAT,SPINDEX)
ELSE
   spSIZE = 9
   call mem_dft_alloc(SPHMAT,spSIZE)
   call mem_dft_alloc(SPINDEX,2)
   SPINDEX(1)=1
   SPINDEX(2)=1
   spsize2 = 2
ENDIF
MMM = BAS%MAXANGMOM*(BAS%MAXANGMOM+1)/2*BAS%MAXANGMOM*(BAS%MAXANGMOM+1)/2
call mem_dft_alloc(LVALUE,MMM,BAS%MAXANGMOM)
call mem_dft_alloc(MVALUE,MMM,BAS%MAXANGMOM)
call mem_dft_alloc(NVALUE,MMM,BAS%MAXANGMOM)
LVALUE=0
MVALUE=0
NVALUE=0
DO SHELLANGMOMA=1,BAS%MAXANGMOM
   KCKTA = SHELLANGMOMA*(SHELLANGMOMA+1)/2
   IJ=0
   DO I = 1,KCKTA
      DO J = 1,I
         IJ=IJ+1
         LVALUE(IJ,SHELLANGMOMA)=SHELLANGMOMA-I
         MVALUE(IJ,SHELLANGMOMA)=I-J
         NVALUE(IJ,SHELLANGMOMA)=J-1
      ENDDO
   ENDDO
ENDDO

BoxMemRequirement = MXBLLEN*MaxNactBast
CALL II_XC_GRID_LOOP(NDMAT,GridObject%NBUFLEN,IT,nbast,maxNactbast,BAS%maxnshell,&
     & BAS%CC,BAS%ushells,BAS%CCSTART,BAS%CCINDEX,TELECTRONS,BAS%CENT,DFTdata,&
     & DFT%DFThri,DMAT,dogga,dolnd,ntypso,BAS%PRIEXPSTART,lupri,&
     & BAS%mxprim,nder,BAS%shellangmom,MMM,BAS%maxangmom,nsob,BAS%nstart,BAS%priexp,&
     & DFT%rhothr,spsize,&
     & sphmat,spsize2,spindex,LVALUE,MVALUE,NVALUE,noOMP,CB,unres,node,&
     & GridObject%Id,BAS%Spherical)

call mem_dft_dealloc(SPHMAT)
call mem_dft_dealloc(SPINDEX)
call mem_dft_dealloc(LVALUE)
call mem_dft_dealloc(MVALUE)
call mem_dft_dealloc(NVALUE)
!
! Test on the number of electrons
!   
ELECTRONS(1:ndmat) = TELECTRONS(1:ndmat)
#ifndef VAR_MPI
IF (DFT%testNelectrons) THEN
  CALL TEST_NELECTRONS(ELECTRONS,NELECTRONS,NDMAT,DFT%DFTELS,UNRES,LUPRI)
ENDIF
#endif

END SUBROUTINE II_DFTINT1

SUBROUTINE TEST_NELECTRONS(ELECTRONS,NELECTRONS,NDMAT,DFTELS,UNRES,LUPRI)
  implicit none
  integer,intent(in) :: NELECTRONS,LUPRI,ndmat
  real(realk),intent(in) :: ELECTRONS(NDMAT),DFTELS
  logical :: UNRES
  !
  integer :: ndmat2,idmat
  real(realk) :: NELE,ERROR
  ndmat2 = ndmat
  IF(UNRES)ndmat2 = ndmat/2
  NELE = REAL(NELECTRONS)
  DO IDMAT=1,NDMAT2
     ERROR  = ELECTRONS(IDMAT) - NELE
     IF (ABS(ERROR) .GT. DFTELS*NELE) THEN
        WRITE (LUPRI,'(4(/2X,A,F14.6),/2X,A)')                              &
             &' Number of electrons from numerical integration:',ELECTRONS(IDMAT), &
             &' Number of electrons from orbital occupations:  ',NELE,&
             &' Error in the number of electrons:              ',ERROR,     &
             &' Error larger than DFTELS (set input):          ',DFTELS,    &
             &' MPI Calculation aborted.'
        CALL LSQUIT('Wrong number of electrons in XC. Calculation aborted.',lupri)
     END IF
  ENDDO
END SUBROUTINE TEST_NELECTRONS

SUBROUTINE II_XC_GRID_LOOP(NDMAT,NBUFLEN,IT,nbast,maxNactbast,maxnshell,&
     & CC,ushells,CCSTART,CCINDEX,TELECTRONS,CENT,DFTdata,&
     & DFThri,DMAT,dogga,dolnd,ntypso,PRIEXPSTART,lupri,&
     & mxprim,nder,shellangmom,MMM,maxangmom,nsob,nstart,priexp,rhothr,spsize,&
     &sphmat,spsize2,spindex,LVALUE,MVALUE,NVALUE,noOMP,CB,unres,node,&
     &GridId,Spherical)
use IIDFTKSMWORK
implicit none
integer,intent(in) :: NDMAT,LUPRI,MAXNSHELL,nbast,NBUFLEN,gridid
integer :: IT
!> Largest number of active orbitals
integer,intent(inout) ::  maxNactbast
!> number of GAOs and geomderiv derivativ GAOs
INTEGER,intent(in)  :: NTYPSO
!> MXPRIM is the total number of (unique) primitive orbitals
INTEGER,intent(in)  :: MXPRIM
!> level of geometrical derivatives
INTEGER,intent(in)  :: NDER
!> the angular momentum for each shell
INTEGER,intent(in)    :: SHELLANGMOM(MAXNSHELL)
!> the maximum angular momentum
INTEGER,intent(in)  :: MAXANGMOM
!> MMM=MAXANGMOM*(MAXANGMOM+1)/2*MAXANGMOM*(MAXANGMOM+1)/2 
INTEGER,intent(in)  :: MMM
!> place in GAO to place LONDON contrib
INTEGER,intent(in)  :: NSOB
real(realk),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
!> the number of unique shells (also different contraction coefficients matrices) 
INTEGER,intent(in)  :: ushells
!> a contraction coefficient matrix for all the unique shells
TYPE(LSMATRIX),intent(in):: CC(ushells)
!> hmmm I do not remember
INTEGER,intent(in) :: CCSTART(ushells)
!> contraction index for given shell, for shellindex => unique shellindex
INTEGER,intent(in) :: CCINDEX(MAXNSHELL)
real(realk),intent(inout) :: TELECTRONS(ndmat)
!> X,Y,Z coordinate for each shell
REAL(REALK),intent(in):: CENT(3,MAXNSHELL)
TYPE(DFTDATATYPE),intent(inout) :: DFTDATA
!> threshold for value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> is it a GGA calc or a LDA calc
LOGICAl,intent(in) :: DOGGA
LOGICAl :: DOMETA = .FALSE. !Should be changed to input (perhaps integer 0 lda, 1 gga, 2 meta or something)
!> shoould OpenMP be deactivated
LOGICAL  :: noOMP
!> do london derivative GAOs?
LOGICAL,intent(in) :: DOLND 
!> the index to start in PRIEXP(MXPRIM) for a given shell index 
INTEGER,intent(in) :: PRIEXPSTART(MAXNSHELL)
!> for a given shell index it gives the corresponding starting orbital
INTEGER,intent(in)    :: NSTART(MAXNSHELL)
!> the unique primitve exponents
REAL(REALK),intent(in) :: PRIEXP(MXPRIM)
!> threshold for value of RHO
REAL(REALK),intent(in) :: RHOTHR
!> size of sphmat
integer,intent(in) :: spsize
!> size of spindex
integer,intent(in) :: spsize2
!> size of spherical transformation matrices
real(realk),intent(in) :: sphmat(spsize)
!> unrestricted calc
logical,intent(in) :: unres,Spherical
integer,intent(in) :: spindex(spsize2)
integer,intent(in) :: LVALUE(MMM,MAXANGMOM),MVALUE(MMM,MAXANGMOM),NVALUE(MMM,MAXANGMOM)
integer(kind=ls_mpik),intent(in) ::  node
!> the NAME OF external SUBROUTINE TO CALL 
EXTERNAL CB 
!INTEGER,intent(in) :: XC_JOB 
!
!real(realk) :: RHOA(MXBLLEN,NDMAT),GRADA(3,MXBLLEN,NDMAT),WEIGHT(NBUFLEN)
!real(realk),target :: COOR(3,NBUFLEN)
!real(realk) :: ACTIVE_DMAT(NBAST*NBAST*NDMAT)
!real(realk) :: GAO(MXBLLEN*NBAST*NTYPSO) 
INTEGER :: NactBAS
INTEGER,pointer :: SHELLBLOCKS(:,:),BLOCKS(:,:)
real(realk),pointer :: RHOA(:,:),GRADA(:,:,:),TAU(:,:),WEIGHT(:)
real(realk),pointer :: COOR(:,:)
real(realk),pointer :: ACTIVE_DMAT(:)
real(realk),pointer :: GAO(:) 
real(realk),pointer :: GAOGMX(:)
!INTEGER,pointer :: SHELLBLOCKS(:,:),BLOCKS(:,:)
!INTEGER :: NactBAS
REAL(REALK) :: GAOMAX
REAL(REALK),pointer :: COOR_pt(:,:),WORK(:)
integer :: XX,NSHELLBLOCKS,NLEN,IPT,NCURLEN,I,IDMAT,WORKLENGTH,WORKLENGTH2
integer :: W1,W2,W3,W4,W5,W6,W7,W8,nthreads,tid,lugrid,idmat1,idmat2
integer :: myBoxMemRequirement
integer,pointer :: INXACT(:)
TYPE(DFTDATATYPE) :: myDFTDATA
character(len=22) :: filename
logical :: grid_exists
#ifdef VAR_OMP
integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
LUGRID=-1
call setNoOMP(noOMP)
call get_quadfilename(filename,nbast,node,GridId)
INQUIRE(file=filename,EXIST=grid_exists)
IF(grid_exists)THEN 
   CALL LSOPEN(LUGRID,filename,'OLD','UNFORMATTED')
   rewind(LUGRID)
ELSE
   print*,'file with filename: ',filename,' does not exist'
   print*,'but it should exist'
   call lsquit('missing file in XC integration',-1)
ENDIF
myBoxMemRequirement = BoxMemRequirement
IF(.NOT.noOMP) call mem_dft_TurnONThread_Memory()
!$OMP PARALLEL IF(.NOT.noOMP) PRIVATE(XX,NSHELLBLOCKS,SHELLBLOCKS,COOR,COOR_pt,WEIGHT,IPT,NCURLEN,&
!$OMP ACTIVE_DMAT,BLOCKS,RHOA,GRADA,TAU,GAO,NLEN,NactBAS,INXACT,I,IDMAT,myDFTdata,WORK,WORKLENGTH,W1,&
!$OMP W2,W3,W4,W5,W6,W7,W8,tid,nthreads,IDMAT1,IDMAT2,WORKLENGTH2,GAOGMX,&
!$OMP GAOMAX) FIRSTPRIVATE(myBoxMemRequirement)& 
!$OMP SHARED(ushells,CC,CCSTART,CCINDEX,CENT,Dmat,Spherical,dometa,&
!$OMP dogga,PRIEXPSTART,nsob,shellangmom,MMM,maxangmom,nstart,priexp,rhothr,spsize,spsize2,sphmat,&
!$OMP spindex,DFTdata,noOMP,NBUFLEN,it,lupri,nbast,maxnshell,ndmat,nder,mxprim,unres,&
!$OMP ntypso,dfthri,dolnd,LVALUE,MVALUE,NVALUE,MaxNactBast,lugrid) REDUCTION(+:TELECTRONS)
IF(.NOT.noOMP) call init_dft_threadmemvar()
#ifdef VAR_OMP
nthreads=OMP_GET_NUM_THREADS()
tid = omp_get_thread_num()
#else
nthreads=1
tid = 0
#endif
!$OMP CRITICAL (initdftDATAblock)
call mem_dft_alloc(SHELLBLOCKS,2,MAXNSHELL)
call mem_dft_alloc(BLOCKS,2,MAXNSHELL)
call mem_dft_alloc(INXACT,maxNactBAST)
call mem_dft_alloc(RHOA,MXBLLEN,NDMAT)
call mem_dft_alloc(GRADA,3,MXBLLEN,NDMAT)
call mem_dft_alloc(TAU,MXBLLEN,NDMAT)
call mem_dft_alloc(WEIGHT,NBUFLEN)
call mem_dft_alloc(COOR,3,NBUFLEN)
call mem_dft_alloc(ACTIVE_DMAT,maxNactBAST*maxNactBAST*NDMAT)
call mem_dft_alloc(GAO,MXBLLEN*maxNactBAST*NTYPSO) 
call mem_dft_alloc(GAOGMX,maxNactBAST)
!mem required in the ksm_worker routines
WORKLENGTH = DFTdata%nWorkNactBastNblen*myBoxMemRequirement &
     &     + DFTdata%nWorkNactBast*maxNactBAST &
     &     + DFTdata%nWorkNactBastNactBast*maxNactBAST*maxNactBAST
!mem required in integral loop
IF(DOGGA) THEN
   WORKLENGTH2 = maxNactBAST*maxNactBAST+maxNactBAST+5*MIN(myBoxMemRequirement,MXBLLEN*maxNactBAST)
ELSE
   WORKLENGTH2 = maxNactBAST*maxNactBAST+maxNactBAST+2*MIN(myBoxMemRequirement,MXBLLEN*maxNactBAST)
ENDIF
IF(WORKLENGTH.EQ.0)THEN
   WORKLENGTH = WORKLENGTH2 
ELSE
   WORKLENGTH = MAX(WORKLENGTH,WORKLENGTH2)
ENDIF
call mem_dft_alloc(WORK,WORKLENGTH)
call initDFTdatatype(myDFTdata)
call copyDFTdata(myDFTdata,DFTdata)
!$OMP END CRITICAL (initdftDATAblock)

DO XX=1+tid,IT,nthreads
!$OMP CRITICAL (gridREAD)
   CALL READ_GRIDPOINTS(LUGRID,NSHELLBLOCKS,SHELLBLOCKS,MAXNSHELL,NBUFLEN,&
        & COOR,WEIGHT,NLEN)
!$OMP END CRITICAL (gridREAD)

!  READ_GRIDPOINTS reads grid data from the grid file.
!
!  The data is read to coor array that will contain
!  the grid point coordinates and to weight that will contain
!  the associated weights. These arrays must be preallocated and have
!  length nbuflen. 
!
!  The information about the basis function shells relevant for 
!  this chunk of grid points is returned
!  via nshellblocks and nshellblocks arguments. We do not
!  return the active shell numbers each other separately. They can be
!  packed instead into blocks of consecutive active shells. Argument
!  nshellblocks will be set to the number of these blocks and shellblocks
!  entries will be filled with the beginnings (shellblocks[1][:]) and
!  ends (shellblocks[2][:]) of the blocks of active shells.
#ifdef VAR_LSDEBUGINT
   IF(NLEN .LE. 0) THEN
      CALL LSQUIT('SOMETHING WRONG IN THE DFT LOOP',lupri)
   ENDIF
#endif
   CALL SHELL_TO_ORB(LUPRI,NSHELLBLOCKS,SHELLBLOCKS,BLOCKS,NSTART,MAXNSHELL,NBAST) !OUTPUT: BLOCKS
   CALL DETERMINE_NACTIVEORB(NactBAS,NSHELLBLOCKS,BLOCKS,MAXNSHELL) !OUTPUT: NactBAS
   CALL BUILD_INXACT(NactBAS,NSHELLBLOCKS,BLOCKS,INXACT(1:NactBas),MAXNSHELL) !OUTPUT: INXACT
   CALL CONSTRUCT_ACTIVE_DMAT(NSHELLBLOCKS,MAXNSHELL,BLOCKS,DMAT,NBAST,ACTIVE_DMAT,NactBAS,NDMAT,INXACT(1:NactBas))
   CALL SHELL_TO_ACTORB(NSHELLBLOCKS,BLOCKS,MAXNSHELL) !CHANGE ORBBLOCKS
   DO IPT = 1, NLEN, MXBLLEN
      NCURLEN=MIN(MXBLLEN,NLEN-IPT+1)
      COOR_pt => COOR(:,IPT:(IPT+NCURLEN-1))
      CALL II_BLGETSOS(LUPRI,NCURLEN,GAO,COOR_pt,&
      &                NSHELLBLOCKS,SHELLBLOCKS,NactBAS,DOLND,DOGGA,DFTHRI,0,&
      &                NTYPSO,NSOB,MMM,MAXANGMOM,spSIZE,SPHMAT,SPINDEX,MAXNSHELL,NSTART,SHELLANGMOM&
      &                ,CENT,PRIEXPSTART,CC,ushells,CCSTART,CCINDEX,MXPRIM,&
      &                PRIEXP,NDER,LVALUE,MVALUE,NVALUE,WORK(1:NactBAS),Spherical)
      ! RETURNS: GAO: evaluated orbitals for a batch of grid points.
      !     GAO(:,:,1) contains orbital values.
      !     GAO(:,:,2:4) contains first geom. derivatives. - if requested.
      !     GAO(:,:,5:10) contains second derivatives - if requested.
      !     GAO(:,:,11:20) contains third derivatives - if requested.
      !     After requested geometric derivatives, london related derivatives
      !     are placed.
      IF(DOMETA) THEN
         W1 = 1; W2 = NactBAS*NactBAS                   !Dred(NactBas*NactBas)
         W3 = NactBAS*NactBAS+1
         W4 = NactBAS*NactBAS+4*NCURLEN*NactBAS         !GAORED(NCURLEN,NactBAS,4)
         W5 = NactBAS*NactBAS+4*NCURLEN*NactBAS+1       
         W6 = NactBAS*NactBAS+8*NCURLEN*NactBAS         !TMP(NCURLEN,NactBAS,4)
         CALL II_GETRHO_BLOCKED_META(LUPRI,ACTIVE_DMAT,NactBAS,GAO,NTYPSO,&
              & NSHELLBLOCKS,MAXNSHELL,BLOCKS,NCURLEN,NDMAT,RHOA,GRADA,TAU,RHOTHR,MXBLLEN,WORK(W1:W2),&
              & GAOGMX,GAOMAX,WORK(W3:W4),WORK(W5:W6),maxNactBAST)
         ! computes rho_a and gradients of rho, Dmat is a density matrix
         ! (it can be a total density and then one will get total density,
         ! or it can be an alpha/beta density    assert(NTYPSO>=NRHO)
      ELSE IF(DOGGA) THEN
         W1 = 1; W2 = NactBAS*NactBAS                   !Dred(NactBas*NactBas)
         W3=  NactBAS*NactBAS+1
         W4 = NactBAS*NactBAS+4*NCURLEN*NactBAS         !GAORED(NCURLEN,NactBAS,4)
         W5 = NactBAS*NactBAS+4*NCURLEN*NactBAS+1       
         W6 = NactBAS*NactBAS+5*NCURLEN*NactBAS         !TMP(NCURLEN,NactBAS)
         CALL II_GETRHO_BLOCKED_GGA(LUPRI,ACTIVE_DMAT,NactBAS,GAO,NTYPSO,&
              & NSHELLBLOCKS,MAXNSHELL,BLOCKS,NCURLEN,NDMAT,RHOA,GRADA,&
              & RHOTHR,MXBLLEN,WORK(W1:W2),GAOGMX,GAOMAX,WORK(W3:W4),&
              & WORK(W5:W6),maxNactBAST)
         ! computes rho_a and gradients of rho, Dmat is a density matrix
         ! (it can be a total density and then one will get total density,
         ! or it can be an alpha/beta density    assert(NTYPSO>=NRHO)
      ELSE
         W1 = 1; W2 = NactBAS*NactBAS                   !DRED(NactBas*NactBas)
         W3=  NactBAS*NactBAS+1
         W4 = NactBAS*NactBAS+NCURLEN*NactBAS           !GAORED(NCURLEN,NactBAS)
         W5 = NactBAS*NactBAS+NCURLEN*NactBAS+1       
         W6 = NactBAS*NactBAS+2*NCURLEN*NactBAS         !TMP(NCURLEN,NactBAS)
         CALL II_GETRHO_BLOCKED_LDA(LUPRI,ACTIVE_DMAT,NactBAS,GAO,NTYPSO,&
              & NSHELLBLOCKS,MAXNSHELL,BLOCKS,NCURLEN,NDMAT,RHOA,RHOTHR,&
              & MXBLLEN,WORK(W1:W2),GAOGMX,GAOMAX,WORK(W3:W4),WORK(W5:W6),&
              & maxNactBAST)
      END IF
      IF(UNRES)THEN
       DO IDMAT = 1,NDMAT/2
        IDMAT1 = 1+(IDMAT-1)*2
        DO I = 1, NCURLEN
         TELECTRONS(IDMAT1)=TELECTRONS(IDMAT1)+WEIGHT(IPT+I-1)*RHOA(I,IDMAT1)
        END DO
        IDMAT2 = 2+(IDMAT-1)*2
        DO I = 1, NCURLEN
         TELECTRONS(IDMAT1)=TELECTRONS(IDMAT1)+WEIGHT(IPT+I-1)*RHOA(I,IDMAT2)
        END DO
       ENDDO
      ELSE
       DO IDMAT = 1,NDMAT
        DO I = 1, NCURLEN
         TELECTRONS(IDMAT)=TELECTRONS(IDMAT)+WEIGHT(IPT+I-1)*RHOA(I,IDMAT) 
        END DO
       ENDDO
      ENDIF
      CALL CB(LUPRI,NCURLEN,NSHELLBLOCKS,BLOCKS(:,:),INXACT(:),NactBas,NBAST,NDMAT,ACTIVE_DMAT(:),NTYPSO,GAO(:),&
           &RHOA(:,:),GRADA(:,:,:),TAU(:,:),MXBLLEN,COOR(:,IPT:IPT+NCURLEN-1),WEIGHT(IPT:IPT+NCURLEN-1),&
           &myDFTDATA,DFTDATA,RHOTHR,DFTHRI,WORK,WORKLENGTH,GAOGMX,GAOMAX,maxNactBAST)
   ENDDO
ENDDO

!$OMP CRITICAL (freeDFTdatablock)
call mem_dft_dealloc(RHOA)
call mem_dft_dealloc(GRADA)
call mem_dft_dealloc(TAU)
call mem_dft_dealloc(WEIGHT)
call mem_dft_dealloc(COOR)
call mem_dft_dealloc(ACTIVE_DMAT)
call mem_dft_dealloc(GAO) 
call mem_dft_dealloc(INXACT)
call mem_dft_dealloc(SHELLBLOCKS)
call mem_dft_dealloc(BLOCKS)
call mem_dft_dealloc(GAOGMX)

call DFTdataReduction(myDFTdata,DFTdata)
call mem_dft_dealloc(WORK)
!$OMP END CRITICAL (freeDFTdatablock)
IF(.NOT.noOMP) call collect_thread_dft_memory()
!$OMP END PARALLEL
IF(.NOT.noOMP) call mem_dft_TurnOffThread_Memory()

CALL LSCLOSE(LUGRID,'KEEP')

END SUBROUTINE II_XC_GRID_LOOP

!> \brief determine the number of active orbitals
!> \author T. Kjaergaard
!> \date 2010
!>
!>  orbitals that contribute to the current gridpoints
SUBROUTINE DETERMINE_NACTIVEORB(NactBAS,NSHELLBLOCKS,BLOCKS,MAXNSHELL)
IMPLICIT NONE
INTEGER,intent(in) :: NSHELLBLOCKS,MAXNSHELL
INTEGER,intent(in) :: BLOCKS(2,MAXNSHELL)
INTEGER,intent(out):: NactBAS
!
INTEGER            :: IBL,IORB

NactBAS = 0
DO IBL = 1, NSHELLBLOCKS
   DO IORB = BLOCKS(1,IBL),BLOCKS(2,IBL) !ORBITALINDEX
      NactBAS = NactBAS + 1
   ENDDO
ENDDO

END SUBROUTINE DETERMINE_NACTIVEORB

!> \brief build indexhandling from active index to orbitalindex
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE BUILD_INXACT(NactBAS,NSHELLBLOCKS,BLOCKS,INXACT,MAXNSHELL) 
IMPLICIT NONE
INTEGER,intent(in) :: NSHELLBLOCKS,MAXNSHELL,NactBAS
INTEGER,intent(in) :: BLOCKS(2,MAXNSHELL)
INTEGER,intent(inout):: INXACT(NactBAS)
!
INTEGER            :: IBL,IORB,IBASIS

IBASIS = 0
DO IBL = 1, NSHELLBLOCKS
   DO IORB = BLOCKS(1,IBL),BLOCKS(2,IBL) !ORBITALINDEX
      IBASIS = IBASIS + 1
      INXACT(IBASIS) = IORB !FOR A GIVEN ACTIVEINDEX - THE CORRESPONDING ORBITALINDEX
   ENDDO
ENDDO

END SUBROUTINE BUILD_INXACT

!> \brief construct active Dmat
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE CONSTRUCT_ACTIVE_DMAT(NSHELLBLOCKS,MAXNSHELL,BLOCKS,DMAT,NBAST,ACTIVE_DMAT,NactBAS,NDMAT,INXACT)
IMPLICIT NONE 
INTEGER,intent(in)     :: NSHELLBLOCKS,MAXNSHELL,NBAST,NactBAS,NDMAT
INTEGER,intent(in)     :: BLOCKS(2,MAXNSHELL)
INTEGER,intent(in)    :: INXACT(NactBAS)
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
REAL(REALK),intent(inout):: ACTIVE_DMAT(NactBAS,NactBAS,NDMAT)
!
INTEGER               :: IMAT,IACT,JACT,I,J
REAL(REALK)           :: DTMP
DO IMAT = 1,NDMAT
   DO JACT = 1,NactBAS
      J = INXACT(JACT)
      DO IACT = 1, NactBAS
         I = INXACT(IACT)
         DTMP = DMAT(I,J,IMAT)
         ACTIVE_DMAT(IACT,JACT,IMAT) = DTMP
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CONSTRUCT_ACTIVE_DMAT

!> \brief determine ntypso and NSOB
!> \author T. Kjaergaard
!> \date 2010
!>
!>  NTYPSO determine how many GAOs, geometrical derivatives and 
!>  london derivatives there is needed. NSOB is the place in GAO
!>  to place mag derivative GAOs 
!>
SUBROUTINE II_SETUPSOS(GEODRV,DOLND,NBAST,NDER,NTYPSO,NSOB)
IMPLICIT NONE
  INTEGER :: NDER,NTYPSO,NBAST
  INTEGER :: NSOB
! index to keep track of where london derivatives fit in GAO
  INTEGER :: GEODRV
!     GEODRV - 0 for just orbital values, 1 for first order orbital
!     derivatives, 2 for laplacian.
  LOGICAL :: DOLND
  !Compute derivatives wrt magnetic field ("london" derivatives).
  !     the computed orbitals are ordered as follows:
  !     (orbital values=O, [dO/dx,dO/dy,dO/dz, [d0/dxx, ..]],
  !      [dO/dBx, dO/dBy, dO/dBz])
  NDER = GEODRV
  IF (NDER.EQ. 0) NTYPSO =  1
  IF (NDER.EQ. 1) NTYPSO =  4
  IF (NDER.EQ. 2) NTYPSO = 10
  IF (NDER.EQ. 3) NTYPSO = 20
  IF (DOLND) THEN
     NTYPSO = NTYPSO + 3
     NSOB   = NTYPSO - 2 
     IF (NDER.GT. 0) THEN
        NTYPSO = NTYPSO + 9
     END IF
  ELSE
     NSOB  = 0
  END IF

END SUBROUTINE II_SETUPSOS

!> \brief evaluates the gaussian atomic orbitals
!> \author T. Kjaergaard
!> \date 2010
!>
!>  RETURNS: GSO: evaluated orbitals for a batch of grid points.
!>     GSO(:,:,1) contains orbital values.
!>     GSO(:,:,2:4) contains first geom. derivatives. - if requested.
!>     GSO(:,:,5:10) contains second derivatives - if requested.
!>     GSO(:,:,11:20) contains third derivatives - if requested.
!>     After requested geometric derivatives, London related derivatives
!>     are placed.
!>  REWRITTEN BY T.KJAERGAARD ORIGINALLY BY  T. Helgaker sep 99, P. Salek 03
!>
SUBROUTINE II_BLGETSOS(LUPRI,NVCLEN,GAO,COOR,NBLCNT,IBLCKS,&
     &        NactBAS,DOLND,DOGGA,DFTHRI,IPRINT,NTYPSO,NSOB,MMM,MAXANGMOM,SIZE,&
     &        SPHMAT,SPINDEX,MAXNSHELL,NSTART,SHELLANGMOM,CENT,&
     &        PRIEXPSTART,CC,ushells,CCSTART,CCINDEX,MXPRIM,PRIEXP,NDER,&
     &        LVALUE,MVALUE,NVALUE,GAOMAX,Spherical)
IMPLICIT NONE
INTEGER,intent(in) :: LUPRI,NVCLEN,NBLCNT,NactBAS,IPRINT,NTYPSO,NSOB,SIZE,MAXANGMOM,MAXNSHELL
INTEGER,intent(in) :: IBLCKS(2,MAXNSHELL),SPINDEX(MAXANGMOM),NSTART(MAXNSHELL),SHELLANGMOM(MAXNSHELL),ushells
INTEGER,intent(in) :: PRIEXPSTART(MAXNSHELL),MXPRIM,NDER,CCSTART(ushells),MMM
integer,intent(in) :: LVALUE(MMM,MAXANGMOM),MVALUE(MMM,MAXANGMOM)
integer,intent(in) :: NVALUE(MMM,MAXANGMOM)
INTEGER,intent(in)     :: CCINDEX(MAXNSHELL)
REAL(REALK),intent(in) :: PRIEXP(MXPRIM)
LOGICAL,intent(in)     :: DOLND, DOGGA,Spherical
real(realk),intent(in) :: CENT(3,MAXNSHELL), COOR(3,NVCLEN)
REAL(REALK),intent(inout) :: GAO(NVCLEN,NactBAS,NTYPSO)
REAL(REALK),intent(inout) :: GAOMAX(NactBAS)
REAL(REALK),intent(in) :: SPHMAT(SIZE)
REAL(REALK),PARAMETER :: D0 = 0.0E0_realk, DTHRS = 20.0E0_realk
TYPE(LSMATRIX),intent(in) :: CC(ushells)
!
REAL(REALK),pointer :: PA(:,:), PA2(:)
REAL(REALK) :: DFTHRI,CENX,CENY,CENZ
INTEGER     :: I,IADR,JSTA,ISHELA,KHKTA,KCKTA,SHELLANGMOMA,IBL,J,oooo,IJ,SPVAL,K,IBL2
call mem_dft_alloc(PA,3,NVCLEN)
call mem_dft_alloc(PA2,NVCLEN)
SPVAL=1
!!#if 1
IF(NDER.GT. 1) THEN
   ! fixme - this won't scale linearly 
   CALL LS_DZERO(GAO(1,1,5),6*NactBAS*NVCLEN)
   IF(NDER.GT. 2)THEN
      CALL LS_DZERO(GAO(1,1,11),10*NactBAS*NVCLEN)
   END IF
END IF
!!!#else
!! this should be activated one I have modified the II_BLGETGA2 routine
!!IF(NDER.GT. 2) THEN
!!   ! fixme - this won't scale linearly 
!!   CALL LS_DZERO(GAO(1,1,11),10*NactBAS*NVCLEN)
!!END IF
!!#endif
IADR = 0
DO IBL = 1, NBLCNT
   DO ISHELA = IBLCKS(1,IBL),IBLCKS(2,IBL)
!      IORB =   NSTART(ISHELA)+1  !ORBITALINDEX NOT USED
      SHELLANGMOMA = SHELLANGMOM(ISHELA)
      KHKTA  = 2*(SHELLANGMOMA-1)+1  
      KCKTA  = SHELLANGMOMA*(SHELLANGMOMA+1)/2  
      SPVAL= KCKTA*KCKTA
      JSTA = PRIEXPSTART(ISHELA)
      CENX = CENT(1,ISHELA)!+ ORIG(1)
      CENY = CENT(2,ISHELA)!+ ORIG(2)
      CENZ = CENT(3,ISHELA)!+ ORIG(3)
!!$      IF (NDER.GT. 1.OR.SHELLANGMOMA .GE. 3.OR.DOLND) THEN !d orbitals or deriv
!!$         call mem_dft_dealloc(LVALUE)
!!$         call mem_dft_dealloc(MVALUE)
!!$         call mem_dft_dealloc(NVALUE)
!!$         call mem_dft_alloc(LVALUE,KCKTA*KCKTA)
!!$         call mem_dft_alloc(MVALUE,KCKTA*KCKTA)
!!$         call mem_dft_alloc(NVALUE,KCKTA*KCKTA)
!!$         SPVAL=KCKTA*KCKTA
!!$         IJ=0
!!$         DO I = 1,KCKTA
!!$            DO J = 1,I
!!$               IJ=IJ+1
!!$               LVALUE(IJ)=SHELLANGMOMA-I
!!$               MVALUE(IJ)=I-J
!!$               NVALUE(IJ)=J-1
!!$            ENDDO
!!$         ENDDO
!!$      ENDIF
      DO I=1,NVCLEN !gridpoints
         PA(1,i) = COOR(1,i)-CENX 
         PA(2,i) = COOR(2,i)-CENY
         PA(3,i) = COOR(3,i)-CENZ
      END DO
      DO I=1,NVCLEN !gridpoints
         PA2(i) = PA(1,i)*PA(1,i)+PA(2,i)*PA(2,i)+PA(3,i)*PA(3,i)
      END DO
      IF (NDER.EQ. 0) THEN
         CALL II_BLGETGAO(LUPRI,NVCLEN,NactBAS,NTYPSO,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
              & SPHMAT(SPINDEX(SHELLANGMOMA):SPINDEX(SHELLANGMOMA)+KHKTA*KCKTA-1),&
              & GAO,IADR,PA,PA2,DFTHRI,CC(CCINDEX(ISHELA)),&
              & CCSTART(CCINDEX(ISHELA)),PRIEXP,LVALUE(1:SPVAL,SHELLANGMOMA),&
              & MVALUE(1:SPVAL,SHELLANGMOMA),NVALUE(1:SPVAL,SHELLANGMOMA),SPVAL,Spherical)
      ELSE IF (NDER.GT. 0) THEN
         CALL II_BLGETGA1(LUPRI,NVCLEN,NactBAS,NTYPSO,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
              & SPHMAT(SPINDEX(SHELLANGMOMA):SPINDEX(SHELLANGMOMA)+KHKTA*KCKTA-1),&
              & GAO,IADR,PA,PA2,DFTHRI,CC(CCINDEX(ISHELA)),&
              & CCSTART(CCINDEX(ISHELA)),PRIEXP,LVALUE(1:SPVAL,SHELLANGMOMA),&
              & MVALUE(1:SPVAL,SHELLANGMOMA),NVALUE(1:SPVAL,SHELLANGMOMA),SPVAL,Spherical)
         IF (NDER.GT. 1) THEN
            CALL II_BLGETGA2(LUPRI,NVCLEN,NactBAS,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
              & SPHMAT(SPINDEX(SHELLANGMOMA):SPINDEX(SHELLANGMOMA)+KHKTA*KCKTA-1),&
              & GAO(:,:,5),GAO(:,:,6),&
              & GAO(:,:,7),GAO(:,:,8),&
              & GAO(:,:,9),GAO(:,:,10),IADR,&
              & PA,PA2,DFTHRI,CC(CCINDEX(ISHELA)),CCSTART(CCINDEX(ISHELA)),PRIEXP,&
              & LVALUE(1:SPVAL,SHELLANGMOMA),MVALUE(1:SPVAL,SHELLANGMOMA),&
              & NVALUE(1:SPVAL,SHELLANGMOMA),SPVAL,Spherical)
            IF (NDER.GT. 2) THEN
!               CALL II_BLGETGA3(NVCLEN,&
!               &       GAO(1,IADR,11),GAO(1,IADR,12),GAO(1,IADR,13),&
!               &       GAO(1,IADR,14),GAO(1,IADR,15),GAO(1,IADR,16),&
!               &       GAO(1,IADR,17),GAO(1,IADR,18),GAO(1,IADR,19),&
!               &       GAO(1,IADR,20),SPHMAT(SPINDEX(SHELLANGMOMA)),PA,PA2,DFTHRI)
            END IF
         END IF
      END IF
      IF (DOLND) THEN
         CALL II_BLGETGB1(NVCLEN,NactBAS,IADR,KHKTA,GAO(:,:,1),PA,&
         &              GAO(:,:,NSOB),GAO(:,:,NSOB+1),GAO(:,:,NSOB+2))
         IF (DOGGA) THEN
            DO I = 1, 3
               CALL II_BLGETGB1(NVCLEN,NactBAS,IADR,KHKTA,GAO(:,:,1+I),PA,&
               & GAO(:,:,NSOB+2+I),GAO(:,:,NSOB+5+I),GAO(:,:,NSOB+8+I))
            END DO
         END IF
      END IF
!!$      IF (NDER.GT. 1.OR.SHELLANGMOMA .GE. 3.OR.DOLND) THEN !d orbitals or deriv      
!!$         call mem_dft_dealloc(LVALUE)
!!$         call mem_dft_dealloc(MVALUE)
!!$         call mem_dft_dealloc(NVALUE)
!!$         call mem_dft_alloc(LVALUE,1)
!!$         call mem_dft_alloc(MVALUE,1)
!!$         call mem_dft_alloc(NVALUE,1)
!!$         SPVAL=1
!!$         LVALUE(1)=0E0_realk
!!$         MVALUE(1)=0E0_realk
!!$         NVALUE(1)=0E0_realk
!!$      ENDIF
      IF(Spherical)THEN
         IADR = IADR + KHKTA !UPDATE ACTIVEORBITALINDEX
      ELSE
         IADR = IADR + KCKTA !UPDATE ACTIVEORBITALINDEX
      ENDIF
   END DO
END DO
!call mem_dft_dealloc(LVALUE)
!call mem_dft_dealloc(MVALUE)
!call mem_dft_dealloc(NVALUE)

!write (lupri,*) 'integrals from II_BLGETSOS',nvclen,NactBAS
!call output(gao(1,1,1),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!write (lupri,*) 'x integrals from II_BLGETSOS'
!call output(gao(1,1,2),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!write (lupri,*) 'y integrals from II_BLGETSOS'
!call output(gao(1,1,3),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!write (lupri,*) 'z integrals from II_BLGETSOS'
!call output(gao(1,1,4),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!if (nder.eq. 2) then
!   WRITE(lupri,*)'COOR(1,:)',(COOR(1,i),I=1,NVCLEN)
!   WRITE(lupri,*)'COOR(2,:)',(COOR(2,i),I=1,NVCLEN)
!   WRITE(lupri,*)'COOR(3,:)',(COOR(3,i),I=1,NVCLEN)
!   write (lupri,*) ' xx integrals from II_BLGETSOS '
!   call output(gao(1,1,5),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' xy integrals from II_BLGETSOS '
!   call output(gao(1,1,6),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' xz integrals from II_BLGETSOS '
!   call output(gao(1,1,7),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' yy integrals from II_BLGETSOS '
!   call output(gao(1,1,8),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' yz integrals from II_BLGETSOS '
!   call output(gao(1,1,9),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!   write (lupri,*) ' zz integrals from II_BLGETSOS '
!   call output(gao(1,1,10),1,nvclen,1,NactBAS,nvclen,NactBAS,1,lupri)
!end if
call mem_dft_dealloc(PA)
call mem_dft_dealloc(PA2)

END SUBROUTINE II_BLGETSOS

!> \brief determines the rho for a LDA type calc
!> \author T. Kjaergaard
!> \date 2010
!>
!>     computes  <o|dmat|o'>  i.e rho_a where dmat is a density matrix.
!>
SUBROUTINE II_GETRHO_BLOCKED_LDA(LUPRI,DMAT,NactBAS,GAO,NTYPSO,NBLOCKS,MAXNSHELL,&
     & BLOCKS,NVCLEN,NDMAT,RHO,RHOTHR,MXBLLEN,DRED,GAOGMX,GAOMAX,GAORED,TMP,maxNactBAST)
IMPLICIT NONE
INTEGER,intent(in)     :: NactBAS,NVCLEN,NBLOCKS,NTYPSO,LUPRI,ndmat,MXBLLEN,MAXNSHELL
REAL(REALK),intent(in) :: DMAT(NactBAS,NactBAS,ndmat),GAO(NVCLEN,NactBAS,NTYPSO),RHOTHR
INTEGER,intent(in)     :: BLOCKS(2,MAXNSHELL),maxNactBAST
real(realk),intent(inout) :: RHO(MXBLLEN,NDMAT)
!tmp
REAL(REALK),intent(inout) :: TMP(NVCLEN,NactBAS)
REAL(REALK),intent(inout) :: GAOGMX(maxNactBAST),GAOMAX
REAL(REALK),intent(inout) :: GAORED(NVCLEN,NactBAS)
REAL(REALK),intent(inout) :: DRED(Nactbas*Nactbas)
!
INTEGER     :: IBL,ISTART,IBLEN,JBL,JSTART,JBLEN
INTEGER     :: IDX,JTOP,JDX,K,I,J,IEND,JEND
REAL(REALK) :: DMAX,GAOMAXTMP
INTEGER     :: NRED,JRED,IRED,idmat,offset
INTEGER,pointer :: INXRED(:)
call mem_dft_alloc(INXRED,NactBAS)
! Set up maximum Gaussian AO elements
GAOMAX = 0.0E0_realk
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      GAOMAXTMP = 0.0E0_realk
      DO K = 1, NVCLEN
         GAOMAXTMP = MAX(GAOMAXTMP,ABS(GAO(K,I,1)))
      ENDDO
      GAOGMX(I) = GAOMAXTMP
      GAOMAX = MAX(GAOMAX,GAOMAXTMP)
   ENDDO
ENDDO

IF(NDMAT.GT.1)THEN
   ! Set up maximum density-matrix elements
   DMAX = 0.0E0_realk
   DO IBL=1, NBLOCKS
    ISTART = BLOCKS(1,IBL)
    IEND = BLOCKS(2,IBL)
    DO JBL=1, IBL  !assume symmetric density matrix
     JSTART = BLOCKS(1,JBL)
     JEND = BLOCKS(2,JBL)
     DO IDMAT=1,NDMAT
      DO J = JSTART, JEND
       DO I = ISTART, IEND
        DMAX = MAX(DMAX,ABS(DMAT(I,J,IDMAT)))
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
ELSE
   ! Set up maximum density-matrix elements
   DMAX = 0.0E0_realk
   DO IBL=1, NBLOCKS
    ISTART = BLOCKS(1,IBL)
    IEND = BLOCKS(2,IBL)
    DO JBL=1, IBL  !assume symmetric density matrix
     JSTART = BLOCKS(1,JBL)
     JEND = BLOCKS(2,JBL)
     DO J = JSTART, JEND
      DO I = ISTART, IEND
       DMAX = MAX(DMAX,ABS(DMAT(I,J,1)))
      ENDDO
     ENDDO
    ENDDO
   ENDDO
ENDIF
   !       Set up reduced Gaussian AO's
NRED = 0
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
         NRED = NRED + 1
         INXRED(NRED) = I
         DO K = 1, NVCLEN
            GAORED(K,NRED) = GAO(K,I,1)
         ENDDO
      ENDIF
   ENDDO
ENDDO

IF (NRED.GT. 0) THEN
   DO IDMAT=1,NDMAT
      ! Set up reduced density-matrix
      DO JRED=1,NRED
         J = INXRED(JRED)
         offset = (JRED-1)*NRED
         DO IRED=1,NRED
            I = INXRED(IRED)
            DRED(IRED+offset) = DMAT(I,J,IDMAT)
            ! DRED(IRED,JRED) = DMAT(I,J,IDMAT)
         ENDDO
      ENDDO
      ! First half-contraction of Gaussian AO with density-matrix
      IF(XCintNoOMP)THEN
         CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED,NVCLEN,&
              &     DRED(1:NRED*NRED),NRED,0.0E0_realk,TMP,NVCLEN)
      ELSE !Use Thread Safe version 
         CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED,NVCLEN,&
              &     DRED(1:NRED*NRED),NRED,0.0E0_realk,TMP,NVCLEN)
      ENDIF
      ! Second half-contraction
      DO K = 1, NVCLEN
         RHO(K,IDMAT)= GAORED(K,1)*TMP(K,1)
      END DO
      DO I = 2, NRED
         DO K = 1, NVCLEN
            RHO(K,IDMAT)    = RHO(K,IDMAT) + GAORED(K,I)*TMP(K,I)
         END DO
      END DO
      ! Hack Severeal functionals does not handle a zero density - so
      ! we set these values explicitly to some small value.
      ! Should instead skip these contributions.
      DO K = 1, NVCLEN
         IF (ABS(RHO(K,IDMAT)).LE. 1.0E-20_realk) RHO(K,IDMAT) = 1.0E-20_realk
      END DO
   ENDDO
ELSE
   DO IDMAT=1,NDMAT
      DO K = 1, NVCLEN
         RHO(K,IDMAT) = 1.0E-20_realk
      END DO
   ENDDO
ENDIF
call mem_dft_dealloc(INXRED)

END SUBROUTINE II_GETRHO_BLOCKED_LDA

!> \brief determines the rho for a GGA type calc
!> \author T. Kjaergaard
!> \date 2010
!>
!>     computes  <o|dmat|o'>
!>     i.e rho_a where dmat is a density matrix (it can be a total
!>     density end then one will get total density, or it can be an
!>     alpha/beta density.
!>     assert(NTYPSO>=NRHO)
!>
SUBROUTINE II_GETRHO_BLOCKED_GGA(LUPRI,DMAT,NactBAS,GAO,NTYPSO,NBLOCKS,MAXNSHELL,BLOCKS,&
     & NVCLEN,NDMAT,RHO,GRAD,RHOTHR,MXBLLEN,DRED,GAOGMX,GAOMAX,GAORED,TMP,maxNactBAST)
IMPLICIT NONE
INTEGER,intent(in)     :: NactBAS,NVCLEN,NBLOCKS,NTYPSO,LUPRI,ndmat,MXBLLEN,MAXNSHELL!,nbast
REAL(REALK),intent(in) :: DMAT(NactBAS,NactBAS,ndmat), GAO(NVCLEN,NactBAS,NTYPSO)
INTEGER,intent(in)     :: BLOCKS(2,MAXNSHELL),MaxNactBast
REAL(REALK),intent(inout) :: RHO(MXBLLEN,NDMAT), GRAD(3,MXBLLEN,NDMAT)
REAL(REALK),intent(inout) :: DRED(NactBAS*NactBAS)
REAL(REALK),intent(inout) :: GAOGMX(maxNactBAST),GAOMAX
REAL(REALK),intent(inout) :: GAORED(NVCLEN,NactBAS,4)
REAL(REALK),intent(inout) :: TMP(NVCLEN,NactBAS)
!
REAL(REALK) :: DMAX,RHOTHR,GAOMAXTMP
INTEGER     :: INXRED(NactBAS)
INTEGER     :: IBL,ISTART,IBLEN,JBL,JSTART,JBLEN,IDX,JTOP,JDX
INTEGER     :: K,I,J,NRED,JRED,IRED,idmat,offset,JEND,IEND

! Set up maximum Gaussian AO elements
GAOMAX = 0.0E0_realk
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      GAOMAXTMP = 0.0E0_realk
      DO K = 1, NVCLEN
         GAOMAXTMP = MAX(GAOMAXTMP,ABS(GAO(K,I,1)))
      ENDDO
      GAOMAX = MAX(GAOMAX,GAOMAXTMP)
      DO J=2,4
         DO K = 1, NVCLEN
            GAOMAXTMP = MAX(GAOMAXTMP,ABS(GAO(K,I,J)))
         ENDDO
      ENDDO
      GAOGMX(I) = GAOMAXTMP
   ENDDO
ENDDO
! Set up maximum density-matrix elements
DMAX = 0.0E0_realk
DO IDMAT=1,NDMAT
 DO IBL=1, NBLOCKS
  DO JBL=1, IBL  !assume symmetric density matrix
    DO J = BLOCKS(1,JBL),BLOCKS(2,JBL)
     DO I = BLOCKS(1,IBL),BLOCKS(2,IBL)
      DMAX = MAX(DMAX,ABS(DMAT(I,J,IDMAT)))
     ENDDO
    ENDDO
  ENDDO
 ENDDO
ENDDO
!Set up reduced Gaussian AO's
NRED = 0
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
         NRED = NRED + 1
         INXRED(NRED) = I
         DO K = 1, NVCLEN
            GAORED(K,NRED,1) = GAO(K,I,1)
            GAORED(K,NRED,2) = GAO(K,I,2)
            GAORED(K,NRED,3) = GAO(K,I,3)
            GAORED(K,NRED,4) = GAO(K,I,4)
         ENDDO
      ENDIF
   ENDDO
ENDDO
IF (NRED.GT. 0) THEN
   DO IDMAT=1,NDMAT
      IF(NRED.EQ.NactBAS)THEN 
         !no screening
         !First half-contraction of Gaussian AO with density-matrix
         IF(XCintNoOMP)THEN
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED,NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP,NVCLEN)
         ELSE !Use Thread Safe version 
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED,NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP,NVCLEN)
         ENDIF
      ELSE
         !Set up reduced density-matrix
         DO JRED=1,NRED
            J = INXRED(JRED)
            offset = (JRED-1)*NRED
            DO IRED=1,NRED
               I = INXRED(IRED)
               DRED(IRED+offset) = DMAT(I,J,IDMAT)
            ENDDO
         ENDDO
         !First half-contraction of Gaussian AO with density-matrix
         IF(XCintNoOMP)THEN
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED,NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP,NVCLEN)
         ELSE !Use Thread Safe version 
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED,NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP,NVCLEN)
         ENDIF
      ENDIF
      !Second half-contraction
      DO K = 1, NVCLEN
         RHO(K,IDMAT)    = GAORED(K,1,1)*TMP(K,1)
         GRAD(1,K,IDMAT) = 2*GAORED(K,1,2)*TMP(K,1)
         GRAD(2,K,IDMAT) = 2*GAORED(K,1,3)*TMP(K,1)
         GRAD(3,K,IDMAT) = 2*GAORED(K,1,4)*TMP(K,1)
      END DO
      DO I = 2, NRED
         DO K = 1, NVCLEN
            RHO(K,IDMAT)    = RHO(K,IDMAT)    +   GAORED(K,I,1)*TMP(K,I)
            GRAD(1,K,IDMAT) = GRAD(1,K,IDMAT) + 2*GAORED(K,I,2)*TMP(K,I)
            GRAD(2,K,IDMAT) = GRAD(2,K,IDMAT) + 2*GAORED(K,I,3)*TMP(K,I)
            GRAD(3,K,IDMAT) = GRAD(3,K,IDMAT) + 2*GAORED(K,I,4)*TMP(K,I)
         END DO
      END DO
   ENDDO
   !Hack Severeal functionals does not handle a zero density - so
   !     we set these values explicitly to some small value.
   !     Should instead skip these contributions.
   DO IDMAT=1,NDMAT   
      DO K = 1, NVCLEN
         IF (ABS(RHO(K,IDMAT)).LE. 1.0E-20_realk) RHO(K,IDMAT) = 1.0E-20_realk
         IF (ABS(GRAD(1,K,IDMAT)).LE. 1.0E-20_realk) GRAD(1,K,IDMAT) = 1.0E-20_realk
         IF (ABS(GRAD(2,K,IDMAT)).LE. 1.0E-20_realk) GRAD(2,K,IDMAT) = 1.0E-20_realk
         IF (ABS(GRAD(3,K,IDMAT)).LE. 1.0E-20_realk) GRAD(3,K,IDMAT) = 1.0E-20_realk
      END DO
   ENDDO
ELSE
   !Hack Severeal functionals does not handle a zero density - so
   !     we set these values explicitly to some small value.
   !     Should instead skip these contributions.
   DO IDMAT=1,NDMAT
      DO K = 1, NVCLEN
         RHO(K,IDMAT) = 1.0E-20_realk
         GRAD(1,K,IDMAT) = 1.0E-20_realk
         GRAD(2,K,IDMAT) = 1.0E-20_realk
         GRAD(3,K,IDMAT) = 1.0E-20_realk
      END DO
   ENDDO
ENDIF
END SUBROUTINE II_GETRHO_BLOCKED_GGA

!> \brief determines the rho for a GGA type calc
!> \author S. _eine
!> \date 2013-03-27
!>
!>     computes       rho(r_i) = sum_ab X_a(r_i) X_b(r_i) D_ab
!>               grad rho(r_i) = sum_ab (grad X_a(r_i) X_b(r_i) + X_a(r_i) grad  X_b(r_i)) D_ab
!>                    tau(r_i) = sum_ab grad X_a(r_i) grad X_b(r_i) D_ab
!>     i.e rho(r_i) where D_ab is a density matrix (it can be a total
!>     density end then one will get total density, or it can be an
!>     alpha/beta density.
!>     assert(NTYPSO>=NRHO)
!>
SUBROUTINE II_GETRHO_BLOCKED_META(LUPRI,DMAT,NactBAS,GAO,NTYPSO,NBLOCKS,MAXNSHELL,BLOCKS,&
     & NVCLEN,NDMAT,RHO,GRAD,TAU,RHOTHR,MXBLLEN,DRED,GAOGMX,GAOMAX,GAORED,TMP,maxNactBAST)
IMPLICIT NONE
INTEGER,intent(in)     :: NactBAS,NVCLEN,NBLOCKS,NTYPSO,LUPRI,ndmat,MXBLLEN,MAXNSHELL!,nbast
REAL(REALK),intent(in) :: DMAT(NactBAS,NactBAS,ndmat), GAO(NVCLEN,NactBAS,NTYPSO)
INTEGER,intent(in)     :: BLOCKS(2,MAXNSHELL),MaxNactBast
REAL(REALK),intent(inout) :: RHO(MXBLLEN,NDMAT), GRAD(3,MXBLLEN,NDMAT), TAU(MXBLLEN,NDMAT)
REAL(REALK),intent(inout) :: DRED(NactBAS*NactBAS)
REAL(REALK),intent(inout) :: GAOGMX(maxNactBAST),GAOMAX
REAL(REALK),intent(inout) :: GAORED(NVCLEN,NactBAS,4)
REAL(REALK),intent(inout) :: TMP(NVCLEN,NactBAS,4)
!
REAL(REALK) :: DMAX,RHOTHR,GAOMAXTMP
INTEGER :: INXRED(NactBAS)
INTEGER     :: IBL,ISTART,IBLEN,JBL,JSTART,JBLEN,IDX,JTOP,JDX
INTEGER     :: K,I,J,NRED,JRED,IRED,idmat,offset,JEND,IEND

! Set up maximum Gaussian AO elements
GAOMAX = 0.0E0_realk
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      GAOMAXTMP = 0.0E0_realk
      DO K = 1, NVCLEN
         GAOMAXTMP = MAX(GAOMAXTMP,ABS(GAO(K,I,1)))
      ENDDO
      GAOMAX = MAX(GAOMAX,GAOMAXTMP)
      DO J=2,4
         DO K = 1, NVCLEN
            GAOMAXTMP = MAX(GAOMAXTMP,ABS(GAO(K,I,J)))
         ENDDO
      ENDDO
      GAOGMX(I) = GAOMAXTMP
   ENDDO
ENDDO

! Set up maximum density-matrix elements
DMAX = 0.0E0_realk
DO IDMAT=1,NDMAT
 DO JBL=1, NBLOCKS
  DO IBL=1, NBLOCKS
    DO J = BLOCKS(1,JBL),BLOCKS(2,JBL)
     DO I = BLOCKS(1,IBL),BLOCKS(2,IBL)
      DMAX = MAX(DMAX,ABS(DMAT(I,J,IDMAT)))
     ENDDO
    ENDDO
  ENDDO
 ENDDO
ENDDO
!Set up reduced Gaussian AO's
NRED = 0
DO IBL=1, NBLOCKS
   DO I = BLOCKS(1,IBL), BLOCKS(2,IBL)
      IF (GAOGMX(I)*GAOMAX*DMAX.GT.RHOTHR) THEN
         NRED = NRED + 1
         INXRED(NRED) = I
         DO K = 1, NVCLEN
            GAORED(K,NRED,1) = GAO(K,I,1)
            GAORED(K,NRED,2) = GAO(K,I,2)
            GAORED(K,NRED,3) = GAO(K,I,3)
            GAORED(K,NRED,4) = GAO(K,I,4)
         ENDDO
      ENDIF
   ENDDO
ENDDO
IF (NRED.GT. 0) THEN
   DO IDMAT=1,NDMAT
      IF(NRED.EQ.NactBAS)THEN 
         !no screening
         !First half-contraction of Gaussian AO with density-matrix
         IF(XCintNoOMP)THEN
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,1),NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP(:,:,1),NVCLEN)
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,2),NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP(:,:,2),NVCLEN)
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,3),NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP(:,:,3),NVCLEN)
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,4),NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP(:,:,4),NVCLEN)
         ELSE !Use Thread Safe version 
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,1),NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP(:,:,1),NVCLEN)
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,2),NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP(:,:,2),NVCLEN)
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,3),NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP(:,:,3),NVCLEN)
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,4),NVCLEN,&
                 &     DMAT(1:NRED,1:NRED,IDMAT),NRED,0.0E0_realk,TMP(:,:,4),NVCLEN)
         ENDIF
      ELSE
         !Set up reduced density-matrix
         DO JRED=1,NRED
            J = INXRED(JRED)
            offset = (JRED-1)*NRED
            DO IRED=1,NRED
               I = INXRED(IRED)
               DRED(IRED+offset) = DMAT(I,J,IDMAT)
            ENDDO
         ENDDO
         !First half-contraction of Gaussian AO with density-matrix
         IF(XCintNoOMP)THEN
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,1),NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP(:,:,1),NVCLEN)
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,2),NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP(:,:,2),NVCLEN)
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,3),NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP(:,:,3),NVCLEN)
            CALL DGEMM("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,4),NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP(:,:,4),NVCLEN)
         ELSE !Use Thread Safe version 
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,1),NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP(:,:,1),NVCLEN)
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,2),NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP(:,:,2),NVCLEN)
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,3),NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP(:,:,3),NVCLEN)
            CALL DGEMM_TS("N","N",NVCLEN,NRED,NRED,1E0_realk,GAORED(:,:,4),NVCLEN,&
                 &     DRED,NRED,0.0E0_realk,TMP(:,:,4),NVCLEN)
         ENDIF
      ENDIF
      !Second half-contraction
      DO K = 1, NVCLEN
         RHO(K,IDMAT)    = GAORED(K,1,1)*TMP(K,1,1)
         GRAD(1,K,IDMAT) = 2*GAORED(K,1,2)*TMP(K,1,1)
         GRAD(2,K,IDMAT) = 2*GAORED(K,1,3)*TMP(K,1,1)
         GRAD(3,K,IDMAT) = 2*GAORED(K,1,4)*TMP(K,1,1)
         TAU(K,IDMAT)    = GAORED(K,1,2)*TMP(K,1,2)+GAORED(K,1,3)*TMP(K,1,3)+GAORED(K,1,4)*TMP(K,1,4)
      END DO
      DO I = 2, NRED
         DO K = 1, NVCLEN
            RHO(K,IDMAT)    = RHO(K,IDMAT)    +   GAORED(K,I,1)*TMP(K,I,1)
            GRAD(1,K,IDMAT) = GRAD(1,K,IDMAT) + 2*GAORED(K,I,2)*TMP(K,I,1)
            GRAD(2,K,IDMAT) = GRAD(2,K,IDMAT) + 2*GAORED(K,I,3)*TMP(K,I,1)
            GRAD(3,K,IDMAT) = GRAD(3,K,IDMAT) + 2*GAORED(K,I,4)*TMP(K,I,1)
            TAU(K,IDMAT)    = TAU(K,IDMAT)    + GAORED(K,I,2)*TMP(K,I,2) + GAORED(K,I,3)*TMP(K,I,3) + GAORED(K,I,4)*TMP(K,I,4)
         END DO
      END DO
   ENDDO
   !Hack Severeal functionals does not handle a zero density - so
   !     we set these values explicitly to some small value.
   !     Should instead skip these contributions.
   DO IDMAT=1,NDMAT   
      DO K = 1, NVCLEN
         IF (ABS(RHO(K,IDMAT)).LE. 1.0E-20_realk) RHO(K,IDMAT) = 1.0E-20_realk
         IF (ABS(GRAD(1,K,IDMAT)).LE. 1.0E-20_realk) GRAD(1,K,IDMAT) = 1.0E-20_realk
         IF (ABS(GRAD(2,K,IDMAT)).LE. 1.0E-20_realk) GRAD(2,K,IDMAT) = 1.0E-20_realk
         IF (ABS(GRAD(3,K,IDMAT)).LE. 1.0E-20_realk) GRAD(3,K,IDMAT) = 1.0E-20_realk
         IF (ABS(TAU(K,IDMAT)).LE. 1.0E-20_realk) TAU(K,IDMAT) = 1.0E-20_realk
      END DO
   ENDDO
ELSE
   !Hack Severeal functionals does not handle a zero density - so
   !     we set these values explicitly to some small value.
   !     Should instead skip these contributions.
   DO IDMAT=1,NDMAT
      DO K = 1, NVCLEN
         RHO(K,IDMAT) = 1.0E-20_realk
         GRAD(1,K,IDMAT) = 1.0E-20_realk
         GRAD(2,K,IDMAT) = 1.0E-20_realk
         GRAD(3,K,IDMAT) = 1.0E-20_realk
         TAU(K,IDMAT) = 1.0E-20_realk
      END DO
   ENDDO
ENDIF
END SUBROUTINE II_GETRHO_BLOCKED_META

!> \brief evaluates the pure GAO(gaussian atomic orbitals) , so no derivatives
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGAO(LUPRI,NVCLEN,NactBAS,NTYPSO,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
                     & CSP,GAO,IADR,PA,PA2,DFTHRI,CC,CCSTART,PRIEXP,LVALUE,MVALUE,&
                     &NVALUE,SPVAL,Spherical)
IMPLICIT NONE
INTEGER,intent(in)        :: NVCLEN,NactBAS,MXPRIM,KHKTA,JSTA,KCKTA,LUPRI,SHELLANGMOMA
INTEGER,intent(in)        :: SPVAL,CCSTART,IADR,NTYPSO,LVALUE(SPVAL),MVALUE(SPVAL)
INTEGER,intent(in)        :: NVALUE(SPVAL)
REAL(REALK),PARAMETER     :: D0 = 0.0E0_realk, D1 = 1.0E0_realk
REAL(REALK),intent(inout) :: GAO(NVCLEN,NactBAS,NTYPSO)
REAL(REALK),intent(in)    :: PA(3,NVCLEN),PA2(NVCLEN),DFTHRI
REAL(REALK),intent(in)    :: CSP(KHKTA,KCKTA),PRIEXP(MXPRIM)
TYPE(LSMATRIX),intent(in) :: CC
!
REAL(REALK)    :: GAZ,GAXX,GAX,GAY,GAYY,GMAX,PA_1,PA_2,PA_3
REAL(REALK)    :: PRICCFVAL,PRIEXPVAL
INTEGER        :: I,J,TEMPI,LVALJ, MVALJ, NVALJ, LVALI, MVALI, NVALI,K
INTEGER        :: ISTART1,ISTART2,ISTART3,ISTART4,ISTART5,ISTART6,ISTART
LOGICAl        :: Spherical
REAL(REALK),pointer  :: GA(:),CINT(:),CINT2(:,:)

call mem_dft_alloc(GA,NVCLEN)

!     loop over primitives
ISTART=IADR
ISTART1=ISTART+1
ISTART2=ISTART+2
ISTART3=ISTART+3
ISTART4=ISTART+4
ISTART5=ISTART+5
ISTART6=ISTART+6
!CALL LS_DZERO(GA,NVCLEN)
I=1
J = JSTA+CCSTART-1+I
PRICCFVAL = CC%elms(I)
PRIEXPVAL = -PRIEXP(J)
DO K = 1, NVCLEN
   GA(K) = PRICCFVAL*EXP(PRIEXPVAL*PA2(K))
END DO
DO I = 2,CC%nrow
   J = JSTA+CCSTART-1+I
   PRICCFVAL = CC%elms(I)
   PRIEXPVAL = -PRIEXP(J)
   DO K = 1, NVCLEN
      GA(K) = GA(K) + PRICCFVAL*EXP(PRIEXPVAL*PA2(K))
   END DO
ENDDO

! screening based on the maximum GA value for the baych of points: GMAX \Andreas Krapp
GMAX = D0
DO K = 1, NVCLEN
   GMAX = MAX(GMAX,ABS(GA(K)))
END DO

!     
!     contracted orbitals
!
IF (SHELLANGMOMA .EQ. 1) THEN      !s orbitals
   IF (GMAX .GT. DFTHRI) THEN 
      DO K = 1, NVCLEN
         GAO(K,ISTART1,1) = GA(K)
      END DO
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1,1),NVCLEN)
   END IF
ELSEIF (SHELLANGMOMA .EQ. 2) THEN !p orbitals
   IF (GMAX .GT. DFTHRI) THEN
      DO K = 1, NVCLEN
         GAO(K,ISTART1,1) = PA(1,K)*GA(K)
         GAO(K,ISTART2,1) = PA(2,K)*GA(K)
         GAO(K,ISTART3,1) = PA(3,K)*GA(K)
      END DO
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1,1),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,1),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,1),NVCLEN)
   END IF
ELSEIF (SHELLANGMOMA .EQ. 3) THEN !d orbitals
   IF (Spherical) THEN
      IF (GMAX .GT. DFTHRI) THEN
         DO K = 1, NVCLEN
            GAO(K,ISTART1,1) = PA(2,K)*PA(1,K)*GA(K)
            GAO(K,ISTART2,1) = PA(2,K)*PA(3,K)*GA(K)
            GAO(K,ISTART3,1) = -0.288675134594813E0_realk*(PA(1,K)*PA(1,K) + PA(2,K)*PA(2,K))*GA(K)&
            &                  + 0.577350269189626E0_realk*PA(3,K)*PA(3,K)*GA(K)
            GAO(K,ISTART4,1) = PA(1,K)*PA(3,K)*GA(K)
            GAO(K,ISTART5,1) = 0.5E0_realk*(PA(1,K)*PA(1,K) - PA(2,K)*PA(2,K))*GA(K)
         END DO
      ELSE 
         CALL LS_DZERO(GAO(1,ISTART1,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART2,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART3,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART4,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART5,1),NVCLEN)
      END IF
   ELSE
      IF (GMAX .GT. DFTHRI) THEN
         DO K = 1, NVCLEN
            PA_1 = PA(1,K)
            PA_2 = PA(2,K)
            PA_3 = PA(3,K)
            GAX  = PA_1*GA(K)
            GAY  = PA_2*GA(K)
            GAZ  = PA_3*GA(K)
            GAO(K,ISTART1,1) = PA_1*GAX
            GAO(K,ISTART2,1) = PA_2*GAX
            GAO(K,ISTART3,1) = PA_3*GAX
            GAO(K,ISTART4,1) = PA_2*GAY
            GAO(K,ISTART5,1) = PA_3*GAY 
            GAO(K,ISTART6,1) = PA_3*GAZ
         END DO
      ELSE
         CALL LS_DZERO(GAO(1,ISTART1,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART2,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART3,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART4,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART5,1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART6,1),NVCLEN)
      END IF
   END IF
ELSEIF (SHELLANGMOMA .EQ. 4) THEN !f orbitals
   IF (Spherical) THEN
      !KCKTA=10
      IF (GMAX .GT. DFTHRI) THEN
         call mem_dft_alloc(CINT2,NVCLEN,KCKTA)
         DO J = 1, KCKTA
            LVALJ = LVALUE(J)
            MVALJ = MVALUE(J)
            NVALJ = NVALUE(J)
            DO K = 1, NVCLEN
               CINT2(K,J) = (PA(1,K)**LVALJ)*(PA(2,K)**MVALJ)&
                    &                 *(PA(3,K)**NVALJ)*GA(K)
            END DO
         ENDDO
         DO K = 1, NVCLEN
            GAO(K,ISTART+1,1) = &
                 & + 0.612372435695794E0_realk*CINT2(K,2) &
                 & - 0.204124145231932E0_realk*CINT2(K,7) 
         ENDDO
         DO K = 1, NVCLEN
            GAO(K,ISTART+2,1) = CINT2(K,5)
         ENDDO
         DO K = 1, NVCLEN
            GAO(K,ISTART+3,1) = &
                 & -0.158113883008419E0_realk*(CINT2(K,2)+CINT2(K,7))&
                 & +0.632455532033676E0_realk*CINT2(K,9)
         ENDDO
         DO K = 1, NVCLEN
            GAO(K,ISTART+4,1) = &
                 & -0.387298334620742E0_realk*(CINT2(K,3)+CINT2(K,8))&
                 & +0.258198889747161E0_realk*CINT2(K,10)
         ENDDO
         DO K = 1, NVCLEN
            GAO(K,ISTART+5,1) = &
                 & -0.158113883008419E0_realk*(CINT2(K,1)+CINT2(K,4))&
                 & +0.632455532033676E0_realk*CINT2(K,6)
         END DO
         DO K = 1, NVCLEN
            GAO(K,ISTART+6,1) = &
                 & 0.5E0_realk*(CINT2(K,3)-CINT2(K,8))
         END DO
         DO K = 1, NVCLEN
            GAO(K,ISTART+7,1) = &
                 & +0.204124145231932E0_realk*CINT2(K,1)&
                 & -0.612372435695794E0_realk*CINT2(K,4)
         END DO
         call mem_dft_dealloc(CINT2)
      ELSE
         DO I = 1, KHKTA
            CALL LS_DZERO(GAO(1,ISTART+I,1),NVCLEN)
         END DO
      END IF
   ELSE
      IF (GMAX .GT. DFTHRI) THEN
         DO I = 1, KCKTA
            LVALI = LVALUE(I)
            MVALI = MVALUE(I)
            NVALI = NVALUE(I)
            TEMPI=ISTART+I
            DO K = 1, NVCLEN
               GAO(K,TEMPI,1) = (PA(1,K)**LVALI)*(PA(2,K)**MVALI)&
                    &                 *(PA(3,K)**NVALI)*GA(K)
            END DO
         END DO
      ELSE
         DO I = 1, KCKTA
            CALL LS_DZERO(GAO(1,ISTART+I,1),NVCLEN)
         END DO
      END IF
   ENDIF
ELSE !higher than f orbitals
   IF (Spherical) THEN
      call mem_dft_alloc(CINT,NVCLEN)
      DO I = 1, KHKTA
         CALL LS_DZERO(GAO(1,ISTART+I,1),NVCLEN)
      END DO
      IF (GMAX .GT. DFTHRI) THEN
         DO J = 1, KCKTA
            LVALJ = LVALUE(J)
            MVALJ = MVALUE(J)
            NVALJ = NVALUE(J)
            DO K = 1, NVCLEN
               CINT(K) = (PA(1,K)**LVALJ)*(PA(2,K)**MVALJ)&
                    &                 *(PA(3,K)**NVALJ)*GA(K)
            END DO
            DO I = 1, KHKTA
               TEMPI=ISTART+I
               DO K = 1, NVCLEN
                  GAO(K,TEMPI,1) = GAO(K,TEMPI,1) + CSP(I,J)*CINT(K)
               END DO
            END DO
         END DO
      END IF
      call mem_dft_dealloc(CINT)
   ELSE
      IF (GMAX .GT. DFTHRI) THEN
         DO I = 1, KCKTA
            LVALI = LVALUE(I)
            MVALI = MVALUE(I)
            NVALI = NVALUE(I)
            TEMPI=ISTART+I
            DO K = 1, NVCLEN
               GAO(K,TEMPI,1) = (PA(1,K)**LVALI)*(PA(2,K)**MVALI)&
                    &                 *(PA(3,K)**NVALI)*GA(K)
            END DO
         END DO
      ELSE
         DO I = 1, KCKTA
            CALL LS_DZERO(GAO(1,ISTART+I,1),NVCLEN)
         END DO
      END IF
   END IF
END IF
call mem_dft_dealloc(GA)

END SUBROUTINE II_BLGETGAO

!> \brief evaluates the GAO(gaussian atomic orbitals) + first geo derivatives
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGA1(LUPRI,NVCLEN,NactBAS,NTYPSO,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
                     & CSP,GAO,IADR,PA,PA2,DFTHRI,CC,CCSTART,PRIEXP,LVALUE,MVALUE,&
                     & NVALUE,SPVAL,Spherical)
 
IMPLICIT NONE
INTEGER,intent(in)        :: NVCLEN,NactBAS,MXPRIM,KHKTA,JSTA,KCKTA,LUPRI,SHELLANGMOMA
INTEGER,intent(in)        :: SPVAL,IADR,NTYPSO,LVALUE(SPVAL),MVALUE(SPVAL),NVALUE(SPVAL)
REAL(REALK),PARAMETER     :: D0 = 0.0E0_realk, D1 = 1.0E0_realk,D2 = 2.0E0_realk
REAL(REALK),intent(inout) :: GAO(NVCLEN,NactBAS,NTYPSO) !NTYPSO=4
REAL(REALK),intent(in)    :: PA(3,NVCLEN),PA2(NVCLEN),DFTHRI
REAL(REALK),intent(in)    :: CSP(KHKTA,KCKTA),PRIEXP(MXPRIM)
LOGICAL,intent(in)        :: Spherical
TYPE(LSMATRIX),intent(in) :: CC
!
INTEGER                   :: I,J,K
INTEGER,PARAMETER         :: I0=1, IX=2, IY=3, IZ=4
REAL(REALK)  :: SPHFAC,FAC,TGX,TGY,TGZ
REAL(REALK)  :: PRICCFVAL,PRIEXPVAL,GMAX,PA_1,PA_2,PA_3,GAVAL,DL,DN,DM
INTEGER      :: ICOMPA,L,M,N,CCSTART,TEMPI,ISTART
INTEGER      :: ISTART1,ISTART2,ISTART3,ISTART4,ISTART5
REAL(REALK), pointer :: CAO(:,:), CAOX(:,:),CAOY(:,:),CAOZ(:,:)
REAL(REALK),pointer  :: GA(:),GU(:),GAX(:),GAY(:),GAZ(:), P0(:)
REAL(REALK),pointer  :: FX(:), FY(:), FZ(:),SPH0(:), SPHX(:),SPHY(:), SPHZ(:)
call mem_dft_alloc(GA,NVCLEN)
call mem_dft_alloc(GU,NVCLEN)
!     loop over primitives
ISTART=IADR
ISTART1=ISTART+1
ISTART2=ISTART+2
ISTART3=ISTART+3
ISTART4=ISTART+4
ISTART5=ISTART+5
!CALL LS_DZERO(GA,NVCLEN)
!CALL LS_DZERO(GU,NVCLEN)
I=1
J = JSTA+CCSTART
PRICCFVAL = CC%elms(I)
PRIEXPVAL = PRIEXP(J)
DO K = 1, NVCLEN
   FAC = PRICCFVAL*EXP(-PRIEXPVAL*PA2(K))
   GA(K) = FAC
   GU(K) = - D2*PRIEXPVAL*FAC
END DO
DO I = 2,CC%nrow
   J = JSTA+CCSTART-1+I
   PRICCFVAL = CC%elms(I)
   PRIEXPVAL = PRIEXP(J)
   DO K = 1, NVCLEN
      FAC   = PRICCFVAL*EXP(-PRIEXPVAL*PA2(K))
      GA(K) = GA(K) + FAC
      GU(K) = GU(K) - D2*PRIEXPVAL*FAC
   END DO
END DO

IF (Spherical) THEN
   call mem_dft_alloc(CAO,NVCLEN,KCKTA)
   call mem_dft_alloc(CAOX,NVCLEN,KCKTA)
   call mem_dft_alloc(CAOY,NVCLEN,KCKTA)
   call mem_dft_alloc(CAOZ,NVCLEN,KCKTA)
END IF

!screening based on the maximum GA value for the batch of points: GMAX
GMAX = D0
DO K = 1, NVCLEN
   GMAX = MAX(GMAX,ABS(GA(K)))
END DO

!s orbitals
IF (SHELLANGMOMA .EQ. 1) THEN
   IF (GMAX .GT. DFTHRI) THEN
      DO K = 1, NVCLEN
         GAO(K,ISTART1,I0) = GA(K)
      END DO
      DO K = 1, NVCLEN
         GAO(K,ISTART1,IX) = PA(1,K)*GU(K) !CHANGE TO PA(K,1-3)
      END DO
      DO K = 1, NVCLEN
         GAO(K,ISTART1,IY) = PA(2,K)*GU(K)
      END DO
      DO K = 1, NVCLEN
         GAO(K,ISTART1,IZ) = PA(3,K)*GU(K)
      END DO  
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1,I0),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IX),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IY),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IZ),NVCLEN)
   END IF
!p orbitals
ELSE IF (SHELLANGMOMA .EQ. 2) THEN
   IF (GMAX .GT. DFTHRI) THEN
      DO K = 1, NVCLEN
         PA_1 = PA(1,K)
         PA_2 = PA(2,K)
         PA_3 = PA(3,K)
         GAVAL = GA(K)
         TGX = PA_1*GU(K)     !CHANGE TO TGX(K)
         TGY = PA_2*GU(K)
         TGZ = PA_3*GU(K)
         GAO(K,ISTART1,I0) = PA_1*GAVAL
         GAO(K,ISTART2,I0) = PA_2*GAVAL
         GAO(K,ISTART3,I0) = PA_3*GAVAL
         GAO(K,ISTART1,IX) = PA_1*TGX + GAVAL
         GAO(K,ISTART2,IX) = PA_2*TGX
         GAO(K,ISTART3,IX) = PA_3*TGX
         GAO(K,ISTART1,IY) = PA_1*TGY
         GAO(K,ISTART2,IY) = PA_2*TGY + GAVAL
         GAO(K,ISTART3,IY) = PA_3*TGY
         GAO(K,ISTART1,IZ) = PA_1*TGZ
         GAO(K,ISTART2,IZ) = PA_2*TGZ
         GAO(K,ISTART3,IZ) = PA_3*TGZ + GAVAL
      END DO
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1,I0),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,I0),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,I0),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IX),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,IX),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,IX),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IY),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,IY),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,IY),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART1,IZ),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2,IZ),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3,IZ),NVCLEN)
   END IF
! d and higher orbitals
ELSE
   IF (GMAX.GT.DFTHRI) THEN
      call mem_dft_alloc(GAX,NVCLEN)
      call mem_dft_alloc(GAY,NVCLEN)
      call mem_dft_alloc(GAZ,NVCLEN)
      call mem_dft_alloc(P0,NVCLEN)
      call mem_dft_alloc(FX,NVCLEN)
      call mem_dft_alloc(FY,NVCLEN)
      call mem_dft_alloc(FZ,NVCLEN)
      DO K = 1, NVCLEN
         FX(K) = PA(1,K)*GU(K)
         FY(K) = PA(2,K)*GU(K)
         FZ(K) = PA(3,K)*GU(K)
      END DO
      DO ICOMPA = 1,KCKTA
         L = LVALUE(ICOMPA)
         M = MVALUE(ICOMPA)
         N = NVALUE(ICOMPA)
         DO K = 1, NVCLEN
            P0(K)  = (PA(1,K)**L)*(PA(2,K)**M)*(PA(3,K)**N) ! make different loops for L=0 or M=0 or N=0 etc ??
         END DO
         DO K = 1, NVCLEN
            GAX(K) = FX(K)*P0(K)
            GAY(K) = FY(K)*P0(K)
            GAZ(K) = FZ(K)*P0(K)
            GAVAL = GA(K)
            PA_1  = PA(1,K)
            PA_2  = PA(2,K)
            PA_3  = PA(3,K)
            IF(L.GT. 0) THEN
               GAX(K) = GAX(K)+L*(PA_1**(L-1))*(PA_2**M)*&
                    &  (PA_3**N)*GAVAL
            ENDIF
            IF(M.GT. 0) THEN
               GAY(K) = GAY(K)+M*(PA_1**L)*(PA_2**(M-1))*&
                    &   (PA_3**N)*GAVAL
            ENDIF
            IF(N.GT. 0) THEN
               GAZ(K) = GAZ(K)+N*(PA_1**L)*(PA_2**M)*&
                    &  (PA_3**(N-1))*GAVAL
            ENDIF
         END DO
         IF (Spherical) THEN
            DO K = 1, NVCLEN
               CAO (K, ICOMPA) = GA(K)*P0(K)
               CAOX(K, ICOMPA) = GAX(K)
               CAOY(K, ICOMPA) = GAY(K)
               CAOZ(K, ICOMPA) = GAZ(K)
            END DO
         ELSE
            TEMPI = ISTART+ICOMPA
            DO K = 1, NVCLEN            
               GAO(K,TEMPI,I0) = GA(K)*P0(K)
            END DO
            DO K = 1, NVCLEN
               GAO(K,TEMPI,IX) = GAX(K)
            END DO
            DO K = 1, NVCLEN
               GAO(K,TEMPI,IY) = GAY(K)
            END DO
            DO K = 1, NVCLEN
               GAO(K,TEMPI,IZ) = GAZ(K)
            END DO
         END IF
      END DO
      call mem_dft_dealloc(GAX)
      call mem_dft_dealloc(GAY)
      call mem_dft_dealloc(GAZ)
      IF(Spherical) THEN
         IF (SHELLANGMOMA.EQ. 3) THEN !d orbitals
            DO K = 1, NVCLEN
               GAO(K,ISTART1,I0) = CSP(1,2)*CAO(K,2)
               GAO(K,ISTART2,I0) = CSP(2,5)*CAO(K,5)
               GAO(K,ISTART3,I0) = CSP(3,1)*CAO(K,1) + CSP(3,4)*CAO(K,4)& 
               &                               + CSP(3,6)*CAO(K,6)
               GAO(K,ISTART4,I0) = CSP(4,3)*CAO(K,3)
               GAO(K,ISTART5,I0) = CSP(5,1)*CAO(K,1) + CSP(5,4)*CAO(K,4)
               GAO(K,ISTART1,IX) = CSP(1,2)*CAOX(K,2)
               GAO(K,ISTART2,IX) = CSP(2,5)*CAOX(K,5)
               GAO(K,ISTART3,IX) = CSP(3,1)*CAOX(K,1) + CSP(3,4)*CAOX(K,4)& 
               &                                + CSP(3,6)*CAOX(K,6)
               GAO(K,ISTART4,IX) = CSP(4,3)*CAOX(K,3)
               GAO(K,ISTART5,IX) = CSP(5,1)*CAOX(K,1) + CSP(5,4)*CAOX(K,4)
               GAO(K,ISTART1,IY) = CSP(1,2)*CAOY(K,2)
               GAO(K,ISTART2,IY) = CSP(2,5)*CAOY(K,5)
               GAO(K,ISTART3,IY) = CSP(3,1)*CAOY(K,1) + CSP(3,4)*CAOY(K,4)&
               &                                + CSP(3,6)*CAOY(K,6)
               GAO(K,ISTART4,IY) = CSP(4,3)*CAOY(K,3)
               GAO(K,ISTART5,IY) = CSP(5,1)*CAOY(K,1) + CSP(5,4)*CAOY(K,4)
               GAO(K,ISTART1,IZ) = CSP(1,2)*CAOZ(K,2)
               GAO(K,ISTART2,IZ) = CSP(2,5)*CAOZ(K,5)
               GAO(K,ISTART3,IZ) = CSP(3,1)*CAOZ(K,1) + CSP(3,4)*CAOZ(K,4)& 
               &                                + CSP(3,6)*CAOZ(K,6)
               GAO(K,ISTART4,IZ) = CSP(4,3)*CAOZ(K,3)
               GAO(K,ISTART5,IZ) = CSP(5,1)*CAOZ(K,1) + CSP(5,4)*CAOZ(K,4)
            END DO
         ELSE
            DO I = 1, KHKTA
               DO K = 1, NVCLEN
                  P0(K) = D0 
                  FX(K) = D0 
                  FY(K) = D0 
                  FZ(K) = D0 
               END DO
               DO J = 1, KCKTA
                  SPHFAC = CSP(I,J)
                  IF (ABS(SPHFAC).GT.D0) THEN
                     DO K = 1, NVCLEN
                        P0(K) = P0(K) + SPHFAC*CAO (K,J)
                        FX(K) = FX(K) + SPHFAC*CAOX(K,J)
                        FY(K) = FY(K) + SPHFAC*CAOY(K,J)
                        FZ(K) = FZ(K) + SPHFAC*CAOZ(K,J)
                     END DO
                  END IF
               END DO
               TEMPI = ISTART+I
               DO K = 1, NVCLEN
                  GAO(K,TEMPI,I0) = P0(K)
                  GAO(K,TEMPI,IX) = FX(K)
                  GAO(K,TEMPI,IY) = FY(K)
                  GAO(K,TEMPI,IZ) = FZ(K)
               END DO
            END DO
         END IF
      END IF
      call mem_dft_dealloc(P0)
      call mem_dft_dealloc(FX)
      call mem_dft_dealloc(FY)
      call mem_dft_dealloc(FZ)
   ELSE
      IF (.NOT.Spherical) THEN
         DO I = 1, KCKTA
            TEMPI = ISTART+I
            CALL LS_DZERO(GAO(1,TEMPI,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,TEMPI,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,TEMPI,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,TEMPI,IZ),NVCLEN)
         END DO
      ELSE
         IF (SHELLANGMOMA.EQ. 3) THEN
            CALL LS_DZERO(GAO(1,ISTART1,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART2,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART3,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART4,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART5,I0),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART1,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART2,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART3,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART4,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART5,IX),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART1,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART2,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART3,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART4,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART5,IY),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART1,IZ),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART2,IZ),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART3,IZ),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART4,IZ),NVCLEN)
            CALL LS_DZERO(GAO(1,ISTART5,IZ),NVCLEN)
         ELSE
            DO I = 1, KHKTA
               TEMPI = ISTART+I
               CALL LS_DZERO(GAO(1,TEMPI,I0),NVCLEN)
               CALL LS_DZERO(GAO(1,TEMPI,IX),NVCLEN)
               CALL LS_DZERO(GAO(1,TEMPI,IY),NVCLEN)
               CALL LS_DZERO(GAO(1,TEMPI,IZ),NVCLEN)
            END DO
         END IF
      END IF
   END IF
END IF

IF (Spherical) THEN
   call mem_dft_dealloc(CAO)
   call mem_dft_dealloc(CAOX)
   call mem_dft_dealloc(CAOY)
   call mem_dft_dealloc(CAOZ)
END IF
call mem_dft_dealloc(GA)
call mem_dft_dealloc(GU)

END SUBROUTINE II_BLGETGA1

!> \brief evaluates the GAO(gaussian atomic orbitals) + first and second geo derivatives
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGA2(LUPRI,NVCLEN,NactBAS,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
                     & CSP,GAOXX,GAOXY,GAOXZ,GAOYY,GAOYZ,GAOZZ,IADR,PA,PA2,&
                     &DFTHRI,CC,CCSTART,PRIEXP,LVALUETK,MVALUETK,NVALUETK,SPVAL,Spherical)
IMPLICIT NONE
INTEGER,intent(in)        :: NVCLEN,MXPRIM,KHKTA,JSTA,KCKTA,LUPRI,SHELLANGMOMA,NactBAS
INTEGER,intent(in)        :: SPVAL,CCSTART,IADR
REAL(REALK),intent(inout) :: GAOXX(NVCLEN,NactBAS), GAOXY(NVCLEN,NactBAS)
REAL(REALK),intent(inout) :: GAOXZ(NVCLEN,NactBAS), GAOYY(NVCLEN,NactBAS)
REAL(REALK),intent(inout) :: GAOYZ(NVCLEN,NactBAS), GAOZZ(NVCLEN,NactBAS)
REAL(REALK),intent(in)    :: PA(3,NVCLEN),PA2(NVCLEN),CSP(KHKTA,KCKTA),PRIEXP(MXPRIM),DFTHRI
INTEGER,intent(in)        :: LVALUETK(SPVAL),MVALUETK(SPVAL),NVALUETK(SPVAL)
LOGICAL,intent(in)        :: Spherical
TYPE(LSMATRIX),intent(in) :: CC
!
INTEGER,PARAMETER         :: I0=1, IX=2, IY=3, IZ=4
REAL(REALK),PARAMETER     :: D0 = 0.0E0_realk, D1 = 1.0E0_realk,D2 = 2.0E0_realk
REAL(REALK)               :: GA,GAZ,GAX,GAY
REAL(REALK)  :: CAOXX(KCKTA),CAOXY(KCKTA),CAOXZ(KCKTA)
REAL(REALK)  :: CAOYY(KCKTA),CAOYZ(KCKTA),CAOZZ(KCKTA)
REAL(REALK)  :: ALPHA,TALPH,TAPAX,TAPAY,TAPAZ,PXD,PYD,PZD,PXM,PYM,PZM
REAL(REALK)  :: PX0,PY0,PZ0,PXP,PYP,PZP,P000,GAXX,GAXY,GAXZ,GAYY,GAYZ,GAZZ
REAL(REALK)  :: SPHXX,SPHXY,SPHXZ,SPHYY,SPHYZ,SPHZZ,SPHFAC,PA_1,PA_2,PA_3,PA2_B
INTEGER      :: I,J,L,M,N,ICOMPA,IPRIMA,IV,TEMPI,ISTART,K
!CHANGE THIS ACCORDING TO NEW BLGETGA2 CHANGED BY ANDREAS
ISTART=IADR
DO IV = 1, NVCLEN
   PA_1=PA(1,IV)
   PA_2=PA(2,IV)
   PA_3=PA(3,IV)
   PA2_B=PA2(IV)
   IF (Spherical) THEN
      DO I=1,KCKTA
         CAOXX(I) = D0
         CAOXY(I) = D0
         CAOXZ(I) = D0
         CAOYY(I) = D0
         CAOYZ(I) = D0
         CAOZZ(I) = D0
      END DO
   END IF
   DO I = 1,CC%nrow
      IPRIMA = JSTA+CCSTART-1+I
      ALPHA = PRIEXP(IPRIMA)
      TALPH = -D2*ALPHA
      TAPAX = TALPH*PA_1
      TAPAY = TALPH*PA_2
      TAPAZ = TALPH*PA_3
      GA = CC%elms(I)*EXP(-ALPHA*PA2(IV))!FAC
      IF (ABS(GA).GT.DFTHRI) THEN
         DO ICOMPA = 1, KCKTA
            L = LVALUETK(ICOMPA)
            M = MVALUETK(ICOMPA)
            N = NVALUETK(ICOMPA)
!                            
            PXD = D0
            PYD = D0
            PZD = D0
            IF (L.GT. 1) PXD = (L*(L-1))*(PA_1**(L-2))
            IF (M.GT. 1) PYD = (M*(M-1))*(PA_2**(M-2))
            IF (N.GT. 1) PZD = (N*(N-1))*(PA_3**(N-2))
            PXM = D0
            PYM = D0
            PZM = D0
            IF (L.GT. 0) PXM = (L)*(PA_1**(L-1))
            IF (M.GT. 0) PYM = (M)*(PA_2**(M-1))
            IF (N.GT. 0) PZM = (N)*(PA_3**(N-1))
            PX0 = PA_1**L
            PY0 = PA_2**M
            PZ0 = PA_3**N
            PXP = TAPAX*PX0
            PYP = TAPAY*PY0
            PZP = TAPAZ*PZ0
            P000 = PX0*PY0*PZ0 
            IF (SHELLANGMOMA.EQ. 1) THEN
               ! s orbitals
               GAXX = TAPAX*TAPAX + TALPH
               GAYY = TAPAY*TAPAY + TALPH
               GAZZ = TAPAZ*TAPAZ + TALPH
               GAXY = TAPAX*TAPAY
               GAXZ = TAPAX*TAPAZ
               GAYZ = TAPAY*TAPAZ
            ELSE IF (SHELLANGMOMA.EQ. 2) THEN
               ! p orbitals
               GAXX = (TAPAX*TAPAX + TALPH*(2*L+1))*P000
               GAYY = (TAPAY*TAPAY + TALPH*(2*M+1))*P000
               GAZZ = (TAPAZ*TAPAZ + TALPH*(2*N+1))*P000
               GAXY = TAPAX*TAPAY*P000+(PXP*PYM+PXM*PYP)*PZ0
               GAXZ = TAPAX*TAPAZ*P000+(PXP*PZM+PXM*PZP)*PY0
               GAYZ = TAPAY*TAPAZ*P000+(PYP*PZM+PYM*PZP)*PX0
            ELSE 
               ! d and higher orbitals
               GAXX = (TAPAX*TAPAX + TALPH*(2*L+1))*P000 + PXD*PY0*PZ0
               GAYY = (TAPAY*TAPAY + TALPH*(2*M+1))*P000 + PX0*PYD*PZ0
               GAZZ = (TAPAZ*TAPAZ + TALPH*(2*N+1))*P000 + PX0*PY0*PZD
               GAXY = TAPAX*TAPAY*P000 + (PXP*PYM+PXM*PYP+PXM*PYM)*PZ0
               GAXZ = TAPAX*TAPAZ*P000 + (PXP*PZM+PXM*PZP+PXM*PZM)*PY0
               GAYZ = TAPAY*TAPAZ*P000 + (PYP*PZM+PYM*PZP+PYM*PZM)*PX0
            END IF
            IF (Spherical.AND.SHELLANGMOMA.GT. 2) THEN
               CAOXX(ICOMPA) = CAOXX(ICOMPA) + GAXX*GA 
               CAOXY(ICOMPA) = CAOXY(ICOMPA) + GAXY*GA 
               CAOXZ(ICOMPA) = CAOXZ(ICOMPA) + GAXZ*GA
               CAOYY(ICOMPA) = CAOYY(ICOMPA) + GAYY*GA
               CAOYZ(ICOMPA) = CAOYZ(ICOMPA) + GAYZ*GA
               CAOZZ(ICOMPA) = CAOZZ(ICOMPA) + GAZZ*GA
            ELSE
               TEMPI = ISTART+ICOMPA
               GAOXX(IV,TEMPI) = GAOXX(IV,TEMPI) + GAXX*GA 
               GAOXY(IV,TEMPI) = GAOXY(IV,TEMPI) + GAXY*GA 
               GAOXZ(IV,TEMPI) = GAOXZ(IV,TEMPI) + GAXZ*GA
               GAOYY(IV,TEMPI) = GAOYY(IV,TEMPI) + GAYY*GA
               GAOYZ(IV,TEMPI) = GAOYZ(IV,TEMPI) + GAYZ*GA
               GAOZZ(IV,TEMPI) = GAOZZ(IV,TEMPI) + GAZZ*GA
            ENDIF
         END DO
      END IF
   END DO
   IF (Spherical.AND.SHELLANGMOMA.GT. 2) THEN
      DO I = 1,KHKTA
         SPHXX = D0 
         SPHXY = D0 
         SPHXZ = D0 
         SPHYY = D0 
         SPHYZ = D0 
         SPHZZ = D0 
         DO J = 1, KCKTA
            SPHFAC = CSP(I,J)
!            IF (ABS(SPHFAC).GT.D0) THEN
               SPHXX = SPHXX + SPHFAC*CAOXX(J)
               SPHXY = SPHXY + SPHFAC*CAOXY(J)
               SPHXZ = SPHXZ + SPHFAC*CAOXZ(J)
               SPHYY = SPHYY + SPHFAC*CAOYY(J)
               SPHYZ = SPHYZ + SPHFAC*CAOYZ(J)
               SPHZZ = SPHZZ + SPHFAC*CAOZZ(J)
!            END IF
         END DO
         TEMPI = ISTART+I
         GAOXX(IV,TEMPI) = SPHXX
         GAOXY(IV,TEMPI) = SPHXY
         GAOXZ(IV,TEMPI) = SPHXZ
         GAOYY(IV,TEMPI) = SPHYY
         GAOYZ(IV,TEMPI) = SPHYZ
         GAOZZ(IV,TEMPI) = SPHZZ
      END DO
   END IF
END DO

END SUBROUTINE II_BLGETGA2

!> \brief evaluates the london derivative GAOs (gaussian atomic orbitals)
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_BLGETGB1(NVCLEN,NactBAS,IADR,KHKTA,GAO,PA,GABX,GABY,GABZ)
IMPLICIT NONE
INTEGER,intent(in) :: NVCLEN,NactBAS,IADR,KHKTA
REAL(REALK),intent(in) :: PA(3,NVCLEN),GAO(NVCLEN,NactBAS)
REAL(REALK),intent(inout) :: GABX(NVCLEN,NactBAS)
REAL(REALK),intent(inout) :: GABZ(NVCLEN,NactBAS),GABY(NVCLEN,NactBAS)
!
REAL(REALK),PARAMETER :: D05 = 0.5E0_realk
INTEGER :: I,J,K
REAL(REALK) :: GA

DO I = 1,KHKTA
 J = IADR+I
 DO K = 1, NVCLEN
    GA = D05*GAO(K,J)
    GABX(K,J) = GA*PA(1,K)
    GABY(K,J) = GA*PA(2,K)
    GABZ(K,J) = GA*PA(3,K)
 END DO
END DO

END SUBROUTINE II_BLGETGB1

!> \brief make blocks that contain the orbitals which contribute to each gridbatch
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE SHELL_TO_ORB(LUPRI,NSHELLBLOCKS,SHELBLOCK,ORBBLOCKS,NSTART,MAXNSHELL,NBAST)
!SHELL_TO_ORB transforms shell block indices to orbital block indices.
IMPLICIT NONE
INTEGER,intent(in)  :: NSHELLBLOCKS,MAXNSHELL,NBAST,LUPRI
INTEGER,intent(inout)  :: ORBBLOCKS(2,MAXNSHELL)
INTEGER,intent(in)  :: SHELBLOCK(2,MAXNSHELL),NSTART(MAXNSHELL)
!
INTEGER  :: ISHELL,I

DO I = 1,NSHELLBLOCKS
   ORBBLOCKS(1,I) = NSTART(SHELBLOCK(1,I))+1
ENDDO

DO I = 1,NSHELLBLOCKS
   IF(SHELBLOCK(2,I) .LT. MAXNSHELL)THEN
      ORBBLOCKS(2,I) = NSTART(SHELBLOCK(2,I)+1)
   ELSE
      ORBBLOCKS(2,I) = NBAST
   ENDIF
!   WRITE(LUPRI,*) '("shell ",2I4," translated to orbital ",2I4)'&
!    &,SHELBLOCK(1,I),SHELBLOCK(2,I),ORBBLOCKS(1,I),ORBBLOCKS(2,I)
ENDDO

END SUBROUTINE SHELL_TO_ORB

!> \brief make blocks that contain the active orbitals which contribute to each gridbatch
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE SHELL_TO_ACTORB(NSHELLBLOCKS,BLOCKS,MAXNSHELL)
!SHELL_TO_ORB transforms shell block indices to orbital block indices.
IMPLICIT NONE
INTEGER,intent(in)  :: NSHELLBLOCKS,MAXNSHELL
INTEGER,intent(inout) :: BLOCKS(2,MAXNSHELL)
!
INTEGER  :: IBASIS,IBL,IORB1,IORB2

IBASIS = 0
DO IBL = 1, NSHELLBLOCKS
   IBASIS = IBASIS + 1
   IORB1 = BLOCKS(1,IBL)
   BLOCKS(1,IBL) = IBASIS
   IORB2 = BLOCKS(2,IBL) !ORBITALINDEX
   IBASIS = IBASIS + (IORB2-IORB1)
   BLOCKS(2,IBL) = IBASIS
ENDDO

END SUBROUTINE SHELL_TO_ACTORB

!> \brief build spherical transformation matrices
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE Build_PRECALCULATED_SPHMAT(LUPRI,MAXANGMOM,SIZE,SPHMAT,SPINDEX) 
use math_fun
IMPLICIT NONE
INTEGER,intent(in)        :: MAXANGMOM,SIZE,LUPRI
INTEGER,intent(inout)     :: SPINDEX(MAXANGMOM+1)
REAL(REALK),intent(inout) :: SPHMAT(SIZE)
!
INTEGER          :: L,I
Real(realk), parameter :: DM1 = -1.0E0_realk, DO = 0.0E0_realk, D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER  :: M1,MADR,MABS,V0, NDER, IOFF,nAngmom
INTEGER  :: EE,FF,GG,BX,BY,BZ,II,JJ,KK,PX,PY,PZ
REAL(realk)  :: FACNRM,FAC1,FAC2, FAC3,FAC4, FACTOR
INTEGER  :: T,U,V,A,B,C,P,Q,R,X,Y,Z,TOT,IADR,AX0,AY0,AZ0,NSIZE
INTEGER  :: M,Ncol,Nrow,INDEX,NDIM,STARTINDEX

IF(MAXANGMOM .LE. 1)CALL LSQUIT('ERROR IN Build_PRECALCULATED_SPHMAT',lupri)
! CALL LS_DZERO(SPHMAT,SIZE) SHOULD BE DONE OUTSIDE

nANGMOM=MAXANGMOM+1 
NSIZE=0
STARTINDEX=1
SPINDEX(1)=1
SPINDEX(2)=1
DO I=3,nANGMOM
   SPINDEX(I)=STARTINDEX
   L = I-1 !angmom
   NRow = 2*L+1
   NCol = (L+1)*(L+2)/2
DO M1 = 0, 2*L 
   M = M1 - L
   IF (L.EQ. 1) THEN
      IF (M .EQ. -1) MADR =  0  
      IF (M .EQ.  0) MADR =  1 
      IF (M .EQ.  1) MADR = -1 
   ELSE
      MADR = M
   END IF
   MABS = ABS(M)
   V0 = 0
   IF (M .LT. 0) V0 = 1 
   FACNRM = D1
   IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(LUPRI,L+MABS)*&
        &FACULT(LUPRI,L-MABS))/(FACULT(LUPRI,L)*(D2**MABS))
   FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
   FACNRM = FACNRM/SQRT(FACUL2(LUPRI,2*L-1))
   DO T = 0, L - MABS, 2
   DO U = 0, T, 2
   DO V = V0, MABS, 2
      !        almost 6.4.48 in the book
      FAC3 = FACNRM*BINOM(LUPRI,L,T/2)*BINOM(LUPRI,L-T/2,MABS+T/2)&
           &                    *BINOM(LUPRI,T/2,U/2)*BINOM(LUPRI,MABS,V)
      DO A = 0, MIN(0,T+MABS-U-V) 
      DO B = 0, MIN(0,U+V)
      DO C = 0, MIN(0,L-T-MABS)
         !           6.4.47 in the book
         DO P = 0, - A, 2
         DO Q = 0, - B, 2
         DO R = 0, - C, 2
            FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
                 &   D2**(-A-B-C-P-Q-R-T)*FAC3
            X = T+MABS-U-V-2*A-P
            Y = U+V-2*B-Q
            Z = L-T-MABS-2*C-R
            TOT = X + Y + Z
            IADR = 1 + (2*L+1)*(NCRT(X,Y,Z)-1) + L + MADR+NSIZE
            SPHMAT(IADR) = SPHMAT(IADR) + FACTOR 
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
   ENDDO
   ENDDO
   ENDDO
ENDDO
   NSIZE= NSIZE+Nrow*Ncol 
   STARTINDEX=STARTINDEX+Nrow*Ncol
ENDDO
END SUBROUTINE BUILD_PRECALCULATED_SPHMAT

!
!     DFTDISP calculates the empirical dispersion correction to the energy and the gradient 
!     as formulated by Grimme in J. Comp. Chem. 2004, 25, 1463. and  J. Comp. Chem. 2006, 27, 1787.
!
!     EDISP = -S6 * sum_i^{N-1}sum_{j=i+1}^N(C_6(ij)/R(ij)^6 * f(R_ij))
!
!     with the damping function f(R_ij) = 1/(1+exp(-d*[R_ij/Rr-1]))
!
!     where N is the number of atoms
!           R_ij is the interatomic distance
!           Rr is the sum of the van-der-Waals radii
!           d is a fixed damping parameter (d=20.0)
!           C_6(ij) = sqrt(C_6(i)*C_6(j)
!           C_6 are fixed, atomic C_6 parameters
!           S_6 is a fixed, functional dependent global scaling factor (defined for BP86, BLYP, PBE, B3LYP, TPSS)
!
!     -------------------------------------------------------------
!     EDISP:  contains the energy correction at output
!     NDERIV: order of derivative wrt nuclei 
!             0: energy only
!             1: first derivative
!             higher: not defined yet
!     -------------------------------------------------------------
!
!     03.2010 Andreas Krapp 
!
SUBROUTINE II_DFTDISP(SETTING,GRAD,DIM1,DIM2,NDERIV,LUPRI,IPRINT)
use ls_util
IMPLICIT NONE
! external integer
INTEGER, INTENT(IN) :: NDERIV, LUPRI, DIM1, DIM2, IPRINT
! internal integer
INTEGER :: SHELL2ATOMA, SHELL2ATOMB, ISCOOA, ISCOOB, ISCOOR, IATOM, IOFF, J, NATOMS
! external real
REAL(REALK), INTENT(INOUT) :: GRAD(DIM1,DIM2)
! internal real
REAL(REALK) :: R0(54), C6(54)
REAL(REALK) :: E, EADD, S6
REAL(REALK) :: CHARGA, CORDAX, CORDAY, CORDAZ, C6A, RvdWA
REAL(REALK) :: CHARGB, CORDBX, CORDBY, CORDBZ, C6B, RvdWB
REAL(REALK) :: RX, RY, RZ, R2, R, RR, R6FAC, C6FAC
REAL(REALK) :: ALPHA, EXPOA, FDMP
REAL(REALK) :: DFAC, GRADX, GRADY, GRADZ
REAL(REALK), pointer :: GRDFT(:)
REAL(realk), PARAMETER ::  CONV1=1.88972612E0_realk   ! convert angstrom to bohr
REAL(realk), PARAMETER ::  CONV2=17.3452771E0_realk   ! convert Joule*nm^6/mol to Bohr^6*hartree
                                       !    1 Bohr = 52.9177 * 10^-3 nm
                                       !    1 Hartree = 2.6255*10^6 Joule/mol
!  external function
!EXTERNAL DISP_FUNCFAC
!  external types
TYPE(LSSETTING),   INTENT(INOUT) :: SETTING
!
! van der Waals radii for the elements H-Xe in Angstrom
! taken from JCC 2006, 27, 1787
!
     DATA &
!         H                He
     & R0/1.001E0_realk,1.012E0_realk, &
!         Li               Be            B              C
     &    0.825E0_realk,1.408E0_realk,1.485E0_realk,1.452E0_realk,&
!         N                O             F              Ne
     &    1.397E0_realk,1.342E0_realk,1.287E0_realk,1.243E0_realk,&
!         Na               Mg            Al             Si   
     &    1.144E0_realk,1.364E0_realk,1.639E0_realk,1.716E0_realk,&
! P     S     Cl    Ar
     &    1.705E0_realk,1.683E0_realk,1.639E0_realk,1.595E0_realk, &
!         K     Ca
     &    1.485E0_realk,1.474E0_realk, &
!         Sc-Zn,
     &    1.562E0_realk,1.562E0_realk,1.562E0_realk,1.562E0_realk,&
     &    1.562E0_realk,1.562E0_realk,1.562E0_realk,1.562E0_realk,1.562E0_realk,1.562E0_realk, &
!         Ga    Ge    As    Se    Br    Kr
     &    1.650E0_realk,1.727E0_realk,1.760E0_realk,1.771E0_realk,1.749E0_realk,1.727E0_realk, &
!         Rb    Sr    
     &    1.628E0_realk,1.606E0_realk, &
!         Y-Cd,
     &    1.639E0_realk,1.639E0_realk,1.639E0_realk,1.639E0_realk,1.639E0_realk,1.639E0_realk, &
     &    1.639E0_realk,1.639E0_realk,1.639E0_realk,1.639E0_realk, &
!         In    Sn    Sb    Te    I     Xe
     &    1.672E0_realk,1.804E0_realk,1.881E0_realk,1.892E0_realk,1.892E0_realk,1.881E0_realk/ 

!
!     C6 parameters for the elements H-Xe in Joule*nm^6/mol
!     taken from JCC 2006, 27, 1787
!
      DATA &
!         H     He
     & C6/0.14E0_realk ,0.08E0_realk , &
!         Li    Be    B     C     N     O     F     Ne
     &    1.61E0_realk ,1.61E0_realk ,3.13E0_realk ,1.75E0_realk ,&
     &    1.23E0_realk ,0.70E0_realk ,0.75E0_realk ,0.63E0_realk ,&
!         Na    Mg    Al    Si    P     S     Cl    Ar
     &    5.71E0_realk ,5.71E0_realk ,10.79E0_realk,9.23E0_realk ,&
     &    7.84E0_realk ,5.57E0_realk ,5.07E0_realk ,4.61E0_realk ,&
!         K     Ca
     &    10.80E0_realk,10.80E0_realk, &
!         Sc-Zn,
     &    10.80E0_realk,10.80E0_realk,10.80E0_realk,10.80E0_realk,&
     &    10.80E0_realk,10.80E0_realk,10.80E0_realk,10.80E0_realk,&
     &    10.80E0_realk,10.80E0_realk, &
!         Ga    Ge    As    Se    Br    Kr
     &    16.99E0_realk,17.10E0_realk,16.37E0_realk,12.64E0_realk,&
     &    12.47E0_realk,12.01E0_realk, &
!         Rb    Sr    
     &    24.67E0_realk,24.67E0_realk, &
!         Y-Cd,
     &    24.67E0_realk,24.67E0_realk,24.67E0_realk,24.67E0_realk,&
     &    24.67E0_realk,24.67E0_realk,24.67E0_realk,24.67E0_realk,&
     &    24.67E0_realk,24.67E0_realk, &
!         In    Sn    Sb    Te    I     Xe
     &    37.32E0_realk,38.71E0_realk,38.44E0_realk,31.74E0_realk,&
     &    31.50E0_realk,29.99E0_realk/


      ! only for DFT 
      IF (.NOT. SETTING%DO_DFT) RETURN

      ! we do not want dispersion correction -> return
      IF (.NOT. SETTING%SCHEME%DFT%DODISP) RETURN

      IF (NDERIV.EQ. 1) THEN
         SETTING%SCHEME%DFT%DISPDONE = .FALSE.
      ELSE
         ! if we have calculated the dispersion correction we can return
         IF (SETTING%SCHEME%DFT%DISPDONE) THEN 
            IF (IPRINT .GE. 1) THEN
               WRITE(LUPRI,'(1X,A39,F24.10)') 'dispersion correction to the KS-energy: ', SETTING%EDISP
               WRITE(LUPRI,*)''
            END IF
            RETURN
         ENDIF
         IF (SETTING%MOLECULE(1)%p%NATOMS.LE. 1) THEN 
            SETTING%SCHEME%DFT%DISPDONE = .FALSE.
            SETTING%EDISP = 0.0E0_realk
            RETURN
         ELSE
            SETTING%SCHEME%DFT%DISPDONE = .TRUE.
         END IF
      END IF

!     error check
      IF (NDERIV.GT. 1 .OR. NDERIV .LT. 0) THEN
         WRITE(LUPRI,'(4X,A62)') &
     &  ' WARNING: dispersion correction only for energies and gradients'
         WRITE(*,*) &
     &  'WARNING: dispersion correction only for energies and gradients'
      END IF

!     Print section
      WRITE(LUPRI,*)''
      WRITE(LUPRI,'(1X,A56)') &
     &   'Add empirical dispersion corr. to the XC-energy/gradient'
      WRITE(LUPRI,'(1X,A25)')'  following S. Grimme,   '
      WRITE(LUPRI,'(1X,A25)')'  JCC 2004, 25, 1463. and'
      WRITE(LUPRI,'(1X,A25)')'  JCC 2006, 27, 1787.    '

!     initialisations and memory for gradient
      NATOMS = SETTING%MOLECULE(1)%p%NATOMS
      E = 0.0E0_realk
      SETTING%EDISP = 0.0E0_realk
      IF (NDERIV.EQ. 1) THEN
         call mem_dft_alloc(GRDFT,DIM1*DIM2)
         CALL LS_DZERO(GRDFT,DIM1*DIM2)
      END IF

!     calculate correction

!
!     Run over nuclei A
!
      DO SHELL2ATOMA = 1, NATOMS-1
         CHARGA = SETTING%MOLECULE(1)%p%ATOM(SHELL2ATOMA)%CHARGE
         IF (ABS(CHARGA-NINT(CHARGA)).GT. 1E-10_realk) THEN
            CALL LSQUIT('Error in DFTDISP. Not implemented for&
                       & non-integer charge!',-1)
         ENDIF
         IF (NINT(CHARGA) .LE. 54 .AND. NINT(CHARGA) .GT. 0) THEN
            ISCOOA = (SHELL2ATOMA-1)*3
            CORDAX = SETTING%MOLECULE(1)%p%ATOM(SHELL2ATOMA)%CENTER(1)
            CORDAY = SETTING%MOLECULE(1)%p%ATOM(SHELL2ATOMA)%CENTER(2)
            CORDAZ = SETTING%MOLECULE(1)%p%ATOM(SHELL2ATOMA)%CENTER(3)
            C6A    = C6(NINT(CHARGA))*CONV2
            RvdWA  = R0(NINT(CHARGA))*CONV1
!
!           Run over nuclei B
!
            DO SHELL2ATOMB =  SHELL2ATOMA+1, NATOMS
               CHARGB = SETTING%MOLECULE(1)%p%ATOM(SHELL2ATOMB)%CHARGE
               IF (ABS(CHARGB-NINT(CHARGB)).GT. 1E-10_realk) THEN
                  CALL LSQUIT('Error in DFTDISP. Not implemented for&
                             & non-integer charge!',-1)
               ENDIF
               IF (NINT(CHARGB).LE. 54 .AND.NINT(CHARGB).GT. 0) THEN
                  ISCOOB = (SHELL2ATOMB-1)*3
                  CORDBX = SETTING%MOLECULE(1)%p%ATOM(SHELL2ATOMB)%CENTER(1)
                  CORDBY = SETTING%MOLECULE(1)%p%ATOM(SHELL2ATOMB)%CENTER(2)
                  CORDBZ = SETTING%MOLECULE(1)%p%ATOM(SHELL2ATOMB)%CENTER(3)
                  C6B    = C6(NINT(CHARGB))*CONV2
                  RvdWB  = R0(NINT(CHARGB))*CONV1

!                 C6 factor for atom pair A-B
                  C6FAC = SQRT(C6A*C6B)

!                 distance R between atoms A and B
!                 and R^6 
                  RX    = CORDAX-CORDBX
                  RY    = CORDAY-CORDBY
                  RZ    = CORDAZ-CORDBZ
                  R2    = RX**2 + RY**2 + RZ**2
                  R6FAC = R2*R2*R2
                  R     = SQRT(R2)

!                 sum of the van der Waals radii RR  
                  RR   = RvdWA + RvdWB

!                 damping function FDMP
                  ALPHA = -20.0*R/RR+20.0E0_realk
                  EXPOA = EXP(ALPHA)
                  FDMP  = 1.0E0_realk/(1.0E0_realk + EXPOA)

!                 dispersion correction contribution
                  EADD  = C6FAC/R6FAC * FDMP
                  E = E + EADD

                  IF (NDERIV.EQ. 1) THEN
!                    derivative wrt to nuclear positions
                     DFAC = (-6.0/R + 20.0/RR * EXPOA * FDMP)*EADD / R
                     GRADX = DFAC * RX
                     GRADY = DFAC * RY
                     GRADZ = DFAC * RZ
                     GRDFT(ISCOOA + 1) = GRDFT( ISCOOA + 1) + GRADX
                     GRDFT(ISCOOA + 2) = GRDFT( ISCOOA + 2) + GRADY
                     GRDFT(ISCOOA + 3) = GRDFT( ISCOOA + 3) + GRADZ
                     GRDFT(ISCOOB + 1) = GRDFT( ISCOOB + 1) - GRADX
                     GRDFT(ISCOOB + 2) = GRDFT( ISCOOB + 2) - GRADY
                     GRDFT(ISCOOB + 3) = GRDFT( ISCOOB + 3) - GRADZ
                  END IF

               ELSE
                  WRITE(LUPRI,'(4X,A42)') 'DISPERSION CORRECTION ONLY FOR ATOMS 1-54.'
                  WRITE(LUPRI,'(4X,A30,I6,A4,I6)') 'ACTUAL CHARGE FOR NUCLEUS ',SHELL2ATOMB,' IS ',NINT(CHARGB)
                  CALL LSQUIT('DISPERSION CORRECTION ONLY FOR ATOMS 1-54.',-1)
               END IF
            END DO
         ELSE 
             WRITE(LUPRI,'(4X,A42)') 'DISPERSION CORRECTION ONLY FOR ATOMS 1-54.'
             WRITE(LUPRI,'(4X,A30,I6,A4,I6)') 'ACTUAL CHARGE FOR NUCLEUS ',SHELL2ATOMA,' IS ',NINT(CHARGA)
             CALL LSQUIT('DISPERSION CORRECTION ONLY FOR ATOMS 1-54.',-1)
         END IF
      END DO

!     final dispersion correction
!     S6 factor is functional dependant
      CALL DISP_FUNCFAC(S6)
      SETTING%EDISP = -S6 * E

!     Add derivative
      IF (NDERIV.EQ. 1) THEN
         IF( (DIM1.NE. 3) .OR. (DIM2.NE.NATOMS)) THEN 
            CALL LSQUIT('ERROR IN DFTDISP WITH GRADIENT',-1)
         ENDIF
         DO IATOM = 1, NATOMS
            ISCOOR = (IATOM-1)*3 
            GRAD(1,IATOM) = GRAD(1,IATOM) - GRDFT(ISCOOR+1)*S6
            GRAD(2,IATOM) = GRAD(2,IATOM) - GRDFT(ISCOOR+2)*S6
            GRAD(3,IATOM) = GRAD(3,IATOM) - GRDFT(ISCOOR+3)*S6
         END DO
         call mem_dft_dealloc(GRDFT)
      END IF

!     Print section
      WRITE(LUPRI,'(1X,A13,F7.2)')'    S6 factor',S6
      WRITE(LUPRI,'(1X,A15)')     '    C6 factors:'
      DO IATOM = 1, NATOMS
         WRITE (LUPRI, '(2X,A6,F7.3)') &
     &          SETTING%MOLECULE(1)%p%ATOM(IATOM)%NAME, C6(NINT(SETTING%MOLECULE(1)%p%ATOM(IATOM)%CHARGE))
      END DO

      IF (IPRINT .GE. 1) THEN
         WRITE(LUPRI,'(1X,A39,F24.10)') 'dispersion correction to the KS-energy: ', SETTING%EDISP
         IF (NDERIV.EQ. 1) THEN
            CALL LSHEADER(LUPRI,'XC-gradient including empir. disp. corr.')
            DO IATOM = 1, NATOMS
               WRITE (LUPRI, '(1X,A6,F17.10,2F24.10)') &
     &            SETTING%MOLECULE(1)%p%ATOM(IATOM)%NAME, (GRAD(J,IATOM),J=1,3)
            END DO
         END IF
         WRITE(LUPRI,*)''
      END IF

   RETURN
END SUBROUTINE II_DFTDISP

SUBROUTINE II_DFTsetFunc(Func,hfweight)
#ifdef VAR_MPI
use infpar_module
use lsmpi_mod
use lsmpi_type
#endif
implicit none
Character(len=80),intent(IN) :: Func
Real(realk),intent(INOUT)    :: hfweight
integer                      :: ierror
CALL DFTsetFunc(Func,hfweight,ierror)
IF(ierror.NE.0)CALL LSQUIT('Unknown Functional',-1)
#ifdef VAR_MPI
!for MPI ne also need to set the functional on the slaves
IF (infpar%mynum.EQ.infpar%master) THEN
  call ls_mpibcast(DFTSETFU,infpar%master,MPI_COMM_LSDALTON)
  call lsmpi_setmasterToSlaveFunc(Func,hfweight)
ELSE
  call lsquit('Error in II_DFTsetFunc. Can only be called from the master',-1)
ENDIF
#endif
END SUBROUTINE II_DFTsetFunc

SUBROUTINE II_DFTaddFunc(Func,GGAfactor)
#ifdef VAR_MPI
use infpar_module
use lsmpi_mod
use lsmpi_type
#endif
implicit none
Character(len=80),intent(IN) :: Func
Real(realk),intent(IN)       :: GGAfactor
!
CALL DFTaddFunc(Func,GGAfactor)
#ifdef VAR_MPI
!for MPI ne also need to set the functional on the slaves
IF (infpar%mynum.EQ.infpar%master) THEN
  call ls_mpibcast(DFTADDFU,infpar%master,MPI_COMM_LSDALTON)
  call lsmpi_setmasterToSlaveFunc(Func,GGAfactor)
ELSE
  call lsquit('Error in II_DFTsetFunc. Can only be called from the master',-1)
ENDIF
#endif
END SUBROUTINE II_DFTaddFunc
END MODULE IIDFTINT
