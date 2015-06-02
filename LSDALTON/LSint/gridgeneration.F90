!> @file
!> Module contains grid generation routines
!> \brief 
!> \author T. Kjaergaard
!> \date 2009 
MODULE gridgenerationmodule
use grid_memory_handling
use gridgenerationboxmodule
!WARNING you must not add memory_handling, all memory goes through 
!grid_memory_handling  module so as to determine the memory used in this module.
use precision
use dft_ld_module
use files
use Fundamental, only: bohr_to_angstrom
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_mod
#endif

SAVE
! some constants for the SSF scheme
real(realk),parameter :: SSF_CUTOFF=0.64E0_realk
real(realk),parameter :: cutoff_inv = 1.0E0_realk/SSF_CUTOFF;
real(realk),parameter :: cutoffinv2 = cutoff_inv*cutoff_inv;
real(realk),parameter :: fac0 = 0.5E0_realk*0.0625E0_realk*cutoff_inv;
real(realk),parameter :: fac1 = 35E0_realk*fac0;
real(realk),parameter :: fac2 = -35E0_realk*cutoffinv2*fac0;
real(realk),parameter :: fac3 = 21E0_realk*cutoffinv2*cutoffinv2*fac0;
real(realk),parameter :: fac4 = -5E0_realk*cutoffinv2*cutoffinv2*cutoffinv2*fac0;
real(realk),parameter :: weight_thr = 1e-15_realk
real(realk),parameter :: compress_thr = 1e-20_realk
INTEGER,parameter   :: NRADPT = 2000
!   Original version from Trond Saue:
!   The below data gives atomic radii in Angstroms and stems from table I of 
!   J.C.Slater: "Atomic Radii in Crystals"
!   J.Chem.Phys. 41(1964) 3199-3204
!   Values for elements (He,Ne,Ar,Kr,Xe,At,Rn,Pb,Fr,Cm,Bk,Cf,Es,Fm,Md,No,Lw) 
!   has been guessed/interpolated
!   bragg_radii(0) is for dummy atom
!   bragg_radii(1) is for Hydrogen
!   bragg_radii(2) is for Helium 
!   ...
real(realk),parameter :: bragg_radii(0:103)=(/ 0.75E0_realk, 0.35E0_realk, 0.35E0_realk, 1.45E0_realk,&
     & 1.05E0_realk,  0.85E0_realk,  0.70E0_realk,  0.65E0_realk,  0.60E0_realk,  0.50E0_realk,  0.45E0_realk,  1.80E0_realk, &
     & 1.50E0_realk,  1.25E0_realk,  1.10E0_realk,  1.00E0_realk,  1.00E0_realk,  1.00E0_realk,  1.00E0_realk,  2.20E0_realk, &
     & 1.80E0_realk,  1.60E0_realk,  1.40E0_realk,  1.35E0_realk,  1.40E0_realk,  1.40E0_realk,  1.40E0_realk,  1.35E0_realk, & 
     & 1.35E0_realk,  1.35E0_realk,  1.35E0_realk,  1.30E0_realk,  1.25E0_realk,  1.15E0_realk,  1.15E0_realk,  1.15E0_realk, &
     & 1.10E0_realk,  2.35E0_realk,  2.00E0_realk,  1.80E0_realk,  1.55E0_realk,  1.45E0_realk,  1.45E0_realk,  1.35E0_realk, &
     & 1.30E0_realk,  1.35E0_realk,  1.40E0_realk,  1.60E0_realk,  1.55E0_realk,  1.55E0_realk,  1.45E0_realk,  1.45E0_realk, &
     & 1.40E0_realk,  1.40E0_realk,  1.40E0_realk,  2.60E0_realk,  2.15E0_realk,  1.95E0_realk,  1.85E0_realk,  1.85E0_realk, &
     & 1.85E0_realk,  1.85E0_realk,  1.85E0_realk,  1.85E0_realk,  1.80E0_realk,  1.75E0_realk,  1.75E0_realk,  1.75E0_realk, &
     & 1.75E0_realk,  1.75E0_realk,  1.75E0_realk,  1.75E0_realk,  1.55E0_realk,  1.45E0_realk,  1.35E0_realk,  1.30E0_realk, &
     & 1.30E0_realk,  1.35E0_realk,  1.35E0_realk,  1.35E0_realk,  1.50E0_realk,  1.90E0_realk,  1.75E0_realk,  1.60E0_realk, &
     & 1.90E0_realk,  1.50E0_realk,  1.50E0_realk,  2.15E0_realk,  2.15E0_realk,  1.95E0_realk,  1.80E0_realk,  1.80E0_realk, &
     & 1.37E0_realk,  1.75E0_realk,  1.75E0_realk,  1.75E0_realk,  1.75E0_realk,  1.75E0_realk,  1.75E0_realk,  1.75E0_realk, &
     & 1.75E0_realk,  1.75E0_realk,  1.75E0_realk,  1.75E0_realk /)
CONTAINS
!> \brief wrapper exchange-correlation integral routine that build basinf.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE GenerateGrid(NBAST,radint,angint,ihardness,iprune,natoms,& 
     & X,Y,Z,Charge,GRIDDONE,SHELL2ATOM,SHELLANGMOM,SHELLNPRIM,MAXANGMOM,&
     & MAXNSHELL,MXPRIM,PRIEXP,PRIEXPSTART,RSHEL,ITERATIONS,TURBO,MaxNbuflen,&
     & RADIALGRID,ZdependenMaxAng,PARTITIONING,nstart,MaxNactBast,LUPRI,&
     & IPRINT,USE_MPI,numnodes,node,GridId,Gridnumnodes)
IMPLICIT NONE
!> How many times do we write (and should we read)  from disk 
INTEGER  :: ITERATIONS
!> the logical unit number for the output file
INTEGER :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER  :: IPRINT
!> number of basis functions
INTEGER  :: NBAST
!> maximum number of shells
INTEGER  :: MAXNSHELL
!> number of atoms
INTEGER  :: natoms
!> the maximum angular momentum
INTEGER  :: MAXANGMOM
!> if the grid is done GRIDDONE=1 else GRIDDONE=0
INTEGER :: GRIDDONE
!> MXPRIM is the total number of (unique) primitive orbitals
INTEGER  :: MXPRIM
!> X coordinate for each atom
REAL(REALK),intent(in):: X(natoms)
!> Y coordinate for each atom
REAL(REALK),intent(in):: Y(natoms)
!> Z coordinate for each atom
REAL(REALK),intent(in):: Z(natoms)
!> charge for each atom used in grid-generation
INTEGER,intent(in)  :: CHARGE(natoms)
!> which atomic center the shell is attached to
INTEGER,intent(in)    :: SHELL2ATOM(MAXNSHELL)
!> the angular momentum for each shell
INTEGER,intent(in)    :: SHELLANGMOM(MAXNSHELL)
!> the number of primitives for each shell
INTEGER,intent(in)    :: SHELLNPRIM(MAXNSHELL)
!> the unique primitve exponents
REAL(REALK),intent(in):: PRIEXP(MXPRIM)
!> the index to start in PRIEXP(MXPRIM) for a given shell index 
INTEGER,intent(in) :: PRIEXPSTART(MAXNSHELL)
!> the radius of each shell 
REAL(REALK),intent(in):: RSHEL(MAXNSHELL)
!> some grid generation specification
REAL(REALK) :: radint
!> some grid generation specification
INTEGER :: angint
!> some grid generation specification
INTEGER :: IHARDNESS
!> some grid generation specification
INTEGER  :: TURBO
!> grid point pruning (true if iprune=1) 
INTEGER :: iprune
!> max length of coordinate buffer -> COOR(3,maxNBUFLEN)
INTEGER,intent(in) :: maxNBUFLEN
!> which radial grid should be used?
INTEGER,intent(in) :: RADIALGRID !(1=GC2,2=LMG,3=TURBO) 
!> should the maximum angular momentum be Z dependent
LOGICAL,intent(in) :: ZdependenMaxAng
!> which partitioning should be used
INTEGER,intent(in) :: PARTITIONING !(1=SSF,2=BECKE,3=BECKEORIG,4=BLOCK,...)
!> maximum number of active orbitals for each box
INTEGER,intent(inout) :: MaxNactBast
!> for a given shell index it gives the corresponding starting orbital
INTEGER,intent(in)    :: NSTART(MAXNSHELL)
!> should we use MPI
LOGICAL,intent(in)    :: USE_MPI
integer(kind=ls_mpik),intent(in)    :: numnodes,node
integer,intent(inout) :: Gridnumnodes
!> if the grid id
INTEGER :: GRIDID
!
!REAL(REALK),pointer :: AA(:,:,:)
INTEGER,pointer     :: nRadialPoints(:),GRIDANG(:,:)
REAL(REALK),pointer :: RADIALPOINTS(:,:),RADIALWEIGHT(:,:)
integer :: totalpoints
integer :: angularpoints ,I,leb_gen_from_point,iang1
!max number of radial points
integer,parameter :: MIN_RAD_PT = 20

IF(GRIDDONE.EQ. 1)RETURN  !Grid already calculated and the grid points are in DALTON.QUAD
!initialise memory counters
call init_gridmemvar()
iang1 = leb_get_from_order(angint)
angularpoints = leb_gen_points(iang1)

IF(iprint.gt. 1)THEN
   WRITE(lupri,'(A,I4,A,I6)')'Angular integration order ',angint,' Angular points ',angularpoints
ENDIF

!mem alloc RADIALPOINTS,RADIALWEIGHT,nRadialPoints = RadialGridPoints,RadialWeights,NumberOfGridPoints 
call mem_grid_alloc(RADIALPOINTS,NRADPT,NATOMS)
call mem_grid_alloc(RADIALWEIGHT,NRADPT,NATOMS)
call mem_grid_alloc(nRadialPoints,NATOMS)   
!Build RADIALPOINTS,RADIALWEIGHT,nRadialPoints = RadialGridPoints,RadialWeights,NumberOfGridPoints 
CALL SET_RADIAL(RADIALPOINTS,RADIALWEIGHT,nRadialPoints,NRADPT,SHELL2ATOM,SHELLANGMOM,SHELLNPRIM,MAXANGMOM,NATOMS,MAXNSHELL,&
        & MXPRIM,PRIEXP,PRIEXPSTART,IPRINT,LUPRI,RADINT,Charge,RADIALGRID,MIN_RAD_PT)
! Computes the angular points for a set of radial points 
! obtained from radial integration scheme.
! 
call mem_grid_alloc(GRIDANG,NRADPT,NATOMS)
CALL SET_ANGULAR(gridang,radint,angint,CHARGE,natoms,NRADPT,nRadialPoints,RADIALPOINTS,RADIALWEIGHT,&
     &           iprune,angularpoints,ZdependenMaxAng)
! Build Atomic and Molecular grids + partitioning
CALL ComputeCoords(totalpoints,nRadialPoints,natoms,RADIALPOINTS,RADIALWEIGHT,NRADPT,GRIDANG,&
     & X,Y,Z,ihardness,MAXNSHELL,RSHEL,SHELL2ATOM,NBAST,maxNBUFLEN,ITERATIONS,&
     & PARTITIONING,Charge,nstart,MaxNactBast,iprint,lupri,USE_MPI,numnodes,node,GridId)

call mem_grid_dealloc(RADIALPOINTS)
call mem_grid_dealloc(RADIALWEIGHT)
call mem_grid_dealloc(nRadialPoints)   
call mem_grid_dealloc(GRIDANG)
Gridnumnodes = numnodes !number of nodes used to construct the grid. 
GRIDDONE=1
#ifdef VAR_MPI
IF (infpar%mynum.EQ.infpar%master) THEN
#endif
call stats_grid_mem(lupri)
#ifdef VAR_MPI
ENDIF
#endif

END SUBROUTINE GENERATEGRID

SUBROUTINE Computecoords(totalpoints,nRadialPoints,Natoms,RADIALPOINTS,RADIALWEIGHT,NRADPT,GRIDANG,&
     & atomcenterX,atomcenterY,atomcenterZ,ihardness,MAXNSHELL,RSHEL,SHELL2ATOM,NBAST,&
     & maxNBUFLEN,ITERATIONS,PARTITIONING,Charge,nstart,FinalMaxNactBast,iprint,lupri,USE_MPI,numnodes,node,GridId)
implicit none
integer,intent(in) :: Natoms,NRADPT,ihardness,iprint,lupri,nbast,maxnbuflen,GridId
!> X coordinate for each atom
REAL(REALK),intent(in):: atomcenterX(natoms)
!> Y coordinate for each atom
REAL(REALK),intent(in):: atomcenterY(natoms)
!> Z coordinate for each atom
REAL(REALK),intent(in):: atomcenterZ(natoms)
integer,intent(inout) :: ITERATIONS !number of times we write to disk 
integer,intent(inout) :: totalpoints !total number of gridpoints
integer,intent(in) :: nRadialPoints(nAtoms)
REAL(REALK),intent(in) :: RADIALPOINTS(NRADPT,NATOMS),RADIALWEIGHT(NRADPT,NATOMS)
integer,intent(in) :: GRIDANG(NRADPT,NATOMS)
!> maximum number of shells
INTEGER  :: MAXNSHELL
!> the radius of each shell (used in grid-generation)
REAL(REALK),intent(in):: RSHEL(MAXNSHELL)
!> which atomic center the shell is attached to
INTEGER,intent(in)    :: SHELL2ATOM(MAXNSHELL)
!> which partitioning should be used
INTEGER,intent(in) :: PARTITIONING !(1=SSF,2=BECKE,3=BECKEORIG,4=BLOCK,...)
!> charge for each atom used in grid-generation
INTEGER,intent(in)  :: CHARGE(natoms)
!> maximum number of active orbitals for each box
INTEGER,intent(inout) :: FinalMaxNactBast
!> for a given shell index it gives the corresponding starting orbital
INTEGER,intent(in)    :: NSTART(MAXNSHELL)
!> should we use MPI
LOGICAL,intent(in)    :: USE_MPI
INTEGER(kind=ls_mpik),intent(in)    :: numnodes,node
!
integer :: j,iatom,dummy,maxNR,N,N2,N3,I,idx,ipoint,iang
integer :: IATOM1,IATOM2,I2,h,igrid,itmp,totalAtomicpointscompressed
integer :: totalAtomicpoints,ix,iy,iz,NactBast,MaxNactBast
real(realk),parameter :: pim4=12.566370614359172E0_realk,D2=2E0_realk !4*pi
real(realk),pointer :: locR(:),WG(:),COOR(:,:)
real(realk),pointer :: WG2(:),COOR2(:,:)
real(realk),pointer :: WG3(:),COOR3(:,:)
real(realk) :: Xatomcoor2,Yatomcoor2,Zatomcoor2
real(realk) :: Xatomcoor,Yatomcoor,Zatomcoor,radial,weight,factor,DX,DY,DZ,MUA
real(realk),pointer :: inverseDistance12(:,:) ,BFAC(:,:)!,TMP(:,:),MU(:),MU2(:),TMP2(:,:)!,VEC(:)
real(realk),pointer :: SSF_ATOMIC_CUTOFF(:)
real(realk),parameter :: D1=1E0_realk,D0=0E0_realk,D05=0.5E0_realk,D3=3E0_realk
real(realk) :: W,RW,fac,center(3),dist
integer :: ikey,npoints2,NSHELLBLOCKS,LUGRID,IT,NGRIDPOINTINBATCH,nthreads,tid
integer :: maxN,MaxAtomicGridpoints,maxGridpoints,AtomicGridpoints,nprocessors
integer :: GlobalmaxAtomicGridpoints,GlobalmaxGridpoints,GlobalmaxN,mynum
integer,pointer :: SHELLBLOCKS(:,:),ATOMIDX(:)
integer :: nx,ny,nz,nkey,mykeysNumber,IMykey,privatetotalpoints
!integer(kind=long),pointer :: key(:),uniquekey(:)
type(bunchpoints),pointer :: keypointer(:)
logical :: unique,postprocess
character(len=22) :: filename
integer,pointer :: npoints(:,:,:),mykeys(:)
logical,pointer :: skip(:),skip2(:)
#ifdef VAR_OMP
integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
type(gridboxtype),pointer :: GridBox
#ifndef VAR_MPI
integer :: infpar
#endif
!logical,pointer :: SHELLLIST(:)

FinalMaxNactBast = 0
ITERATIONS=0
call mem_grid_alloc(inverseDistance12,NATOMS,NATOMS)
!This subroutine actually the inversedistance 1/r12
CALL determine_distance12(inverseDistance12,NATOMS,atomcenterX,atomcenterY,atomcenterZ)

IF(PARTITIONING.EQ. 2.OR.PARTITIONING.EQ. 4)THEN !needed for becke preprocess
   call mem_grid_alloc(Bfac,NATOMS,NATOMS)
   CALL Determine_BraggFac(NATOMS,BFAC,CHARGE)
ENDIF
IF(PARTITIONING.EQ. 5.OR.PARTITIONING.EQ. 1)THEN !needed for block-SSF postprocess
   call mem_grid_alloc(SSF_ATOMIC_CUTOFF,NATOMS)
   CALL Determine_SSF_ATOMIC_CUTOFF(NATOMS,SSF_ATOMIC_CUTOFF,&
        & atomcenterX,atomcenterY,atomcenterZ)
ENDIF
IF(PARTITIONING.EQ. 4.OR.PARTITIONING.EQ. 5)THEN
   postprocess=.TRUE.
ELSE
   postprocess=.FALSE.
ENDIF
IF(NATOMS.EQ. 1)postprocess=.FALSE.
! This can be improved by using atomic types
! we find the maximum gridpoints for a single atom 
! and the maximum number of gridpoints in total
nthreads=1
tid=0
call DetermineMaxGridPoints(GlobalmaxAtomicGridpoints,GlobalmaxGridpoints,&
     & GlobalmaxN,NATOMS,nRadialPoints,GRIDANG,NRADPT,tid,nthreads)
call mem_grid_alloc(WG2,GlobalmaxGridpoints)
call mem_grid_alloc(COOR2,3,GlobalmaxGridpoints)
totalpoints=0
IF(postprocess)THEN
   call mem_grid_alloc(ATOMIDX,GlobalmaxGridpoints)   
ENDIF
IF(PARTITIONING.EQ. 5)call mem_grid_alloc(skip2,GlobalmaxGridpoints)   

call mem_grid_TurnONThread_Memory()
!$OMP PARALLEL DEFAULT(none) PRIVATE(IATOM,Xatomcoor,Yatomcoor,Zatomcoor,idx,&
!$OMP IPOINT,iang,radial,weight,factor,N,locR,N2,N3,COOR,WG,totalAtomicpoints,&
!$OMP totalAtomicpointscompressed,igrid,ITMP,tid,nthreads,maxAtomicGridpoints,&
!$OMP maxGridpoints,I,maxN,dummy,skip,IMyKey,Ikey,npoints2,ix,iy,iz,center,&
!$OMP NSHELLBLOCKS,SHELLBLOCKS,w,rw,fac,NactBast,privatetotalpoints,WG3,&
!$OMP COOR3,MaxNactBast) SHARED(COOR2,WG2,USE_MPI,&
!$OMP atomcenterX,atomcenterY,atomcenterZ,totalpoints,inverseDistance12,ihardness,&
!$OMP bfac,GRIDANG,NRADPT,LUPRI,RADIALPOINTS,RADIALWEIGHT,PARTITIONING,IPRINT,&
!$OMP ATOMIDX,postprocess,skip2,SSF_ATOMIC_CUTOFF,nx,ny,nz,npoints,GridBox,&
!$OMP nkey,keypointer,LUGRID,ITERATIONS,filename,nprocessors,numnodes,node,mynum,mykeys,&
!$OMP mykeysnumber,MAXNSHELL,RSHEL,SHELL2ATOM,NBAST,GlobalmaxGridpoints,&
!$OMP NATOMS,nRadialPoints,maxNBUFLEN,nstart,infpar,GridId,FinalMaxNactBast)
call init_grid_threadmemvar()
#ifdef VAR_OMP
nthreads=OMP_GET_NUM_THREADS()
tid = omp_get_thread_num()
#ifdef VAR_MPI
IF (infpar%mynum.EQ.infpar%master) THEN
#endif
!$OMP MASTER
WRITE(lupri,'(4X,A,I3,A)')'This is an OpenMP Gridgeneration calculation using',&
     &omp_get_num_threads(),' threads.'
!$OMP END MASTER
#ifdef VAR_MPI
ENDIF
#endif
#else
nthreads=1
tid=0
#endif
! This can be improved by using atomic types
! we find the maximum gridpoints for a single atom 
! and the maximum number of gridpoints in total
! in case of OpenMP (OMP) these numbers is the values for each thread 
call DetermineMaxGridPoints(maxAtomicGridpoints,maxGridpoints,maxN,&
     & NATOMS,nRadialPoints,GRIDANG,NRADPT,tid,nthreads)
call mem_grid_alloc(locR,4*maxN)
call mem_grid_alloc(COOR,3,MaxAtomicGridpoints)
call mem_grid_alloc(WG,MaxAtomicGridpoints)
IF(PARTITIONING.EQ. 5)call mem_grid_alloc(skip,MaxAtomicGridpoints)   
!fixme add gridgeneration memory consumption!!

DO IATOM=1+tid,NATOMS,nthreads
   Xatomcoor = atomcenterX(IATOM)
   Yatomcoor = atomcenterY(IATOM)
   Zatomcoor = atomcenterZ(IATOM)
   idx = 0
   DO IPOINT = 1,nRadialPoints(IATOM)
      iang = GRIDANG(IPOINT,IATOM)
      radial = RADIALPOINTS(IPOINT,IATOM)
      weight = RADIALWEIGHT(IPOINT,IATOM)
      FACTOR = pim4*weight
      IF(iang.EQ. 1)THEN
         N=14
         call ld0014(locR(1:14),locR(15:28),locR(29:42),locR(43:56),dummy)
      ELSEIF(iang.EQ. 2)THEN
         N=38
         call ld0038(locR(1:38),locR(39:76),locR(77:114),locR(115:152),dummy)
      ELSEIF(iang.EQ. 3)THEN
         N=50
         call ld0050(locR(1:50),locR(51:100),locR(101:150),locR(151:200),dummy)
      ELSEIF(iang.EQ. 4)THEN
         N=86
         call ld0086(locR(1:86),locR(87:172),locR(173:258),locR(259:344),dummy)
      ELSEIF(iang.EQ. 5)THEN
         N=110
         call ld0110(locR(1:110),locR(111:220),locR(221:330),locR(331:440),dummy)
      ELSEIF(iang.EQ. 6)THEN
         N=146
         call ld0146(locR(1:146),locR(147:292),locR(293:438),locR(439:584),dummy)
      ELSEIF(iang.EQ. 7)THEN
         N=170
         call ld0170(locR(1:170),locR(171:340),locR(341:510),locR(511:680),dummy)
      ELSEIF(iang.EQ. 8)THEN
         N=194
         call ld0194(locR(1:194),locR(195:388),locR(389:582),locR(583:776),dummy)
      ELSEIF(iang.EQ. 9)THEN
         N=230
         call ld0230(locR(1:230),locR(231:460),locR(461:690),locR(691:920),dummy)
      ELSEIF(iang.EQ. 10)THEN
         N=266
         call ld0266(locR(1:266),locR(267:532),locR(533:798),locR(799:1064),dummy)
      ELSEIF(iang.EQ. 11)THEN
         N=302
         call ld0302(locR(1:302),locR(303:604),locR(605:906),locR(907:1208),dummy)
      ELSEIF(iang.EQ. 12)THEN
         N=350
         call ld0350(locR(1:350),locR(351:700),locR(701:1050),locR(1051:1400),dummy)
      ELSEIF(iang.EQ. 13)THEN
         N=434
         call ld0434(locR(1:434),locR(435:868),locR(869:1302),locR(1303:1736),dummy)
      ELSEIF(iang.EQ. 14)THEN
         N=590
         call ld0590(locR(1:590),locR(591:1180),locR(1181:1770),locR(1771:2360),dummy)
      ELSEIF(iang.EQ. 15)THEN
         N=770
         call ld0770(locR(1:770),locR(771:1540),locR(1541:2310),locR(2311:3080),dummy)
      ELSEIF(iang.EQ. 16)THEN
         N=974
         call ld0974(locR(1:974),locR(975:1948),locR(1949:2922),locR(2923:3896),dummy)
      ELSEIF(iang.EQ. 17)THEN
         N=1202
         call ld1202(locR(1:1202),locR(1203:2404),locR(2405:3606),locR(3607:4808),dummy)
      ELSEIF(iang.EQ. 18)THEN
         N=1454
         call ld1454(locR(1:1454),locR(1455:2908),locR(2909:4362),locR(4363:5816),dummy)
      ELSE
         CALL LSQUIT('wrong case in ld, maximum ANGINT exceeded',-1)
      ENDIF
      N2=2*N
      DO I=1,N
         COOR(1,I+idx) = locR(I)*radial    + Xatomcoor
         COOR(2,I+idx) = locR(N+I)*radial  + Yatomcoor
         COOR(3,I+idx) = locR(N2+I)*radial + Zatomcoor
      ENDDO
      N3=3*N
      DO I=1,N
         WG(I+idx) = locR(N3+I)*FACTOR
      ENDDO
      IF(NATOMS.GT. 1)THEN
         IF(PARTITIONING.EQ. 1)THEN !(.SSF)
            IF(radial .GE. SSF_atomic_cutoff(IATOM))THEN 
               ! Multiply current sphere/shell by space partitioning weights as
               ! described in SSF article. (Chem. Phys. Lett. 1996, 257, 213)
               CALL SSFPartitioning(N,NATOMS,atomcenterX,atomcenterY,atomcenterZ,&
                    & COOR,idx,inverseDistance12,ihardness,MaxAtomicGridpoints,WG,IATOM)
            ENDIF
         ELSEIF(PARTITIONING.EQ. 2)THEN !becke preprocess  (.BECKE)
            ! Becke scheme with atomic size correction
            CALL BeckePartitioning(N,NATOMS,atomcenterX,atomcenterY,atomcenterZ,&
                 & COOR,idx,inverseDistance12,bfac,ihardness,MaxAtomicGridpoints,WG,IATOM)
         ELSEIF(PARTITIONING.EQ. 3)THEN 
            ! Becke scheme without atomic size correction
            CALL BeckeOriginalPartitioning(N,NATOMS,atomcenterX,atomcenterY,&
                 & atomcenterZ,COOR,idx,inverseDistance12,ihardness,&
                 & MaxAtomicGridpoints,WG,IATOM)         
         ELSEIF(PARTITIONING.EQ. 5)THEN !block ssf 
            IF(radial .GE. SSF_atomic_cutoff(IATOM))THEN
               DO I=1,N
                  skip(I+idx)=.FALSE.
               ENDDO
            ELSE
               DO I=1,N
                  skip(I+idx)=.TRUE.
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      idx = idx+N
   ENDDO
   totalAtomicpoints = idx
   !compress points
   totalAtomicpointscompressed = totalAtomicpoints
   ! we find the first element which is lower than threshold
   igrid = 0
   compressloop: DO I=1,totalAtomicpoints
      IF(ABS(WG(I)).LT.compress_thr)THEN
         igrid = I
         EXIT compressloop
      ENDIF
   ENDDO compressloop
   if(igrid.gt. 0)THEN
      ! until igrid all elements are larger
      ITMP=igrid-1
      IF(postprocess.AND.PARTITIONING.EQ. 5)THEN
         DO I=igrid,totalAtomicpoints
            IF(ABS(WG(I)).GT.compress_thr)THEN
               ITMP=ITMP + 1
               COOR(1,iTMP) = COOR(1,I)
               COOR(2,iTMP) = COOR(2,I)
               COOR(3,iTMP) = COOR(3,I)
               WG(iTMP) = WG(I)
               skip(iTMP) = skip(I)
            ENDIF
         ENDDO
      ELSE
         DO I=igrid,totalAtomicpoints
            IF(ABS(WG(I)).GT.compress_thr)THEN
               ITMP=ITMP + 1
               COOR(1,iTMP) = COOR(1,I)
               COOR(2,iTMP) = COOR(2,I)
               COOR(3,iTMP) = COOR(3,I)
               WG(iTMP) = WG(I)
            ENDIF
         ENDDO
      ENDIF
      totalAtomicpointscompressed = ITMP
   endif
   IF(iprint.gt. 1)THEN
!$OMP CRITICAL (GRIDWRITE)
      WRITE(Lupri,'(A,I5,A,I6,A,I6,A,I5,A)')'Atom: ',IATOM,' points=',&
           & totalAtomicpointscompressed,' compressed from',&
           & totalAtomicpoints,'(',nRadialPoints(IATOM),' radial)'
!$OMP END CRITICAL (GRIDWRITE)
   ENDIF
!$OMP CRITICAL 
   privatetotalpoints = totalpoints
   totalpoints = totalpoints + totalAtomicpointscompressed   
!$OMP END CRITICAL 
! We only write once to each element so this does not need to be critical
! But we do need to know where in the list to write to 
   DO I=1,totalAtomicpointscompressed
      COOR2(1,privatetotalpoints + I) = COOR(1,I)
      COOR2(2,privatetotalpoints + I) = COOR(2,I)
      COOR2(3,privatetotalpoints + I) = COOR(3,I)
   ENDDO
   DO I=1,totalAtomicpointscompressed
      WG2(privatetotalpoints + I) = WG(I)
   ENDDO
   IF(postprocess.AND.PARTITIONING.EQ. 5)THEN
      DO I=1,totalAtomicpointscompressed
         skip2(privatetotalpoints + I) = skip(I)
      ENDDO
   ENDIF
   IF(postprocess)THEN
      DO I=privatetotalpoints+1,privatetotalpoints+totalAtomicpointscompressed
         ATOMIDX(I) = IATOM
      ENDDO
   ENDIF
ENDDO

call mem_grid_dealloc(locR)
call mem_grid_dealloc(COOR)
call mem_grid_dealloc(WG)
IF(PARTITIONING.EQ. 5)call mem_grid_dealloc(skip)

!the calculation of totalpoints must be done before the master starts work
!so we need a barrier 
!$OMP BARRIER

!$OMP MASTER

WRITE(Lupri,'(A,I12)')'Total Number of grid points:',totalpoints

!build the Keypointer and determines nkey
call BuildBoxes(gridBox,COOR2,GlobalmaxGridpoints,totalpoints,RSHEL,&
     & MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,atomcenterX,atomcenterY,atomcenterZ,&
     & NSTART,Maxnbuflen,nkey,keypointer,iprint,lupri)
IF(postprocess)then
   call AddAtomicidxToBox(keypointer,nkey,ATOMIDX,GlobalmaxGridpoints)
ENDIF
!We now write the gridpoints to file(s)
!In case of MPI we use one file for each processor with a unique name
LUGRID=-1
call get_quadfilename(filename,nbast,node,GridId)
CALL LSOPEN(LUGRID,filename,'NEW','UNFORMATTED')

#ifdef VAR_MPI
IF(USE_MPI)THEN 
   nprocessors = numnodes
   mynum = node
ELSE
   nprocessors = 1
   mynum = 0
ENDIF
#else
nprocessors = 1
mynum = 0
#endif

!MPI determination of which keys/boxes it should do
!so that OpenMP can be done on the loop later
call mem_grid_alloc(mykeys,nkey)
mykeysNumber=0
DO Ikey=1+mynum,nkey,nprocessors
   mykeysNumber=mykeysNumber+1
   mykeys(mykeysNumber)=Ikey  
ENDDO

!$OMP END MASTER

!the remaining calculation require stuff calculated by the master
!so we need a barrier 
!$OMP BARRIER

MaxNactBAST = 0
call mem_grid_alloc(SHELLBLOCKS,2,MAXNSHELL)
!call mem_grid_alloc(SHELLLIST,MAXNSHELL)
!DO I=1,MAXNSHELL
!   SHELLLIST(I) = .FALSE.
!ENDDO
! Another OpenMP loop 
call mem_grid_alloc(WG3,maxNBUFLEN)
call mem_grid_alloc(COOR3,3,maxNBUFLEN)
DO IMykey=1+tid,mykeysNumber,nthreads
   Ikey = mykeys(IMykey)
   npoints2=keypointer(Ikey)%npoints
!   IF(npoints2.EQ. 0)CYCLE !this should no longer be necesarry 
   CALL GRID_GETBLOCKS(RSHEL,MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,atomcenterX,atomcenterY,&
        & atomcenterZ,NSHELLBLOCKS,SHELLBLOCKS,NactBAST,Nstart,keypointer(Ikey)%Center,&
        & keypointer(Ikey)%CELL_SIZE,LUPRI)
!   CALL GRID_BUILDSHELLLIST(SHELLLIST,MAXNSHELL,NSHELLBLOCKS,SHELLBLOCKS,LUPRI)
   MaxNactBAST = MAX(MaxNactBAST,NactBAST)
   IF(NSHELLBLOCKS.EQ. 0)CYCLE
   IF(postprocess.AND.PARTITIONING.EQ. 5)THEN 
      CALL BlockSSFPartitioning(keypointer(Ikey)%Center,NATOMS,atomcenterX,atomcenterY,&
           & atomcenterZ,COOR2,inverseDistance12,ihardness,&
           & GlobalmaxGridpoints,npoints2,keypointer(Ikey)%points(1:npoints2),&
           & keypointer(Ikey)%atom_idx(1:npoints2),WG2,maxNBUFLEN,&
           & ITERATIONS,RADIALPOINTS,NRADPT,weight_thr,LUGRID,NSHELLBLOCKS,&
           & SHELLBLOCKS(1:2,1:NSHELLBLOCKS),keypointer(Ikey)%CELL_SIZE,skip2,COOR3,WG3)
   ELSEIF(postprocess.AND.PARTITIONING.EQ. 4)THEN 
      CALL BlockPartitioning(keypointer(Ikey)%Center,NATOMS,atomcenterX,atomcenterY,&
           & atomcenterZ,COOR2,inverseDistance12,bfac,ihardness,&
           & GlobalmaxGridpoints,npoints2,keypointer(Ikey)%points(1:npoints2),&
           & keypointer(Ikey)%atom_idx(1:npoints2),WG2,maxNBUFLEN,&
           & ITERATIONS,RADIALPOINTS,NRADPT,weight_thr,LUGRID,NSHELLBLOCKS,&
           & SHELLBLOCKS(1:2,1:NSHELLBLOCKS),keypointer(Ikey)%CELL_SIZE,COOR3,WG3)
   ELSE
      ! if npoints greater than NBUFLEN loop over npoints in batches of maxNBUFLEN
      if(npoints2.GT.maxNBUFLEN)THEN
         CALL WRITE_COORD2(npoints2,ITERATIONS,keypointer(Ikey)%points(1:npoints2),&
              & COOR2,WG2,GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS(1:2,1:NSHELLBLOCKS),&
              & maxNBUFLEN,LUGRID,COOR3,WG3)
      else
         CALL WRITE_COORD1(npoints2,ITERATIONS,keypointer(Ikey)%points(1:npoints2),&
              & COOR2,WG2,GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS(1:2,1:NSHELLBLOCKS),&
              & maxNBUFLEN,LUGRID,COOR3,WG3)
      endif
   ENDIF
   ! In principle I should be able to call mem_dealloc here, but 
   ! since this is a Shared variable which is allocated by the master 
   ! the memory booking gets fucked up if someone other than master 
   ! deallocates it.
   deallocate(keypointer(Ikey)%points)
   IF(postprocess)THEN
      deallocate(keypointer(ikey)%atom_idx)
   ENDIF
ENDDO
call mem_grid_dealloc(WG3)
call mem_grid_dealloc(COOR3)

!call SHELLLIST_ORBBLOCK(LUPRI,SHELLLIST,MAXNSHELL,SHELLBLOCKS)

call mem_grid_dealloc(SHELLBLOCKS)
!call mem_grid_dealloc(SHELLLIST)
!$OMP CRITICAL
FinalMaxNactBast = MAX(FinalMaxNactBast,MaxNactBast)
!$OMP END CRITICAL
call collect_thread_grid_memory()
!$OMP END PARALLEL
call mem_grid_TurnOffThread_Memory()

call mem_grid_dealloc(mykeys)
deallocate(keypointer) 

IF(postprocess)THEN
   call mem_grid_dealloc(ATOMIDX)   
ENDIF
call mem_grid_dealloc(inverseDistance12)
IF(PARTITIONING.EQ. 2.OR.PARTITIONING.EQ. 4)THEN !needed for becke preprocess
   call mem_grid_dealloc(Bfac)
ENDIF
IF(PARTITIONING.EQ. 5.OR.PARTITIONING.EQ. 1)THEN
   call mem_grid_dealloc(SSF_ATOMIC_CUTOFF)
ENDIF

call mem_grid_dealloc(WG2)
call mem_grid_dealloc(COOR2)
IF(PARTITIONING.EQ. 5)THEN
   call mem_grid_dealloc(skip2)   
ENDIF

CALL LSCLOSE(LUGRID,'KEEP')

END SUBROUTINE Computecoords

SUBROUTINE WRITE_COORD1(npoints2,ITERATIONS,points,COOR2,WG2,&
     & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,maxNBUFLEN,&
     & LUGRID,COOR3,WG3)
  implicit none
  INTEGER,intent(in) :: npoints2,GlobalmaxGridpoints,NSHELLBLOCKS,maxNBUFLEN,LUGRID
  INTEGER,intent(inout) :: ITERATIONS
  integer,intent(in) :: points(npoints2)
  REAL(REALK),intent(in) :: COOR2(3,GlobalmaxGridpoints),WG2(GlobalmaxGridpoints)
  integer,intent(in) :: SHELLBLOCKS(2,NSHELLBLOCKS)
!
  real(realk),intent(inout) :: COOR3(3,maxNBUFLEN),WG3(maxNBUFLEN)
  Integer :: I,Ipoint
  DO Ipoint=1,npoints2
     I = points(ipoint)
     COOR3(1,Ipoint)=COOR2(1,I)
     COOR3(2,Ipoint)=COOR2(2,I)
     COOR3(3,Ipoint)=COOR2(3,I)
     WG3(Ipoint)=WG2(I)
  ENDDO
!$OMP CRITICAL
  ITERATIONS=ITERATIONS+1
  WRITE(LUGRID) npoints2
  WRITE(LUGRID) COOR3(1:3,1:npoints2)
  WRITE(LUGRID) WG3(1:npoints2)
  WRITE(LUGRID) NSHELLBLOCKS
  WRITE(LUGRID) SHELLBLOCKS(1:2,1:NSHELLBLOCKS)
!$OMP END CRITICAL
END SUBROUTINE WRITE_COORD1

SUBROUTINE WRITE_COORD2(npoints2,ITERATIONS,points,COOR2,WG2,&
     & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,maxNBUFLEN,LUGRID,COOR3,WG3)
  implicit none
  INTEGER,intent(in) :: npoints2,GlobalmaxGridpoints,NSHELLBLOCKS,maxNBUFLEN,LUGRID
  INTEGER,intent(inout) :: ITERATIONS
  integer,intent(in) :: points(npoints2)
  REAL(REALK),intent(in) :: COOR2(3,GlobalmaxGridpoints),WG2(GlobalmaxGridpoints)
  integer,intent(in) :: SHELLBLOCKS(2,NSHELLBLOCKS)
  real(realk),intent(inout) :: COOR3(3,maxNBUFLEN),WG3(maxNBUFLEN)
!
  Integer :: I,Ipoint,NGRIDPOINTINBATCH,IT
  
  DO Ipoint=1,npoints2,maxNBUFLEN
     NGRIDPOINTINBATCH=MIN(maxNBUFLEN,npoints2-ipoint+1)
     DO IT=0,NGRIDPOINTINBATCH-1
        I = points(ipoint+IT)
        COOR3(1,IT+1)=COOR2(1,I)
        COOR3(2,IT+1)=COOR2(2,I)
        COOR3(3,IT+1)=COOR2(3,I)
        WG3(IT+1)=WG2(I)
     ENDDO
!$OMP CRITICAL
     ITERATIONS=ITERATIONS+1
     WRITE(LUGRID) NGRIDPOINTINBATCH
     WRITE(LUGRID) COOR3(1:3,1:NGRIDPOINTINBATCH)
     WRITE(LUGRID) WG3(1:NGRIDPOINTINBATCH)
     WRITE(LUGRID) NSHELLBLOCKS
     WRITE(LUGRID) SHELLBLOCKS(1:2,1:NSHELLBLOCKS)
!$OMP END CRITICAL
  ENDDO
END SUBROUTINE WRITE_COORD2

SUBROUTINE WRITE_COORD3(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
     & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,ATOM_IDX,weight_thr,&
     & maxNBUFLEN,LUGRID,COOR3,WG3)
  implicit none
  INTEGER,intent(in) :: npoints,GlobalmaxGridpoints,NSHELLBLOCKS,LUGRID,NATOMS,maxNBUFLEN
  INTEGER,intent(inout) :: ITERATIONS
  integer,intent(in) :: points(npoints),atom_idx(npoints)
  REAL(REALK),intent(in) :: COOR2(3,GlobalmaxGridpoints),WG2(GlobalmaxGridpoints)
  REAL(REALK),intent(in) :: weight_thr
  integer,intent(in) :: SHELLBLOCKS(2,NSHELLBLOCKS)
  real(realk),intent(in) :: MU(npoints),TMP2(npoints,NATOMS)
  real(realk),intent(inout) :: COOR3(3,maxNBUFLEN),WG3(maxNBUFLEN)
  real(realk) :: FACTOR
  Integer :: I,Ipoint,I3,IATOM
  I=0
  DO Ipoint=1,npoints
     IATOM = atom_idx(Ipoint)
     FACTOR = TMP2(Ipoint,IATOM)/MU(Ipoint) 
     IF(FACTOR.GT.weight_thr)THEN
        I=I+1
        I3 = points(ipoint)
        WG3(I) = WG2(I3)*FACTOR
        COOR3(1,I) = COOR2(1,I3)
        COOR3(2,I) = COOR2(2,I3)
        COOR3(3,I) = COOR2(3,I3)
     ENDIF
  ENDDO
  IF(I.GT. 0)THEN
!$OMP CRITICAL
     ITERATIONS=ITERATIONS+1
     WRITE(LUGRID) I
     WRITE(LUGRID) COOR3(1:3,1:I)
     WRITE(LUGRID) WG3(1:I)
     WRITE(LUGRID) NSHELLBLOCKS
     WRITE(LUGRID) SHELLBLOCKS(1:2,1:NSHELLBLOCKS)
!$OMP END CRITICAL
  ENDIF
END SUBROUTINE WRITE_COORD3

SUBROUTINE WRITE_COORD4(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
     &GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,ATOM_IDX,weight_thr,&
     &maxNBUFLEN,LUGRID,COOR3,WG3)
  implicit none
  INTEGER,intent(in) :: npoints,GlobalmaxGridpoints,NSHELLBLOCKS,LUGRID,NATOMS
  integer,intent(in) :: points(npoints),atom_idx(npoints),maxNBUFLEN
  INTEGER,intent(inout) :: ITERATIONS
  REAL(REALK),intent(in) :: COOR2(3,GlobalmaxGridpoints),WG2(GlobalmaxGridpoints)
  integer,intent(in) :: SHELLBLOCKS(2,NSHELLBLOCKS)
  real(realk),intent(in) :: MU(npoints),TMP2(npoints,NATOMS),weight_thr
!
  real(realk),intent(inout) :: COOR3(3,maxNBUFLEN),WG3(maxNBUFLEN)
  real(realk) :: FACTOR
  Integer :: I,Ipoint,I3,NGRIDPOINTINBATCH,IT,IATOM
  !FIXME CHANGE THIS SO 
!  DO Ipoint=1,npoints
!
!        IF(FACTOR.GT.weight_thr)THEN
!           I=I+1
!   IF(I.EQ.maxNBUFLEN)call WRITEBUFFER
!
!  ENDDO
!  call WRITEBUFFER
  DO Ipoint=1,npoints,maxNBUFLEN
     NGRIDPOINTINBATCH=MIN(maxNBUFLEN,npoints-ipoint+1)
     I=0
     DO IT=0,NGRIDPOINTINBATCH-1
        IATOM = atom_idx(Ipoint+IT)
        FACTOR = TMP2(Ipoint+IT,IATOM)/MU(Ipoint+IT) 
        IF(FACTOR.GT.weight_thr)THEN
           I=I+1
           I3=points(ipoint+IT)
           WG3(I) = WG2(I3)*FACTOR
           COOR3(1,I) = COOR2(1,I3)
           COOR3(2,I) = COOR2(2,I3)
           COOR3(3,I) = COOR2(3,I3)
        ENDIF
     ENDDO
     IF(I.GT. 0)THEN
!$OMP CRITICAL
        ITERATIONS=ITERATIONS+1
        WRITE(LUGRID) I
        WRITE(LUGRID) COOR3(1:3,1:I)
        WRITE(LUGRID) WG3(1:I)
        WRITE(LUGRID) NSHELLBLOCKS
        WRITE(LUGRID) SHELLBLOCKS(1:2,1:NSHELLBLOCKS)
!$OMP END CRITICAL
     ENDIF
  ENDDO
END SUBROUTINE WRITE_COORD4

SUBROUTINE WRITE_COORD5(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
     & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,ATOM_IDX,weight_thr,maxNBUFLEN,LUGRID,&
     & relevantpoints,nrelevantpoints,nWeight1points,Weight1points,COOR3,WG3)
  implicit none
  INTEGER,intent(in) :: npoints,GlobalmaxGridpoints,NSHELLBLOCKS,LUGRID,NATOMS,maxNBUFLEN
  INTEGER,intent(inout) :: ITERATIONS
  integer,intent(in) :: points(npoints),atom_idx(npoints)
  REAL(REALK),intent(in) :: COOR2(3,GlobalmaxGridpoints),WG2(GlobalmaxGridpoints)
  REAL(REALK),intent(in) :: weight_thr
  integer,intent(in) :: SHELLBLOCKS(2,NSHELLBLOCKS)
  real(realk),intent(in) :: MU(npoints),TMP2(npoints,NATOMS)
  integer,intent(in) :: nrelevantpoints,nWeight1points
  integer,intent(in) :: relevantpoints(nrelevantpoints)
  integer,intent(in) :: Weight1points(nWeight1points)
  real(realk),intent(inout) :: COOR3(3,maxNBUFLEN),WG3(maxNBUFLEN)
!
  real(realk) :: FACTOR
  Integer :: I,Ipoint,I3,IATOM,Ipoint_rel
  I=0
  DO Ipoint_rel=1,nrelevantpoints!npoints
     Ipoint = relevantpoints(Ipoint_rel)
     IATOM = atom_idx(Ipoint)
     FACTOR = TMP2(Ipoint_rel,IATOM)/MU(Ipoint_rel) 
     IF(FACTOR.GT.weight_thr)THEN
        I=I+1
        I3 = points(ipoint)
        WG3(I) = WG2(I3)*FACTOR
        COOR3(1,I) = COOR2(1,I3)
        COOR3(2,I) = COOR2(2,I3)
        COOR3(3,I) = COOR2(3,I3)
     ENDIF
  ENDDO
  DO Ipoint_rel=1,nWeight1points
     Ipoint = Weight1points(Ipoint_rel)
     IATOM = atom_idx(Ipoint)
     I=I+1
     I3 = points(ipoint)
     WG3(I) = WG2(I3)
     COOR3(1,I) = COOR2(1,I3)
     COOR3(2,I) = COOR2(2,I3)
     COOR3(3,I) = COOR2(3,I3)
  ENDDO
  IF(I.GT. 0)THEN
!$OMP CRITICAL
     ITERATIONS=ITERATIONS+1
     WRITE(LUGRID) I
     WRITE(LUGRID) COOR3(1:3,1:I)
     WRITE(LUGRID) WG3(1:I)
     WRITE(LUGRID) NSHELLBLOCKS
     WRITE(LUGRID) SHELLBLOCKS(1:2,1:NSHELLBLOCKS)
!$OMP END CRITICAL
  ENDIF
END SUBROUTINE WRITE_COORD5

SUBROUTINE WRITE_COORD6(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
     & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,ATOM_IDX,weight_thr,maxNBUFLEN,LUGRID,&
     & relevantpoints,nrelevantpoints,nWeight1points,Weight1points,COOR3,WG3)
  implicit none
  INTEGER,intent(in) :: npoints,GlobalmaxGridpoints,NSHELLBLOCKS,LUGRID,NATOMS
  integer,intent(in) :: points(npoints),atom_idx(npoints),maxNBUFLEN
  INTEGER,intent(inout) :: ITERATIONS
  REAL(REALK),intent(in) :: COOR2(3,GlobalmaxGridpoints),WG2(GlobalmaxGridpoints)
  integer,intent(in) :: SHELLBLOCKS(2,NSHELLBLOCKS)
  real(realk),intent(in) :: MU(npoints),TMP2(npoints,NATOMS),weight_thr
  integer,intent(in) :: nrelevantpoints,nWeight1points
  integer,intent(in) :: relevantpoints(nrelevantpoints)
  integer,intent(in) :: Weight1points(nWeight1points)
  real(realk),intent(inout) :: COOR3(3,maxNBUFLEN),WG3(maxNBUFLEN)
!
  real(realk) :: FACTOR
  Integer :: I,Ipoint,I3,NGRIDPOINTINBATCH,IT,IATOM,Ipoint_rel
  
  DO Ipoint_rel=1,nrelevantpoints,maxNBUFLEN
     NGRIDPOINTINBATCH=MIN(maxNBUFLEN,nrelevantpoints-ipoint_rel+1)
     I=0
     DO IT=0,NGRIDPOINTINBATCH-1
        Ipoint = relevantpoints(Ipoint_rel+IT)
        IATOM = atom_idx(Ipoint)
        FACTOR = TMP2(Ipoint_rel+IT,IATOM)/MU(Ipoint_rel+IT) 
        IF(FACTOR.GT.weight_thr)THEN
           I=I+1
           I3=points(ipoint)
           WG3(I) = WG2(I3)*FACTOR
           COOR3(1,I) = COOR2(1,I3)
           COOR3(2,I) = COOR2(2,I3)
           COOR3(3,I) = COOR2(3,I3)
        ENDIF
     ENDDO
     IF(I.GT. 0)THEN
!$OMP CRITICAL
        ITERATIONS=ITERATIONS+1
        WRITE(LUGRID) I
        WRITE(LUGRID) COOR3(1:3,1:I)
        WRITE(LUGRID) WG3(1:I)
        WRITE(LUGRID) NSHELLBLOCKS
        WRITE(LUGRID) SHELLBLOCKS(1:2,1:NSHELLBLOCKS)
!$OMP END CRITICAL
     ENDIF
  ENDDO
  DO Ipoint_rel=1,nWeight1points,maxNBUFLEN
     NGRIDPOINTINBATCH=MIN(maxNBUFLEN,nWeight1points-ipoint_rel+1)
     I=0
     DO IT=0,NGRIDPOINTINBATCH-1
        Ipoint = Weight1points(Ipoint_rel+IT)
        IATOM = atom_idx(Ipoint)
        I=I+1
        I3=points(ipoint)
        WG3(I) = WG2(I3)
        COOR3(1,I) = COOR2(1,I3)
        COOR3(2,I) = COOR2(2,I3)
        COOR3(3,I) = COOR2(3,I3)
     ENDDO
     IF(I.GT. 0)THEN
!$OMP CRITICAL
        ITERATIONS=ITERATIONS+1
        WRITE(LUGRID) I
        WRITE(LUGRID) COOR3(1:3,1:I)
        WRITE(LUGRID) WG3(1:I)
        WRITE(LUGRID) NSHELLBLOCKS
        WRITE(LUGRID) SHELLBLOCKS(1:2,1:NSHELLBLOCKS)
!$OMP END CRITICAL
     ENDIF
  ENDDO
END SUBROUTINE WRITE_COORD6

SUBROUTINE READ_GRIDPOINTS(LUGRID,NSHELL,SHELLBLOCKS,MAXNSHELL,NBUFLEN,COOR,WEIGHT,NLEN)
use precision
IMPLICIT NONE
integer,intent(in)    :: LUGRID,NBUFLEN,MAXNSHELL
integer,intent(inout) :: NLEN,NSHELL
INTEGER,intent(inout) :: SHELLBLOCKS(2,MAXNSHELL)
REAL(REALK),intent(inout) :: COOR(3,NBUFLEN),WEIGHT(NBUFLEN)
!
INTEGER :: I,layer
READ(LUGRID) NLEN
READ(LUGRID) COOR(1:3,1:NLEN)
READ(LUGRID) WEIGHT(1:NLEN)
READ(LUGRID) NSHELL
READ(LUGRID) SHELLBLOCKS(1:2,1:NSHELL)

END SUBROUTINE READ_GRIDPOINTS

!!$SUBROUTINE DETERMINE_NACTIVEORB_FROM_SHELL(NactBAST,NSHELLBLOCKS,SHELLBLOCKS,MAXNSHELL,NSTART,NBAST)
!!$implicit none
!!$INTEGER :: NactBAST,NSHELLBLOCKS,MAXNSHELL,NBAST
!!$INTEGER :: NSTART(MAXNSHELL),SHELLBLOCKS(2,MAXNSHELL)
!!$!
!!$INTEGER :: ORBBLOCKS(2,MAXNSHELL),I,IBL,IORB
!!$
!!$DO I = 1,NSHELLBLOCKS
!!$   ORBBLOCKS(1,I) = NSTART(SHELLBLOCKS(1,I))+1
!!$ENDDO
!!$
!!$DO I = 1,NSHELLBLOCKS
!!$   IF(SHELLBLOCKS(2,I) .LT. MAXNSHELL)THEN
!!$      ORBBLOCKS(2,I) = NSTART(SHELLBLOCKS(2,I)+1)
!!$   ELSE
!!$      ORBBLOCKS(2,I) = NBAST
!!$   ENDIF
!!$ENDDO
!!$
!!$NactBAST = 0
!!$DO IBL = 1, NSHELLBLOCKS
!!$   DO IORB = ORBBLOCKS(1,IBL),ORBBLOCKS(2,IBL)
!!$      NactBAST = NactBAST + 1
!!$   ENDDO
!!$ENDDO
!!$
!!$END SUBROUTINE DETERMINE_NACTIVEORB_FROM_SHELL

SUBROUTINE determine_distance12(inverseDistance12,NATOMS,atomcenterX,atomcenterY,atomcenterZ)
IMPLICIT NONE
INTEGER,intent(in) :: NATOMS
real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS)
real(realk),intent(inout) :: inverseDistance12(NATOMS,NATOMS)
!
Integer :: IATOM1,IATOM2
real(realk) :: Xatomcoor,Yatomcoor,Zatomcoor,DX,DY,DZ
real(realk),parameter :: D1=1E0_realk,D0=0E0_realk
 
DO IATOM1=1,NATOMS
   Xatomcoor = atomcenterX(IATOM1)
   Yatomcoor = atomcenterY(IATOM1)
   Zatomcoor = atomcenterZ(IATOM1)   
   DO IATOM2=1,IATOM1-1
      inverseDistance12(IATOM2,IATOM1) = inverseDistance12(IATOM1,IATOM2)
   ENDDO
   inverseDistance12(IATOM1,IATOM1) = D0
   DO IATOM2=IATOM1+1,NATOMS
      DX = Xatomcoor-atomcenterX(IATOM2)
      DY = Yatomcoor-atomcenterY(IATOM2)
      DZ = Zatomcoor-atomcenterZ(IATOM2)
      inverseDistance12(IATOM2,IATOM1) = D1/SQRT(DX*DX+DY*DY+DZ*DZ)
   ENDDO
ENDDO
END SUBROUTINE DETERMINE_DISTANCE12

SUBROUTINE determine_SSF_ATOMIC_CUTOFF(NATOMS,SSF_ATOMIC_CUTOFF,&
     & atomcenterX,atomcenterY,atomcenterZ)
IMPLICIT NONE
INTEGER,intent(in) :: NATOMS
real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS)
real(realk),intent(inout) :: SSF_ATOMIC_CUTOFF(NATOMS)
!
Integer :: IATOM1,IATOM2
real(realk) :: Xatomcoor,Yatomcoor,Zatomcoor,DX,DY,DZ,dist,dist2
real(realk),parameter :: D1=1E0_realk,D0=0E0_realk,D05=0.5E0_realk
 
DO IATOM1=1,NATOMS
   Xatomcoor = atomcenterX(IATOM1)
   Yatomcoor = atomcenterY(IATOM1)
   Zatomcoor = atomcenterZ(IATOM1)   
   dist=HUGE(1E0_realk)!1.0E+10_realk
   DO IATOM2=1,NATOMS
      IF(IATOM2.EQ.IATOM1)CYCLE
      DX = Xatomcoor-atomcenterX(IATOM2)
      DY = Yatomcoor-atomcenterY(IATOM2)
      DZ = Zatomcoor-atomcenterZ(IATOM2)
      dist2 = SQRT(DX*DX+DY*DY+DZ*DZ)
      IF(dist .GE. dist2) dist = dist2
   ENDDO
   SSF_ATOMIC_CUTOFF(IATOM1) = D05 * (D1-SSF_CUTOFF)*dist
ENDDO
END SUBROUTINE DETERMINE_SSF_ATOMIC_CUTOFF

!based on jcp vol 88, page 2547 (1988)
Subroutine BeckeOriginalPartitioning(N,NATOMS,atomcenterX,atomcenterY,&
     & atomcenterZ,COOR,idx,inverseDistance12,ihardness,MaxAtomicGridpoints,WG,IATOM)
  implicit none
  integer,intent(in)     :: N,NATOMS,idx,ihardness,IATOM,MaxAtomicGridpoints
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS)
  real(realk),intent(in) :: COOR(3,MaxAtomicGridpoints)
  real(realk),intent(in) :: inverseDistance12(NATOMS,NATOMS)
  real(realk),intent(inout) :: WG(MaxAtomicGridpoints)
  !
  real(realk),parameter :: D1=1E0_realk,D05=0.5E0_realk,D3=3E0_realk
  real(realk),pointer :: TMP(:,:),TMP2(:,:),MU(:),MU2(:)
  real(realk) :: dist
  real(realk) :: Xatomcoor2,Yatomcoor2,Zatomcoor2,DX,DY,DZ,MUA
  integer     :: IATOM2,I2,IATOM1,h
  call mem_grid_alloc(TMP,N,NATOMS)
  call mem_grid_alloc(TMP2,N,NATOMS)
  call mem_grid_alloc(MU,N)
  call mem_grid_alloc(MU2,N)
  CALL DetermineGridDistance(TMP,atomcenterX,atomcenterY,atomcenterZ,&
       & NATOMS,N,COOR,MaxAtomicGridpoints,idx)
  call ls_SetToOne(TMP2,NATOMS*N)
  DO IATOM1=2,NATOMS
     DO IATOM2=1,IATOM1-1
        dist = inverseDistance12(IATOM2,IATOM1)
        DO I2=1,N
           MUA = (TMP(I2,IATOM1)-TMP(I2,IATOM2))*dist
           MU(I2)=MUA
           MU2(I2)=MUA*MUA
        ENDDO
        DO h=1,ihardness-1                     !Eq. 20 of jcp 88, 2547
           DO I2=1,N
              mu(I2) = D05*mu(I2)*(D3-mu2(I2)) !Eq. 19 of jcp 88, 2547
              mu2(I2) = mu(I2)*mu(I2)
           ENDDO
        ENDDO
        h=ihardness
        DO I2=1,N
           mu(I2) = D05*mu(I2)*(D3-mu2(I2))    !Eq. 19 of jcp 88, 2547
        ENDDO
        DO I2=1,N
           TMP2(I2,IATOM1) = TMP2(I2,IATOM1)*D05*(D1-mu(I2)) !Eq. 20 of jcp 88, 2547
           TMP2(I2,IATOM2) = TMP2(I2,IATOM2)*D05*(D1+mu(I2)) 
        ENDDO
     ENDDO
  ENDDO
  DO I2=1,N
     MU(I2)=TMP2(I2,1)
  ENDDO
  DO IATOM1=2,NATOMS
     DO I2=1,N
        MU(I2)=MU(I2) + TMP2(I2,IATOM1)
     ENDDO
  ENDDO
  DO I2=1,N
     WG(I2+idx) = WG(I2+idx)*TMP2(I2,IATOM)/MU(I2) 
  ENDDO
  call mem_grid_dealloc(TMP)
  call mem_grid_dealloc(TMP2)
  call mem_grid_dealloc(MU)
  call mem_grid_dealloc(MU2)
END Subroutine BECKEORIGINALPARTITIONING

!based on jcp vol 88, page 2547 (1988)
Subroutine BeckePartitioning(N,NATOMS,atomcenterX,atomcenterY,&
     & atomcenterZ,COOR,idx,inverseDistance12,BFAC,ihardness,MaxAtomicGridpoints,WG,IATOM)
  implicit none
  integer,intent(in)     :: N,NATOMS,idx,ihardness,IATOM,MaxAtomicGridpoints
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS)
  real(realk),intent(in) :: COOR(3,MaxAtomicGridpoints)
  real(realk),intent(in) :: inverseDistance12(NATOMS,NATOMS),BFAC(NATOMS,NATOMS)
  real(realk),intent(inout) :: WG(MaxAtomicGridpoints)
  !
  real(realk),parameter :: D1=1E0_realk,D05=0.5E0_realk,D3=3E0_realk
  real(realk),pointer :: TMP(:,:),TMP2(:,:),MU(:),MU2(:)
  real(realk) :: Xatomcoor2,Yatomcoor2,Zatomcoor2,DX,DY,DZ,MUA,B,dist
  integer     :: IATOM2,I2,IATOM1,h
  call mem_grid_alloc(TMP,N,NATOMS)
  call mem_grid_alloc(TMP2,N,NATOMS)
  call mem_grid_alloc(MU,N)
  call mem_grid_alloc(MU2,N)

  CALL DetermineGridDistance(TMP,atomcenterX,atomcenterY,atomcenterZ,&
       & NATOMS,N,COOR,MaxAtomicGridpoints,idx)
  call ls_SetToOne(TMP2,NATOMS*N)
  DO IATOM1=2,NATOMS
     DO IATOM2=1,IATOM1-1
        dist = inverseDistance12(IATOM2,IATOM1)
        B = BFAC(IATOM2,IATOM1)
        DO I2=1,N
           MUA = (TMP(I2,IATOM1)-TMP(I2,IATOM2))*dist
           MUA = MUA + B*(D1-MUA*MUA);
           MU(I2)=MUA
           MU2(I2)=MUA*MUA
        ENDDO
        DO h=1,ihardness-1                     !Eq. 20 of jcp 88, 2547
           DO I2=1,N
              mu(I2) = D05*mu(I2)*(D3-mu2(I2)) !Eq. 19 of jcp 88, 2547
              mu2(I2) = mu(I2)*mu(I2)
           ENDDO
        ENDDO
        h=ihardness
        DO I2=1,N
           mu(I2) = D05*mu(I2)*(D3-mu2(I2))    !Eq. 19 of jcp 88, 2547
        ENDDO
        DO I2=1,N
           TMP2(I2,IATOM1) = TMP2(I2,IATOM1)*D05*(D1-mu(I2)) !Eq. 20 of jcp 88, 2547
           TMP2(I2,IATOM2) = TMP2(I2,IATOM2)*D05*(D1+mu(I2)) 
        ENDDO
     ENDDO
  ENDDO
  DO I2=1,N
     MU(I2)=TMP2(I2,1)
  ENDDO
  DO IATOM1=2,NATOMS
     DO I2=1,N
        MU(I2)=MU(I2) + TMP2(I2,IATOM1)
     ENDDO
  ENDDO
  DO I2=1,N
     WG(I2+idx) = WG(I2+idx)*TMP2(I2,IATOM)/MU(I2) 
  ENDDO
  call mem_grid_dealloc(TMP)
  call mem_grid_dealloc(TMP2)
  call mem_grid_dealloc(MU)
  call mem_grid_dealloc(MU2)
END Subroutine BECKEPARTITIONING

!based on Chem. Phys. Lett. 1996, 257, 213
Subroutine SSFPartitioning(N,NATOMS,atomcenterX,atomcenterY,&
     & atomcenterZ,COOR,idx,inverseDistance12,ihardness,MaxAtomicGridpoints,WG,IATOM)
  implicit none
  integer,intent(in)     :: N,NATOMS,idx,ihardness,IATOM,MaxAtomicGridpoints
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS)
  real(realk),intent(in) :: COOR(3,MaxAtomicGridpoints)
  real(realk),intent(in) :: inverseDistance12(NATOMS,NATOMS)
  real(realk),intent(inout) :: WG(MaxAtomicGridpoints)
  !
  real(realk),parameter :: D1=1E0_realk,D05=0.5E0_realk,D3=3E0_realk,D2=2E0_realk,D0=0E0_realk
! some constants for the SSF scheme
  real(realk) :: dist,radius,Ria1,W,V
  real(realk) :: Xatomcoor2,Yatomcoor2,Zatomcoor2,DX,DY,DZ,MUA,MUA2,GMU
  integer     :: IATOM2,I2,IATOM1,h,I,ATOMINDEX
  Integer     :: nRELEVANT_ATOMS,IATOM1_REL,Iatom2_REL
  integer,pointer     :: RELEVANT_ATOMS(:)
  real(realk),pointer :: TMP(:,:),P(:)
! Testing 1
  Logical,pointer     :: zero(:)
  call mem_grid_alloc(TMP,N,NATOMS)
  call mem_grid_alloc(P,MAX(N,NATOMS))
  call mem_grid_alloc(RELEVANT_ATOMS,NATOMS)
  call mem_grid_alloc(zero,NATOMS)

  DX = atomcenterX(IATOM) - COOR(1,1+idx)
  DY = atomcenterY(IATOM) - COOR(2,1+idx)
  DZ = atomcenterZ(IATOM) - COOR(3,1+idx)
  RADIUS=SQRT(DX*DX + DY*DY + DZ*DZ) !radius of the sphere 
  ! we right now do not make use of the fact that the outermost points come first, 
  ! so that the list of relevant atoms becomes shorter in the course of evaluation
  ! and we can save the list of relevant atoms and reduce the list instead.
  ! We now determine the relevant atoms
  ! First up to IATOM
  I=0
  DO IATOM1=1,IATOM-1
     dist = inverseDistance12(IATOM1,IATOM)
     IF(D2*RADIUS*dist.GT.D1-SSF_CUTOFF)THEN
        I=I+1
        RELEVANT_ATOMS(I) = IATOM1
     ENDIF
  ENDDO
  !We always include IATOM itself as a relevant atom
  I=I+1
  ATOMINDEX = I !Index for when in the list of relevant atoms
                !the atoms are .LT.IATOM 
  RELEVANT_ATOMS(I) = IATOM
  ! Now we continue through the rest of the atoms
  DO IATOM1=IATOM+1,NATOMS
     dist = inverseDistance12(IATOM1,IATOM)
     IF(D2*RADIUS*dist.GT.D1-SSF_CUTOFF)THEN
        I=I+1
        RELEVANT_ATOMS(I) = IATOM1
     ENDIF
  ENDDO
  nRELEVANT_ATOMS = I
  IF(nRELEVANT_ATOMS.EQ. 1)RETURN !no need to do anymore work
  
! compute confocal ellipical coordinates for the batch of points to
! be processed. 
  DO IATOM2_REL=1,nRELEVANT_ATOMS
     IATOM2 = RELEVANT_ATOMS(IATOM2_REL)
     Xatomcoor2 = -atomcenterX(IATOM2)
     Yatomcoor2 = -atomcenterY(IATOM2)
     Zatomcoor2 = -atomcenterZ(IATOM2)
     DO I2 = 1,N
        DX = Xatomcoor2 + COOR(1,I2+idx) !=-(Xatomcoor2 - COOR(1,I2+idx))
        DY = Yatomcoor2 + COOR(2,I2+idx) !=-(Yatomcoor2 - COOR(2,I2+idx))
        DZ = Zatomcoor2 + COOR(3,I2+idx) !=-(Zatomcoor2 - COOR(3,I2+idx))
        TMP(I2,IATOM2_REL)=SQRT(DX*DX + DY*DY + DZ*DZ)
     ENDDO
  ENDDO

!   evaluate the partition functions

!  This is grid point driven 
!  seems to be the most efficient one, 
!  due to very efficient "screening" \Andreas Krapp
GRIDLOOP: DO I2=1,N
  DO IATOM1_REL=1,nRELEVANT_ATOMS
     P(IATOM1_REL) = D1
     zero(IATOM1_REL) = .FALSE.
  ENDDO

  !we first sort out all values with zero weight
  !first until IATOM (not including IATOM)
  DO IATOM1_REL=1,ATOMINDEX-1
     IATOM1 = RELEVANT_ATOMS(IATOM1_REL)
     dist = inverseDistance12(IATOM1,IATOM)
     MUA = (RADIUS - TMP(I2,IATOM1_REL))*dist
     IF(MUA .GE. SSF_CUTOFF)THEN
        WG(idx+I2)=D0
        CYCLE GRIDLOOP !GO TO I2
     ENDIF
  ENDDO
  !then from IATOM (not including IATOM)
  DO IATOM1_REL=ATOMINDEX+1,nRELEVANT_ATOMS
     IATOM1 = RELEVANT_ATOMS(IATOM1_REL)
     dist = inverseDistance12(IATOM1,IATOM)
     MUA = (RADIUS - TMP(I2,IATOM1_REL))*dist
     IF(MUA .GE. SSF_CUTOFF)THEN
        WG(idx+I2)=D0
        CYCLE GRIDLOOP !GO TO I2
     ENDIF
  ENDDO

  !now for the rest of points
  DO IATOM1_REL=1,nRELEVANT_ATOMS
     IATOM1 = RELEVANT_ATOMS(IATOM1_REL)
     Ria1 = TMP(I2,IATOM1_REL)
     DO IATOM2_REL=1,IATOM1_REL-1
        IF (zero(IATOM1_REL).AND.zero(IATOM2_REL)) CYCLE
        IATOM2 = RELEVANT_ATOMS(IATOM2_REL)
        dist = inverseDistance12(IATOM2,IATOM1)
        MUA = (Ria1 - TMP(I2,IATOM2_REL))*dist
        IF(MUA.LE.-SSF_CUTOFF)THEN
           P(IATOM2_REL) = D0
           zero(IATOM2_REL) = .TRUE.
        ELSEIF(MUA.GE.SSF_CUTOFF)THEN
           P(IATOM1_REL) = D0
           zero(IATOM1_REL) = .TRUE.
        ELSE
           MUA2= MUA*MUA
           GMU=MUA*(fac1+MUA2*(fac2+MUA2*(fac3+fac4*MUA2)));
           P(IATOM1_REL) = P(IATOM1_REL) * (D05-GMU)
           P(IATOM2_REL) = P(IATOM2_REL) * (D05+GMU)
        ENDIF
     ENDDO
  ENDDO
  IF(ABS(P(ATOMINDEX)).LT. 1.0E-30_realk)THEN
     WG(idx+I2) = D0
     W = D0
  ELSE
     !compute weight normalization factors
     V = D0
     DO IATOM2_REL=1,nRELEVANT_ATOMS
        V = V + P(IATOM2_REL)
     ENDDO
     W = P(ATOMINDEX)/V
     WG(I2+idx) = WG(I2+idx)*W
  ENDIF
ENDDO GRIDLOOP

call mem_grid_dealloc(TMP)
call mem_grid_dealloc(P)
call mem_grid_dealloc(RELEVANT_ATOMS)
call mem_grid_dealloc(zero)
END Subroutine SSFPARTITIONING

!based on jcp vol 171, page 2915 (2004)
Subroutine BlockPartitioning(Center,NATOMS,atomcenterX,atomcenterY,&
     & atomcenterZ,COOR2,inverseDistance12,bfac,ihardness,GlobalmaxGridpoints,&
     & npoints,points,atom_idx,WG2,maxNBUFLEN,ITERATIONS,&
     & RADIALPOINTS,NRADPT,weight_thr,LUGRID,NSHELLBLOCKS,SHELLBLOCKS,CELL_SIZE,COOR3,WG3)
  implicit none
  real(realk),intent(in) :: center(3)
  integer,intent(in)     :: NATOMS,ihardness,npoints,GlobalmaxGridpoints
  integer,intent(in)     :: points(npoints),atom_idx(npoints),NRADPT,LUGRID
  integer,intent(in)     :: NSHELLBLOCKS,SHELLBLOCKS(2,NSHELLBLOCKS),maxNBUFLEN
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS)
  real(realk),intent(in) :: inverseDistance12(NATOMS,NATOMS),BFAC(NATOMS,NATOMS)
  real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
  real(realk),intent(in) :: WG2(GlobalmaxGridpoints)
  real(realk),intent(inout) :: COOR3(3,maxNBUFLEN)
  real(realk),intent(inout) :: WG3(maxNBUFLEN)
  integer,intent(inout)   :: ITERATIONS
  REAL(REALK),intent(in) :: RADIALPOINTS(NRADPT,NATOMS),CELL_SIZE,weight_thr
  !
  real(realk),parameter :: D1=1E0_realk,D05=0.5E0_realk,D3=3E0_realk
  real(realk) :: dist,Xatomcoor,Yatomcoor,Zatomcoor,DX,DY,DZ,DIST2,R
  real(realk) :: X,Y,Z,B,FACTOR
  real(realk) :: MU1,MUTMP1(4),MUTMP2(4),MUTMP3(4),MUTMP4(4)
  logical     :: modnpoints
  integer     :: relevant_atoms(NATOMS),IATOM,h,nRELEVANT_ATOMS
  integer     :: IATOM_REL,I,JATOM,JATOM_REL,Ipoint,npoints2
  integer     :: NGRIDPOINTINBATCH,IT,I3
  real(realk),pointer :: RIJ(:,:),TMP2(:,:),MU(:),MU2(:)
  IF(ihardness.EQ.3)THEN !unroll the hardness loop
     !determine relevant atoms
     nRELEVANT_ATOMS=0
     DO IATOM=1,NATOMS
        DX = atomcenterX(IATOM) - CENTER(1)
        DY = atomcenterY(IATOM) - CENTER(2)
        DZ = atomcenterZ(IATOM) - CENTER(3)
        DIST2=DX*DX + DY*DY + DZ*DZ  
        r = RADIALPOINTS(1,IATOM) + CELL_SIZE;
        IF(r*r>dist2)THEN
           nRELEVANT_ATOMS=nRELEVANT_ATOMS+1
           relevant_atoms(nRELEVANT_ATOMS) = IATOM
        ENDIF
     ENDDO
     IF(nRELEVANT_ATOMS.EQ. 0)CALL LSQUIT('CRITICAL ERROR QUIT TK1253',-1)
     IF(nRELEVANT_ATOMS.GT. 1)THEN  
        call mem_grid_alloc(MU,npoints)
!        call mem_grid_alloc(MU2,npoints)     
        call mem_grid_alloc(RIJ,npoints,nRELEVANT_ATOMS)
        call mem_grid_alloc(TMP2,npoints,NATOMS)     
        call build_rij(nRELEVANT_ATOMS,relevant_atoms,NATOMS,atomcenterX,&
             & atomcenterY,atomcenterZ,points,npoints,COOR2,GlobalmaxGridpoints,RIJ)
        call ls_SetToOne(TMP2,NATOMS*npoints)
        npoints2 = (npoints/4)*4          !unroll with 4
        modnpoints = (mod(npoints,4).GT.0)
        DO IATOM_REL=2,nRELEVANT_ATOMS
           IATOM=relevant_atoms(IATOM_REL) 
           DO JATOM_REL = 1,IATOM_REL-1
              JATOM = relevant_atoms(JATOM_REL) 
              dist = inverseDistance12(JATOM,IATOM)
              B = BFAC(JATOM,IATOM)
              DO Ipoint=1,npoints2,4
                 MUTMP1(1) = (RIJ(Ipoint,IATOM_REL)-RIJ(Ipoint,JATOM_REL))*dist
                 MUTMP1(2) = (RIJ(Ipoint+1,IATOM_REL)-RIJ(Ipoint+1,JATOM_REL))*dist
                 MUTMP1(3) = (RIJ(Ipoint+2,IATOM_REL)-RIJ(Ipoint+2,JATOM_REL))*dist
                 MUTMP1(4) = (RIJ(Ipoint+3,IATOM_REL)-RIJ(Ipoint+3,JATOM_REL))*dist

                 MUTMP2(1) = MUTMP1(1)+B*(D1-MUTMP1(1)*MUTMP1(1))
                 MUTMP2(2) = MUTMP1(2)+B*(D1-MUTMP1(2)*MUTMP1(2))
                 MUTMP2(3) = MUTMP1(3)+B*(D1-MUTMP1(3)*MUTMP1(3))
                 MUTMP2(4) = MUTMP1(4)+B*(D1-MUTMP1(4)*MUTMP1(4))

                 MUTMP3(1) = D05*MUTMP2(1)*(D3-MUTMP2(1)*MUTMP2(1)) 
                 MUTMP3(2) = D05*MUTMP2(2)*(D3-MUTMP2(2)*MUTMP2(2)) 
                 MUTMP3(3) = D05*MUTMP2(3)*(D3-MUTMP2(3)*MUTMP2(3)) 
                 MUTMP3(4) = D05*MUTMP2(4)*(D3-MUTMP2(4)*MUTMP2(4)) 

                 MUTMP4(1) = D05*MUTMP3(1)*(D3-MUTMP3(1)*MUTMP3(1))
                 MUTMP4(2) = D05*MUTMP3(2)*(D3-MUTMP3(2)*MUTMP3(2))
                 MUTMP4(3) = D05*MUTMP3(3)*(D3-MUTMP3(3)*MUTMP3(3))
                 MUTMP4(4) = D05*MUTMP3(4)*(D3-MUTMP3(4)*MUTMP3(4))

                 mu(Ipoint)   = D05*MUTMP4(1)*(D3-MUTMP4(1)*MUTMP4(1))
                 mu(Ipoint+1) = D05*MUTMP4(2)*(D3-MUTMP4(2)*MUTMP4(2))
                 mu(Ipoint+2) = D05*MUTMP4(3)*(D3-MUTMP4(3)*MUTMP4(3))
                 mu(Ipoint+3) = D05*MUTMP4(4)*(D3-MUTMP4(4)*MUTMP4(4))
              ENDDO
              IF(modnpoints)THEN
                 DO Ipoint=npoints2+1,npoints
                    MUTMP1(1) = (RIJ(Ipoint,IATOM_REL)-RIJ(Ipoint,JATOM_REL))*dist
                    MUTMP2(1) = MUTMP1(1)+B*(D1-MUTMP1(1)*MUTMP1(1))
                    MUTMP3(1) = D05*MUTMP2(1)*(D3-MUTMP2(1)*MUTMP2(1)) 
                    MUTMP4(1) = D05*MUTMP3(1)*(D3-MUTMP3(1)*MUTMP3(1))
                    mu(Ipoint)   = D05*MUTMP4(1)*(D3-MUTMP4(1)*MUTMP4(1))
                 ENDDO
              ENDIF
              DO Ipoint=1,npoints
                 TMP2(Ipoint,IATOM) = TMP2(Ipoint,IATOM)*D05*(D1-mu(Ipoint)) 
                 TMP2(Ipoint,JATOM) = TMP2(Ipoint,JATOM)*D05*(D1+mu(Ipoint)) 
              ENDDO
           ENDDO
        ENDDO
        call mem_grid_dealloc(RIJ)
        IATOM = relevant_atoms(1) 
        DO Ipoint=1,npoints
           MU(Ipoint)=TMP2(Ipoint,IATOM)
        ENDDO
        DO IATOM_REL=2,nRELEVANT_ATOMS
           IATOM = relevant_atoms(IATOM_REL) 
           DO Ipoint=1,npoints
              MU(Ipoint)=MU(Ipoint) + TMP2(Ipoint,IATOM)
           ENDDO
        ENDDO
        if(npoints.GT.maxNBUFLEN)THEN
           CALL WRITE_COORD4(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
                & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,atom_idx,weight_thr,&
                & maxNBUFLEN,LUGRID,COOR3,WG3)
        ELSE
           CALL WRITE_COORD3(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
                & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,atom_idx,weight_thr,&
                & maxNBUFLEN,LUGRID,COOR3,WG3)
        ENDIF
        call mem_grid_dealloc(TMP2) 
        call mem_grid_dealloc(MU)
!        call mem_grid_dealloc(MU2)         
     ELSE
        if(npoints.GT.maxNBUFLEN)THEN
           CALL WRITE_COORD2(npoints,ITERATIONS,points(1:npoints),&
                & COOR2,WG2,GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,&
                & maxNBUFLEN,LUGRID,COOR3,WG3)
        else
           CALL WRITE_COORD1(npoints,ITERATIONS,points(1:npoints),&
                & COOR2,WG2,GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,&
                & maxNBUFLEN,LUGRID,COOR3,WG3)
        endif
     ENDIF

  ELSE
     !determine relevant atoms
     nRELEVANT_ATOMS=0
     DO IATOM=1,NATOMS
        DX = atomcenterX(IATOM) - CENTER(1)
        DY = atomcenterY(IATOM) - CENTER(2)
        DZ = atomcenterZ(IATOM) - CENTER(3)
        DIST2=DX*DX + DY*DY + DZ*DZ  
        r = RADIALPOINTS(1,IATOM) + CELL_SIZE;
        IF(r*r>dist2)THEN
           nRELEVANT_ATOMS=nRELEVANT_ATOMS+1
           relevant_atoms(nRELEVANT_ATOMS) = IATOM
        ENDIF
     ENDDO
     IF(nRELEVANT_ATOMS.EQ. 0)CALL LSQUIT('CRITICAL ERROR QUIT TK1253',-1)
     IF(nRELEVANT_ATOMS.GT. 1)THEN  
        call mem_grid_alloc(MU,npoints)
        call mem_grid_alloc(MU2,npoints)     
        call mem_grid_alloc(RIJ,npoints,nRELEVANT_ATOMS)
        call mem_grid_alloc(TMP2,npoints,NATOMS)     
        call build_rij(nRELEVANT_ATOMS,relevant_atoms,NATOMS,atomcenterX,&
             & atomcenterY,atomcenterZ,points,npoints,COOR2,GlobalmaxGridpoints,RIJ)
        call ls_SetToOne(TMP2,NATOMS*npoints)
        DO IATOM_REL=2,nRELEVANT_ATOMS
           IATOM=relevant_atoms(IATOM_REL) 
           DO JATOM_REL = 1,IATOM_REL-1
              JATOM = relevant_atoms(JATOM_REL) 
              dist = inverseDistance12(JATOM,IATOM)
              B = BFAC(JATOM,IATOM)
              DO Ipoint=1,npoints
                 MU1=(RIJ(Ipoint,IATOM_REL)-RIJ(Ipoint,JATOM_REL))*dist
                 MU1=MU1+B*(D1-MU1*MU1)
                 MU(Ipoint) = MU1
                 MU2(Ipoint) = MU1*MU1
              ENDDO
              DO h=1,ihardness-1              
                 DO Ipoint=1,npoints
                    mu(Ipoint) = D05*mu(Ipoint)*(D3-mu2(Ipoint)) 
                    mu2(Ipoint) = mu(Ipoint)*mu(Ipoint)
                 ENDDO
              ENDDO
              !        h=ihardness
              DO Ipoint=1,npoints
                 mu(Ipoint) = D05*mu(Ipoint)*(D3-mu2(Ipoint))    
              ENDDO
              DO Ipoint=1,npoints
                 TMP2(Ipoint,IATOM) = TMP2(Ipoint,IATOM)*D05*(D1-mu(Ipoint)) 
                 TMP2(Ipoint,JATOM) = TMP2(Ipoint,JATOM)*D05*(D1+mu(Ipoint)) 
              ENDDO
           ENDDO
        ENDDO
        call mem_grid_dealloc(RIJ)
        IATOM = relevant_atoms(1) 
        DO Ipoint=1,npoints
           MU(Ipoint)=TMP2(Ipoint,IATOM)
        ENDDO
        DO IATOM_REL=2,nRELEVANT_ATOMS
           IATOM = relevant_atoms(IATOM_REL) 
           DO Ipoint=1,npoints
              MU(Ipoint)=MU(Ipoint) + TMP2(Ipoint,IATOM)
           ENDDO
        ENDDO
        if(npoints.GT.maxNBUFLEN)THEN
           CALL WRITE_COORD4(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
                & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,atom_idx,weight_thr,&
                & maxNBUFLEN,LUGRID,COOR3,WG3)
        ELSE
           CALL WRITE_COORD3(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
                & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,atom_idx,weight_thr,&
                & maxNBUFLEN,LUGRID,COOR3,WG3)
        ENDIF
        call mem_grid_dealloc(TMP2) 
        call mem_grid_dealloc(MU)
        call mem_grid_dealloc(MU2)         
     ELSE
        if(npoints.GT.maxNBUFLEN)THEN
           CALL WRITE_COORD2(npoints,ITERATIONS,points(1:npoints),&
                & COOR2,WG2,GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,&
                & maxNBUFLEN,LUGRID,COOR3,WG3)
        else
           CALL WRITE_COORD1(npoints,ITERATIONS,points(1:npoints),&
                & COOR2,WG2,GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,&
                & maxNBUFLEN,LUGRID,COOR3,WG3)
        endif
     ENDIF
  ENDIF
END Subroutine BLOCKPARTITIONING

! SSF scheme for grid weight evaluation 
! combined with a blockwise handling of grid points 
! the later be described in JCP, 2004, 171, 2915.
!
! the later allows to reduce the list of relevant atom pairs 
! considerably 
! \Andreas Krapp, 07/2010
!
Subroutine BlockSSFPartitioning(Center,NATOMS,atomcenterX,atomcenterY,&
     & atomcenterZ,COOR2,inverseDistance12,ihardness,GlobalmaxGridpoints,&
     & npoints,points,atom_idx,WG2,maxNBUFLEN,ITERATIONS,&
     & RADIALPOINTS,NRADPT,weight_thr,LUGRID,NSHELLBLOCKS,SHELLBLOCKS,CELL_SIZE,skip,COOR3,WG3)
  implicit none
  real(realk),intent(in) :: center(3)
  integer,intent(in)     :: NATOMS,ihardness,npoints,GlobalmaxGridpoints
  integer,intent(in)     :: points(npoints),atom_idx(npoints),NRADPT,LUGRID
  integer,intent(in)     :: NSHELLBLOCKS,SHELLBLOCKS(2,NSHELLBLOCKS),maxNBUFLEN
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS)
  real(realk),intent(in) :: inverseDistance12(NATOMS,NATOMS)
  real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
  real(realk),intent(in) :: WG2(GlobalmaxGridpoints)
  real(realk),intent(inout) :: COOR3(3,maxNBUFLEN)
  real(realk),intent(inout) :: WG3(maxNBUFLEN)
  integer,intent(inout)   :: ITERATIONS
  REAL(REALK),intent(in) :: RADIALPOINTS(NRADPT,NATOMS),CELL_SIZE,weight_thr
  logical,intent(in) :: skip(GlobalmaxGridpoints)
  !
  real(realk),parameter :: D1=1E0_realk,D05=0.5E0_realk,D3=3E0_realk,D0=0E0_realk
  real(realk) :: dist,Xatomcoor,Yatomcoor,Zatomcoor,DX,DY,DZ,DIST2,R
  real(realk) :: X,Y,Z,FACTOR,MU1,MU2,GMU,RIJA
  integer     :: IATOM,h,nRELEVANT_ATOMS,IATOM_REL,I,JATOM,JATOM_REL,Ipoint
  integer     :: NGRIDPOINTINBATCH,IT,I3,I2,nWeight1points,nrelevantpoints
  logical     :: addpoint
  integer,pointer     :: relevantpoints(:),relevant_atoms(:),Weight1points(:),TMPINDEX(:)
  real(realk),pointer :: RIJ(:,:),TMP2(:,:),MU(:)
  logical,pointer     :: zero(:,:)
  
  call mem_grid_alloc(relevant_atoms,NATOMS)
  call mem_grid_alloc(TMPINDEX,NATOMS)
  !determine relevant atoms
  nRELEVANT_ATOMS=0
  DO IATOM=1,NATOMS
     DX = atomcenterX(IATOM) - CENTER(1)
     DY = atomcenterY(IATOM) - CENTER(2)
     DZ = atomcenterZ(IATOM) - CENTER(3)
     DIST2=DX*DX + DY*DY + DZ*DZ  
     r = RADIALPOINTS(1,IATOM) + CELL_SIZE
     IF(r*r>dist2)THEN
        nRELEVANT_ATOMS=nRELEVANT_ATOMS+1
        relevant_atoms(nRELEVANT_ATOMS) = IATOM
        TMPINDEX(IATOM) = nRELEVANT_ATOMS
     ENDIF
  ENDDO
  IF(nRELEVANT_ATOMS.EQ. 0)CALL LSQUIT('CRITICAL ERROR QUIT TK1253',-1)
  IF(nRELEVANT_ATOMS.GT. 1)THEN       
     call mem_grid_alloc(relevantpoints,npoints)
     call mem_grid_alloc(Weight1points,npoints)
     call mem_grid_alloc(zero,NATOMS,npoints)
     call mem_grid_alloc(RIJ,nRELEVANT_ATOMS,npoints)
     call build_rij_inv(nRELEVANT_ATOMS,relevant_atoms,NATOMS,atomcenterX,&
          & atomcenterY,atomcenterZ,points,npoints,COOR2,GlobalmaxGridpoints,RIJ)
     I=0
     I2=0
     DO Ipoint=1,npoints
        IF(skip(points(Ipoint)))THEN
           I=I+1
           Weight1points(I)=Ipoint
        ELSE
           !Collect points with nonzero weights
           ADDPOINT=.TRUE.
           IATOM = atom_idx(Ipoint)
           IATOM_REL = TMPINDEX(IATOM)
           RIJA = RIJ(IATOM_REL,Ipoint)
           ATOMLOOP: DO JATOM_REL=1,nRELEVANT_ATOMS
              IF(JATOM_REL.EQ.IATOM_REL)CYCLE
              JATOM=relevant_atoms(JATOM_REL) 
              dist= inverseDistance12(JATOM,IATOM)
              MU1 = (RIJA-RIJ(JATOM_REL,Ipoint))*dist
              IF(MU1.GE.SSF_CUTOFF)THEN
                 ADDPOINT=.FALSE.
                 EXIT ATOMLOOP
              ENDIF
           ENDDO ATOMLOOP
           IF(ADDPOINT)THEN
              I2 = I2 + 1
              relevantpoints(I2) = Ipoint
           ENDIF
        ENDIF
     ENDDO
     nWeight1points = I
     nrelevantpoints=I2
     !build only for relevant gridpoints
     call mem_grid_dealloc(RIJ)
     call mem_grid_alloc(RIJ,nRELEVANT_ATOMS,nrelevantpoints)
     call build_rij2_inv(nRELEVANT_ATOMS,relevant_atoms,NATOMS,atomcenterX,&
          & atomcenterY,atomcenterZ,points,npoints,COOR2,GlobalmaxGridpoints,&
          & relevantpoints(1:nrelevantpoints),nrelevantpoints,RIJ)
     ! Keep track of zero elements   
     DO iatom=1,nRELEVANT_ATOMS
       DO ipoint=1,nrelevantpoints
         zero(iatom,ipoint) = .FALSE.
       ENDDO
     ENDDO
     !compute weights
     call mem_grid_alloc(TMP2,npoints,NATOMS)
     call ls_SetToOne(TMP2,npoints*NATOMS)
     DO Ipoint=1,nrelevantpoints!npoints
       DO IATOM_REL=1,nRELEVANT_ATOMS
         IF (zero(iatom_rel,ipoint)) CYCLE
         RIJA = RIJ(IATOM_REL,Ipoint)
         IATOM=relevant_atoms(IATOM_REL) 
         DO JATOM_REL = 1,nRELEVANT_ATOMS
            IF (IATOM_REL.EQ.JATOM_REL) CYCLE
 
            JATOM = relevant_atoms(JATOM_REL) 
            dist = inverseDistance12(JATOM,IATOM)
            MU1=(RIJA-RIJ(JATOM_REL,Ipoint))*dist
            IF(MU1.LE.-SSF_CUTOFF)THEN
               TMP2(Ipoint,JATOM) = D0
               zero(jatom_rel,ipoint) = .TRUE.
            ELSEIF(MU1.GE.SSF_CUTOFF)THEN
               TMP2(Ipoint,IATOM) = D0
               zero(iatom_rel,ipoint) = .TRUE.
               EXIT
            ELSE
               MU2=MU1*MU1
               GMU=MU1*(fac1+MU2*(fac2+MU2*(fac3+fac4*MU2)))
               TMP2(Ipoint,IATOM) = TMP2(Ipoint,IATOM) * (D05-GMU)
            ENDIF
         ENDDO
       ENDDO
     ENDDO
     call mem_grid_dealloc(RIJ)
     call mem_grid_dealloc(zero)
     call mem_grid_alloc(MU,npoints)
     IATOM_REL=1
     IATOM=relevant_atoms(IATOM_REL) 
     DO Ipoint=1,nrelevantpoints!npoints
        MU(Ipoint)=TMP2(Ipoint,IATOM)
     ENDDO
     DO IATOM_REL=2,nRELEVANT_ATOMS
        IATOM=relevant_atoms(IATOM_REL) 
        DO Ipoint=1,nrelevantpoints!npoints
           MU(Ipoint)=MU(Ipoint) + TMP2(Ipoint,IATOM)
        ENDDO
     ENDDO
     if(npoints.GT.maxNBUFLEN)THEN
        CALL WRITE_COORD6(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
             & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,atom_idx,weight_thr,&
             & maxNBUFLEN,LUGRID,relevantpoints(1:nrelevantpoints),&
             & nrelevantpoints,nWeight1points,Weight1points(1:nWeight1points),COOR3,WG3)
     ELSE
        CALL WRITE_COORD5(npoints,ITERATIONS,points,COOR2,WG2,MU,TMP2,NATOMS,&
             & GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,atom_idx,weight_thr,&
             & maxNBUFLEN,LUGRID,relevantpoints(1:nrelevantpoints),&
             & nrelevantpoints,nWeight1points,Weight1points(1:nWeight1points),COOR3,WG3)
     ENDIF
     call mem_grid_dealloc(relevantpoints)
     call mem_grid_dealloc(MU)
     call mem_grid_dealloc(TMP2)
     call mem_grid_dealloc(Weight1points)
  ELSE
     if(npoints.GT.maxNBUFLEN)THEN
        CALL WRITE_COORD2(npoints,ITERATIONS,points(1:npoints),&
             & COOR2,WG2,GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,maxNBUFLEN,LUGRID,COOR3,WG3)
     else
        CALL WRITE_COORD1(npoints,ITERATIONS,points(1:npoints),&
             & COOR2,WG2,GlobalmaxGridpoints,NSHELLBLOCKS,SHELLBLOCKS,maxNBUFLEN,LUGRID,COOR3,WG3)
     endif
  ENDIF
  call mem_grid_dealloc(relevant_atoms)
  call mem_grid_dealloc(TMPINDEX)
END Subroutine BLOCKSSFPARTITIONING

subroutine build_rij(nRELEVANT_ATOMS,relevant_atoms,NATOMS,atomcenterX,&
     & atomcenterY,atomcenterZ,points,npoints,COOR2,GlobalmaxGridpoints,RIJ)
  IMPLICIT NONE
  INTEGER,intent(in) :: nRELEVANT_ATOMS,npoints,NATOMS,GlobalmaxGridpoints
  INTEGER,intent(in) :: relevant_atoms(NATOMS),points(npoints)
  real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS)
  real(realk),intent(in) :: atomcenterZ(NATOMS)
  real(realk),intent(inout) :: RIJ(npoints,nRELEVANT_ATOMS)
!
  INTEGER :: IATOM_REL,IATOm,Ipoint,I
  REAL(REALK) :: X,Y,Z,DX,DY,DZ

  DO IATOM_REL=1,nRELEVANT_ATOMS
     IATOM=relevant_atoms(IATOM_REL) 
     X=-atomcenterX(IATOM)
     Y=-atomcenterY(IATOM)
     Z=-atomcenterZ(IATOM)
     DO Ipoint=1,npoints
        I = points(ipoint)
        DX=COOR2(1,I)+X
        DY=COOR2(2,I)+Y
        DZ=COOR2(3,I)+Z
        RIJ(Ipoint,IATOM_REL) = SQRT(DX*DX+DY*DY+DZ*DZ)        
     ENDDO
  ENDDO
END subroutine BUILD_RIJ

subroutine build_rij_inv(nRELEVANT_ATOMS,relevant_atoms,NATOMS,atomcenterX,&
     & atomcenterY,atomcenterZ,points,npoints,COOR2,GlobalmaxGridpoints,RIJ)
  IMPLICIT NONE
  INTEGER,intent(in) :: nRELEVANT_ATOMS,npoints,NATOMS,GlobalmaxGridpoints
  INTEGER,intent(in) :: relevant_atoms(NATOMS),points(npoints)
  real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS)
  real(realk),intent(in) :: atomcenterZ(NATOMS)
  real(realk),intent(inout) :: RIJ(nRELEVANT_ATOMS,npoints)
!
  INTEGER :: IATOM_REL,IATOm,Ipoint,I
  REAL(REALK) :: X,Y,Z,DX,DY,DZ

  DO Ipoint=1,npoints
     I = points(ipoint)
     X=-COOR2(1,I)
     Y=-COOR2(2,I)
     Z=-COOR2(3,I)
     DO IATOM_REL=1,nRELEVANT_ATOMS
        IATOM=relevant_atoms(IATOM_REL) 
        DX=atomcenterX(IATOM)+X
        DY=atomcenterY(IATOM)+Y
        DZ=atomcenterZ(IATOM)+Z
        RIJ(IATOM_REL,Ipoint) = SQRT(DX*DX+DY*DY+DZ*DZ)        
     ENDDO
  ENDDO
END subroutine BUILD_RIJ_inv

subroutine build_rij2(nRELEVANT_ATOMS,relevant_atoms,NATOMS,atomcenterX,&
     & atomcenterY,atomcenterZ,points,npoints,COOR2,GlobalmaxGridpoints,&
     & relevantpoints,nrelevantpoints,RIJ)
  IMPLICIT NONE
  INTEGER,intent(in) :: nRELEVANT_ATOMS,npoints,NATOMS,GlobalmaxGridpoints
  INTEGER,intent(in) :: relevant_atoms(NATOMS),points(npoints)
  INTEGER,intent(in) :: nrelevantpoints
  real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS)
  real(realk),intent(in) :: atomcenterZ(NATOMS)
  real(realk),intent(inout) :: RIJ(nrelevantpoints,nRELEVANT_ATOMS)
  integer,intent(in) :: relevantpoints(nrelevantpoints)
!
  INTEGER :: IATOM_REL,IATOm,Ipoint,I,Irelpoint
  REAL(REALK) :: X,Y,Z,DX,DY,DZ

  DO IATOM_REL=1,nRELEVANT_ATOMS
     IATOM=relevant_atoms(IATOM_REL) 
     X=-atomcenterX(IATOM)
     Y=-atomcenterY(IATOM)
     Z=-atomcenterZ(IATOM)
     DO Ipoint=1,nrelevantpoints
        Irelpoint = relevantpoints(Ipoint)
        I = points(Irelpoint)
        DX=COOR2(1,I)+X
        DY=COOR2(2,I)+Y
        DZ=COOR2(3,I)+Z
        RIJ(Ipoint,IATOM_REL) = SQRT(DX*DX+DY*DY+DZ*DZ)        
     ENDDO
  ENDDO
END subroutine BUILD_RIJ2

subroutine build_rij2_inv(nRELEVANT_ATOMS,relevant_atoms,NATOMS,atomcenterX,&
     & atomcenterY,atomcenterZ,points,npoints,COOR2,GlobalmaxGridpoints,&
     & relevantpoints,nrelevantpoints,RIJ)
  IMPLICIT NONE
  INTEGER,intent(in) :: nRELEVANT_ATOMS,npoints,NATOMS,GlobalmaxGridpoints
  INTEGER,intent(in) :: relevant_atoms(NATOMS),points(npoints)
  INTEGER,intent(in) :: nrelevantpoints
  real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS)
  real(realk),intent(in) :: atomcenterZ(NATOMS)
  real(realk),intent(inout) :: RIJ(nRELEVANT_ATOMS,nrelevantpoints)
  integer,intent(in) :: relevantpoints(nrelevantpoints)
!
  INTEGER :: IATOM_REL,IATOm,Ipoint,I,Irelpoint
  REAL(REALK) :: X,Y,Z,DX,DY,DZ

  DO Ipoint=1,nrelevantpoints
     Irelpoint = relevantpoints(Ipoint)
     I = points(Irelpoint)
     X=-COOR2(1,I)
     Y=-COOR2(2,I)
     Z=-COOR2(3,I)
     DO IATOM_REL=1,nRELEVANT_ATOMS
        IATOM=relevant_atoms(IATOM_REL) 
        DX=atomcenterX(IATOM)+X
        DY=atomcenterY(IATOM)+Y
        DZ=atomcenterZ(IATOM)+Z
        RIJ(IATOM_REL,Ipoint) = SQRT(DX*DX+DY*DY+DZ*DZ)        
     ENDDO
  ENDDO
END subroutine BUILD_RIJ2_INV

Subroutine DetermineGridDistance(TMP,atomcenterX,atomcenterY,atomcenterZ,&
     & NATOMS,N,COOR,MaxAtomicGridpoints,idx)
  implicit none
  Integer,intent(in) :: NATOMS,N,idx,MaxAtomicGridpoints
  real(realk),intent(in) :: atomcenterX(NATOMS),atomcenterY(NATOMS)
  real(realk),intent(in) :: atomcenterZ(NATOMS),COOR(3,MaxAtomicGridpoints)
  real(realk),intent(inout) :: TMP(N,NATOMS)
  !
  INTEGER :: IATOM2,I2
  real(realk) :: Xatomcoor2,Yatomcoor2,Zatomcoor2,DX,DY,DZ
  
  DO IATOM2=1,NATOMS
     Xatomcoor2 = -atomcenterX(IATOM2)
     Yatomcoor2 = -atomcenterY(IATOM2)
     Zatomcoor2 = -atomcenterZ(IATOM2)
     DO I2 = 1,N
        DX = Xatomcoor2 + COOR(1,I2+idx) !=-(Xatomcoor2 - COOR(1,I2+idx))
        DY = Yatomcoor2 + COOR(2,I2+idx) !=-(Yatomcoor2 - COOR(2,I2+idx))
        DZ = Zatomcoor2 + COOR(3,I2+idx) !=-(Zatomcoor2 - COOR(3,I2+idx))
        TMP(I2,IATOM2)=SQRT(DX*DX + DY*DY + DZ*DZ)
     ENDDO
  ENDDO
end Subroutine DetermineGridDistance

SUBROUTINE GRID_GETBLOCKS(RSHEL,MAXNSHELL,SHELL2ATOM,NBAST,NATOMS&
     & ,X,Y,Z,NSHELLBLOCKS,SHELLBLOCKS,NactBast,NSTART,CENTER,CELL_SIZE,LUPRI)
use precision
IMPLICIT NONE
  !     get blocks of active SHELLS in box 
  !
  !     RSHEL2 - precomputed shell extents (squared).
  !     NSHELLBLOCKS (output) - number of active blocks
  !     SHELLBLOCKS (output) - pairs of (startindex, stopindex)
INTEGER,intent(in)    :: NBAST,MAXNSHELL,NATOMS,LUPRI
INTEGER,intent(in)    :: SHELL2ATOM(MAXNSHELL),NSTART(MAXNSHELL)
REAL(REALK),intent(in):: X(NATOMS),Y(NATOMS),Z(NATOMS),RSHEL(MAXNSHELL)
REAL(REALK),intent(in):: CELL_SIZE,CENTER(3)
INTEGER,intent(inout) :: NactBast
INTEGER,intent(inout) :: SHELLBLOCKS(2,MAXNSHELL),NSHELLBLOCKS
!
INTEGER               :: ISHLEN,IPREV,ICENT,ISHELA,I,J
REAL(REALK)           :: DST,PX,PY,PZ,RSHEL1,LengthX,LengthY,LengthZ
INTEGER               :: ORBBLOCKS(2,MAXNSHELL),IBL,IORB
LOGICAL               :: XVERIFY,YVERIFY,ZVERIFY
REAL(REALK),parameter :: D05=0.5E0_realk
NSHELLBLOCKS = 0
IPREV = -1111
DO ISHELA=1,MAXNSHELL
   ICENT = SHELL2ATOM(ISHELA)
   PX = center(1)-X(ICENT)
   PY = center(2)-Y(ICENT)
   PZ = center(3)-Z(ICENT)
   DST = SQRT(PX*PX + PY*PY + PZ*PZ)
   IF(DST.LE.RSHEL(ISHELA)+CELL_SIZE) THEN
      !       accepted...
      IF(ISHELA.NE.IPREV+1) THEN
         NSHELLBLOCKS = NSHELLBLOCKS + 1
         SHELLBLOCKS(1,NSHELLBLOCKS) = ISHELA
      END IF
      IPREV = ISHELA
      SHELLBLOCKS(2,NSHELLBLOCKS) = ISHELA
!      write(lupri,'(A,2I3,A,3F6.2,A,F6.2,A)')'shell',ISHELA,ICENT,' at ',&
!      &X(ICENT),Y(ICENT),Z(ICENT),' rad. ',&
!      & SQRT(RSHEL(ISHELA)),' accepted'
   ELSE
!      write(lupri,'(A,2I3,A,3F6.2,A,F6.2,A)')'shell',ISHELA,ICENT,' at ',&
!           &X(ICENT),Y(ICENT),Z(ICENT),' rad. ',&
!           & SQRT(RSHEL(ISHELA)),' rejected'
   END IF
END DO !shellloop
!write(lupri,'(A,3F15.5,A,I3)') 'cell at:',CENTER(1),CENTER(2),CENTER(3),' blocks:',NSHELLBLOCKS
!write(lupri,'(8(A,I3,I3,A)/,(8(A,I3,I3,A)))') (('[',SHELLBLOCKS(1,J),SHELLBLOCKS(2,J),']'),J=1,NSHELLBLOCKS)
IF(NSHELLBLOCKS.GT.0)THEN
   DO I = 1,NSHELLBLOCKS
      ORBBLOCKS(1,I) = NSTART(SHELLBLOCKS(1,I))+1
   ENDDO
   DO I = 1,NSHELLBLOCKS-1
      ORBBLOCKS(2,I) = NSTART(SHELLBLOCKS(2,I)+1)
   ENDDO
   I = NSHELLBLOCKS
   IF(SHELLBLOCKS(2,I) .LT. MAXNSHELL)THEN
      ORBBLOCKS(2,I) = NSTART(SHELLBLOCKS(2,I)+1)
   ELSE
      ORBBLOCKS(2,I) = NBAST
   ENDIF
   NactBAST = 0
   DO IBL = 1, NSHELLBLOCKS
      NactBAST = NactBAST + ORBBLOCKS(2,IBL) - ORBBLOCKS(1,IBL) + 1
      !   DO IORB = ORBBLOCKS(1,IBL),ORBBLOCKS(2,IBL)
      !      NactBAST = NactBAST + 1
      !   ENDDO
   ENDDO
ELSE
   ORBBLOCKS = 0
   NactBAST = 0
ENDIF

END SUBROUTINE GRID_GETBLOCKS

subroutine SHELLLIST_ORBBLOCK(LUPRI,SHELLLIST,MAXNSHELL,SHELLBLOCKS)
implicit none
integer,intent(in)    :: LUPRI,MAXNSHELL
logical,intent(in)    :: SHELLLIST(MAXNSHELL)
integer,intent(inout) :: SHELLBLOCKS(2,MAXNSHELL)
!
INTEGER :: I,N
LOGICAL :: INSIDEBLOCK

INSIDEBLOCK=.FALSE.
I=0
N=0
blockloop: DO 
   IF(INSIDEBLOCK)THEN
      insideloop: DO
         I=I+1
         IF(.NOT.SHELLLIST(I))THEN
            INSIDEBLOCK=.FALSE.
            SHELLBLOCKS(2,N) = I-1
            exit insideloop
         ENDIF
         IF(I.EQ.MAXNSHELL)THEN
            INSIDEBLOCK=.FALSE.
            SHELLBLOCKS(2,N) = I
            exit blockloop
         ENDIF
      ENDDO insideloop
   ELSE
      outsideloop: DO
         I=I+1
         IF(SHELLLIST(I))THEN
            N=N+1
            INSIDEBLOCK=.TRUE.
            SHELLBLOCKS(1,N) = I
            exit outsideloop
         ENDIF
         IF(I.EQ.MAXNSHELL)THEN
            exit blockloop
         ENDIF
      ENDDO outsideloop
   ENDIF
ENDDO blockloop

DO I=1,N
   WRITE(6,*)'ACTIVE SHELL BLOCKS',SHELLBLOCKS(1,I),SHELLBLOCKS(2,I)
ENDDO

end subroutine SHELLLIST_ORBBLOCK

SUBROUTINE GRID_BUILDSHELLLIST(SHELLLIST,MAXNSHELL,NSHELLBLOCKS,SHELLBLOCKS,LUPRI)
implicit none
INTEGER,intent(in) :: NSHELLBLOCKS,MAXNSHELL,LUPRI
INTEGER,intent(in) :: SHELLBLOCKS(2,MAXNSHELL)
LOGICAL,intent(inout) :: SHELLLIST(MAXNSHELL)
!
INTEGER :: I,IBL
DO I=1,NSHELLBLOCKS
   DO IBL = SHELLBLOCKS(1,I),SHELLBLOCKS(2,I)
      SHELLLIST(IBL) = .TRUE.
   ENDDO
ENDDO

end SUBROUTINE GRID_BUILDSHELLLIST

SUBROUTINE Determine_BraggFac(NATOMS,BFAC,CHARGE)
implicit none
Integer,intent(in) :: NATOMS
Integer,intent(in) :: CHARGE(NATOMS)
REAL(REALK),intent(inout) :: BFAC(NATOMS,NATOMS) 
!
integer :: IATOM1,IATOM2,ICHARGE1,ICHARGE2
real(realk) :: B1,B2,CHI,temp
real(realk),parameter :: D1=1E0_realk,D05=0.5E0_realk

DO IATOM1=1,NATOMS
   ICHARGE1 = CHARGE(IATOM1)
   B1 = bragg_radii(ICHARGE1)
   DO IATOM2=1,NATOMS
      ICHARGE2 = CHARGE(IATOM2)
      B2 = bragg_radii(ICHARGE2)
      CHI = B1/B2
      temp = (CHI-D1)/(CHI+D1)  !Eq. A6 in jcp vol 88, page 2547
      temp = temp/(temp*temp-1) !Eq. A5 in jcp vol 88, page 2547
      IF(temp.GT.D05)THEN
         temp=D05
      ELSEIF(temp.LT.-D05)THEN
         temp = -D05
      ENDIF
      BFAC(IATOM2,IATOM1) = temp
   ENDDO
ENDDO

END SUBROUTINE DETERMINE_BRAGGFAC

SUBROUTINE SET_ANGULAR(GRIDANG,radint,angint,charge,natoms,NRADPT,nRadialPoints,RADIALPOINTS,RADIALWEIGHT,&
     &                 iprune,angularpoints,ZdependentMaxAng)
implicit none
integer,intent(in)  :: angint,natoms,iprune,NRADPT,angularpoints
!> charge for each atom used in grid-generation
INTEGER,intent(in)  :: CHARGE(natoms),nRadialPoints(natoms)
LOGICAL,intent(in)  :: ZdependentMaxAng
real(realk),intent(in) :: radint
integer,intent(inout) :: GRIDANG(NRADPT,NATOMS)
REAL(REALK),intent(in) :: RADIALPOINTS(NRADPT,NATOMS),RADIALWEIGHT(NRADPT,NATOMS)
!
REAL(REALK),parameter :: D1=1E0_realk,D2=2E0_realk
REAL(REALK) :: INVERSEBOHR,RBRAGG
INTEGER :: atom_maxang,iatom,icharge,I,iang,n_points_optimal,iang1
INVERSEBOHR = D1/(D2*bohr_to_angstrom)
DO IATOM=1,NATOMS
   atom_maxang = angint
   icharge = charge(IATOM)
   rbragg = bragg_radii(icharge)*INVERSEBOHR
   IF(ZdependentMaxAng)THEN 
      !this is the maxang scheme from turbomole
      IF(icharge.LE. 2)THEN
         atom_maxang = atom_maxang-6
         IF(angint.EQ. 47) atom_maxang = atom_maxang-6
      ENDIF
   ENDIF
!   thís pruning, follows Murray, Handy, Laming (MolPhys.1993,78,997)
   IF(IPRUNE.EQ. 1)THEN
      DO I=1,nRadialPoints(IATOM)
         iang = leb_get_from_order(atom_maxang)
         IF(RADIALPOINTS(I,IATOM).LT.RBRAGG)THEN
            n_points_optimal=INT(RADIALPOINTS(I,IATOM)*angularpoints/RBRAGG)
            iang1 = leb_get_from_point(n_points_optimal)
            IF(iang1 .LT. iang)THEN
               iang = iang1
            ENDIF
         ENDIF
         GRIDANG(I,IATOM)=iang
      ENDDO      
   ELSE
      DO I=1,nRadialPoints(IATOM)
         GRIDANG(I,IATOM) = leb_get_from_order(atom_maxang)
      ENDDO
   ENDIF
ENDDO
END SUBROUTINE SET_ANGULAR
   
INTEGER FUNCTION leb_get_from_order(angint)
integer,intent(in) :: angint
! I am unsure about leb_gen_poly_order(1) and leb_gen_poly_order(2)  and leb_gen_poly_order(3) 
integer,parameter :: leb_gen_poly_order(18) = (/  5, 9,11,15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53,  59,  64 /)
integer,parameter :: leb_gen_point(18)      = (/  14,38,50,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454 /)
integer :: I
leb_get_from_order = 0
DO I=1,18 !size of leb_gen_poly_order
   if(leb_gen_poly_order(I).GE.angint)THEN
      leb_get_from_order = I
      RETURN
   ENDIF
ENDDO
CALL LSQUIT('leb_get_from_order error ',-1)
END FUNCTION LEB_GET_FROM_ORDER

INTEGER FUNCTION leb_get_from_point(point)
integer,intent(in) :: point
! I am unsure about leb_gen_poly_order(1) and leb_gen_poly_order(2)  and leb_gen_poly_order(3) 
integer,parameter :: leb_gen_poly_order(18) = (/  5, 9,11,15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53,  59,  64 /)
integer,parameter :: leb_gen_point(18)      = (/  14,38,50,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454 /)
integer :: I
leb_get_from_point = 1
DO I=18,1,-1 !size of leb_gen_point
   if(point.GE.leb_gen_point(I))THEN
      leb_get_from_point = I
      RETURN
   ENDIF
ENDDO
END FUNCTION LEB_GET_FROM_POINT

INTEGER FUNCTION  leb_gen_points(I)
integer,intent(in) :: I
integer,parameter :: leb_gen_point(18)      = (/  14,38,50,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454 /)

leb_gen_points = leb_gen_point(I)

END FUNCTION LEB_GEN_POINTS

SUBROUTINE SET_RADIAL(RADIALPOINTS,RADIALWEIGHT,nRadialPoints,NRADPT,SHELL2ATOM,SHELLANGMOM,SHELLNPRIM,MAXANGMOM,NATOMS,MAXNSHELL,&
        & MXPRIM,PRIEXP,PRIEXPSTART,IPRINT,LUPRI,RADINT,Charge,RADIALGRID,MIN_RAD_PT)
implicit none
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)  :: IPRINT
!> maximum number of shells
INTEGER,intent(in)  :: MAXNSHELL
!> number of atoms
INTEGER,intent(in)  :: natoms
!> the maximum angular momentum
INTEGER,intent(in)  :: MAXANGMOM
!> MXPRIM is the total number of (unique) primitive orbitals
INTEGER,intent(in)  :: MXPRIM
!> which atomic center the shell is attached to
INTEGER,intent(in)    :: SHELL2ATOM(MAXNSHELL)
!> the angular momentum for each shell
INTEGER,intent(in)    :: SHELLANGMOM(MAXNSHELL)
!> the number of primitives for each shell
INTEGER,intent(in)    :: SHELLNPRIM(MAXNSHELL)
!> the unique primitve exponents
REAL(REALK),intent(in):: PRIEXP(MXPRIM)
!> the index to start in PRIEXP(MXPRIM) for a given shell index 
INTEGER,intent(in) :: PRIEXPSTART(MAXNSHELL)
!> some grid generation specification
REAL(REALK),intent(in) :: radint
!> charge for each atom used in grid-generation
INTEGER,intent(in)  :: CHARGE(natoms)
!> which radial grid should be used?
INTEGER,intent(in) :: RADIALGRID
!> minimum radial grid points 
INTEGER,intent(in) :: MIN_RAD_PT
!> 
INTEGER :: NRADPT
!>
INTEGER,intent(inout)     :: nRadialPoints(NATOMS)                 
!>
REAL(REALK),intent(inout) :: RADIALPOINTS(NRADPT,NATOMS)
!>
REAL(REALK),intent(inout) :: RADIALWEIGHT(NRADPT,NATOMS) 
!
REAL(REALK),pointer :: AA(:,:,:)
INTEGER,pointer     :: NUCORB(:,:,:)
INTEGER :: IATOM

!USE ATOMTYPE IN THIS ROUTINE AS THIS IS THE SAME FOR ALL TYPES OF ATOMS.
SELECT CASE(RADIALGRID)
CASE(1) !GC2
   DO IATOM=1,NATOMS
      CALL GRID_RADGC2(CHARGE(IATOM),RADIALPOINTS(:,IATOM),RADIALWEIGHT(:,IATOM),nRadialPoints(IATOM),&
           & RADINT,NRADPT,MAXANGMOM,IPRINT,LUPRI)
   ENDDO
CASE(2) !LMG
   call mem_grid_alloc(AA,2,MAXANGMOM,NATOMS)
   call LS_DZERO(AA,2*MAXANGMOM*NATOMS)
   call mem_grid_alloc(NUCORB,MAXANGMOM,2,NATOMS)   
   call LS_IZERO(NUCORB,MAXANGMOM*2*NATOMS)
   !Build NUCORB and AA
   CALL GRID_NUCBAS(NUCORB,AA,SHELL2ATOM,SHELLANGMOM,SHELLNPRIM,MAXANGMOM,NATOMS,MAXNSHELL,&
        & MXPRIM,PRIEXP,PRIEXPSTART,LUPRI,IPRINT)
   !Build RADIALPOINTS,RADIALWEIGHT,nRadialPoints
   DO IATOM=1,NATOMS
      CALL GRID_RADLMG(RADIALPOINTS(:,IATOM),RADIALWEIGHT(:,IATOM),nRadialPoints(IATOM),RADINT,&
           & NRADPT,NUCORB(:,:,IATOM),AA(1:2,1:MAXANGMOM,IATOM),MAXANGMOM,IPRINT,LUPRI)
   ENDDO
   call mem_grid_dealloc(AA)
   call mem_grid_dealloc(NUCORB)
CASE(3) !TURBO
   DO IATOM=1,NATOMS
      CALL GRID_RADTURBO(CHARGE(IATOM),RADIALPOINTS(:,IATOM),RADIALWEIGHT(:,IATOM),nRadialPoints(IATOM),&
           & RADINT,NRADPT,MAXANGMOM,IPRINT,LUPRI,MIN_RAD_PT)
   ENDDO
END SELECT
!VERIFY THAT RADIALPOINTS HAS THE LARGEST VALUE FIRST
!CALL MAXLOC(RADIALPOINTS) IF not eq 1 call lsquit()
END SUBROUTINE SET_RADIAL

SUBROUTINE GRID_NUCBAS(NUCORB,AA,SHELL2ATOM,SHELLANGMOM,SHELLNPRIM,MAXANGMOM,NATOMS,MAXNSHELL,&
     & MXPRIM,PRIEXP,PRIEXPSTART,LUPRI,IPRINT)
  use precision
  IMPLICIT NONE
  !  Extract basis information for all atoms
  !  Written by T.Saue March 12 2001
  INTEGER,intent(in) :: MAXANGMOM,NATOMS,MAXNSHELL,IPRINT,MXPRIM
  INTEGER,intent(in) :: SHELL2ATOM(MAXNSHELL),SHELLANGMOM(MAXNSHELL),SHELLNPRIM(MAXNSHELL),PRIEXPSTART(MAXNSHELL)
  REAL(REALK),intent(in) :: PRIEXP(MXPRIM)
  REAL(REALK),intent(inout):: AA(2,MAXANGMOM,NATOMS)
  INTEGER,intent(inout)    :: NUCORB(MAXANGMOM,2,NATOMS)
  !
  REAL(REALK)           :: A,DUMMY
  REAL(REALK),PARAMETER :: D0=0.0E0_realk
  INTEGER  :: NDIM,JCENT,JC,JLVAL,ISHELL,ICENT,IC,ILVAL,IPRIM,JPRIM,I,IEXP,LL,L
  INTEGER  :: NPRIM,LUPRI
  NDIM = 2*MAXANGMOM*NATOMS
  CALL LS_IZERO(NUCORB,NDIM)
  DUMMY=1.0E+20_realk
  JCENT = 0
  JPRIM = -1
  JC = -1
  JLVAL = -1
  DO ISHELL = 1,MAXNSHELL
     ICENT = SHELL2ATOM(ISHELL)
     IF(ICENT.NE.JCENT) THEN
        JCENT = ICENT
        JLVAL = 0
     ENDIF
     IC = 1
     IF(IC.NE.JC) THEN
        JC    = IC
        JLVAL = 0
     ENDIF
     ILVAL = SHELLANGMOM(ISHELL)
     IF(ILVAL.NE.JLVAL) THEN
        JLVAL = ILVAL
        NUCORB(ILVAL,IC,ICENT) = 0
        AA(1,ILVAL,ICENT)=D0
        AA(2,ILVAL,ICENT)=DUMMY
     ENDIF
     NUCORB(ILVAL,IC,ICENT)=NUCORB(ILVAL,IC,ICENT)+1
     IPRIM = PRIEXPSTART(ISHELL)
     IF(IPRIM.NE.JPRIM) THEN
        JPRIM = IPRIM
        NPRIM = SHELLNPRIM(ISHELL)
        DO IEXP = 1,NPRIM
           A=PRIEXP(IPRIM+IEXP)            
           AA(1,ILVAL,ICENT)=MAX(AA(1,ILVAL,ICENT),A)
           AA(2,ILVAL,ICENT)=MIN(AA(2,ILVAL,ICENT),A)
        ENDDO
     ENDIF
  ENDDO

END SUBROUTINE GRID_NUCBAS

SUBROUTINE GRID_RADLMG(RADIALPOINTS,RADIALWEIGHT,nRadialPoints,RADINT,NRADPT,NUCORB,AA,MAXANGMOM,IPRINT,LUPRI)
  use precision
  IMPLICIT NONE
  INTEGER,intent(in) :: NRADPT,IPRINT,MAXANGMOM,LUPRI
  INTEGER,intent(in) :: NUCORB(MAXANGMOM,2)
  REAL(REALK),intent(in) :: AA(2,MAXANGMOM),RADINT
  INTEGER,intent(inout) :: nRadialPoints
  REAL(REALK),intent(inout) :: RADIALPOINTS(NRADPT),RADIALWEIGHT(NRADPT)
  !
  REAL(REALK), PARAMETER  :: D0 = 0E0_realk, D1 = 1E0_realk,D2 = 2E0_realk, D3 = 3E0_realk
  INTEGER     :: LL,L,NBAS,IR,I
  REAL(REALK) :: AH,H,EPH,RL,RH,AL,HTMP,RHTMP,GRDC,DUMMY
  REAL(realk) :: DISERR90,OUTERR90
  DUMMY=HUGE(AH)
  ! Grid spacing to H and inner grid point to AH
  nRadialPoints = 0
  H  = DUMMY
  AH = 0E0_realk
  DO LL = 1,MAXANGMOM
     L = LL-1
     NBAS=NUCORB(LL,1)+NUCORB(LL,2)
     IF(NBAS.GT. 0) THEN
        HTMP = DISERR90(RADINT,L) !function
        H = MIN(H,HTMP)
     ENDIF
     AH = MAX(AH,AA(1,LL))
  ENDDO
  IF(AH .EQ. 0E0_realk) RETURN
  EPH = EXP(H)
  AH = D2*AH
  RL = ((1.9E0_realk+LOG(RADINT))/D3)-(LOG(AH)/D2)
  RL = EXP(RL)
  RH = D0
  DO LL = 1,MAXANGMOM
     L = LL-1
     AL=DUMMY
     IF(NUCORB(LL,1).GT. 0) AL=AA(2,LL)
     IF(NUCORB(LL,2).GT. 0) AL=MIN(AL,D0)
     IF(AL.LT.DUMMY) THEN
        AL = AL+AL
        RHTMP = OUTERR90(AL,L,RADINT) !function          
        RH=MAX(RH,RHTMP)
     ENDIF
  ENDDO
  GRDC = RL/(EPH-D1)
  nRadialPoints = NINT(LOG(D1+(RH/GRDC))/H)
!  WRITE(LUPRI,'(A,I5)') 'Number of points:',nRadialPoints
  IF(nRadialPoints.GT.NRADPT) CALL LSQUIT('GRID_RADLMG: Too many radial points.',-1)
  RADIALPOINTS(nRadialPoints) = RL
  RADIALWEIGHT(nRadialPoints)  = (RL+GRDC)*RL*RL*H
  DO IR = nRadialPoints-1,1,-1
     RADIALPOINTS(IR) = (RADIALPOINTS(IR+1)+GRDC)*EPH-GRDC
     RADIALWEIGHT(IR) = (RADIALPOINTS(IR)+GRDC)*RADIALPOINTS(IR)*RADIALPOINTS(IR)*H
  ENDDO
END SUBROUTINE GRID_RADLMG

SUBROUTINE GRID_RADTURBO(CHARGE,RADIALPOINTS,RADIALWEIGHT,nRadialPoints,RADINT,NRADPT,MAXANGMOM,IPRINT,LUPRI,MIN_RAD_PT)
  use precision
  IMPLICIT NONE
  INTEGER,intent(in) :: NRADPT,IPRINT,MAXANGMOM,LUPRI,CHARGE,MIN_RAD_PT
  REAL(REALK),intent(in) :: RADINT
  INTEGER,intent(inout) :: nRadialPoints
  REAL(REALK),intent(inout) :: RADIALPOINTS(NRADPT),RADIALWEIGHT(NRADPT)
  !
  REAL(REALK), PARAMETER  :: D0 = 0E0_realk, D1 = 1E0_realk,D2 = 2E0_realk, D3 = 3E0_realk
  REAL(REALK), PARAMETER  :: PI=3.1415926535897932384626E0_realk
  REAL(REALK), PARAMETER  :: D06=0.6E0_realk,D04=0.4E0_realk,D05=0.5E0_realk
  REAL(REALK), PARAMETER  :: A=1E0_realk
  REAL(REALK) :: zeta,rfac,PiOverN,ANGL,X,s,w,aPlusX06,logAPlus1Over1MinusX,r,rdiff
  INTEGER :: I
  !H-Kr
  real(realk),parameter  :: zetas(36) = (/ 0.8E0_realk,  0.9E0_realk,  1.8E0_realk,  1.4E0_realk,  1.3E0_realk,&
       & 1.1E0_realk,  0.9E0_realk,  0.9E0_realk,  0.9E0_realk,  0.9E0_realk, 1.4E0_realk,  1.3E0_realk,  1.3E0_realk,&
       & 1.2E0_realk,  1.1E0_realk, 1.0E0_realk,  1.0E0_realk,  1.0E0_realk,  1.5E0_realk,  1.4E0_realk, 1.3E0_realk, &
       & 1.2E0_realk,  1.2E0_realk,  1.2E0_realk,  1.2E0_realk, 1.2E0_realk,  1.2E0_realk,  1.1E0_realk,  1.1E0_realk, &
       & 1.1E0_realk, 1.1E0_realk,  1.0E0_realk,  0.9E0_realk,  0.9E0_realk,  0.9E0_realk, 0.9E0_realk/)

  IF(charge.GE. 1.AND.charge.LE.size(zetas))THEN
     zeta = zetas(charge)
  ELSE
     zeta = 0.9E0_realk
  ENDIF
  rfac = zeta/LOG(D2)

  call GRID_TURBO_RADIALPOINTS(charge,RADINT,nRadialPoints,MIN_RAD_PT)
  IF(nRadialPoints.GT.NRADPT)CALL LSQUIT('something wrong in GRID_RADTURBO, nRadialPoints gt than allowed',lupri)
  PiOverN = PI/nRadialPoints
  !Radial points
  DO I=1,nRadialPoints
     ANGL = (I-1+D05)*PiOverN
     X = COS(angl)
     s = SIN(angl)
     w = PiOverN*s
     aPlusX06 = (A+X)**D06
     logAPlus1Over1MinusX = LOG((A+D1)/(D1-X))
     r = rfac*aPlusX06*logAPlus1Over1MinusX
     rdiff = rfac*(aPlusX06/(D1-X) + D06*logAPlus1Over1MinusX/((A+X)**D04))
     RADIALWEIGHT(I)=w*rdiff*r*r
     RADIALPOINTS(I) = r
  ENDDO

END SUBROUTINE GRID_RADTURBO

subroutine GRID_TURBO_RADIALPOINTS(Z,thrl,nRadialPoints,MIN_RAD_PT)
  implicit none
  integer,intent(out)     :: nRadialPoints
  integer,intent(in)     :: Z,MIN_RAD_PT
  real(realk),intent(in) :: thrl
  !
  integer :: ta,accuracy_correction,z_correction
  real(realk),parameter :: D5=5E0_realk,D3=3E0_realk,D05=0.5E0_realk
  real(realk) :: TMP
  
  IF(Z.LE. 2)THEN
     ta=0
  ELSEIF(Z.LE. 10)THEN
     ta=1
  ELSEIF(Z.LE. 18)THEN
     ta=2
  ELSEIF(Z.LE. 36)THEN
     ta=3
  ELSEIF(Z.LE. 54)THEN
     ta=4
  ELSEIF(Z.LE. 86)THEN
     ta=5
  ELSE
     ta=6
  ENDIF
  TMP = (-log10(thrl)-D5)*D3
  IF(TMP.LT. 0E0_realk)THEN
     CALL LSQUIT('Negative threshold in GRID_TURBO_RADIALPOINTS',-1)
  ENDIF
  accuracy_correction=NINT(TMP)
  IF(accuracy_correction.LT. 0)accuracy_correction=0
  z_correction = ta*5
  nRadialPoints = MIN_RAD_PT + accuracy_correction + z_correction
end subroutine GRID_TURBO_RADIALPOINTS

subroutine get_quadfilename(filename,nbast,node,gridid)
implicit none
character(len=22) :: filename
integer,intent(in) :: nbast,gridid
integer(kind=ls_mpik),intent(in) :: node
filename(1:11)='DALTON.QUAD'
#ifdef VAR_MPI 
filename(12:16)=Char(node/10000+48)//Char(mod(node,10000)/1000+48)&
   &//Char(mod(mod(node,10000),1000)/100+48)&
   &//Char(mod(mod(mod(node,10000),1000),100)/10+48)&
   &//Char(mod(mod(mod(mod(node,10000),1000),100),10)+48)
#else
filename(12:16)='00000'
#endif
filename(17:21)=Char(nbast/10000+48)//Char(mod(nbast,10000)/1000+48)&
   &//Char(mod(mod(nbast,10000),1000)/100+48)&
   &//Char(mod(mod(mod(nbast,10000),1000),100)/10+48)&
   &//Char(mod(mod(mod(mod(nbast,10000),1000),100),10)+48)
filename(22:22)=Char(GridId+48)
end subroutine get_quadfilename

SUBROUTINE GRID_RADGC2(CHARGE,RADIALPOINTS,RADIALWEIGHT,nRadialPoints,RADINT,NRADPT,MAXANGMOM,IPRINT,LUPRI)
  use precision
  IMPLICIT NONE
  INTEGER,intent(in) :: NRADPT,IPRINT,MAXANGMOM,LUPRI,CHARGE
  REAL(REALK),intent(in) :: RADINT
  INTEGER,intent(inout) :: nRadialPoints
  REAL(REALK),intent(inout) :: RADIALPOINTS(NRADPT),RADIALWEIGHT(NRADPT)
  !
  REAL(REALK), PARAMETER  :: D0 = 0E0_realk, D1 = 1E0_realk,D2 = 2E0_realk, D3 = 3E0_realk
  REAL(REALK), PARAMETER  :: PI=3.1415926535897932384626E0_realk,D16=16E0_realk
  REAL(REALK), PARAMETER  :: PI2=2E0_realk/3.1415926535897932384626E0_realk
  REAL(REALK), PARAMETER  :: SFAC=2E0_realk/3E0_realk
  REAL(REALK), PARAMETER  :: RFAC=1.442695040888963387E0_realk
  REAL(REALK) :: W,X,sinangl2,sinangl,ANGL,WFAC,N_PI,N_ONE
  INTEGER :: I

  call GRID_GC2_RADIALPOINTS(charge,RADINT,nRadialPoints)
  IF(nRadialPoints.GT.NRADPT)CALL LSQUIT('something wrong in GRID_RADGC2, nRadialPoints gt than allowed',lupri)
  N_ONE=nRadialPoints+D1
  N_PI=PI/N_ONE
  WFAC=D16/(D3*N_ONE)
  !Radial points
  DO I=1,nRadialPoints
     X = (nRadialPoints-1-2*(I-1))/N_ONE
     ANGL = N_PI*I
     sinangl = SIN(angl)
     sinangl2 = sinangl*sinangl
     X = X + PI2*(D1+SFAC*sinangl2)*COS(angl)*sinangl;
     RADIALPOINTS(I) = RFAC*LOG(D2/(D1-X))
     W = WFAC*sinangl2*sinangl2;
     RADIALWEIGHT(I)=W*RFAC/(D1-X)*RADIALPOINTS(I)*RADIALPOINTS(I)
     !transformation factor accumulated in weight
  ENDDO

END SUBROUTINE GRID_RADGC2

subroutine GRID_GC2_RADIALPOINTS(Z,thrl,nRadialPoints)
  implicit none
  integer,intent(out)     :: nRadialPoints
  integer,intent(in)     :: Z
  real(realk),intent(in) :: thrl
  !
  integer :: ta
  integer,parameter :: MIN_RAD_PT = 20
  real(realk) :: D32=3.2E0_realk,D18=1.8E0_realk,D4=4E0_realk
  
  IF(Z.LE. 2)THEN
     ta=0
  ELSEIF(Z.LE. 10)THEN
     ta=1
  ELSEIF(Z.LE. 18)THEN
     ta=2
  ELSEIF(Z.LE. 36)THEN
     ta=3
  ELSEIF(Z.LE. 54)THEN
     ta=4
  ELSEIF(Z.LE. 86)THEN
     ta=5
  ELSE
     ta=6
  ENDIF
  nRadialPoints = NINT(-D32*(D18*LOG10(thrl)-ta*D4))
  IF(nRadialPoints.LE.MIN_RAD_PT)nRadialPoints=MIN_RAD_PT
  !ri = rint( -3.2*(1.8*log10(thrl)-ta*4.0 ) );
  !return ri>MIN_RAD_PT ? ri : MIN_RAD_PT;
end subroutine GRID_GC2_RADIALPOINTS

SUBROUTINE DetermineMaxGridPoints(maxAtomicGridpoints,maxGridpoints,maxN,NATOMS,nRadialPoints,GRIDANG,NRADPT,tid,nthreads)
IMPLICIT NONE
INTEGER,intent(out) :: maxAtomicGridpoints,maxGridpoints,maxN
INTEGER,intent(in)  :: NATOMS,NRADPT
INTEGER,intent(in)  :: nRadialPoints(NATOMS),GRIDANG(NRADPT,NATOMS)
INTEGER,intent(in)  :: tid,nthreads
!
INTEGER :: IATOM,AtomicGridpoints,IPOINT,iang,N
integer,parameter :: Npoints(18) = (/  14,38,50,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454 /)
maxAtomicGridpoints = 0
maxGridpoints = 0
maxN = 0
DO IATOM=1+tid,NATOMS,nthreads
   AtomicGridpoints = 0
   DO IPOINT = 1,nRadialPoints(IATOM)
      iang = GRIDANG(IPOINT,IATOM)
      N=Npoints(iang)
      maxN = MAX(maxN,N)
      AtomicGridpoints = AtomicGridpoints+N
   ENDDO
   MaxAtomicGridpoints = MAX(MaxAtomicGridpoints,AtomicGridpoints)
   maxGridpoints = maxGridpoints+AtomicGridpoints
ENDDO
END SUBROUTINE DETERMINEMAXGRIDPOINTS

END MODULE GRIDGENERATIONMODULE
