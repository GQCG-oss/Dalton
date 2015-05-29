!> @file
!> Module contains grid generation routines
!> \brief 
!> \author T. Kjaergaard
!> \date 2009 
MODULE THCgridgenerationmodule
use gridgenerationmodule
use memory_handling
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
integer :: nTHCGridPoints
real(realk),pointer :: THCCOOR(:,:) !(3,nTHCGridPoints)
real(realk),pointer :: THCWG(:)

CONTAINS
!> \brief wrapper exchange-correlation integral routine that build basinf.
!> \author T. Kjaergaard
!> \date 2010
  SUBROUTINE THCGenerateGrid(NBAST,radint,angmin,angint,ihardness,iprune,natoms,& 
       & X,Y,Z,Charge,SHELL2ATOM,SHELLANGMOM,SHELLNPRIM,MAXANGMOM,&
       & MAXNSHELL,MXPRIM,PRIEXP,PRIEXPSTART,RSHEL,TURBO,&
       & RADIALGRID,ZdependenMaxAng,PARTITIONING,LUPRI,IPRINT,NGRID)
    IMPLICIT NONE
    !> number of gridpoints
    INTEGER,intent(inout) :: NGRID    
    !> the logical unit number for the output file
    INTEGER,intent(in) :: LUPRI
    !> the printlevel integer, determining how much output should be generated
    INTEGER,intent(in)  :: IPRINT
    !> number of basis functions
    INTEGER,intent(in)  :: NBAST
    !> maximum number of shells
    INTEGER,intent(in)  :: MAXNSHELL
    !> number of atoms
    INTEGER,intent(in)  :: natoms
    !> the maximum angular momentum
    INTEGER,intent(in)  :: MAXANGMOM
    !> MXPRIM is the total number of (unique) primitive orbitals
    INTEGER,intent(in)  :: MXPRIM
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
    REAL(REALK),intent(in) :: radint
    !> some grid generation specification
    INTEGER,intent(in) :: angmin
    !> some grid generation specification
    INTEGER,intent(in) :: angint
    !> some grid generation specification
    INTEGER,intent(in) :: IHARDNESS
    !> some grid generation specification
    INTEGER,intent(in)  :: TURBO
    !> grid point pruning (true if iprune=1) 
    INTEGER,intent(in) :: iprune
    !> which radial grid should be used?
    INTEGER,intent(in) :: RADIALGRID !(1=GC2,2=LMG,3=TURBO) 
    !> should the maximum angular momentum be Z dependent
    LOGICAL,intent(in) :: ZdependenMaxAng
    !> which partitioning should be used
    INTEGER,intent(in) :: PARTITIONING !(1=SSF,2=BECKE,3=BECKEORIG,4=BLOCK,...)
    !
    INTEGER,pointer     :: nRadialPoints(:),GRIDANG(:,:)
    REAL(REALK),pointer :: RADIALPOINTS(:,:),RADIALWEIGHT(:,:)
    integer :: totalpoints
    integer :: angularpoints ,I,leb_gen_from_point,iang1

    call init_gridmemvar()
    iang1 = leb_get_from_order(angint)
    angularpoints = leb_gen_points(iang1)
    !mem alloc RADIALPOINTS,RADIALWEIGHT,nRadialPoints=RadialGridPoints,RadialWeights,NumberOfGridPoints 
    call mem_alloc(RADIALPOINTS,NRADPT,NATOMS)
    call mem_alloc(RADIALWEIGHT,NRADPT,NATOMS)
    call mem_alloc(nRadialPoints,NATOMS)   
    !Build RADIALPOINTS,RADIALWEIGHT,nRadialPoints = RadialGridPoints,RadialWeights,NumberOfGridPoints 
    CALL SET_RADIAL(RADIALPOINTS,RADIALWEIGHT,nRadialPoints,NRADPT,SHELL2ATOM,SHELLANGMOM,&
         & SHELLNPRIM,MAXANGMOM,NATOMS,MAXNSHELL,MXPRIM,PRIEXP,PRIEXPSTART,IPRINT,LUPRI,RADINT,&
         & Charge,RADIALGRID)
    ! Computes the angular points for a set of radial points 
    ! obtained from radial integration scheme.
    call mem_alloc(GRIDANG,NRADPT,NATOMS)
    CALL SET_ANGULAR(gridang,radint,angmin,angint,CHARGE,natoms,NRADPT,nRadialPoints,RADIALPOINTS,&
         & RADIALWEIGHT,iprune,angularpoints,ZdependenMaxAng)
    ! Build Atomic and Molecular grids + partitioning
    CALL ComputeCoordsTHC(NGRID,totalpoints,nRadialPoints,natoms,RADIALPOINTS,RADIALWEIGHT,NRADPT,GRIDANG,&
         & X,Y,Z,ihardness,MAXNSHELL,RSHEL,SHELL2ATOM,NBAST,PARTITIONING,Charge,iprint,lupri)
    call mem_dealloc(RADIALPOINTS)
    call mem_dealloc(RADIALWEIGHT)
    call mem_dealloc(nRadialPoints)   
    call mem_dealloc(GRIDANG)
    call stats_grid_mem(lupri)
  END SUBROUTINE THCGENERATEGRID

  SUBROUTINE ComputecoordsTHC(NGRID,totalpoints,nRadialPoints,Natoms,RADIALPOINTS,RADIALWEIGHT,NRADPT,&
       & GRIDANG,atomcenterX,atomcenterY,atomcenterZ,ihardness,MAXNSHELL,RSHEL,SHELL2ATOM,NBAST,&
       & PARTITIONING,Charge,iprint,lupri)
    implicit none
    !> number of gridpoints
    INTEGER,intent(inout) :: NGRID    
    integer,intent(in) :: Natoms,NRADPT,ihardness,iprint,lupri,nbast
    !> X coordinate for each atom
    REAL(REALK),intent(in):: atomcenterX(natoms)
    !> Y coordinate for each atom
    REAL(REALK),intent(in):: atomcenterY(natoms)
    !> Z coordinate for each atom
    REAL(REALK),intent(in):: atomcenterZ(natoms)
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
    real(realk),pointer :: inverseDistance12(:,:) ,BFAC(:,:)
    real(realk),pointer :: SSF_ATOMIC_CUTOFF(:)
    real(realk),parameter :: D1=1E0_realk,D0=0E0_realk,D05=0.5E0_realk,D3=3E0_realk
    real(realk) :: W,RW,fac,center(3),dist
    integer :: ikey,npoints2,NSHELLBLOCKS,LUGRID,IT,NGRIDPOINTINBATCH,nthreads,tid
    integer :: maxN,MaxAtomicGridpoints,maxGridpoints,AtomicGridpoints,nprocessors
    integer :: GlobalmaxAtomicGridpoints,GlobalmaxGridpoints,mynum
    integer,pointer :: SHELLBLOCKS(:,:),ATOMIDX(:)
    integer :: nx,ny,nz,nkey,mykeysNumber,IMykey,privatetotalpoints
    logical :: unique,postprocess

    call mem_alloc(inverseDistance12,NATOMS,NATOMS)
    !This subroutine actually the inversedistance 1/r12
    CALL determine_distance12(inverseDistance12,NATOMS,atomcenterX,atomcenterY,atomcenterZ)

    IF(PARTITIONING.EQ. 2.OR.PARTITIONING.EQ. 4)THEN !needed for becke preprocess
       call mem_alloc(Bfac,NATOMS,NATOMS)
       CALL Determine_BraggFac(NATOMS,BFAC,CHARGE)
    ENDIF
    IF(PARTITIONING.EQ. 5.OR.PARTITIONING.EQ. 1)THEN !needed for block-SSF postprocess
       call mem_alloc(SSF_ATOMIC_CUTOFF,NATOMS)
       CALL Determine_SSF_ATOMIC_CUTOFF(NATOMS,SSF_ATOMIC_CUTOFF,&
            & atomcenterX,atomcenterY,atomcenterZ)
    ENDIF
    postprocess=.FALSE.
    IF(PARTITIONING.EQ. 4.OR.PARTITIONING.EQ. 5)THEN
       call lsquit('block post process deactivated',-1)
    ENDIF

    nthreads=1
    tid=0
    call DetermineMaxGridPointsA(GlobalmaxGridpoints,&
         & NATOMS,nRadialPoints,GRIDANG,NRADPT)
    call mem_alloc(WG2,GlobalmaxGridpoints)
    call mem_alloc(COOR2,3,GlobalmaxGridpoints)
    totalpoints=0

    call DetermineMaxGridPointsB(maxAtomicGridpoints,maxGridpoints,maxN,&
         & NATOMS,nRadialPoints,GRIDANG,NRADPT)

    call mem_alloc(locR,4*maxN)
    call mem_alloc(COOR,3,MaxAtomicGridpoints)
    call mem_alloc(WG,MaxAtomicGridpoints)

    print*,'THC: PARTITIONING',PARTITIONING

    DO IATOM=1,NATOMS
       Xatomcoor = atomcenterX(IATOM)
       Yatomcoor = atomcenterY(IATOM)
       Zatomcoor = atomcenterZ(IATOM)
!       print*,'THC: IATOM',IATOM,'atomcenterX(IATOM)',atomcenterX(IATOM)
!       print*,'THC: atomcenterY(IATOM)',atomcenterY(IATOM)
!       print*,'THC: atomcenterZ(IATOM)',atomcenterZ(IATOM)
       idx = 0
       DO IPOINT = 1,nRadialPoints(IATOM)
          iang = GRIDANG(IPOINT,IATOM)
          radial = RADIALPOINTS(IPOINT,IATOM)
          weight = RADIALWEIGHT(IPOINT,IATOM)
          FACTOR = pim4*weight

!          print*,'THC: IPOINT',IPOINT,'GRIDANG(IPOINT,IATOM)',GRIDANG(IPOINT,IATOM)
!          print*,'THC: RADIALPOINTS(IPOINT,IATOM)',RADIALPOINTS(IPOINT,IATOM)
!          print*,'THC: RADIALWEIGHT(IPOINT,IATOM)',RADIALWEIGHT(IPOINT,IATOM)
!          print*,'THC: pim4*weight',pim4*weight

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
          print*,'THC: ELSE VERSION'
          DO I=igrid,totalAtomicpoints
             IF(ABS(WG(I)).GT.compress_thr)THEN
                ITMP=ITMP + 1
                COOR(1,iTMP) = COOR(1,I)
                COOR(2,iTMP) = COOR(2,I)
                COOR(3,iTMP) = COOR(3,I)
                WG(iTMP) = WG(I)
             ENDIF
          ENDDO
          totalAtomicpointscompressed = ITMP
       endif
       !   IF(iprint.gt. 1)THEN
!!$OMP CRITICAL (GRIDWRITE)
       WRITE(Lupri,'(A,I5,A,I6,A,I6,A,I5,A)')'Atom: ',IATOM,' points=',&
            & totalAtomicpointscompressed,' compressed from',&
            & totalAtomicpoints,'(',nRadialPoints(IATOM),' radial)'
!!$OMP END CRITICAL (GRIDWRITE)
       !   ENDIF
!!$OMP CRITICAL 
       privatetotalpoints = totalpoints
       totalpoints = totalpoints + totalAtomicpointscompressed   
!!$OMP END CRITICAL 
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
       print*,'THC: postprocess',postprocess
    ENDDO

    call mem_dealloc(locR)
    call mem_dealloc(COOR)
    call mem_dealloc(WG)

    WRITE(Lupri,'(A,I12)')'Total Number of grid points:',totalpoints

    call mem_dealloc(inverseDistance12)
    IF(PARTITIONING.EQ. 2.OR.PARTITIONING.EQ. 4)THEN !needed for becke preprocess
       call mem_dealloc(Bfac)
    ENDIF
    IF(PARTITIONING.EQ. 5.OR.PARTITIONING.EQ. 1)THEN
       call mem_dealloc(SSF_ATOMIC_CUTOFF)
    ENDIF

    nTHCGridPoints = totalpoints
    call mem_alloc(THCCOOR,3,totalpoints)
    call mem_alloc(THCWG,totalpoints)
    DO I=1,totalpoints
       THCCOOR(1,I) = COOR2(1,I)
       THCCOOR(2,I) = COOR2(2,I)
       THCCOOR(3,I) = COOR2(3,I)
       THCWG(I) = WG2(I)
    ENDDO
    call mem_dealloc(WG2)
    call mem_dealloc(COOR2)
    NGRID = nTHCGridPoints
  END SUBROUTINE ComputecoordsTHC

SUBROUTINE DetermineMaxGridPointsB(maxAtomicGridpoints,maxGridpoints,maxN,&
     & NATOMS,nRadialPoints,GRIDANG,NRADPT)
IMPLICIT NONE
INTEGER,intent(out) :: maxAtomicGridpoints,maxGridpoints,maxN
INTEGER,intent(in)  :: NATOMS,NRADPT
INTEGER,intent(in)  :: nRadialPoints(NATOMS),GRIDANG(NRADPT,NATOMS)
!
INTEGER :: IATOM,AtomicGridpoints,IPOINT,iang,N
INTEGER,parameter :: Npoints(18) = (/ 14,38,50,86,110,146,170,&
     &194,230,266,302,350,434,590,770,974,1202,1454/)
maxAtomicGridpoints = 0
maxGridpoints = 0
maxN = 0
DO IATOM=1,NATOMS
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
END SUBROUTINE DETERMINEMAXGRIDPOINTSB

SUBROUTINE DetermineMaxGridPointsA(maxGridpoints,NATOMS,nRadialPoints,GRIDANG,NRADPT)
IMPLICIT NONE
INTEGER,intent(out) :: maxGridpoints
INTEGER,intent(in)  :: NATOMS,NRADPT
INTEGER,intent(in)  :: nRadialPoints(NATOMS),GRIDANG(NRADPT,NATOMS)
!
INTEGER :: IATOM,AtomicGridpoints,IPOINT,iang
INTEGER,parameter :: Npoints(18) = (/ 14,38,50,86,110,146,170,&
     & 194,230,266,302,350,434,590,770,974,1202,1454/)
maxGridpoints = 0
DO IATOM=1,NATOMS
   AtomicGridpoints = 0
   DO IPOINT = 1,nRadialPoints(IATOM)
      iang = GRIDANG(IPOINT,IATOM)
      AtomicGridpoints = AtomicGridpoints+Npoints(iang)
   ENDDO
   maxGridpoints = maxGridpoints+AtomicGridpoints
ENDDO
END SUBROUTINE DETERMINEMAXGRIDPOINTSA

END MODULE THCGRIDGENERATIONMODULE
