!> @file
!> Module containing main exchange-correlation integral driver, and routines to evaluate AOs and electron-densities
!> \brief 
!> \author T. Kjaergaard
!> \date 2009 
MODULE IIABSVALINT
use gridgenerationmodule
use dft_memory_handling
use LSparameters
!use memory_handling
!WARNING you must not add memory_handling, all memory goes through 
!grid_memory_handling  module so as to determine the memory used in this module.
use precision
use TYPEDEF
use dft_type
use dft_typetype
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
integer :: ABSVAL_MXBLLEN 
private 
public :: II_ABSVALINT
CONTAINS
!> \brief wrapper exchange-correlation integral routine that build basinf.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_ABSVALINT(LUPRI,IPRINT,SETTING,CMAT1,CMAT2,NBAST,NMO1,NMO2,&
     & ABSVALOVERLAP,USE_MPI,DFTHRI,SameCmat)
use BUILDAOBATCH
IMPLICIT NONE
INTEGER,intent(in)     :: LUPRI,IPRINT,NBAST,NMO1,nMO2
REAL(REALK),intent(in) :: CMAT1(NBAST,NMO1)
REAL(REALK),intent(in) :: CMAT2(NBAST,NMO2),DFTHRI
REAL(REALK),intent(inout) :: ABSVALOVERLAP(NMO1,NMO2)
TYPE(LSSETTING) :: SETTING
LOGICAL         :: USE_MPI !use MPI ?
LOGICAL,intent(in) :: SameCmat !Cmat1 and Cmat2, same?
TYPE(BASINF)  :: BAS
integer :: igrid,GRIDDONE
REAL(REALK),parameter :: D0=0E0_realk
INTEGER       :: maxNactbast,GRIDITERATIONS
igrid = SETTING%scheme%DFT%igrid
GRIDDONE = SETTING%scheme%DFT%GridObject(igrid)%GRIDDONE
maxNactbast = dft_maxNactbast(igrid)
GRIDITERATIONS = dft_GRIDITERATIONS(igrid)
CALL BUILD_BASINF(LUPRI,IPRINT,BAS,SETTING,GRIDDONE,.FALSE.)
CALL II_ABSVALINT1(LUPRI,IPRINT,CMAT1,CMAT2,NBAST,NMO1,nMO2,BAS,&
     &ABSVALOVERLAP,SETTING%SCHEME%noOMP,USE_MPI,setting%numnodes,setting%node,&
     &DFTHRI,SETTING%scheme%DFT%GridObject(igrid),maxNactbast,&
     &GRIDITERATIONS,SameCmat)
CALL FREE_BASINF(BAS)
dft_maxNactbast(igrid) = maxNactbast
dft_GRIDITERATIONS(igrid) = GRIDITERATIONS
END SUBROUTINE II_ABSVALINT

SUBROUTINE SET_ABSVAL_MXBLLEN(NBAST)
implicit none
integer :: NBAST
ABSVAL_MXBLLEN=NBAST
END SUBROUTINE SET_ABSVAL_MXBLLEN

SUBROUTINE II_ABSVALINT1(LUPRI,IPRINT,CMAT1,CMAT2,NBAST,NMO1,nMO2,BAS,&
     &ABSVALOVERLAP,noOMP,USE_MPI,numnodes,node,DFTHRI,GridObject,&
     &maxNactbast,GRIDITERATIONS,SameCmat) 
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in)  :: IPRINT
!> number of basis functions
INTEGER,intent(in)  :: NBAST
!> number of MO coeff
INTEGER,intent(in)  :: NMO1,NMO2
!> the MO coeff matrix 1
REAL(REALK),intent(in) :: CMAT1(NBAST,NMO1) !ao,mo
!> the MO coeff matrix 2
REAL(REALK),intent(in) :: CMAT2(NBAST,NMO2) !ao,mo
!> Basis-set information
TYPE(BASINF),intent(INOUT)  :: BAS
!> thresholds
REAL(REALK),intent(in) :: DFTHRI
!> same Cmat1 and Cmat2
logical,intent(in) :: SameCmat
!> 
real(realk),intent(inout) :: ABSVALOVERLAP(NMO1,NMO2)
!> shoould OpenMP be deactivated
LOGICAL,intent(in)  :: noOMP
!> should we use MPI
LOGICAL,intent(in) :: USE_MPI
!> MPI info
INTEGER(kind=ls_mpik),intent(in) :: numnodes,node
!> Grid parameters
type(Griditem),intent(INOUT) :: GridObject
integer :: maxNactBAST,GRIDITERATIONS
!
integer :: NPOINTS,NLEN,NSHELLBLOCKS,NROW,NCOL
INTEGER :: iprune,L_prev,L_curr,IPT,spSIZE,L,NCURLEN,I,J
INTEGER,pointer     :: SPINDEX(:)
INTEGER   :: KCKTA,KHKTA,SHELLANGMOMA,IT,IBUF,IBUF_PREV,IDUM,ILEN,K,NactBAS,NRED,NRED2,IDMAT
INTEGER   :: spsize2,IJ
REAL(REALK) :: ERROR,TS,TE
REAL(REALK),pointer :: SPHMAT(:)
LOGICAL     :: CHECKELS,SETIT,LDUM,PRINTTIM
integer,pointer :: LVALUE(:,:),MVALUE(:,:),NVALUE(:,:)
integer :: mmm

IT=0
SETIT=.FALSE.
IF(GridObject%GRIDDONE .EQ. 0)THEN
   GRIDITERATIONS=0
   SETIT=.TRUE.
   GridObject%NBUFLEN=0
   maxNactBAST = 0
ENDIF
!pruning: per default on
IPRUNE = 1
IF (GridObject%NOPRUN) IPRUNE = 0
IF(GridObject%GRIDDONE.EQ. 0)THEN
   PRINTTIM=.TRUE.
ELSE
   PRINTTIM=.FALSE.
ENDIF

GridObject%NBUFLEN=1024

!call SET_ABSVAL_MXBLLEN(NBAST)
call SET_ABSVAL_MXBLLEN(128)
!GridObject%NBUFLEN=NBAST
BoxMemRequirement = 5*NBAST*NBAST
CALL LSTIMER('START',TS,TE,LUPRI)
CALL GenerateGrid(NBAST,GridObject%radint,GridObject%angint,GridObject%HRDNES,iprune,BAS%natoms,& 
     & BAS%X,BAS%Y,BAS%Z,BAS%Charge,GridObject%GRIDDONE,BAS%SHELL2ATOM,BAS%SHELLANGMOM,BAS%SHELLNPRIM,BAS%MAXANGMOM,&
     & BAS%MAXNSHELL,BAS%MXPRIM,BAS%PRIEXP,BAS%PRIEXPSTART,BAS%RSHEL,IT,GridObject%TURBO,GridObject%nbuflen,&
     & GridObject%RADIALGRID,GridObject%ZdependenMaxAng,GridObject%PARTITIONING,BAS%nstart,MaxNactBast,LUPRI,&
     & IPRINT,USE_MPI,numnodes,node,GridObject%Id,GridObject%numnodes)
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

IF(BAS%MAXANGMOM .GE. 3)THEN
   spSIZE=0
   DO L=2,BAS%MAXANGMOM-1
      spSIZE=spSIZE+(2*L+1)*(L+1)*(L+2)/2
   ENDDO
   call mem_dft_alloc(SPHMAT,spSIZE)
   CALL LS_DZERO(SPHMAT,spSIZE)
   call mem_dft_alloc(SPINDEX,BAS%MAXANGMOM)
   spsize2 = BAS%maxangmom
   CALL Build_PRECALCULATED_SPHMAT_ABSVAL(LUPRI,BAS%MAXANGMOM-1,spSIZE,SPHMAT,SPINDEX)
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

BoxMemRequirement = ABSVAL_MXBLLEN*MaxNactBast
CALL II_ABSVAL_GRID_LOOP(GridObject%NBUFLEN,IT,nbast,nmo1,nmo2,maxNactbast,&
     & BAS%maxnshell,BAS%CC,BAS%ushells,BAS%CCSTART,BAS%CCINDEX,&
     & BAS%CENT,ABSVALOVERLAP,DFThri,CMAT1,CMAT2,BAS%PRIEXPSTART,lupri,&
     & BAS%mxprim,BAS%shellangmom,MMM,BAS%maxangmom,BAS%nstart,BAS%priexp,&
     & spsize,sphmat,spsize2,spindex,LVALUE,MVALUE,NVALUE,noOMP,&
     & node,GridObject%Id,SameCmat)
call mem_dft_dealloc(SPHMAT)
call mem_dft_dealloc(SPINDEX)
call mem_dft_dealloc(LVALUE)
call mem_dft_dealloc(MVALUE)
call mem_dft_dealloc(NVALUE)
END SUBROUTINE II_ABSVALINT1

SUBROUTINE II_ABSVAL_GRID_LOOP(NBUFLEN,IT,nbast,nmo1,nmo2,maxNactbast,maxnshell,&
     & CC,ushells,CCSTART,CCINDEX,CENT,ABSVALOVERLAP,&
     & DFThri,CMAT1,CMAT2,PRIEXPSTART,lupri,&
     & mxprim,shellangmom,MMM,maxangmom,nstart,priexp,spsize,&
     & sphmat,spsize2,spindex,LVALUE,MVALUE,NVALUE,noOMP,node,&
     & GridId,SameCmat)
implicit none
integer,intent(in) :: LUPRI,MAXNSHELL,nbast,NBUFLEN,gridid,nmo1,nmo2
integer :: IT
!> Largest number of active orbitals
integer,intent(inout) ::  maxNactbast
!> MXPRIM is the total number of (unique) primitive orbitals
INTEGER,intent(in)  :: MXPRIM
!> the angular momentum for each shell
INTEGER,intent(in)    :: SHELLANGMOM(MAXNSHELL)
!> the maximum angular momentum
INTEGER,intent(in)  :: MAXANGMOM
!> MMM=MAXANGMOM*(MAXANGMOM+1)/2*MAXANGMOM*(MAXANGMOM+1)/2 
INTEGER,intent(in)  :: MMM
real(realk),intent(in) :: CMAT1(NBAST,NMO1)
real(realk),intent(in) :: CMAT2(NBAST,NMO2)
logical,intent(in)     :: SameCmat
!> the number of unique shells (also different contraction coefficients matrices) 
INTEGER,intent(in)  :: ushells
!> a contraction coefficient matrix for all the unique shells
TYPE(LSMATRIX),intent(in):: CC(ushells)
!> hmmm I do not remember
INTEGER,intent(in) :: CCSTART(ushells)
!> contraction index for given shell, for shellindex => unique shellindex
INTEGER,intent(in) :: CCINDEX(MAXNSHELL)
!> X,Y,Z coordinate for each shell
REAL(REALK),intent(in):: CENT(3,MAXNSHELL)
real(realk),intent(inout) :: ABSVALOVERLAP(NMO1,NMO2)
!> threshold for value of GAOs
REAL(REALK),intent(in) :: DFTHRI
!> shoould OpenMP be deactivated
LOGICAL  :: noOMP
!> the index to start in PRIEXP(MXPRIM) for a given shell index 
INTEGER,intent(in) :: PRIEXPSTART(MAXNSHELL)
!> for a given shell index it gives the corresponding starting orbital
INTEGER,intent(in)    :: NSTART(MAXNSHELL)
!> the unique primitve exponents
REAL(REALK),intent(in) :: PRIEXP(MXPRIM)
!> size of sphmat
integer,intent(in) :: spsize
!> size of spindex
integer,intent(in) :: spsize2
!> size of spherical transformation matrices
real(realk),intent(in) :: sphmat(spsize)
!> use the new grid?
integer,intent(in) :: spindex(spsize2)
integer,intent(in) :: LVALUE(MMM,MAXANGMOM),MVALUE(MMM,MAXANGMOM),NVALUE(MMM,MAXANGMOM)
integer(kind=ls_mpik),intent(in) ::  node
!
INTEGER :: NactBAS
INTEGER,pointer :: SHELLBLOCKS(:,:),BLOCKS(:,:)
real(realk),pointer :: WEIGHT(:)
real(realk),pointer :: COOR(:,:)
real(realk),pointer :: ACTIVE_CMAT1(:)
real(realk),pointer :: ACTIVE_CMAT2(:)
real(realk),pointer :: GAO(:),TMP(:),TMP2(:)
REAL(REALK),pointer :: COOR_pt(:,:)
integer :: XX,NSHELLBLOCKS,NLEN,IPT,NCURLEN,I
integer :: W1,W2,W3,W4,W5,W6,W7,W8,nthreads,tid,lugrid,idmat1,idmat2
integer :: myBoxMemRequirement
real(realk),pointer :: myABSVALOVERLAP(:,:)
character(len=22) :: filename
logical :: grid_exists
#ifdef VAR_OMP
integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif

LUGRID=-1
call get_quadfilename(filename,nbast,node,GridId)
INQUIRE(file=filename,EXIST=grid_exists)
IF(grid_exists)THEN 
   CALL LSOPEN(LUGRID,filename,'OLD','UNFORMATTED')
ELSE
   print*,'file with filename: ',filename,' does not exist'
   print*,'but it should exist'
   call lsquit('missing file in XC integration',-1)
ENDIF
IF(.NOT.noOMP) call mem_dft_TurnONThread_Memory()
!$OMP PARALLEL IF(.NOT.noOMP) PRIVATE(XX,NSHELLBLOCKS,SHELLBLOCKS,COOR,COOR_pt,WEIGHT,IPT,NCURLEN,&
!$OMP ACTIVE_CMAT1,ACTIVE_CMAT2,BLOCKS,GAO,NLEN,NactBAS,I,myABSVALOVERLAP,&
!$OMP tid,nthreads,TMP,TMP2) SHARED(ushells,CC,CCSTART,CCINDEX,CENT,Cmat1,Cmat2,&
!$OMP PRIEXPSTART,shellangmom,MMM,maxangmom,nstart,priexp,spsize,spsize2,sphmat,&
!$OMP spindex,ABSVALOVERLAP,noOMP,NBUFLEN,it,lupri,nbast,nmo1,nmo2,maxnshell,mxprim,&
!$OMP dfthri,LVALUE,MVALUE,NVALUE,MaxNactBast,lugrid,SameCmat)
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
call mem_dft_alloc(WEIGHT,NBUFLEN)
call mem_dft_alloc(COOR,3,NBUFLEN)
call mem_dft_alloc(ACTIVE_CMAT1,maxNactBAST*NMO1)
IF(.NOT.SameCmat)THEN
   call mem_dft_alloc(ACTIVE_CMAT2,maxNactBAST*NMO2)
ENDIF
call mem_dft_alloc(GAO,ABSVAL_MXBLLEN*maxNactBAST) 
call mem_dft_alloc(TMP,ABSVAL_MXBLLEN*NMO1) 
call mem_dft_alloc(TMP2,ABSVAL_MXBLLEN*NMO2) 
call mem_dft_alloc(myABSVALOVERLAP,NMO1,NMO2)
call ls_dzero(myABSVALOVERLAP,NMO1*NMO2)
!$OMP END CRITICAL (initdftDATAblock)

DO XX=1+tid,IT,nthreads
!$OMP CRITICAL (reaquaREAD)
   CALL READ_GRIDPOINTS(LUGRID,NSHELLBLOCKS,SHELLBLOCKS,MAXNSHELL,NBUFLEN,COOR,WEIGHT,NLEN)
!$OMP END CRITICAL (reaquaREAD)
#ifdef VAR_LSDEBUGINT
   IF(NLEN .LE. 0) THEN
      CALL LSQUIT('SOMETHING WRONG IN THE DFT LOOP',lupri)
   ENDIF
#endif
   CALL SHELL_TO_ABSVAL_ORB(LUPRI,NSHELLBLOCKS,SHELLBLOCKS,BLOCKS,NSTART,MAXNSHELL,NBAST) !OUTPUT: BLOCKS
   CALL DETERMINE_ABSVAL_NACTIVEORB(NactBAS,NSHELLBLOCKS,BLOCKS,MAXNSHELL) !OUTPUT: NactBAS
   IF(SameCmat)THEN
      CALL CONSTRUCT_ACTIVE_CMAT(NSHELLBLOCKS,MAXNSHELL,BLOCKS,CMAT1,NBAST,NMO1,ACTIVE_CMAT1,NactBAS)
   ELSE
      CALL CONSTRUCT_ACTIVE_CMAT12(NSHELLBLOCKS,MAXNSHELL,BLOCKS,CMAT1,CMAT2,NBAST,&
           & NMO1,NMO2,ACTIVE_CMAT1,ACTIVE_CMAT2,NactBAS)
   ENDIF
   CALL SHELL_TO_ABSVAL_ACTORB(NSHELLBLOCKS,BLOCKS,MAXNSHELL) !CHANGE ORBBLOCKS
   !CHANGE TO NO LOOP MAYBE
   DO IPT = 1, NLEN, ABSVAL_MXBLLEN
      NCURLEN=MIN(ABSVAL_MXBLLEN,NLEN-IPT+1)
      COOR_pt => COOR(:,IPT:(IPT+NCURLEN-1))
      CALL II_ABSVAL_BLGETSOS(LUPRI,NCURLEN,GAO(1:NCURLEN*NactBas),COOR_pt,&
      &                NSHELLBLOCKS,SHELLBLOCKS,NactBAS,DFTHRI,0,&
      &                MMM,MAXANGMOM,spSIZE,SPHMAT,SPINDEX,MAXNSHELL,NSTART,SHELLANGMOM&
      &                ,CENT,PRIEXPSTART,CC,ushells,CCSTART,CCINDEX,MXPRIM,&
      &                PRIEXP,LVALUE,MVALUE,NVALUE)
      IF(SameCmat)THEN
         CALL ABSVALOVERLAP_WORKER(LUPRI,NCURLEN,NactBas,NBAST,NMO1,&
              & WEIGHT(IPT:IPT+NCURLEN-1),GAO(1:NCURLEN*NactBas),ACTIVE_CMAT1,myABSVALOVERLAP,&
              & TMP(1:NCURLEN*NMO1),TMP2(1:NCURLEN*NMO1))  
      ELSE
         CALL ABSVALOVERLAP_WORKER12(LUPRI,NCURLEN,NactBas,NBAST,NMO1,NMO2,&
              & WEIGHT(IPT:IPT+NCURLEN-1),GAO(1:NCURLEN*NactBas),&
              & ACTIVE_CMAT1,ACTIVE_CMAT2,myABSVALOVERLAP,&
              & TMP(1:NCURLEN*NMO1),TMP2(1:NCURLEN*NMO2))  
      ENDIF
   ENDDO
ENDDO

!$OMP CRITICAL (freeDFTdatablock)
call mem_dft_dealloc(TMP) 
call mem_dft_dealloc(TMP2) 
call mem_dft_dealloc(WEIGHT)
call mem_dft_dealloc(COOR)
call mem_dft_dealloc(ACTIVE_CMAT1)
IF(.NOT.SameCmat)THEN
   call mem_dft_dealloc(ACTIVE_CMAT2)
ENDIF
call mem_dft_dealloc(GAO) 
call mem_dft_dealloc(SHELLBLOCKS)
call mem_dft_dealloc(BLOCKS)
CALL DAXPY(nmo1*nmo2,1E0_realk,myABSVALOVERLAP,1,ABSVALOVERLAP,1)
call mem_dft_dealloc(myABSVALOVERLAP)
!$OMP END CRITICAL (freeDFTdatablock)
IF(.NOT.noOMP) call collect_thread_dft_memory()
!$OMP END PARALLEL
IF(.NOT.noOMP) call mem_dft_TurnOffThread_Memory()

CALL LSCLOSE(LUGRID,'KEEP')

END SUBROUTINE II_ABSVAL_GRID_LOOP

SUBROUTINE ABSVALOVERLAP_WORKER(LUPRI,NBLEN,Nactbast,NBAST,&
     &  NMO,WEIGHT,GAOS,CMAT,ABSVALOVERLAP,TMP,TMP2)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of MO coef
INTEGER,intent(in) :: NMO
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: WEIGHT(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST)
!> The active AO to MO coeff matrix
REAL(REALK),intent(inout) :: CMAT(NACTBAST,NMO)
!> The absolute valued Overlap matrix in MO
REAL(REALK),intent(inout) :: ABSVALOVERLAP(NMO,NMO)
! TEMPORARY MEM FROM WORK
REAL(REALK),intent(inout) :: TMP(NBLEN,NMO) 
REAL(REALK),intent(inout) :: TMP2(NBLEN,NMO)
!
Integer :: J,K,IBL,IACT
!MO Values on grid point 
!TMP(NBLEN,NMO) = GAO(NBLEN,NACTBAST)*CMAT(NACTBAST,NMO)
CALL DGEMM('N','N',NBLEN,NMO,NACTBAST,1.0E0_realk,&
     &             GAOS,NBLEN,CMAT,NACTBAST,0.0E0_realk,&
     &             TMP,NBLEN)

! First half-contraction of MO's with grid point weight
DO J=1,NMO
   DO K=1, NBLEN
      TMP(K,J) = ABS(TMP(K,J))
   ENDDO
   DO K=1, NBLEN
      TMP2(K,J)= WEIGHT(K)*TMP(K,J)
   ENDDO
ENDDO
!  Second half-contraction of MO's with grid point weight
!ABSVALOVERLAP(NMO,NMO) =+ op(TMP)(NBLEN,NMO)*TMP2(NBLEN,NMO)
CALL DGEMM('T','N',NMO,NMO,NBLEN,1.0E0_realk,&
     &             TMP,NBLEN,TMP2,NBLEN,1.0E0_realk,&
     &             ABSVALOVERLAP,NMO)

END SUBROUTINE ABSVALOVERLAP_WORKER

SUBROUTINE ABSVALOVERLAP_WORKER12(LUPRI,NBLEN,Nactbast,NBAST,&
     &  NMO1,NMO2,WEIGHT,GAOS,CMAT1,CMAT2,ABSVALOVERLAP,TMP,TMP2)
IMPLICIT NONE
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the number of gridpoints
INTEGER,intent(in) :: NBLEN
!> Number of active basis functions
INTEGER,intent(in) :: Nactbast
!> Number of basis functions
INTEGER,intent(in) :: Nbast
!> Number of MO coef
INTEGER,intent(in) :: NMO1,NMO2
!> coefficient (combination of functional derivative)
REAL(REALK),intent(in) :: WEIGHT(NBLEN)
!> gaussian atomic orbitals
REAL(REALK),intent(in) :: GAOS(NBLEN,NACTBAST)
!> The active AO to MO coeff matrix
REAL(REALK),intent(inout) :: CMAT1(NACTBAST,NMO1)
!> The active AO to MO coeff matrix
REAL(REALK),intent(inout) :: CMAT2(NACTBAST,NMO2)
!> The absolute valued Overlap matrix in MO
REAL(REALK),intent(inout) :: ABSVALOVERLAP(NMO1,NMO2)
! TEMPORARY MEM FROM WORK
REAL(REALK),intent(inout) :: TMP(NBLEN,NMO1) 
REAL(REALK),intent(inout) :: TMP2(NBLEN,NMO2)
!
Integer :: J,K,IBL,IACT
!MO Values on grid point 
!TMP(NBLEN,NMO1) = GAO(NBLEN,NACTBAST)*CMAT1(NACTBAST,NMO1)
CALL DGEMM('N','N',NBLEN,NMO1,NACTBAST,1.0E0_realk,&
     &             GAOS,NBLEN,CMAT1,NACTBAST,0.0E0_realk,&
     &             TMP,NBLEN)

! First half-contraction of MO's with grid point weight
DO J=1,NMO1
   DO K=1, NBLEN
      TMP(K,J) = WEIGHT(K)*ABS(TMP(K,J))
   ENDDO
ENDDO
!Second MO Values on grid point 
!TMP2(NBLEN,NMO2) = GAO(NBLEN,NACTBAST)*CMAT2(NACTBAST,NMO2)
CALL DGEMM('N','N',NBLEN,NMO2,NACTBAST,1.0E0_realk,&
     &             GAOS,NBLEN,CMAT2,NACTBAST,0.0E0_realk,&
     &             TMP2,NBLEN)
!take absolute value
DO J=1,NMO2
   DO K=1, NBLEN
      TMP2(K,J) = ABS(TMP2(K,J))
   ENDDO
ENDDO

!  Second half-contraction of MO's with grid point weight
!ABSVALOVERLAP(NMO1,NMO2) =+ op(TMP)(NBLEN,NMO1)*TMP2(NBLEN,NMO2)
CALL DGEMM('T','N',NMO1,NMO2,NBLEN,1.0E0_realk,&
     &             TMP,NBLEN,TMP2,NBLEN,1.0E0_realk,&
     &             ABSVALOVERLAP,NMO1)
END SUBROUTINE ABSVALOVERLAP_WORKER12

!> \brief determine the number of active orbitals
!> \author T. Kjaergaard
!> \date 2010
!>
!>  orbitals that contribute to the current gridpoints
SUBROUTINE DETERMINE_ABSVAL_NACTIVEORB(NactBAS,NSHELLBLOCKS,BLOCKS,MAXNSHELL)
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

END SUBROUTINE DETERMINE_ABSVAL_NACTIVEORB

!> \brief construct active Dmat
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE CONSTRUCT_ACTIVE_CMAT(NSHELLBLOCKS,MAXNSHELL,BLOCKS,CMAT,NBAST,NMO,ACTIVE_CMAT,NactBAS)
IMPLICIT NONE 
INTEGER,intent(in)     :: NSHELLBLOCKS,MAXNSHELL,NBAST,NactBAS,NMO
INTEGER,intent(in)     :: BLOCKS(2,MAXNSHELL)
REAL(REALK),intent(in) :: CMAT(NBAST,NMO)
REAL(REALK),intent(inout):: ACTIVE_CMAT(NactBAS,NMO)
!
INTEGER               :: IBL,IORB,J,IBASIS
IF(NBAST.EQ.NactBAS)THEN
   DO J = 1, NMO
      DO IBASIS = 1,NBAST
         ACTIVE_CMAT(IBASIS,J) = CMAT(IBASIS,J)
      ENDDO
   ENDDO
ELSE
   DO J = 1, NMO
      IBASIS = 0
      DO IBL = 1, NSHELLBLOCKS
         DO IORB = BLOCKS(1,IBL),BLOCKS(2,IBL) !ORBITALINDEX
            IBASIS = IBASIS + 1
            ACTIVE_CMAT(IBASIS,J) = CMAT(IORB,J)
         ENDDO
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE CONSTRUCT_ACTIVE_CMAT

!> \brief construct active Dmat
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE CONSTRUCT_ACTIVE_CMAT12(NSHELLBLOCKS,MAXNSHELL,BLOCKS,&
     & CMAT1,CMAT2,NBAST,NMO1,NMO2,ACTIVE_CMAT1,ACTIVE_CMAT2,NactBAS)
IMPLICIT NONE 
INTEGER,intent(in)     :: NSHELLBLOCKS,MAXNSHELL,NBAST,NactBAS,NMO1,NMO2
INTEGER,intent(in)     :: BLOCKS(2,MAXNSHELL)
REAL(REALK),intent(in) :: CMAT1(NBAST,NMO1)
REAL(REALK),intent(in) :: CMAT2(NBAST,NMO2)
REAL(REALK),intent(inout):: ACTIVE_CMAT1(NactBAS,NMO1)
REAL(REALK),intent(inout):: ACTIVE_CMAT2(NactBAS,NMO2)
!
INTEGER               :: IBL,IORB,J,IBASIS
IF(NBAST.EQ.NactBAS)THEN
   DO J = 1, NMO1
      DO IBASIS = 1,NBAST
         ACTIVE_CMAT1(IBASIS,J) = CMAT1(IBASIS,J)
      ENDDO
   ENDDO
   DO J = 1, NMO2
      DO IBASIS = 1,NBAST
         ACTIVE_CMAT2(IBASIS,J) = CMAT2(IBASIS,J)
      ENDDO
   ENDDO
ELSE
   DO J = 1, NMO1
      IBASIS = 0
      DO IBL = 1, NSHELLBLOCKS
         DO IORB = BLOCKS(1,IBL),BLOCKS(2,IBL) !ORBITALINDEX
            IBASIS = IBASIS + 1
            ACTIVE_CMAT1(IBASIS,J) = CMAT1(IORB,J)
         ENDDO
      ENDDO
   ENDDO
   DO J = 1, NMO2
      IBASIS = 0
      DO IBL = 1, NSHELLBLOCKS
         DO IORB = BLOCKS(1,IBL),BLOCKS(2,IBL) !ORBITALINDEX
            IBASIS = IBASIS + 1
            ACTIVE_CMAT2(IBASIS,J) = CMAT2(IORB,J)
         ENDDO
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE CONSTRUCT_ACTIVE_CMAT12

!> \brief evaluates the gaussian atomic orbitals
!> \author T. Kjaergaard
!> \date 2010
!>
!>  RETURNS: GSO: evaluated orbitals for a batch of grid points.
!>     GSO(:,:) contains orbital values.
!>     After requested geometric derivatives, London related derivatives
!>     are placed.
!>  REWRITTEN BY T.KJAERGAARD ORIGINALLY BY  T. Helgaker sep 99, P. Salek 03
!>
SUBROUTINE II_ABSVAL_BLGETSOS(LUPRI,NVCLEN,GAO,COOR,NBLCNT,IBLCKS,&
     &        NactBAS,DFTHRI,IPRINT,MMM,MAXANGMOM,SIZE,&
     &        SPHMAT,SPINDEX,MAXNSHELL,NSTART,SHELLANGMOM,CENT,&
     &        PRIEXPSTART,CC,ushells,CCSTART,CCINDEX,MXPRIM,PRIEXP,&
     &        LVALUE,MVALUE,NVALUE)
IMPLICIT NONE
INTEGER,intent(in) :: LUPRI,NVCLEN,NBLCNT,NactBAS,IPRINT,SIZE,MAXANGMOM,MAXNSHELL
INTEGER,intent(in) :: IBLCKS(2,MAXNSHELL),SPINDEX(MAXANGMOM),NSTART(MAXNSHELL),SHELLANGMOM(MAXNSHELL),ushells
INTEGER,intent(in) :: PRIEXPSTART(MAXNSHELL),MXPRIM,CCSTART(ushells),MMM
integer,intent(in) :: LVALUE(MMM,MAXANGMOM),MVALUE(MMM,MAXANGMOM)
integer,intent(in) :: NVALUE(MMM,MAXANGMOM)
INTEGER,intent(in)     :: CCINDEX(MAXNSHELL)
REAL(REALK),intent(in) :: PRIEXP(MXPRIM)
real(realk),intent(in) :: CENT(3,MAXNSHELL), COOR(3,NVCLEN)
REAL(REALK),intent(inout) :: GAO(NVCLEN,NactBAS)
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
IADR = 0
DO IBL = 1, NBLCNT
   DO ISHELA = IBLCKS(1,IBL),IBLCKS(2,IBL) !active Orbitals
      SHELLANGMOMA = SHELLANGMOM(ISHELA)
      KHKTA  = 2*(SHELLANGMOMA-1)+1  
      KCKTA  = SHELLANGMOMA*(SHELLANGMOMA+1)/2  
      SPVAL= KCKTA*KCKTA
      JSTA = PRIEXPSTART(ISHELA)
      CENX = CENT(1,ISHELA)!+ ORIG(1)
      CENY = CENT(2,ISHELA)!+ ORIG(2)
      CENZ = CENT(3,ISHELA)!+ ORIG(3)
      DO I=1,NVCLEN !gridpoints
         PA(1,i) = COOR(1,i)-CENX 
         PA(2,i) = COOR(2,i)-CENY
         PA(3,i) = COOR(3,i)-CENZ
      END DO
      DO I=1,NVCLEN !gridpoints
         PA2(i) = PA(1,i)*PA(1,i)+PA(2,i)*PA(2,i)+PA(3,i)*PA(3,i)
      END DO
      CALL II_ABSVAL_BLGETGAO(LUPRI,NVCLEN,NactBAS,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
           & SPHMAT(SPINDEX(SHELLANGMOMA):SPINDEX(SHELLANGMOMA)+KHKTA*KCKTA-1),&
           & GAO,IADR,PA,PA2,DFTHRI,CC(CCINDEX(ISHELA)),&
           & CCSTART(CCINDEX(ISHELA)),PRIEXP,LVALUE(1:SPVAL,SHELLANGMOMA),&
           & MVALUE(1:SPVAL,SHELLANGMOMA),NVALUE(1:SPVAL,SHELLANGMOMA),SPVAL,.TRUE.)
      IADR = IADR + KHKTA !UPDATE ACTIVEORBITALINDEX
   END DO
END DO

call mem_dft_dealloc(PA)
call mem_dft_dealloc(PA2)

END SUBROUTINE II_ABSVAL_BLGETSOS

!> \brief evaluates the pure GAO(gaussian atomic orbitals) , so no derivatives
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE II_ABSVAL_BLGETGAO(LUPRI,NVCLEN,NactBAS,SHELLANGMOMA,KHKTA,KCKTA,MXPRIM,JSTA,&
                     & CSP,GAO,IADR,PA,PA2,DFTHRI,CC,CCSTART,PRIEXP,LVALUE,MVALUE,&
                     &NVALUE,SPVAL,SPHRA)
IMPLICIT NONE
INTEGER,intent(in)        :: NVCLEN,NactBAS,MXPRIM,KHKTA,JSTA,KCKTA,LUPRI,SHELLANGMOMA
INTEGER,intent(in)        :: SPVAL,CCSTART,IADR,LVALUE(SPVAL),MVALUE(SPVAL)
INTEGER,intent(in)        :: NVALUE(SPVAL)
REAL(REALK),PARAMETER     :: D0 = 0.0E0_realk, D1 = 1.0E0_realk
REAL(REALK),intent(inout) :: GAO(NVCLEN,NactBAS)
REAL(REALK),intent(in)    :: PA(3,NVCLEN),PA2(NVCLEN),DFTHRI
REAL(REALK),intent(in)    :: CSP(KHKTA,KCKTA),PRIEXP(MXPRIM)
TYPE(LSMATRIX),intent(in) :: CC
!
REAL(REALK)    :: GAZ,GAXX,GAX,GAY,GAYY,GMAX,PA_1,PA_2,PA_3
REAL(REALK)    :: PRICCFVAL,PRIEXPVAL
INTEGER        :: I,J,TEMPI,LVALJ, MVALJ, NVALJ, LVALI, MVALI, NVALI,K
INTEGER        :: ISTART1,ISTART2,ISTART3,ISTART4,ISTART5,ISTART6,ISTART
LOGICAl        :: SPHRA
REAL(REALK),pointer  :: GA(:),CINT(:)

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
         GAO(K,ISTART1) = GA(K)
      END DO
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1),NVCLEN)
   END IF
ELSEIF (SHELLANGMOMA .EQ. 2) THEN !p orbitals
   IF (GMAX .GT. DFTHRI) THEN
      DO K = 1, NVCLEN
         GAO(K,ISTART1) = PA(1,K)*GA(K)
         GAO(K,ISTART2) = PA(2,K)*GA(K)
         GAO(K,ISTART3) = PA(3,K)*GA(K)
      END DO
   ELSE
      CALL LS_DZERO(GAO(1,ISTART1),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART2),NVCLEN)
      CALL LS_DZERO(GAO(1,ISTART3),NVCLEN)
   END IF
ELSEIF (SHELLANGMOMA .EQ. 3) THEN !d orbitals
   IF (SPHRA) THEN
      IF (GMAX .GT. DFTHRI) THEN
         DO K = 1, NVCLEN
            PA_1 = PA(1,K)
            PA_2 = PA(2,K)
            PA_3 = PA(3,K)
            GAX  = PA_1*GA(K)
            GAY  = PA_2*GA(K)
            GAZ  = PA_3*GA(K)
            GAXX = PA_1*GAX
            GAYY = PA_2*GAY
            GAO(K,ISTART1) = CSP(1,2)*PA_2*GAX
            GAO(K,ISTART2) = CSP(2,5)*PA_2*GAZ
            GAO(K,ISTART3) = CSP(3,1)*GAXX + CSP(3,4)*GAYY&
            &                  + CSP(3,6)*PA_3*GAZ
            GAO(K,ISTART4) = CSP(4,3)*PA_1*GAZ
            GAO(K,ISTART5) = CSP(5,1)*GAXX + CSP(5,4)*GAYY
         END DO
      ELSE 
         CALL LS_DZERO(GAO(1,ISTART1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART2),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART3),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART4),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART5),NVCLEN)
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
            GAO(K,ISTART1) = PA_1*GAX
            GAO(K,ISTART2) = PA_2*GAX
            GAO(K,ISTART3) = PA_3*GAX
            GAO(K,ISTART4) = PA_2*GAY
            GAO(K,ISTART5) = PA_3*GAY 
            GAO(K,ISTART6) = PA_3*GAZ
         END DO
      ELSE
         CALL LS_DZERO(GAO(1,ISTART1),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART2),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART3),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART4),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART5),NVCLEN)
         CALL LS_DZERO(GAO(1,ISTART6),NVCLEN)
      END IF
   END IF
ELSE !higher than d orbitals
   IF (SPHRA) THEN
      call mem_dft_alloc(CINT,NVCLEN)
      DO I = 1, KHKTA
         CALL LS_DZERO(GAO(1,ISTART+I),NVCLEN)
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
            !              do a dgemm here?
            DO I = 1, KHKTA
               TEMPI=ISTART+I
               DO K = 1, NVCLEN
                  GAO(K,TEMPI) = GAO(K,TEMPI) + CSP(I,J)*CINT(K)
               END DO
            END DO
         END DO
      END IF
      call mem_dft_dealloc(CINT)
   ELSE
      IF (GMAX .GT. DFTHRI) THEN
         DO I = 1, KHKTA
            LVALI = LVALUE(I)
            MVALI = MVALUE(I)
            NVALI = NVALUE(I)
            TEMPI=ISTART+I
            DO K = 1, NVCLEN
               GAO(K,TEMPI) = (PA(1,K)**LVALI)*(PA(2,K)**MVALI)&
                    &                 *(PA(3,K)**NVALI)*GA(K)
            END DO
         END DO
      END IF
   END IF
END IF
call mem_dft_dealloc(GA)

END SUBROUTINE II_ABSVAL_BLGETGAO

!> \brief make blocks that contain the orbitals which contribute to each gridbatch
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE SHELL_TO_ABSVAL_ORB(LUPRI,NSHELLBLOCKS,SHELBLOCK,ORBBLOCKS,NSTART,MAXNSHELL,NBAST)
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

END SUBROUTINE SHELL_TO_ABSVAL_ORB

!> \brief make blocks that contain the active orbitals which contribute to each gridbatch
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE SHELL_TO_ABSVAL_ACTORB(NSHELLBLOCKS,BLOCKS,MAXNSHELL)
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

END SUBROUTINE SHELL_TO_ABSVAL_ACTORB

!> \brief build spherical transformation matrices
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE Build_PRECALCULATED_SPHMAT_ABSVAL(LUPRI,MAXANGMOM,SIZE,SPHMAT,SPINDEX) 
use math_fun
IMPLICIT NONE
INTEGER,intent(in)        :: MAXANGMOM,SIZE,LUPRI
INTEGER,intent(inout)     :: SPINDEX(MAXANGMOM+1)
REAL(REALK),intent(inout) :: SPHMAT(SIZE)
!
INTEGER          :: L,I
Real(realk), parameter :: DM1 = -1.0E0_realk, DO = 0.0E0_realk, D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER  :: M1,MADR,MABS,V0, IOFF,nAngmom
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
END SUBROUTINE BUILD_PRECALCULATED_SPHMAT_ABSVAL
END MODULE IIABSVALINT
