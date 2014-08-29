!> @file
!> Module contains grid generation routines
!> \brief 
!> \author T. Kjaergaard
!> \date 2009 
MODULE gridgenerationboxmodule
use grid_memory_handling
!WARNING you must not add memory_handling, all memory goes through 
!grid_memory_handling  module so as to determine the memory used in this module.
use precision
use files
use dft_typetype
use lstiming
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_mod
#endif

SAVE
type gridboxtype
REAL(REALK) :: CELL_SIZE
REAL(REALK) :: maxR(3)
REAL(REALK) :: minR(3)
INTEGER :: npoints
INTEGER,pointer :: points(:)
TYPE(gridboxtype), POINTER :: next  
TYPE(gridboxtype), POINTER :: prev  
end type gridboxtype

TYPE bunchpoints
REAL(REALK) :: Center(3)
REAL(REALK) :: CELL_SIZE
integer,pointer :: atom_idx(:)
integer,pointer :: points(:)
INTEGER         :: npoints
END TYPE bunchpoints
CONTAINS
Subroutine BuildBoxes(gridBox,COOR2,GlobalmaxGridpoints,totalpoints,RSHEL,&
     & MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,atomcenterX,atomcenterY,atomcenterZ,&
     & NSTART,maxnbuflen,nBoxes,keypointer,iprint,lupri)
implicit none
INTEGER,intent(in)     ::  NBAST,MAXNSHELL,NATOMS,LUPRI,IPRINT
integer,intent(in)     :: GlobalmaxGridpoints,totalpoints,maxnbuflen
type(gridboxtype),pointer :: GridBox 
real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
INTEGER,intent(in)     ::  SHELL2ATOM(MAXNSHELL)
REAL(REALK),intent(in) ::  atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS),RSHEL(MAXNSHELL)
INTEGER,intent(in)     :: NSTART(MAXNSHELL)
INTEGER,intent(inout)  :: nBoxes
type(bunchpoints),pointer :: keypointer(:)
!
real(realk) :: minR(3),maxR(3),Lengthinv,center1(3),center2(3),length(3),split
integer :: maxDim,nBox1,nBox2,I,ix,NactBast,ip,MaxnBox
real(realk),parameter :: D05=0.5E0_realk,D8=8E0_realk
real(realk) :: TS,TE,Volumen1
real(realk) :: gridLengthX,gridLengthY,gridLengthZ
type(gridboxtype),pointer :: TMPBOX
integer,pointer :: IBLCKS(:,:),tmp(:)

CALL LSTIMER('START',TS,TE,LUPRI)

minR(1) = +1E+20_realk 
minR(2) = +1E+20_realk 
minR(3) = +1E+20_realk 
maxR(1) = -1E+20_realk 
maxR(2) = -1E+20_realk 
maxR(3) = -1E+20_realk 
DO I=1,totalpoints
   IF(COOR2(1,I).LT.minR(1))minR(1) = COOR2(1,I)
   IF(COOR2(1,I).GT.maxR(1))maxR(1) = COOR2(1,I)
   IF(COOR2(2,I).LT.minR(2))minR(2) = COOR2(2,I)
   IF(COOR2(2,I).GT.maxR(2))maxR(2) = COOR2(2,I)
   IF(COOR2(3,I).LT.minR(3))minR(3) = COOR2(3,I)
   IF(COOR2(3,I).GT.maxR(3))maxR(3) = COOR2(3,I)
ENDDO

!WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'X boxdimension [',minR(1),',',maxR(1),']'
!WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'Y boxdimension [',minR(2),',',maxR(2),']'
!WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'Z boxdimension [',minR(3),',',maxR(3),']'

!CENTER(1) = minR(1) + (maxR(1)-minR(1))*D05
!CENTER(2) = minR(1) + (maxR(1)-minR(1))*D05
!CENTER(3) = minR(1) + (maxR(1)-minR(1))*D05

maxDim = getmaxDim(minR,maxR)
Length(1) = (maxR(1)-minR(1))
Length(2) = (maxR(2)-minR(2))
Length(3) = (maxR(3)-minR(3))
!Lengthinv = 1/(maxR(maxDim)-minR(maxDim))
split = minR(maxDim)+Length(maxdim)*d05
nBox1=0
nBox2=0
DO I=1,totalpoints
   if(COOR2(maxDim,I).GT.split)THEN
      nBox1 = nBox1 +1
   else
      nBox2 = nBox2 +1
   endif
ENDDO

allocate(gridBox)
gridBox%npoints=nBox1
call mem_grid_alloc(gridBox%points,nBox1)
nullify(gridBox%next)
nullify(gridBox%prev)

allocate(tmpBox)
tmpBox%npoints=nBox2
call mem_grid_alloc(tmpBox%points,nBox2)
nullify(tmpBox%next)
nullify(tmpBox%prev)

gridBox%minR(1) = minR(1)
gridBox%minR(2) = minR(2)
gridBox%minR(3) = minR(3)
gridBox%maxR(1) = maxR(1)
gridBox%maxR(2) = maxR(2)
gridBox%maxR(3) = maxR(3)
gridBox%minR(maxDim) = minR(maxDim)+Length(maxdim)*d05

tmpBox%minR(1) = minR(1)
tmpBox%minR(2) = minR(2)
tmpBox%minR(3) = minR(3)
tmpBox%maxR(1) = maxR(1)
tmpBox%maxR(2) = maxR(2)
tmpBox%maxR(3) = maxR(3)
tmpBox%maxR(maxDim) = minR(maxDim)+Length(maxdim)*d05

gridLengthX = gridBox%maxR(1)-gridBox%minR(1)
gridLengthY = gridBox%maxR(2)-gridBox%minR(2)
gridLengthZ = gridBox%maxR(3)-gridBox%minR(3)
gridBox%CELL_SIZE = SQRT(gridLengthX*gridLengthX+gridLengthY*gridLengthY+gridLengthZ*gridLengthZ)*D05
gridBox%CELL_SIZE = SQRT(gridLengthX*gridLengthX+gridLengthY*gridLengthY+gridLengthZ*gridLengthZ)*D05
gridLengthX = tmpbox%maxR(1)-tmpbox%minR(1)
gridLengthY = tmpbox%maxR(2)-tmpbox%minR(2)
gridLengthZ = tmpbox%maxR(3)-tmpbox%minR(3)
tmpbox%CELL_SIZE = SQRT(gridLengthX*gridLengthX+gridLengthY*gridLengthY+gridLengthZ*gridLengthZ)*D05
tmpbox%CELL_SIZE = SQRT(gridLengthX*gridLengthX+gridLengthY*gridLengthY+gridLengthZ*gridLengthZ)*D05

nBox1=0
nBox2=0
DO I=1,totalpoints
   if(COOR2(maxDim,I).GT.split)THEN
      nBox1 = nBox1 + 1
      gridBox%points(nBox1) = I
   else
      nBox2 = nBox2 + 1
      tmpBox%points(nBox2) = I
   endif
ENDDO

gridBox%next => TmpBox
TmpBox%prev => gridBox
call mem_grid_alloc(IBLCKS,2,MAXNSHELL)
call BOX_NactOrbGETBLOCKS(RSHEL,MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,&
     & atomcenterX,atomcenterY,atomcenterZ,NactBast,NSTART,minR,maxR,LUPRI,IBLCKS)
call BOX_NactOrbGETBLOCKS(RSHEL,MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,&
     & atomcenterX,atomcenterY,atomcenterZ,NactBast,NSTART,GridBox%minR,GridBox%maxR,LUPRI,IBLCKS)
call BOX_NactOrbGETBLOCKS(RSHEL,MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,&
     & atomcenterX,atomcenterY,atomcenterZ,NactBast,NSTART,TmpBox%minR,TmpBox%maxR,LUPRI,IBLCKS)
Volumen1 = (maxR(1)-minR(1))*(maxR(2)-minR(2))*(maxR(3)-minR(3))

MaxnBox = MAX(gridBox%npoints,tmpBox%npoints)
call mem_grid_alloc(tmp,MaxnBox)
call split_box(GridBox%next,COOR2,GlobalmaxGridpoints,totalpoints,RSHEL,&
     & MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,atomcenterX,atomcenterY,atomcenterZ,&
     & NSTART,maxnbuflen,iprint,lupri,IBLCKS,tmp,MaxnBox)
call split_box(GridBox,COOR2,GlobalmaxGridpoints,totalpoints,RSHEL,&
     & MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,atomcenterX,atomcenterY,atomcenterZ,&
     & NSTART,maxnbuflen,iprint,lupri,IBLCKS,tmp,MaxnBox)
call mem_grid_dealloc(tmp)
call mem_grid_dealloc(IBLCKS)

nBoxes=1
call determine_nboxes(GridBox,nBoxes,lupri)

allocate(keypointer(nBoxes))
call build_keypointer(keypointer,nBoxes,GridBox,lupri)

!WRITE(LUPRI,*)'PRINT BOXES'
!call print_box(Gridbox,.TRUE.,lupri)

call free_box(Gridbox)

#ifdef VAR_MPI
IF (infpar%mynum.EQ.infpar%master) THEN
#endif
CALL LSTIMER('Boxify',TS,TE,LUPRI)
#ifdef VAR_MPI
ENDIF
#endif

end Subroutine BuildBoxes

function GetmaxDim(minR,maxR)
implicit none
real(realk) :: minR(3),maxR(3)
integer :: GetMaxdim
!
real(realk) :: range,dist
range = maxR(1) - minR(1)
Dist = range
GetmaxDim = 1
range = maxR(2) - minR(2)
IF(range.GT.Dist)THEN
   Dist = range
   GetmaxDim = 2
ENDIF
range = maxR(3) - minR(3)
IF(range.GT.Dist)GetmaxDim = 3
end function GetmaxDim

subroutine build_keypointer(keypointer,nBoxes,GridBox,lupri)
implicit none
integer,intent(in) :: nBoxes
type(bunchpoints)  :: keypointer(nBoxes)
type(gridboxtype),pointer :: gridbox
integer,intent(in) :: lupri
!
real(realk),parameter :: D05=0.5E0_realk
integer :: I,IP
I=1
keypointer(I)%npoints = GridBox%npoints
keypointer(I)%CELL_SIZE = GridBox%CELL_SIZE
keypointer(I)%Center(1) = GridBox%minR(1)+(GridBox%maxR(1)-GridBox%minR(1))*d05
keypointer(I)%Center(2) = GridBox%minR(2)+(GridBox%maxR(2)-GridBox%minR(2))*d05
keypointer(I)%Center(3) = GridBox%minR(3)+(GridBox%maxR(3)-GridBox%minR(3))*d05
!cannot use mem_alloc due to OpenMP mem bookkeeping
allocate(keypointer(I)%points(GridBox%npoints))
do IP=1,GridBox%npoints
   keypointer(I)%points(IP) = gridbox%points(IP)      
enddo
if(associated(gridbox%next))then
   call add_keypointer(keypointer,nBoxes,GridBox%next,I,lupri)
endif
IF(I.NE.nBoxes)call lsquit('Internal grid error. ',lupri)

end subroutine build_keypointer

recursive subroutine add_keypointer(keypointer,nBoxes,GridBox,I,lupri)
implicit none
integer,intent(inout) :: I
integer,intent(in) :: nBoxes
type(bunchpoints) :: keypointer(nBoxes)
type(gridboxtype),pointer :: gridbox
integer,intent(in) :: lupri
!
real(realk),parameter :: D05=0.5E0_realk
integer :: IP

if(GridBox%npoints.GT. 0)then
   I=I+1
   keypointer(I)%npoints = GridBox%npoints
   keypointer(I)%CELL_SIZE = GridBox%CELL_SIZE
   keypointer(I)%Center(1) = GridBox%minR(1)+(GridBox%maxR(1)-GridBox%minR(1))*d05
   keypointer(I)%Center(2) = GridBox%minR(2)+(GridBox%maxR(2)-GridBox%minR(2))*d05
   keypointer(I)%Center(3) = GridBox%minR(3)+(GridBox%maxR(3)-GridBox%minR(3))*d05
   !cannot use mem_alloc due to OpenMP mem bookkeeping  
   allocate(keypointer(I)%points(GridBox%npoints))
   do IP=1,GridBox%npoints
      keypointer(I)%points(IP) = gridbox%points(IP)      
   enddo
endif
if(associated(gridbox%next))then
   call add_keypointer(keypointer,nBoxes,GridBox%next,I,lupri)
endif

end subroutine add_keypointer

recursive subroutine determine_nboxes(gridbox,nBoxes,lupri)
implicit none
type(gridboxtype),pointer :: gridbox
integer,intent(inout) :: nBoxes
type(gridboxtype),pointer :: tmp
integer :: lupri
!
if(associated(gridbox%next))then
   if(gridbox%next%npoints.GT. 0)nBoxes =  nBoxes + 1
   call determine_nboxes(Gridbox%next,nBoxes,lupri)
endif
end subroutine determine_nboxes

recursive subroutine free_box(gridbox)
implicit none
type(gridboxtype),pointer :: gridbox

if(associated(gridbox%next))then
   call free_box(GridBox%next)
   call mem_grid_dealloc(GridBox%points)
   deallocate(GridBox)
else
   call mem_grid_dealloc(GridBox%points)
   deallocate(GridBox)
endif

end subroutine free_box

subroutine print_box2(gridbox,COOR2,GlobalmaxGridpoints,lupri)
implicit none
type(gridboxtype),pointer :: gridbox
integer,intent(in)     :: GlobalmaxGridpoints
real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
integer :: lupri
!
integer :: IP,I
if(gridbox%npoints.GT. 0)then
   WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'X Gridboxdimension [',GridBox%minR(1),',',GridBox%maxR(1),']'
   WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'Y Gridboxdimension [',GridBox%minR(2),',',GridBox%maxR(2),']'
   WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'Z Gridboxdimension [',GridBox%minR(3),',',GridBox%maxR(3),']'
   write(lupri,'(A,I9)')'GridBox%npoints',gridbox%npoints

   DO Ip=1,GridBox%npoints
      I = GridBox%points(Ip)
      WRITE(lupri,'(A,F12.5,A,F12.5,A,F12.5,A)')'COOR=(',COOR2(1,I),',',COOR2(2,I),',',COOR2(3,I),')'
   enddo
endif
end subroutine print_box2

recursive subroutine print_box(gridbox,full,lupri)
implicit none
type(gridboxtype),pointer :: gridbox
logical :: full
integer :: lupri
!
if(gridbox%npoints.GT. 0)then
   WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'X Gridboxdimension [',GridBox%minR(1),',',GridBox%maxR(1),']'
   WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'Y Gridboxdimension [',GridBox%minR(2),',',GridBox%maxR(2),']'
   WRITE(lupri,'(A,ES16.8,A,ES16.8,A)')'Z Gridboxdimension [',GridBox%minR(3),',',GridBox%maxR(3),']'
   write(lupri,'(A,F18.6)')'GridBox%CELL_SIZE               ',gridbox%CELL_SIZE
   write(lupri,'(A,I9)')'GridBox%npoints',gridbox%npoints
endif
if(associated(GridBox%points))then
   WRITE(lupri,'(A)')'GridBox%points(1:GridBox%npoints)'
   WRITE(lupri,'(5X,5I9/,(5X,5I9))')GridBox%points(1:GridBox%npoints)
endif
IF(full.AND.associated(GridBox%next))THEN
   call print_box(Gridbox%next,full,lupri)
ENDIF
end subroutine print_box

SUBROUTINE BOX_NactOrbGETBLOCKS(RSHEL,MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,X,Y,Z,NactBast,NSTART,minR,maxR,LUPRI,IBLCKS)
use precision
IMPLICIT NONE
INTEGER,intent(in)    ::  NBAST,MAXNSHELL,NATOMS,LUPRI
INTEGER,intent(in)    ::  SHELL2ATOM(MAXNSHELL),NSTART(MAXNSHELL)
REAL(REALK),intent(in)::  X(NATOMS),Y(NATOMS),Z(NATOMS),RSHEL(MAXNSHELL)
REAL(REALK),intent(in) ::  minR(3),maxR(3)
INTEGER,intent(inout)  ::  NactBast
INTEGER,intent(inout) ::  IBLCKS(2,MAXNSHELL)
!
INTEGER     ::  NBLCNT
INTEGER     ::  ISHLEN,IPREV,ICENT,ISHELA,I,J
REAL(REALK) ::  RSHEL1,center(3),LengthX,LengthY,LengthZ,CELLDG,PX,PY,PZ,DST
INTEGER     :: IBL,IORB
LOGICAL     :: XVERIFY,YVERIFY,ZVERIFY
real(realk),parameter :: D05=0.5E0_realk
!each box can save which ISHELAs that contribute ... 
!so that this does not need to be reevaluated.
NBLCNT = 0
IPREV = -1111
LengthX = maxR(1)-minR(1)
LengthY = maxR(2)-minR(2)
LengthZ = maxR(3)-minR(3)
Center(1) = minR(1)+LengthX*d05
Center(2) = minR(2)+LengthY*d05
Center(3) = minR(3)+LengthZ*d05
CELLDG = SQRT(LengthX*LengthX+LengthY*LengthY+LengthZ*LengthZ)*d05
DO ISHELA=1,MAXNSHELL
   ICENT = SHELL2ATOM(ISHELA)
   PX = center(1)-X(ICENT)
   PY = center(2)-Y(ICENT)
   PZ = center(3)-Z(ICENT)
   DST = SQRT(PX*PX + PY*PY + PZ*PZ)
   IF(DST.LE.RSHEL(ISHELA)+CELLDG) THEN
      !       accepted...
      IF(ISHELA.NE.IPREV+1) THEN
         NBLCNT = NBLCNT + 1
         IBLCKS(1,NBLCNT) = ISHELA
      END IF
      IPREV = ISHELA
      IBLCKS(2,NBLCNT) = ISHELA
   END IF
END DO 

IF(NBLCNT.GT.0)THEN
!!$   DO I = 1,NBLCNT
!!$      ORBBLOCKS(1,I) = NSTART(IBLCKS(1,I))+1
!!$   ENDDO
!!$   
!!$   DO I = 1,NBLCNT-1
!!$      ORBBLOCKS(2,I) = NSTART(IBLCKS(2,I)+1)
!!$   ENDDO
!!$   I=NBLCNT
!!$   IF(IBLCKS(2,I) .LT. MAXNSHELL)THEN
!!$      ORBBLOCKS(2,I) = NSTART(IBLCKS(2,I)+1)
!!$   ELSE
!!$      ORBBLOCKS(2,I) = NBAST
!!$   ENDIF
!!$
!!$   NactBAST = 0
!!$   DO IBL = 1, NBLCNT
!!$      NactBAST = NactBAST + ORBBLOCKS(2,IBL) - ORBBLOCKS(1,IBL) + 1
!!$   ENDDO
   NactBAST = 0
   DO I = 1,NBLCNT-1
      NactBAST = NactBAST + NSTART(IBLCKS(2,I)+1) - NSTART(IBLCKS(1,I))
   ENDDO   
   I=NBLCNT
   IF(IBLCKS(2,I) .LT. MAXNSHELL)THEN
      NactBAST = NactBAST + NSTART(IBLCKS(2,I)+1) - NSTART(IBLCKS(1,I))
   ELSE
      NactBAST = NactBAST + NBAST - NSTART(IBLCKS(1,I))
   ENDIF
ELSE
   NactBAST = 0
ENDIF
END SUBROUTINE BOX_NACTORBGETBLOCKS

recursive subroutine split_box(gridBox,COOR2,GlobalmaxGridpoints,totalpoints,&
     & RSHEL,MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,atomcenterX,atomcenterY,&
     & atomcenterZ,NSTART,maxnbuflen,iprint,lupri,IBLCKS,tmp,maxnBox)
implicit none
INTEGER,intent(in)     ::  NBAST,MAXNSHELL,NATOMS,LUPRI,IPRINT,MaxnBox
integer,intent(in)     :: GlobalmaxGridpoints,totalpoints,maxnbuflen
type(gridboxtype),pointer :: GridBox 
real(realk),intent(in) :: COOR2(3,GlobalmaxGridpoints)
INTEGER,intent(in)     ::  SHELL2ATOM(MAXNSHELL)
REAL(REALK),intent(in) ::  atomcenterX(NATOMS),atomcenterY(NATOMS),atomcenterZ(NATOMS),RSHEL(MAXNSHELL)
INTEGER,intent(in)     :: NSTART(MAXNSHELL)
INTEGER,intent(inout) :: IBLCKS(2,MAXNSHELL)
integer,intent(inout) :: tmp(maxnBox)
!
real(realk) :: Lengthinv,center1(3),center2(3),length(3),split,minR(3),maxR(3)
integer :: maxDim,nBox1,nBox2,I,ix,NactBast,ip
real(realk),parameter :: D05=0.5E0_realk,D8=8E0_realk
real(realk) :: Volumen1,gridLengthX,gridLengthY,gridLengthZ
type(gridboxtype),pointer :: TMPBOX

!FIXME we prefer boxes with GridBox%npoints and NactBast a multiple of 4 or better 8 - as much as possible. 

IF(GridBox%npoints.LT.maxnbuflen)THEN
   !determine number of active orbitals in this box 
   !if it is too big - which means that the memory 
   !requirements are too big we devide again
   call BOX_NactOrbGETBLOCKS(RSHEL,MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,&
        & atomcenterX,atomcenterY,atomcenterZ,NactBast,NSTART,&
        & GridBox%minR,GridBox%maxR,LUPRI,IBLCKS)
   IF(NactBAST.EQ.0)THEN
      RETURN
   ENDIF
   IF(MAX(MXBLLEN,GridBox%npoints,NactBast)*NactBast.LT.BoxMemRequirement)THEN
!      print*,'BoxMemRequirement',BoxMemRequirement,'NBAST',NBAST
!      print*,'MXBLLEN*NactBast',MXBLLEN*NactBast,'NactBast',NactBast
!      print*,'MAX(MXBLLEN,GridBox%npoints,NactBast)*NactBast',MAX(MXBLLEN,GridBox%npoints,NactBast)*NactBast
      RETURN
!   ELSE
!      continue to divide
!      print*,'BoxMemRequirement',BoxMemRequirement,'NBAST',NBAST,100000000
!      print*,'MXBLLEN*NactBast',MXBLLEN*NactBast,'NactBast',NactBast
   ENDIF
ENDIF
minR(1) = gridBox%minR(1)
minR(2) = gridBox%minR(2)
minR(3) = gridBox%minR(3)
maxR(1) = gridBox%maxR(1)
maxR(2) = gridBox%maxR(2)
maxR(3) = gridBox%maxR(3)
maxDim = getmaxDim(minR,maxR)
Length(1) = (maxR(1)-minR(1))
Length(2) = (maxR(2)-minR(2))
Length(3) = (maxR(3)-minR(3))
!Lengthinv = 1/(maxR(maxDim)-minR(maxDim))
split = minR(maxDim)+Length(maxdim)*d05
nBox1=0
nBox2=0
DO Ip=1,GridBox%npoints
   I = GridBox%points(Ip)
   if(COOR2(maxDim,I).GT.split)THEN
      nBox1 = nBox1 +1
   else
      nBox2 = nBox2 +1
   endif
ENDDO

allocate(tmpBox)
tmpBox%npoints=nBox2
call mem_grid_alloc(tmpBox%points,nBox2)
nullify(tmpBox%next)
nullify(tmpBox%prev)

gridBox%minR(1) = minR(1)
gridBox%minR(2) = minR(2)
gridBox%minR(3) = minR(3)
gridBox%maxR(1) = maxR(1)
gridBox%maxR(2) = maxR(2)
gridBox%maxR(3) = maxR(3)
gridBox%minR(maxDim) = minR(maxDim)+Length(maxdim)*d05

tmpBox%minR(1) = minR(1)
tmpBox%minR(2) = minR(2)
tmpBox%minR(3) = minR(3)
tmpBox%maxR(1) = maxR(1)
tmpBox%maxR(2) = maxR(2)
tmpBox%maxR(3) = maxR(3)
tmpBox%maxR(maxDim) = minR(maxDim)+Length(maxdim)*d05

gridLengthX = gridBox%maxR(1)-gridBox%minR(1)
gridLengthY = gridBox%maxR(2)-gridBox%minR(2)
gridLengthZ = gridBox%maxR(3)-gridBox%minR(3)
gridBox%CELL_SIZE = SQRT(gridLengthX*gridLengthX+gridLengthY*gridLengthY+gridLengthZ*gridLengthZ)*D05
gridBox%CELL_SIZE = SQRT(gridLengthX*gridLengthX+gridLengthY*gridLengthY+gridLengthZ*gridLengthZ)*D05
gridLengthX = tmpbox%maxR(1)-tmpbox%minR(1)
gridLengthY = tmpbox%maxR(2)-tmpbox%minR(2)
gridLengthZ = tmpbox%maxR(3)-tmpbox%minR(3)
tmpbox%CELL_SIZE = SQRT(gridLengthX*gridLengthX+gridLengthY*gridLengthY+gridLengthZ*gridLengthZ)*D05
tmpbox%CELL_SIZE = SQRT(gridLengthX*gridLengthX+gridLengthY*gridLengthY+gridLengthZ*gridLengthZ)*D05

nBox1=0
nBox2=0
DO Ip=1,GridBox%npoints
   I = GridBox%points(Ip)
   if(COOR2(maxDim,I).GT.split)THEN
      nBox1 = nBox1 +1
      tmp(nBox1)=I
   else
      nBox2 = nBox2 +1
      tmpBox%points(nBox2) = I
   endif
ENDDO

call mem_grid_dealloc(GridBox%points)
call mem_grid_alloc(GridBox%points,nBox1)
GridBox%npoints=nBox1
do I=1,nBox1
   gridBox%points(I) = tmp(I)
enddo

if(associated(gridBox%next))then
   !not the last box in the chain
   TmpBox%next => gridBox%next
endif
GridBox%next => TmpBox
TmpBox%prev => GridBox

call split_box(GridBox%next,COOR2,GlobalmaxGridpoints,totalpoints,RSHEL,&
     & MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,atomcenterX,atomcenterY,atomcenterZ,&
     & NSTART,maxnbuflen,iprint,lupri,IBLCKS,tmp,MaxnBox)

call split_box(GridBox,COOR2,GlobalmaxGridpoints,totalpoints,RSHEL,&
     & MAXNSHELL,SHELL2ATOM,NBAST,NATOMS,atomcenterX,atomcenterY,atomcenterZ,&
     & NSTART,maxnbuflen,iprint,lupri,IBLCKS,tmp,MaxnBox)

end subroutine split_box

subroutine AddAtomicidxToBox(keypointer,nkey,ATOMIDX,GlobalmaxGridpoints)
implicit none
integer,intent(in) :: nkey,GlobalmaxGridpoints
integer,intent(in) :: ATOMIDX(GlobalmaxGridpoints)   
type(bunchpoints),intent(inout) :: keypointer(nkey)
!
integer :: Ikey,IP,I,npoints2

do Ikey=1,nkey
   npoints2=keypointer(Ikey)%npoints
   !cannot use mem_alloc due to OpenMP mem bookkeeping
   allocate(keypointer(Ikey)%atom_idx(npoints2))
   do IP = 1,npoints2
      I = keypointer(Ikey)%points(IP)
      keypointer(Ikey)%atom_idx(IP) = ATOMIDX(I)
   enddo
enddo

end subroutine AddAtomicidxToBox

end MODULE gridgenerationboxmodule
