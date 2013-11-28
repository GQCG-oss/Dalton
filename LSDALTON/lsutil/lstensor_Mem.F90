!> @file
!> Contains the orbital and overlap types, and subroutines to set these up
MODULE lstensorMem
use precision
use memory_handling

type lstmemchainitem
integer :: index
!sizes
integer(kind=long) :: lstmem_INTWORKLENGTH
integer(kind=long) :: lstmem_INTSWORKLENGTH
integer(kind=long) :: lstmem_REALKWORKLENGTH
!last unused element 
integer(kind=long) :: LSTMEM_LINTWORK
integer(kind=long) :: LSTMEM_LINTSWORK
integer(kind=long) :: LSTMEM_LREALKWORK
integer,pointer :: LSTMEM_INTWORK(:)
integer(kind=short),pointer :: LSTMEM_INTSWORK(:)
real(realk),pointer :: LSTMEM_REALKWORK(:)
type(lstmemchainitem),pointer :: next
end type lstmemchainitem

TYPE(lstmemchainitem),pointer :: lstmem_chain
TYPE(lstmemchainitem),pointer :: current
TYPE(lstmemchainitem),pointer :: last
integer                       :: nlstmem
integer :: lstindex
!$OMP THREADPRIVATE(lstmem_chain,nlstmem,current,lstindex,last)

INTERFACE mem_LSTpointer_alloc
  MODULE PROCEDURE real_setLSTpointer_alloc, int_setLSTpointer_alloc, &
       & ints_setLSTpointer_alloc
END INTERFACE
CONTAINS
subroutine lstmem_init()
implicit none
lstindex= 0
nlstmem = 0
nullify(lstmem_chain)
allocate(lstmem_chain)
nullify(lstmem_chain%LSTMEM_INTWORK)
nullify(lstmem_chain%LSTMEM_INTSWORK)
nullify(lstmem_chain%LSTMEM_REALKWORK)
allocate(lstmem_chain%next)
nullify(lstmem_chain%next%LSTMEM_INTWORK)
nullify(lstmem_chain%next%LSTMEM_INTSWORK)
nullify(lstmem_chain%next%LSTMEM_REALKWORK)
nullify(lstmem_chain%next%next)
nullify(current)
nullify(last)
lstmem_chain%index = 0
end subroutine lstmem_init

!subroutine mpicopy_lstensorMem(master)
!  use lsmpi_type 
!implicit none
!integer :: master
!call LS_MPI_BUFFER(current%LSTMEM_INTWORK,current%lstmem_INTWORKLENGTH,master)
!call LS_MPI_BUFFER(current%LSTMEM_INTSWORK,current%lstmem_INTSWORKLENGTH,master)
!call LS_MPI_BUFFER(current%LSTMEM_REALKWORK,current%lstmem_REALKWORKLENGTH,master)
!end subroutine mpicopy_lstensorMem

subroutine lstmem_free()
implicit none
IF(associated(lstmem_chain%next))THEN
   call lstmem_rec_free(lstmem_chain%next)
   deallocate(lstmem_chain%next)
   nullify(lstmem_chain%next)   
ENDIF
IF(associated(lstmem_chain%LSTMEM_INTWORK))THEN
   call lsquit('error',-1)
ENDIF
IF(associated(lstmem_chain%LSTMEM_INTSWORK))THEN
   call lsquit('error',-1)
ENDIF
IF(associated(lstmem_chain%LSTMEM_REALKWORK))THEN
   call lsquit('error',-1)
ENDIF
deallocate(lstmem_chain)
nullify(lstmem_chain) 
nullify(current)
nullify(last)
lstindex=0
nlstmem = 0
end subroutine lstmem_free

recursive subroutine lstmem_rec_free(chitm)
implicit none
type(lstmemchainitem),pointer :: chitm
integer(kind=long) :: nsize
IF(associated(chitm%next))THEN
   call lstmem_rec_free(chitm%next)
   deallocate(chitm%next)
   nullify(chitm%next)   
ELSE
   nsize=0
   IF(associated(chitm%LSTMEM_INTWORK))THEN
      nsize = nsize + size(chitm%LSTMEM_INTWORK)*mem_intsize
      call mem_dealloc(chitm%LSTMEM_INTWORK)
   ENDIF
   IF(associated(chitm%LSTMEM_INTSWORK))THEN
      nsize = nsize + size(chitm%LSTMEM_INTSWORK)*mem_shortintsize
      call mem_dealloc(chitm%LSTMEM_INTSWORK)
   ENDIF
   IF(associated(chitm%LSTMEM_REALKWORK))THEN
      nsize = nsize + size(chitm%LSTMEM_REALKWORK)*mem_realsize
      call mem_dealloc(chitm%LSTMEM_REALKWORK)
   ENDIF
   if(nsize.GT.0)call mem_deallocated_mem_lstensor(nsize)
ENDIF
end subroutine lstmem_rec_free

subroutine init_lstensorMem(maxInt,maxRealk,maxIntS,index)
implicit none
integer,intent(out) :: index
integer(kind=long),intent(in) :: maxInt,maxRealk,maxIntS
integer :: n,i
index=1
if(nlstmem.EQ.0)then
   lstindex = 1
   !first lstensors
   call init_lstensorMemchainitem(maxInt,maxRealk,maxIntS,lstmem_chain%next)
   nlstmem = 1
   lstmem_chain%next%index=lstindex
   current => lstmem_chain%next
   allocate(lstmem_chain%next%next)
   nullify(lstmem_chain%next%next%LSTMEM_INTWORK)
   nullify(lstmem_chain%next%next%LSTMEM_INTSWORK)
   nullify(lstmem_chain%next%next%LSTMEM_REALKWORK)
   nullify(lstmem_chain%next%next%next)
   last => lstmem_chain%next%next
else
   lstindex = lstindex + 1
   call init_lstensorMemchainitem(maxInt,maxRealk,maxIntS,last)
   last%index = lstindex
   index = lstindex
   current => last
   IF(associated(current%next))call lsquit('error init_lstensorMem',-1)

   allocate(current%next)
   nullify(current%next%LSTMEM_INTWORK)
   nullify(current%next%LSTMEM_INTSWORK)
   nullify(current%next%LSTMEM_REALKWORK)
   nullify(current%next%next)
   last => current%next
   nlstmem = nlstmem+1
endif
end subroutine init_lstensorMem

subroutine zero_lstensorMem()
implicit none
integer :: n
n=SIZE(current%LSTMEM_REALKWORK)
call ls_dzero(current%LSTMEM_REALKWORK,n)
n=SIZE(current%LSTMEM_INTWORK)
CALL LS_IZERO(current%LSTMEM_INTWORK,n)
n=SIZE(current%LSTMEM_INTSWORK)
call ls_sizero(current%LSTMEM_INTSWORK,n)
end subroutine zero_lstensorMem

subroutine zero_lstensorMemindex(index)
implicit none
integer :: n,i,index
if(lstmem_chain%next%index.EQ.index)then
   call zero_lstmemchainitem(lstmem_chain%next)
else
   IF(associated(lstmem_chain%next%next))then
      IF(current%index.EQ.index)then
         call zero_lstmemchainitem(current)
      else
         call zero_rec_lstensorMemindex(index,lstmem_chain%next%next)
      endif
   endif
endif
end subroutine zero_lstensorMemindex

recursive subroutine zero_rec_lstensorMemindex(index,chit)
implicit none
integer,intent(in) :: index
type(lstmemchainitem),pointer :: chit

if(chit%index.EQ.index)then
   call zero_lstmemchainitem(chit)
else
   IF(associated(chit%next))then
      call zero_rec_lstensorMemindex(index,chit%next)
   ELSE
      call lsquit('did not find index zero_rec_lstensorMemindex',-1)
   endif
endif
end subroutine zero_rec_lstensorMemindex

subroutine zero_lstmemchainitem(chit)
type(lstmemchainitem),pointer :: chit
!
integer :: n
n=SIZE(chit%LSTMEM_REALKWORK)
call ls_dzero(chit%LSTMEM_REALKWORK,n)
!n=SIZE(chit%LSTMEM_INTWORK)
!CALL LS_IZERO(chit%LSTMEM_INTWORK,n)
n=SIZE(chit%LSTMEM_INTSWORK)
call ls_sizero(chit%LSTMEM_INTSWORK,n)
end subroutine zero_lstmemchainitem

subroutine set_lstmemrealkbufferpointer(index,buffer,nbuffer)
implicit none
integer,intent(in)  :: index
integer(kind=long)  :: nbuffer
real(realk),pointer :: buffer(:)

if(lstmem_chain%next%index.EQ.index)then
   nbuffer = SIZE(lstmem_chain%LSTMEM_REALKWORK)
   buffer => lstmem_chain%LSTMEM_REALKWORK
else
   IF(associated(lstmem_chain%next%next))then
      call set_rec_lstmemrealkbufferpointer(index,buffer,nbuffer,lstmem_chain%next%next)
   endif
endif
end subroutine set_lstmemrealkbufferpointer

recursive subroutine set_rec_lstmemrealkbufferpointer(index,buffer,nbuffer,chit)
implicit none
integer,intent(in) :: index
integer(kind=long)  :: nbuffer
real(realk),pointer :: buffer(:)
type(lstmemchainitem),pointer :: chit

if(chit%index.EQ.index)then
   nbuffer = SIZE(chit%LSTMEM_REALKWORK)
   buffer => chit%LSTMEM_REALKWORK
else
   IF(associated(chit%next))then
      call set_rec_lstmemrealkbufferpointer(index,buffer,nbuffer,chit%next)
   ELSE
      call lsquit('did not find index set_rec_lstmemrealkbufferpointer',-1)
   endif
endif
end subroutine set_rec_lstmemrealkbufferpointer

subroutine print_lstmemrealkbufferpointer(index)
implicit none
integer,intent(in)  :: index

if(lstmem_chain%next%index.EQ.index)then
   print*,'size(lstmembuffer):',size(lstmem_chain%next%LSTMEM_REALKWORK)
   print*,'lstmembuffer:',lstmem_chain%next%LSTMEM_REALKWORK
else
   IF(associated(lstmem_chain%next%next))then
      call print_rec_lstmemrealkbufferpointer(index,lstmem_chain%next%next)
   endif
endif
end subroutine print_lstmemrealkbufferpointer

recursive subroutine print_rec_lstmemrealkbufferpointer(index,chit)
implicit none
integer,intent(in) :: index
type(lstmemchainitem),pointer :: chit

if(chit%index.EQ.index)then
   print*,'size(lstmembuffer):',size(chit%LSTMEM_REALKWORK)
   print*,'lstmembuffer:',chit%LSTMEM_REALKWORK
else
   IF(associated(chit%next))then
      call print_rec_lstmemrealkbufferpointer(index,chit%next)
   ELSE
      call lsquit('did not find index print_rec_lstmemrealkbufferpointer',-1)
   endif
endif
end subroutine print_rec_lstmemrealkbufferpointer

subroutine retrieve_lstMemVal(index,maxInt,maxRealk,maxIntS)
implicit none
integer,intent(in) :: index
integer(kind=long),intent(out) :: maxInt,maxRealk,maxIntS
if(lstmem_chain%next%index.EQ.index)then
   maxInt = SIZE(lstmem_chain%next%LSTMEM_INTWORK)
   maxRealk = SIZE(lstmem_chain%next%LSTMEM_REALKWORK)
   maxIntS = SIZE(lstmem_chain%next%LSTMEM_INTSWORK)
   current => lstmem_chain%next
else
   IF(associated(lstmem_chain%next%next))then
      call retrieve_rec_lstMemVal(index,maxInt,maxRealk,maxIntS,lstmem_chain%next%next)
   endif
endif
end subroutine retrieve_lstMemVal

recursive subroutine retrieve_rec_lstMemVal(index,maxInt,maxRealk,maxIntS,chit)
implicit none
integer,intent(in) :: index
integer(kind=long),intent(out) :: maxInt,maxRealk,maxIntS
type(lstmemchainitem),pointer :: chit

if(chit%index.EQ.index)then
   maxInt = SIZE(chit%LSTMEM_INTWORK)
   maxRealk = SIZE(chit%LSTMEM_REALKWORK)
   maxIntS = SIZE(chit%LSTMEM_INTSWORK)
   current => chit
else
   IF(associated(chit%next))then
      call retrieve_rec_lstMemVal(index,maxInt,maxRealk,maxIntS,chit%next)
   ELSE
      call lsquit('did not find index retrieve_rec_lstMemVal',-1)
   endif
endif
end subroutine retrieve_rec_lstMemVal

subroutine copy_lstensorMemToCurrent(index)
implicit none
integer,intent(in) :: index
integer(kind=long) :: maxInt,maxRealk,maxIntS
if(lstmem_chain%next%index.EQ.index)then
   maxInt = SIZE(current%LSTMEM_INTWORK)
   maxRealk = SIZE(current%LSTMEM_REALKWORK)
   maxIntS = SIZE(current%LSTMEM_INTSWORK)
#if VAR_LSDEBUGINT
   IF(maxInt.NE.SIZE(lstmem_chain%next%LSTMEM_INTWORK))&
        & call lsquit('dim1 mismatch lstememcopy',-1)
   IF(maxIntS.NE.SIZE(lstmem_chain%next%LSTMEM_INTSWORK))&
        & call lsquit('dim2 mismatch lstememcopy',-1)
   IF(maxRealk.NE.SIZE(lstmem_chain%next%LSTMEM_REALKWORK))&
        & call lsquit('dim3 mismatch lstememcopy',-1)
#endif
   call LS_DCOPY(maxRealk,lstmem_chain%next%LSTMEM_REALKWORK,current%LSTMEM_REALKWORK)
   call LS_ICOPY(maxInt,lstmem_chain%next%LSTMEM_INTWORK,current%LSTMEM_INTWORK)
   call LS_SICOPY(maxIntS,lstmem_chain%next%LSTMEM_INTSWORK,current%LSTMEM_INTSWORK)
else
   IF(associated(lstmem_chain%next%next))then
      call copy_rec_lstensorMemToCurrent(index,lstmem_chain%next%next)
   endif
endif
end subroutine copy_lstensorMemToCurrent

recursive subroutine copy_rec_lstensorMemToCurrent(index,chit)
implicit none
integer,intent(in) :: index
integer(kind=long) :: maxInt,maxRealk,maxIntS
type(lstmemchainitem),pointer :: chit
if(chit%index.EQ.index)then
   maxInt = SIZE(chit%LSTMEM_INTWORK)
   maxRealk = SIZE(chit%LSTMEM_REALKWORK)
   maxIntS = SIZE(chit%LSTMEM_INTSWORK)
#if VAR_LSDEBUGINT
   IF(maxInt.NE.SIZE(current%LSTMEM_INTWORK))&
        & call lsquit('dim1 mismatch lstememcopy',-1)
   IF(maxIntS.NE.SIZE(current%LSTMEM_INTSWORK))&
        & call lsquit('dim2 mismatch lstememcopy',-1)
   IF(maxRealk.NE.SIZE(current%LSTMEM_REALKWORK))&
        & call lsquit('dim3 mismatch lstememcopy',-1)
#endif
   call LS_DCOPY(maxRealk,chit%LSTMEM_REALKWORK,current%LSTMEM_REALKWORK)
   call LS_ICOPY(maxInt,chit%LSTMEM_INTWORK,current%LSTMEM_INTWORK)
   call LS_SICOPY(maxIntS,chit%LSTMEM_INTSWORK,current%LSTMEM_INTSWORK)
else
   IF(associated(chit%next))then
      call copy_rec_lstensorMemToCurrent(index,chit%next)
   ELSE
      call lsquit('did not find index copy_rec_lstensorMemToCurrent',-1)
   endif
endif
end subroutine copy_rec_lstensorMemToCurrent

subroutine init_lstensorMemchainitem(maxInt,maxRealk,maxIntS,chainitm)
implicit none
integer(kind=long),intent(in) :: maxInt,maxRealk,maxIntS
type(lstmemchainitem) :: chainitm
integer(kind=long) :: maxInt2,maxRealk2,maxIntS2,nsize
IF(maxInt.EQ.0)then
   maxint2 = 1
ELSE
   maxint2 = maxInt
ENDIF
IF(maxIntS.EQ.0)then
   maxintS2 = 1
ELSE
   maxintS2 = maxIntS
ENDIF
IF(maxRealk.EQ.0)then
   maxRealk2 = 1
ELSE
   maxRealk2 = maxRealk
ENDIF
chainitm%lstmem_INTWORKLENGTH = maxInt2
chainitm%lstmem_INTSWORKLENGTH = maxIntS2
chainitm%lstmem_REALKWORKLENGTH = maxRealk2
call mem_alloc(chainitm%LSTMEM_INTWORK,maxInt2)
call mem_alloc(chainitm%LSTMEM_INTSWORK,maxIntS2)
call mem_alloc(chainitm%LSTMEM_REALKWORK,maxRealk2)
nsize=maxInt2*mem_intsize+maxIntS2*mem_shortintsize+maxRealk2*mem_realsize
call mem_allocated_mem_lstensor(nsize)

chainitm%LSTMEM_LINTWORK = 1
chainitm%LSTMEM_LINTSWORK = 1
chainitm%LSTMEM_LREALKWORK = 1
end subroutine init_lstensorMemchainitem

subroutine free_lstensorMem(index)
implicit none
integer,intent(in)  :: index
logical :: found
type(lstmemchainitem),pointer :: next

if(lstmem_chain%next%index.EQ.index)then
   IF (ASSOCIATED(lstmem_chain%next%next)) THEN
     next => lstmem_chain%next%next
     call free_lstensorMemchainitem(lstmem_chain%next)
     deallocate(lstmem_chain%next)
     lstmem_chain%next => next
    ELSE
     call free_lstensorMemchainitem(lstmem_chain%next)
     deallocate(lstmem_chain%next)
     nullify(lstmem_chain%next)
    ENDIF
else
   IF(associated(lstmem_chain%next%next))then
      found=.FALSE.
      call free_rec_lstensorMem(index,lstmem_chain%next,found)
      IF(found)then
         ! Remove freed lstmemchainitem from chain
         next => lstmem_chain%next%next
         deallocate(lstmem_chain%next)
         lstmem_chain%next=> next
      endif
   else
      call lsquit('index not found A1',-1)
   endif
endif
nlstmem = nlstmem-1
IF(nlstmem.EQ.0)then
   IF(associated(lstmem_chain%next%next))THEN
      call lstmem_rec_free(lstmem_chain%next%next)
      deallocate(lstmem_chain%next%next)
      nullify(lstmem_chain%next%next)   
   ENDIF
   nullify(lstmem_chain%next%next)
   nullify(current)
   nullify(last)
   lstindex=0
   nlstmem = 0
endif
end subroutine free_lstensorMem

recursive subroutine free_rec_lstensorMem(index,chit,found)
implicit none
integer,intent(in)  :: index
type(lstmemchainitem),pointer :: chit,next
logical,intent(inout) :: found
if(chit%index.EQ.index)then
   call free_lstensorMemchainitem(chit)
   found=.TRUE.
else
   IF(associated(chit%next))then
      found=.false.
      call free_rec_lstensorMem(index,chit%next,found)
      IF(found)then
         ! Remove freed lstmemchainitem from chain
         next => chit%next%next
         deallocate(chit%next)
         chit%next => next
         found=.FALSE.
      endif
   else
      call lsquit('index not found A2',-1)
   endif
endif
end subroutine free_rec_lstensorMem

subroutine free_lstensorMemchainitem(chainitm)
implicit none
type(lstmemchainitem),pointer :: chainitm
integer(kind=long) :: nsize
chainitm%lstmem_INTWORKLENGTH = 0
chainitm%lstmem_INTSWORKLENGTH = 0
chainitm%lstmem_REALKWORKLENGTH = 0
nsize=size(chainitm%LSTMEM_INTWORK)*mem_intsize&
     & + size(chainitm%LSTMEM_INTSWORK)*mem_shortintsize&
     & + size(chainitm%LSTMEM_REALKWORK)*mem_realsize
call mem_deallocated_mem_lstensor(nsize)

call mem_dealloc(chainitm%LSTMEM_INTWORK)
call mem_dealloc(chainitm%LSTMEM_INTSWORK)
call mem_dealloc(chainitm%LSTMEM_REALKWORK)

nullify(chainitm%LSTMEM_INTWORK)
nullify(chainitm%LSTMEM_INTSWORK)
nullify(chainitm%LSTMEM_REALKWORK)
chainitm%index = 0
chainitm%LSTMEM_LINTWORK = 1
chainitm%LSTMEM_LINTSWORK = 1
chainitm%LSTMEM_LREALKWORK = 1
end subroutine free_lstensorMemchainitem

subroutine real_setLSTpointer_alloc(A,n)
implicit none
integer,intent(in)  :: n
REAL(REALK),pointer :: A(:)
integer(kind=long) :: LREALKWORK,WORKLENGTH
nullify(A)
LREALKWORK = current%LSTMEM_LREALKWORK
WORKLENGTH = current%lstmem_REALKWORKLENGTH
IF(LREALKWORK+n-1.GT.WORKLENGTH)then
   print*,'n               ',n
   print*,'LREALKWORK      ',LREALKWORK
   print*,'REALKWORKLENGTH ',WORKLENGTH
   print*,'size            ',size(current%LSTMEM_REALKWORK)
   call lsquit('Real Work buffer full in lstensor',-1)
ENDIF
A => current%LSTMEM_REALKWORK(LREALKWORK:LREALKWORK+n-1)
LREALKWORK = LREALKWORK+n
current%LSTMEM_LREALKWORK = LREALKWORK
end subroutine real_setLSTpointer_alloc

subroutine int_setLSTpointer_alloc(A,n)
implicit none
integer,intent(in)  :: n
integer,pointer :: A(:)
integer(kind=long) :: LINTWORK,WORKLENGTH
nullify(A)
LINTWORK = current%LSTMEM_LINTWORK
WORKLENGTH = current%lstmem_INTWORKLENGTH
IF(LINTWORK+n-1.GT.WORKLENGTH)then
   print*,'n               ',n
   print*,'LINTWORK      ',LINTWORK
   print*,'INTWORKLENGTH ',WORKLENGTH
   print*,'size            ',size(current%LSTMEM_INTWORK)
   call lsquit('Int Work buffer full in lstensor',-1)
ENDIF
A => current%LSTMEM_INTWORK(LINTWORK:LINTWORK+n-1)
LINTWORK = LINTWORK+n
current%LSTMEM_LINTWORK = LINTWORK
end subroutine int_setLSTpointer_alloc

subroutine ints_setLSTpointer_alloc(A,n)
implicit none
integer,intent(in)  :: n
integer(kind=short),pointer :: A(:)
integer(kind=long) :: LINTSWORK,WORKLENGTH
nullify(A)
LINTSWORK = current%LSTMEM_LINTSWORK
WORKLENGTH = current%lstmem_INTSWORKLENGTH
IF(LINTSWORK+n-1.GT.WORKLENGTH)then
   print*,'n               ',n
   print*,'LINTWORK      ',LINTSWORK
   print*,'INTSWORKLENGTH ',WORKLENGTH
   print*,'size            ',size(current%LSTMEM_INTSWORK)
   call lsquit('int(short) Work buffer full in lstensor',-1)
ENDIF
A => current%LSTMEM_INTSWORK(LINTSWORK:LINTSWORK+n-1)
LINTSWORK = LINTSWORK+n
current%LSTMEM_LINTSWORK = LINTSWORK
end subroutine ints_setLSTpointer_alloc

END MODULE LstensorMem
