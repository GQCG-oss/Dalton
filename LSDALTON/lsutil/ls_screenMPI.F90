!> @file 
!> Contain the IOITEM which keeps track of files written to disk
MODULE screen_modMPI
use screen_mod
use lsmpi_type
use lsmpi_op

#ifdef VAR_MPI
CONTAINS
SUBROUTINE mpicopy_screen(Slave,Master)
implicit none
logical :: slave
integer(kind=ls_mpik) :: master
!
logical :: assStart

IF(SLAVE)call screen_init()
assStart = associated(SCREENFROMLSSCREEN%start)
call LS_MPI_BUFFER(assStart,Master)
IF(assStart)THEN
   IF(SLAVE)THEN
      allocate(SCREENFROMLSSCREEN%start)
      nullify(SCREENFROMLSSCREEN%start%next)
   ENDIF
   call mpicopy_screenchain(SCREENFROMLSSCREEN,SCREENFROMLSSCREEN%start,slave,Master)
ELSE
   IF(SLAVE)nullify(SCREENFROMLSSCREEN%start)
ENDIF
end SUBROUTINE mpicopy_screen

recursive SUBROUTINE mpicopy_screenchain(screen,screenchain,Slave,Master)
implicit none
type(screenitem) :: screen
type(screenchainitem),pointer :: ScreenChain
logical :: slave
integer(kind=ls_mpik) :: master
!
type(LSTENSOR),pointer :: GAB  
logical :: assStart

call LS_MPI_BUFFER(ScreenChain%filename,80,master)
IF(SLAVE)THEN
   nullify(ScreenChain%LST%p)
   allocate(ScreenChain%LST%p)   
ENDIF
call mpicopy_lstensor(ScreenChain%LST%p,slave,master)
assStart = associated(Screenchain%next)
call LS_MPI_BUFFER(assStart,Master)
IF(assStart)THEN
   IF(SLAVE)THEN
      allocate(ScreenChain%next)
      nullify(ScreenChain%next%next)
   ENDIF
   call mpicopy_screenchain(SCREEN,ScreenChain%next,Slave,Master)
ELSE
   IF(SLAVE)THEN
      nullify(Screenchain%next)
      SCREEN%end => Screenchain
   ENDIF
ENDIF
end SUBROUTINE mpicopy_screenchain
#else

contains

!Added to avoid "has no symbols" linking warning
subroutine screen_modMPI_void()
end subroutine screen_modMPI_void
#endif

end MODULE screen_modMPI
