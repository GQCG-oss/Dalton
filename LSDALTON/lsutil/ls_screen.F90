!> @file 
!> Contain the IOITEM which keeps track of files written to disk
MODULE screen_mod
use precision
use LSTENSOR_OPERATIONSMOD
use io
use io_type
use lsmpi_type
#ifdef VAR_MPI
use lsmpi_param
#endif
use lsparameters

TYPE SCREENCHAINITEM
Character(80) :: filename
TYPE(LSTENSORp) :: LST
TYPE(SCREENCHAINITEM),pointer :: next
END TYPE SCREENCHAINITEM

TYPE SCREENITEM
TYPE(SCREENCHAINITEM),pointer :: start
TYPE(SCREENCHAINITEM),pointer :: end
END TYPE SCREENITEM

TYPE DECSCREENITEM
TYPE(LSTENSOR),pointer :: MasterGabRHS
TYPE(LSTENSOR),pointer :: MasterGabLHS
TYPE(LSTENSORp),pointer :: batchGab(:,:)
TYPE(LSTENSORp),pointer :: batchGabKLHS(:)
TYPE(LSTENSORp),pointer :: batchGabKRHS(:)
END TYPE DECSCREENITEM

TYPE(SCREENITEM),save :: screenFromlsScreen
CONTAINS
!> \brief initialise the SCREENitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param SCREENFROMLSSCREEN the SCREENITEM
SUBROUTINE screen_init()!(SCREEN)
  implicit none
!  TYPE(SCREENITEM)  :: SCREEN
  nullify(SCREENFROMLSSCREEN%start)
  nullify(SCREENFROMLSSCREEN%end)
#ifdef VAR_MPI
  IF(infpar%mynum.EQ.infpar%master)THEN
     call ls_mpibcast(IISCREENINIT,infpar%master,MPI_COMM_LSDALTON)
  ENDIF
#endif
END SUBROUTINE screen_init

SUBROUTINE nullify_decscreen(SCREEN)
  implicit none
  TYPE(DECSCREENITEM)  :: SCREEN
  nullify(SCREEN%masterGabRHS)
  nullify(SCREEN%masterGabLHS)
  nullify(SCREEN%batchGab)
  nullify(SCREEN%batchGabKLHS)
  nullify(SCREEN%batchGabKRHS)
END SUBROUTINE nullify_decscreen

SUBROUTINE null_decscreen_and_associate_MasterGab_LHS(GAB,SCREEN)
  implicit none
  type(LSTENSOR),pointer :: GAB
  TYPE(DECSCREENITEM)  :: SCREEN
  nullify(SCREEN%masterGabLHS)
  allocate(SCREEN%masterGabLHS)
  call lstensor_nullify(SCREEN%masterGabLHS)
  GAB => SCREEN%masterGabLHS
End SUBROUTINE null_decscreen_and_associate_MasterGab_LHS

SUBROUTINE null_decscreen_and_associate_MasterGab_RHS(GAB,SCREEN)
  implicit none
  type(LSTENSOR),pointer :: GAB
  TYPE(DECSCREENITEM)  :: SCREEN
  nullify(SCREEN%masterGabRHS)
  allocate(SCREEN%masterGabRHS)
  call lstensor_nullify(SCREEN%masterGabRHS)
  GAB => SCREEN%masterGabRHS
  nullify(SCREEN%batchGab)
  nullify(SCREEN%batchGabKLHS)
  nullify(SCREEN%batchGabKRHS)
End SUBROUTINE null_decscreen_and_associate_MasterGab_RHS

SUBROUTINE init_decscreen_batch(ndimA,ndimG,SCREEN)
  implicit none
  integer :: ndimA,ndimG
  TYPE(DECSCREENITEM)  :: SCREEN
!
  integer :: ndim,I,J
  nullify(SCREEN%batchGab)
  allocate(SCREEN%batchGab(ndimA,ndimG))
  do J=1,ndimG
     do I=1,ndimA
        nullify(SCREEN%batchGab(I,J)%p)
     enddo
  enddo
  nullify(SCREEN%batchGabKLHS)
  allocate(SCREEN%batchGabKLHS(ndimA))
  do I=1,ndimA
     nullify(SCREEN%batchGabKLHS(I)%p)
  enddo
  nullify(SCREEN%batchGabKRHS)
  allocate(SCREEN%batchGabKRHS(ndimG))
  do J=1,ndimG
     nullify(SCREEN%batchGabKRHS(J)%p)
  enddo
END SUBROUTINE init_decscreen_batch

SUBROUTINE free_decscreen(SCREEN)
  implicit none
  TYPE(DECSCREENITEM)  :: SCREEN
!
  integer :: I,J,n1,n2
  if(associated(SCREEN%masterGabLHS))THEN
     call lstensor_free(SCREEN%masterGabLHS)
     deallocate(SCREEN%masterGabLHS)
     nullify(SCREEN%masterGabLHS)
  ENDIF
  if(associated(SCREEN%masterGabRHS))THEN
     call lstensor_free(SCREEN%masterGabRHS)
     deallocate(SCREEN%masterGabRHS)
     nullify(SCREEN%masterGabRHS)
  ENDIF
  if(associated(SCREEN%batchGab))THEN
     n2 = size(SCREEN%batchGab,2)
     n1 = size(SCREEN%batchGab,1)
!     IF(n1.EQ.1.AND.n2.EQ.1)THEN
!        IF(SCREEN%batchGab(1,1)%screenTensor)THEN
!           call lstensor_free(SCREEN%batchGab(1,1))
!        ENDIF
!     ELSE
        do J=1,n2
           do I=1,n1
              if(associated(SCREEN%batchGab(I,J)%p))THEN
                 call lstensor_free(SCREEN%batchGab(I,J)%p)
                 deallocate(SCREEN%batchGab(I,J)%p)
                 nullify(SCREEN%batchGab(I,J)%p)
              endif
           enddo
        enddo
!     ENDIF
     deallocate(SCREEN%batchGab)
     nullify(SCREEN%batchGab)
  ENDIF  
  if(associated(SCREEN%batchGabKLHS))THEN
     n1 = size(SCREEN%batchGabKLHS)
     do I=1,n1
        if(associated(SCREEN%batchGabKLHS(I)%p))THEN
           call lstensor_free(SCREEN%batchGabKLHS(I)%p)
           deallocate(SCREEN%batchGabKLHS(I)%p)
           nullify(SCREEN%batchGabKLHS(I)%p)
        endif
     enddo
     deallocate(SCREEN%batchGabKLHS)
     nullify(SCREEN%batchGabKLHS)
  ENDIF  
  if(associated(SCREEN%batchGabKRHS))THEN
     n1 = size(SCREEN%batchGabKRHS)
     do I=1,n1
        if(associated(SCREEN%batchGabKRHS(I)%p))THEN
           call lstensor_free(SCREEN%batchGabKRHS(I)%p)
           deallocate(SCREEN%batchGabKRHS(I)%p)
           nullify(SCREEN%batchGabKRHS(I)%p)
        endif
     enddo
     deallocate(SCREEN%batchGabKRHS)
     nullify(SCREEN%batchGabKRHS)
  ENDIF  
END SUBROUTINE free_decscreen

SUBROUTINE decscreen_associateMaster_RHS(GAB,SCREEN)
implicit none 
type(LSTENSOR),pointer :: GAB
TYPE(DECSCREENITEM),intent(in)  :: SCREEN
GAB => SCREEN%masterGabRHS
END SUBROUTINE DECSCREEN_ASSOCIATEMASTER_RHS

SUBROUTINE decscreen_associateMaster_LHS(GAB,SCREEN)
implicit none 
type(LSTENSOR),pointer :: GAB
TYPE(DECSCREENITEM),intent(in)  :: SCREEN
GAB => SCREEN%masterGabLHS
END SUBROUTINE DECSCREEN_ASSOCIATEMASTER_LHS

!SUBROUTINE copy_and_alloc_screenitem(SCREENNEW,SCREENOLD)
!  implicit none
!  TYPE(SCREENITEM)  :: SCREENNEW,SCREENOLD
!  call screen_init(SCREENNEW)
!  IF(associated(SCREENOLD%start))THEN
!     allocate(SCREENNEW%start)
!     nullify(SCREENNEW%start%next)
!     call copy_screenchain(SCREENNEW,SCREENNEW%start,SCREENOLD%start)
!  ENDIF
!END SUBROUTINE copy_and_alloc_screenitem

!!$recursive subroutine copy_screenchain(SCREENNEW,ScreenChainNew,ScreenChainOld)
!!$  implicit none  
!!$  TYPE(SCREENITEM)  :: SCREENNEW
!!$  type(screenchainitem),target ::ScreenChainNEW
!!$  type(screenchainitem) :: ScreenChainOLD
!!$  type(LSTENSOR),pointer :: GAB  
!!$
!!$  ScreenChainNew%filename = ScreenChainOld%filename
!!$  nullify(ScreenChainNew%LST%p)
!!$  allocate(ScreenChainNew%LST%p)   
!!$  call copy_lstensor_to_lstensor(ScreenChainOld%LST%p,ScreenChainNew%LST%p)
!!$  IF(associated(ScreenchainOld%next))THEN
!!$     allocate(ScreenChainNew%next)
!!$     nullify(ScreenChainNew%next%next)
!!$     call copy_screenchain(SCREENNEW,ScreenChainNew%next,ScreenChainOld%next)
!!$  ELSE
!!$     SCREENNEW%end => ScreenchainNew
!!$  ENDIF
!!$
!!$end subroutine copy_screenchain

SUBROUTINE screen_add_associate_item(GAB,FILENAME)
implicit none
type(LSTENSOR),pointer :: GAB
!TYPE(SCREENITEM)  :: SCREEN
Character(80),intent(in) :: Filename

IF(.NOT.associated(SCREENFROMLSSCREEN%start))THEN
   allocate(SCREENFROMLSSCREEN%start)
   SCREENFROMLSSCREEN%start%filename = Filename
   nullify(SCREENFROMLSSCREEN%start%LST%p)
   allocate(SCREENFROMLSSCREEN%start%LST%p)
   GAB => SCREENFROMLSSCREEN%start%LST%p
   SCREENFROMLSSCREEN%end => SCREENFROMLSSCREEN%start
   nullify(SCREENFROMLSSCREEN%start%next)
ELSE
   allocate(SCREENFROMLSSCREEN%end%next)
   nullify(SCREENFROMLSSCREEN%end%next%next)
   SCREENFROMLSSCREEN%end%next%filename = Filename
   nullify(SCREENFROMLSSCREEN%end%next%LST%p)
   allocate(SCREENFROMLSSCREEN%end%next%LST%p)
   GAB => SCREENFROMLSSCREEN%end%next%LST%p
   SCREENFROMLSSCREEN%end => SCREENFROMLSSCREEN%end%next
ENDIF
END SUBROUTINE screen_add_associate_item

SUBROUTINE determine_lst_in_screenlist(Filename,FoundInMem,IO)!SCREEN)
implicit none 
Character(80),intent(in)       :: Filename
logical,intent(out) :: FoundInMem
TYPE(IOITEM),intent(in)  :: IO
!TYPE(SCREENITEM),intent(in)  :: SCREEN
!
integer :: I
FoundInMem=.FALSE.
IF(associated(SCREENFROMLSSCREEN%start))THEN
   call FindFilename(SCREENFROMLSSCREEN%start,Filename,FoundInMem)
ENDIF
END SUBROUTINE DETERMINE_LST_IN_SCREENLIST

recursive subroutine FindFilename(ScreenChain,Filename,FoundInMem)
implicit none
type(screenchainitem) ::ScreenChain
Character(80),intent(in) :: Filename
logical,intent(inout) :: FoundInMem

IF(Screenchain%Filename.EQ.Filename)THEN
   FoundInMem=.TRUE.
ELSE
   IF(associated(SCREENchain%next))THEN
      call FindFilename(SCREENchain%next,Filename,FoundInMem)
   ENDIF
ENDIF
end subroutine FindFilename

SUBROUTINE screen_associate(GAB,Filename,FoundInMem)!,SCREEN)
implicit none 
type(LSTENSOR),pointer :: GAB
Character(80),intent(in)       :: Filename
logical,intent(out) :: FoundInMem
!TYPE(SCREENITEM),intent(in)  :: SCREEN
!
integer :: I
call associateGAB(GAB,SCREENFROMLSSCREEN%start,Filename,FoundInMem)
IF(.NOT.FoundInMem)THEN
   CALL LSQUIT('Error in SCREEN_ASSOCIATE did not find Gab in screenitem',-1)
ENDIF
END SUBROUTINE SCREEN_ASSOCIATE

recursive subroutine associateGAB(GAB,ScreenChain,Filename,FoundInMem)
implicit none
type(LSTENSOR),pointer :: GAB
type(screenchainitem) ::ScreenChain
Character(80),intent(in)       :: Filename
logical,intent(inout) :: FoundInMem

IF(Screenchain%Filename.EQ.Filename)THEN
   FoundInMem=.TRUE.
   GAB => Screenchain%LST%p
   RETURN
ELSE
   IF(associated(SCREENchain%next))THEN
      call associateGAB(GAB,ScreenChain%next,Filename,FoundInMem)
   ENDIF
ENDIF
end subroutine AssociateGAB

!> \brief free the SCREENitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param SCREEN the SCREENITEM
SUBROUTINE screen_free()!(SCREEN)
implicit none
!TYPE(SCREENITEM)  :: SCREEN
IF(associated(screenFromlsScreen%start))THEN
   call screenchain_free(ScreenFromlsScreen%start)
   deallocate(ScreenFromlsScreen%start)
ENDIF
nullify(ScreenFromlsScreen%start)
nullify(ScreenFromlsScreen%end)
#ifdef VAR_MPI
  IF(infpar%mynum.EQ.infpar%master)THEN
     call ls_mpibcast(IISCREENFREE,infpar%master,MPI_COMM_LSDALTON)
  ENDIF
#endif
END SUBROUTINE screen_free

!> \brief free the SCREENitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param SCREEN the SCREENITEM
recursive SUBROUTINE screenchain_free(SCREENCHAIN)
implicit none
TYPE(SCREENCHAINITEM)  :: SCREENCHAIN
if(associated(screenchain%next))then
   call screenchain_free(Screenchain%next)
   deallocate(Screenchain%next)
endif
call lstensor_free(SCREENchain%LST%p)
deallocate(Screenchain%LST%p)
nullify(Screenchain%LST%p)
nullify(Screenchain%next)

END SUBROUTINE screenchain_free

SUBROUTINE PRINT_SCREENITEM(LUPRI)
implicit none
!TYPE(SCREENITEM) :: SCREEN
INTEGER :: LUPRI

WRITE(lupri,*)'Screening matrices in SCREENITEM'
IF(associated(SCREENFROMLSSCREEN%start))THEN
   CALL PRINT_SCREENCHAINITEM(SCREENFROMLSSCREEN%start,LUPRI)
ENDIF

END SUBROUTINE PRINT_SCREENITEM

recursive SUBROUTINE PRINT_SCREENCHAINITEM(SCREENCHAIN,LUPRI)
implicit none
TYPE(SCREENCHAINITEM) :: SCREENCHAIN
INTEGER :: LUPRI
WRITE(lupri,'(A,A80)')'  ',SCREENCHAIN%filename
call lstensor_Print(SCREENCHAIN%LST%p,LUPRI)
IF(associated(SCREENCHAIN%next))THEN
   CALL PRINT_SCREENCHAINITEM(SCREENCHAIN%next,LUPRI)
ENDIF
END SUBROUTINE PRINT_SCREENCHAINITEM

end MODULE screen_mod
