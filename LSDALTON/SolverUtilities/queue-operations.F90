!> @file
!> Contains queue type definition module.

!> \brief Contains operations for type(modFIFO).
!> \author Stinne Host
!> \date 2007
MODULE queue_ops
use files
use queue_module
use Matrix_Operations

contains
!> \brief Initialize a type(modFIFO).
!> \author Stinne Host
!> \date 2007
!> \param queue The type(modFIFO) to be initialized
!> \param qsize The maximum of matrices to store in each of the two queues, including expansion points
!> \param rowdim Row size of matrices to be stored (rowdim x coldim)
!> \param coldim Column size of matrices to be stored (rowdim x coldim)
!> \param disk True if vectors in queue should be kept on disk instead of in core
   SUBROUTINE modFIFO_init(queue,qsize,rowdim,coldim,disk)
      implicit none
      TYPE(modFIFO) :: queue
      integer, intent(in) :: qsize 
      integer, intent(in) :: rowdim 
      integer, intent(in) :: coldim 
      logical, intent(in) :: disk
      integer :: i, j
         queue%disk = disk
         queue%queuesize = qsize
         queue%offset = 0  
         if (queue%disk) then
            allocate(queue%iDarray(qsize-1))
            allocate(queue%iFarray(qsize-1))
            do i = 1, qsize-1
               nullify(queue%iDarray(i)%p)
               nullify(queue%iFarray(i)%p)
            enddo
            nullify(queue%iD_exp)
            nullify(queue%iF_exp)

            call mat_init(queue%tempD,rowdim,coldim)
            call mat_init(queue%tempF,rowdim,coldim)
         else
            allocate(queue%Darray(qsize-1))
            allocate(queue%Farray(qsize-1))
            do i = 1, qsize-1
               nullify(queue%Darray(i)%p)
               nullify(queue%Farray(i)%p)
            enddo
            nullify(queue%D_exp)
            nullify(queue%F_exp)
         endif
         queue%qcounter = 1
   END SUBROUTINE modFIFO_init

!> \brief Free a type(modFIFO).
!> \author Stinne Host
!> \date 2007
!> \param queue The type(modFIFO) to be free'd
   SUBROUTINE modFIFO_free(queue)
      implicit none
      TYPE(modFIFO) :: queue
      integer :: i
         if (queue%disk) then
            if(associated(queue%iD_exp)) call lsclose(queue%iD_exp,'DELETE')
            if(associated(queue%iF_exp)) call lsclose(queue%iF_exp,'DELETE')

            do i = 1, queue%offset
               call lsclose(queue%iDarray(i)%p,'DELETE')
               call lsclose(queue%iFarray(i)%p,'DELETE')
            enddo

            deallocate(queue%iDarray)
            deallocate(queue%iFarray)

            call mat_free(queue%tempD)
            call mat_free(queue%tempF)
         else
            if(associated(queue%D_exp)) call mat_free(queue%D_exp)
            if(associated(queue%F_exp)) call mat_free(queue%F_exp)

            do i = 1, queue%offset
               call mat_free(queue%Darray(i)%p)
               call mat_free(queue%Farray(i)%p)
            enddo

            deallocate(queue%Darray)
            deallocate(queue%Farray)
         endif
   END SUBROUTINE modFIFO_free

!> \brief Branches out to the proper add-to-queue routine (disk or memory).
!> \author Stinne Host
!> \date May 2010
!> \param queue The type(modFIFO) to add to
!> \param F Matrix to be added to the Farray queue (usually Fock matrix, but could be anything)
!> \param D Matrix to be added to the Darray queue (usually density matrix, but could be anything)
!> \param energy Only referenced if exppoint=true. Energy (or other info) for expansion point to be saved.
!> \param exppoint True if the added pair of matrices should be used as expansion point
   subroutine add_to_modFIFO(queue, F, D, energy, exppoint)
   implicit none
      type(modFIFO), target :: queue
      type(matrix), intent(in) :: F, D
      logical, intent(in) :: exppoint
      real(realk), intent(in) :: energy

      if (queue%disk) then
         call add_to_modFIFO_disk(queue, F, D, energy, exppoint)
      else
         call add_to_modFIFO_memory(queue, F, D, energy, exppoint)
      endif

   end subroutine add_to_modFIFO

!> \brief Add a matrix to a type(modFIFO) kept in memory.
!> \author Stinne Host
!> \date 2007
!>
!> If the queue is not full, all previous data is kept. If the queue is full,
!> the oldest pair of matrices (F,D) is discarded to make room for the new
!> matrices. If the new matrix is not to be used as an expansion point, the old
!> expansion point is left untouched. If the new matrix is to be used as an expansion
!> point, the old expansion point is moved to the "regular" queue as the newest
!> member, discarding a pair (F,D) of old matrices if necessary.
!> 
!> \param queue The type(modFIFO) to add to
!> \param F Matrix to be added to the Farray queue (usually Fock matrix, but could be anything)
!> \param D Matrix to be added to the Darray queue (usually density matrix, but could be anything)
!> \param energy Only referenced if exppoint=true. Energy (or other info) for expansion point to be saved.
!> \param exppoint True if the added pair of matrices should be used as expansion point
   subroutine add_to_modFIFO_memory(queue, F, D, energy, exppoint)
   implicit none
      type(modFIFO), target :: queue
      type(matrix), intent(in) :: F, D
      type(matrix), pointer :: Fpointer, Dpointer  
      type(matrix), pointer :: Dplaceinqueue, Fplaceinqueue
      logical, intent(in) :: exppoint
      logical :: addtoqueue
      real(realk), intent(in) :: energy
      integer :: i, rowdim, coldim

      rowdim = D%nrow
      coldim = D%ncol
      if (queue%qcounter == 200) queue%qcounter = 1

      call mat_init(queue%Fmatrices(queue%qcounter),rowdim,coldim)
      call mat_init(queue%Dmatrices(queue%qcounter),rowdim,coldim)
      call mat_assign(queue%Fmatrices(queue%qcounter),F)
      call mat_assign(queue%Dmatrices(queue%qcounter),D)

      Fpointer=>queue%Fmatrices(queue%qcounter)
      Dpointer=>queue%Dmatrices(queue%qcounter)

      queue%qcounter = queue%qcounter + 1

      addtoqueue = .false.
      if (exppoint) then !if new expansion point, move existing expansion point to queue
         if(associated(queue%D_exp)) then !if it is not the first addition, 
                                          !move current expansion point to queue
            Dplaceinqueue => queue%D_exp
            Fplaceinqueue => queue%F_exp
            addtoqueue = .true.
         endif
         queue%D_exp => Dpointer
         queue%F_exp => Fpointer
         queue%energy = energy
      else
         addtoqueue = .true.
         Dplaceinqueue => Dpointer
         Fplaceinqueue => Fpointer
      endif

      if (addtoqueue) then
         if (queue%offset < queue%queuesize-1) then
            queue%offset = queue%offset + 1
            queue%Darray(queue%offset)%p => Dplaceinqueue
            queue%Farray(queue%offset)%p => Fplaceinqueue
         else
            call mat_free(queue%Darray(1)%p)
            call mat_free(queue%Farray(1)%p)
            do i = 1, queue%queuesize-2
               queue%Darray(i)%p => queue%Darray(i+1)%p
               queue%Farray(i)%p => queue%Farray(i+1)%p
            enddo
            queue%Darray(queue%offset)%p => Dplaceinqueue
            queue%Farray(queue%offset)%p => Fplaceinqueue
         endif
      endif
   end subroutine add_to_modFIFO_memory

!> \brief Add a matrix to a type(modFIFO) kept on disk.
!> \author Stinne Host
!> \date May 2010
!>
!> See add_to_modFIFO_memory above. Instead of in memory, matrices are kept
!> in separate files.
!> 
!> \param queue The type(modFIFO) to add to
!> \param F Matrix to be added to the Farray queue (usually Fock matrix, but could be anything)
!> \param D Matrix to be added to the Darray queue (usually density matrix, but could be anything)
!> \param energy Only referenced if exppoint=true. Energy (or other info) for expansion point to be saved.
!> \param exppoint True if the added pair of matrices should be used as expansion point
   subroutine add_to_modFIFO_disk(queue, F, D, energy, exppoint)
   implicit none
      type(modFIFO), target :: queue
      type(matrix), intent(in) :: F, D
      integer, pointer :: iFpointer, iDpointer  
      integer, pointer :: iDplaceinqueue, iFplaceinqueue
      logical, intent(in) :: exppoint
      logical :: addtoqueue, fileexists,OnMaster
      real(realk), intent(in) :: energy
      integer :: i, j, rowdim, coldim
      integer, target :: lunF, lunD
      character(len=80)  :: internalFile 
      character(len=8)   :: filenameD, filenameF !number allowed to be up to 3 digits

      rowdim = D%nrow
      coldim = D%ncol
      if (queue%qcounter == 200) queue%qcounter = 1

      !Construct filenames:
      !====================
      i = 0
      do
         i = i + 1
         if (i == 100) call lsquit('Something wrong in add_to_modFIFO_disk',-1)
         !Convert integer to character via internal file:
         write(internalFile,*) queue%qcounter
         !Left-justify the string:
         internalFile = adjustl(internalFile)
         !Construct filenames:
         filenameD = 'D' // trim(internalFile) // '.mat'
         filenameF = 'F' // trim(internalFile) // '.mat'
         !Check if filename already exists (possible if several queues are in use at the same time):
         INQUIRE(file=filenameD,EXIST=fileexists)
         if (fileexists) then
            queue%qcounter = queue%qcounter + 1
         else
            exit
         endif
      enddo

      !call mat_init(queue%Fmatrices(queue%qcounter),rowdim,coldim)
      !call mat_init(queue%Dmatrices(queue%qcounter),rowdim,coldim)
      lunF = -1 ; lunD = -1
      call lsopen(lunF,filenameF,'new','UNFORMATTED')
      call lsopen(lunD,filenameD,'new','UNFORMATTED')

      !Write matrices to disk:
      !=======================
      OnMaster = .TRUE.
      call mat_write_to_disk(lunF,F,OnMaster)
      call mat_write_to_disk(lunD,D,OnMaster)

      queue%iFmatrices(queue%qcounter) = lunF
      queue%iDmatrices(queue%qcounter) = lunD

      !Fpointer=>queue%Fmatrices(queue%qcounter)
      !Dpointer=>queue%Dmatrices(queue%qcounter)
      iFpointer=>queue%iFmatrices(queue%qcounter)
      iDpointer=>queue%iDmatrices(queue%qcounter)

      queue%qcounter = queue%qcounter + 1

      addtoqueue = .false.
      if (exppoint) then !if new expansion point, move existing expansion point to queue
         if(associated(queue%iD_exp)) then !if it is not the first addition, 
                                           !move current expansion point to queue
            iDplaceinqueue => queue%iD_exp
            iFplaceinqueue => queue%iF_exp
            addtoqueue = .true.
         endif
         queue%iD_exp => iDpointer
         queue%iF_exp => iFpointer
         queue%energy = energy
      else
         addtoqueue = .true.
         iDplaceinqueue => iDpointer
         iFplaceinqueue => iFpointer
      endif

      if (addtoqueue) then
         if (queue%offset < queue%queuesize-1) then
            queue%offset = queue%offset + 1
            queue%iDarray(queue%offset)%p => iDplaceinqueue
            queue%iFarray(queue%offset)%p => iFplaceinqueue
         else
            !call mat_free(queue%Darray(1)%p)
            !call mat_free(queue%Farray(1)%p)
            call lsclose(queue%iDarray(1)%p,'DELETE')
            call lsclose(queue%iFarray(1)%p,'DELETE')
            do i = 1, queue%queuesize-2
               queue%iDarray(i)%p => queue%iDarray(i+1)%p
               queue%iFarray(i)%p => queue%iFarray(i+1)%p
            enddo
            queue%iDarray(queue%offset)%p => iDplaceinqueue
            queue%iFarray(queue%offset)%p => iFplaceinqueue
         endif
      endif
   end subroutine add_to_modFIFO_disk

!> \brief Branch out to the proper remove-from-queue routine.
!> \author Stinne Host
!> \date May 2010
!> \param nremove Number of pairs of matrices (F,D) to remove
!> \param queue The type(modFIFO) to remove from
   subroutine remove_from_modFIFO(nremove, queue)
   implicit none
      integer, intent(in) :: nremove
      type(modFIFO), intent(inout) :: queue

      if (queue%disk) then
         call remove_from_modFIFO_disk(nremove, queue)
      else
         call remove_from_modFIFO_memory(nremove, queue)
      endif

   end subroutine remove_from_modFIFO

!> \brief Remove the nremove oldest pairs (F,D) in a type(modFIFO) on disk.
!> \author Stinne Host
!> \date May 2010
!> \param nremove Number of pairs of matrices (F,D) to remove
!> \param queue The type(modFIFO) to remove from
   subroutine remove_from_modFIFO_disk(nremove, queue)
   implicit none
      integer, intent(in) :: nremove
      type(modFIFO), intent(inout) :: queue
      integer :: i, nkeep

      nkeep = queue%offset - nremove !Number of matrices to keep
      if (nkeep < 0) then
         call lsquit('Cannot remove that many matrices from queue!',-1)
      endif

      !Free matrices which we want to discard
      do i = 1, nremove
         call lsclose(queue%iDarray(i)%p, 'DELETE')
         call lsclose(queue%iFarray(i)%p, 'DELETE')
         !call mat_free(queue%Darray(i)%p)
         !call mat_free(queue%Farray(i)%p)
      enddo

      !Move pointers
      do i = 1, nkeep
         queue%iDarray(i)%p => queue%iDarray(i+nremove)%p
         queue%iFarray(i)%p => queue%iFarray(i+nremove)%p
      enddo

      !Reset offset
      queue%offset = nkeep
   end subroutine remove_from_modFIFO_disk

!> \brief Remove the nremove oldest pairs (F,D) in a type(modFIFO) in memory.
!> \author Stinne Host
!> \date 2007
!> \param nremove Number of pairs of matrices (F,D) to remove
!> \param queue The type(modFIFO) to remove from
   subroutine remove_from_modFIFO_memory(nremove, queue)
   implicit none
      integer, intent(in) :: nremove
      type(modFIFO), intent(inout) :: queue
      integer :: i, nkeep

      nkeep = queue%offset - nremove !Number of matrices to keep
      if (nkeep < 0) then
         call lsquit('Cannot remove that many matrices from queue!',-1)
      endif

      !Free matrices which we want to discard
      do i = 1, nremove
         call mat_free(queue%Darray(i)%p)
         call mat_free(queue%Farray(i)%p)
      enddo

      !Move pointers
      do i = 1, nkeep
         queue%Darray(i)%p => queue%Darray(i+nremove)%p
         queue%Farray(i)%p => queue%Farray(i+nremove)%p
      enddo

      !Reset offset
      queue%offset = nkeep
   end subroutine remove_from_modFIFO_memory

!> \brief Change expansion point in a type(modFIFO). Discard old exp. point, don't touch the queue (mostly for debugging).
!> \author Stinne Host
!> \date 2007
!> \param queue The type(modFIFO) which should have a new expansion point
!> \param Fpointer Pointer to the new F expansion point
!> \param Dpointer Pointer to the new D expansion point
!> \param energy Energy for the new expansion point
   subroutine change_exp_modFIFO(queue, Fpointer, Dpointer,energy) 
   implicit none
      type(modFIFO) :: queue
      type(matrix), pointer :: Fpointer, Dpointer
      real(realk), intent(in) :: energy

      if (queue%disk) then
         call lsquit('change_exp_modFIFO not implemented for queue on disk',-1)
      else
         queue%D_exp => Dpointer
         queue%F_exp => Fpointer
         queue%energy = energy
      endif

   end subroutine change_exp_modFIFO

!> \brief Retrieve a pair of matrices (F,D) from a type(modFIFO). 
!> \author Stinne Host
!> \date May 2010
!> \param queue The type(modFIFO) to remove from
!> \param i The number of the desired pair in queue. If i=0, retrieve expansion point
!> \param Fpointer Pointer to retrieved F matrix
!> \param Dpointer Pointer to retrieved D matrix
   subroutine get_from_modFIFO(queue, i, Fpointer, Dpointer)
   !FIXME: It should be possible to keep the queue on disk
   implicit none
      type(modFIFO), intent(inout) :: queue
      type(matrix), pointer :: Fpointer, Dpointer
      integer, intent(in) :: i

      if (queue%disk) then
         call get_from_modFIFO_disk(queue, i, Fpointer, Dpointer)
      else
         call get_from_modFIFO_memory(queue, i, Fpointer, Dpointer)
      endif

   end subroutine get_from_modFIFO

!> \brief Retrieve a pair of matrices (F,D) from a type(modFIFO) on disk. 
!> \author Stinne Host
!> \date 2007
!> \param queue The type(modFIFO) to retrieve from
!> \param i The number of the desired pair in queue. If i=0, retrieve expansion point
!> \param Fpointer Pointer to retrieved F matrix
!> \param Dpointer Pointer to retrieved D matrix
   subroutine get_from_modFIFO_disk(queue, i, Fpointer, Dpointer)
   !FIXME: It should be possible to keep the queue on disk
   implicit none
      type(modFIFO), intent(inout), target :: queue
      type(matrix), pointer :: Fpointer, Dpointer
      integer, intent(in) :: i
      logical :: OnMaster
      OnMaster = .TRUE.
      if (i == 0) then !Get the expansion point
         rewind(queue%iF_exp)
         rewind(queue%iD_exp)
         call mat_read_from_disk(queue%iF_exp,queue%tempF,OnMaster)
         call mat_read_from_disk(queue%iD_exp,queue%tempD,OnMaster)
         Fpointer => queue%tempF
         Dpointer => queue%tempD
         !Fpointer => queue%F_exp
         !Dpointer => queue%D_exp
      else if (i >= 1 .and. i <= queue%queuesize-1) then
         rewind(queue%iFarray(i)%p)
         rewind(queue%iDarray(i)%p)
         call mat_read_from_disk(queue%iFarray(i)%p,queue%tempF,OnMaster)
         call mat_read_from_disk(queue%iDarray(i)%p,queue%tempD,OnMaster)
         Fpointer => queue%tempF
         Dpointer => queue%tempD
         !Fpointer => queue%Farray(i)%p
         !Dpointer => queue%Darray(i)%p
      else
         write(mat_lu,*) 'Index out of bounds in get_from_modFIFO:'
         write(mat_lu,*) 'Requested index:', i
         write(mat_lu,*) 'queuesize-1:', queue%queuesize-1
         CALL lsQUIT('Index out of bounds in get_from_modFIFO',-1)
      endif

   end subroutine get_from_modFIFO_disk

!> \brief Retrieve a pair of matrices (F,D) from a type(modFIFO) in memory. 
!> \author Stinne Host
!> \date 2007
!> \param queue The type(modFIFO) to retrieve from
!> \param i The number of the desired pair in queue. If i=0, retrieve expansion point
!> \param Fpointer Pointer to retrieved F matrix
!> \param Dpointer Pointer to retrieved D matrix
   subroutine get_from_modFIFO_memory(queue, i, Fpointer, Dpointer)
   !FIXME: It should be possible to keep the queue on disk
   implicit none
      type(modFIFO), intent(in) :: queue
      type(matrix), pointer :: Fpointer, Dpointer
      integer, intent(in) :: i

      if (i == 0) then !Get the expansion point
         Fpointer => queue%F_exp
         Dpointer => queue%D_exp
      else if (i >= 1 .and. i <= queue%queuesize-1) then
         Fpointer => queue%Farray(i)%p
         Dpointer => queue%Darray(i)%p
      else
         write(mat_lu,*) 'Index out of bounds in get_from_modFIFO:'
         write(mat_lu,*) 'Requested index:', i
         write(mat_lu,*) 'queuesize-1:', queue%queuesize-1
         CALL lsQUIT('Index out of bounds in get_from_modFIFO',-1)
      endif

   end subroutine get_from_modFIFO_memory

!> \brief Memory saving routine. Dump a type(modFIFO) to disk. 
!> 
!> The dumped type(modFIFO) is not free'd with free_mod_FIFO, but all the memory
!> is free'd. This is because the intention would usually be to use the type(modFIFO)
!> later, otherwise it could just be free'd instead of dumped.
!> 
!> \author Stinne Host
!> \date 2009
!> \param queue The type(modFIFO) to dump to disk
!> \param lun Logical unit number of file to dump to
!> \param ndim Number of basis functions/size of matrices
   SUBROUTINE fifoqueue_on_disk(queue,lun,ndim)
      implicit none
      type(modFIFO), intent(inout) :: queue
      integer, intent(out) :: lun,ndim
      integer :: i,j
      logical :: OnMaster
      OnMaster=.TRUE.
      if (queue%disk) then
         call lsquit('Programming error, fifoqueue_on_disk should not be called queue%disk=true',-1)
      endif
      lun = -1
      CALL lsOPEN(lun,'queue','unknown','UNFORMATTED')
      write(lun) queue%offset
      do i = 1, queue%offset
         ndim = queue%Farray(1)%p%nrow
         call mat_write_to_disk(lun,queue%Farray(i)%p,OnMaster)
         call mat_write_to_disk(lun,queue%Darray(i)%p,OnMaster)
      enddo
      call remove_from_modFIFO(queue%offset, queue)
   END SUBROUTINE fifoqueue_on_disk

!> \brief Memory saving routine. Read a type(modFIFO) to disk. 
!> 
!> The type(modFIFO) is not initialized with modFIFO_init. This should be done
!> before this routine is called. This is because the intention would usually be to read a type(modFIFO)
!> that has been used earlier in the program.
!> 
!> \author Stinne Host
!> \date 2009
!> \param queue The type(modFIFO) to read from disk
!> \param lun Logical unit number of file to read from
!> \param ndim Number of basis functions/size of matrices
   SUBROUTINE fifoqueue_from_disk(queue,lun,ndim)
      implicit none
      type(modFIFO), intent(inout) :: queue
      integer, intent(inout) :: lun,ndim
      type(matrix) :: D, F
      real(realk) :: energy
      integer :: i, qsize, offset
      logical :: OnMaster

      if (queue%disk) then
         call lsquit('Programming error, fifoqueue_from_disk should not be called queue%disk=true',-1)
      endif

      rewind(lun)
      read(lun) offset

      if (offset > 0) then
         call mat_init(D,ndim,ndim)
         call mat_init(F,ndim,ndim)
      endif
      OnMaster=.TRUE.
      do i = 1, offset
         call mat_read_from_disk(lun,F,OnMaster)
         call mat_read_from_disk(lun,D,OnMaster)
         call add_to_modFIFO(queue, F, D, 0.0E0_realk, .false.)
      enddo

      CALL lsCLOSE(lun,'DELETE')

      if (offset > 0) then
         call mat_free(D)
         call mat_free(F)
      endif
   END SUBROUTINE fifoqueue_from_disk

END MODULE queue_ops
