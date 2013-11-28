MODULE queue_ops4
use queue4_module
use Matrix_Operations

contains

   SUBROUTINE modFIFO_init4(queue,qsize)
   !FIXME: It should be possible to keep the queue on disk
      implicit none
      TYPE(modFIFO4) :: queue
      integer, intent(in) :: qsize 
      integer :: i
         queue%queuesize = qsize
         queue%offset = 0  
         allocate(queue%xarray(qsize-1))
         allocate(queue%parray(qsize-1))
         allocate(queue%Aparray(qsize-1))
         allocate(queue%Apparray(qsize-1))
         do i = 1, qsize-1
            nullify(queue%xarray(i)%p)
            nullify(queue%parray(i)%p)
            nullify(queue%Aparray(i)%p)
            nullify(queue%Apparray(i)%p)
         enddo
         nullify(queue%x_exp)
         nullify(queue%p_exp)
         nullify(queue%Ap_exp)
         nullify(queue%App_exp)
         !allocate(queue%bmatrices(200))
         !allocate(queue%Ematrices(200))
         !allocate(queue%Smatrices(200))
         queue%qcounter = 1
   END SUBROUTINE modFIFO_init4

   SUBROUTINE modFIFO_free4(queue)
      implicit none
      TYPE(modFIFO4) :: queue
      integer :: i
         if(associated(queue%x_exp)) call mat_free(queue%x_exp)
         if(associated(queue%p_exp)) call mat_free(queue%p_exp)
         if(associated(queue%Ap_exp)) call mat_free(queue%Ap_exp)
         if(associated(queue%App_exp)) call mat_free(queue%App_exp)

         do i = 1, queue%offset
            call mat_free(queue%xarray(i)%p)
            call mat_free(queue%parray(i)%p)
            call mat_free(queue%Aparray(i)%p)
            call mat_free(queue%Apparray(i)%p)
         enddo

         deallocate(queue%xarray)
         deallocate(queue%parray)
         deallocate(queue%Aparray)
         deallocate(queue%Apparray)
   END SUBROUTINE modFIFO_free4

   subroutine add_to_modFIFO4(queue, x, p, Ap, App, energy, exppoint)
   implicit none
      type(modFIFO4), target :: queue
      type(matrix), intent(in) :: x, p, Ap, App
      type(matrix), pointer :: xpointer, ppointer, Appointer, Apppointer  
      type(matrix), pointer :: xplaceinqueue,pplaceinqueue,Applaceinqueue, Appplaceinqueue
      logical, intent(in) :: exppoint
      logical :: addtoqueue
      real(realk), intent(in) :: energy
      integer :: i, ndim

      ndim = p%nrow
      if (queue%qcounter == 500) queue%qcounter = 1

      call mat_init(queue%xmatrices(queue%qcounter),ndim,ndim)
      call mat_init(queue%pmatrices(queue%qcounter),ndim,ndim)
      call mat_init(queue%Apmatrices(queue%qcounter),ndim,ndim)
      call mat_init(queue%Appmatrices(queue%qcounter),ndim,ndim)
      queue%xmatrices(queue%qcounter) = x
      queue%pmatrices(queue%qcounter) = p
      queue%Apmatrices(queue%qcounter) = Ap
      queue%Appmatrices(queue%qcounter) = App

      xpointer=>queue%xmatrices(queue%qcounter)
      ppointer=>queue%pmatrices(queue%qcounter)
      Appointer=>queue%Apmatrices(queue%qcounter)
      Apppointer=>queue%Appmatrices(queue%qcounter)

      queue%qcounter = queue%qcounter + 1

      addtoqueue = .false.
      if (exppoint) then !if new expansion point, move existing expansion point to queue
         if(associated(queue%x_exp)) then !if it is not the first addition, 
                                          !move current expansion point to queue
            xplaceinqueue => queue%x_exp
            pplaceinqueue => queue%p_exp
            Applaceinqueue => queue%Ap_exp
            Appplaceinqueue => queue%App_exp
            addtoqueue = .true.
         endif
         queue%x_exp => xpointer
         queue%p_exp => ppointer
         queue%Ap_exp => Appointer
         queue%App_exp => Apppointer
         queue%energy = energy
      else
         addtoqueue = .true.
         xplaceinqueue => xpointer
         pplaceinqueue => ppointer
         Applaceinqueue => Appointer
         Appplaceinqueue => Apppointer
      endif

      if (addtoqueue) then
         if (queue%offset < queue%queuesize-1) then
            queue%offset = queue%offset + 1
            queue%xarray(queue%offset)%p => xplaceinqueue
            queue%parray(queue%offset)%p => pplaceinqueue
            queue%Aparray(queue%offset)%p => Applaceinqueue
            queue%Apparray(queue%offset)%p => Appplaceinqueue
         else
            call mat_free(queue%xarray(1)%p)
            call mat_free(queue%parray(1)%p)
            call mat_free(queue%Aparray(1)%p)
            call mat_free(queue%Apparray(1)%p)
            do i = 1, queue%queuesize-2
               queue%xarray(i)%p => queue%xarray(i+1)%p
               queue%parray(i)%p => queue%parray(i+1)%p
               queue%Aparray(i)%p => queue%Aparray(i+1)%p
               queue%Apparray(i)%p => queue%Apparray(i+1)%p
            enddo
            queue%xarray(queue%offset)%p => xplaceinqueue
            queue%parray(queue%offset)%p => pplaceinqueue
            queue%Aparray(queue%offset)%p => Applaceinqueue
            queue%Apparray(queue%offset)%p => Appplaceinqueue
         endif
      endif
   end subroutine add_to_modFIFO4

   subroutine remove_from_modFIFO4(nremove, queue)
   !Remove nremove oldest densities in queue.
   implicit none
      integer, intent(in) :: nremove
      type(modFIFO4) :: queue
      !type(matrix), pointer :: bpointer, Epointer, Spointer
      !type(matrix), pointer :: bplaceinqueue, Eplaceinqueue, Splaceinqueue
      !logical, intent(in) :: exppoint
      !logical :: addtoqueue
      !real(realk), intent(in) :: energy
      integer :: i, ndim, nkeep
                                                                                                                                                      
      ndim = queue%x_exp%nrow
      nkeep = queue%offset - nremove !Number of matrices to keep

      !Free matrices which we want to discard
      do i = 1, nremove
         call mat_free(queue%xarray(i)%p)
         call mat_free(queue%parray(i)%p)
         call mat_free(queue%Aparray(i)%p)
         call mat_free(queue%Apparray(i)%p)
      enddo

      !Move pointers
      do i = 1, nkeep
         queue%xarray(i)%p => queue%xarray(i+nremove)%p
         queue%parray(i)%p => queue%parray(i+nremove)%p
         queue%Aparray(i)%p => queue%Aparray(i+nremove)%p
         queue%Apparray(i)%p => queue%Apparray(i+nremove)%p
      enddo

      !Reset offset
      queue%offset = nkeep
   end subroutine remove_from_modFIFO4

   subroutine change_exp_modFIFO4(queue, xpointer, ppointer, Appointer,Apppointer,energy) !For debug: Don't touch queue, change only exp point
   implicit none
      type(modFIFO4) :: queue
      type(matrix), pointer :: xpointer, ppointer, Appointer, Apppointer
      real(realk), intent(in) :: energy

      queue%x_exp => xpointer
      queue%p_exp => ppointer
      queue%Ap_exp => Appointer
      queue%App_exp => Apppointer
      queue%energy = energy

   end subroutine change_exp_modFIFO4

   subroutine get_from_modFIFO4(queue, i, xpointer, ppointer, Appointer, Apppointer)
   !FIXME: It should be possible to keep the queue on disk
   implicit none
      type(modFIFO4), intent(in) :: queue
      type(matrix), pointer :: xpointer,ppointer, Appointer, Apppointer
      integer, intent(in) :: i

      if (i == 0) then !Get the expansion point
         xpointer => queue%x_exp
         ppointer => queue%p_exp
         Appointer => queue%Ap_exp
         Apppointer => queue%App_exp
      else if (i >= 1 .and. i <= queue%queuesize-1) then
         xpointer => queue%xarray(i)%p
         ppointer => queue%parray(i)%p
         Appointer => queue%Aparray(i)%p
         Apppointer => queue%Apparray(i)%p
      else
         write(mat_lu,*) 'Index out of bounds in get_from_modFIFO4:'
         write(mat_lu,*) 'Requested index:', i
         write(mat_lu,*) 'queuesize-1:', queue%queuesize-1
         CALL LSQUIT('Index out of bounds in get_from_modFIFO4', mat_lu)
      endif

   end subroutine get_from_modFIFO4

END MODULE queue_ops4
