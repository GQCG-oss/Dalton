!> @file
!> Memory manager for three dimensional arrays
module array3_memory_manager

  use precision
  use LSTIMING!,only:lstimer
  use memory_handling!, only: mem_alloc, mem_dealloc
  use dec_typedef_module

  !> Allocated memory
  real(realk) :: array3_allocated_memory = 0.0E0_realk
  !> Deallocated memory
  real(realk) :: array3_deallocated_memory = 0.0E0_realk
  !> Currently allocated memory
  real(realk) :: array3_memory_in_use = 0.0E0_realk
  !> Max allocated memory
  real(realk) :: array3_max_memory = 0.0E0_realk

public :: memory_allocate_3d
public :: memory_deallocate_3d
public :: print_memory_currents_3d
private
  contains

  !> \brief Allocate memory for 3d arrays with memory statistics
  !> \author Janus Eriksen (after 4d template by Marcin Ziolkowski)
  !> \param mat Three dimensional pointer to be allocated
  !> \param dims Dimensions of indices
    subroutine memory_allocate_3d(mat,dims)
      implicit none
      real(realk), pointer :: mat(:,:,:)
      integer, dimension(3), intent(in) :: dims
      real(realk) :: vector_size
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      call memory_deallocate_3d(mat)
      vector_size = dble(dims(1))*dble(dims(2))*dble(dims(3))*realk
      call mem_alloc(mat,dims(1),dims(2),dims(3))

!$OMP CRITICAL
      array3_allocated_memory = array3_allocated_memory + vector_size
      array3_memory_in_use = array3_memory_in_use + vector_size
      array3_max_memory = max(array3_max_memory,array3_memory_in_use)
!$OMP END CRITICAL

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)

    end subroutine memory_allocate_3d

  !> \brief Deallocate memory for 3d arrays with memory statistics
  !> \author Janus Eriksen (after 4d template by Marcin Ziolkowski)
  !> \param mat Three dimensional pointer to be deallocated
    subroutine memory_deallocate_3d(mat)
      implicit none
      real(realk), pointer :: mat(:,:,:)
      real(realk) :: vector_size
      real(realk) :: dim1,dim2,dim3
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      if(associated(mat)) then
         dim1 = dble(size(mat(:,1,1)))
         dim2 = dble(size(mat(1,:,1)))
         dim3 = dble(size(mat(1,1,:)))
         vector_size = dim1*dim2*dim3*realk
         call mem_dealloc(mat)
         mat => null()
!$OMP CRITICAL
         array3_deallocated_memory = array3_deallocated_memory + vector_size
         array3_memory_in_use = array3_memory_in_use - vector_size
!$OMP END CRITICAL
      end if

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)


    end subroutine memory_deallocate_3d

  !> \brief Print statistics of array3 objects
  !> \author Marcin Ziolkowski (modified by Janus Eriksen)
  !> \param output File unit for output 
    subroutine print_memory_statistics_3d(output)

      implicit none
      integer, intent(in) :: output

      write(DECinfo%output,'(/,a)')    '  Array3 memory statistics    '
      write(DECinfo%output,'(a)')      ' =================================================='
      write(DECinfo%output,'(a,f12.2,a)') ' Allocated memory    : ',array3_allocated_memory/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Deallocated memory  : ',array3_deallocated_memory/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Alloc-dealloc mem   : ', &
           (array3_allocated_memory-array3_deallocated_memory)/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Memory in use       : ',array3_memory_in_use/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Max memory in use   : ',array3_max_memory/2**20,' MB'

      write(DECinfo%output,'(a,/,a)') '  Time ', &
           ' ======='

    end subroutine print_memory_statistics_3d

  !> \brief Print currenly used memory and max allocated memory so far (for array3)
  !> \author Marcin Ziolkowski (modified by Janus Eriksen)
  !> \param output File unit for output
    subroutine print_memory_currents_3d(output)
  
      implicit none
      integer, intent(in) :: output
      write(DECinfo%output,'(a,g12.4,a)') ' Allocated memory for array3   :',&
           & array3_allocated_memory/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Memory in use for array3      :',&
           & array3_memory_in_use/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Max memory in use for array3  :',&
           & array3_max_memory/(1.0E9_realk),' GB'
  
    end subroutine print_memory_currents_3d

end module array3_memory_manager
