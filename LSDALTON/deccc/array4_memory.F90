!> @file
!> Memory manager for four dimensional arrays
module array4_memory_manager

  use precision
  use LSTIMING!,only:lstimer
  use memory_handling!, only: mem_alloc,mem_dealloc
  use dec_typedef_module

  !> Allocated memory
  real(realk) :: array4_alloc_memory = 0.0E0_realk
  !> Deallocated memory
  real(realk) :: array4_dealloc_memory = 0.0E0_realk
  !> Currently allocated memory
  real(realk) :: array4_memory_in_use = 0.0E0_realk
  !> Max allocated memory
  real(realk) :: array4_max_memory = 0.0E0_realk

public :: memory_allocate_4d
public :: memory_deallocate_4d
public :: print_memory_currents4
private
  contains

  !> \brief Allocate memory for 4d arrays with memory statistics
  !> \author Marcin Ziolkowski
  !> \param mat Four dimensional pointer to be allocated
  !> \param dims Dimensions of indeces
    subroutine memory_allocate_4d(mat,dims)
      implicit none
      real(realk), pointer :: mat(:,:,:,:)
      integer, dimension(4), intent(in) :: dims
      real(realk) :: vector_size
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      call memory_deallocate_4d(mat)
      vector_size = dble(dims(1))*dble(dims(2))*dble(dims(3))*dble(dims(4))*realk

      call mem_alloc(mat,dims(1),dims(2),dims(3),dims(4))

!$OMP CRITICAL
      array4_alloc_memory = array4_alloc_memory + vector_size
      array4_memory_in_use = array4_memory_in_use + vector_size
      array4_max_memory = max(array4_max_memory,array4_memory_in_use)
!$OMP END CRITICAL

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)

    end subroutine memory_allocate_4d

  !> \brief Deallocate memory for 4d arrays with memory statistics
  !> \author Marcin Ziolkowski
  !> \param mat Four dimensional pointer to be deallocated
    subroutine memory_deallocate_4d(mat)
      implicit none
      real(realk), pointer :: mat(:,:,:,:)
      real(realk) :: vector_size
      real(realk) :: dim1,dim2,dim3,dim4
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      if(associated(mat)) then
         dim1 = dble(size(mat(:,1,1,1)))
         dim2 = dble(size(mat(1,:,1,1)))
         dim3 = dble(size(mat(1,1,:,1)))
         dim4 = dble(size(mat(1,1,1,:)))
         vector_size = dim1*dim2*dim3*dim4*realk
         call mem_dealloc(mat)
         mat => null()
!$OMP CRITICAL
         array4_dealloc_memory = array4_dealloc_memory + vector_size
         array4_memory_in_use = array4_memory_in_use - vector_size
!$OMP END CRITICAL
      end if

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)


    end subroutine memory_deallocate_4d

  !> \brief Print statistics of array4 objects
  !> \author Marcin Ziolkowski
  !> \param output File unit for output 
    subroutine print_memory_statistics4(output)

      implicit none
      integer, intent(in) :: output

      write(DECinfo%output,'(/,a)')    '  Array4 memory statistics    '
      write(DECinfo%output,'(a)')      ' =================================================='
      write(DECinfo%output,'(a,f12.2,a)') ' Allocated memory    : ',array4_alloc_memory/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Deallocated memory  : ',array4_dealloc_memory/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Alloc-dealloc mem   : ', &
           (array4_alloc_memory-array4_dealloc_memory)/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Memory in use       : ',array4_memory_in_use/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Max memory in use   : ',array4_max_memory/2**20,' MB'

      write(DECinfo%output,'(a,/,a)') '  Time ', &
           ' ======='

    end subroutine print_memory_statistics4
  
  !> \brief Print currenly used memory and max allocated memory so far
  !> \author Marcin Ziolkowski
  !> \param output File unit for output
  subroutine print_memory_currents4(output)

    implicit none
    integer, intent(in) :: output
    write(DECinfo%output,'(a,g12.4,a)') ' Allocated memory for array4   :',&
         & array4_alloc_memory/(1.0E9_realk),' GB'
    write(DECinfo%output,'(a,g12.4,a)') ' Memory in use for array4      :',&
         & array4_memory_in_use/(1.0E9_realk),' GB'
    write(DECinfo%output,'(a,g12.4,a)') ' Max memory in use for array4  :',&
         & array4_max_memory/(1.0E9_realk),' GB'

  end subroutine print_memory_currents4

end module array4_memory_manager
