!@file This file only contains a workaround for an issue with the
!current cray compiler, i.e.
! 1) if too much memory is allocated using the
!    simple operator "=" does not work. This module hides the actual size of
!    the arrays from the compiler use VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN

module dec_workarounds_module
  use precision

#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
  abstract interface
    subroutine subblock_op(drain,source,nel,scal1,scal2)
    import
    implicit none
    integer(kind=8), intent(in) :: nel
    real(realk), intent(in), optional :: scal1,scal2
    real(realk), intent(inout) :: drain(nel)
    real(realk), intent(in) :: source(nel)
    end subroutine subblock_op
  end interface
#endif
#endif


  contains
  
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
  subroutine assign_in_subblocks(drain,op,source,nel,scal1,scal2,gpu)
    implicit none
    integer(kind=8), intent(in) :: nel
    character, intent(in) :: op
    real(realk),intent(inout) :: drain(:)
    real(realk),intent(in)    :: source(:)
    real(realk), intent(in), optional :: scal1,scal2
    logical, intent(in), optional :: gpu
    integer(kind=8) :: block,i
    integer, parameter :: k = 50000000

    if (present(gpu) .and. (gpu)) then
       gpu = .true.
    else
       gpu = .false.
    endif
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
    procedure(subblock_op), pointer :: op_blocks

    select case(op)
    case("+")
      op_blocks => add_subblock
    case("=")
      op_blocks => copy_subblock
    case("-")
      op_blocks => subtract_subblock
    case default
      call lsquit("ERROR(in_subblocks): wrong choice of op",-1)
    end select


    do i = 1, nel, k
      block = k
      if( (nel-i)<k .and. mod(nel-i+1,k)/=0 ) block=mod(nel,k)
      !if(infpar%mynum == 0 ) print *,"begin",i,"end",i+block-1,"length",block,"of",nel
      call op_blocks(drain(i:i+block-1),source(i:i+block-1),block,gpu,scal1=scal1,scal2=scal2)
    enddo
#else
    print *,"ERROR(assign_in_subblocks): you cannot use this workaround&
    & without a compiler that is able to understand fortran 2003 function pointers"
    stop 1
#endif
  end subroutine assign_in_subblocks

  subroutine copy_subblock(drain,source,nel,gpu,scal1,scal2)
    implicit none
    integer(kind=8), intent(in) :: nel
    logical, intent(in) :: gpu
    real(realk), intent(in), optional :: scal1,scal2
    real(realk), intent(inout) :: drain(nel)
    real(realk), intent(in) :: source(nel)
    if(present(scal1))then
      print *,"ERROR(copy_subblock):wrong input, the scal2 only exists for being&
      & able to use function pointers"
      stop 1
    endif

    if (gpu) then

       if(present(scal2))then
!$acc kernels present(drain,source)
         drain = scal2 * source
!$acc end kernels
       else
!$acc kernels present(drain,source)
         drain = source
!$acc end kernels
       endif

    else

       if(present(scal2))then
         !OMP WORKSHARE
         drain = scal2 * source
         !OMP WORKSHARE
       else
         drain = source
       endif

    endif

  end subroutine copy_subblock
  subroutine add_subblock(drain,source,nel,gpu,scal1,scal2)
    implicit none
    integer(kind=8), intent(in) :: nel
    logical, intent(in) :: gpu
    real(realk), intent(in), optional :: scal1,scal2
    real(realk), intent(inout) :: drain(nel)
    real(realk), intent(in) :: source(nel)

    if (gpu) then

       if(present(scal1).and.present(scal2))then
!$acc kernels present(drain,source)
         drain = scal1 * drain + scal2 * source
!$acc end kernels
       elseif(present(scal1))then
!$acc kernels present(drain,source)
         drain = scal1 * drain + source
!$acc end kernels
       elseif(present(scal2))then
!$acc kernels present(drain,source)
         drain = drain + scal2 * source
!$acc end kernels
       else
!$acc kernels present(drain,source)
         drain = drain + source
!$acc end kernels
       endif

    else

       if(present(scal1).and.present(scal2))then
         drain = scal1 * drain + scal2 * source
       elseif(present(scal1))then
         drain = scal1 * drain + source
       elseif(present(scal2))then
         drain = drain + scal2 * source
       else
         drain = drain + source
       endif

    endif
  end subroutine add_subblock
  subroutine subtract_subblock(drain,source,nel,gpu,scal1,scal2)
    implicit none
    integer(kind=8), intent(in) :: nel
    logical, intent(in) :: gpu
    real(realk), intent(in), optional :: scal1,scal2
    real(realk), intent(inout) :: drain(nel)
    real(realk), intent(in) :: source(nel)

    if (gpu) then

       if(present(scal1).and.present(scal2))then
!$acc kernels present(drain,source)
         drain = scal1 * drain - scal2 * source
!$acc end kernels
       elseif(present(scal1))then
!$acc kernels present(drain,source)
         drain = scal1 * drain - source
!$acc end kernels
       elseif(present(scal2))then
!$acc kernels present(drain,source)
         drain = drain - scal2 * source
!$acc end kernels
       else
!$acc kernels present(drain,source)
         drain = drain - source
!$acc end kernels
       endif

    else

       if(present(scal1).and.present(scal2))then
         drain = scal1 * drain - scal2 * source
       elseif(present(scal1))then
         drain = scal1 * drain - source
       elseif(present(scal2))then
         drain = drain - scal2 * source
       else
         drain = drain - source
       endif

    endif

  end subroutine subtract_subblock
#endif
  subroutine dec_workarounds_dummy
    implicit none
  end subroutine dec_workarounds_dummy
end module dec_workarounds_module
