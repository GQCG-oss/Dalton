module ls_pcm_write
   
implicit none

public init_host_writer
public host_writer

private

integer :: global_print_unit = -1

contains

subroutine init_host_writer(print_unit)

   integer :: print_unit

   if (global_print_unit .eq. -1) then
     global_print_unit = print_unit
   else
     write(*, *) "Printing to stdout"
   end if

end subroutine init_host_writer                

subroutine host_writer(message, message_length) bind(c, name='host_writer')
    
   use, intrinsic :: iso_c_binding

   integer(c_size_t), intent(in) :: message_length
   character(kind=c_char)        :: message(message_length)

   if (global_print_unit .eq. -1) then
     write(*, '(1000A)') message
   else
     write(global_print_unit, '(1000A)') message
     flush(global_print_unit)
   end if

end subroutine host_writer
   
end module
