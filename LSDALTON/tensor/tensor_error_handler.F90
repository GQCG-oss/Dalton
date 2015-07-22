module tensor_error_handler
   use tensor_parameters_and_counters

   public :: tensor_status_quit
   private

   interface tensor_status_quit
      module procedure tensor_status_quit_std,&
                     & tensor_status_quit_long
   end interface tensor_status_quit

   contains

   subroutine tensor_status_quit_std(msg,stat)
      integer(kind=tensor_standard_int), intent(in) :: stat
      include "error_output.inc"
   end subroutine tensor_status_quit_std
   subroutine tensor_status_quit_long(msg,stat)
      integer(kind=tensor_long_int), intent(in) :: stat
      include "error_output.inc"
   end subroutine tensor_status_quit_long

end module tensor_error_handler
