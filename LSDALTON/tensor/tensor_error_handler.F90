module tensor_error_handler
   use tensor_parameters_and_counters
   contains

   subroutine tensor_status_quit(msg,stat)
      character*(*),intent(in) :: msg
      integer, intent(in) :: stat
      print *, "Tensor status quit has been called with:",stat
      stop -1
   end subroutine tensor_status_quit

end module tensor_error_handler
