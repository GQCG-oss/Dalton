      SUBROUTINE fde_quit(text)
         character(len=60) :: text
! using the dalton quit, to handle parallel executaion etc
         call quit(text)
      END

