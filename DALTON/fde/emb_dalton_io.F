c modify here for using GPOPEN and GPCLOSE
c the optional arguments there should use a 1-length character string

      SUBROUTINE fde_qccode_file_open(file_unit,file_name)

         integer, intent(in) :: file_unit
         character(*), intent(in) :: file_name

         if (file_unit.eq.-1) then
             call gpopen(file_unit,file_name,'UNKNOWN','SEQUENTIAL',  
     &                   'FORMATTED',0,.FALSE.)
         else
             call fde_quit('unknown action in fde_qccode_file_open')
         end if
      END

      SUBROUTINE fde_qccode_file_close(file_unit)
         integer, intent(in) :: file_unit

         call GPCLOSE(file_unit,'KEEP')
      END

