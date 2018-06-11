!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_input

   use fde_types
   use fde_cfg
   use fde_data
! aspg, 23/06/2017
! before the interface to xcfun/xcint is properly set up, anything having to do
! with functional derivatives will be disabled
!  use fde_nadd_derv
!  use fde_evaluators
   use fde_evaluators_dalton
   
   private

   public fde_input_init
   public fde_input_validate_and_postprocess
   public fde_input_show_setup

   public skip_all_imports
   
   logical, save :: skip_all_imports = .false.

   contains

!   ----------------------------------------------------------------------------
      subroutine fde_input_init(log_unit,log_name)
!   ----------------------------------------------------------------------------
         integer          :: log_unit
         character(len=*) :: log_name
         type(fde_files)  :: ftmp

         call fde_initialize_cfg
         
         call fde_get_files_info(ftmp)
         ftmp%logfile%unit = log_unit
         ftmp%logfile%name = log_name
         call fde_set_files_info(ftmp)
         
         call fde_set_qccode_fileops(.true.)
         call fde_qccode_data_interface
      end subroutine fde_input_init

      
!   ----------------------------------------------------------------------------
      subroutine fde_input_validate_and_postprocess
!   ----------------------------------------------------------------------------
         type(fde_export)  :: etmp
         type(fde_import)  :: itmp
         type(fde_files)   :: ftmp
         integer :: pl

      call fde_get_files_info(ftmp)
      call fde_get_export_info(etmp)
      call fde_get_import_info(itmp)

      if (.not.skip_all_imports) then

      if (itmp%im_vemb) then
         if (itmp%im_update_vemb) then
            write (ftmp%logfile%unit,*) '.UPDATE and .EMBPOT cannot be used together!'
            call fde_quit('conflicting keywords selected')
         endif
      else
         itmp%im_update_vemb = .true.
         itmp%im_frozen      = .true.
      endif
 
! we should arrive here if the update of the embedding potential during the
! scf is asked for, or if fde response        
      if (itmp%im_frozen) then
         itmp%im_n        = .true.
         itmp%im_gn       = .true.
         itmp%im_coulomb  = .true.
      endif

      if (etmp%do_grid_out) then
         etmp%ex_n        = .true.
         etmp%ex_coulomb  = .true.
         etmp%ex_gn       = .true.
      endif

      endif ! .not.skip_all_imports

      if (itmp%im_vemb) then
         call fde_import_static      
      endif

! aspg, 23/06/2017
! before the interface to xcfun/xcint is properly set up, anything having to do
! with functional derivatives will be disabled
#if 0
      if (itmp%im_frozen) then
         call fde_import_frozen
         call fde_initialize_nadd_functionals         
      endif
#endif

      if (etmp%do_grid_out) then
         call fde_import_gridout
      endif
      
! saving information on 
      call fde_set_files_info(ftmp)
      call fde_set_export_info(etmp)
      call fde_set_import_info(itmp)

     end subroutine  


      subroutine fde_input_show_setup
         integer :: unit
         character(len=60) :: string
         type(fde_export)  :: etmp
         type(fde_import)  :: itmp
         type(fde_files)   :: ftmp
         integer :: pl

         call fde_get_print_level(pl)

         call fde_get_files_info(ftmp)
         call fde_get_export_info(etmp)
         call fde_get_import_info(itmp)

         unit = ftmp%logfile%unit
         
         write(unit,'(/A,I5)')  ' * FDE print level               : ',pl

! aspg, 23/06/2017
! before the interface to xcfun/xcint is properly set up, anything having to do
! with functional derivatives will be disabled
#if 0
         if (itmp%im_frozen) call fde_print_nadd_functionals(unit)
#endif

         if (itmp%im_vemb) then
            string = ftmp%embpot%name
            
            write(unit,'(A,A60)')  ' * FDE Potential read from file : ',string
            write(unit,'(/3X,2A,/3X,2A/)') &
            'Enviroment effects included via the ',          &
            'fixed potential method described in:',          &
            'A.S.P. Gomes, C. R. Jacob and L. Visscher, ', &
            'PCCP 10 (2008) 5353-5362.'

         else
            write(unit,'(A/)')  ' * FDE Potential generated from frozen, active densities'
            call fde_test_frozen_density
         endif

         if (itmp%im_frozen) then
            string = ftmp%frozen%name
            write(unit,'(A,A60)') &
                ' * FDE Density (and gradient) from frozen subsystems read from file: ',string
         endif

         if (etmp%do_grid_out) then
            string = ftmp%export%name
            write(unit,'(2A)') &
            ' * FDE Gridfile with updated density written to:',string

            select case(etmp%ex_level)
               case('DHF','MP2')
                  write(unit,*)                               &
                  'Outputted density will be taken from: ',trim(etmp%ex_level)
                  
               case default
                  write(unit,*) 'Input for outputted density not a &
     &  x  valid calculation, HF will be used instead. Given:',trim(etmp%ex_level)
                  etmp%ex_level = 'DHF'
                  
            end select
         endif
                  
      end subroutine fde_input_show_setup

end module fde_input
