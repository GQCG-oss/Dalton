!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_cfg
! this module holds all options for the embedding code, and routines to set/access
! this information

   type fde_export
      character(len=4) :: ex_level
      character(len=3) :: ex_format
      logical          :: ex_n
      logical          :: ex_gn
      logical          :: ex_coulomb
      logical          :: do_grid_out
   end type
   
   type fde_import
      character(len=4) :: im_level
      character(len=3) :: im_format
      logical          :: im_n
      logical          :: im_gn
      logical          :: im_coulomb
      logical          :: im_vemb
      logical          :: im_frozen
      logical          :: im_update_vemb
   end type

   type fde_file
      character(len=60) :: name
      integer           :: unit
   end type
   
   type fde_files
      type(fde_file) :: embpot
      type(fde_file) :: frozen
      type(fde_file) :: export
      type(fde_file) :: test
      type(fde_file) :: logfile
   end type

   public fde_export
   public fde_import
   public fde_file
   public fde_files

   public fde_initialize_cfg
   
   public fde_set_print_level
   public fde_set_import_info
   public fde_set_export_info
   public fde_set_files_info
   public fde_set_qccode_fileops
   public fde_set_screening
   public fde_set_parallel

   public fde_get_print_level
   public fde_get_import_info
   public fde_get_export_info
   public fde_get_files_info
   public fde_get_qccode_fileops
   public fde_get_screening
   public fde_get_parallel
   
   private

   integer, save      :: fde_print_level   = 0
   real(kind=8), save :: fde_cfg_screening = 1.0d-9
   logical, save      :: fde_is_parallel   = .false.
   logical, save      :: fde_use_qccode_fileops = .true.

   type(fde_export), save :: fde_ex_control
   type(fde_import), save :: fde_im_control
   type(fde_files),  save :: fde_f

   contains
   
!   ----------------------------------------------------------------------------
      subroutine fde_set_screening(screening_value)
!   ----------------------------------------------------------------------------
         real(kind=8), intent(in) :: screening_value
         fde_cfg_screening = screening_value
      end subroutine fde_set_screening

!   ----------------------------------------------------------------------------
      subroutine fde_get_screening(screening_value)
!   ----------------------------------------------------------------------------
         real(kind=8), intent(out) :: screening_value
         screening_value = fde_cfg_screening 
      end subroutine fde_get_screening

!   ----------------------------------------------------------------------------
      subroutine fde_set_qccode_fileops(fops_value)
!   ----------------------------------------------------------------------------
         logical, intent(in) :: fops_value 
         fde_use_qccode_fileops = fops_value
      end subroutine fde_set_qccode_fileops

!   ----------------------------------------------------------------------------
      subroutine fde_get_qccode_fileops(fops_value)
!   ----------------------------------------------------------------------------
         logical, intent(out) :: fops_value
         fops_value = fde_use_qccode_fileops
      end subroutine fde_get_qccode_fileops

!   ----------------------------------------------------------------------------
      subroutine fde_set_parallel(par_value)
!   ----------------------------------------------------------------------------
         logical, intent(in) :: par_value
         fde_is_parallel = par_value
      end subroutine fde_set_parallel

!   ----------------------------------------------------------------------------
      subroutine fde_get_parallel(par_value)
!   ----------------------------------------------------------------------------
         logical, intent(out) :: par_value
         par_value = fde_is_parallel 
      end subroutine fde_get_parallel

!   ----------------------------------------------------------------------------
      subroutine fde_set_print_level(level)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: level
         fde_print_level = level
      end subroutine fde_set_print_level

!   ----------------------------------------------------------------------------
      subroutine fde_get_print_level(level)
!   ----------------------------------------------------------------------------
         integer, intent(out) :: level
         level = fde_print_level 
      end subroutine fde_get_print_level

!   ----------------------------------------------------------------------------
      subroutine fde_initialize_cfg
!   ----------------------------------------------------------------------------
         call fde_initialize_files(fde_f)
         call fde_initialize_export(fde_ex_control)
         call fde_initialize_import(fde_im_control)
      end subroutine fde_initialize_cfg
         
!   ----------------------------------------------------------------------------
      subroutine fde_initialize_files(f)
!   ----------------------------------------------------------------------------
         type(fde_files) :: f
         
! file units are initialized to -1, and should be modified to point to a proper
! unit in each case. codes such as DALTON have routines to open/close files and 
! dynamically attribute unit, other codes do not...
         f%embpot%name = 'EMBPOT'
         f%embpot%unit = -1

         f%frozen%name = 'FRZDNS'
         f%frozen%unit = -1

         f%export%name = 'GRIDOUT'
         f%export%unit = -1

         f%test%name = 'TESTFILE'
         f%test%unit = -1

! this is the output file name, which is set to embed.out but can be modified whenever
! another evaluator code (DALTON, DIRAC etc) is used, through a call to the fde_CODE_init()  
! routine
         f%logfile%name = 'EMBED.OUT'
         f%logfile%unit = -1         
      end subroutine fde_initialize_files
      
!   ----------------------------------------------------------------------------
      subroutine fde_initialize_export(ec)
!   ----------------------------------------------------------------------------
         type(fde_export) :: ec
         
         ec%ex_level     = 'DHF'
         ec%ex_format    = 'XML'
         ec%ex_n         = .false.
         ec%ex_gn        = .false.
         ec%ex_coulomb   = .false.
         ec%do_grid_out  = .false.
      end subroutine fde_initialize_export

!   ----------------------------------------------------------------------------
      subroutine fde_initialize_import(ic)
!   ----------------------------------------------------------------------------
         type(fde_import) :: ic
        
         ic%im_level     = 'DHF'
         ic%im_format    = 'TXT'
         ic%im_n         = .false.
         ic%im_gn        = .false.
         ic%im_coulomb   = .false.
         ic%im_vemb      = .false.
         ic%im_frozen    = .false.
         ic%im_update_vemb=.false.
      end subroutine fde_initialize_import

!   ----------------------------------------------------------------------------
      subroutine fde_set_export_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_export), intent(in) :: tmp
         
         fde_ex_control%ex_level     = tmp%ex_level
         fde_ex_control%ex_format    = tmp%ex_format
         fde_ex_control%ex_n         = tmp%ex_n
         fde_ex_control%ex_gn        = tmp%ex_gn 
         fde_ex_control%ex_coulomb   = tmp%ex_coulomb
         fde_ex_control%do_grid_out  = tmp%do_grid_out
      end subroutine fde_set_export_info

!   ----------------------------------------------------------------------------
      subroutine fde_set_import_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_import), intent(in) :: tmp
         
         fde_im_control%im_level     = tmp%im_level
         fde_im_control%im_format    = tmp%im_format
         fde_im_control%im_n         = tmp%im_n
         fde_im_control%im_gn        = tmp%im_gn 
         fde_im_control%im_coulomb   = tmp%im_coulomb
         fde_im_control%im_vemb      = tmp%im_vemb
         fde_im_control%im_frozen    = tmp%im_frozen
         fde_im_control%im_update_vemb=tmp%im_update_vemb
      end subroutine fde_set_import_info

!   ----------------------------------------------------------------------------
      subroutine fde_set_files_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_files), intent(in) :: tmp
         
         fde_f%embpot%name = tmp%embpot%name
         fde_f%embpot%unit = tmp%embpot%unit

         fde_f%frozen%name = tmp%frozen%name
         fde_f%frozen%unit = tmp%frozen%unit

         fde_f%export%name = tmp%export%name
         fde_f%export%unit = tmp%export%unit

         fde_f%test%name = tmp%test%name
         fde_f%test%unit = tmp%test%unit 

         fde_f%logfile%name = tmp%logfile%name
         fde_f%logfile%unit = tmp%logfile%unit         
      end subroutine fde_set_files_info
      
!   ----------------------------------------------------------------------------
      subroutine fde_get_export_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_export), intent(out) :: tmp
         
         tmp%ex_level   = fde_ex_control%ex_level
         tmp%ex_format  = fde_ex_control%ex_format
         tmp%ex_n       = fde_ex_control%ex_n         
         tmp%ex_gn      = fde_ex_control%ex_gn           
         tmp%ex_coulomb = fde_ex_control%ex_coulomb
         tmp%do_grid_out= fde_ex_control%do_grid_out    
      end subroutine fde_get_export_info

!   ----------------------------------------------------------------------------
      subroutine fde_get_import_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_import), intent(out) :: tmp
         
         tmp%im_level       = fde_im_control%im_level       
         tmp%im_format      = fde_im_control%im_format      
         tmp%im_n           = fde_im_control%im_n           
         tmp%im_gn          = fde_im_control%im_gn           
         tmp%im_coulomb     = fde_im_control%im_coulomb     
         tmp%im_vemb        = fde_im_control%im_vemb        
         tmp%im_frozen      = fde_im_control%im_frozen      
         tmp%im_update_vemb = fde_im_control%im_update_vemb 
      end subroutine fde_get_import_info
      
!   ----------------------------------------------------------------------------
      subroutine fde_get_files_info(tmp)
!   ----------------------------------------------------------------------------
         type(fde_files), intent(out) :: tmp
         
         tmp%embpot%name = fde_f%embpot%name
         tmp%embpot%unit = fde_f%embpot%unit

         tmp%frozen%name = fde_f%frozen%name
         tmp%frozen%unit = fde_f%frozen%unit

         tmp%export%name = fde_f%export%name
         tmp%export%unit = fde_f%export%unit

         tmp%test%name = fde_f%test%name
         tmp%test%unit = fde_f%test%unit 

         tmp%logfile%name = fde_f%logfile%name
         tmp%logfile%unit = fde_f%logfile%unit         
      end subroutine fde_get_files_info
      

end module fde_cfg
