!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
!...  This file contains the data used to generate cube files.
!
!...  2012-03-09, Bin Gao
!...  * first version

#include "gen1int_host.h"

!> \brief data module of cube files
!> \author Bin Gao
!> \date 2012-03-09
module gen1int_cube

  implicit none

  logical, save, public :: do_density_cube = .false.                !electron density cube file
  logical, save, public :: do_homo_cube = .false.                   !HOMO cube file
  logical, save, public :: do_lumo_cube = .false.                   !LUMO cube file
  logical, save, public :: do_mo_cube = .false.                     !MO cube files
  integer, save, public :: num_cube_mo = 0                          !number of MOs to generate
  integer, save, allocatable, public :: idx_cube_mo(:)              !indices of MOs
  character(MAX_LEN_STR), save, public :: cube_format = "GAUSSIAN"  !format of cube file
  real(REALK), save, public :: cube_origin(3) = 0.0_8           !origin of cube file
  real(REALK), save, public :: cube_increment(3,3) = 0.0_8      !increments of cube file
  integer, save, public :: cube_num_inc(3) = 0                      !number of increments of cube file

end module gen1int_cube
