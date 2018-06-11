!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_mod

   use fde_types
   use fde_cfg
   use fde_data
   use fde_io
   use fde_export_data
! still need to sort out xcfun interface
!  use fde_evaluators
!  for the same problem we need to specifically add
   use fde_evaluators_dalton
   use fde_input
   use fde_input_dalton

end module fde_mod
