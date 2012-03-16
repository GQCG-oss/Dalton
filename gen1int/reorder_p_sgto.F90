!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
!...  This file reorders the p-shell contracted real solid-harmonic Gaussians
!...  in Dalton's order.
!
!...  2011-10-02, Bin Gao
!...  * first version

#include "xkind.h"

  !> \brief reorders the p-shell contracted real solid-harmonic Gaussians in Dalton's order
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param dim_bra_sgto is the dimension of SGTOs on bra center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_opt is the number of operators
  !> \param gen_ints contains the contracted integrals from Gen1Int
  subroutine reorder_p_sgto(dim_bra_sgto, num_contr_ket, num_opt, gen_ints)
    implicit none
    integer, intent(in) :: dim_bra_sgto
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(inout) :: gen_ints(dim_bra_sgto,3,num_contr_ket,num_opt)
    real(REALK), allocatable :: pshell_ints(:)  !temporary integrals
    integer icontr, iopt                        !incremental recorders
    integer ierr                                !error information
    ! Dalton's order of SGTOs: px(1), py(-1), pz(0),
    ! while those in Gen1Int is: py(-1), pz(0), px(1)
    allocate(pshell_ints(dim_bra_sgto), stat=ierr)
    if (ierr/=0) stop "reorder_p_sgto>> failed to allocate pshell_ints!"
    do iopt = 1, num_opt
      do icontr = 1, num_contr_ket
        pshell_ints = gen_ints(:,3,icontr,iopt)
        gen_ints(:,3,icontr,iopt) = gen_ints(:,2,icontr,iopt)
        gen_ints(:,2,icontr,iopt) = gen_ints(:,1,icontr,iopt)
        gen_ints(:,1,icontr,iopt) = pshell_ints
      end do
    end do
    deallocate(pshell_ints)
    return
  end subroutine reorder_p_sgto
