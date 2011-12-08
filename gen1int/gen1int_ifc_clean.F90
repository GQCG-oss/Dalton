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
!...  This file cleans up the data used in Gen1Int interface after all calculations.
!
!...  2011-10-02, Bin Gao
!...  * first version

  !> \brief cleans up the data used in Gen1Int interface after all calculations
  !> \author Bin Gao
  !> \date 2011-10-02
  subroutine gen1int_ifc_clean
    ! AO sub-shells
    use gen1int_shell
    implicit none
#ifdef BUILD_GEN1INT
    if (.not. shell_init) return
    call QENTER("gen1int_ifc_clean")
    ! cleans up the AO sub-shells
    call gen1int_shell_clean(num_ao_shells, ao_shells)
    deallocate(ao_shells)
    num_ao_shells = 0
    call QEXIT("gen1int_ifc_clean")
    return
#else
    call QUIT("gen1int_ifc_clean>> Gen1Int is not installed!")
#endif
  end subroutine gen1int_ifc_clean
