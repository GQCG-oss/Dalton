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
!...  This file initializes the data used in Gen1Int interface before
!...  any calculation.
!
!...  2011-10-03, Bin Gao
!...  * rewrites in Fortran 90 format
!
!...  2010-12-06, Bin Gao
!...  * first version

#include "xkind.h"

  !> \brief initializes the data used in Gen1Int interface before any calculation
  !> \author Bin Gao
  !> \date 2010-12-06
  subroutine gen1int_dal_init
    ! AO sub-shells
    use gen1int_shell
#ifdef BUILD_GEN1INT
#include "implicit.h"
#ifdef GEN1INT_DEBUG
#include "priunit.h"
#endif
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
!-#include "nuclei.h"
!-#include "onecom.h"
#include "aovec.h"
!-#include "symmet.h"
!-#include "symind.h"
#include "primit.h"
!-#include "lmns.h"
#include "shells.h"
#include "orgcom.h"
!-#include "ibtfun.h"
!-#include "dummy.h"
    integer ishell, jshell       !incremental recorders over sub-shells
    integer strt_prim, end_prim  !incremental recorders of start and end indices of primitive GTOs
    integer end_contr            !incremental recorder of end indices of contractions
    real(REALK), allocatable :: contr_coef(:,:)
                                 !contraction coefficients
    integer icontr               !incremental recorder over contractions
    integer ierr                 !error information
    if (shell_init) return
    call QENTER("gen1int_dal_init")
    ! gets the number of AO sub-shells
    num_ao_shells = 1
    do ishell = 1, KMAX-1
      ! found a new AO sub-shell
      if (JSTRT(ishell)/=JSTRT(ishell+1)) num_ao_shells = num_ao_shells+1
    end do
    ! initializes the AO sub-shells
    allocate(ao_shells(num_ao_shells), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_dal_init>> failed to allocate ao_shells!")
    ! the first AO sub-shells
    strt_prim = JSTRT(1)+1
    end_prim = JSTRT(1)+NUCO(1)
    end_contr = NUMCF(1)+NRCO(1)-1
    allocate(contr_coef(NUMCF(1):end_contr,strt_prim:end_prim), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_dal_init>> failed to allocate contr_coef!")
    do icontr = NUMCF(1), end_contr
      contr_coef(icontr,strt_prim:end_prim) = PRICCF(strt_prim:end_prim,icontr)
    end do
    call gen1int_shell_set(is_sgto=SPHR(1), idx_cent=NCENT(1),            &
           coord_cent=CENT(1,1:3,1), ang_num=NHKT(1)-1, num_prim=NUCO(1), &
           exponents=PRIEXP(strt_prim:end_prim), num_contr=NRCO(1),       &
           contr_coef=contr_coef, ao_shell=ao_shells(1))
    deallocate(contr_coef)
    ! other AO sub-shells
    kshell = 1
    do ishell = 1, KMAX-1
      ! found a new AO sub-shell
      if (JSTRT(ishell)/=JSTRT(ishell+1)) then
        jshell = ishell+1
        strt_prim = JSTRT(jshell)+1
        end_prim = JSTRT(jshell)+NUCO(jshell)
        end_contr = NUMCF(jshell)+NRCO(jshell)-1
        allocate(contr_coef(NUMCF(jshell):end_contr,strt_prim:end_prim), stat=ierr)
        if (ierr/=0) call QUIT("gen1int_dal_init>> failed to allocate contr_coef!")
        do icontr = NUMCF(jshell), end_contr
          contr_coef(icontr,strt_prim:end_prim) = PRICCF(strt_prim:end_prim,icontr)
        end do
        kshell = kshell+1
        call gen1int_shell_set(is_sgto=SPHR(jshell), idx_cent=NCENT(jshell), &
               coord_cent=CENT(jshell,1:3,1), ang_num=NHKT(jshell)-1,        &
               num_prim=NUCO(jshell), exponents=PRIEXP(strt_prim:end_prim),  &
               num_contr=NRCO(jshell), contr_coef=contr_coef, ao_shell=ao_shells(kshell))
        deallocate(contr_coef)
      end if
    end do
#ifdef GEN1INT_DEBUG
    call gen1int_shell_dump(num_ao_shells, ao_shells, LUPRI)
#endif
    shell_init = .true.
    call QEXIT("gen1int_dal_init")
    return
#else
    call QUIT("Gen1Int is not installed!")
#endif
  end subroutine gen1int_dal_init
