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
!...  2011-10-06, Bin Gao
!...  * rewrites based on \fn(ORBPRO) subroutine, getting the unnormalized
!...    contraction coefficients
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
  !> \param num_atom_type is the number of atomic types
  !> \param num_sym_atom contains the number of symmetry independent centers of atomic types
  !> \param ang_numbers contains the angular momentum (1=s, 2=p, 3=d, ...)
  !> \param num_cgto contains the number of CGTOs in the AO blocks for an angular momentum
  !> \param num_prim contains the number of uncontracted functions
  !> \param num_contr contains the number of contracted functions
  !> \param exponents contains the exponents of primitive shells
  !> \param ucontr_coefs contains the unnormalized contraction coefficients
  subroutine gen1int_ifc_init(num_atom_type, KATOM, num_sym_atom, ang_numbers,    &
                              NBLCK, KANG, num_cgto, KBLOCK, num_prim, num_contr, &
                              KPRIM, exponents, ucontr_coefs)
    ! AO sub-shells
    use gen1int_shell
    implicit none
    integer, intent(in) :: num_atom_type
    integer, intent(in) :: KATOM
    integer, intent(in) :: num_sym_atom(KATOM)
    integer, intent(in) :: ang_numbers(KATOM)
    integer, intent(in) :: NBLCK(KATOM)
    integer, intent(in) :: KANG
    integer, intent(in) :: num_cgto(KANG,KATOM)
    integer, intent(in) :: KBLOCK
    integer, intent(in) :: num_prim(KBLOCK)
    integer, intent(in) :: num_contr(KBLOCK)
    integer, intent(in) :: KPRIM
    real(REALK), intent(in) :: exponents(KPRIM,KBLOCK)
    real(REALK), intent(in) :: ucontr_coefs(KPRIM,KPRIM,KBLOCK)
#ifdef BUILD_GEN1INT
#include "mxcent.h"
#include "maxaqn.h"
#include "ccom.h"
#include "nuclei.h"
#ifdef GEN1INT_DEBUG
#include "priunit.h"
#endif
    integer IDX_CENT   !index of symmetry independent center
    integer IDX_BLOCK  !
    integer ITYP       !incremental recorder over number of atomic types
    integer ICENT      !incremental recorder over number of symmetry independent centers
    integer IANG       !incremental recorder over angular momentum
    integer KBCH       !
    logical is_sgto    !if SGTOs
    integer ishell     !incremental recorder over AO sub-shells
    integer ang_num    !angular number
    real(REALK), allocatable :: contr_coef(:,:)  !contraction coefficients
    integer icontr, iprim                        !incremental recorder over contractions
    integer ierr                                 !error information
    ! backs if already initialized the basis sets
    if (shell_init) return
    call QENTER("gen1int_ifc_init")
    ! gets the number of AO sub-shells and kind of GTOs
    num_ao_shells = 0
    is_sgto = .false.
    ! number of atomic types
    do ITYP = 1, num_atom_type
      ! number of symmetry independent centers of this type
      do ICENT = 1, num_sym_atom(ITYP)
        ! angular momentum 1=s, 2=p, 3=d, etc.
        do IANG = 1, ang_numbers(ITYP)
          num_ao_shells = num_ao_shells+1
          if (SPH(IANG)) is_sgto = .true.
        end do
       end do
    end do
    ! initializes the AO sub-shells
    allocate(ao_shells(num_ao_shells), stat=ierr)
    if (ierr/=0) call QUIT("gen1int_ifc_init>> failed to allocate ao_shells!")
    IDX_CENT = 0
    IDX_BLOCK = 0
    ishell = 0
    ! number of atomic types
    do ITYP = 1, num_atom_type
      ! number of symmetry independent centers of this type
      do ICENT = 1, num_sym_atom(ITYP)
        IDX_CENT = IDX_CENT+1
        KBCH = IDX_BLOCK
        ! angular momentum 1=s, 2=p, 3=d, etc.
        do IANG = 1, ang_numbers(ITYP)
          ! next block
          KBCH = KBCH+1
          ! gets the contraction coefficients
          allocate(contr_coef(num_contr(KBCH),num_prim(KBCH)), stat=ierr)
          if (ierr/=0) call QUIT("gen1int_ifc_init>> failed to allocate contr_coef!")
          do iprim = 1, num_prim(KBCH)
            do icontr = 1, num_contr(KBCH)
              contr_coef(icontr,iprim) = ucontr_coefs(iprim,icontr,KBCH)
            end do
          end do
          ! normalizes the contraction coefficients
          ang_num = IANG-1
          !-if (SPH(IANG)) then
          if (is_sgto) then
            call norm_contr_sgto(ang_num, num_prim(KBCH),          &
                                 exponents(1:num_prim(KBCH),KBCH), &
                                 num_contr(KBCH), contr_coef)
          else
            call norm_contr_cgto(ang_num, num_prim(KBCH),          &
                                 exponents(1:num_prim(KBCH),KBCH), &
                                 num_contr(KBCH), contr_coef)
          end if
          !-ISTBNU(IDX_CENT)  !stabiliser: basic sym. op. that do not move center
          ishell = ishell+1
          if (ishell>1) then
            !-call gen1int_shell_set(is_sgto=SPH(IANG), idx_cent=IDX_CENT,  &
            call gen1int_shell_set(is_sgto=is_sgto, idx_cent=IDX_CENT, &
                   coord_cent=CORD(1:3,IDX_CENT), ang_num=ang_num,     &
                   num_prim=num_prim(KBCH),                            &
                   exponents=exponents(1:num_prim(KBCH),KBCH),         &
                   num_contr=num_contr(KBCH), contr_coef=contr_coef,   &
                   last_shell=ao_shells(ishell-1), ao_shell=ao_shells(ishell))
          ! sets the first AO sub-shell
          else
            !-call gen1int_shell_set(is_sgto=SPH(IANG), idx_cent=IDX_CENT,  &
            call gen1int_shell_set(is_sgto=is_sgto, idx_cent=IDX_CENT, &
                   coord_cent=CORD(1:3,IDX_CENT), ang_num=ang_num,     &
                   num_prim=num_prim(KBCH),                            &
                   exponents=exponents(1:num_prim(KBCH),KBCH),         &
                   num_contr=num_contr(KBCH), contr_coef=contr_coef,   &
                   ao_shell=ao_shells(ishell))
          end if
          deallocate(contr_coef)
          ! skips other CGTOs in this AO block for this angular momentum
          KBCH = KBCH+num_cgto(IANG,ITYP)-1
        end do
      end do
      IDX_BLOCK = IDX_BLOCK+NBLCK(ITYP)
    end do
#ifdef GEN1INT_DEBUG
    call gen1int_shell_dump(num_ao_shells, ao_shells, LUPRI)
#endif
    shell_init = .true.
    call QEXIT("gen1int_ifc_init")
    return
#else
    call QUIT("gen1int_ifc_init>> Gen1Int is not installed!")
#endif
  end subroutine gen1int_ifc_init
