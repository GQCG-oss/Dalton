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
!...  This file takes care the atomic orbital (AO) sub-shells used in Gen1Int interface.
!
!...  2011-10-02, Bin Gao
!...  * first version

#include "stdout.h"
#include "xkind.h"

!> \brief defines the AO sub-shells used in Gen1Int interface and corresponding subroutines
!> \author Bin Gao
!> \date 2011-10-02
module gen1int_shell
  implicit none
  ! contracted Cartesian or spherical GTOs, or an atomic orbital (AO) sub-shell
  ! here, an AO sub-shell contains the orbitals with the same exponents, and the
  ! same/different contraction coefficients; for instance 2p and 3p sub-shells
  type, public :: sub_shell_t
    private
    ! if spherical GTOs
    logical is_sgto
    ! index of atomic center where this AO sub-shell locates
    integer idx_cent
    ! coordinates of the atomic center
    real(REALK) coord_cent(3)
    ! angular number
    integer ang_num
    ! number of primitive GTOs
    integer num_prim
    ! exponents of primitive GTOs
    real(REALK), allocatable :: exponents(:)
    ! number of contractions
    integer num_contr
    ! contraction coefficients
    real(REALK), allocatable :: contr_coef(:,:)
    ! number of atomic orbitals
    integer num_ao
    !-Dalton uses the same sequence of SGTOs as those in Gen1Int except for p-shell
    !-which could be done by hand coding
    !-! magnetic numbers if spherical GTOs
    !-integer, allocatable :: mag_num(:)
    ! Cartesian powers if Cartesian GTOs
    integer, allocatable :: powers(:,:)
    ! base index of the atomic orbitals in this sub-shell, i.e., the index of the
    ! first atomic orbital in this shell minus 1
    integer base_idx
  end type sub_shell_t
  ! contracted integrals between two sub-shells
  type, public :: shell_int_t
    ! base index of the atomic orbitals on bra center
    integer base_idx_bra
    ! base index of the atomic orbitals on ket center
    integer base_idx_ket
    ! number of orbitals on bra center
    integer num_orb_bra
    ! number of orbitals on ket center
    integer num_orb_ket
    ! number of operators
    integer num_opt
    ! number of total geometric derivatives
    integer num_geo_derv
    ! contracted integrals between sub-shells on bra and ket centers
    real(REALK), allocatable :: contr_ints(:,:,:,:)
  end type shell_int_t

  ! for Dalton
  !
  ! if the AO sub-shells from Dalton initialized
  logical, public, save :: shell_init = .false.
  ! number of AO sub-shells from Dalton
  integer, public, save :: num_ao_shells = 0
  ! AO sub-shells from Dalton
  type(sub_shell_t), public, allocatable, save :: ao_shells(:)

  public :: gen1int_shell_set
  public :: gen1int_shell_dims
  public :: gen1int_shell_idx_orb
  public :: gen1int_shell_dump
  public :: gen1int_shell_clean
  public :: gen1int_shell_prop
  public :: gen1int_shell_int_clean
  public :: gen1int_shell_int_tri_diag
  public :: gen1int_shell_int_tri_off
  public :: gen1int_shell_int_square
  public :: gen1int_shell_tr_tri_diag
  public :: gen1int_shell_tr_tri_off
  public :: gen1int_shell_tr_square

  private :: gen1int_shell_carmom
  private :: gen1int_shell_nucpot
  ! for Dalton p-shell spherical GTOs
  private :: dal_reorder_p_sgto

  contains

  !> \brief sets an AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param is_sgto indicates if spherical GTOs
  !> \param idx_cent is the index of atomic center where this AO sub-shell locates
  !> \param coord_cent contains the coordinates of the atomic center
  !> \param ang_num is the angular number
  !> \param num_prim is the number of primitive GTOs
  !> \param exponents contains the exponents of primitive GTOs
  !> \param num_contr is the number of contractions
  !> \param contr_coef contains the contraction coefficients
  !> \param powers contains the Cartesian powers if Cartesian GTOs
  !> \param last_shell is the last sub-shell before this AO sub-shell
  !> \return ao_shell is the initialized AO sub-shell
  subroutine gen1int_shell_set(is_sgto, idx_cent, coord_cent, ang_num,     &
                               num_prim, exponents, num_contr, contr_coef, &
                               powers, last_shell, ao_shell)
    implicit none
    logical, intent(in) :: is_sgto
    integer, intent(in) :: idx_cent
    real(REALK), intent(in) :: coord_cent(3)
    integer, intent(in) :: ang_num
    integer, intent(in) :: num_prim
    real(REALK), intent(in) :: exponents(num_prim)
    integer, intent(in) :: num_contr
    real(REALK), intent(in) :: contr_coef(num_contr,num_prim)
    integer, optional, intent(in) :: powers(3,(ang_num+1)*(ang_num+2)/2)
    type(sub_shell_t), optional, intent(in) :: last_shell
    type(sub_shell_t), intent(inout) :: ao_shell
    integer iao         !incremental recorder
    integer xpow, ypow  !incremental recorders over xy powers
    integer ierr        !error information
    ao_shell%is_sgto = is_sgto
    ao_shell%idx_cent = idx_cent
    ao_shell%coord_cent = coord_cent
    ao_shell%ang_num = ang_num
    ao_shell%num_prim = num_prim
    allocate(ao_shell%exponents(num_prim), stat=ierr)
    if (ierr/=0) stop "gen1int_shell_set>> failed to allocate exponents!"
    ao_shell%exponents = exponents
    ao_shell%num_contr = num_contr
    allocate(ao_shell%contr_coef(num_contr,num_prim), stat=ierr)
    if (ierr/=0) stop "gen1int_shell_set>> failed to allocate contr_coef!"
    ao_shell%contr_coef = contr_coef
    ! spherical GTOs
    if (is_sgto) then
      ao_shell%num_ao = 2*ang_num+1
      !-allocate(ao_shell%mag_num(ao_shell%num_ao), stat=ierr)
      !-if (ierr/=0) stop "gen1int_shell_set>> failed to allocate mag_num!"
      !-if (present(mag_num)) then
      !-  ao_shell%mag_num = mag_num
      !-else
      !-  ! Dalton's order of SGTOs
      !-  if (ang_num==1) then
      !-    ao_shell%mag_num(1) = 1   !px
      !-    ao_shell%mag_num(2) = -1  !py
      !-    ao_shell%mag_num(3) = 0   !pz
      !-  else
      !-    do iao = -ang_num, ang_num
      !-      ao_shell%mag_num(ang_num+iao+1) = iao
      !-    end do
      !-  end if
      !-end if
    ! Cartesian GTOs
    else
      ao_shell%num_ao = (ang_num+1)*(ang_num+2)/2
      allocate(ao_shell%powers(3,ao_shell%num_ao), stat=ierr)
      if (ierr/=0) stop "gen1int_shell_set>> failed to allocate powers!"
      if (present(powers)) then
        ao_shell%powers = powers
      ! Dalton's order of CGTOs, for instance, dxx, dxy, dxz, dyy, dyz, dzz
      else
        iao = 0
        do xpow = ang_num, 0, -1
          do ypow = ang_num-xpow, 0, -1
            iao = iao+1
            ao_shell%powers(1,iao) = xpow
            ao_shell%powers(2,iao) = ypow
            ao_shell%powers(3,iao) = ang_num-(xpow+ypow)
          end do
        end do
      end if  
    end if
    if (present(last_shell)) then
      ao_shell%base_idx = last_shell%base_idx &
                        + last_shell%num_ao*last_shell%num_contr
    ! this is the first sub-shell
    else
      ao_shell%base_idx = 0
    end if
    return
  end subroutine gen1int_shell_set

  !> \brief gets dimensions of an AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param ao_shell is the AO sub-shell
  !> \return num_ao is the number of atomic orbitals
  !> \return num_contr is the number of contractions
  subroutine gen1int_shell_dims(ao_shell, num_ao, num_contr)
    implicit none
    type(sub_shell_t), intent(in) :: ao_shell
    integer, intent(out) :: num_ao
    integer, intent(out) :: num_contr
    num_ao = ao_shell%num_ao
    num_contr = ao_shell%num_contr
    return
  end subroutine gen1int_shell_dims

  !> \brief gets the index of the first and last orbtials for the given AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-06
  !> \param ao_shell is the AO sub-shell
  !> \return idx_first is the index of the first orbital
  !> \return idx_last is the index of the last orbital
  subroutine gen1int_shell_idx_orb(ao_shell, idx_first, idx_last)
    implicit none
    type(sub_shell_t), intent(in) :: ao_shell
    integer, intent(out) :: idx_first
    integer, intent(out) :: idx_last
    idx_first = ao_shell%base_idx+1
    idx_last = ao_shell%base_idx+ao_shell%num_ao*ao_shell%num_contr
    return
  end subroutine gen1int_shell_idx_orb

  !> \brief dumps information of several AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param num_shell is the number of AO sub-shells
  !> \param ao_shells are the AO sub-shells to dump
  !> \param io_std is the IO unit to dump
  subroutine gen1int_shell_dump(num_shell, ao_shells, io_std)
    implicit none
    integer, intent(in) :: num_shell
    type(sub_shell_t), intent(in) :: ao_shells(*)
    integer, intent(in) :: io_std
    integer ishell, iprim  !incremental recorder
    write(io_std,100) "number of AO sub-shells", num_shell
    do ishell = 1, num_shell
      write(io_std,100) "sub-shell", ishell
      write(io_std,100) "index of atomic center of the sub-shell", &
                         ao_shells(ishell)%idx_cent
      write(io_std,110) "coordinates of the atomic center", &
                         ao_shells(ishell)%coord_cent
      write(io_std,100) "angular number", ao_shells(ishell)%ang_num
      write(io_std,100) "number of primitive GTOs", ao_shells(ishell)%num_prim
      write(io_std,100) "number of contractions", ao_shells(ishell)%num_contr
      write(io_std,100) "    exponents      coefficients"
      do iprim = 1, ao_shells(ishell)%num_prim
        write(io_std,120)                     &
          ao_shells(ishell)%exponents(iprim), &
          ao_shells(ishell)%contr_coef(1:min(4,ao_shells(ishell)%num_contr),iprim)
        write(io_std,130) &
          ao_shells(ishell)%contr_coef(5:ao_shells(ishell)%num_contr,iprim)
      end do
      if (ao_shells(ishell)%is_sgto) then
        write(io_std,100) "SGTOs used"
        !-write(io_std,140) ao_shells(ishell)%mag_num
      else
        write(io_std,100) "CGTOs used"
        write(io_std,150) ao_shells(ishell)%powers
      end if
      write(io_std,100) "base index of the orbitals", ao_shells(ishell)%base_idx
    end do
    return
100 format("gen1int_shell_dump>> "A,2I8)
110 format("gen1int_shell_dump>> "A,3F16.8)
120 format("gen1int_shell_dump>> ",5Es16.8)
130 format("gen1int_shell_dump>> ",16X,4Es16.8)
140 format("gen1int_shell_dump>> magnetic numbers>> ",8I4)
150 format("gen1int_shell_dump>> Cartesian powers>> ",4(3I4,2X))
  end subroutine gen1int_shell_dump

  !> \brief cleans up several AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param num_shell is the number of AO sub-shells
  !> \param ao_shells are the AO sub-shells to clean up
  subroutine gen1int_shell_clean(num_shell, ao_shells)
    implicit none
    integer, intent(in) :: num_shell
    type(sub_shell_t), intent(inout) :: ao_shells(*)
    integer ishell  !incremental recorder
    do ishell = 1, num_shell
      deallocate(ao_shells(ishell)%exponents)
      deallocate(ao_shells(ishell)%contr_coef)
      if (ao_shells(ishell)%is_sgto) then
        !-deallocate(ao_shells(ishell)%mag_num)
      else
        deallocate(ao_shells(ishell)%powers)
      end if
    end do
    return
  end subroutine gen1int_shell_clean

  !> \brief calculates contracted property integrals between two AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-07
  !> \param prop_name is the name of property integrals to calculated
  !> \param bra_shell is the sub-shell on bra center
  !> \param ket_shell is the sub-shell on ket center
  !> \param is_lao indicates if using London atomic orbitals
  !> \param order_mom is the order of Cartesian multipole moments, only used for
  !>        integrals "CARMOM"
  !> \param num_cents is the number of geometric differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \param num_geo_derv is the number of total geometric derivatives
  !> \return prop_int contains the contracted property integrals between two AO sub-shells
  subroutine gen1int_shell_prop(prop_name, bra_shell, ket_shell, is_lao,    &
                                order_mom, num_cents, idx_cent, order_cent, &
                                num_geo_derv, prop_int)
    implicit none
    character*(*), intent(in) :: prop_name
    type(sub_shell_t), intent(in) :: bra_shell
    type(sub_shell_t), intent(in) :: ket_shell
    logical, intent(in) :: is_lao
    integer, intent(in) :: order_mom
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: num_geo_derv
    type(shell_int_t), intent(inout) :: prop_int
    ! origins in Dalton
#include "orgcom.h"
    ! uses \var(MXCENT)
#include "mxcent.h"
    ! coordinates and charges of atoms
#include "nuclei.h"
    integer, parameter :: IDX_DIPORG = -1        !index of dipole origin
    integer icent                                !incremental recorder over atomic centers
    real(REALK), allocatable :: tmp_ints(:,:,:)  !temporary contracted integrals
    integer ierr                                 !error information
    ! sets the base index of orbitals on bra and ket centers
    prop_int%base_idx_bra = bra_shell%base_idx
    prop_int%base_idx_ket = ket_shell%base_idx
    ! sets the dimensions of contracted integrals between two sub-shells
    prop_int%num_orb_bra = bra_shell%num_ao*bra_shell%num_contr
    prop_int%num_orb_ket = ket_shell%num_ao*ket_shell%num_contr
    prop_int%num_geo_derv = num_geo_derv
    ! different property integrals, please use alphabetical order!
    ! the following keywords and \var(prop_int%num_opt) should be
    ! consistent with subroutine \fn(gen1int_prop_attr)!!
    select case(trim(prop_name))
    ! Cartesian multipole integrals
    case("CARMOM")
      prop_int%num_opt = (order_mom+1)*(order_mom+2)/2
      allocate(prop_int%contr_ints(prop_int%num_orb_bra,prop_int%num_orb_ket, &
                                   prop_int%num_opt,prop_int%num_geo_derv), stat=ierr)
      if (ierr/=0) stop "gen1int_shell_prop>> failed to allocate contr_ints!"
      call gen1int_shell_carmom(bra_shell, ket_shell, is_lao, 0,          &
                                IDX_DIPORG, DIPORG, 1.0_REALK, order_mom, &
                                0, 0, 0, 0, 0, 0, 0, 0, 0,                &
                                num_cents, idx_cent, order_cent,          &
                                bra_shell%num_ao, bra_shell%num_contr,    &
                                ket_shell%num_ao, ket_shell%num_contr,    &
                                prop_int%num_opt*prop_int%num_geo_derv,   &
                                prop_int%contr_ints)
    ! dipole length integrals
    case("DIPLEN")
      prop_int%num_opt = 3
      allocate(prop_int%contr_ints(prop_int%num_orb_bra,prop_int%num_orb_ket, &
                                   prop_int%num_opt,prop_int%num_geo_derv), stat=ierr)
      if (ierr/=0) stop "gen1int_shell_prop>> failed to allocate contr_ints!"
      call gen1int_shell_carmom(bra_shell, ket_shell, is_lao, 0,        &
                                IDX_DIPORG, DIPORG, 1.0_REALK, 1,       &
                                0, 0, 0, 0, 0, 0, 0, 0, 0,              &
                                num_cents, idx_cent, order_cent,        &
                                bra_shell%num_ao, bra_shell%num_contr,  &
                                ket_shell%num_ao, ket_shell%num_contr,  &
                                prop_int%num_opt*prop_int%num_geo_derv, &
                                prop_int%contr_ints)
    ! kinetic energy integrals
    case("KINENERG")
      prop_int%num_opt = 6
      allocate(prop_int%contr_ints(prop_int%num_orb_bra,prop_int%num_orb_ket, &
                                   prop_int%num_opt,prop_int%num_geo_derv), stat=ierr)
      if (ierr/=0) stop "gen1int_shell_prop>> failed to allocate contr_ints!"
      call gen1int_shell_carmom(bra_shell, ket_shell, is_lao, 2,        &
                                IDX_DIPORG, DIPORG, -0.5_REALK, 0,      &
                                0, 0, 0, 0, 0, 0, 0, 0, 0,              &
                                num_cents, idx_cent, order_cent,        &
                                bra_shell%num_ao, bra_shell%num_contr,  &
                                ket_shell%num_ao, ket_shell%num_contr,  &
                                prop_int%num_opt*prop_int%num_geo_derv, &
                                prop_int%contr_ints)
      ! sums xx, yy and zz components
      prop_int%contr_ints(:,:,1,:) = prop_int%contr_ints(:,:,1,:) &
                                   + prop_int%contr_ints(:,:,3,:) &
                                   + prop_int%contr_ints(:,:,6,:)
      prop_int%num_opt = 1
    ! overlap integrals
    case("OVERLAP")
      prop_int%num_opt = 1
      allocate(prop_int%contr_ints(prop_int%num_orb_bra,prop_int%num_orb_ket, &
                                   prop_int%num_opt,prop_int%num_geo_derv), stat=ierr)
      if (ierr/=0) stop "gen1int_shell_prop>> failed to allocate contr_ints!"
      call gen1int_shell_carmom(bra_shell, ket_shell, is_lao, 0,        &
                                IDX_DIPORG, DIPORG, 1.0_REALK, 0,       &
                                0, 0, 0, 0, 0, 0, 0, 0, 0,              &
                                num_cents, idx_cent, order_cent,        &
                                bra_shell%num_ao, bra_shell%num_contr,  &
                                ket_shell%num_ao, ket_shell%num_contr,  &
                                prop_int%num_opt*prop_int%num_geo_derv, &
                                prop_int%contr_ints)
    ! one-electron potential energy integrals
    case("POTENERG")
      prop_int%num_opt = 1
      allocate(prop_int%contr_ints(prop_int%num_orb_bra,prop_int%num_orb_ket, &
                                   prop_int%num_opt,prop_int%num_geo_derv), stat=ierr)
      if (ierr/=0) stop "gen1int_shell_prop>> failed to allocate contr_ints!"
      ! the first atomic center
      call gen1int_shell_nucpot(bra_shell, ket_shell, is_lao, 0,       &
                                1, CORD(:,1), IDX_DIPORG,              &
                                DIPORG, CHARGE(1), 0,                  &
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          &
                                num_cents, idx_cent, order_cent,       &
                                bra_shell%num_ao, bra_shell%num_contr, &
                                ket_shell%num_ao, ket_shell%num_contr, &
                                prop_int%num_geo_derv, prop_int%contr_ints)
      ! allocates memory for temporary contracted integrals
      allocate(tmp_ints(prop_int%num_orb_bra,prop_int%num_orb_ket, &
                        prop_int%num_geo_derv), stat=ierr)
      if (ierr/=0) stop "gen1int_shell_prop>> failed to allocate tmp_ints!"
      ! other atomic centers
      do icent = 2, NUCDEP
        call gen1int_shell_nucpot(bra_shell, ket_shell, is_lao, 0,       &
                                  icent, CORD(:,icent),                  &
                                  IDX_DIPORG, DIPORG, CHARGE(icent), 0,  &
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,          &
                                  num_cents, idx_cent, order_cent,       &
                                  bra_shell%num_ao, bra_shell%num_contr, &
                                  ket_shell%num_ao, ket_shell%num_contr, &
                                  prop_int%num_geo_derv, tmp_ints)
        prop_int%contr_ints(:,:,1,:) = prop_int%contr_ints(:,:,1,:)+tmp_ints
      end do
      deallocate(tmp_ints)
    case default
      write(STDOUT,999) trim(prop_name)
      stop
    end select
    return
999 format("gen1int_shell_prop>> ",A," is not implemented!")
  end subroutine gen1int_shell_prop

  !> \brief cleans the contracted property integrals between two AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-07
  !> \param prop_int contains the contracted property integrals between two AO sub-shells
  subroutine gen1int_shell_int_clean(prop_int)
    implicit none
    type(shell_int_t), intent(inout) :: prop_int    
    deallocate(prop_int%contr_ints)
    return
  end subroutine gen1int_shell_int_clean

  !> \brief gets the contracted property integrals between two AO sub-shells back,
  !>        and writes the integrals on file is required, the integrals are belong
  !>        to the diagonal parts in triangular format
  !> \author Bin Gao
  !> \date 2011-10-07
  !> \param prop_int contains the contracted property integrals between two AO sub-shells
  !> \param base_geo_derv is the base address of current total geometric derivatives in the
  !>        list of all total geometric derivatives
  !> \param is_triang indicates if returning integral matrices are in triangular format,
  !>        otherwise square format
  !> \param num_ao_orb is the total number of atomic orbitals
  !> \param wrt_int indicates if writing integrals on file
  !> \param dim_int is the dimension of integral matrices
  !> \return val_ints contains the returned integrals
  subroutine gen1int_shell_int_tri_diag(prop_int, base_geo_derv, is_triang, &
                                        num_ao_orb, wrt_int, dim_int, val_ints)
    implicit none
    type(shell_int_t), intent(in) :: prop_int
    integer, intent(in) :: base_geo_derv
    logical, intent(in) :: is_triang
    integer, intent(in) :: num_ao_orb
    logical, intent(in) :: wrt_int
    integer, intent(in) :: dim_int
    real(REALK), intent(inout) :: val_ints(dim_int,*)
    integer addr_opt_derv            !address of operators and derivatives
    integer iderv, iopt, iket, ibra  !incremental recorders
    integer addr_ket                 !address of orbitals on ket center in returned integrals
    addr_opt_derv = base_geo_derv*prop_int%num_opt
    ! triangular format
    if (is_triang) then
      do iderv = 1, prop_int%num_geo_derv
        do iopt = 1, prop_int%num_opt
          addr_opt_derv = addr_opt_derv+1
          addr_ket = (prop_int%base_idx_ket-1)*prop_int%base_idx_ket/2
          do iket = 1, prop_int%num_orb_ket
            addr_ket = addr_ket+iket+prop_int%base_idx_ket-1  !=(iket+base_idx_ket-1)*(iket+base_idx_ket)/2
            do ibra = 1, iket  !only upper and diagonal parts returned
              val_ints(addr_ket+ibra+prop_int%base_idx_bra,addr_opt_derv) &
                = prop_int%contr_ints(ibra,iket,iopt,iderv)
            end do
          end do
        end do
      end do
  !FIXME
      ! the path, dimension, and kind of integral matrices should already be created
      if (wrt_int) then
      end if
    ! square format
    else
      do iderv = 1, prop_int%num_geo_derv
        do iopt = 1, prop_int%num_opt
          addr_opt_derv = addr_opt_derv+1
          do iket = 1, prop_int%num_orb_ket
            addr_ket = prop_int%base_idx_ket*num_ao_orb
            do ibra = 1, prop_int%num_orb_bra
              val_ints(ibra+prop_int%base_idx_bra+addr_ket,addr_opt_derv) &
                = prop_int%contr_ints(ibra,iket,iopt,iderv)
            end do
          end do
        end do
      end do
!FIXME
      ! the path, dimension, and kind of integral matrices should already be created
      if (wrt_int) then
      end if
    end if
    return
  end subroutine gen1int_shell_int_tri_diag

  !> \brief gets the contracted property integrals between two AO sub-shells back,
  !>        and writes the integrals on file is required, the integrals are belong
  !>        to the off-diagonal parts in triangular format
  !> \author Bin Gao
  !> \date 2011-10-07
  !> \param prop_int contains the contracted property integrals between two AO sub-shells
  !> \param base_geo_derv is the base address of current total geometric derivatives in the
  !>        list of all total geometric derivatives
  !> \param is_triang indicates if returning integral matrices are in triangular format,
  !>        otherwise square format
  !> \param sym_int indicates if the integrals matrices are symmetric, otherwise anti-symmetric
  !> \param num_ao_orb is the total number of atomic orbitals
  !> \param wrt_int indicates if writing integrals on file
  !> \param dim_int is the dimension of integral matrices
  !> \return val_ints contains the returned integrals
  subroutine gen1int_shell_int_tri_off(prop_int, base_geo_derv, is_triang,    &
                                       sym_int, num_ao_orb, wrt_int, dim_int, &
                                       val_ints)
    implicit none
    type(shell_int_t), intent(in) :: prop_int
    integer, intent(in) :: base_geo_derv
    logical, intent(in) :: is_triang
    logical, intent(in) :: sym_int
    integer, intent(in) :: num_ao_orb
    logical, intent(in) :: wrt_int
    integer, intent(in) :: dim_int
    real(REALK), intent(inout) :: val_ints(dim_int,*)
    integer addr_opt_derv            !address of operators and derivatives
    integer iderv, iopt, iket, ibra  !incremental recorders
    integer addr_bra                 !address of orbitals on bra center in returned integrals
    integer addr_ket                 !address of orbitals on ket center in returned integrals
    addr_opt_derv = base_geo_derv*prop_int%num_opt
    ! triangular format
    if (is_triang) then
      do iderv = 1, prop_int%num_geo_derv
        do iopt = 1, prop_int%num_opt
          addr_opt_derv = addr_opt_derv+1
          addr_ket = (prop_int%base_idx_ket-1)*prop_int%base_idx_ket/2
          do iket = 1, prop_int%num_orb_ket
            addr_ket = addr_ket+iket+prop_int%base_idx_ket-1  !=(iket+base_idx_ket-1)*(iket+base_idx_ket)/2
            do ibra = 1, prop_int%num_orb_bra
              val_ints(addr_ket+ibra+prop_int%base_idx_bra,addr_opt_derv) &
                = prop_int%contr_ints(ibra,iket,iopt,iderv)
            end do
          end do
        end do
      end do
!FIXME
      ! the path, dimension, and kind of integral matrices should already be created
      if (wrt_int) then
      end if
    ! square format
    else
      ! symmetric matrices
      if (sym_int) then
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do iket = 1, prop_int%num_orb_ket
              addr_ket = prop_int%base_idx_ket*num_ao_orb
              do ibra = 1, prop_int%num_orb_bra
                val_ints(ibra+prop_int%base_idx_bra+addr_ket,addr_opt_derv) &
                  = prop_int%contr_ints(ibra,iket,iopt,iderv)
              end do
            end do
            do ibra = 1, prop_int%num_orb_bra
              addr_bra = prop_int%base_idx_bra*num_ao_orb
              do iket = 1, prop_int%num_orb_ket
                val_ints(iket+prop_int%base_idx_ket+addr_bra,addr_opt_derv) &
                  = prop_int%contr_ints(ibra,iket,iopt,iderv)
              end do
            end do
          end do
        end do
!FIXME
        ! the path, dimension, and kind of integral matrices should already be created
        if (wrt_int) then
        end if
      ! anti-symmetric matrices
      else
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do iket = 1, prop_int%num_orb_ket
              addr_ket = prop_int%base_idx_ket*num_ao_orb
              do ibra = 1, prop_int%num_orb_bra
                val_ints(ibra+prop_int%base_idx_bra+addr_ket,addr_opt_derv) &
                  = prop_int%contr_ints(ibra,iket,iopt,iderv)
              end do
            end do
            do ibra = 1, prop_int%num_orb_bra
              addr_bra = prop_int%base_idx_bra*num_ao_orb
              do iket = 1, prop_int%num_orb_ket
                val_ints(iket+prop_int%base_idx_ket+addr_bra,addr_opt_derv) &
                  = -prop_int%contr_ints(ibra,iket,iopt,iderv)
              end do
            end do
          end do
        end do
!FIXME
        ! the path, dimension, and kind of integral matrices should already be created
        if (wrt_int) then
        end if
      end if
    end if
    return
  end subroutine gen1int_shell_int_tri_off

  !> \brief gets the contracted property integrals between two AO sub-shells back,
  !>        and writes the integrals on file is required, the integrals are parts
  !>        of square matrix
  !> \author Bin Gao
  !> \date 2011-10-07
  !> \param prop_int contains the contracted property integrals between two AO sub-shells
  !> \param base_geo_derv is the base address of current total geometric derivatives in the
  !>        list of all total geometric derivatives
  !> \param num_ao_orb is the total number of atomic orbitals
  !> \param wrt_int indicates if writing integrals on file
  !> \return val_ints contains the returned integrals
  subroutine gen1int_shell_int_square(prop_int, base_geo_derv, num_ao_orb, &
                                      wrt_int, val_ints)
    implicit none
    type(shell_int_t), intent(in) :: prop_int
    integer, intent(in) :: base_geo_derv
    integer, intent(in) :: num_ao_orb
    logical, intent(in) :: wrt_int
    real(REALK), intent(inout) :: val_ints(num_ao_orb,num_ao_orb,*)
    integer addr_opt_derv            !address of operators and derivatives
    integer iderv, iopt, iket, ibra  !incremental recorders
    integer addr_ket                 !address of orbitals on ket center in returned integrals
    addr_opt_derv = base_geo_derv*prop_int%num_opt
    do iderv = 1, prop_int%num_geo_derv
      do iopt = 1, prop_int%num_opt
        addr_opt_derv = addr_opt_derv+1
        do iket = 1, prop_int%num_orb_ket
          addr_ket = iket+prop_int%base_idx_ket
          do ibra = 1, prop_int%num_orb_bra
            val_ints(ibra+prop_int%base_idx_bra,addr_ket,addr_opt_derv) &
              = prop_int%contr_ints(ibra,iket,iopt,iderv)
          end do
        end do
      end do
    end do
!FIXME
    ! the path, dimension, and kind of integral matrices should already be created
    if (wrt_int) then
    end if
    return
  end subroutine gen1int_shell_int_square

  !> \brief calculates the expectation values and writes them on file is required,
  !>        the involved integrals are belong to the diagonal parts in triangular format
  !> \author Bin Gao
  !> \date 2011-10-11
  !> \param prop_int contains the contracted property integrals between two AO sub-shells
  !> \param base_geo_derv is the base address of current total geometric derivatives in the
  !>        list of all total geometric derivatives
  !> \param sym_int indicates if the integrals matrices are symmetric, otherwise anti-symmetric
  !> \param num_ao_orb is the total number of atomic orbitals
  !> \param kind_dens indicates if the kind of AO density matrices, 1 for symmetric, -1 for
  !>        anti-symmetric, others for square
  !> \param dim_dens is the dimension of AO density matrices
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \return val_expt contains the calculated expectation values
  subroutine gen1int_shell_tr_tri_diag(prop_int, base_geo_derv, sym_int, &
                                       num_ao_orb, kind_dens, dim_dens,  &
                                       num_dens, ao_dens, val_expt)
    implicit none
    type(shell_int_t), intent(in) :: prop_int
    integer, intent(in) :: base_geo_derv
    logical, intent(in) :: sym_int
    integer, intent(in) :: num_ao_orb
    integer, intent(in) :: kind_dens
    integer, intent(in) :: dim_dens
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: ao_dens(dim_dens,num_dens)
    real(REALK), intent(inout) :: val_expt(num_dens,*)
    ! zero expectation values
    if ((kind_dens==1 .and. .not.sym_int) .and. (kind_dens==-1 .and. sym_int)) then
      return
    else
      call gen1int_shell_tr_square(prop_int, base_geo_derv, num_ao_orb,    &
                                   kind_dens, dim_dens, num_dens, ao_dens, &
                                   val_expt)
    end if
    return
  end subroutine gen1int_shell_tr_tri_diag

  !> \brief calculates the expectation values and writes them on file is required,
  !>        the involved integrals are belong to the off-diagonal parts in triangular format
  !> \author Bin Gao
  !> \date 2011-10-11
  !> \param prop_int contains the contracted property integrals between two AO sub-shells
  !> \param base_geo_derv is the base address of current total geometric derivatives in the
  !>        list of all total geometric derivatives
  !> \param sym_int indicates if the integrals matrices are symmetric, otherwise anti-symmetric
  !> \param num_ao_orb is the total number of atomic orbitals
  !> \param kind_dens indicates if the kind of AO density matrices, 1 for symmetric, -1 for
  !>        anti-symmetric, others for square
  !> \param dim_dens is the dimension of AO density matrices
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \return val_expt contains the calculated expectation values
  subroutine gen1int_shell_tr_tri_off(prop_int, base_geo_derv, sym_int, &
                                      num_ao_orb, kind_dens, dim_dens,  &
                                      num_dens, ao_dens, val_expt)
    implicit none
    type(shell_int_t), intent(inout) :: prop_int
    integer, intent(in) :: base_geo_derv
    logical, intent(in) :: sym_int
    integer, intent(in) :: num_ao_orb
    integer, intent(in) :: kind_dens
    integer, intent(in) :: dim_dens
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: ao_dens(dim_dens,num_dens)
    real(REALK), intent(inout) :: val_expt(num_dens,*)
    real(REALK), allocatable :: tmp_expt(:,:,:)  !temporary expectation values
    integer addr_opt_derv                        !address of operators and derivatives
    integer iderv, iopt, idens                   !incremental recorders
    integer ierr                                 !error information
    ! zero expectation values
    if ((kind_dens==1 .and. .not.sym_int) .and. (kind_dens==-1 .and. sym_int)) then
      return
    else
      ! both AO density and integral matrices are symmetric or anti-symmetric
      if ((kind_dens==1 .and. sym_int) .or. (kind_dens==-1 .and. .not.sym_int)) then
        allocate(tmp_expt(num_dens,prop_int%num_opt,prop_int%num_geo_derv), &
                 stat=ierr)
        if (ierr/=0) stop "gen1int_shell_tr_tri_off>> failed to allocate tmp_expt!"
        call gen1int_shell_tr_square(prop_int, 0, num_ao_orb, kind_dens, &
                                     dim_dens, num_dens, ao_dens, tmp_expt)
        addr_opt_derv = base_geo_derv*prop_int%num_opt
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do idens = 1, num_dens
              val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                + 2.0_REALK*tmp_expt(idens,iopt,iderv)
            end do
          end do
        end do
        deallocate(tmp_expt)
      else
        call gen1int_shell_tr_square(prop_int, base_geo_derv, num_ao_orb,    &
                                     kind_dens, dim_dens, num_dens, ao_dens, &
                                     val_expt)
        ! takes the opposite sign for anti-symmetric integral matrices
        if (.not.sym_int) prop_int%contr_ints = -prop_int%contr_ints
        call gen1int_shell_tr_square(prop_int, base_geo_derv, num_ao_orb,    &
                                     kind_dens, dim_dens, num_dens, ao_dens, &
                                     val_expt)
      end if
    end if
    return
  end subroutine gen1int_shell_tr_tri_off

  !> \brief calculates the expectation values and writes them on file is required,
  !>        the involved integrals are parts of square matrix
  !> \author Bin Gao
  !> \date 2011-10-11
  !> \param prop_int contains the contracted property integrals between two AO sub-shells
  !> \param base_geo_derv is the base address of current total geometric derivatives in the
  !>        list of all total geometric derivatives
  !> \param num_ao_orb is the total number of atomic orbitals
  !> \param kind_dens indicates if the kind of AO density matrices, 1 for symmetric, -1 for
  !>        anti-symmetric, others for square
  !> \param dim_dens is the dimension of AO density matrices
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \return val_expt contains the calculated expectation values
  subroutine gen1int_shell_tr_square(prop_int, base_geo_derv, num_ao_orb,    &
                                     kind_dens, dim_dens, num_dens, ao_dens, &
                                     val_expt)
    implicit none
    type(shell_int_t), intent(in) :: prop_int
    integer, intent(in) :: base_geo_derv
    integer, intent(in) :: num_ao_orb
    integer, intent(in) :: kind_dens
    integer, intent(in) :: dim_dens
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: ao_dens(dim_dens,num_dens)
    real(REALK), intent(inout) :: val_expt(num_dens,*)
    integer addr_opt_derv                   !address of operators and derivatives
    integer iderv, iopt, idens, iket, ibra  !incremental recorders
    integer addr_ket                        !address of orbitals on ket center in AO density matrices
    integer addr_bra                        !address of orbitals on bra center in AO density matrices
    addr_opt_derv = base_geo_derv*prop_int%num_opt
    select case(kind_dens)
    ! symmetric matrices (column- or ket-major, upper and diagonal parts)
    case(1)
      ! diagonal part integral matrices
      if (prop_int%base_idx_bra==prop_int%base_idx_ket) then
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do idens = 1, num_dens
              ! lower and diagonal parts
              addr_bra = (prop_int%base_idx_bra-1)*prop_int%base_idx_bra/2
              do ibra = 1, prop_int%num_orb_bra
                addr_bra = addr_bra+ibra+prop_int%base_idx_bra-1  !=(ibra+base_idx_bra-1)*(ibra+base_idx_bra)/2
                do iket = 1, ibra
                  val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                    + prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                    * ao_dens(addr_bra+iket+prop_int%base_idx_ket,idens)
                end do
              end do
              ! upper part
              addr_ket = (prop_int%base_idx_ket-1)*prop_int%base_idx_ket/2
              do iket = 1, prop_int%num_orb_ket
                addr_ket = addr_ket+iket+prop_int%base_idx_ket-1  !=(iket+base_idx_ket-1)*(iket+base_idx_ket)/2
                do ibra = 1, iket-1
                  val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                    + prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                    * ao_dens(addr_ket+ibra+prop_int%base_idx_bra,idens)
                end do
              end do
            end do
          end do
        end do
      ! upper part integral matrices, needs lower part AO density matrices
      else if (prop_int%base_idx_bra<prop_int%base_idx_ket) then
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do idens = 1, num_dens
              addr_ket = (prop_int%base_idx_ket-1)*prop_int%base_idx_ket/2
              do iket = 1, prop_int%num_orb_ket
                addr_ket = addr_ket+iket+prop_int%base_idx_ket-1  !=(iket+base_idx_ket-1)*(iket+base_idx_ket)/2
                do ibra = 1, prop_int%num_orb_bra
                  val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                    + prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                    * ao_dens(addr_ket+ibra+prop_int%base_idx_bra,idens)
                end do
              end do
            end do
          end do
        end do
      ! lower part integral matrices, needs upper part AO density matrices,
      ! or the lower part AO density matrices
      else
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do idens = 1, num_dens
              addr_bra = (prop_int%base_idx_bra-1)*prop_int%base_idx_bra/2
              do ibra = 1, prop_int%num_orb_bra
                addr_bra = addr_bra+ibra+prop_int%base_idx_bra-1  !=(ibra+base_idx_bra-1)*(ibra+base_idx_bra)/2
                do iket = 1, prop_int%num_orb_ket
                  val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                    + prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                    * ao_dens(addr_bra+iket+prop_int%base_idx_ket,idens)
                end do
              end do
            end do
          end do
        end do
      end if
    ! anti-symmetric matrices (column- or ket-major, upper and diagonal parts)
    case(-1)
      ! diagonal part integral matrices
      if (prop_int%base_idx_bra==prop_int%base_idx_ket) then
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do idens = 1, num_dens
              ! lower and diagonal parts
              addr_bra = (prop_int%base_idx_bra-1)*prop_int%base_idx_bra/2
              do ibra = 1, prop_int%num_orb_bra
                addr_bra = addr_bra+ibra+prop_int%base_idx_bra-1  !=(ibra+base_idx_bra-1)*(ibra+base_idx_bra)/2
                do iket = 1, ibra
                  val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                    - prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                    * ao_dens(addr_bra+iket+prop_int%base_idx_ket,idens)
                end do
              end do
              ! upper part
              addr_ket = (prop_int%base_idx_ket-1)*prop_int%base_idx_ket/2
              do iket = 1, prop_int%num_orb_ket
                addr_ket = addr_ket+iket+prop_int%base_idx_ket-1  !=(iket+base_idx_ket-1)*(iket+base_idx_ket)/2
                do ibra = 1, iket-1
                  val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                    + prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                    * ao_dens(addr_ket+ibra+prop_int%base_idx_bra,idens)
                end do
              end do
            end do
          end do
        end do
      ! upper part integral matrices, needs lower part AO density matrices
      else if (prop_int%base_idx_bra<prop_int%base_idx_ket) then
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do idens = 1, num_dens
              addr_ket = (prop_int%base_idx_ket-1)*prop_int%base_idx_ket/2
              do iket = 1, prop_int%num_orb_ket
                addr_ket = addr_ket+iket+prop_int%base_idx_ket-1  !=(iket+base_idx_ket-1)*(iket+base_idx_ket)/2
                do ibra = 1, prop_int%num_orb_bra
                  val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                    + prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                    * ao_dens(addr_ket+ibra+prop_int%base_idx_bra,idens)
                end do
              end do
            end do
          end do
        end do
      ! lower part integral matrices, needs upper part AO density matrices,
      ! or the opposite of lower part AO density matrices
      else
        do iderv = 1, prop_int%num_geo_derv
          do iopt = 1, prop_int%num_opt
            addr_opt_derv = addr_opt_derv+1
            do idens = 1, num_dens
              addr_bra = (prop_int%base_idx_bra-1)*prop_int%base_idx_bra/2
              do ibra = 1, prop_int%num_orb_bra
                addr_bra = addr_bra+ibra+prop_int%base_idx_bra-1  !=(ibra+base_idx_bra-1)*(ibra+base_idx_bra)/2
                do iket = 1, prop_int%num_orb_ket
                  val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                    - prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                    * ao_dens(addr_bra+iket+prop_int%base_idx_ket,idens)
                end do
              end do
            end do
          end do
        end do
      end if
    ! square matrices (column- or ket-major)
    case default
      do iderv = 1, prop_int%num_geo_derv
        do iopt = 1, prop_int%num_opt
          addr_opt_derv = addr_opt_derv+1
          do idens = 1, num_dens
            do ibra = 1, prop_int%num_orb_bra
              addr_bra = prop_int%base_idx_bra*num_ao_orb
              do iket = 1, prop_int%num_orb_ket
                val_expt(idens,addr_opt_derv) = val_expt(idens,addr_opt_derv) &
                  + prop_int%contr_ints(ibra,iket,iopt,iderv)                 &
                  * ao_dens(iket+prop_int%base_idx_ket+addr_bra,idens)
              end do
            end do
          end do
        end do
      end do
    end select
    return
  end subroutine gen1int_shell_tr_square

  !> \brief calculates the Cartesian multipole moment integrals using contracted GTOs
  !> \author Bin Gao
  !> \date 2011-10-05
  !> \param bra_shell is the sub-shell on bra center
  !> \param ket_shell is the sub-shell on ket center
  !> \param is_lao indicates if using London atomic orbitals
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scaling constant for Cartesian multipole moments
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param num_cents is the number of geometric differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \param num_ao_bra is the number of AOs on bra center
  !> \param num_contr_bra is the number of contractions on bra center
  !> \param num_ao_ket is the number of AOs on ket center
  !> \param num_contr_ket is the number of contractions on ket center
  !> \param num_opt is the number of operators including derivatives
  !> \return contr_ints contains the calculated contracted integrals
  subroutine gen1int_shell_carmom(bra_shell, ket_shell, is_lao, order_elec,         &
                                  idx_diporg, dipole_origin, scal_const, order_mom, &
                                  order_mag_bra, order_mag_ket, order_mag_total,    &
                                  order_ram_bra, order_ram_ket, order_ram_total,    &
                                  order_geo_bra, order_geo_ket, order_geo_mom,      &
                                  num_cents, idx_cent, order_cent, num_ao_bra,      &
                                  num_contr_bra, num_ao_ket, num_contr_ket,         &
                                  num_opt, contr_ints)
    implicit none
    type(sub_shell_t), intent(in) :: bra_shell
    type(sub_shell_t), intent(in) :: ket_shell
    logical, intent(in) :: is_lao
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: num_ao_bra
    integer, intent(in) :: num_contr_bra
    integer, intent(in) :: num_ao_ket
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_ao_bra,num_contr_bra, &
                                           num_ao_ket,num_contr_ket,num_opt)
    ! contracted integrals from Gen1Int
    real(REALK), allocatable :: gen_ints(:,:,:,:,:)
    !error information
    integer ierr
#ifdef BUILD_GEN1INT
    ! London atomic orbitals
    if (is_lao) then
      stop "gen1int_shell_carmom>> LAO is not implemented!"
      ! SGTOs (Dalton does not use mixed Cartesian and spherical GTOs)
      if (bra_shell%is_sgto .or. ket_shell%is_sgto) then
      ! CGTOs
      else
      end if
    else
      ! SGTOs (Dalton does not use mixed Cartesian and spherical GTOs)
      if (bra_shell%is_sgto .or. ket_shell%is_sgto) then
        call contr_sgto_carmom(bra_shell%idx_cent, bra_shell%coord_cent, &
                               bra_shell%ang_num, bra_shell%num_prim,    &
                               bra_shell%exponents, bra_shell%num_contr, &
                               bra_shell%contr_coef,                     &
                               ket_shell%idx_cent, ket_shell%coord_cent, &
                               ket_shell%ang_num, ket_shell%num_prim,    &
                               ket_shell%exponents, ket_shell%num_contr, &
                               ket_shell%contr_coef,                     &
                               order_elec, idx_diporg, dipole_origin,    &
                               scal_const, order_mom, order_geo_bra,     &
                               order_geo_ket, order_geo_mom, num_cents,  &
                               idx_cent, order_cent, num_ao_bra,         &
                               num_ao_ket, num_opt, contr_ints)
        ! reorders p-shell
        if (bra_shell%ang_num==1)                                            &
          call dal_reorder_p_sgto(1, num_contr_bra*num_ao_ket*num_contr_ket, &
                                  num_opt, contr_ints)
        if (ket_shell%ang_num==1)                           &
          call dal_reorder_p_sgto(num_ao_bra*num_contr_bra, &
                                  num_contr_ket, num_opt, contr_ints)
      ! CGTOs
      else
        ! s-shell
        if (bra_shell%ang_num==0 .and. ket_shell%ang_num==0) then
          call contr_cgto_carmom(bra_shell%idx_cent, bra_shell%coord_cent, &
                                 bra_shell%ang_num, bra_shell%num_prim,    &
                                 bra_shell%exponents, bra_shell%num_contr, &
                                 bra_shell%contr_coef,                     &
                                 ket_shell%idx_cent, ket_shell%coord_cent, &
                                 ket_shell%ang_num, ket_shell%num_prim,    &
                                 ket_shell%exponents, ket_shell%num_contr, &
                                 ket_shell%contr_coef,                     &
                                 order_elec, idx_diporg, dipole_origin,    &
                                 scal_const, order_mom, order_geo_bra,     &
                                 order_geo_ket, order_geo_mom, num_cents,  &
                                 idx_cent, order_cent, num_ao_bra,         &
                                 num_ao_ket, num_opt, contr_ints)
        else
          ! contracted integrals from Gen1Int
          allocate(gen_ints(num_ao_bra,num_contr_bra, &
                            num_ao_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0) stop "gen1int_shell_carmom>> failed to allocate gen_ints!"
          call contr_cgto_carmom(bra_shell%idx_cent, bra_shell%coord_cent, &
                                 bra_shell%ang_num, bra_shell%num_prim,    &
                                 bra_shell%exponents, bra_shell%num_contr, &
                                 bra_shell%contr_coef,                     &
                                 ket_shell%idx_cent, ket_shell%coord_cent, &
                                 ket_shell%ang_num, ket_shell%num_prim,    &
                                 ket_shell%exponents, ket_shell%num_contr, &
                                 ket_shell%contr_coef,                     &
                                 order_elec, idx_diporg, dipole_origin,    &
                                 scal_const, order_mom, order_geo_bra,     &
                                 order_geo_ket, order_geo_mom, num_cents,  &
                                 idx_cent, order_cent, num_ao_bra,         &
                                 num_ao_ket, num_opt, gen_ints)
          ! reorder the integrals
          call reorder_cgto_ints(bra_shell%ang_num, num_ao_bra, bra_shell%powers, &
                                 ket_shell%ang_num, num_ao_ket, ket_shell%powers, &
                                 num_contr_bra, num_contr_ket, num_opt, gen_ints, &
                                 contr_ints)
          deallocate(gen_ints)
        end if
      end if
    end if
    return
#else
    stop "gen1int_shell_carmom>> Gen1Int is not installed!"
#endif
  end subroutine gen1int_shell_carmom

  !> \brief calculates the nuclear attraction potential (with Cartesian multipole
  !>        moment) integrals using contracted GTOs
  !> \author Bin Gao
  !> \date 2011-10-05
  !> \param bra_shell is the sub-shell on bra center
  !> \param ket_shell is the sub-shell on ket center
  !> \param is_lao indicates if using London atomic orbitals
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_nucorg is the atomic center of nuclear potential origin (<1 for non-atomic center)
  !> \param nucpot_origin contains the coordinates of nuclear potential origin
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scaling constant for Cartesian multipole moments
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_nuc is the order of geometric derivatives on nuclear attraction potential origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param num_cents is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the order of derivatives of the differentiated centers
  !> \param num_ao_bra is the number of AOs on bra center
  !> \param num_contr_bra is the number of contractions on bra center
  !> \param num_ao_ket is the number of AOs on ket center
  !> \param num_contr_ket is the number of contractions on ket center
  !> \param num_opt is the number of operators including derivatives
  !> \return contr_ints contains the calculated contracted integrals
  subroutine gen1int_shell_nucpot(bra_shell, ket_shell, is_lao, order_elec,       &
                                  idx_nucorg, nucpot_origin, idx_diporg,          &
                                  dipole_origin, scal_const, order_mom,           &
                                  order_mag_bra, order_mag_ket, order_mag_total,  &
                                  order_ram_bra, order_ram_ket, order_ram_total,  &
                                  order_geo_bra, order_geo_ket, order_geo_nuc,    &
                                  order_geo_mom, num_cents, idx_cent, order_cent, &
                                  num_ao_bra, num_contr_bra, num_ao_ket,          &
                                  num_contr_ket, num_opt, contr_ints)
    implicit none
    type(sub_shell_t), intent(in) :: bra_shell
    type(sub_shell_t), intent(in) :: ket_shell
    logical, intent(in) :: is_lao
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_nucorg
    real(REALK), intent(in) :: nucpot_origin(3)
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_geo_nuc
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: num_ao_bra
    integer, intent(in) :: num_contr_bra
    integer, intent(in) :: num_ao_ket
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_ao_bra,num_contr_bra, &
                                           num_ao_ket,num_contr_ket,num_opt)
    ! contracted integrals from Gen1Int
    real(REALK), allocatable :: gen_ints(:,:,:,:,:)
    !error information
    integer ierr
#ifdef BUILD_GEN1INT
    ! London atomic orbitals
    if (is_lao) then
      stop "gen1int_shell_nucpot>> LAO is not implemented!"
      ! SGTOs (Dalton does not use mixed Cartesian and spherical GTOs)
      if (bra_shell%is_sgto .or. ket_shell%is_sgto) then
      ! CGTOs
      else
      end if
    else
      ! SGTOs (Dalton does not use mixed Cartesian and spherical GTOs)
      if (bra_shell%is_sgto .or. ket_shell%is_sgto) then
        call contr_sgto_nucpot(bra_shell%idx_cent, bra_shell%coord_cent, &
                               bra_shell%ang_num, bra_shell%num_prim,    &
                               bra_shell%exponents, bra_shell%num_contr, &
                               bra_shell%contr_coef,                     &
                               ket_shell%idx_cent, ket_shell%coord_cent, &
                               ket_shell%ang_num, ket_shell%num_prim,    &
                               ket_shell%exponents, ket_shell%num_contr, &
                               ket_shell%contr_coef,                     &
                               order_elec, idx_nucorg, nucpot_origin,    &
                               idx_diporg, dipole_origin, scal_const,    &
                               order_mom, order_geo_bra, order_geo_ket,  &
                               order_geo_nuc, order_geo_mom, num_cents,  &
                               idx_cent, order_cent, num_ao_bra,         &
                               num_ao_ket, num_opt, contr_ints)
        ! reorders p-shell
        if (bra_shell%ang_num==1)                                            &
          call dal_reorder_p_sgto(1, num_contr_bra*num_ao_ket*num_contr_ket, &
                                     num_opt, contr_ints)
        if (ket_shell%ang_num==1)                           &
          call dal_reorder_p_sgto(num_ao_bra*num_contr_bra, &
                                  num_contr_ket, num_opt, contr_ints)
      ! CGTOs
      else
        ! s-shell
        if (bra_shell%ang_num==0 .and. ket_shell%ang_num==0) then
          call contr_cgto_nucpot(bra_shell%idx_cent, bra_shell%coord_cent, &
                                 bra_shell%ang_num, bra_shell%num_prim,    &
                                 bra_shell%exponents, bra_shell%num_contr, &
                                 bra_shell%contr_coef,                     &
                                 ket_shell%idx_cent, ket_shell%coord_cent, &
                                 ket_shell%ang_num, ket_shell%num_prim,    &
                                 ket_shell%exponents, ket_shell%num_contr, &
                                 ket_shell%contr_coef,                     &
                                 order_elec, idx_nucorg, nucpot_origin,    &
                                 idx_diporg, dipole_origin, scal_const,    &
                                 order_mom, order_geo_bra, order_geo_ket,  &
                                 order_geo_nuc, order_geo_mom, num_cents,  &
                                 idx_cent, order_cent, num_ao_bra,         &
                                 num_ao_ket, num_opt, contr_ints)
        else
          ! contracted integrals from Gen1Int
          allocate(gen_ints(num_ao_bra,num_contr_bra, &
                            num_ao_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0) stop "gen1int_shell_nucpot>> failed to allocate gen_ints!"
          call contr_cgto_nucpot(bra_shell%idx_cent, bra_shell%coord_cent, &
                                 bra_shell%ang_num, bra_shell%num_prim,    &
                                 bra_shell%exponents, bra_shell%num_contr, &
                                 bra_shell%contr_coef,                     &
                                 ket_shell%idx_cent, ket_shell%coord_cent, &
                                 ket_shell%ang_num, ket_shell%num_prim,    &
                                 ket_shell%exponents, ket_shell%num_contr, &
                                 ket_shell%contr_coef,                     &
                                 order_elec, idx_nucorg, nucpot_origin,    &
                                 idx_diporg, dipole_origin, scal_const,    &
                                 order_mom, order_geo_bra, order_geo_ket,  &
                                 order_geo_nuc, order_geo_mom, num_cents,  &
                                 idx_cent, order_cent, num_ao_bra,         &
                                 num_ao_ket, num_opt, gen_ints)
          ! reorder the integrals
          call reorder_cgto_ints(bra_shell%ang_num, num_ao_bra, bra_shell%powers, &
                                 ket_shell%ang_num, num_ao_ket, ket_shell%powers, &
                                 num_contr_bra, num_contr_ket, num_opt, gen_ints, &
                                 contr_ints)
          deallocate(gen_ints)
        end if
      end if
    end if
    return
#else
    stop "gen1int_shell_carmom>> Gen1Int is not installed!"
#endif
  end subroutine gen1int_shell_nucpot

  !> \brief reorders the p-shell contracted real solid-harmonic Gaussians in Dalton's
  !>        order on bra or ket center
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param dim_bra_sgto is the dimension of SGTOs on bra center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_opt is the number of operators
  !> \param gen_ints contains the contracted integrals from Gen1Int
  subroutine dal_reorder_p_sgto(dim_bra_sgto, num_contr_ket, num_opt, gen_ints)
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
    if (ierr/=0) stop "dal_reorder_p_sgto>> failed to allocate pshell_ints!"
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
  end subroutine 

end module gen1int_shell
