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
  end type sub_shell_t

  ! if the AO sub-shells from Dalton initialized
  logical, public, save :: shell_init = .false.
  ! number of AO sub-shells from Dalton
  integer, public, save :: num_ao_shells = 0
  ! AO sub-shells from Dalton
  type(sub_shell_t), public, allocatable, save :: ao_shells(:)

  public :: gen1int_shell_set
  public :: gen1int_shell_dims
  public :: gen1int_shell_num_orb
  public :: gen1int_shell_dump
  public :: gen1int_shell_clean
  public :: gen1int_shell_carmom
  public :: gen1int_shell_nucpot

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
  !> \return ao_shell is the initialized AO sub-shell
  subroutine gen1int_shell_set(is_sgto, idx_cent, coord_cent, ang_num,    &
                              num_prim, exponents, num_contr, contr_coef, &
                              powers, ao_shell)
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

  !> \brief gets the number of orbtials for the given AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-06
  !> \param num_shell is the number of AO sub-shells
  !> \param ao_shells are the AO sub-shells
  !> \return num_orb is the number of orbitals
  subroutine gen1int_shell_num_orb(num_shell, ao_shells, num_orb)
    implicit none
    integer, intent(in) :: num_shell
    type(sub_shell_t), intent(in) :: ao_shells(*)
    integer, intent(out) :: num_orb
    integer ishell  !incremental recorder
    num_orb = ao_shells(1)%num_ao*ao_shells(1)%num_contr
    do ishell = 2, num_shell
      num_orb = num_orb+ao_shells(ishell)%num_ao*ao_shells(ishell)%num_contr
    end do
    return
  end subroutine gen1int_shell_num_orb

  !> \brief dumps information of several AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param num_shell is the number of AO sub-shells
  !> \param ao_shells are the AO sub-shells to dump
  !> \param io_unit is the IO unit to dump
  subroutine gen1int_shell_dump(num_shell, ao_shells, io_unit)
    implicit none
    integer, intent(in) :: num_shell
    type(sub_shell_t), intent(in) :: ao_shells(*)
    integer, intent(in) :: io_unit
    integer ishell, iprim  !incremental recorder
    write(io_unit,100) "number of AO sub-shells", num_shell
    do ishell = 1, num_shell
      write(io_unit,100) "sub-shell", ishell
      write(io_unit,100) "index of atomic center of the sub-shell", &
                         ao_shells(ishell)%idx_cent
      write(io_unit,110) "coordinates of the atomic center", &
                         ao_shells(ishell)%coord_cent
      write(io_unit,100) "angular number", ao_shells(ishell)%ang_num
      write(io_unit,100) "number of primitive GTOs", ao_shells(ishell)%num_prim
      write(io_unit,100) "number of contractions", ao_shells(ishell)%num_contr
      write(io_unit,100) "    exponents      coefficients"
      do iprim = 1, ao_shells(ishell)%num_prim
        write(io_unit,120)                    &
          ao_shells(ishell)%exponents(iprim), &
          ao_shells(ishell)%contr_coef(1:min(4,ao_shells(ishell)%num_contr),iprim)
        write(io_unit,130) &
          ao_shells(ishell)%contr_coef(5:ao_shells(ishell)%num_contr,iprim)
      end do
      if (ao_shells(ishell)%is_sgto) then
        write(io_unit,100) "SGTOs used"
        !-write(io_unit,140) ao_shells(ishell)%mag_num
      else
        write(io_unit,100) "CGTOs used"
        write(io_unit,150) ao_shells(ishell)%powers
      end if
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

  !> \brief calculates the Cartesian multipole moment integrals using contracted GTOs
  !> \author Bin Gao
  !> \date 2010-07-28
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
    ! London atomic orbitals
    if (is_lao) then
      stop "gen1int_shell_carmom>> LAO is not implemented!"
      ! SGTOs
      if (bra_shell%is_sgto .or. ket_shell%is_sgto) then
      ! CGTOs
      else
      end if
    else
      ! SGTOs
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
        if (bra_shell%ang_num==1)                                               &
          call reorder_dalton_pshell(1, num_contr_bra*num_ao_ket*num_contr_ket, &
                                     num_opt, contr_ints)
        if (ket_shell%ang_num==1)                              &
          call reorder_dalton_pshell(num_ao_bra*num_contr_bra, &
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
  end subroutine gen1int_shell_carmom

  !> \brief calculates the nuclear attraction potential (with Cartesian multipole
  !>        moment) integrals using contracted GTOs
  !> \author Bin Gao
  !> \date 2010-07-28
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
    ! London atomic orbitals
    if (is_lao) then
      stop "gen1int_shell_nucpot>> LAO is not implemented!"
      ! SGTOs
      if (bra_shell%is_sgto .or. ket_shell%is_sgto) then
      ! CGTOs
      else
      end if
    else
      ! SGTOs
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
        if (bra_shell%ang_num==1)                                               &
          call reorder_dalton_pshell(1, num_contr_bra*num_ao_ket*num_contr_ket, &
                                     num_opt, contr_ints)
        if (ket_shell%ang_num==1)                              &
          call reorder_dalton_pshell(num_ao_bra*num_contr_bra, &
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
  end subroutine gen1int_shell_nucpot

  !> \brief reorders the p-shell contracted real solid-harmonic Gaussians in Dalton's
  !>        order on bra or ket center
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param dim_bra_sgto is the dimension of SGTOs on bra center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_opt is the number of operators
  !> \param gen_ints contains the contracted integrals from Gen1Int
  subroutine reorder_dalton_pshell(dim_bra_sgto, num_contr_ket, num_opt, gen_ints)
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
    if (ierr/=0) stop "reorder_dalton_pshell>> failed to allocate pshell_ints!"
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
