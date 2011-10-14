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
!...  This file returns the number of operators (including different derivatives)
!...  and total geometric derivatives, kind of integral matrices according to given
!...  property name and the information of total geometric derivatives.
!
!...  2011-10-08, Bin Gao
!...  * first version

#include "stdout.h"

  !> \brief returns the number of operators (including different derivatives)
  !>        and total geometric derivatives, kind of integral matrices according
  !>        to given property name and the information of total geometric derivatives
  !> \note this file should be consistent with subroutine \fn(gen1int_shell_prop)
  !>       in gen1int_shell.F90!!
  !> \author Bin Gao
  !> \date 2011-10-08
  !> \param prop_name is the name of property integrals to calculated
  !> \param is_lao indicates if using London atomic orbitals
  !> \param order_mom is the order of Cartesian multipole moments, only used for
  !>        integrals "CARMOM"
  !> \param max_num_cent is the maximum number of geometric differentiated centers
  !> \param order_geo_total is the order of total geometric derivatives
  !> \param num_atoms is the number of atoms
  !> \return num_opt is the number of operators (including different derivatives)
  !>         and total geometric derivatives
  !> \return kind_int indicates if the kind of integral matrices, 1 for symmetric, -1 for
  !>         anti-symmetric, others for square
  subroutine gen1int_prop_attr(prop_name, is_lao, order_mom,             &
                               max_num_cent, order_geo_total, num_atoms, &
                               num_opt, kind_int)
    implicit none
    character*(*), intent(in) :: prop_name
    logical, intent(in) :: is_lao
    integer, intent(in) :: order_mom
    integer, intent(in) :: max_num_cent
    integer, intent(in) :: order_geo_total
    integer, intent(in) :: num_atoms
    integer, intent(out) :: num_opt
    integer, intent(out) :: kind_int
    ! number of total geometric derivatives
    integer num_geo_derv
#ifdef BUILD_GEN1INT
    ! gets the number of total geometric derivatives (up to four-center total geometric derivatives)
    if (max_num_cent>0 .and. order_geo_total>0) then
      call geom_total_num_derv(order_geo_total, max_num_cent, &
                               num_atoms, num_geo_derv)
    else
      num_geo_derv = 1
    end if
    ! different property integrals, please use alphabetical order!
    ! the following keywords and \var(num_opt) should be consistent
    ! with subroutine \fn(gen1int_shell_prop) in gen1int_shell.F90!!
    select case(trim(prop_name))
    ! electronic angular momentum around the nuclei
    case("ANGLON")
      num_opt = 3*num_geo_derv
      kind_int = 0
    ! Cartesian multipole integrals
    case("CARMOM")
      num_opt = (order_mom+1)*(order_mom+2)*num_geo_derv/2
      kind_int = 1
    ! dipole length integrals
    case("DIPLEN")
      num_opt = 3*num_geo_derv
      kind_int = 1
    ! kinetic energy integrals
    case("KINENERG")
      num_opt = num_geo_derv
      kind_int = 1
    ! London orbital contribution to nuclear shielding tensor integrals
    case("NSTLON")
      num_opt = 9*num_atoms
      kind_int = 0
    ! nuclear shielding integrals without London orbital contribution
    case("NSTNOL")
      num_opt = 9*num_atoms
      kind_int = 0
    ! overlap integrals
    case("OVERLAP")
      num_opt = num_geo_derv
      kind_int = 1
    ! one-electron potential energy integrals
    case("POTENERG")
      num_opt = num_geo_derv
      kind_int = 1
    ! paramagnetic spin-orbit integrals
    case("PSO")
      num_opt = 3*num_atoms
      kind_int = -1
    case default
      write(STDOUT,999) trim(prop_name)
      stop
    end select
    return
999 format("gen1int_prop_attr>> ",A," is not implemented!")
#else
    stop "gen1int_prop_attr>> Gen1Int is not installed!"
#endif
  end subroutine gen1int_prop_attr
