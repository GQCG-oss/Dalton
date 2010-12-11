C  gen1int: compute one-electron integrals using rotational London atomic-orbitals
C  Copyright 2009, 2010 Bin Gao, and Andreas Thorvaldsen
C
C  gen1int is free software: you can redistribute it and/or modify
C  it under the terms of the GNU Lesser General Public License as published by
C  the Free Software Foundation, either version 3 of the License, or
C  (at your option) any later version.
C
C  gen1int is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C  GNU Lesser General Public License for more details.
C
C  You should have received a copy of the GNU Lesser General Public License
C  along with gen1int. If not, see <http://www.gnu.org/licenses/>.
C
C  This file takes care of ECP using Gen1Int from Dalton.
C
C  2010-12-06, Bin Gao
C  * first version

C>    \brief reads ECP information
C>    \detail modified files include
C>            (1) abacus/herrdn.F  -- adding interface of calling gen1int_ecpinp()
C>    \author Bin Gao
C>    \date 2010-12-06
C>    \param idx_atom is the index of ECP atom in all atom list
C>    \param anum_atom is the atomic number of this ECP atom
C>    \param ecp_name is the ECP name
      subroutine gen1int_ecpinp(idx_atom, anum_atom, ecp_name)
      implicit none
      integer idx_atom
      integer anum_atom
      character*(*) ecp_name
#ifdef BUILD_GEN1INT
#include "gen1int.h"
      return
#else
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_ecpinp
