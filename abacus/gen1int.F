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
C  This file provides the interface of calling Gen1Int from Dalton.
C
C  2010-07-28, Bin Gao
C  * first version

C> \brief calculates integrals using Cartesian Gaussian type orbitals
C> \author Bin Gao
C> \date 2010-07-28
C> \param 
C> \return 
      subroutine gen1int_cgto()
#include<implicit.h>
#ifdef BUILD_GEN1INT
#else
#endif
1000  format(A)
      end subroutine gen1int_cgto

C> \brief calculates integrals using real solid-harmonic Gaussian type orbitals
C> \author Bin Gao
C> \date 2010-07-28
C> \param 
C> \return 
      subroutine gen1int_sgto()
#include<implicit.h>
#ifdef BUILD_GEN1INT
#else
#endif
1000  format(A)
      end subroutine gen1int_sgto
