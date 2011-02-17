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
C  This file sets AO blocks used for Gen1Int.
C
C  2010-12-06, Bin Gao
C  * first version

C>    \brief sets up AO-blocks
C>    \author Bin Gao
C>    \date 2010-12-06
      subroutine gen1int_block
#ifdef BUILD_GEN1INT
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "shells.h"
#include "gen1int.h"
C     incremental recorder over shells
      integer ibra
      call QENTER('gen1int_block')
C     initializes the number of AO-blocks
      NUM_BLOCK = 1
      ibra = 1
      IDX_SHELL(NUM_BLOCK) = ibra
      do ibra = 1, KMAX-1
C       finds a new AO-block
        if (JSTRT(ibra).ne.JSTRT(ibra+1)) then
          NUM_BLOCK = NUM_BLOCK+1
          if (NUM_BLOCK.gt.MAX_NBLOCK) then
            write(LUPRI,1000)
     &        "Maximum number of AO-blocks (MAX_NBLOCK)", MAX_NBLOCK
            call QUIT("Too many AO-blocks")
          else
            IDX_SHELL(NUM_BLOCK) = ibra+1
          end if
        end if
      end do
#ifdef GEN1INT_DEBUG
      write(LUPRI,9000) "Number of AO-blocks:", NUM_BLOCK
      write(LUPRI,9001) IDX_SHELL(1:NUM_BLOCK)
#endif
      call QEXIT('gen1int_block')
      return
1000  format("GEN1INT_ERROR>> ",A,I8)
#ifdef GEN1INT_DEBUG
9000  format("GEN1INT_DEBUG>> ",A,I8)
9001  format("GEN1INT_DEBUG>> IDX_SHELL:",5I8)
#endif
#else
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_block
