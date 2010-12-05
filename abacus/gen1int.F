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

C>    \brief calculates integrals using Cartesian or real solid-harmonic Gaussian type orbitals
C>    \detail modified files include
C>            (1) include/gnrinf.h -- adding GENINT
C>            (2) abacus/dalgnr.F  -- adding interface of calling gen1int(), such as
C>                (2.1) modifying GENINT
C>                (2.2) calling gen1int()
C>                (2.3) adding '*GEN1IN' into TABDIR
C>    \author Bin Gao
C>    \date 2010-07-28
C>    \param WORK is the workspace of Dalton
C>    \param LWORK is the dimension of the workspace
      subroutine gen1int(WORK, LWORK)
#ifdef BUILD_GEN1INT
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C#include "mxcent.h"
#include "maxorb.h"
Cinclude "maxaqn.h"
#include "gnrinf.h"
#include "shells.h"
C     workspace of Dalton
      integer LWORK
      dimension WORK(LWORK)
C     name of property integrals to calculated
      character*15 prop_name
C     start and end of the read in keyword
      integer start_prop, end_prop
C     start of comment
      integer start_pound, start_bang
C     incremental recorder over shells
      integer ibra, iket
      call QENTER('gen1int')
      call GETTIM(TIMHER, WALHER)
      call FLSHFO(LUPRI)
      call TITLER('Starting Gen1Int', ' ', 200)
      if (GENINT) then
C       (1) processing DALTON.INP
        call TITLER('Processing DALTON.INP', '*', 118)
        if (LUCMD.lt.0) then
          call GPOPEN(LUCMD, 'DALTON.INP', 'OLD', ' ', 'FORMATTED',
     &                IDUMMY, .FALSE.)
        end if
        rewind(LUCMD)
C       tries to find keywords of Gen1Int
        do while (.true.)
          read(LUCMD, '(A)', end=100) prop_name
          call UPCASE(prop_name)
          if (index(trim(prop_name), '*GEN1INT').gt.0) exit
        end do
        do while (.true.)
          read(LUCMD, '(A)', end=100) prop_name
          call UPCASE(prop_name)
          if (index(trim(prop_name), '*END').gt.0) exit
          start_prop = index(prop_name, '.')
          end_prop = len_trim(prop_name)
          start_pound = index(prop_name, '#')
          start_bang = index(prop_name, '!')
          if (start_pound.gt.0) end_prop = min(end_prop, start_pound)
          if (start_bang.gt.0) end_prop = min(end_prop, start_bang)
C         reads in a property integral to calculate
          if (start_prop.gt.0 .and. start_prop.lt.end_prop) then
C           (2) doing calculations now ...
            call TITLER('Calculating '//prop_name(start_prop:end_prop),
     &                  '*', 118)
C           loops over shells of bra
            do ibra = 1, KMAX
C             loops over shells of ket
              do iket = 1, KMAX
C                print *, SPHR(ibra), SPHR(iket)
                select case(prop_name(start_prop:end_prop))
                end select
C                if (DOCART) then
C                  call contr_cgto
C                else
C                  call contr_sgto
C                end if
              end do
            end do
          end if
        end do
      end if
100   call GPCLOSE(LUCMD, 'KEEP')
      call GETTIM(TEND, WEND)
      TIMHER = TEND-TIMHER
      WALHER = WEND-WALHER
      call TIMTXT('>>>> Total CPU  time used in Gen1Int:',
     &            TIMHER, LUPRI)
      call TIMTXT('>>>> Total wall time used in Gen1Int:',
     &            WALHER, LUPRI)

      call TITLER('End of Gen1Int', ' ', 200)
      call FLSHFO(LUPRI)
      call QEXIT('gen1int')
      return
#else
      implicit none
C     workspace of Dalton
      integer LWORK
      dimension WORK(LWORK)
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int
