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

#ifdef BUILD_GEN1INT
C>    \brief sets up AO-blocks
C>    \author Bin Gao
C>    \date 2010-12-06
      subroutine gen1int_set_block
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "shells.h"
#include "gen1int.h"
C     incremental recorder over shells
      integer ibra
      call QENTER('gen1int_set_block')
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
      call QEXIT('gen1int_set_block')
1000  format("GEN1INT_ERROR>> ",A,I8)
#ifdef GEN1INT_DEBUG
9000  format("GEN1INT_DEBUG>> ",A,I8)
9001  format("GEN1INT_DEBUG>> IDX_SHELL:",5I8)
#endif
      return
      end subroutine gen1int_set_block
#endif

C>    \brief calculates integrals using Cartesian or real solid-harmonic Gaussian type orbitals
C>    \detail modified files include
C>            (1) abacus/dalgnr.F  -- adding interface of calling gen1int_prop(), such as
C>                (a) including "gen1int.h"
C>                (a) modifying RUN_GENINT
C>                (b) calling gen1int_prop()
C>                (c) adding '*GEN1IN' into TABDIR
C>            (2) abacus/abaopt.F -- adding interface of calling gen1int_prop()
C>                and including "gen1int.h"
C>    \author Bin Gao
C>    \date 2010-07-28
C>    \param WORK is the workspace of Dalton
C>    \param LWORK is the dimension of the workspace
      subroutine gen1int_prop(WORK, LWORK)
#ifdef BUILD_GEN1INT
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"
#include "lmns.h"
#include "gen1int.h"
C     workspace of Dalton
      integer LWORK
      dimension WORK(LWORK)
C     name of property integrals to calculated
      character*15 prop_name
C     start and end of the read in keyword
      integer start_prop, end_prop
C     start of comment
      integer start_pound, start_bang
C     start of primitives
      integer st_bra_prim, st_ket_prim
C     end of primitives
      integer end_bra_prim, end_ket_prim
C     start of contractions
      integer st_bra_contr, st_ket_contr
C     end of contractions
      integer end_bra_contr, end_ket_contr
C     incremental recorders over AO-blocks
      integer ibra, iket
      call QENTER('gen1int_prop')
      call GETTIM(TIMHER, WALHER)
      call FLSHFO(LUPRI)
      call TITLER('Starting Gen1Int', ' ', 200)
      if (RUN_GENINT) then
C       processes DALTON.INP
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
#ifdef GEN1INT_DEBUG
        write(LUPRI,9000) "Number of shells:", KMAX
#endif
C       sets up AO-blocks
        call gen1int_set_block
C       tries to find property integrals to calculate
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
C           doing calculations now ...
            call TITLER('Calculating '//prop_name(start_prop:end_prop),
     &                  '*', 118)
C           loops over blocks of bra
            do ibra = 1, NUM_BLOCK
C             generates the Cartesian powers
              call LMNVAL(NHKT(IDX_SHELL(ibra)), KCKT(IDX_SHELL(ibra)),
     &                    LVALUA, MVALUA, NVALUA)
              st_bra_prim = JSTRT(IDX_SHELL(ibra))+1
              end_bra_prim = JSTRT(IDX_SHELL(ibra))
     &                     + NUCO(IDX_SHELL(ibra))
              st_bra_contr = NUMCF(IDX_SHELL(ibra))
              end_bra_contr = NUMCF(IDX_SHELL(ibra))
     &                      + NRCO(IDX_SHELL(ibra))-1
#ifdef GEN1INT_DEBUG
              write(LUPRI,9000) "Center of bra:", NCENT(IDX_SHELL(ibra))
              write(LUPRI,9100) "Coordinates of bra",
     &                          CENT(IDX_SHELL(ibra),:,1)
              write(LUPRI,9000) "Angular number of bra:",
     &                          NHKT(IDX_SHELL(ibra))-1
              write(LUPRI,9000) "Number of bra prim.:",
     &                          NUCO(IDX_SHELL(ibra))
              write(LUPRI,9000) "Number of bra contr.:",
     &                          NRCO(IDX_SHELL(ibra))
              write(LUPRI,9000) "Number of bra Cart.:",
     &                          KCKT(IDX_SHELL(ibra))
              write(LUPRI,9000) "Start of bra prim.:", st_bra_prim
              write(LUPRI,9000) "End of bra prim.:", end_bra_prim
              write(LUPRI,9200) PRIEXP(st_bra_prim:end_bra_prim)
              write(LUPRI,9000) "Start of bra contr.:", st_bra_contr
              write(LUPRI,9000) "End of bra contr.:", end_bra_contr
              write(LUPRI,9300) PRICCF(st_bra_prim:end_bra_prim,
     &                                 st_bra_contr:end_bra_contr)
#endif
C             loops over shells of ket
              do iket = 1, NUM_BLOCK
C               generates the Cartesian powers
                call LMNVAL(NHKT(IDX_SHELL(iket)),
     &                      KCKT(IDX_SHELL(iket)),
     &                      LVALUB, MVALUB, NVALUB)
                st_ket_prim = JSTRT(IDX_SHELL(iket))+1
                end_ket_prim = JSTRT(IDX_SHELL(iket))
     &                       + NUCO(IDX_SHELL(iket))
                st_ket_contr = NUMCF(IDX_SHELL(iket))
                end_ket_contr = NUMCF(IDX_SHELL(iket))
     &                        + NRCO(IDX_SHELL(iket))-1
                select case(prop_name(start_prop:end_prop))
C (a) contr_cgto_carmom: Cartesian multipole moments
C (b) contr_cgto_delta:  Delta function
C (c) contr_cgto_nucpot: nuclear attraction potential
C (d) contr_cgto_isdpot: inverse square distance potential
C (e) contr_cgto_gaupot: Gaussian charge potential
C (f) contr_cgto_dso:    diamagnetic spin-orbit coupling
C (g) contr_cgto_ecp:    effective core potential
C (h) contr_cgto_mcp1:   model core potential (Version 1)
C     &    ' Center of mass  (bohr):', (CMXYZ(I),I=1,3),
C     &    ' Operator center (bohr):', (ORIGIN(I),I=1,3),
C     &    ' Gauge origin    (bohr):', (GAGORG(I),I=1,3),
C     &    ' Dipole origin   (bohr):', (DIPORG(I),I=1,3)
                case('.OVERLAP')
                  call contr_cgto_carmom(NCENT(IDX_SHELL(ibra)),
     &              CENT(IDX_SHELL(ibra),:,1),
     &              NHKT(IDX_SHELL(ibra))-1, NUCO(IDX_SHELL(ibra)),
     &              PRIEXP(st_bra_prim:end_bra_prim),
     &              NRCO(IDX_SHELL(ibra)),
     &              PRICCF(st_bra_prim:end_bra_prim,
     &                     st_bra_contr:end_bra_contr),
     &              NCENT(IDX_SHELL(iket)), CENT(IDX_SHELL(iket),:,1),
     &              NHKT(IDX_SHELL(iket))-1, NUCO(IDX_SHELL(iket)),
     &              PRIEXP(st_ket_prim:end_ket_prim),
     &              NRCO(IDX_SHELL(iket)),
     &              PRICCF(st_ket_prim:end_ket_prim,
     &                     st_ket_contr:end_ket_contr), DIPORG)
                end select
C                if (SPHR(IDX_SHELL(iket))) then
C                  call cgto_to_sgto()
C                endif
C                if (SPHR(IDX_SHELL(ibra))) then
C                  call cgto_to_sgto()
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
      call QEXIT('gen1int_prop')
      return
#ifdef GEN1INT_DEBUG
9000  format("GEN1INT_DEBUG>> ",A,I8)
9100  format("GEN1INT_DEBUG>> ",A,3F16.8)
9200  format("GEN1INT_DEBUG>> Exponents:",5Es16.8)
9300  format("GEN1INT_DEBUG>> Contr. Coef.:",5Es16.8)
#endif
#else
      implicit none
C     workspace of Dalton
      integer LWORK
      double precision WORK(LWORK)
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_prop
