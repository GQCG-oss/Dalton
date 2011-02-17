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
C  2011-02-13, Bin Gao
C  * the subroutine does not read information from input anymore, it will only
C    be called by other subroutines
C
C  2010-07-28, Bin Gao
C  * first version

C>    \brief calculates integrals/expectation values using Cartesian or
C>           real solid-harmonic Gaussians
C>    \detail modified files include
C>            (1) abacus/herrdn.F: some modifications for processing different ECPs,
C>                                 adds the interface of calling gen1int_block.
C>            (2) abacus/her1car.F: removes the interface of calling old ECP subroutine.
C>            (3) abacus/her1drv.F: adds the interface of calling gen1int_wrapper().
C>            (4) include/gen1int.h: saves AO block information for Gen1Int
C>    \author Bin Gao
C>    \date 2010-07-28
      subroutine gen1int_wrapper(prop_name, order_geo, order_mag,
     &                           vals_int, propty, vals_expect,
     &                           len_work, work)
#include "implicit.h"
C     name of property integrals to calculated
      character*(*) prop_name
C     order of geometric derivatives
      integer order_geo
C     order of magnetic derivaives
      integer order_mag
C     integrals
      dimension vals_int(*)
C     if calculating expectation values
      logical propty
C     expectation values
      dimension vals_expect(*)
C     workspace of Dalton
      integer len_work
      dimension work(len_work)
#ifdef BUILD_GEN1INT
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "inforb.h"
#include "onecom.h"
#include "aovec.h"
#include "symmet.h"
#include "primit.h"
#include "lmns.h"
#include "shells.h"
#include "gen1int.h"
#include "ecpinf.h"
#include "ibtfun.h"
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
C     Dalton stuff
      call QENTER('gen1int_wrapper')
      call GETTIM(TIMHER, WALHER)
      call FLSHFO(LUPRI)
      call TITLER('Gen1Int calculates '//trim(prop_name), '*', 120)
C     sets up the order of geometric derivatives
C      IA0000 =  1
C      IA0X00 =  2
C      IA0Y00 =  3
C      IA0Z00 =  4
C      IAXX00 =  5
C      IAXY00 =  6
C      IAXZ00 =  7
C      IAYY00 =  8
C      IAYZ00 =  9

C     loops over blocks of bra
      do ibra = 1, NUM_BLOCK

C       generates the Cartesian powers
        call LMNVAL(NHKT(IDX_SHELL(ibra)), KCKT(IDX_SHELL(ibra)),
     &              LVALUA, MVALUA, NVALUA)
        st_bra_prim = JSTRT(IDX_SHELL(ibra))+1
        end_bra_prim = JSTRT(IDX_SHELL(ibra))
     &               + NUCO(IDX_SHELL(ibra))
        st_bra_contr = NUMCF(IDX_SHELL(ibra))
        end_bra_contr = NUMCF(IDX_SHELL(ibra))
     &                + NRCO(IDX_SHELL(ibra))-1
C       computes symmetry integral pointers for contributions from this block
C       note that at present this assumes all components from a shell are included

#ifdef GEN1INT_DEBUG
        write(LUPRI,9000) "Center of bra:", NCENT(IDX_SHELL(ibra))
        write(LUPRI,9100) "Coordinates of bra",
     &                    CENT(IDX_SHELL(ibra),:,1)
        write(LUPRI,9000) "Angular number of bra:",
     &                    NHKT(IDX_SHELL(ibra))-1
        write(LUPRI,9000) "Number of bra prim.:",
     &                    NUCO(IDX_SHELL(ibra))
        write(LUPRI,9000) "Number of bra contr.:",
     &                    NRCO(IDX_SHELL(ibra))
        write(LUPRI,9000) "Number of bra Cart.:",
     &                    KCKT(IDX_SHELL(ibra))
        write(LUPRI,9000) "Start of bra prim.:", st_bra_prim
        write(LUPRI,9000) "End of bra prim.:", end_bra_prim
        write(LUPRI,9200) PRIEXP(st_bra_prim:end_bra_prim)
        write(LUPRI,9000) "Start of bra contr.:", st_bra_contr
        write(LUPRI,9000) "End of bra contr.:", end_bra_contr
        write(LUPRI,9300) PRICCF(st_bra_prim:end_bra_prim,
     &                           st_bra_contr:end_bra_contr)
#endif
C FIXME: lower triangular?
C       loops over shells of ket
        do iket = 1, NUM_BLOCK

C FIXME: extra calculations, ECP with symmetry
C         loops over symmetry operations

          do isymop = 0, MAXOPR
            if (IBTAND(isymop,MAB).eq.0) then

C             generates the Cartesian powers
              call LMNVAL(NHKT(IDX_SHELL(iket)),
     &                    KCKT(IDX_SHELL(iket)),
     &                    LVALUB, MVALUB, NVALUB)
              st_ket_prim = JSTRT(IDX_SHELL(iket))+1
              end_ket_prim = JSTRT(IDX_SHELL(iket))
     &                     + NUCO(IDX_SHELL(iket))
              st_ket_contr = NUMCF(IDX_SHELL(iket))
              end_ket_contr = NUMCF(IDX_SHELL(iket))
     &                      + NRCO(IDX_SHELL(iket))-1
C             computes symmetry integral pointers for contributions from this block
C             note that at present this assumes all components from a shell are included

              select case(trim(prop_name))
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
     &            CENT(IDX_SHELL(ibra),:,1),
     &            NHKT(IDX_SHELL(ibra))-1, NUCO(IDX_SHELL(ibra)),
     &            PRIEXP(st_bra_prim:end_bra_prim),
     &            NRCO(IDX_SHELL(ibra)),
     &            PRICCF(st_bra_prim:end_bra_prim,
     &                   st_bra_contr:end_bra_contr),
     &            NCENT(IDX_SHELL(iket)), CENT(IDX_SHELL(iket),:,1),
     &            NHKT(IDX_SHELL(iket))-1, NUCO(IDX_SHELL(iket)),
     &            PRIEXP(st_ket_prim:end_ket_prim),
     &            NRCO(IDX_SHELL(iket)),
     &            PRICCF(st_ket_prim:end_ket_prim,
     &                   st_ket_contr:end_ket_contr), DIPORG)
              end select
C             transforms to spherical harmonic basis
C              if (SPHR(IDX_SHELL(iket))) then
C                call cgto_to_sgto()
C              endif
C              if (SPHR(IDX_SHELL(ibra))) then
C                call cgto_to_sgto()
C              end if
C              or
C              call SPHRM1()
C             computes expectation values
C              if (propty) then
C                collects density and Fock elements
C              end if
C FIXME: if undifferentiated, first- and second-order derivatives are calculated seperately
C             transforms integrals to SO basis
C              FULMAT = .TRUE.
C              ANTI   = .FALSE.
C              call SYM1S() or SYM1N() by MAXREP(irreps of the differentiation operator)
C FIXME:
C             writes differentiated integrals on file?
C              if (propty .and. (SECDER .OR. DIFDIP .OR. DIFQDP)) then
C              end if
            end if
          end do
        end do
      end do
C     Dalton stuff
      call GETTIM(TEND, WEND)
      TIMHER = TEND-TIMHER
      WALHER = WEND-WALHER
      write(LUPRI,'()')
      call TIMTXT('>>>> Total CPU  time used in Gen1Int:',
     &            TIMHER, LUPRI)
      call TIMTXT('>>>> Total wall time used in Gen1Int:',
     &            WALHER, LUPRI)
      call TITLER('End of Gen1Int', '*', 120)
      call FLSHFO(LUPRI)
      call QEXIT('gen1int_wrapper')
      return
#ifdef GEN1INT_DEBUG
9000  format("GEN1INT_DEBUG>> ",A,I8)
9100  format("GEN1INT_DEBUG>> ",A,3F16.8)
9200  format("GEN1INT_DEBUG>> Exponents:",5Es16.8)
9300  format("GEN1INT_DEBUG>> Contr. Coef.:",5Es16.8)
#endif
#else
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_wrapper
