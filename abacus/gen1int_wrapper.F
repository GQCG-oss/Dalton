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
     &                           doint, len_work, work)
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
C     if doing calculations
      logical doint(2,2)
C     workspace of Dalton
      integer len_work
      dimension work(len_work)
#ifdef BUILD_GEN1INT
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "inforb.h"
#include "onecom.h"
#include "aovec.h"
#include "symmet.h"
#include "symind.h"
#include "primit.h"
#include "lmns.h"
#include "shells.h"
#include "gen1int.h"
#include "ecpinf.h"
#include "orgcom.h"
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

      IDENA = 0
      call IZERO(ISOFRA, 8)
C     loops over blocks of bra
      do ibra = 1, NUM_BLOCK
        ISHELA = IDX_SHELL(ibra)
        NHKTA = NHKT(ISHELA)        !angular quantum number (s=1,p=2,d=3 etc.)
        KHKTA = KHKT(ISHELA)        !number of spherical (Cartesian) components
        KCKTA = KCKT(ISHELA)        !number of Cartesian components
        ICA   = LCLASS(ISHELA)      !class: large component (1), small or Huckel(2), density fitting (0)
        SPHRA = SPHR(ISHELA)        !if spherical Gaussians
        call LMNVAL(NHKTA, KCKTA, LVALUA, MVALUA, NVALUA)
        NCENTA = NCENT(ISHELA)      !index of symmetry independent center
        ICENTA = NUCNUM(NCENTA,1)   !total number of nuclei
        MULA   = ISTBAO(ISHELA)     !stabiliser: basic sym. op. that do not move center
        MULTA  = MULT(MULA)         !multiplicity of center
        NUCA   = NUCO(ISHELA)       !number of uncontracted functions
        NUMCFA = NUMCF(ISHELA)      !index of shell in AO-block
        NRCA   = NRCO(ISHELA)       !number of contracted functions
C FIXME: JSTA, NBCH, NUMCF
        JSTA   = JSTRT(ISHELA)      !
        CORAX  = CENT(ISHELA,1,1)   !coordinates of center
        CORAY  = CENT(ISHELA,2,1)
        CORAZ  = CENT(ISHELA,3,1)
        IDENB0 = 0

        st_bra_prim = JSTA+1
        end_bra_prim = JSTA+NUCA
        st_bra_contr = NUMCFA
        end_bra_contr = NUMCFA+NRCA-1
C       computes symmetry integral pointers for contributions from this block
C       note that at present this assumes all components from a shell are included
        do I = 1, 8
          ISOFRB(I) = 0
C         loops over maximum number of Cartesian powers
          do ixyz = 1, MXAQN
            INDFA(I,ixyz) = -10 000 000
          end do
        end do
C FIXME: INDFA used in SYM1S, SYM1N, SYM1NC
        do NA = 1, KHKTA
          do irep = 0, MAXREP
            if (IBTAND(MULA,IBTXOR(irep,ISYMAO(NHKTA,NA))).eq.0) then
              ISOFRA(irep+1) = ISOFRA(irep+1)+1
              INDFA(irep+1,NA) = ISOFRA(irep+1)
            end if
          end do
        end do
#ifdef GEN1INT_DEBUG

        write(LUPRI,'(A,I4)') 'IA address offsets for shell ', ISHELA
        do NA = 1,KHKTA
          write(LUPRI,'(8(1X,I5))') (INDFA(I,NA), I = 1, MAXREP+1)
        end do

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

          ISHELB = IDX_SHELL(iket)
          LDIAG = ISHELA.eq.ISHELB        !
          NHKTB = NHKT(ISHELB)
          KHKTB = KHKT(ISHELB)
          KCKTB = KCKT(ISHELB)
          ICB   = LCLASS(ISHELB)
          SPHRB = SPHR(ISHELB)
          CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
          NCENTB = NCENT(ISHELB)
          NHKTAB = NHKTA + NHKTB
          MULB   = ISTBAO(ISHELB)
          MULTB  = MULT(MULB)
          NUCB   = NUCO(ISHELB)
          NUMCFB = NUMCF(ISHELB)
          NRCB   = NRCO(ISHELB)
          JSTB   = JSTRT(ISHELB)
          CORBX0 = CENT(ISHELB,1,1)
          CORBY0 = CENT(ISHELB,2,1)
          CORBZ0 = CENT(ISHELB,3,1)
          KHKTAB = KHKTA*KHKTB
          KCKTAB = KCKTA*KCKTB
          MAB    = IBTOR(MULA,MULB)       !
          KAB    = IBTAND(MULA,MULB)      !
          HKAB   = FMULT(KAB)             !symmetry factor, which determines how many different integrals
C
          SPHRAB = SPHRA .OR. SPHRB
C         generates the Cartesian powers
          call LMNVAL(NHKTB, KCKTB, LVALUB, MVALUB, NVALUB)

          st_ket_prim = JSTB+1
          end_ket_prim = JSTB+NUCB
          st_ket_contr = NUMCFB
          end_ket_contr = NUMCFB+NRCB-1

C         computes symmetry integral pointers for contributions from this block
C         note that at present this assumes all components from a shell are included
          do I = 1, 8
            do ixyz = 1, MXAQN
              INDFB(I,ixyz) = -10 000 000
            end do
          end do
          do NB = 1, KHKTB
            do irep = 0, MAXREP
              if (IBTAND(MULB,IBTXOR(irep,ISYMAO(NHKTB,NB))).eq.0) then
                ISOFRB(irep+1) = ISOFRB(irep+1)+1
                INDFB(irep+1,NB) = ISOFRB(irep+1)
              end if
            end do
          end do

          if (.not.doint(ICA,ICB)) cycle

C FIXME: extra calculations, ECP with symmetry
C         loops over symmetry operations

          IDENB = IDENB0 - KHKTB
          do isymop = 0, MAXOPR
            if (IBTAND(isymop,MAB).eq.0) then

              IDENB  = IDENB + KHKTB
              ICENTB = NUCNUM(NCENTB,ISYMOP+1)
              ONECEN = ICENTA .eq. ICENTB
              SIGNBX = PT(IBTAND(ISYMAX(1,1),ISYMOP))
              SIGNBY = PT(IBTAND(ISYMAX(2,1),ISYMOP))
              SIGNBZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))
              CORBX  = SIGNBX*CORBX0
              CORBY  = SIGNBY*CORBY0
              CORBZ  = SIGNBZ*CORBZ0

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
          IDENB0 = IDENB0+KHKTB*MULTB
        end do
        IDENA = IDENA+KHKTA*MULTA
      end do
C     cleans
      call IZERO(ISOFRA, 8)
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
