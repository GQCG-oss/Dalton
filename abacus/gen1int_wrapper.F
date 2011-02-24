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
     &                           getint, vals_int, propty, vals_expect,
     &                           doint, len_work, work)
#include "implicit.h"
C     name of property integrals to calculated
      character*(*) prop_name
C     order of geometric derivatives
      integer order_geo
C     order of magnetic derivaives
      integer order_mag
C     if getting integrals back
      logical getint
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
C     incremental recorders over spherical (Cartesian) components
      integer icomp, jcomp
C     incremental recorders over contractions
      integer icontr, jcontr
C     incremental recorders over AO-shells
      integer ishell, jshell
C     incremental recorders over AO-blocks
      integer iblock, jblock
C     Dalton stuff
      call QENTER("gen1int_wrapper")
      call GETTIM(TIMHER, WALHER)
      call HEADER("Gen1Int calculates "//trim(prop_name),-1)
C     sets up the order of geometric derivatives
C-      IA0000 =  1
C-      IA0X00 =  2
C-      IA0Y00 =  3
C-      IA0Z00 =  4
C-      IAXX00 =  5
C-      IAXY00 =  6
C-      IAXZ00 =  7
C-      IAYY00 =  8
C-      IAYZ00 =  9

      IDENA = 0
      call IZERO(ISOFRA, 8)
C     loops over blocks of bra
      do iblock = 1, NUM_BLOCK
        ishell = IDX_SHELL(iblock)
        NHKTA = NHKT(ishell)        !angular quantum number (s=1,p=2,d=3 etc.)
        KHKTA = KHKT(ishell)        !number of spherical (Cartesian) components
        KCKTA = KCKT(ishell)        !number of Cartesian components
        ICA   = LCLASS(ishell)      !class: large component (1), small or Huckel(2), density fitting (0)
        SPHRA = SPHR(ishell)        !if spherical Gaussians
        call LMNVAL(NHKTA, KCKTA, LVALUA, MVALUA, NVALUA)
        NCENTA = NCENT(ishell)      !index of symmetry independent center
        ICENTA = NUCNUM(NCENTA,1)   !total number of nuclei
        MULA   = ISTBAO(ishell)     !stabiliser: basic sym. op. that do not move center
        MULTA  = MULT(MULA)         !multiplicity of center
        NUCA   = NUCO(ishell)       !number of uncontracted functions
        NUMCFA = NUMCF(ishell)      !index of shell in AO-block
        NRCA   = NRCO(ishell)       !number of contracted functions
C FIXME: JSTA, NBCH, NUMCF
        JSTA   = JSTRT(ishell)      !start index of the AO-block
        CORAX  = CENT(ishell,1,1)   !coordinates of center
        CORAY  = CENT(ishell,2,1)
        CORAZ  = CENT(ishell,3,1)
        IDENB0 = 0

        st_bra_prim = JSTA+1
        end_bra_prim = JSTA+NUCA
        st_bra_contr = NUMCFA
        end_bra_contr = NUMCFA+NRCA-1
C       computes symmetry integral pointers for contributions from this block
C       note that at present this assumes all components from a shell are included

        call IZERO(ISOFRB, 8)
        do icontr = 1, NRCA
          do icomp = 1, MXAQN
            INDFA(:,icomp,icontr) = -10 000 000
          end do
        end do
        do icontr = 1, NRCA
          do icomp = 1, KHKTA
            do irep = 0, MAXREP
              if (IBTAND(MULA,IBTXOR(irep,ISYMAO(NHKTB,icomp)))
     &            .eq.0) then
                ISOFRA(irep+1) = ISOFRA(irep+1)+1
                INDFA(irep+1,icomp,icontr) = ISOFRA(irep+1)
              end if
            end do
          end do
        end do

#ifdef GEN1INT_DEBUG
        write(LUPRI,9000) "Center of bra:", NCENTA
        write(LUPRI,9000) "ICENTA: ", ICENTA
        write(LUPRI,9000) "Class of bra: ", ICA
        write(LUPRI,9100) "Coordinates of bra:", CORAX, CORAY, CORAZ
        write(LUPRI,9000) "Angular number of bra:", NHKTA
        write(LUPRI,9000) "Number of bra prim.:", NUCA
        write(LUPRI,9000) "Number of bra contr.:", NRCA
        if (SPHRA) then
          write(LUPRI,9000) "Number of bra spher. components: ", KHKTA
        else
          write(LUPRI,9000) "Number of bra Cart. components: ", KHKTA
        end if
        write(LUPRI,9000) "Number of bra Cart.:", KCKTA
        write(LUPRI,9000) "Start of bra prim.:", st_bra_prim
        write(LUPRI,9000) "End of bra prim.:", end_bra_prim
        write(LUPRI,9200) PRIEXP(st_bra_prim:end_bra_prim)
        write(LUPRI,9000) "Start of bra contr.:", st_bra_contr
        write(LUPRI,9000) "End of bra contr.:", end_bra_contr
        write(LUPRI,9300) PRICCF(st_bra_prim:end_bra_prim,
     &                           st_bra_contr:end_bra_contr)
        do icontr = 1, NRCA
          write(LUPRI,9000) "IA address offsets for bra contraction: ",
     &                      icontr
          do icomp = 1, KHKTA
            write(LUPRI,'(8(1X,I8))')
     &        (INDFA(irep,icomp,icontr),irep=1,MAXREP+1)
          end do
        end do
#endif
C FIXME: lower triangular?
C       loops over shells of ket
        do jblock = 1, NUM_BLOCK

          jshell = IDX_SHELL(jblock)
          LDIAG = ishell.eq.jshell        !lower triangle (.true.) or full square block
          NHKTB = NHKT(jshell)            !angular quantum number (s=1,p=2,d=3 etc.)
          KHKTB = KHKT(jshell)            !number of spherical (Cartesian) components
          KCKTB = KCKT(jshell)            !number of Cartesian components
          ICB   = LCLASS(jshell)          !class: large component (1), small or Huckel(2), density fitting (0)
          SPHRB = SPHR(jshell)            !if spherical Gaussians
          CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
          NCENTB = NCENT(jshell)          !index of symmetry independent center
          NHKTAB = NHKTA+NHKTB
          MULB   = ISTBAO(jshell)         !stabiliser: basic sym. op. that do not move center
          MULTB  = MULT(MULB)             !multiplicity of center
          NUCB   = NUCO(jshell)           !number of uncontracted functions
          NUMCFB = NUMCF(jshell)          !index of shell in AO-block
          NRCB   = NRCO(jshell)           !number of contracted functions
          JSTB   = JSTRT(jshell)          !start index of the AO-block
          CORBX0 = CENT(jshell,1,1)       !coordinates of center
          CORBY0 = CENT(jshell,2,1)
          CORBZ0 = CENT(jshell,3,1)
          KHKTAB = KHKTA*KHKTB
          KCKTAB = KCKTA*KCKTB
          MAB    = IBTOR(MULA,MULB)       !stabliser
          KAB    = IBTAND(MULA,MULB)      !multiplicity to multiplied with the integral
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
          do jcontr = 1, NRCB
            do jcomp = 1, MXAQN
              INDFB(:,jcomp,jcontr) = -10 000 000
            end do
          end do
          do jcontr = 1, NRCB
            do jcomp = 1, KHKTB
              do irep = 0, MAXREP
                if (IBTAND(MULB,IBTXOR(irep,ISYMAO(NHKTB,jcomp)))
     &              .eq.0) then
                  ISOFRB(irep+1) = ISOFRB(irep+1)+1
                  INDFB(irep+1,jcomp,jcontr) = ISOFRB(irep+1)
                end if
              end do
            end do
          end do

          if (.not.doint(ICA,ICB)) cycle

C         loops over symmetry operations

          IDENB = IDENB0-KHKTB
          do isymop = 0, MAXOPR
            if (IBTAND(isymop,MAB).eq.0) then

              IDENB  = IDENB+KHKTB
              ICENTB = NUCNUM(NCENTB,isymop+1)   !total number of nuclei? FIXME
              ONECEN = ICENTA.eq.ICENTB
              SIGNBX = PT(IBTAND(ISYMAX(1,1),isymop))
              SIGNBY = PT(IBTAND(ISYMAX(2,1),isymop))
              SIGNBZ = PT(IBTAND(ISYMAX(3,1),isymop))
              CORBX  = SIGNBX*CORBX0
              CORBY  = SIGNBY*CORBY0
              CORBZ  = SIGNBZ*CORBZ0

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
              select case(trim(prop_name))
              case('.OVERLAP')
C-                call contr_cgto_carmom(NCENTA, CORAX, CORAY, CORAZ,
C-     &            NHKTA-1, NUCA, PRIEXP(st_bra_prim:end_bra_prim),
C-     &            NRCA, PRICCF(st_bra_prim:end_bra_prim,
C-     &                         st_bra_contr:end_bra_contr),
C-     &            NCENTB, CORBX, CORBY, CORBZ,
C-     &            NHKTB-1, NUCB, PRIEXP(st_ket_prim:end_ket_prim),
C-     &            NRCB, PRICCF(st_ket_prim:end_ket_prim,
C-     &                         st_ket_contr:end_ket_contr), DIPORG)
              case('.ECP')
C               loops over number of all ECP centers
C-                do ITYP = 1, NTYECP
C-                  do INONT = 1, NECP(ITYP)
C-                    IATOM = INDECP(ITYP)+INONT-1  !
C-                    KCRS = ITYP                   !
C-                    LCRU = LCR(ITYP)              !
C-                    XCI = CORD(1,IATOM)           !
C-                    YCI = CORD(2,IATOM)
C-                    ZCI = CORD(3,IATOM)
C-                    do jsymop = 0, MAXOPR
C-                      if (IBTAND(jsymop,ISTBNU(IATOM)).eq.0) then
C-                        XC = PT(IBTAND(ISYMAX(1,1),jsymop))*XCI
C-                        YC = PT(IBTAND(ISYMAX(2,1),jsymop))*YCI
C-                        ZC = PT(IBTAND(ISYMAX(3,1),jsymop))*ZCI
C-C                        call PSEUD1(WORK(KGOUT),L1+1,L2+1,EXP1,EXP2)
C-C                        call PSEUD2(WORK(KGOUT),L1+1,L2+1,EXP1,EXP2)
C-                      end if
C-                    end do
C-                  end do
C-                end do
              end select
C             transforms to spherical harmonic basis
C              if (SPHRA) then
C                call cgto_to_sgto()
C              endif
C              if (SPHRB) then
C                call cgto_to_sgto()
C              end if
C             computes expectation values
C-             if (propty) then
C-               collects density and Fock elements
C-               MAXCMP = 0
C-               DO 200 IORBA = IDENA + 1, IDENA + KHKTA
C-                  DO 250 IORBB = IDENB + 1, IDENB + KHKTB
C-                     FAC = D1
C-                     IF (ONECEN .AND. LDIAG .AND. IORBB.NE.IORBA) FAC = DP5
C-                     MAXCMP = MAXCMP + 1
C-                     IORBAB = ITRI(IORBA,IORBB)
C-                     DSHELL(MAXCMP) = FAC*DENMAT(IORBAB)
C-                     FSHELL(MAXCMP) = FAC*FOCMAT(IORBAB)
C-  250             CONTINUE
C-  200          CONTINUE
C-             end if
C FIXME: if undifferentiated, first- and second-order derivatives are calculated seperately
C             transforms integrals to SO basis
C              FULMAT = .true.
C              ANTI   = .false.
C              call SYM1S(ADER,STHMAT(1,2),ISYMOP,MULA,MULB,NHKTA,NHKTB,
C     &                   KHKTA,KHKTB,HKAB,LDIAG,FULMAT,DUMMY,IDUMMY,
C     &                   IPRINT)
C FIXME:
C             writes differentiated integrals on file for the right-hand side
C-             if (getint) then
C-               vals_int() = work()
C-               call WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
C-    &                      ISYMOP,NATOMC,IPRINT,PCM)
C-             end if
            end if
          end do  !loops over symmetry operations
          IDENB0 = IDENB0+KHKTB*MULTB
        end do  !loops over shells of ket
        IDENA = IDENA+KHKTA*MULTA
      end do  !loops over shells of bra
C     cleans
      call IZERO(ISOFRA, 8)
C     Dalton stuff
      call GETTIM(TEND, WEND)
      TIMHER = TEND-TIMHER
      WALHER = WEND-WALHER
      write(LUPRI,'()')
      call TIMTXT(">>>> Total CPU  time used in Gen1Int:",
     &            TIMHER, LUPRI)
      call TIMTXT(">>>> Total wall time used in Gen1Int:",
     &            WALHER, LUPRI)
      call HEADER("End of Gen1Int",-1)
      call QEXIT("gen1int_wrapper")
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
