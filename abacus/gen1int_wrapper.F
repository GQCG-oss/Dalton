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
C  * this file does not read information from input anymore, it will only
C    be called by other subroutines
C
C  2010-07-28, Bin Gao
C  * first version

C>    \brief main driver of calling Gen1Int by Dalton
C>    \detail modified files include
C>            (1) abacus/herrdn.F: some modifications for processing different ECPs,
C>                                 adds the interface of calling gen1int_block.
C>            (2) abacus/her1car.F: removes the interface of calling old ECP subroutine.
C>            (3) abacus/her1drv.F: adds the interface of calling gen1int_driver().
C>            (4) include/gen1int.h: saves AO block information for Gen1Int
C>    \author Bin Gao
C>    \date 2010-07-28
C>    \param prop_name is the name of property integrals to calculated
C>    \param order_geo is the order of geometric derivatives
C>    \param order_mag is the order of magnetic derivaives
C>    \param get_int indicates if getting integrals back
C>    \param wrt_int indicates if writing integrals to file
C>    \param do_exp indicates if calculating expectation values
C>    \param ndens is the number of AO density matrices
C>    \param ao_dens contains the AO density matrices
C>    \param get_exp indicates if getting expectation values back
C>    \param wrt_exp indicates if writing expectation values to file
C>    \param len_work is the length of Dalton workspace
C>    \param dal_work is the Dalton workspace
C>    \param level_print is the level of print
C>    \return vals_int contains the calculated integrals
C>    \return vals_expect contains the calculated expectation values
      subroutine gen1int_driver(prop_name, order_geo, order_mag,
     &                          get_int, wrt_int, vals_int,
     &                          do_exp, ndens, ao_dens,
     &                          get_exp, wrt_exp, vals_expect,
     &                          len_work, dal_work, level_print)
#include "implicit.h"
#include "inforb.h"
      character*(*) prop_name
      integer order_geo
      integer order_mag
      logical get_int
      logical wrt_int
      dimension vals_int(*)
      logical do_exp
      integer ndens
      dimension ao_dens(NNBASX,ndens)
      dimension vals_expect(*)
      integer len_work
      dimension dal_work(len_work)
      integer level_print
#ifdef BUILD_GEN1INT
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
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
#include "dummy.h"
C     incremental recorder over different paths from subroutine \fn(geom_deriv_cent_tree)
      integer ipath
C     number of geometric derivatives per path
      integer num_dgeo
C     position of geometric derivatives
      integer pos_dgeo
C     maximum number of centers in the integral
      parameter(MAX_NCENT = 4)
C     depth of atom to visit
      integer visit_depth
C     selected atom nodes
      integer node_list(MAX_NCENT)
C     weights of the selected atom nodes
      integer wt_node(MAX_NCENT)
C     indices of atomic centers
      integer idx_cent(MAX_NCENT)
C     order of derivatives of the corresponding atomic centers
      integer order_cent(MAX_NCENT)
C     number of geometric differentiated centers
      integer num_cents
C     number of different paths from subroutine \fn(geom_deriv_cent_tree)
      integer num_paths
C     Dalton stuff
      call QENTER("gen1int_driver")
      call GETTIM(TIMHER, WALHER)
      call HEADER("Gen1Int calculates "//trim(prop_name),-1)
C     dumps to check
      write(LUPRI,1000) "Order of geometric derivatives:", order_geo
      write(LUPRI,1000) "Order of magnetic derivaives:", order_mag
      if (get_int) then
        write(LUPRI,1000) "Return integrals."
      else if (wrt_int) then
        write(LUPRI,1000) "Write integrals on file."
      else
        write(LUPRI,1000) "Neither return nor write integrals."
      end if
      if (do_exp) then
        write(LUPRI,1000) "Calculate expectation values."
        write(LUPRI,1000) "Number of atomic-orbital density matrices:",
     &                    ndens
        if (get_exp) then
          write(LUPRI,1000) "Return expectation values."
        else if (wrt_exp) then
          write(LUPRI,1000) "Write expectation values on file."
        else
          write(LUPRI,1000)
     &      "Neither return nor write expectation values."
        end if
      end if
      write(LUPRI,1001) "Length of Dalton workspace:", len_work
      write(LUPRI,"()")
C     number of centers in the operator
      select case(trim(prop_name))
C     two-center operators
      case('.DSO')
        num_cents = 4
C     operators without any center
      case('.ANGLON','.ANGMOM','.CARMOM','.DIPLEN','.DIPVEL','.DPTOVL',
     &     '.DSUSNL','.KINENE','.MASSVE','.OVERLAP','.S1MAG','.SECMOM',
     &     '.SQHDOL','.SQHDOR','.THETA','.THIRDM','.2NDMM','.3RDMM')
        num_cents = 2
C     one-center operators
      case default
        num_cents = 3
      end select
C     calculates geometric derivatives
      if (order_geo>0) then
C       computes the number of different paths from subroutine \fn(geom_deriv_cent_tree)
C       FIXME: NUCDEP is the number of atoms?
        call geom_deriv_cent_npath(NUCDEP, order_geo,
     &                             num_cents, num_paths)
        write(LUPRI,1001)
     &    "Number of compositions of differentiated centers:", num_paths
C       the first path (composition of centers), geometric derivatives
C       with respect to the first atom
        visit_depth = order_geo
        node_list = 1
        wt_node = 1
        idx_cent = 0
        idx_cent(1) = 1
        order_cent = 0
        order_cent(1) = order_geo
        pos_dgeo = 1
C       calculates the geometric and magnetic derivatives
        call gen1int_prop(prop_name, pos_dgeo, wt_node(order_geo), 
     &                    idx_cent, order_cent, order_mag,
     &                    get_int, wrt_int, vals_int,
     &                    do_exp, ndens, ao_dens,
     &                    get_exp, wrt_exp, vals_expect,
     &                    len_work, dal_work, level_print)
C       loops over other paths
        do ipath = 2, num_paths
C         computes the position of current geometric derivatives
          call geom_deriv_cent_ndgeo(wt_node(order_geo),
     &                               order_cent(1:wt_node(order_geo)),
     &                               num_dgeo)
          pos_dgeo = pos_dgeo+num_dgeo
C         generates the differentiated centers and their orders
          call geom_deriv_cent_tree(NUCDEP, order_geo, num_cents,
     &                              visit_depth, node_list, wt_node,
     &                              idx_cent, order_cent)
C         calculates the geometric and magnetic derivatives
          call gen1int_prop(prop_name, pos_dgeo, wt_node(order_geo), 
     &                      idx_cent, order_cent, order_mag,
     &                      get_int, wrt_int, vals_int,
     &                      do_exp, ndens, ao_dens,
     &                      get_exp, wrt_exp, vals_expect,
     &                      len_work, dal_work, level_print)
        end do
C     no geometric derivatives
      else
        call gen1int_prop(prop_name, 1, 0, idx_cent, order_cent,
     &                    order_mag, get_int, wrt_int, vals_int,
     &                    do_exp, ndens, ao_dens,
     &                    get_exp, wrt_exp, vals_expect,
     &                    len_work, dal_work, level_print)
      end if
C     Dalton stuff
      call GETTIM(TEND, WEND)
      TIMHER = TEND-TIMHER
      WALHER = WEND-WALHER
      write(LUPRI,"()")
      call TIMTXT(">>>> Total CPU  time used in Gen1Int:",
     &            TIMHER, LUPRI)
      call TIMTXT(">>>> Total wall time used in Gen1Int:",
     &            WALHER, LUPRI)
      call HEADER("End of Gen1Int",-1)
      call QEXIT("gen1int_driver")
      return
1000  format("GEN1INT_INFO>> ",A,I4)
1001  format("GEN1INT_INFO>> ",A,I12)
#else
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_driver

C>    \brief calculates integrals/expectation values using Cartesian or
C>           real solid-harmonic Gaussians
C>    \author Bin Gao
C>    \date 2010-07-28
C>    \param prop_name is the name of property integrals to calculated
C>    \param order_geo is the order of geometric derivatives
C>    \param order_mag is the order of magnetic derivaives
C>    \param get_int indicates if getting integrals back
C>    \param wrt_int indicates if writing integrals to file
C>    \param do_exp indicates if calculating expectation values
C>    \param ndens is the number of AO density matrices
C>    \param ao_dens contains the AO density matrices
C>    \param get_exp indicates if getting expectation values back
C>    \param wrt_exp indicates if writing expectation values to file
C>    \param len_work is the length of Dalton workspace
C>    \param dal_work is the Dalton workspace
C>    \param level_print is the level of print
C>    \return vals_int contains the calculated integrals
C>    \return vals_expect contains the calculated expectation values
      subroutine gen1int_prop(prop_name, order_geo, order_mag,
     &                        get_int, wrt_int, vals_int,
     &                        do_exp, ndens, ao_dens,
     &                        get_exp, wrt_exp, vals_expect,
     &                        len_work, dal_work, level_print)
#include "implicit.h"
#include "inforb.h"
      character*(*) prop_name
      integer order_geo
      integer order_mag
      logical get_int
      logical wrt_int
      dimension vals_int(*)
      logical do_exp
      integer ndens
      dimension ao_dens(NNBASX,ndens)
      dimension vals_expect(*)
      integer len_work
      dimension dal_work(len_work)
      integer level_print
#ifdef BUILD_GEN1INT
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
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
#include "dummy.h"
C     symmetry integral pointers for AO-blocks
      integer BLOCK_INDFA(8,MXAQN,MXCONT), BLOCK_INDFB(8,MXAQN,MXCONT)
C     incremental recorders over spherical (Cartesian) components
      integer icomp, jcomp
C     incremental recorders over contractions
      integer icontr, jcontr
C     incremental recorders over AO-shells
      integer ishell, jshell
C     incremental recorders over AO-blocks
      integer iblock, jblock
C     incremental recorder over different paths from subroutine \fn(geom_deriv_cent_tree)
      integer ipath
C     if in the format of full matrix
      logical full_matrix
C     if anti-symmetric matrix
      logical anti_matrix
C     number of AO shells in the AO block
      integer nbra_shell, nket_shell
C     maximum number of centers in the integral
      parameter(MAX_NCENT = 4)
C     depth of atom to visit
      integer visit_depth
C     selected atom nodes
      integer node_list(MAX_NCENT)
C     weights of the selected atom nodes
      integer wt_node(MAX_NCENT)
C     indices of atomic centers
      integer idx_cent(MAX_NCENT)
C     order of derivatives of the corresponding atomic centers
      integer order_cent(MAX_NCENT)
C     number of geometric differentiated centers
      integer num_cents
C     number of different paths from subroutine \fn(geom_deriv_cent_tree)
      integer num_paths
C     number of geometric derivatives
      integer num_dgeo
C     number of magnetic derivatives
      integer num_dmag
C     number of derivatives
      integer num_deriv

C     number of the magnetic derivatives
      num_dmag = (order_mag+1)*(order_mag+2)/2

      num_deriv = num_dgeo*num_dmag

C     used for SYM1S and SYM1N(FIXME)
C-      full_matrix = .true.
C-      anti_matrix = .false.

C     FIXME
      IDENA = 0
      call IZERO(ISOFRA, 8)
C     loops over blocks of bra
      do iblock = 1, NUM_BLOCK
        ishell = IDX_SHELL(iblock)
        NHKTA = NHKT(ishell)        !angular quantum number (s=1,p=2,d=3 etc.)
        KHKTA = KHKT(ishell)        !number of spherical (Cartesian) components
        KCKTA = KCKT(ishell)        !number of Cartesian components
        SPHRA = SPHR(ishell)        !if spherical Gaussians
C       generates the Cartesian powers
        call LMNVAL(NHKTA, KCKTA, LVALUA, MVALUA, NVALUA)
        NCENTA = NCENT(ishell)      !index of symmetry independent center
        ICENTA = IPTNUC(NCENTA,1)   !index of this atomic center(FIXME)
        MULA   = ISTBAO(ishell)     !stabiliser: basic sym. op. that do not move center
        MULTA  = MULT(MULA)         !multiplicity of center
        NUCA   = NUCO(ishell)       !number of uncontracted functions
        NUMCFA = NUMCF(ishell)      !index of shell in AO-block
        NRCA   = NRCO(ishell)       !number of contracted functions
        JSTA   = JSTRT(ishell)      !start index of the AO-block
        CORAX  = CENT(ishell,1,1)   !coordinates of center
        CORAY  = CENT(ishell,2,1)
        CORAZ  = CENT(ishell,3,1)
C
        IDENB0 = 0
C       computes symmetry integral pointers for contributions from this block
C       note that at present this assumes all components from a shell are included
        call IZERO(ISOFRB, 8)
        do icontr = 1, NRCA
          do icomp = 1, MXAQN
            BLOCK_INDFA(:,icomp,icontr) = -10 000 000
          end do
        end do
        do icontr = 1, NRCA
          do icomp = 1, KHKTA
            do irep = 0, MAXREP
              if (IBTAND(MULA,IBTXOR(irep,ISYMAO(NHKTA,icomp)))
     &            .eq.0) then
                ISOFRA(irep+1) = ISOFRA(irep+1)+1
                BLOCK_INDFA(irep+1,icomp,icontr) = ISOFRA(irep+1)
              end if
            end do
          end do
        end do
C       number of the AO shells
        nbra_shell = KCKTA*NRCA
#ifdef GEN1INT_DEBUG
        write(LUPRI,9000) "Index of symmetry-independent bra center:",
     &                    NCENTA
        write(LUPRI,9000) "Index of atomic bra center: ", ICENTA
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
        write(LUPRI,9000) "Start of bra prim.:", JSTA+1
        write(LUPRI,9000) "End of bra prim.:", JSTA+NUCA
        write(LUPRI,9200) PRIEXP(JSTA+1:JSTA+NUCA)
        write(LUPRI,9000) "Start of bra contr.:", NUMCFA
        write(LUPRI,9000) "End of bra contr.:", NUMCFA+NRCA-1
        write(LUPRI,9300) PRICCF(JSTA+1:JSTA+NUCA,NUMCFA:NUMCFA+NRCA-1)
        do icontr = 1, NRCA
          write(LUPRI,9000) "IA address offsets for bra contraction: ",
     &                      icontr
          do icomp = 1, KHKTA
            write(LUPRI,'(8(1X,I8))')
     &        (BLOCK_INDFA(irep,icomp,icontr),irep=1,MAXREP+1)
          end do
        end do
#endif
C       lower triangular(FIXME)
C       loops over shells of ket
        do jblock = 1, iblock
          jshell = IDX_SHELL(jblock)
          NHKTB = NHKT(jshell)            !angular quantum number (s=1,p=2,d=3 etc.)
          KHKTB = KHKT(jshell)            !number of spherical (Cartesian) components
          KCKTB = KCKT(jshell)            !number of Cartesian components
          SPHRB = SPHR(jshell)            !if spherical Gaussians
          CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
          NCENTB = NCENT(jshell)          !index of symmetry independent center
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
C         generates the Cartesian powers
          call LMNVAL(NHKTB, KCKTB, LVALUB, MVALUB, NVALUB)
C         computes symmetry integral pointers for contributions from this block
C         note that at present this assumes all components from a shell are included
          do jcontr = 1, NRCB
            do jcomp = 1, MXAQN
              BLOCK_INDFB(:,jcomp,jcontr) = -10 000 000
            end do
          end do
          do jcontr = 1, NRCB
            do jcomp = 1, KHKTB
              do irep = 0, MAXREP
                if (IBTAND(MULB,IBTXOR(irep,ISYMAO(NHKTB,jcomp)))
     &              .eq.0) then
                  ISOFRB(irep+1) = ISOFRB(irep+1)+1
                  BLOCK_INDFB(irep+1,jcomp,jcontr) = ISOFRB(irep+1)
                end if
              end do
            end do
          end do
C         number of the AO shells
          nket_shell = KCKTB*NRCB
C         loops over symmetry operations(FIXME)
          IDENB = IDENB0-KHKTB
          do isymop = 0, MAXOPR
            if (IBTAND(isymop,MAB).eq.0) then
C             FIXME
              IDENB  = IDENB+KHKTB
              ICENTB = IPTNUC(NCENTB,isymop+1)   !index of this atomic center(FIXME)
              ONECEN = ICENTA.eq.ICENTB

              SIGNBX = PT(IBTAND(ISYMAX(1,1),isymop))
              SIGNBY = PT(IBTAND(ISYMAX(2,1),isymop))
              SIGNBZ = PT(IBTAND(ISYMAX(3,1),isymop))
              CORBX  = SIGNBX*CORBX0
              CORBY  = SIGNBY*CORBY0
              CORBZ  = SIGNBZ*CORBZ0

C             different property integrals using Gen1Int
C             we may also use the following origins
C             (1) center of mass : CMXYZ
C             (2) operator center: ORIGIN
C             (3) gauge origin   : GAGORG
C             (4) dipole origin  : DIPORG
              select case(trim(prop_name))
              case('.OVERLAP')
C               computes the contracted integrals
                call contr_cgto_carmom(ICENTA, (/CORAX,CORAY,CORAZ/),
     &            NHKTA-1, NUCA, PRIEXP(JSTA+1:JSTA+NUCA),
     &            NRCA, PRICCF(JSTA+1:JSTA+NUCA,NUMCFA:NUMCFA+NRCA-1),
     &            ICENTB, (/CORBX,CORBY,CORBZ/), NHKTB-1,
     &            NUCB, PRIEXP(JSTB+1:JSTB+NUCB), NRCB,
     &            PRICCF(JSTB+1:JSTB+NUCB,NUMCFB:NUMCFB+NRCB-1),
     &            DIPORG, 0, dal_work)
              case('.ECP')
C               loops over number of all ECP centers
                do ITYP = 1, NTYECP
                  do INONT = 1, NECP(ITYP)
                    IATOM = INDECP(ITYP)+INONT-1  !
C FIXME                    KCRS = ITYP                   !
C FIXME                    LCRU = LCR(ITYP)              !
                    XCI = CORD(1,IATOM)           !coordinates of ECP center
                    YCI = CORD(2,IATOM)
                    ZCI = CORD(3,IATOM)
                    do jsymop = 0, MAXOPR
                      if (IBTAND(jsymop,ISTBNU(IATOM)).eq.0) then
                        XC = PT(IBTAND(ISYMAX(1,1),jsymop))*XCI
                        YC = PT(IBTAND(ISYMAX(2,1),jsymop))*YCI
                        ZC = PT(IBTAND(ISYMAX(3,1),jsymop))*ZCI
C                        call PSEUD1(WORK(KGOUT),L1+1,L2+1,EXP1,EXP2)
C                        call PSEUD2(WORK(KGOUT),L1+1,L2+1,EXP1,EXP2)
                      end if
                    end do
                  end do
                end do
              end select
C             transforms to spherical Gaussian integrals
              if (SPHRB) then
                call cgto_to_sgto(NHKTB-1, (/nbra_shell,KCKTB,
     &                            NRCB*num_deriv/), dal_work)
              end if
              if (SPHRA) then
                call cgto_to_sgto(NHKTA-1, (/1,KCKTA,
     &                            NRCA*nket_shell*num_deriv/),
     &                            dal_work)
              endif
C
              do icontr = 1, NRCA
                INDFA = BLOCK_INDFA(:,:,icontr)
                do jcontr = 1, NRCB
                  LDIAG = icontr.eq.jcontr  !lower triangle (.true.) or full square block
                  INDFB = BLOCK_INDFB(:,:,jcontr)
C                 computes expectation values
                  if (do_exp) then
     &              call gen1int_expect(ndens, ao_dens,
     &                                  get_exp, wrt_exp, vals_expect,
     &                                  len_work, dal_work, level_print)
                  end if
C                 FIXME
                  if (get_int) then
C                  transforms integrals to SO basis
C-                     call SYM1S(ADER, vals_int, isymop, MULA, MULB, NHKTA, NHKTB,
C-     &                          KHKTA, KHKTB, HKAB, LDIAG, FULMAT, DUMMY, IDUMMY,
C-     &                          level_print)
                  else if (wrt_int) then
C                     writes differentiated integrals on file for the right-hand side
C-                    call WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
C-    &                           ISYMOP,NATOMC,IPRINT,PCM)
                  end if
                end do
              end do
            end if
          end do  !loops over symmetry operations
          IDENB0 = IDENB0+NUCB*KHKTB*MULTB
        end do  !loops over shells of ket
        IDENA = IDENA+NUCA*KHKTA*MULTA
      end do  !loops over shells of bra
C     cleans
      call IZERO(ISOFRA, 8)
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
      end subroutine gen1int_prop

C>    \brief computes the expectation values
C>    \author Bin Gao
C>    \date 2010-12-06
      subroutine gen1int_expect(ndens, ao_dens,
     &                          get_exp, wrt_exp, vals_expect,
     &                          len_work, dal_work, level_print)
C     collects density elements
      do iden = 1, ndens
        MAXCMP = 0
        do IORBA = IDENA+1, IDENA+KHKTA
           do IORBB = IDENB+1, IDENB+KHKTB
              MAXCMP = MAXCMP+1
              MAXAB = max(IORBA,IORBB)
              MINAB = min(IORBA,IORBB)
              IORBAB = MAXAB*(MAXAB-1)/2+MINAB
              if (ONECEN .and. LDIAG .and. IORBB.ne.IORBA) then
                dal_work(MAXCMP) = 5.0D-01*ao_dens(IORBAB,iden)
              else
                dal_work(MAXCMP) = ao_dens(IORBAB,iden)
              end if
           end do
        end do
      end do
      return
      end subroutine gen1int_expect

C>    \brief sets up AO-blocks used for Gen1Int
C>    \author Bin Gao
C>    \date 2010-12-06
      subroutine gen1int_block
#ifdef BUILD_GEN1INT
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "shells.h"
#include "aovec.h"
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
