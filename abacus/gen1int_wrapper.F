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
!...  This file provides the interface of calling Gen1Int from Dalton.
!...
!...  2011-02-13, Bin Gao
!...  * this file does not read information from input anymore, it will only
!...    be called by other subroutines
!...
!...  2010-07-28, Bin Gao
!...  * first version

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
C>    \param is_lao indicates if using London atomic orbitals
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
     &                          is_lao, get_int, wrt_int, vals_int,
     &                          do_exp, ndens, ao_dens,
     &                          get_exp, wrt_exp, vals_expect,
     &                          len_work, dal_work, level_print)
#include "implicit.h"
#include "inforb.h"
      character*(*) prop_name
      integer order_geo
      integer order_mag
      logical is_lao
      logical get_int
      logical wrt_int
      dimension vals_int(*)
      logical do_exp
      integer ndens
      dimension ao_dens(NNBASX,ndens)
      logical get_exp
      logical wrt_exp
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
C     number of geometric differentiated centers
      integer num_cents
C     number of different paths from subroutine \fn(geom_deriv_cent_tree)
      integer num_paths
C     depth of atom to visit
      integer visit_depth
C     selected atom nodes
      integer node_list(MAX_ORDER_GEO)
C     weights of the selected atom nodes
      integer wt_node(MAX_ORDER_GEO)
C     incremental recorder over different paths from subroutine \fn(geom_deriv_cent_tree)
      integer ipath
C     number of geometric derivatives per path
      integer num_dgeo
C     start position of geometric derivatives per path
      integer start_geo
C     incremental recorder over centers
      integer icent
C     Dalton stuff
      call QENTER("gen1int_driver")
      call GETTIM(TIMHER, WALHER)
      call HEADER("Gen1Int calculates "//trim(prop_name),-1)
C     checks the validity of orders
      if (order_geo>MAX_ORDER_GEO) then
        write(LUPRI,1000) "Order of geometric derivatives:", order_geo
        write(LUPRI,1000)
     &    "Maximum order of geometric derivatives (MAX_ORDER_GEO):",
     &    MAX_ORDER_GEO
        call QUIT("Increase MAX_ORDER_GEO in gen1int.h")
      end if
      if (order_mag>MAX_ORDER_MAG) then
        write(LUPRI,1000) "Order of magnetic derivatives:", order_mag
        write(LUPRI,1000)
     &    "Maximum order of magnetic derivatives (MAX_ORDER_MAG):",
     &    MAX_ORDER_MAG
        call QUIT("Increase MAX_ORDER_MAG in gen1int.h")
      end if
C     dumps to check
      if (level_print>=5) then
        write(LUPRI,1000) "Order of geometric derivatives:", order_geo
        write(LUPRI,1000) "Order of magnetic derivaives:", order_mag
        if (is_lao) write(LUPRI,1000) "Using London atomic orbitals"
        if (get_int) write(LUPRI,1000) "Return integrals"
        if (wrt_int) write(LUPRI,1000) "Write integrals on file"
        if (do_exp) then
          write(LUPRI,1000) "Calculate expectation values"
          write(LUPRI,1000) "Number of AO density matrices:", ndens
          if (get_exp) write(LUPRI,1000) "Return expectation values"
          if (wrt_exp)
     &      write(LUPRI,1000) "Write expectation values on file"
        end if
        write(LUPRI,1001) "Length of Dalton workspace:", len_work
      end if
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
        if (level_print>=10) write(LUPRI,1001)
     &    "Number of compositions of differentiated centers:", num_paths
C       the first path (composition of centers), geometric derivatives
C       with respect to the first atom
        start_geo = 1
        call geom_deriv_cent_tree_init(order_geo, num_cents,
     &                                 visit_depth, node_list, wt_node,
     &                                 idx_cent, order_cent, num_dgeo)
C       dumps to check
        if (level_print>=20) then
          write(LUPRI,1002) ipath, start_geo, wt_node(order_geo),
     &                      (idx_cent(icent),"(",order_cent(icent),
     &                       ")",icent=1,wt_node(order_geo))
        end if
C       calculates the geometric and magnetic derivatives
        call gen1int_prop(prop_name, start_geo,
     &                    wt_node(order_geo), order_mag, is_lao,
     &                    get_int, wrt_int, vals_int,
     &                    do_exp, ndens, ao_dens,
     &                    get_exp, wrt_exp, vals_expect,
     &                    len_work, dal_work, level_print)
C       loops over other paths
        do ipath = 2, num_paths
C         computes the position of current geometric derivatives
          start_geo = start_geo+num_dgeo
C         generates the differentiated centers and their orders
          call geom_deriv_cent_tree(NUCDEP, order_geo, num_cents,
     &                              visit_depth, node_list, wt_node,
     &                              idx_cent, order_cent, num_dgeo)
C         dumps to check
          if (level_print>=20) then
            write(LUPRI,1002) ipath, start_geo, wt_node(order_geo),
     &                        (idx_cent(icent),"(",order_cent(icent),
     &                         ")",icent=1,wt_node(order_geo))
          end if
C         calculates the geometric and magnetic derivatives
          call gen1int_prop(prop_name, start_geo,
     &                      wt_node(order_geo), order_mag, is_lao,
     &                      get_int, wrt_int, vals_int,
     &                      do_exp, ndens, ao_dens,
     &                      get_exp, wrt_exp, vals_expect,
     &                      len_work, dal_work, level_print)
        end do
        if (level_print>=10) then
          write(LUPRI,1001) "Current number of geometric derivatives:",
     &                      start_geo+num_dgeo-1
        end if
C     no geometric derivatives
      else
        call gen1int_prop(prop_name, 1, 0, order_mag, is_lao,
     &                    get_int, wrt_int, vals_int,
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
1000  format("gen1int_driver>> ",A,I4)
1001  format("gen1int_driver>> ",A,I12)
1002  format("gen1int_driver>> ","Path",I8,", Start",I10,
     &       ", N_{cent}",I3,", Atom(Order)",20(I3,A,I2,A))
#else
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_driver

C>    \brief calculates integrals/expectation values using Cartesian or
C>           real solid-harmonic Gaussians
C>    \author Bin Gao
C>    \date 2010-07-28
C>    \param prop_name is the name of property integrals to calculated
C>    \param start_geo is the start position of geometric derivatives
C>    \param num_cents is the number of differetiated centers
C>    \param order_mag is the order of magnetic derivaives
C>    \param is_lao indicates if using London atomic orbitals
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
      subroutine gen1int_prop(prop_name, start_geo,
     &                        num_cents, order_mag, is_lao,
     &                        get_int, wrt_int, vals_int,
     &                        do_exp, ndens, ao_dens,
     &                        get_exp, wrt_exp, vals_expect,
     &                        len_work, dal_work, level_print)
#include "implicit.h"
#include "inforb.h"
      character*(*) prop_name
      integer start_geo
      integer num_cents
      integer order_mag
      logical is_lao
      logical get_int
      logical wrt_int
      dimension vals_int(*)
      logical do_exp
      integer ndens
      dimension ao_dens(NNBASX,ndens)
      logical get_exp
      logical wrt_exp
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
C     if in the format of full matrix
      logical full_matrix
C     if anti-symmetric matrix
      logical anti_matrix
C     number of AOs in the AO block
      integer nbra_ao, nket_ao
C     number of contracted Gaussians in the AO block
      integer nbra_cgto, nket_cgto
C     number of geometric derivatives
      integer num_dgeo
C     number of magnetic derivatives
      integer num_dmag
C     number of derivatives
      integer num_deriv
C     number of operators
      integer num_ops
C     number of AO integrals
      integer num_ao_ints
C     number of contracted integrals
      integer num_cints
C     Cartesian powers
      integer power_bra(MXAQN,3), power_ket(MXAQN,3)
C     computes the number of geometric derivatives
      if (num_cents==0) then
        num_dgeo = 1
      else
        call geom_deriv_cent_ndgeo(num_cents, order_cent(1:num_cents),
     &                             num_dgeo)
      end if
C     computes the number of the magnetic derivatives
      num_dmag = (order_mag+1)*(order_mag+2)/2
C     computes the number of derivatives
      num_deriv = num_dgeo*num_dmag
C     used for SYM1S and SYM1N(FIXME)
C-      full_matrix = .true.
C-      anti_matrix = .false.
C     checks if the loop works or not(FIXME)
      IDENA = 0
      call IZERO(ISOFRA, 8)
C     loops over blocks of bra
      do iblock = 1, NUM_BLOCK
        ishell = IDX_SHELL(iblock)
        NHKTA = NHKT(ishell)        !angular quantum number (s=1, p=2, d=3 etc.)
        KHKTA = KHKT(ishell)        !number of spherical (Cartesian) components
        KCKTA = KCKT(ishell)        !number of Cartesian components
        SPHRA = SPHR(ishell)        !if spherical Gaussians
C       generates the Cartesian powers
        call LMNVAL(NHKTA, KCKTA, LVALUA, MVALUA, NVALUA)
        power_bra(:,1) = LVALUA
        power_bra(:,2) = MVALUA
        power_bra(:,3) = NVALUA
        NCENTA = NCENT(ishell)      !index of symmetry independent center
        ICENTA = IPTNUC(NCENTA,0)   !index of this atomic center(FIXME)
        MULA   = ISTBAO(ishell)     !stabiliser: basic sym. op. that do not move center
        MULTA  = MULT(MULA)         !multiplicity of center
        NUCA   = NUCO(ishell)       !number of uncontracted functions
        NUMCFA = NUMCF(ishell)      !index of shell in AO-block
        NRCA   = NRCO(ishell)       !number of contracted functions
        JSTA   = JSTRT(ishell)      !start index of the AO-block
        CORAX  = CENT(ishell,1,1)   !coordinates of center
        CORAY  = CENT(ishell,2,1)
        CORAZ  = CENT(ishell,3,1)
C       checks if the loop works or not(FIXME)
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
C       number of AOs
        nbra_ao = KHKTA*NRCA
C       number of the contracted Gaussians
        nbra_cgto = KCKTA*NRCA
#ifdef GEN1INT_DEBUG
        write(LUPRI,9000)
     &    "Index of symmetry-independent bra center:", NCENTA
        write(LUPRI,9000) "Index of atomic bra center:", ICENTA
        write(LUPRI,9100) "Coordinates of bra:", CORAX, CORAY, CORAZ
        write(LUPRI,9000) "Angular number of bra:", NHKTA
        write(LUPRI,9000) "Number of bra prim.:", NUCA
        write(LUPRI,9000) "Number of bra contr.:", NRCA
        if (SPHRA) then
          write(LUPRI,9000) "Number of bra spher. components:", KHKTA
        else
          write(LUPRI,9000) "Number of bra Cart. components:", KHKTA
        end if
        write(LUPRI,9000) "Number of bra Cart.:", KCKTA
        write(LUPRI,9000) "Start of bra prim.:", JSTA+1
        write(LUPRI,9000) "End of bra prim.:", JSTA+NUCA
        write(LUPRI,9200) PRIEXP(JSTA+1:JSTA+NUCA)
        write(LUPRI,9000) "Start of bra contr.:", NUMCFA
        write(LUPRI,9000) "End of bra contr.:", NUMCFA+NRCA-1
        write(LUPRI,9300) PRICCF(JSTA+1:JSTA+NUCA,NUMCFA:NUMCFA+NRCA-1)
        do icontr = 1, NRCA
          write(LUPRI,9000)
     &      "IA address offsets for bra contraction:", icontr
          do icomp = 1, KHKTA
            write(LUPRI,"(8(1X,I8))")
     &        (BLOCK_INDFA(irep,icomp,icontr),irep=1,MAXREP+1)
          end do
        end do
#endif
C       lower triangular(FIXME)
C       loops over shells of ket
        do jblock = 1, iblock
          jshell = IDX_SHELL(jblock)
          NHKTB = NHKT(jshell)            !angular quantum number (s=1, p=2, d=3 etc.)
          KHKTB = KHKT(jshell)            !number of spherical (Cartesian) components
          KCKTB = KCKT(jshell)            !number of Cartesian components
          SPHRB = SPHR(jshell)            !if spherical Gaussians
C         generates the Cartesian powers
          CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
          power_ket(:,1) = LVALUB
          power_ket(:,2) = MVALUB
          power_ket(:,3) = NVALUB
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
C         number of AOs
          nket_ao = KHKTB*NRCB
C         number of the contracted Gaussians
          nket_cgto = KCKTB*NRCB
C         number of AO integrals
          num_ao_ints = nbra_ao*nket_ao*num_deriv
C         loops over irreducible representations(FIXME)
          IDENB = IDENB0-KHKTB
          do irep = 0, MAXREP
            if (IBTAND(irep,MAB).eq.0) then
C             FIXME
              IDENB  = IDENB+KHKTB
              ICENTB = IPTNUC(NCENTB,irep)   !index of this atomic center(FIXME)
              ONECEN = ICENTA.eq.ICENTB      !if the bra and ket centers coincide
              SIGNBX = PT(IBTAND(ISYMAX(1,1),irep))
              SIGNBY = PT(IBTAND(ISYMAX(2,1),irep))
              SIGNBZ = PT(IBTAND(ISYMAX(3,1),irep))
              CORBX  = SIGNBX*CORBX0         !coordinates of ket center
              CORBY  = SIGNBY*CORBY0
              CORBZ  = SIGNBZ*CORBZ0
#ifdef GEN1INT_DEBUG
              write(LUPRI,9000)
     &          "Index of symmetry-independent ket center:", NCENTB
              write(LUPRI,9000) "Index of atomic ket center: ", ICENTB
              write(LUPRI,9100)
     &          "Coordinates of ket:", CORBX, CORBY, CORBZ
              write(LUPRI,9000) "Angular number of ket:", NHKTB
              write(LUPRI,9000) "Number of ket prim.:", NUCB
              write(LUPRI,9000) "Number of ket contr.:", NRCB
              if (SPHRB) then
                write(LUPRI,9000)
     &            "Number of ket spher. components:", KHKTB
              else
                write(LUPRI,9000)
     &            "Number of ket Cart. components:", KHKTB
              end if
              write(LUPRI,9000) "Number of ket Cart.:", KCKTB
              write(LUPRI,9000) "Start of ket prim.:", JSTB+1
              write(LUPRI,9000) "End of ket prim.:", JSTB+NUCB
              write(LUPRI,9200) PRIEXP(JSTB+1:JSTB+NUCB)
              write(LUPRI,9000) "Start of ket contr.:", NUMCFB
              write(LUPRI,9000) "End of ket contr.:", NUMCFB+NRCB-1
              write(LUPRI,9300)
     &          PRICCF(JSTB+1:JSTB+NUCB,NUMCFB:NUMCFB+NRCB-1)
              do jcontr = 1, NRCB
                write(LUPRI,9000)
     &            "IB address offsets for ket contraction:", jcontr
                do jcomp = 1, KHKTB
                  write(LUPRI,"(8(1X,I8))")
     &              (BLOCK_INDFB(jrep,jcomp,jcontr),jrep=1,MAXREP+1)
                end do
              end do
#endif
C             different property integrals using Gen1Int
C             we may also use the following origins
C               (1) center of mass : CMXYZ
C               (2) operator center: ORIGIN
C               (3) gauge origin   : GAGORG
C               (4) dipole origin  : DIPORG
C             London atomic orbitals at zero field
              if (is_lao) then
C               computes the contracted integrals
                select case(trim(prop_name))
                case default
                  call QUIT("Property"//trim(prop_name)//
     &                      " is not implemented")
                end select
C             normal Gaussians
              else
C               computes the contracted integrals
                select case(trim(prop_name))
                case('.OVERLAP')
                  num_ops = 1
C                 number of contracted integrals
                  num_cints = nbra_cgto*nket_cgto*num_deriv
                  call contr_cgto_carmom(ICENTA, (/CORAX,CORAY,CORAZ/),
     &              NHKTA-1, NUCA, PRIEXP(JSTA+1:JSTA+NUCA),
     &              NRCA, PRICCF(JSTA+1:JSTA+NUCA,NUMCFA:NUMCFA+NRCA-1),
     &              ICENTB, (/CORBX,CORBY,CORBZ/), NHKTB-1,
     &              NUCB, PRIEXP(JSTB+1:JSTB+NUCB), NRCB,
     &              PRICCF(JSTB+1:JSTB+NUCB,NUMCFB:NUMCFB+NRCB-1),
     &              DIPORG, 0, 0, num_cents,
     &              idx_cent(1:2), order_cent(1:2), num_cints, dal_work)
                case('.ECP')
C                 loops over number of all ECP centers
                  do ITYP = 1, NTYECP
                    do INONT = 1, NECP(ITYP)
                      IATOM = INDECP(ITYP)+INONT-1  !
C FIXME                      KCRS = ITYP                   !
C FIXME                      LCRU = LCR(ITYP)              !
                      XCI = CORD(1,IATOM)           !coordinates of ECP center
                      YCI = CORD(2,IATOM)
                      ZCI = CORD(3,IATOM)
                      do jrep = 0, MAXREP
                        if (IBTAND(jrep,ISTBNU(IATOM)).eq.0) then
                          XC = PT(IBTAND(ISYMAX(1,1),jrep))*XCI
                          YC = PT(IBTAND(ISYMAX(2,1),jrep))*YCI
                          ZC = PT(IBTAND(ISYMAX(3,1),jrep))*ZCI
C                          call PSEUD1(WORK(KGOUT),L1+1,L2+1,EXP1,EXP2)
C                          call PSEUD2(WORK(KGOUT),L1+1,L2+1,EXP1,EXP2)
                        end if
                      end do
                    end do
                  end do
                  num_ops = 1
                case default
                  call QUIT("Property"//trim(prop_name)//
     &                      " is not implemented")
                end select
              end if
C             transforms to spherical Gaussian integrals
              if (SPHRA.and.SPHRB) then
                call cgto_to_sgto(NHKTB-1, (/nbra_cgto,KCKTB,
     &                            NRCB*num_deriv*num_ops/), dal_work)
                call cgto_to_sgto(NHKTA-1, (/1,KCKTA,
     &                            NRCA*nket_cgto*num_deriv*num_ops/),
     &                            dal_work)
C             reorders the Cartesian Gaussian integrals
              else
                call reorder_cgto(NHKTB, power_ket,
     &                            (/nbra_cgto,KCKTB,
     &                              NRCB*num_deriv*num_ops/),
     &                            dal_work(1:num_cints),
     &                            dal_work(num_cints+1:
     &                                     num_cints+num_cints))
                call reorder_cgto(NHKTA, power_bra,
     &                            (/1,KCKTA,
     &                              NRCA*nbra_cgto*num_deriv*num_ops/),
     &                            dal_work(num_cints+1:
     &                                     num_cints+num_cints),
     &                            dal_work(1:num_cints))
              end if
C             post-processes the integrals or computes the expectation values
              do icontr = 1, NRCA
                INDFA = BLOCK_INDFA(:,:,icontr)
                do jcontr = 1, NRCB
                  LDIAG = icontr.eq.jcontr  !lower triangle (.true.) or full square block
                  INDFB = BLOCK_INDFB(:,:,jcontr)
C                 computes the expectation values
                  if (do_exp) then
                    call gen1int_expect(nbra_ao, nket_ao, num_dmag,
     &                                  num_dgeo, start_geo,
     &                                  dal_work(1:num_ao_ints),
     &                                  ndens, ao_dens, get_exp,
     &                                  wrt_exp, vals_expect, len_work,
     &                                  dal_work(num_ao_ints+1:),
     &                                  level_print)
                  end if
C                 FIXME
                  if (get_int) then
C                   transforms integrals to SO basis
C-                     call SYM1S(ADER, vals_int, irep, MULA, MULB, NHKTA, NHKTB,
C-     &                          KHKTA, KHKTB, HKAB, LDIAG, FULMAT, DUMMY, IDUMMY,
C-     &                          level_print)
C                   writes differentiated integrals on file for the right-hand side
                    if (wrt_int) then
C-                    call WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
C-     &                          ISYMOP,NATOMC,IPRINT,PCM)
                    end if
                  else if (wrt_int) then
C                   transforms integrals to SO basis
C-                    call SYM1S(ADER, vals_int, irep, MULA, MULB, NHKTA, NHKTB,
C-     &                         KHKTA, KHKTB, HKAB, LDIAG, FULMAT, DUMMY, IDUMMY,
C-     &                         level_print)
C                   writes differentiated integrals on file for the right-hand side
C-                     call WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
C-    &                            ISYMOP,NATOMC,IPRINT,PCM)
                  end if
                end do
              end do
            end if
          end do  !loops over irreducible representations
          IDENB0 = IDENB0+NUCB*KHKTB*MULTB
        end do  !loops over shells of ket
        IDENA = IDENA+NUCA*KHKTA*MULTA
      end do  !loops over shells of bra
C     cleans
      call IZERO(ISOFRA, 8)
      return
#ifdef GEN1INT_DEBUG
9000  format("gen1int_prop>> ",A,I8)
9100  format("gen1int_prop>> ",A,3F16.8)
9200  format("gen1int_prop>> Exponents:",5Es16.8)
9300  format("gen1int_prop>> Contr. Coef.:",5Es16.8)
#endif
#else
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_prop

C>    \brief computes the expectation values
C>    \author Bin Gao
C>    \date 2011-03-24
C>    \param nbra_ao is the number of AOs on bra center
C>    \param nket_ao is the number of AOs on ket center
C>    \param num_dmag is the number of magnetic derivatives
C>    \param num_dgeo is the number of geometric derivatives
C>    \param start_geo is the start position of geometric derivatives
C>    \param vals_int contains the calculated integrals
C>    \param ndens is the number of AO density matrices
C>    \param ao_dens contains the AO density matrices
C>    \param get_exp indicates if getting expectation values back
C>    \param wrt_exp indicates if writing expectation values to file
C>    \param len_work is the length of Dalton workspace
C>    \param dal_work is the Dalton workspace
C>    \param level_print is the level of print
C>    \return vals_expect contains the calculated expectation values
      subroutine gen1int_expect(nbra_ao, nket_ao, num_dmag, num_dgeo,
     &                          start_geo, vals_int, ndens, ao_dens,
     &                          get_exp, wrt_exp, vals_expect,
     &                          len_work, dal_work, level_print)
#include "implicit.h"
#include "inforb.h"
      integer nbra_ao
      integer nket_ao
      integer num_dmag
      integer num_dgeo
      integer start_geo
      dimension vals_int(nbra_ao,nket_ao,num_dmag,num_dgeo)
      integer ndens
      dimension ao_dens(NNBASX,ndens)
      logical get_exp
      logical wrt_exp
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
C-   abacus/her1car.F
C     collects density elements(FIXME)
C-      do iden = 1, ndens
C-        MAXCMP = 0
C-        do IORBA = IDENA+1, IDENA+KHKTA
C-           do IORBB = IDENB+1, IDENB+KHKTB
C-              MAXCMP = MAXCMP+1
C-              MAXAB = max(IORBA,IORBB)
C-              MINAB = min(IORBA,IORBB)
C-              IORBAB = MAXAB*(MAXAB-1)/2+MINAB
C-              if (ONECEN .and. LDIAG .and. IORBB.ne.IORBA) then
C-                dal_work(MAXCMP) = 5.0D-01*ao_dens(IORBAB,iden)
C-              else
C-                dal_work(MAXCMP) = ao_dens(IORBAB,iden)
C-              end if
C-           end do
C-        end do
C-      end do
      return
#else
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_expect

C      subroutine gen1int_fwrite
C        GETONE
C      end subroutine gen1int_fwrite

C      subroutine gen1int_fread
C        GETONE
C      end subroutine gen1int_fread

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
      write(LUPRI,"()")
      write(LUPRI,9000) "Number of AO-blocks:", NUM_BLOCK
      write(LUPRI,9001) IDX_SHELL(1:NUM_BLOCK)
#endif
      call QEXIT('gen1int_block')
      return
1000  format("gen1int_block>> ",A,I8)
#ifdef GEN1INT_DEBUG
9000  format("gen1int_block>> ",A,I8)
9001  format("gen1int_block>> IDX_SHELL:",5I8)
#endif
#else
      call QUIT("Gen1Int is not installed")
#endif
      end subroutine gen1int_block
