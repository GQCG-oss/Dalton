!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!

C----------------------------------------------------------------------------
      SUBROUTINE DFTINI
C----------------------------------------------------------------------------
c     DFTINI
C     intialize DFTCOM common block
C------------------------------------------------------------------------------
#include "implicit.h"
#include "mxcent.h"
      PARAMETER (D0=0D0, D1=1D0)
#include "gnrinf.h"
#include "dftcom.h"
#include "dftacb.h"
#include "dftd.h"
      IPRDFT = IPRUSR
      HFXFAC = D1
      HFXMU  = D0
      HFXSET = .FALSE. !hjaaj: HFXFAC has not been set (i.e. redefined) by user
      WDFTX  = D0
      WDFTC  = D0
      WDFTL  = D0
      WDFTB  = D0
      WDFTMP = D0
      DFTHR0 = 1.0D-9
      DFTHRL = 1.0D-10
      DFTHRI = 2.0D-12
      DFTELS = 1.0D-3
      RADINT = 1.0D-13
      ANGINT = 35
      ANGMIN = 15
c     LEBMIN is in the common block but is used only internally in make_dftgrid
      DFTADD = .FALSE.
      DFTGRID_DONE = .FALSE.
      DFTGRID_DONE_OLD = .FALSE.
      DFTRUN = .FALSE.
      DFTPOT = .FALSE.
      DFTORD = .FALSE.
      DFTASC = .FALSE.
      DFTHES = .FALSE.
      DFTHRS = .FALSE.
      NOPRUN = .FALSE.
C
      DFTIPT = 1.0D-20
      DFTBR1 = 1.0D-20
      DFTBR2 = 1.0D-20
C
C     Common block variables used when SR-DFT hybrids.
C
      SRDFTRUN = .FALSE.
      IF (ERFEXP(1)) THEN
C        default scaling factor for ERFGAU, as proposed by Savin & Toulouse
C        (COPFAC becomes approx. 1/3.375)
         COPFAC   = 1.0D0+6.0D0*SQRT(3.0D0)
         COPFAC   = 1.0D0/SQRT(COPFAC)
      ELSE IF (ERFEXP(2)) THEN
         COPFAC = 1.0D0/2.25D0 ! mu=1 for erf is mu=2.25 for new erfexp(2)
         ! COPFAC*MU is used when the erf functionals are used for new erfexp(2)
      ELSE
         COPFAC   = D1
      ENDIF

      IWINT    = 0
      SRXFUN   = '      '
      SRCFUN   = '      '
      DOSRX_LDA = .FALSE.
      DOSRX_GGA = .FALSE.
      DOSRX_WIB = .FALSE.
      DOSRGGA2 = .FALSE.
      DOSRBCK  = .FALSE.
      DOHFEXCH = .FALSE.

C JT 12-02-05 beg
      ISJT     = .FALSE.
      DOSRX_PBEHSE = .FALSE.
      DOSRX_PBETCS = .FALSE.
      DOSRC_PBETCS = .FALSE.
      DOSRX_PBERI  = .FALSE.
      DOSRC_PBERI  = .FALSE.
C JT 12-02-05 end
      DOSRC_LYPRI  = .FALSE.

C JT 11-08-09 beg
      DOSRX_PBEGWS = .FALSE.
      DOSRC_PBEGWS = .FALSE.
C JT 11-08-09 end

C JT 17-08-09 beg
      DOSRX_LDA_S = .FALSE.
      DOSRC_LDA_S = .FALSE.

C EDH 10-02-16
      DOSRX_PBEGWS_S = .FALSE.
      DOSRC_PBEGWS_S = .FALSE.

C SR correlation functionals
      DOSRC_LDA = .FALSE.
      DOSRC_GGA = .FALSE.
      DOSRC_WIB = .FALSE.
      DOSRC_MULOCAL(0:3) = .FALSE.
      DOSRLYPT = .FALSE.

C MNP 30.01.12
      SRCMULOFAC = .FALSE.
      DSLOCALFAC = .FALSE.
      DOSRC_MULOC_GGA = .FALSE.
      DOSRC_MULOD_GGA = .FALSE.
      DOSRC_MULOE_GGA = .FALSE.
      DOSRC_PBELO = .FALSE.
      HEAVISIDE_PVALUE = 4.0d0

CAMT --- DFT-D correction initialize
      DODFTD     = .FALSE.
      DFTD_TEST  = .FALSE.
      L_INP_D2PAR= .FALSE.
      L_INP_D3PAR= .FALSE.
      DO_DFTD2   = .FALSE.
      DO_DFTD3   = .FALSE.
      DO_BJDAMP  = .FALSE.
      DO_3BODY   = .FALSE.
      D2_s6_inp  = 0.0d0
      D2_alp_inp = 0.0d0
      D2_rs6_inp = 0.0d0
      D3_s6_inp  = 0.0d0
      D3_alp_inp = 0.0d0
      D3_rs6_inp = 0.0d0
      D3_rs18_inp= 0.0d0
      D3_s18_inp = 0.0d0

C DFTVXC Init
      LDFTVXC = .FALSE.
C DFTAC Init
      LGRAC = .FALSE.
      LLIN  = .FALSE.
      LTAN  = .FALSE.
      DOLB94 = .FALSE.

      RETURN
      END
c
      SUBROUTINE WRITE_SOI(SOINTS,WRK,LWRK)
c
c     Pawel Salek, 2003-03-18.
c     Write DFT SO integrals to AOPROPER.
c
#include "implicit.h"
#include "dummy.h"
#include "inforb.h"
#include "inftap.h"
      DIMENSION SOINTS(N2BASX,3)
c
      CHARACTER*1 CXYZ(3)
      DATA CXYZ /'X','Y','Z'/
      CHARACTER*8 RTNLBL(2)
      IF(LWRK.LE.NNBASX) CALL QUIT('NO memory in WRITE_SOI')
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)='ANTISYMM'
      CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      DO I=1,3
         CALL DGETAP(NBAST,SOINTS(1,I),WRK)
         CALL WRTPRO(WRK,NNBASX,CXYZ(I)//'1DFT-SO',RTNLBL)
      END DO
      CALL GPCLOSE(LUPROP,'KEEP')
      RETURN
      END
c
      SUBROUTINE GET_MAXL_NUCIND(maxl,nuniqat)
#include "implicit.h"
#include "maxaqn.h"
#include "ccom.h"
#include "mxcent.h"
#include "nuclei.h"
      maxl = NHTYP
      nuniqat = NUCIND
      END

      FUNCTION ishell_cnt()
#include "implicit.h"
#include "maxorb.h"
#include "shells.h"
      ishell_cnt = KMAX
      END
c
      SUBROUTINE get_grid_paras(LDFTGRID_DONE,Pradint,Iangmin,Iangint)
#include "implicit.h"
#include "dftcom.h"
      LOGICAL LDFTGRID_DONE
      LDFTGRID_DONE = DFTGRID_DONE
      PRADINT = RADINT
      IANGMIN = ANGMIN
      IANGINT = ANGINT
      END
c
      SUBROUTINE set_grid_done
#include "implicit.h"
#include "dftcom.h"
      DFTGRID_DONE = .TRUE.
      END
      SUBROUTINE get_no_atoms(natom_CNT)
c return total number of atoms.
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
      NATOM_CNT = NATOMS
      END
      SUBROUTINE GET_ATOM_BY_ICENT(ICENT,ZCNT,NAT,NDEG,X,Y,Z)
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
      DIMENSION X(8), Y(8), Z(8)
c

csymmetry is not completely implemented here yet.
      ZCNT   = CHARGE(ICENT)
      MULCNT = ISTBNU(ICENT)
      NDEG   = MULT(MULCNT)
      NAT = 0
      IF(NAMN(ICENT)(1:2).NE.'Gh') THEN
         DO ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
               NAT = NAT + 1
               X(NAT) = PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
               Y(NAT) = PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
               Z(NAT) = PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
            ENDIF
         ENDDO
      ENDIF
      END
c
      SUBROUTINE FORT_WRT(str, len)
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(len) str
      WRITE(LUPRI,'(1X,A)') str
      END
      INTEGER FUNCTION ISETKSYMOP(IKSYM)
#include "implicit.h"
#include "wrkrsp.h"
      ISETKSYMOP = KSYMOP
      KSYMOP = IKSYM
      END
      SUBROUTINE LRAO2MO(CMO,KSYMOP,FAOMAT,FMAT,WORK,LWORK)
#include "implicit.h"
#include "inforb.h"
      DIMENSION CMO(*), FAOMAT(*), FMAT(*),WORK(LWORK)
      PARAMETER (D2 = 2D0)
      KDA1 = 1
      KLST = KDA1 + N2ORBX
      LWRK = LWORK - KLST + 1
      IF(KLST.GT.LWORK) CALL QUIT('not enough memory in LRAO2MO')
      CALL DZERO(WORK(KDA1),N2ORBX)
      DO ISYM = 1, NSYM
         JSYM  = MULD2H(ISYM,KSYMOP)
         NORBI = NORB(ISYM)
         NORBJ = NORB(JSYM)
         IF (NORBI.GT.0 .AND. NORBJ.GT.0) THEN
            CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &                 FAOMAT,NBAS,NBAST,WORK(KDA1),NORB,
     &                 NORBT,WORK(KLST),LWRK)
         END IF
      END DO
      CALL DAXPY(N2ORBX,D2,WORK(KDA1),1,FMAT,1)
      END
      SUBROUTINE dft_get_ao_dens_mat(cmo,dmat,work,lwork)
c
c     dft_get_ao_dens_mat:
c     compute AO density matrix dmat using provided CMO.
c     dmat must be previously allocated.
c
#include "implicit.h"
#include "dummy.h"
#include "infvar.h"
      DIMENSION CMO(*), DMAT(*),WORK(LWORK)
      JKEEP = JWOPSY
      JWOPSY = 1
      CALL FCKDEN(1,0,DMAT,VDUMMY,CMO,VDUMMY,WORK,LWORK)
      JWOPSY = JKEEP
      END
      SUBROUTINE dft_get_ao_dens_matab(cmo,dmata,dmatb,work,lwork)
c
c dft_get_ao_dens_matab:
c   compute AO active and inactive density matrices using provided CMO.
c
#include "implicit.h"
#include "inforb.h"
      DIMENSION CMO(*), DMATA(*),DMATB(*),WORK(LWORK)
      KUDV = 1
      KLAST = KUDV + N2ASHX
      IF(KLAST.GT.LWORK) CALL QUIT('no mem in dft_get_ao_dens_matab')
      LWRK = LWORK - KLAST + 1
      CALL DUNIT(WORK(KUDV),NASHT)
C     ... high spin, i.e. active density matrix UDV is a unit matrix.
      CALL GTDMSO(WORK(KUDV),CMO,DMATA,DMATB,WORK(KLAST))
      END
c
      SUBROUTINE OUTMAT(A,IRLOW,IRHI,ICLOW,ICHI,IROWS,ICOLS)
#include "implicit.h"
#include "priunit.h"
      CALL OUTPUT(A,IRLOW,IRHI,ICLOW,ICHI,IROWS,ICOLS,0,LUPRI)
      END
c
      SUBROUTINE dftsethf(W)
#include "implicit.h"
#include "dftcom.h"
      HFXFAC = W
      END
C
      FUNCTION dftgethf()
#include "implicit.h"
#include "dftcom.h"
      dftgethf = HFXFAC
      END

      SUBROUTINE dftsetmp2(W)
#include "implicit.h"
#include "dftcom.h"
      WDFTMP = W
      END
C
      FUNCTION dftgetmp2()
#include "implicit.h"
#include "dftcom.h"
      dftgetmp2 = WDFTMP
      END
c
      SUBROUTINE dftsetcam(W,BETA)
#include "implicit.h"
#include "dftcom.h"
      HFXMU = W
      HFXATT = BETA
      END
c
      SUBROUTINE DFT_GET_THRESHOLDS(TDFTI,TDFT0)
#include "implicit.h"
#include "dftcom.h"
      TDFTI = DFTHRI
      TDFT0 = DFTHR0
      END
#ifdef VAR_MPI
      SUBROUTINE DFTLRSYNC
#include "implicit.h"
#include "mpif.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
#include "infvar.h"
#include "wrkrsp.h"
      CALL MPI_Bcast(NASHT, 1,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NISHT, 1,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(JWOPSY,1,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(ICMO,  8,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(IORB,  8,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(IASH,  8,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NASH,  8,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NISH,  8,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NORB,  8,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(KSYMOP,8,  my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(MULD2H,8*8,my_MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!     print *, 'done with sync mytid',mytid
      END
#endif
      SUBROUTINE GETSHELLSCNT(nshells)
#include "implicit.h"
#include "maxorb.h"
#include "shells.h"
      nshells = kmax
      END
      SUBROUTINE GETSHELLNO(ISHELL,ICNT,IL,X,Y,Z,MXCONTR,
     &                      COEFF, ALPHA)
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
c     RETURNS INFORMATION ABOUT SHELL ISHELL
c
      DIMENSION COEFF(MXCONTR), ALPHA(MXCONTR)
#include "inforb.h"
#include "onecom.h"
#include "lmns.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "primit.h"
#include "sphtrm.h"
#include "orgcom.h"
      NHKTA  = NHKT(ISHELL)
      KHKTA  = KHKT(ISHELL)
      KCKTA  = KCKT(ISHELL)
      SPHRA  = SPHR(ISHELL)
      NUMCFA = NUMCF(ISHELL)
      JSTA   = JSTRT(ISHELL)
      NUCA   = NUCO(ISHELL)
c        dump shell:
      ICNT = 0
      DO IC = JSTA +1, JSTA + NUCA
         IF(PRICCF(IC,NUMCFA).NE.0D0) ICNT = ICNT +1
      END DO
      IL = NHKTA - 1
      X  = CENT(ISHELL,1,1)
      Y  = CENT(ISHELL,2,1)
      Z  = CENT(ISHELL,3,1)
      IF(ICNT.LE.MXCONTR) THEN
         IDX = 1
         DO IC = JSTA +1, JSTA + NUCA
            IF(PRICCF(IC,NUMCFA).NE.0D0) THEN
               COEFF(IDX) = PRICCF(IC, NUMCFA)
               ALPHA(IDX) = PRIEXP(IC)
               IDX = IDX + 1
            END IF
         END DO
      END IF
      END
