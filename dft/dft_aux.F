C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
c     DFTINI
C     intialize DFTCOM common block
c /* deck dftini */
      SUBROUTINE DFTINI
#include <implicit.h>
#include <mxcent.h>
      PARAMETER (D0=0D0, D1=1D0)
#include <dftcom.h>
#include <dftacb.h>
      HFXFAC = D1
      WDFTX  = D0
      WDFTC  = D0
      WDFTL  = D0
      WDFTB  = D0
      DFTHR0 = 1.0D-9
      DFTHRL = 1.0D-10
      DFTHRI = 1.0D-13 
      DFTELS = 1.0D-3
      RADINT = 1.0D-13
      ANGINT = 35
      ANGMIN = 15
c     LEBMIN is in the common block but is used only internall in dftgrd
      DFTADD = .FALSE.
      GRDONE = .FALSE.
      DFTRUN = .FALSE.
      DFTPOT = .FALSE.
      DFTORD = .FALSE.
      DFTASC = .FALSE.
      DFTHES = .FALSE.
      DFTHRS = .FALSE.
      NOPRUN = .FALSE.
      DOVWN3 = .FALSE.
      DOVWNI = .TRUE.
C
      DFTIPT = 1.0D-20
      DFTBR1 = 1.0D-20
      DFTBR2 = 1.0D-20
C
      RETURN
      END
c
      SUBROUTINE WRITE_SOI(SOINTS,WRK,LWRK)
c
c     Pawel Salek, 2003-03-18.
c     Write DFT SO integrals to AOPROPER.
c
#include <implicit.h>
#include <inforb.h>
#include <inftap.h>
      DIMENSION SOINTS(N2BASX,3)
c
      CHARACTER*1 CXYZ(3)
      DATA CXYZ /'X','Y','Z'/
      CHARACTER*8 CDUMMY,RTNLBL(2)
      IF(LWRK.LE.NNBASX) CALL QUIT('NO memory in WRITE_SOI')
      CALL GETDAT(RTNLBL(1),CDUMMY)
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
#include <implicit.h>
#include <maxaqn.h>
#include <ccom.h>
#include <mxcent.h>
#include <nuclei.h>
      maxl = NHTYP
      nuniqat = NUCIND
      END

      FUNCTION ishell_cnt()
#include <implicit.h>
#include <maxorb.h>
#include <shells.h>
      ishell_cnt = KMAX
      END
c
      SUBROUTINE get_grid_paras(Lgrdone, Pradint, Iangmin, Iangint)
#include <implicit.h>
#include <dftcom.h>
      LOGICAL LGRDONE
      LGRDONE = GRDONE
      PRADINT = RADINT
      IANGMIN = ANGMIN
      IANGINT = ANGINT
      END
c
      SUBROUTINE set_grid_done
#include <implicit.h>
#include <dftcom.h>
      GRDONE = .TRUE.
      END
      SUBROUTINE get_no_atoms(natom_CNT)
c return total number of atoms.
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
      NATOM_CNT = NATOMS
      END
      SUBROUTINE GET_ATOM_BY_ICENT(ICENT,ZCNT,NAT,NDEG,X,Y,Z)
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
      DIMENSION X(8), Y(8), Z(8)
c
#include <ibtfun.h>
csymmetry is not completely implemented here yet.
      ZCNT   = CHARGE(ICENT)
      MULCNT = ISTBNU(ICENT)
      NDEG   = MULT(MULCNT)
      NAT = 0
      IF(NAMN(ICENT)(1:2).NE.'Gh') THEN
         DO ISYMOP = 0, MAXOPR         
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               NAT = NAT + 1
               X(NAT) = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
               Y(NAT) = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
               Z(NAT) = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
            ENDIF
         ENDDO
      ENDIF
      END
c
      SUBROUTINE FORT_WRT(str, len)
#include <implicit.h>
#include <priunit.h>
      CHARACTER str(len)
      WRITE(LUPRI,*) str
      END
      SUBROUTINE LRAO2MO(CMO,KSYMOP,FAOMAT,FMAT,WORK,LWORK)
#include <implicit.h>
#include <inforb.h>
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
c     dft_get_ao_dens_mat:
c     compute AO density matrix dmat using provided CMO.
c     dmat must be previously allocated.
c
      SUBROUTINE dft_get_ao_dens_mat(cmo,dmat,work,lwork)
#include <implicit.h>
#include <infvar.h>      
      DIMENSION CMO(*), DMAT(*),WORK(LWORK)
      JKEEP = JWOPSY
      JWOPSY = 1
      CALL FCKDEN(1,0,DMAT,DUMMY1,CMO,DUMMY2,WORK,LWORK)
      JWOPSY = JKEEP
      END
c
c dft_get_ao_dens_matab:
c   compute AO active and inactive density matrices using provided CMO.
c
      SUBROUTINE dft_get_ao_dens_matab(cmo,dmata,dmatb,work,lwork)
#include <implicit.h>
#include <inforb.h>
      DIMENSION CMO(*), DMATA(*),DMATB(*),WORK(LWORK)
      KUNIT = 1
      KLAST = KUNIT + N2ASHX
      IF(KLAST.GT.LWORK) CALL QUIT('no mem in dft_get_ao_dens_matab')
      LWRK = LWORK - KLAST + 1
      CALL DUNIT(WORK(KUNIT),NASHT)
      CALL GTDMSO(WORK(KUNIT),CMO,DMATA,DMATB,WORK(KLAST))
      END
c
      SUBROUTINE OUTMAT(A,IRLOW,IRHI,ICLOW,ICHI,IROWS,ICOLS)
#include <implicit.h>
#include <priunit.h>
      CALL OUTPUT(A,IRLOW,IRHI,ICLOW,ICHI,IROWS,ICOLS,0,LUPRI)
      END
c
      SUBROUTINE dftsethf(W)
#include <implicit.h>
#include <dftcom.h>
      HFXFAC = W
      END
C     
      FUNCTION dftgethf()
#include <implicit.h>
#include <dftcom.h>
      dftgethf = HFXFAC
      END
c
      SUBROUTINE DFT_GET_THRESHOLDS(TDFTI,TDFT0)
#include <implicit.h>
#include <dftcom.h>
      TDFTI = DFTHRI
      TDFT0 = DFTHR0
      END
      SUBROUTINE GET_KSYMOP(KS)
#include <implicit.h>
#include <wrkrsp.h>
      KS = KSYMOP
      END
#ifdef VAR_MPI
      SUBROUTINE DFTLRSYNC
#include <implicit.h>
#include <mpif.h>
#include <inforb.h>
#include <infvar.h>
#include <wrkrsp.h>
      CALL MPI_Bcast(NASHT, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NISHT, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(JWOPSY,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(ICMO,  8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(IORB,  8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(IORB,  8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(IASH,  8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NASH,  8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NISH,  8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NORB,  8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(KSYMOP,8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(MULD2H,8*8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      END
#endif
