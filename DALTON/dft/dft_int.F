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
      SUBROUTINE DFTINT(DMAT,NDMAT,NGEODRV,DOLND,WORK,LWORK,CB,
     &                  CBDATA,ELECTRONS)
#include "implicit.h"
C
C Used
C  dftcom.h : IPRDFT, ?
C
#include "priunit.h"
#include "inforb.h"
#include "mxcent.h"
#include "dftcom.h"
#include "nuclei.h"
#include "dftacb.h"
#include "maxorb.h"
      DIMENSION DMAT(NBAST,NBAST,NDMAT), WORK(LWORK), CBDATA(*)
      EXTERNAL CB
      LOGICAL DOLND
C
      CALL QENTER('DFTINT')
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK
      LGAO   = MXBLLEN*NBAST*22
      LDMGAO = MXBLLEN*NBAST
      CALL MEMGET2('REAL','   GAO',KGAO  ,LGAO  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL',' DMGAO',KDMGAO,LDMGAO,WORK,KFREE,LFREE)
CAMT  Allocations for DFT-AC terms
      LDAOx  = MXBLLEN*NBAST
      LDAOy  = MXBLLEN*NBAST
      LDAOz  = MXBLLEN*NBAST
      LHESA  = 6*MXBLLEN
      LDST   = NATOMS
      LVFA   = MXBLLEN
      LXCPOT = MXBLLEN
      CALL MEMGET2('REAL','  DAOx',KDAOx ,LDAOx ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','  DAOy',KDAOy ,LDAOy ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','  DAOz',KDAOz ,LDAOz ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','  HESA',KHESA ,LHESA ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','   DST',KDST  ,LDST  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','   VFA',KVFA  ,LVFA  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL',' XCPOT',KXCPOT,LXCPOT,WORK,KFREE,LFREE)
CAMT
C
      CALL DFTINT_1(DMAT,NDMAT,NGEODRV,DOLND,CB,CBDATA,ELECTRONS,
     &            WORK(KGAO),WORK(KDMGAO),WORK(KDAOx),
     &            WORK(KDAOy),WORK(KDAOz),WORK(KHESA),WORK(KDST),
     &            WORK(KVFA),WORK(KXCPOT),WORK,KFREE,LFREE,IPRDFT)
C
      CALL QEXIT('DFTINT')
      END
      SUBROUTINE DFTINT_1(DMAT,NDMAT,NGEODRV,DOLND,CB,CBDATA,ELECTRONS,
     &                    GAO,DMAGAO,DAOx,DAOy,DAOz,HESA,DST,VFA,XCPOT,
     &                    WORK,KFREE,LFREE,IPRINT)

C
C     P. Salek and T. Helgaker, October 2003
C
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "inforb.h"
#include "shells.h"
#include "mxcent.h"
#include "dftacb.h"
#include "dftcom.h"
#include "dftbrhs.h"
#include "abainf.h"
#include "nuclei.h"
#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
#endif
C     choose reasonably large. Exceeding this limit means that boxes are
C     too large.
      PARAMETER(NBUFLEN=99000)
C
      DIMENSION DMAT(NBAST,NBAST,NDMAT), WORK(*), CBDATA(*)
      EXTERNAL CB
      LOGICAL DOLND
      DIMENSION DMAGAO(MXBLLEN,NBAST),GAO(MXBLLEN,NBAST,22)
c     for Hessians
      DIMENSION DAOx(MXBLLEN,NBAST),DAOy(MXBLLEN,NBAST)
      DIMENSION DAOz(MXBLLEN,NBAST)
      DIMENSION HESA(6,MXBLLEN)
C     For AC
      DIMENSION DST(NATOMS), VFA(MXBLLEN)
      DIMENSION XCPOT(MXBLLEN)
C
      DIMENSION COOR(3,NBUFLEN), WEIGHT(NBUFLEN)
      DIMENSION NSHLBLCK(2,KMAX), NBLCNT(8), NBLOCKS(2,NBAST,8)
      DIMENSION RHOA(MXBLLEN), GRADA(3,MXBLLEN)
      DIMENSION RHOB(MXBLLEN), GRADB(3,MXBLLEN)
C
      LOGICAL DOGGA, DFT_ISGGA
      EXTERNAL DFT_ISGGA
      DIMENSION IORIDX(KMAX,2,8)

#ifdef ENABLE_XCINT
#define EXPORT_GRID
#endif
#ifdef EXPORT_GRID
      integer              :: nr_points_export
      integer, parameter   :: max_export_grid_length = 1000000
      integer, parameter   :: lunit_export = 137
      real(8), allocatable :: rx(:)
      real(8), allocatable :: ry(:)
      real(8), allocatable :: rz(:)
      real(8), allocatable :: rw(:)
      logical              :: export_file_exists
      integer              :: irank, ierr, nr_proc
#ifdef VAR_MPI
      integer              :: status_array(MPI_STATUS_SIZE)
#endif
#endif /* EXPORT_GRID */

      ELECTRONS = 0.D0
      DOGGA = DFT_ISGGA() ! C code
      IGEODRV = NGEODRV
      IF (DOGGA) IGEODRV = IGEODRV + 1
CAMT Require higher derivatives for explicit potential construction
      IF (LDFTVXC) THEN
        IF (IGEODRV.LT.2) IGEODRV = 2
      ENDIF
      CALL SETUPSOS(IGEODRV,DOLND,IDUM1,IDUM2)
      CALL OPNQUA(NBAST,DMAT,WORK(KFREE),LFREE,IPRINT) ! C code
      NPOINTS = 0
      CALL CONSTRUCT_IORIDX(IORIDX)

#ifdef EXPORT_GRID
      irank   = 0
      nr_proc = 1
#ifdef VAR_MPI
      call mpi_comm_rank(mpi_comm_world, irank, ierr)
      call mpi_comm_size(mpi_comm_world, nr_proc, ierr)
#endif
      if (irank == 0) then
         inquire(file='num_grid', exist=export_file_exists)
      end if
#ifdef VAR_MPI
      if (nr_proc > 1) then
         call mpi_bcast(export_file_exists, 1, mpi_logical,
     &                  0, mpi_comm_world, ierr)
      end if
#endif
      if (.not. export_file_exists) then
         nr_points_export = 0
         allocate(rx(max_export_grid_length))
         allocate(ry(max_export_grid_length))
         allocate(rz(max_export_grid_length))
         allocate(rw(max_export_grid_length))
      end if
#endif /* EXPORT_GRID */

 100  CONTINUE
         CALL REAQUA(NSHELL,NSHLBLCK,NBUFLEN,COOR,WEIGHT,NLEN) ! C code
         if (nlen .le. 0) then
           go to 200
         end if
         NPOINTS = NPOINTS + NLEN

#ifdef EXPORT_GRID
         if (.not. export_file_exists) then
            if ((nr_points_export + nlen) > max_export_grid_length) then
               print *, 'ERROR in grid export, array too small'
               call quit('ERROR in grid export, array too small')
            end if
            do i = 1, nlen
               rx(nr_points_export + i) = coor(1, i)
               ry(nr_points_export + i) = coor(2, i)
               rz(nr_points_export + i) = coor(3, i)
               rw(nr_points_export + i) = weight(i)
               if (dabs(rw(nr_points_export + i)) < 1.0e-40) then
                  rw(nr_points_export + i) = 0.0d0
              end if
            enddo
            nr_points_export = nr_points_export + nlen
         end if
#endif /* EXPORT_GRID */

         DO IPT = 1, NLEN, MXBLLEN
            NCURLEN=MIN(MXBLLEN,NLEN-IPT+1)
            CALL BLGETSOS(NCURLEN,GAO,COOR(1,IPT),NSHELL,NSHLBLCK,
     &                    WORK(KFREE),LFREE,NBAST,DOLND,DOGGA,DFTHRI,0)
            CALL SHLTOORB(NSHELL,NSHLBLCK,NBLCNT,NBLOCKS,IORIDX)
            IF(DOGGA.OR.DFTASC) THEN
              IF (LDFTVXC.AND.(.NOT.DOVB)) THEN
               CALL GETRHO_BLOCKED_HES(DMAT,GAO,NBLCNT,NBLOCKS,
     &              NSHELL,DMAGAO,DAOx,DAOy,DAOz,NCURLEN,
     &              RHOA,GRADA,HESA)
               IF (NDMAT .EQ. 2) THEN
                 CALL QUIT('DFTAC / VXC NYI FOR NDMAT.GT.1')
               ENDIF
              ELSE
               CALL GETRHO_BLOCKED_GGA(DMAT,GAO,NBLCNT,NBLOCKS,NSHELL,
     &                                 DMAGAO,NCURLEN,RHOA,GRADA)
               IF (NDMAT .EQ. 2) CALL GETRHO_BLOCKED_GGA(DMAT(1,1,2),
     &                                 GAO,NBLCNT,NBLOCKS,NSHELL,
     &                                 DMAGAO,NCURLEN,RHOB,GRADB)
              ENDIF
            ELSE
               CALL GETRHO_BLOCKED_LDA(DMAT,GAO,NBLCNT,NBLOCKS,NSHELL,
     &                                 DMAGAO,NCURLEN,RHOA)
               IF (NDMAT .EQ. 2) CALL GETRHO_BLOCKED_LDA(DMAT(1,1,2),
     &                                 GAO,NBLCNT,NBLOCKS,NSHELL,
     &                                 DMAGAO,NCURLEN,RHOB)
            END IF


            IF (NDMAT .EQ. 2) THEN
              DO I = 1, NCURLEN
                 ELECTRONS = ELECTRONS + WEIGHT(IPT+I-1)
     &                       * (RHOA(I)+RHOB(I))
              END DO
            ELSE
              DO I = 1, NCURLEN
                 ELECTRONS = ELECTRONS + WEIGHT(IPT+I-1)*RHOA(I)
              END DO
            ENDIF

            IF ((LDFTVXC.AND.DOGGA).AND.(.NOT.DOVB)) THEN
              IF (NDMAT .EQ. 2) THEN
                 CALL QUIT('DFTAC / VXC NYI FOR NDMAT.GT.1')
              ELSE
                CALL CB(NCURLEN,NBLCNT,NBLOCKS,NSHELL,GAO,RHOA,
     &                  GRADA,HESA,DST,VFA,XCPOT,COOR(1,IPT),
     &                  WEIGHT(IPT),CBDATA)
              ENDIF
            ELSE
              IF (NDMAT .EQ. 2) THEN
                CALL CB(NCURLEN,NBLCNT,NBLOCKS,NSHELL,GAO,RHOA,RHOB,
     &                  GRADA,GRADB,COOR(1,IPT),WEIGHT(IPT),CBDATA)
              ELSE
                 CALL CB(NCURLEN,NBLCNT,NBLOCKS,NSHELL,GAO,RHOA,GRADA,
     &                   DST,VFA,XCPOT,COOR(1,IPT),WEIGHT(IPT),CBDATA)
              END IF
            ENDIF
         END DO
         GO TO 100
 200     CONTINUE

#ifdef EXPORT_GRID
      if (.not. export_file_exists) then
         if (irank == 0) then
            open(lunit_export,
     &           form   = 'formatted',
     &           status = 'new',
     &           access = 'sequential',
     &           file   = 'num_grid')
            rewind lunit_export
            write(lunit_export, *) nr_proc

            ! write points from master
            write(lunit_export, *) nr_points_export
            do i = 1, nr_points_export
               write(lunit_export, '(e20.12)') rx(i)
               write(lunit_export, '(e20.12)') ry(i)
               write(lunit_export, '(e20.12)') rz(i)
               write(lunit_export, '(e20.12)') rw(i)
            enddo
         end if

#ifdef VAR_MPI
         ! write points from workers
         if (nr_proc > 1) then
            do iworker = 1, nr_proc-1
               if (irank == iworker) then
                  call mpi_send(nr_points_export, 1, mpi_integer,
     &                          0, 92, mpi_comm_world, ierr)
                  call mpi_send(rx, nr_points_export, mpi_real8,
     &                          0, 93, mpi_comm_world, ierr)
                  call mpi_send(ry, nr_points_export, mpi_real8,
     &                          0, 94, mpi_comm_world, ierr)
                  call mpi_send(rz, nr_points_export, mpi_real8,
     &                          0, 95, mpi_comm_world, ierr)
                  call mpi_send(rw, nr_points_export, mpi_real8,
     &                          0, 96, mpi_comm_world, ierr)
               end if
               if (irank == 0) then
                  call mpi_recv(nr_points_export, 1, mpi_integer,
     &                          iworker, 92, mpi_comm_world,
     &                          status_array, ierr)
                  call mpi_recv(rx, nr_points_export, mpi_real8,
     &                          iworker, 93, mpi_comm_world,
     &                          status_array, ierr)
                  call mpi_recv(ry, nr_points_export, mpi_real8,
     &                          iworker, 94, mpi_comm_world,
     &                          status_array, ierr)
                  call mpi_recv(rz, nr_points_export, mpi_real8,
     &                          iworker, 95, mpi_comm_world,
     &                          status_array, ierr)
                  call mpi_recv(rw, nr_points_export, mpi_real8,
     &                          iworker, 96, mpi_comm_world,
     &                          status_array, ierr)

                  write(lunit_export, *) nr_points_export
                  do i = 1, nr_points_export
                     write(lunit_export, '(e20.12)') rx(i)
                     write(lunit_export, '(e20.12)') ry(i)
                     write(lunit_export, '(e20.12)') rz(i)
                     write(lunit_export, '(e20.12)') rw(i)
                  enddo
               end if
            end do
         end if
#endif

         if (irank == 0) then
            close(lunit_export, status = 'keep')
         end if

         deallocate(rx)
         deallocate(ry)
         deallocate(rz)
         deallocate(rw)
      end if
#endif /* EXPORT_GRID */

      CALL CLSQUA !C code
C
C     Test on the number of electrons
C
C     WRITE(2,*) ' node ', mynum,' knows of ', electrons
      CALL  DFTINTCOLLECT(ELECTRONS)
      IF (MYNUM.EQ.0) THEN
        ELCTRX = FLOAT(2*NISHT + NASHT)
        ERROR  = ELECTRONS - ELCTRX
        IF (ABS(ERROR) .GT. DFTELS*ELCTRX) THEN
          WRITE (LUPRI,'(4(/2X,A,F14.6),/2X,A)')
     &    ' Number of electrons from numerical integration:',ELECTRONS,
     &    ' Number of electrons from orbital occupations:  ',ELCTRX,
     &    ' Error in the number of electrons:              ',ERROR,
     &    ' Error larger than DFTELS (set input):          ',DFTELS,
     &    ' Calculation aborted.'
          CALL QUIT
     &    ('Wrong number of electrons in DFTINT. Calculation aborted.')
        END IF
      END IF
      RETURN
      END
      SUBROUTINE CONSTRUCT_IORIDX(IORIDX)
#include "implicit.h"
#include "maxorb.h"
#include "shells.h"
#include "inforb.h"
c KMAX in this context is an upper limit for the number of blocks.
      DIMENSION IORIDX(2,KMAX,NSYM)
c
c     ISHELL contains a shell index for given basis function.
c
      DIMENSION ISHELL(NBAST)
c
c     CONSTRUCT_IROIDX construct a list of indices containing
c     for a given symmetry ranges of basis sets that given shell
c     contributes to.
c     if given shell does not contribute to any orbital in given
c     symmetry, then the corresponding iordix elements are [0,-1]
c     and the shell should be skipped.
c
#include "pincom.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "symmet.h"
c

c       Just to find out errors in the code.
      ISTRA = 1
      DO IREPA = 0, MAXREP
         NORBA = NAOS(IREPA+1)
         DO I = ISTRA,ISTRA + NORBA - 1
            ISHEL = IAND(ISHFT(IPIND(I),-16),65535)
            ISHELL(I) = ISHEL
         END DO
         ISTRA = ISTRA + NORBA
      END DO
      DO ISYM = 1, NSYM
         DO I = 1, KMAX
            IORIDX(1,I,ISYM) =  0
            IORIDX(2,I,ISYM) = -1
         END DO
         DO I = IBAS(ISYM)+1, IBAS(ISYM) + NBAS(ISYM)
            ISHEL = ISHELL(I)
            IF(IORIDX(1,ISHEL,ISYM).LE.0) IORIDX(1,ISHEL,ISYM) = I
            IORIDX(2,ISHEL,ISYM) = I
         END DO
      END DO
      END
c     ===============================================================
c     MPI-related routines.
c     They are responsible for data distribution and collection.
c     ===============================================================
      SUBROUTINE DFTINTBCAST
#ifdef VAR_MPI
c     Equivalent of lifesupport from dft-qr: synchronizes all data needed
c     for evaluation of basis functions.  One probably wants to call it
c     as early as possible so that slaves can take right decisions
c     regarding which functionals to run etc.
c
c     NORBT, N2ORBX,NOCCT, NVIRT are needed for....?
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
#include "mxcent.h"
#include "mpif.h"
c
#include "infpar.h"
#include "inforb.h"
#include "lmns.h"
#include "primit.h"
#include "nuclei.h"
#include "onecom.h"
#include "pincom.h"
#include "shells.h"
#include "dftcom.h"
#include "sphtrm.h"
#include "symmet.h"
#include "xyzpow.h"
c     sync inforb:
!     write(lupri,*) ' arrived in DFTINTBCAST', mytid
      CALL MPI_Bcast(muld2h, 8*8,   my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(nbas,   8,     my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(nsym,   1,     my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NBAST, 1, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(N2BASX, 1, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
c     /lmns/
      CALL MPI_Bcast(lvalua, MXAQN, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(mvalua, MXAQN, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(nvalua, MXAQN, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
c     /ONECOM/
      CALL MPI_Bcast(jsta,   1,     my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(nuca,   1,     my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
c     /PINCOM/
      CALL MPI_Bcast(ipind,  MXCORB, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
c     /PRIMIT/
      NPRICCF = MXCONT*MXPRIM
      CALL MPI_Bcast(priccf, NPRICCF,MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(priexp, MXPRIM, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ierr)
c     /SHELLS/
      CALL MPI_Bcast(cent, MXSHEL*2*3, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(istbao, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(jstrt, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(kckt, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(khkt, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(kmax, 1, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(kstrt, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(nhkt, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(nuco, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(numcf, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(sphr, MXSHEL, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
c     /SPHTRM/
      CALL MPI_Bcast(csp, NCSP, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(ispadr, MXQN,       my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      NPTCNT = 3*MXCENT*8*2
      CALL MPI_Bcast(iptcnt, NPTCNT,     my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(iptsym, MXCORB*8,   my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
c     /SYMMET/
      CALL MPI_Bcast(isymao, MXAQN*MXQN, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(isymax, 3*2, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(maxopr, 1, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(maxrep, 1, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(naos, 8, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(pt, 8, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ie)
c     /XYZPOW/
      CALL MPI_Bcast(istep, MXAQNM, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(mval, MXAQNM, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(nval, MXAQNM, my_MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
c     /DFTCOM/
      CALL MPI_Bcast(DFTADD,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ie)

!     sync thresholds
      call mpi_bcast(dfthr0, 1, mpi_double_precision,
     &                       0, mpi_comm_world, ie)
      call mpi_bcast(dfthrl, 1, mpi_double_precision,
     &                       0, mpi_comm_world, ie)
      call mpi_bcast(dfthri, 1, mpi_double_precision,
     &                       0, mpi_comm_world, ie)

c     /DFTASC/
      CALL MPI_Bcast(DFTASC,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(DOMPOLE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(DOLB94,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(LDFTVXC,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ie)

c     sync functionals, too.
!     print *, 'calling DFTFUNCSYNC for mytid',mytid, NODTOT
      CALL DFTFUNCSYNC(MYNUM,NODTOT)
c
#endif
      END
c
      SUBROUTINE DFTINTCOLLECT(ELECTRONS)
#include "implicit.h"
#ifdef VAR_MPI
c
#include "mpif.h"
#include "maxorb.h"
#include "infpar.h"
      A = ELECTRONS
!     print *, 'A is for mytid',mytid
      CALL MPI_Reduce(A,ELECTRONS,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)
#endif
      END


      SUBROUTINE GETRHO_BLOCKED_HES(DMAT,GAO,NBLOCKS,IBLOCKS,
     *                              LDAIB,TMP,TMPx,TMPy,TMPz,NVCLEN,
     *                              RHO,GRAD,HES)
#include "implicit.h"
#include "inforb.h"
#include "dftinf.h"
      DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,NTYPSO)
      DIMENSION NBLOCKS(NSYM),IBLOCKS(2,LDAIB,NSYM)
      DIMENSION RHO(NVCLEN), GRAD(3,NVCLEN),TMP(NVCLEN,NBAST)
      DIMENSION TMPx(NVCLEN,NBAST),TMPy(NVCLEN,NBAST)
      DIMENSION TMPz(NVCLEN,NBAST),HES(6,NVCLEN)
#include "priunit.h"
c
      CALL zeroorbs(TMP,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)
      CALL zeroorbs(TMPx,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)
      CALL zeroorbs(TMPy,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)
      CALL zeroorbs(TMPz,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)

        DO ISYM = 1, NSYM
        DO IBL=1, NBLOCKS(ISYM)
          ISTART=IBLOCKS(1,IBL,ISYM)
          ILEN=IBLOCKS(2,IBL,ISYM)-ISTART+1
           DO JBL=1, NBLOCKS(ISYM)
           JSTART=IBLOCKS(1,JBL,ISYM)
           JLEN=IBLOCKS(2,JBL,ISYM)-JSTART+1
           call dgemm('N','N',NVCLEN,JLEN,ILEN,1.0d0,
     *            gao(1,ISTART,1),NVCLEN,
     *            DMAT(ISTART,JSTART),NBAST,1.0d0,
     *            TMP(1,JSTART),NVCLEN)

           call dgemm('N','N',NVCLEN,JLEN,ILEN,1.0d0,
     *            gao(1,ISTART,2),NVCLEN,
     *            DMAT(ISTART,JSTART),NBAST,1.0d0,
     *            TMPx(1,JSTART),NVCLEN)

           call dgemm('N','N',NVCLEN,JLEN,ILEN,1.0d0,
     *            gao(1,ISTART,3),NVCLEN,
     *            DMAT(ISTART,JSTART),NBAST,1.0d0,
     *            TMPy(1,JSTART),NVCLEN)

           call dgemm('N','N',NVCLEN,JLEN,ILEN,1.0d0,
     *            gao(1,ISTART,4),NVCLEN,
     *            DMAT(ISTART,JSTART),NBAST,1.0d0,
     *            TMPz(1,JSTART),NVCLEN)

      END DO
      END DO
      END DO
      call dzero(RHO, NVCLEN)
      call dzero(GRAD, 3*NVCLEN)
      call dzero(HES, 6*NVCLEN)
      DO ISYM = 1, NSYM
      DO IBL=1, NBLOCKS(ISYM)
      DO IDX=IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
         DO K = 1, NVCLEN
            RHO(K)    = RHO(K)    + GAO(K,IDX,1)*TMP(K,IDX)
            GRAD(1,K) = GRAD(1,K) + 2*GAO(K,IDX,2)*TMP(K,IDX)
            GRAD(2,K) = GRAD(2,K) + 2*GAO(K,IDX,3)*TMP(K,IDX)
            GRAD(3,K) = GRAD(3,K) + 2*GAO(K,IDX,4)*TMP(K,IDX)
            HES(1,K)  = HES(1,K)  + 2*GAO(K,IDX,5)*TMP(K,IDX)
     *                 + 2*GAO(K,IDX,2)*TMPx(K,IDX)
            HES(2,K)  = HES(2,K)  + 2*GAO(K,IDX,6)*TMP(K,IDX)
     *        + GAO(K,IDX,2)*TMPy(K,IDX)+GAO(K,IDX,3)*TMPx(K,IDX)
            HES(3,K)  = HES(3,K)  + 2*GAO(K,IDX,7)*TMP(K,IDX)
     *        + GAO(K,IDX,2)*TMPz(K,IDX)+GAO(K,IDX,4)*TMPx(K,IDX)
            HES(4,K)  = HES(4,K)  + 2*GAO(K,IDX,8)*TMP(K,IDX)
     *                 + 2*GAO(K,IDX,3)*TMPy(K,IDX)
            HES(5,K)  = HES(5,K)  + 2*GAO(K,IDX,9)*TMP(K,IDX)
     *        + GAO(K,IDX,3)*TMPz(K,IDX)+GAO(K,IDX,4)*TMPy(K,IDX)
            HES(6,K)  = HES(6,K)  + 2*GAO(K,IDX,10)*TMP(K,IDX)
     *                 + 2*GAO(K,IDX,4)*TMPz(K,IDX)

         END DO
      END DO
      END DO
      END DO
      END




