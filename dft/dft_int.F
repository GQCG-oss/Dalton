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
      SUBROUTINE DFTINT(DMAT,NDMAT,NGEODRV,DOLND,WORK,LWORK,CB,
     &                  CBDATA,ELECTRONS)
#include <implicit.h>      
#include <inforb.h>
#include <dftcom.h>
#include <priunit.h>
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
C
      CALL DFTINT_1(DMAT,NDMAT,NGEODRV,DOLND,CB,CBDATA,ELECTRONS,
     &            WORK(KGAO),WORK(KDMGAO),WORK,KFREE,LFREE)
C
      CALL QEXIT('DFTINT')
      END
      SUBROUTINE DFTINT_1(DMAT,NDMAT,NGEODRV,DOLND,CB,CBDATA,ELECTRONS,
     &                  GAO,DMAGAO,WORK,KFREE,LFREE)

C
C     P. Salek and T. Helgaker, October 2003
C
#include <implicit.h>      
#include <priunit.h>
C
#include <maxorb.h>
#include <inforb.h>
#include <shells.h>
#include <dftcom.h>
#ifdef VAR_MPI
#include <infpar.h>
#include <mpif.h>
#endif
C     choose reasonably large. Exceeding this limit means that boxes are
C     too large.
      PARAMETER(NBUFLEN=99000)
C
      DIMENSION DMAT(NBAST,NBAST,NDMAT), WORK(*), CBDATA(*)
      EXTERNAL CB
      LOGICAL DOLND
      DIMENSION DMAGAO(MXBLLEN,NBAST),GAO(MXBLLEN,NBAST,22)
C
      DIMENSION COOR(3,NBUFLEN), WEIGHT(NBUFLEN)
      DIMENSION NSHLBLCK(2,KMAX), NBLCNT(8), NBLOCKS(2,NBAST,8)
      DIMENSION RHOA(MXBLLEN), GRADA(3,MXBLLEN)
      DIMENSION RHOB(MXBLLEN), GRADB(3,MXBLLEN) 
C
      LOGICAL DOGGA, DFT_ISGGA
      EXTERNAL DFT_ISGGA
      DIMENSION IORIDX(KMAX,2,8)
C     
      ELECTRONS = 0.D0
      DOGGA = DFT_ISGGA() ! C code
      IGEODRV = NGEODRV
      IF (DOGGA) IGEODRV = IGEODRV + 1
      CALL SETUPSOS(IGEODRV,DOLND,IDUM1,IDUM2)
      CALL OPNQUA(NBAST,DMAT,WORK(KFREE),LFREE) ! C code
      NPOINTS = 0
      CALL CONSTRUCT_IORIDX(IORIDX)

!#define EXPORT_GRID
#ifdef EXPORT_GRID
      lunit_export = 137
      open(lunit_export,
     &     form   = 'formatted',
     &     status = 'unknown',
     &     access = 'sequential',
     &     file   =  'grid.export')

      rewind lunit_export
#endif /* ifdef EXPORT_GRID */

 100  CONTINUE
         CALL REAQUA(NSHELL,NSHLBLCK,NBUFLEN,COOR,WEIGHT,NLEN) ! C code
         if (nlen .le. 0) then
#ifdef EXPORT_GRID
!          there are no more points
!          write some negative integer
!          so that DIRAC knows that it's over
           write(lunit_export, *) -1
#endif /* ifdef EXPORT_GRID */
           go to 200
         end if
         NPOINTS = NPOINTS + NLEN

#ifdef EXPORT_GRID
         write(lunit_export, *) nlen
         do i = 1, nlen
           write(lunit_export, '(4d20.12)') coor(1, i),
     &                                      coor(2, i),
     &                                      coor(3, i),
     &                                      weight(i)
         enddo
#endif /* ifdef EXPORT_GRID */

         DO IPT = 1, NLEN, MXBLLEN
            NCURLEN=MIN(MXBLLEN,NLEN-IPT+1)
            CALL BLGETSOS(NCURLEN,GAO,COOR(1,IPT),NSHELL,NSHLBLCK,
     &                    WORK(KFREE),LFREE,NBAST,DOLND,DOGGA,DFTHRI,0)
            CALL SHLTOORB(NSHELL,NSHLBLCK,NBLCNT,NBLOCKS,IORIDX)
            IF(DOGGA) THEN
               CALL GETRHO_BLOCKED_GGA(DMAT,GAO,NBLCNT,NBLOCKS,NSHELL,
     &                                 DMAGAO,NCURLEN,RHOA,GRADA)
               IF (NDMAT .EQ. 2) CALL GETRHO_BLOCKED_GGA(DMAT(1,1,2),
     &                                 GAO,NBLCNT,NBLOCKS,NSHELL,
     &                                 DMAGAO,NCURLEN,RHOB,GRADB)
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
     &                     * (RHOA(I)+RHOB(I))
              END DO
              CALL CB(NCURLEN,NBLCNT,NBLOCKS,NSHELL,GAO,RHOA,RHOB, 
     &                GRADA,GRADB,COOR(1,IPT),WEIGHT(IPT),CBDATA)
            ELSE 
               DO I = 1, NCURLEN
                  ELECTRONS = ELECTRONS + WEIGHT(IPT+I-1)*RHOA(I)
               END DO
               CALL CB(NCURLEN,NBLCNT,NBLOCKS,NSHELL,GAO,RHOA,GRADA,
     &                 COOR(1,IPT),WEIGHT(IPT),CBDATA)
            END IF 
         END DO
         GO TO 100
 200     CONTINUE

#ifdef EXPORT_GRID
      close(lunit_export, status = 'keep')
#endif /* ifdef EXPORT_GRID */

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
#include <implicit.h>
#include <maxorb.h>
#include <shells.h>
#include <inforb.h>
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
#include <pincom.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <symmet.h>
c
#include <ibtfun.h>
c       Just to find out errors in the code.
      ISTRA = 1
      DO IREPA = 0, MAXREP
         NORBA = NAOS(IREPA+1)
         DO I = ISTRA,ISTRA + NORBA - 1
            ISHEL = IBTAND(IBTSHR(IPIND(I),16),65535)
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
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <maxorb.h>
#include <mxcent.h>
#include <mpif.h>
c
#include <infpar.h>
#include <inforb.h>
#include <lmns.h>
#include <primit.h>
#include <nuclei.h>
#include <onecom.h>
#include <pincom.h>
#include <shells.h>
#include <dftcom.h>
#include <sphtrm.h>
#include <symmet.h>
#include <xyzpow.h>
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

c     sync functionals, too.
!     print *, 'calling DFTFUNCSYNC for mytid',mytid, NODTOT
      CALL DFTFUNCSYNC(MYNUM,NODTOT)
c      
#endif
      END
c
      SUBROUTINE DFTINTCOLLECT(ELECTRONS)
#ifdef VAR_MPI
c
#include <implicit.h>
#include <mpif.h>
#include <maxorb.h>
#include <infpar.h>
      A = ELECTRONS
!     print *, 'A is for mytid',mytid
      CALL MPI_Reduce(A,ELECTRONS,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)
#endif
      END
