c     DFTINT: perform numerical integration.
c     DMAT - NDMAT density matrices.
c     CB - routine to be called for each batch of points.
c     CBDATA - data to be passed to the CB routine as the extra argument.
c     
c     MPI model:
c     The integrator itself is responsible for distributing data
c     related to orbital and functional evaluation and collecting 
c     integrated charge. Nothing else.
c
c     All the MPI-related routines are at the bottom of this file.
c
      SUBROUTINE DFTINT(DMAT,NDMAT,NGEODRV,DOLND,WRK,LWRK,CB,
     &                  CBDATA,ELECTRONS)
#include <implicit.h>      
#include <inforb.h>
#include <dftcom.h>
#include <priunit.h>
      DIMENSION DMAT(NBAST,NBAST,NDMAT), WRK(LWRK), CBDATA(*)
      EXTERNAL CB
      LOGICAL DOLND
C
      KGAO   = 1
      KDMGAO = KGAO   + MXBLLEN*NBAST*22 
      KFREE  = KDMGAO + MXBLLEN*NBAST
      IF (KFREE.GT.LWRK) CALL STOPIT("DFTINT",'DFTINT',KFREE,LWRK)
      LFREE = LWRK - KFREE + 1
C
      CALL DFTIN1(DMAT,NDMAT,NGEODRV,DOLND,WRK(KFREE),LFREE,CB,
     &            CBDATA,ELECTRONS,WRK(KGAO),WRK(KDMGAO))
C
      END
      SUBROUTINE DFTIN1(DMAT,NDMAT,NGEODRV,DOLND,WORK,LWORK,CB,
     &                  CBDATA,ELECTRONS,GAO,DMAGAO)

C
C     P. Salek and T. Helgaker, October 2003
C
#include <implicit.h>      
#include <inforb.h>
#include <maxorb.h>
#include <shells.h>
#include <dftinf.h>
#include <dftcom.h>
#include <priunit.h>
#ifdef VAR_MPI
#include <infpar.h>
#include <mpif.h>
#endif
C     choose reasonably large. Exceeding this limit means that boxes are
C     too large.
      PARAMETER(NBUFLEN=50000)
C
      DIMENSION DMAT(NBAST,NBAST,NDMAT), WORK(LWORK), CBDATA(*)
      EXTERNAL CB
      LOGICAL DOLND
      DIMENSION DMAGAO(MXBLLEN,NBAST),GAO(MXBLLEN,NBAST,22)
C
      DIMENSION COOR(3,NBUFLEN), WEIGHT(NBUFLEN)
      DIMENSION NSHLBLCK(2,KMAX), NBLCNT(8), NBLOCKS(2,NBAST,8)
      DIMENSION RHOA(MXBLLEN), GRADA(3,MXBLLEN)
C
      LOGICAL DOGGA, DFT_ISGGA, CHECKELS
      EXTERNAL DFT_ISGGA
      DIMENSION IORIDX(KMAX,2,8)
C     
      ELECTRONS = 0.D0
      DOGGA = DFT_ISGGA() ! C code
      IGEODRV = NGEODRV
      IF (DOGGA) IGEODRV = IGEODRV + 1
      CALL SETUPSOS(IGEODRV,DOLND,NTYPSO,IDUMMY)
      CALL OPNQUA(NBAST,DMAT,WORK,LWORK) ! C code
      NPOINTS = 0
      CALL CONSTRUCT_IORIDX(IORIDX)
 100  CONTINUE
         CALL REAQUA(NSHELL,NSHLBLCK,NBUFLEN,COOR,WEIGHT,NLEN) ! C code
         IF(NLEN .LE. 0) GO TO 200
         NPOINTS = NPOINTS + NLEN
         DO IPT = 1, NLEN, MXBLLEN
            NCURLEN=MIN(MXBLLEN,NLEN-IPT+1)
            CALL BLGETSOS(NCURLEN,GAO,COOR(1,IPT),NSHELL,NSHLBLCK,
     &                    WORK,LWORK,NBAST,DOLND,DOGGA,DFTHRI,0)
            CALL SHLTOORB(NSHELL,NSHLBLCK,NBLCNT,NBLOCKS,IORIDX)
            IF(DOGGA) THEN
               CALL GETRHO_BLOCKED_GGA(DMAT,GAO,NBLCNT,NBLOCKS,NSHELL,
     &                                 DMAGAO,NCURLEN,RHOA,GRADA)
            ELSE
               CALL GETRHO_BLOCKED_LDA(DMAT,GAO,NBLCNT,NBLOCKS,NSHELL,
     &                                 DMAGAO,NCURLEN,RHOA)
            END IF
            DO I = 1, NCURLEN
               ELECTRONS = ELECTRONS + WEIGHT(IPT+I-1)*RHOA(I)
            END DO
            CALL CB(NCURLEN,NBLCNT,NBLOCKS,NSHELL,GAO,RHOA,GRADA,
     &              COOR(1,IPT),WEIGHT(IPT),CBDATA)
         END DO
         GO TO 100
 200     CONTINUE
      CALL CLSQUA !C code
C
C     Test on the number of electrons
C   
      CALL  DFTINTCOLLECT(ELECTRONS)
#ifdef VAR_MPI
      CHECKELS = MYNUM.EQ.0
#else
      CHECKELS = .TRUE.
#endif
      ELCTRX = FLOAT(2*NISHT + NASHT)
      ERROR  = ELECTRONS - ELCTRX
      IF (CHECKELS.AND.ABS(ERROR) .GT. DFTELS) THEN
         print *,"ERROR=",ERROR,ELECTRONS
         WRITE (LUPRI,'(4(/2X,A,F14.6),/2X,A)')
     &   ' Number of electrons from numerical integration:',ELECTRONS,
     &   ' Number of electrons from orbial occupations:   ',ELCTRX,
     &   ' Error in the number of electrons:              ',ERROR,
     &   ' Error larger than DFTELS (set input):          ',DFTELS,
     &   ' Calculation aborted.'
         CALL QUIT
     &    ('Wrong number of electrons in DFTINT. Calculation aborted.')
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
#ifdef VAR_MPI
      SUBROUTINE DFTINTBCAST
c     Equivalent of lifesupport from dft-qr: synchronizes all data needed
c     for evaluation of basis functions.  One probably wants to call it
c     as early as possible so that slaves can take right decisions
c     regarding which functionals to run etc.
c     
c     NORBT, N2ORBX,NOCCT, NVIRT are needed for....?
#include <implicit.h>
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
#include <sphtrm.h>
#include <symmet.h>
#include <xyzpow.h>
c     sync inforb:
      CALL MPI_Bcast(muld2h, 8*8,   MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(nbas,   8,     MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(nsym,   1,     MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NBAST, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(N2BASX, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
c     /lmns/
      CALL MPI_Bcast(lvalua, MXAQN, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(mvalua, MXAQN, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(nvalua, MXAQN, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
c     /ONECOM/
      CALL MPI_Bcast(jsta,   1,     MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(nuca,   1,     MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
c     /PINCOM/
      CALL MPI_Bcast(ipind,  MXCORB, MPI_INTEGER,
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
      CALL MPI_Bcast(istbao, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(jstrt, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(kckt, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(khkt, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(kmax, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(kstrt, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(nhkt, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(nuco, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(numcf, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(sphr, MXSHEL, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
c     /SPHTRM/
      CALL MPI_Bcast(csp, NCSP, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(ispadr, MXQN,       MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      NPTCNT = 3*MXCENT*8*2
      CALL MPI_Bcast(iptcnt, NPTCNT,     MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(iptsym, MXCORB*8,   MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
c     /SYMMET/
      CALL MPI_Bcast(isymao, MXAQN*MXQN, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(isymax, 3*2, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(maxopr, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(maxrep, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(naos, 8, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(pt, 8, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ie)
c     /XYZPOW/
      CALL MPI_Bcast(istep, MXAQNM, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(mval, MXAQNM, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
      CALL MPI_Bcast(nval, MXAQNM, MPI_INTEGER,
     &               0, MPI_COMM_WORLD,ie)
c
c     sync functionals, too.
      CALL DFTFUNCSYNC(MYNUM,NODTOT)
c      
      END
#endif
c
      SUBROUTINE DFTINTCOLLECT(ELECTRONS)
#ifdef VAR_MPI
c
#include <implicit.h>
#include <mpif.h>
#include <maxorb.h>
#include <infpar.h>
      A = ELECTRONS
      CALL MPI_Reduce(A,ELECTRONS,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)
#endif
      END
