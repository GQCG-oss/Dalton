c     DFTINT: perform numerical integration.
c     DMAT - NDMAT density matrices.
c     CB - routine to be called for each batch of points.
c     CBDATA - data to be passed to the CB routine as the extra argument.
c     
      SUBROUTINE DFTINT(DMAT,NDMAT,NEEDLAP,DOLND,WORK,LWORK,CB,CBDATA,
     &                  ELECTRONS)
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
C     choose reasonably large - so that loop unrolling gives speedup 
C     but small compared with a cache size max block length MXBLLEN.
      PARAMETER(MXBLLEN = 100)
C     choose reasonably large.
      PARAMETER(NBUFLEN=50000)
C
      DIMENSION DMAT(NBAST,NBAST,NDMAT), WORK(LWORK), CBDATA(*)
      EXTERNAL CB
      LOGICAL NEEDLAP, DOLND
C
      DIMENSION DMAGAO(MXBLLEN,NBAST),GAO(MXBLLEN,NBAST,10)
      DIMENSION COOR(3,NBUFLEN), WEIGHT(NBUFLEN)
      DIMENSION NSHLBLCK(2,KMAX),NBLCNT(8), NBLOCKS(2,NBAST,8)
      DIMENSION RHOA(MXBLLEN), GRADA(3,MXBLLEN)
C
      LOGICAL DOGGA, DFT_ISGGA
      EXTERNAL DFT_ISGGA
      DIMENSION IORIDX(KMAX,2,8)
C     
      ELECTRONS = 0.D0
      DOGGA = DFT_ISGGA() ! C code
      IGEODRV = 0
      IF (DOGGA) IGEODRV = IGEODRV + 1
      IF (NEEDLAP) IGEODRV = IGEODRV + 1
      CALL SETUPSOS(IGEODRV,DOLND,NTYPSO,IDUMMY)
      CALL OPNQUA(NBAST,WORK,LWORK) ! C code
      NPOINTS = 0
      CALL CONSTRUCT_IORIDX(IORIDX)
 100  CONTINUE
         CALL REAQUA(NSHELL,NSHLBLCK,NBUFLEN,COOR,WEIGHT,NLEN) ! C code
         IF(NLEN .LE. 0) GO TO 200
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
      ELCTRX = FLOAT(2*NISHT + NASHT)
      ERROR  = ELECTRONS - ELCTRX
      IF (ABS(ERROR) .GT. DFTELS) THEN
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
