c     DFTKSM - compute exchange-correlation contribution to the KS
c     matrix.
      SUBROUTINE DFTKSMb(DMAT,FKSM,EDFTY,WORK,LWORK,IPRFCK)
C
C     P. Salek and T. Helgaker October 2003
C
#include <implicit.h>
#include <inforb.h>
#include <priunit.h>
      DIMENSION DMAT(NBAST,NBAST),FKSM(NBAST,NBAST),WORK(LWORK)
      EXTERNAL DFTKSMGGA,DFTKSMLDA
      LOGICAL DOGGA, DFT_ISGGA
      EXTERNAL DFT_ISGGA
C      
      KENER = 1
      KKSM  = KENER + 1
      KLST  = KKSM  + NBAST*NBAST
      LWRK  = LWORK - KLST +1
      IF(KLST.GT.LWORK) CALL QUIT('NOMEM IN DFTKSM')
      CALL KICK_KSM_SLAVES_ALIVE(NBAST,DMAT,IPRFCK)
      DOGGA = DFT_ISGGA()
      CALL DZERO(WORK(KENER),NBAST*NBAST+1)
      IF(DOGGA) THEN
         CALL DFTINT(DMAT,1,0,.FALSE.,WORK(KLST),LWRK,
     &               DFTKSMGGA,WORK(KENER),ELE)
      ELSE
         CALL DFTINT(DMAT,1,0,.FALSE.,WORK(KLST),LWRK,
     &               DFTKSMLDA,WORK(KENER),ELE)
C
      END IF
      CALL KSMCOLLECT(WORK(KENER),WORK(KLST),LWORK)
      DO I = 1, NBAST
         DO J = 1, NBAST
            AVERAG = 0.5*(WORK(KKSM+(J-1)+(I-1)*NBAST)
     &                  + WORK(KKSM+(I-1)+(J-1)*NBAST))
            WORK(KKSM+(J-1)+(I-1)*NBAST) = AVERAG
            WORK(KKSM+(I-1)+(J-1)*NBAST) = AVERAG
         END DO
      END DO
      EDFTY = WORK(KENER)
#if 0
      DO I =1, NBAST
         WRITE(LUPRI,*) (WORK(KKSM+(I-1)+(J-1)*NBAST),J=1,NBAST)
      END DO
#endif
      CALL DAXPY(N2BASX,1D0,WORK(KKSM),1, FKSM,1)
#if 0
      DO I =1, NBAST
         WRITE(LUPRI,*) (FKSM(I,J),J=1,NBAST)
      END DO
#endif
      IF(IPRFCK.GE.0) WRITE(LUPRI,'(A,2F20.14)')
     &     'Kohn-Sham electrons/energy :', ELE, WORK(KENER)
      RETURN
      END
C
      SUBROUTINE DFTKSMLDA(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                     RHOA,GRADA,COORD,WGHT,CBDATA)
C
C     P. Salek and T. Helgaker oct 2003
C
#include <implicit.h>
      PARAMETER (D2 = 2.0D0)
#include <inforb.h>
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &     RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &     CBDATA(*)
      EXTERNAL DFTENE
C
#include <dftinf.h>
c
      DIMENSION VXC(NBLEN),VX(5)
      PARAMETER (DUMMY = 0D0)
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      KENERG = 1
      KKOHN  = 2
      DO IPNT = 1, NBLEN
         CALL DFTPTF0(RHOA(IPNT),DUMMY,WGHT(IPNT),VX)
         CBDATA(KENERG) = CBDATA(KENERG)
     &        + DFTENE(RHOA(IPNT),DUMMY)*WGHT(IPNT)
         VXC(IPNT) = 2D0*VX(1)
      END DO
      CALL DISTRLDA(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN,
     &     VXC,GAO,CbDATA(KKOHN))
      RETURN
      END
c
      SUBROUTINE DFTKSMGGA(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                     RHOA,GRADA,COORD,WGHT,CBDATA)
C
C     P. Salek and T. Helgaker oct 2003
C
#include <implicit.h>
#include <inforb.h>
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &     RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &     CBDATA(*)
      EXTERNAL DFTENE
C
#include <dftinf.h>
c
      DIMENSION VXC(2,NBLEN),VX(5)
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      KENERG = 1
      KKOHN  = 2
      DO IPNT = 1, NBLEN
         GRD = SQRT(GRADA(1,IPNT)**2+GRADA(2,IPNT)**2+GRADA(3,IPNT)**2)
         CALL DFTPTF0(RHOA(IPNT),GRD,WGHT(IPNT),VX)
         CBDATA(KENERG) = CBDATA(KENERG)
     &        + DFTENE(RHOA(IPNT),GRD)*WGHT(IPNT)
         VXC(1,IPNT) = VX(1)
         IF(GRD.GT.1D-40) THEN
            VXC(2,IPNT) = VX(2)/(0.5*GRD)
         ELSE
            VXC(2,IPNT) = 0D0
         END IF
      END DO
      CALL DISTRGGA(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN,
     &              VXC,GAO,GRADA,CbDATA(KKOHN))
      RETURN
      END
c     
c     DISTRLDA - helper subroutine to distribute given Omega with
c     given set of rho-dependent coefficients.
c
      SUBROUTINE DISTRLDA(NBLEN,NBLCNT,NBLOCKS,LDAIB,IBLSTART,IBLEND,
     &                    COEF,GAOS,EXCMAT)
#include <implicit.h>
c
c     NBLEN - number of grid points in the batch.
c     NBLCNT - number of "active" blocks in each symmetry.
c     NBLOCKS(:,) - start and stop indexes for each blocks in each symmetry.
c
#include <inforb.h>
c
      DIMENSION NBLCNT(NSYM),NBLOCKS(2,LDAIB,NSYM)
      DIMENSION GAOS(NBLEN,NBAST), COEF(NBLEN)
      DIMENSION EXCMAT(NBAST,NBAST)
c
      DIMENSION TMP(NBLEN,NBAST)
c
      DO ISYM = 1, NSYM
         DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO K = IBLSTART, IBLEND
                  TMP(K,J) = GAOS(K,J)*coef(K)
               END DO
            END DO
         END DO
         
         DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO IBL = 1, NBLCNT(ISYM)
                  ITOP = MIN(J-1,NBLOCKS(2,IBL,ISYM))
                  DO I = NBLOCKS(1,IBL,ISYM),ITOP
                     DO K = IBLSTART, IBLEND
                        EXCMAT(I,J) = EXCMAT(I,J) + GAOS(k,I)*TMP(k,J)
                     END DO
                  END DO
               END DO
               DO K = IBLSTART, IBLEND
                  EXCMAT(J,J) = EXCMAT(J,J) + 0.5*GAOS(k,J)*TMP(k,J)
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
c     
c     DISTRGGA - helper subroutine to distribute given Omega with
c     given set of rho and grad-dependent coefficients.
c
      SUBROUTINE DISTRGGA(NBLEN,NBLCNT,NBLOCKS,LDAIB,IBLSTART,IBLEND,
     &                    COEF,GAOS,GRAD,EXCMAT)
#include <implicit.h>
c
c     NBLEN - number of grid points in the batch.
c     NBLCNT - number of "active" blocks in each symmetry.
c     NBLOCKS(:,) - start and stop indexes for each blocks in each symmetry.
c
#include <inforb.h>
c
      DIMENSION NBLCNT(NSYM),NBLOCKS(2,LDAIB,NSYM)
      DIMENSION GAOS(NBLEN,NBAST,4), COEF(2,NBLEN),GRAD(3,NBLEN)
      DIMENSION EXCMAT(NBAST,NBAST)
c
      DIMENSION TMP(NBLEN,NBAST)
c
      DO ISYM = 1, NSYM
         DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO K = IBLSTART, IBLEND
                  TMP(K,J) =
     &                   coef(1,K)* GAOS(K,J,1)
     &                 + coef(2,K)*(GAOS(K,J,2)*GRAD(1,K)+
     &                              GAOS(K,J,3)*GRAD(2,K)+
     &                              GAOS(K,J,4)*GRAD(3,K))
               END DO
            END DO
         END DO
c         
         DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO IBL = 1, NBLCNT(ISYM)
                  DO I = NBLOCKS(1,IBL,ISYM),NBLOCKS(2,IBL,ISYM)
                     DO K = IBLSTART, IBLEND
                        EXCMAT(I,J) = EXCMAT(I,J) + GAOS(K,I,1)*TMP(K,J)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
c
      SUBROUTINE KICK_KSM_SLAVES_ALIVE(NBAST,DMAT,IPRINT)
#include <implicit.h>
c     do nothing when running serial code.
#ifdef VAR_MPI
#include <maxorb.h>
#include <infpar.h>
#include <mpif.h>
      DIMENSION DMAT(NBAST,NBAST)
      IF(MYNUM.EQ.MASTER) THEN
c        IPRTYP=4 is registered for KSM Slave driver.
         IPRTYP = 4
      WRITE(40+MYNUM,*)
     &     "master blows the horn!",mynum,master
         CALL MPI_Bcast(iprtyp,1, MPI_INTEGER,
     &                  MASTER, MPI_COMM_WORLD,IERR)
         CALL MPI_Bcast(iprint,1, MPI_INTEGER,
     &                  MASTER, MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL KSMSYNC(DMAT)
      END IF
      WRITE(40+MYNUM,*)
     &     "Slave checked status!",mynum,master
#endif
      END
c
#ifdef VAR_MPI
      SUBROUTINE DFT_KSMSLAVE(WRK,LWRK,IPRINT)
#include <implicit.h>
#include <maxorb.h>
#include <infpar.h>
#include <inforb.h>
      DIMENSION WRK(LWRK)
c
      CALL DFTINTBCAST
      KDMAT = 1
      KFKSM = KDMAT + N2BASX
      KFREE = KFKSM + N2BASX
      IF(KFREE.GT.LWRK) CALL STOPIT('DFT_KSMSLAVE','',KFREE,LWRK)
      LFREE = LWRK - KFREE + 1
      CALL KSMSYNC(WRK(KDMAT))
      CALL DZERO(WRK(KFKSM),N2BASX)
      CALL DFTKSMb(WRK(KDMAT),WRK(KFKSM),EDFTY,
     &             WRK(KFREE),LFREE,IPRINT)
      END
c      
      SUBROUTINE KSMSYNC(DMAT)
#include <implicit.h>
#include <mpif.h>
#include <priunit.h>
#include <inforb.h>
      CALL MPI_Bcast(DMAT,NBAST*NBAST,MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ierr)
      END
#endif
      SUBROUTINE KSMCOLLECT(FMAT,WRK,LFREE)
#ifdef VAR_MPI
c
#include <implicit.h>
#include <mpif.h>
#include <maxorb.h>
#include <infpar.h>
#include <inforb.h>
      DIMENSION FMAT(N2BASX)
      DIMENSION WRK(LFREE)
      IF(LFREE.LE.N2BASX) CALL STOPIT('KSMCOLLECT','',LFREE,N2BASX)
      CALL DCOPY(N2BASX+1,FMAT,1,WRK,1)
      CALL MPI_Reduce(WRK,FMAT,N2BASX+1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)
#endif
      END
