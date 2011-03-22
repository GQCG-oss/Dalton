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
      SUBROUTINE DFTKSMb(DMAT,FKSM,EDFTY,WORK,LWORK,IPRFCK)
C
C     P. Salek and T. Helgaker October 2003
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <infpar.h>
#include <inforb.h>
      DIMENSION DMAT(NBAST,NBAST),FKSM(NBAST,NBAST),WORK(LWORK)
      EXTERNAL DFTKSMGGA,DFTKSMLDA
      LOGICAL DOGGA, DFT_ISGGA
      EXTERNAL DFT_ISGGA
#ifdef MOD_XCFUN
      external xcfun_KSMLDA
      external xcfun_KSMGGA
#endif
      COMMON /DKSMPRIV/ ENERGY_DFTKSM
C      
      CALL QENTER('DFTKSMb')
C
      IF (NODTOT.gt.0) CALL KICK_KSM_SLAVES_ALIVE(NBAST,DMAT,IPRFCK)
      DOGGA = DFT_ISGGA()
      ENERGY_DFTKSM = 0.D0
      IF(DOGGA) THEN
#ifdef MOD_XCFUN
         CALL DFTINT(DMAT,1,0,.FALSE.,WORK,LWORK,
     &               xcfun_KSMGGA,FKSM,ELE)
#else
         CALL DFTINT(DMAT,1,0,.FALSE.,WORK,LWORK,
     &               DFTKSMGGA,FKSM,ELE)
#endif
      ELSE
         DO I = 1, NBAST
            FKSM(I,I) = FKSM(I,I)*2.0D0
         END DO
#ifdef MOD_XCFUN
         CALL DFTINT(DMAT,1,0,.FALSE.,WORK,LWORK,
     &               xcfun_KSMLDA,FKSM,ELE)
#else
         CALL DFTINT(DMAT,1,0,.FALSE.,WORK,LWORK,
     &               DFTKSMLDA,FKSM,ELE)
#endif
C
      END IF
      IF (NODTOT.gt.0) CALL KSMCOLLECT(FKSM,ENERGY_DFTKSM,WORK,LWORK)
      DO I = 1, NBAST
         DO J = 1, I-1
            AVERAG = 0.5*(FKSM(J,I) + FKSM(I,J))
            FKSM(J,I) = AVERAG
            FKSM(I,J) = AVERAG
         END DO
         IF(.NOT.DOGGA) FKSM(I,I) = FKSM(I,I)*0.5D0
      END DO
      EDFTY = ENERGY_DFTKSM
      ELE_ERROR = ELE - 2.0D0*NRHFT
      IF (ABS(ELE_ERROR) .gt. 1.0D-8 .or. IPRFCK .ge. 0)
     &   WRITE(LUPRI,'(T7,A,F20.12,F15.10,1P,D12.2)')
     &     'K-S energy, electrons, error :', EDFTY, ELE, ELE_ERROR
      CALL QEXIT('DFTKSMb')
      RETURN
      END
C
      SUBROUTINE DFTKSMLDA(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                     RHOA,GRADA,COORD,WGHT,FKSM)
C
C     P. Salek and T. Helgaker oct 2003
C
#include <implicit.h>
      PARAMETER (D2 = 2.0D0)
#include <inforb.h>
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &     RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &     FKSM(*)
      EXTERNAL DFTENE
      COMMON /DKSMPRIV/ ENERGY_DFTKSM
C
      DIMENSION VXC(NBLEN),VX(5)
      PARAMETER (DUMMY = 0D0)
c     CALL QENTER('DFTKSMLDA')
c     .. called often - thus QENTER is too much overhead
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      DO IPNT = 1, NBLEN
         CALL DFTPTF0(RHOA(IPNT),DUMMY,WGHT(IPNT),VX)
         ENERGY_DFTKSM = ENERGY_DFTKSM
     &        + DFTENE(RHOA(IPNT),DUMMY)*WGHT(IPNT)
         VXC(IPNT) = 2D0*VX(1)
      END DO
      CALL DISTLDAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN,
     &     VXC,GAO,FKSM)
c     CALL QEXIT('DFTKSMLDA')
      RETURN
      END
c
      SUBROUTINE DFTKSMGGA(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                     RHOA,GRADA,COORD,WGHT,FKSM)
C
C     P. Salek and T. Helgaker oct 2003
C
#include <implicit.h>
#include <inforb.h>
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &     RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &     FKSM(*)
      EXTERNAL DFTENE
C
      DIMENSION VXC(2,NBLEN),VX(5)
      COMMON /DKSMPRIV/ ENERGY_DFTKSM
c     CALL QENTER('DFTKSMGGA')
c     .. called often - thus QENTER is too much overhead
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      DO IPNT = 1, NBLEN
         GRD = SQRT(GRADA(1,IPNT)**2+GRADA(2,IPNT)**2+GRADA(3,IPNT)**2)
         CALL DFTPTF0(RHOA(IPNT),GRD,WGHT(IPNT),VX)
         ENERGY_DFTKSM = ENERGY_DFTKSM
     &        + DFTENE(RHOA(IPNT),GRD)*WGHT(IPNT)
         VXC(1,IPNT) = VX(1)
         IF(GRD.GT.1D-40) THEN
            VXC(2,IPNT) = VX(2)/(0.5*GRD)
         ELSE
            VXC(2,IPNT) = 0D0
         END IF
      END DO
      CALL DISTGGAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN,
     &              VXC,GAO,GRADA,FKSM)
c     CALL QEXIT('DFTKSMGGA')
C
      RETURN
      END
c     ================================================================
c     DISTLDAI AND DISTGGAI are matrix element distribution routines
c     using inline versions of the loops. This make sense for very fast
c     compilers and slow BLAS routines.
c     
c     DISTRLDAI - helper subroutine to distribute given Omega with
c     given set of rho-dependent coefficients.
c
      SUBROUTINE DISTLDAI(NBLEN,NBLCNT,NBLOCKS,LDAIB,IBLSTART,IBLEND,
     &                    COEF,GAOS,EXCMAT)
c
c     NBLEN - number of grid points in the batch.
c     NBLCNT - number of "active" blocks in each symmetry.
c     NBLOCKS(:,) - start and stop indexes for each blocks in each symmetry.
c
#include <implicit.h>
#include <inforb.h>
c
      DIMENSION NBLCNT(NSYM),NBLOCKS(2,LDAIB,NSYM)
      DIMENSION GAOS(NBLEN,NBAST), COEF(NBLEN)
      DIMENSION EXCMAT(NBAST,NBAST)
c
      DIMENSION TMP(NBLEN)
c
      DO ISYM = 1, NSYM         
         DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO K = IBLSTART, IBLEND
                  TMP(K) = GAOS(K,J)*coef(K)
               END DO
               DO IBL = 1, NBLCNT(ISYM)
                  ITOP = MIN(J,NBLOCKS(2,IBL,ISYM))
                  DO I = NBLOCKS(1,IBL,ISYM),ITOP
                     DO K = IBLSTART, IBLEND
                        EXCMAT(I,J) = EXCMAT(I,J) + GAOS(k,I)*TMP(k)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
c     
c     DISTGGAI - helper subroutine to distribute given Omega with
c     given set of rho and grad-dependent coefficients.
c
      SUBROUTINE DISTGGAI(NBLEN,NBLCNT,NBLOCKS,LDAIB,IBLSTART,IBLEND,
     &                    COEF,GAOS,GRAD,EXCMAT)
c
c     NBLEN - number of grid points in the batch.
c     NBLCNT - number of "active" blocks in each symmetry.
c     NBLOCKS(:,) - start and stop indexes for each blocks in each symmetry.
c
#include <implicit.h>
#include <inforb.h>
c
      DIMENSION NBLCNT(NSYM),NBLOCKS(2,LDAIB,NSYM)
      DIMENSION GAOS(NBLEN,NBAST,4), COEF(2,NBLEN),GRAD(3,NBLEN)
      DIMENSION EXCMAT(NBAST,NBAST)
c
      DIMENSION TMP(NBLEN)
c
      DO ISYM = 1, NSYM
         DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO K = IBLSTART, IBLEND
                  TMP(K) =
     &                   coef(1,K)* GAOS(K,J,1)
     &                 + coef(2,K)*(GAOS(K,J,2)*GRAD(1,K)+
     &                              GAOS(K,J,3)*GRAD(2,K)+
     &                              GAOS(K,J,4)*GRAD(3,K))
               END DO
               DO IBL = 1, NBLCNT(ISYM)
                  DO I = NBLOCKS(1,IBL,ISYM),NBLOCKS(2,IBL,ISYM)
                     DO K = IBLSTART, IBLEND
                        EXCMAT(I,J) = EXCMAT(I,J) + GAOS(K,I,1)*TMP(K)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
c     ================================================================
c     DISTLDAB AND DISTGGAB are matrix element distribution routines
c     using BLAS routines for the kernel operations. This make sense for
c     hardware-optimized BLAS routines only that are able to compensate
c     for the call overhead.
c     
      SUBROUTINE DISTLDAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,IBLSTART,IBLEND,
     &                    COEF,GAOS,EXCMAT)
c
c     NBLEN - number of grid points in the batch.
c     NBLCNT - number of "active" blocks in each symmetry.
c     NBLOCKS(:,) - start and stop indexes for each blocks in each symmetry.
c
#include <implicit.h>
#include <inforb.h>
c
      DIMENSION NBLCNT(NSYM),NBLOCKS(2,LDAIB,NSYM)
      DIMENSION GAOS(NBLEN,NBAST), COEF(NBLEN)
      DIMENSION EXCMAT(NBAST,NBAST)
c
      DIMENSION TMP(NBLEN,NBAST)
c
      NK=IBLEND-IBLSTART+1
      DO ISYM = 1, NSYM         
         DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO K = IBLSTART, IBLEND
                  TMP(K,J) = GAOS(K,J)*coef(K)
               END DO
            END DO
         END DO
         DO JBL = 1, NBLCNT(ISYM)
            JSTART=NBLOCKS(1,JBL,ISYM)
            JLEN=NBLOCKS(2,JBL,ISYM)-JSTART+1
            DO IBL = 1, JBL-1
               ISTART=NBLOCKS(1,IBL,ISYM)
               ILEN=NBLOCKS(2,IBL,ISYM)-ISTART+1
               call dgemm('T','N',ILEN,JLEN,NK,1.0d0,
     &              GAOS(IBLSTART,ISTART),NBLEN,
     &              TMP(IBLSTART,JSTART),NBLEN,1.0d0,
     &              EXCMAT(ISTART,JSTART),NBAST)
            END DO
c           the diagonal block which always exists (functions always
c           overlap with themselves without diagonal)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO I = NBLOCKS(1,JBL,ISYM), J
                  DO K = IBLSTART, IBLEND
                     EXCMAT(I,J) = EXCMAT(I,J) + GAOS(k,I)*TMP(k,J)
                  END DO
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
c
c     DISTGGAB - helper subroutine to distribute given Omega with
c     given set of rho and grad-dependent coefficients.
c
      SUBROUTINE DISTGGAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,IBLSTART,IBLEND,
     &                    COEF,GAOS,GRAD,EXCMAT)
c
c     NBLEN - number of grid points in the batch.
c     NBLCNT - number of "active" blocks in each symmetry.
c     NBLOCKS(:,) - start and stop indexes for each blocks in each symmetry.
c
#include <implicit.h>
#include <inforb.h>
c
      DIMENSION NBLCNT(NSYM),NBLOCKS(2,LDAIB,NSYM)
      DIMENSION GAOS(NBLEN,NBAST,4), COEF(2,NBLEN),GRAD(3,NBLEN)
      DIMENSION EXCMAT(NBAST,NBAST)
c
      DIMENSION TMP(NBLEN,NBAST)
c
      NK=IBLEND-IBLSTART+1
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
         DO JBL = 1, NBLCNT(ISYM)
            JSTART=NBLOCKS(1,JBL,ISYM)
            JLEN=NBLOCKS(2,JBL,ISYM)-JSTART+1
            DO IBL = 1, NBLCNT(ISYM)
               ISTART=NBLOCKS(1,IBL,ISYM)
               ILEN=NBLOCKS(2,IBL,ISYM)-ISTART+1
               call dgemm('T','N',ILEN,JLEN,NK,1.0d0,
     &              GAOS(IBLSTART,ISTART,1),NBLEN,
     &              TMP(IBLSTART,JSTART),NBLEN,1.0d0,
     &              EXCMAT(ISTART,JSTART),NBAST)
            END DO
         END DO
      END DO
      RETURN
      END
c     ================================================================
c     Parallelization-related part of the exchange-correlation
c     evaluation module.
c
      SUBROUTINE KICK_KSM_SLAVES_ALIVE(NBAST,DMAT,IPRINT)
c     do nothing when running serial code.
#include <implicit.h>
      DIMENSION DMAT(NBAST,NBAST)
C
#ifdef VAR_MPI
#include <maxorb.h>
#include <infpar.h>
#include <mpif.h>
      IF(MYNUM.EQ.MASTER) THEN
c        IPRTYP=4 is registered for KSM Slave driver.
         IPRTYP = 4
         CALL MPI_Bcast(iprtyp,1, my_MPI_INTEGER,
     &                  MASTER, MPI_COMM_WORLD,IERR)
         CALL MPI_Bcast(iprint,1, my_MPI_INTEGER,
     &                  MASTER, MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL KSMSYNC(DMAT)
      END IF
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
      CALL QENTER('DFT_KSMSLAVE')
c
      CALL DFTINTBCAST
      KDMAT = 1
      KFKSM = KDMAT + N2BASX
      KFREE = KFKSM + N2BASX
      IF(KFREE.GT.LWRK) CALL STOPIT('DFT_KSMSLAVE',' ',KFREE,LWRK)
      LFREE = LWRK - KFREE + 1
      CALL KSMSYNC(WRK(KDMAT))
      CALL DZERO(WRK(KFKSM),N2BASX)
      CALL DFTKSMb(WRK(KDMAT),WRK(KFKSM),EDFTY,
     &             WRK(KFREE),LFREE,IPRINT)
      CALL QEXIT('DFT_KSMSLAVE')
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
      SUBROUTINE KSMCOLLECT(FMAT,ENERGY,WRK,LFREE)
c
#include <implicit.h>
      DIMENSION FMAT(*), WRK(LFREE)
#ifdef VAR_MPI
c
#include <mpif.h>
#include <maxorb.h>
#include <infpar.h>
#include <inforb.h>
      IF(LFREE.LE.N2BASX) CALL STOPIT('KSMCOLLECT',' ',LFREE,N2BASX)
      ETMP = ENERGY
      CALL MPI_Reduce(ETMP,ENERGY,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)
      CALL DCOPY(N2BASX,FMAT,1,WRK,1)
      CALL MPI_Reduce(WRK,FMAT,N2BASX,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)
      
#endif
      END
