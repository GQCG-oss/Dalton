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
      SUBROUTINE DFTKSMb(DMAT,FKSM,EDFTY,WORK,LWORK,IPRFCK)
C
C     P. Salek and T. Helgaker October 2003
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
#include "dftcom.h"
#include "mxcent.h"
#include "dftacb.h"
      DIMENSION DMAT(NBAST,NBAST),FKSM(NBAST,NBAST),WORK(LWORK)
      EXTERNAL DFTKSMGGA,DFTKSMLDA, DFTKSMHES
      LOGICAL DOGGA, DFT_ISGGA
      EXTERNAL DFT_ISGGA
      COMMON /DKSMPRIV/ ENERGY_DFTKSM
C
      CALL QENTER('DFTKSMb')
C
      IF (NODTOT.gt.0) CALL KICK_KSM_SLAVES_ALIVE(NBAST,DMAT,IPRFCK)
      DOGGA = DFT_ISGGA()
      ENERGY_DFTKSM = 0.D0
      IF(DOGGA) THEN
CAMT If explicit v_xc use DFTKSMHES
        IF (LDFTVXC) THEN
         CALL DFTINT(DMAT,1,1,.FALSE.,WORK,LWORK,
     &               DFTKSMHES,FKSM,ELE)
        ELSE
         CALL DFTINT(DMAT,1,0,.FALSE.,WORK,LWORK,
     &               DFTKSMGGA,FKSM,ELE)
        ENDIF
      ELSE
         DO I = 1, NBAST
            FKSM(I,I) = FKSM(I,I)*2.0D0
         END DO
         CALL DFTINT(DMAT,1,0,.FALSE.,WORK,LWORK,
     &               DFTKSMLDA,FKSM,ELE)
C
      END IF
      IF (NODTOT.gt.0) CALL KSMCOLLECT(FKSM,ENERGY_DFTKSM,WORK,LWORK)
      DO I = 1, NBAST
         DO J = 1, I-1
            AVERAG = 0.5d0*(FKSM(J,I) + FKSM(I,J))
            FKSM(J,I) = AVERAG
            FKSM(I,J) = AVERAG
         END DO
         IF(.NOT.DOGGA) FKSM(I,I) = FKSM(I,I)*0.5D0
      END DO
      EDFTY = ENERGY_DFTKSM
      ELE_ERROR = ELE - 2.0D0*NRHFT
      IF (ABS(ELE_ERROR) .gt. 1.0D-8 .or. IPRFCK .ge. 0) THEN
        IF ((.NOT.SLAVE) .OR. (NODTOT.LE.20)) THEN
          WRITE(LUPRI,'(T7,A,F20.12,F15.10,1P,D12.2)')
     &     'K-S energy, electrons, error :', EDFTY, ELE, ELE_ERROR
        END IF
      END IF
      CALL QEXIT('DFTKSMb')
      RETURN
      END
C
      SUBROUTINE DFTKSMLDA(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                     RHOA,GRADA,DST,VFA,XCPOT,
     &                     COORD,WGHT,FKSM)
C
C     P. Salek and T. Helgaker oct 2003
C
#include "implicit.h"
      PARAMETER (D2 = 2.0D0)
#include "inforb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dftcom.h"
#include "dftacb.h"
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &     RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &     FKSM(*), DST(NATOMS), VFA(NBLEN), XCPOT(NBLEN)
      EXTERNAL DFTENE
      COMMON /DKSMPRIV/ ENERGY_DFTKSM
C
      DIMENSION VXC(NBLEN),VX(5)
      PARAMETER (DUMMY = 0.D0)
c     CALL QENTER('DFTKSMLDA')
c     .. called often - thus QENTER is too much overhead
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      DO IPNT = 1, NBLEN
         CALL DFTPTF0(RHOA(IPNT),DUMMY,WGHT(IPNT),VX)

         IF (DFTASC) THEN
           IF (LGRAC) GRDAC = 0.5d0*DSQRT(GRADA(1,IPNT)*GRADA(1,IPNT)
     &                      + GRADA(2,IPNT)*GRADA(2,IPNT)
     &                      + GRADA(3,IPNT)*GRADA(3,IPNT))
           IF (LGRAC.OR.DOLB94) THEN
                IF (RHOA(IPNT).GT.1.0d-13) THEN
                  RHO43 = (0.5d0*RHOA(IPNT))**(4.0d0/3.0d0)
                ELSE
                  RHO43 = 1.0d-13
                ENDIF
           ENDIF
           IF (DOLB94) THEN
             IF (RHOA(IPNT).GT.1.0d-13) THEN
               RHO13 = (0.5d0*RHOA(IPNT))**(1.0d0/3.0d0)
             ELSE
               RHO13 = 1.0d-13
             ENDIF
           ENDIF
           IF (WGHT(IPNT).NE.0.0d0) VA = VX(1) / WGHT(IPNT)
           IF (WGHT(IPNT).EQ.0.0d0) VA = 0.0d0
           CALL DFTAC(VA,VFA(IPNT),DST,COORD(1,IPNT),COORD(2,IPNT),
     &                COORD(3,IPNT),GRDAC,RHO43,RHO13)
           VX(1) = VA * WGHT(IPNT)
         ENDIF

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
     &                     RHOA,GRADA,DST,VFA,XCPOT,
     &                     COORD,WGHT,FKSM)
C
C     P. Salek and T. Helgaker oct 2003
C
#include "implicit.h"
#include "inforb.h"
#include "mxcent.h"
#include "nuclei.h"
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &     RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &     FKSM(*), DST(NATOMS), VFA(NBLEN), XCPOT(NBLEN)
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
            VXC(2,IPNT) = 2.0d0*VX(2)/GRD
         ELSE
            VXC(2,IPNT) = 0.D0
         END IF
      END DO
      CALL DISTGGAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN,
     &              VXC,GAO,GRADA,FKSM)
c     CALL QEXIT('DFTKSMGGA')
C
      RETURN
      END

      SUBROUTINE DFTKSMHES(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &            RHOA,GRADA,HESA,DST,VFA,XCPOT,
     &            COORD,WGHT,FKSM)
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "dftcom.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dftacb.h"
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &     RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &     FKSM(*),HESA(6,NBLEN),DST(NATOMS),VFA(NBLEN),XCPOT(NBLEN)
      EXTERNAL DFTENE
      COMMON /DKSMPRIV/ ENERGY_DFTKSM
#include "dftinf.h"
c
      DIMENSION VXC(NBLEN),VX(5)
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      KENERG = 1
      KKOHN  = 2
      DO IPNT = 1, NBLEN

         GRD = DSQRT(GRADA(1,IPNT)**2+GRADA(2,IPNT)**2+GRADA(3,IPNT)**2)
         IF (RHOA(IPNT).GT.1.0d-13) THEN
           call dftptfh(RHOA(IPNT),GRADA(1,IPNT),HESA(1,IPNT),
     &                  WGHT(IPNT),VX)
         ELSE
           VX(1)=0.0d0
           VX(2)=0.0d0
         ENDIF

         IF (DFTASC) THEN
           IF (LGRAC) GRDAC = 0.5d0*DSQRT(GRADA(1,IPNT)*GRADA(1,IPNT)
     &                      + GRADA(2,IPNT)*GRADA(2,IPNT)
     &                      + GRADA(3,IPNT)*GRADA(3,IPNT))
           IF (LGRAC.OR.DOLB94) THEN
                IF (RHOA(IPNT).GT.1.0d-13) THEN
                  RHO43 = (0.5d0*RHOA(IPNT))**(4.0d0/3.0d0)
                ELSE
                  RHO43 = 1.0d-13
                ENDIF
           ENDIF

           IF (DOLB94) THEN
             IF (RHOA(IPNT).GT.1.0d-13) THEN
               RHO13 = (0.5d0*RHOA(IPNT))**(1.0d0/3.0d0)
             ELSE
               RHO13 = 1.0d-13
             ENDIF
           ENDIF

           IF (WGHT(IPNT).NE.0.0d0) VA = VX(1) / WGHT(IPNT)
           IF (WGHT(IPNT).EQ.0.0d0) VA = 0.0d0
           CALL DFTAC(VA,VFA(IPNT),DST,COORD(1,IPNT),COORD(2,IPNT),
     &                COORD(3,IPNT),GRDAC,RHO43,RHO13)
           VX(1) = VA * WGHT(IPNT)
         ENDIF

         ENERGY_DFTKSM = ENERGY_DFTKSM
     &        + DFTENE(RHOA(IPNT),GRD)*WGHT(IPNT)
C
         VXC(IPNT) = VX(1)
         XCPOT(IPNT) = VX(1) / WGHT(IPNT)
      END DO
      CALL DISTRHES(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN,
     &              VXC,GAO,FKSM)
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
#include "implicit.h"
#include "inforb.h"
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
#include "implicit.h"
#include "inforb.h"
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
#include "implicit.h"
#include "inforb.h"
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
#include "implicit.h"
#include "inforb.h"
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

c
c     DISTRHES - helper subroutine to distribute given Omega with
c     given set of rho-dependent coefficients from calculations with
C     Hessian
c
      SUBROUTINE DISTRHES(NBLEN,NBLCNT,NBLOCKS,LDAIB,
     *                    IBLSTART,IBLEND,
     &                    COEF,GAOS,EXCMAT)
#include "implicit.h"
c
c     NBLEN - number of grid points in the batch.
c     NBLCNT - number of "active" blocks in each symmetry.
c     NBLOCKS(:,) - start and stop indexes for each blocks in each
c     symmetry.
c
#include "inforb.h"
c
      DIMENSION NBLCNT(NSYM),NBLOCKS(2,LDAIB,NSYM)
      DIMENSION GAOS(NBLEN,NBAST), COEF(NBLEN)
      DIMENSION EXCMAT(NBAST,NBAST)
c
      DIMENSION TMP(NBLEN,NBAST)

      INTEGER NK

      NK=IBLEND-IBLSTART+1
      if (nk .gt. NBLEN) call quit("NBLEN < NK in DISTRHES")
      DO ISYM = 1, NSYM
         DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO K = IBLSTART, IBLEND
                  TMP(K,J) =
     &                   coef(K)* GAOS(K,J)
               END DO
            END DO
         END DO
         DO IBL = 1, NBLCNT(ISYM)
           ISTART=NBLOCKS(1,IBL,ISYM)
           ILEN=NBLOCKS(2,IBL,ISYM)-ISTART+1
           DO JBL = 1, NBLCNT(ISYM)
            JSTART=NBLOCKS(1,JBL,ISYM)
            JLEN=NBLOCKS(2,JBL,ISYM)-JSTART+1
             call dgemm('T','N',ILEN,JLEN,NK,1.0d0,
     *                   GAOS(IBLSTART,ISTART),NBLEN,
     *                   TMP(IBLSTART,JSTART),NBLEN,1.0d0,
     *                   EXCMAT(ISTART,JSTART),NBAST)
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
#include "implicit.h"
      DIMENSION DMAT(NBAST,NBAST)
C
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
C defined parallel calculation types
#include "iprtyp.h"
      IF(MYNUM.EQ.MASTER) THEN
c        DFT_KSM_WORK is registered for KSM Slave driver.
         IPRTYP = DFT_KSM_WORK
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
#include "implicit.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
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
#include "implicit.h"
#include "mpif.h"
#include "priunit.h"
#include "inforb.h"
      REAL*8 DMAT(NBAST,NBAST)
      CALL MPI_Bcast(DMAT,NBAST*NBAST,MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD,ierr)
      END
#endif
      SUBROUTINE KSMCOLLECT(FMAT,ENERGY,WRK,LFREE)
c
#include "implicit.h"
      DIMENSION FMAT(*), WRK(LFREE)
#ifdef VAR_MPI
c
#include "mpif.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
      IF(LFREE.LE.N2BASX) CALL STOPIT('KSMCOLLECT',' ',LFREE,N2BASX)
      ETMP = ENERGY
      CALL MPI_Reduce(ETMP,ENERGY,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)
      CALL DCOPY(N2BASX,FMAT,1,WRK,1)
      CALL MPI_Reduce(WRK,FMAT,N2BASX,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,IERR)

#endif
      END
