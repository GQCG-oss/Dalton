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

c /* Deck dftlnd */
      SUBROUTINE DFTBRHS(FD,WORK,LWORK,IPRINT)
C
C     T. Helgaker Jan 2004 
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, DP5 = 0.5D0)
#include <inforb.h>
#include "dftbrhs.h"
      DIMENSION FD(N2BASX,3), WORK(LWORK)
#include <dftcom.h>
      EXTERNAL DFTBRHS1
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      DOVB = DFT_ISGGA() .AND. .NOT.DFTPOT
C
      KDMAT = 1
      KFMAT = KDMAT +   N2BASX 
      KLAST = KFMAT + 3*N2BASX 
      LWRK  = LWORK - KLAST + 1
      IF(KLAST.GT.LWORK) CALL QUIT('NOMEM IN DFTBRHS')
      CALL DFTDNS(WORK(KDMAT),WORK(KLAST),LWRK,0)
      CALL DZERO(WORK(KFMAT),3*N2BASX)
      IF(DOVB) THEN
         NGEODRV = 1
      ELSE
         NGEODRV = 0
      END IF
      CALL DFTINT(WORK(KDMAT),1,NGEODRV,.TRUE.,WORK(KLAST),LWRK,
     &            DFTBRHS1,WORK(KFMAT),ELE)
      DO K = 1, 3
         DO I = 1, NBAST
         DO J = 1, I - 1
            IJADR = KFMAT + N2BASX*(K-1) + NBAST*(J-1) + I - 1
            JIADR = KFMAT + N2BASX*(K-1) + NBAST*(I-1) + J - 1
            AVERAG = DP5*(WORK(IJADR) - WORK(JIADR))
            WORK(IJADR) =   AVERAG
            WORK(JIADR) = - AVERAG
         END DO
         END DO
      END DO
      CALL DAXPY(3*N2BASX,D1,WORK(KFMAT),1,FD,1)
C
      IF (IPRINT.GT.10) THEN
         CALL HEADER('DMAT in DFTBRHS ',-1)
         CALL OUTPUT(WORK(KDMAT),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('x component of FMAT in DFTBRHS ',-1)
         CALL OUTPUT(WORK(KFMAT          ),1,NBAST,1,NBAST,NBAST,NBAST,
     &               1,LUPRI)
         CALL HEADER('y component of FMAT in DFTBRHS ',-1)
         CALL OUTPUT(WORK(KFMAT +  N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &               1,LUPRI)
         CALL HEADER('z component of FMAT in DFTBRHS ',-1)
         CALL OUTPUT(WORK(KFMAT +2*N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &               1,LUPRI)
      END IF
      RETURN
      END
C
      SUBROUTINE DFTBRHS1(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                    COORD,WGHT,FMAT)
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          FMAT(NBAST,NBAST,3)
#include <dftinf.h>
C
      CALL DFTBRHS2(NBLEN,NBLCNT,NBLOCKS,LDAIB,
     &              GAO(1,1,1),GAO(1,1,2),GAO(1,1,NSOB),GAO(1,1,NSOB1),
     &              RHOA,GRADA,COORD,WGHT,FMAT)
C
      END
      SUBROUTINE DFTBRHS2(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,GAO1,GAB1,GAB2,
     &                    RHOA,GRADA,COORD,WGHT,FMAT)
C
C     T. Helgaker sep 99/ oct 00 / feb 01 / jan 04
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C     differentiated with respect to magnetic field
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
C
#include "dftbrhs.h"
#include <inforb.h>
#include <nuclei.h>
#include <dfterg.h>
#include <energy.h>
#include <dftcom.h>
#include <dftinf.h>
#include <orgcom.h>
#include <symmet.h>
C
      DIMENSION COORD(3,NBLEN),WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          FMAT(NBAST,NBAST,3)
      DIMENSION GAO(NBLEN,NBAST), GAO1(NBLEN,NBAST,3), 
     &          GAB1(NBLEN,NBAST,3), GAB2(NBLEN,NBAST,3,3) 
      DIMENSION VXC(NBLEN), VXB(NBLEN), VX(5), GRAVB(3)
C
      IF (.NOT.DOVB) THEN
         DO I = 1, NBLEN
            CALL DFTPTF0(RHOA(I),D0,WGHT(I),VX)
            VXC(I) = D2*VX(1)
         END DO
         DO N = 1, NBLEN
            DO K = 1, 3
               IF (K.EQ.1) THEN
                  KY = 2
                  KZ = 3
               ELSE IF (K.EQ.2) THEN
                  KY = 3
                  KZ = 1
               ELSE
                  KY = 1
                  KZ = 2
               END IF 
               G10Z = VXC(N)*(COORD(KY,N) - ORIGIN(KY))
               G10Y = VXC(N)*(COORD(KZ,N) - ORIGIN(KZ))
               DO JSYM = 1, NSYM
                  ISYM = MULD2H(JSYM,ISYMAX(K,2) + 1) 
                  DO JBLB = 1, NBLCNT(JSYM)
                  DO JB = NBLOCKS(1,JBLB,JSYM), NBLOCKS(2,JBLB,JSYM)
                     FYZ = G10Y*GAB1(N,JB,KY) - G10Z*GAB1(N,JB,KZ)
                     DO IBLB = 1, NBLCNT(ISYM)
                     DO IB = NBLOCKS(1,IBLB,ISYM), NBLOCKS(2,IBLB,ISYM)
                        FMAT(IB,JB,K) = FMAT(IB,JB,K) + FYZ*GAO(N,IB)
                     END DO
                     END DO
                  END DO
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO I = 1, NBLEN
            GRDNRM = SQRT(GRADA(1,I)**2 + GRADA(2,I)**2 + GRADA(3,I)**2)
            CALL DFTPTF0(RHOA(I),GRDNRM,WGHT(I),VX)
            VXC(I) = D2*VX(1)
            VXB(I) = D2*VX(2)/GRDNRM
         END DO
         DO N = 1, NBLEN
            GRAVB(1) = VXB(N)*GRADA(1,N)
            GRAVB(2) = VXB(N)*GRADA(2,N)
            GRAVB(3) = VXB(N)*GRADA(3,N)
            DO K = 1, 3
               IF (K.EQ.1) THEN
                  KY = 2
                  KZ = 3
               ELSE IF (K.EQ.2) THEN
                  KY = 3
                  KZ = 1
               ELSE
                  KY = 1
                  KZ = 2
               END IF 
               PY = COORD(KY,N) - ORIGIN(KY)
               PZ = COORD(KZ,N) - ORIGIN(KZ)
               G10Z = VXC(N)*PY + GRAVB(KY)
               G10Y = VXC(N)*PZ + GRAVB(KZ)
               G21Z = PY*GRAVB(1)
               G21Y = PZ*GRAVB(1)
               G22Z = PY*GRAVB(2)
               G22Y = PZ*GRAVB(2)
               G23Z = PY*GRAVB(3)
               G23Y = PZ*GRAVB(3)
               DO JSYM = 1, NSYM
                  ISYM = MULD2H(JSYM,ISYMAX(K,2) + 1)
                  DO JBLB = 1, NBLCNT(JSYM)
                  DO JB = NBLOCKS(1,JBLB,JSYM), NBLOCKS(2,JBLB,JSYM)
                     FYZ = G10Y*GAB1(N,JB,  KY) - G10Z*GAB1(N,JB,  KZ)
     &                   + G21Y*GAB2(N,JB,1,KY) - G21Z*GAB2(N,JB,1,KZ)
     &                   + G22Y*GAB2(N,JB,2,KY) - G22Z*GAB2(N,JB,2,KZ)
     &                   + G23Y*GAB2(N,JB,3,KY) - G23Z*GAB2(N,JB,3,KZ)
                     GD = GRAVB(1)*GAO1(N,JB,1) + GRAVB(2)*GAO1(N,JB,2) 
     &                                          + GRAVB(3)*GAO1(N,JB,3)
                     FY = PY*GD
                     FZ = PZ*GD
                     DO IBLB = 1, NBLCNT(ISYM)
                     DO IB = NBLOCKS(1,IBLB,ISYM), NBLOCKS(2,IBLB,ISYM)
                        FMAT(IB,JB,K) = FMAT(IB,JB,K) + FYZ*GAO(N,IB)
     &                                                + FY*GAB1(N,IB,KZ) 
     &                                                - FZ*GAB1(N,IB,KY)
                     END DO
                     END DO
                  END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
      END
      SUBROUTINE DFTMAG(EXCMAT,COORX,COORY,COORZ,GAO,GAO1,GAB1,
     &                  GAB2,VXC,VXB,RH,DOGGA,FROMVX)
C
C     T. Helgaker sep 99/ oct 00 / feb 01
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C     differentiated with respect to magnetic field
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
      PARAMETER (D2 = 2.0D0)
C
#include <inforb.h>
#include <nuclei.h>
#include <dfterg.h>
#include <energy.h>
#include <dftcom.h>
#include <dftinf.h>
#include <orgcom.h>
#include <symmet.h>
C
      LOGICAL FROMVX, DOGGA
      DIMENSION GAO(NBAST), GAO1(NBAST,3), 
     &          GAB1(NBAST,3), GAB2(NBAST,3,3), 
     &          EXCMAT(NBAST,NBAST,3), RH(3)
C
      POX = COORX - ORIGIN(1)
      POY = COORY - ORIGIN(2)
      POZ = COORZ - ORIGIN(3)
C
      DO K = 1, 3
         KSYM = ISYMAX(K,2) + 1
         IF (K.EQ.1) THEN
            KY = 2
            KZ = 3
            PY = POY
            PZ = POZ
            RY = RH(2)
            RZ = RH(3)
         ELSE IF (K.EQ.2) THEN
            KY = 3
            KZ = 1
            PY = POZ
            PZ = POX
            RY = RH(3)
            RZ = RH(1)
         ELSE
            KY = 1
            KZ = 2
            PY = POX
            PZ = POY
            RY = RH(1)
            RZ = RH(2)
         END IF 
         IF (.NOT.(DOGGA .AND. .NOT.FROMVX)) THEN
            G10Z = D2*VXC*PY
            G10Y = D2*VXC*PZ
            DO JSYM = 1, NSYM
               ISYM = MULD2H(JSYM,KSYM) 
               JSTR = IBAS(JSYM) + 1
               JEND = IBAS(JSYM) + NBAS(JSYM)
               ISTR = IBAS(ISYM) + 1
               IEND = IBAS(ISYM) + NBAS(ISYM)
               DO J = JSTR, JEND
                  FYZ = G10Y*GAB1(J,KY) - G10Z*GAB1(J,KZ)
                  DO I = ISTR, IEND
                     EXCMAT(I,J,K) = EXCMAT(I,J,K) + FYZ*GAO(I)
                  END DO
               END DO
            END DO
         ELSE 
            G10Z = D2*VXC*PY + VXB*RY
            G10Y = D2*VXC*PZ + VXB*RZ
            G21Z = VXB*PY*RH(1)
            G21Y = VXB*PZ*RH(1)
            G22Z = VXB*PY*RH(2)
            G22Y = VXB*PZ*RH(2)
            G23Z = VXB*PY*RH(3)
            G23Y = VXB*PZ*RH(3)
            IF (NSYM.EQ.1) THEN
               DO J = 1, NBAST
                  FYZ = G10Y*GAB1(J,  KY) - G10Z*GAB1(J,  KZ)
     &                + G21Y*GAB2(J,1,KY) - G21Z*GAB2(J,1,KZ)
     &                + G22Y*GAB2(J,2,KY) - G22Z*GAB2(J,2,KZ)
     &                + G23Y*GAB2(J,3,KY) - G23Z*GAB2(J,3,KZ)
                  GD = VXB*(RH(1)*GAO1(J,1)+RH(2)*GAO1(J,2)
     $                 +RH(3)*GAO1(J,3))
                  FY = PY*GD
                  FZ = PZ*GD
                  DO I = 1, NBAST 
                     EXCMAT(I,J,K) = EXCMAT(I,J,K) + FYZ*GAO(I)
     &                             + FY*GAB1(I,KZ) - FZ*GAB1(I,KY)
                  END DO
               END DO
            ELSE
               DO JSYM = 1, NSYM
                  ISYM = MULD2H(JSYM,KSYM) 
                  JSTR = IBAS(JSYM) + 1
                  JEND = IBAS(JSYM) + NBAS(JSYM)
                  ISTR = IBAS(ISYM) + 1
                  IEND = IBAS(ISYM) + NBAS(ISYM)
                  DO J = JSTR, JEND
                     FYZ = G10Y*GAB1(J,  KY) - G10Z*GAB1(J,  KZ)
     &                   + G21Y*GAB2(J,1,KY) - G21Z*GAB2(J,1,KZ)
     &                   + G22Y*GAB2(J,2,KY) - G22Z*GAB2(J,2,KZ)
     &                   + G23Y*GAB2(J,3,KY) - G23Z*GAB2(J,3,KZ)
                     GD = VXB*(RH(1)*GAO1(J,1) + RH(2)*GAO1(J,2) 
     &                                       + RH(3)*GAO1(J,3))
                     FY = PY*GD
                     FZ = PZ*GD
                     DO I = ISTR, IEND 
                        EXCMAT(I,J,K) = EXCMAT(I,J,K) + FYZ*GAO(I)
     &                                + FY*GAB1(I,KZ) - FZ*GAB1(I,KY)
                     END DO
                  END DO
               END DO
            END IF
         END IF
      END DO
      RETURN
      END
