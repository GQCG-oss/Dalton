      SUBROUTINE DFTMOLGRAD(WORK,LWORK,IPRINT)
C
C     T. Helgaker and P. Salek, October 2003
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <inforb.h>
#include <priunit.h>
#include "dftmolgrad.h"
      DIMENSION WORK(LWORK)
#include <energy.h>
      EXTERNAL DFTMOLGRD 
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      DOGGA = DFT_ISGGA()
C      
      CALL DZERO(GRADFT,MXCOOR)
C
      KDMAT = 1
      KLAST = KDMAT + NBAST*NBAST
      LWRK  = LWORK - KLAST +1
      IF(KLAST.GT.LWORK) CALL QUIT('NOMEM IN DFTGRAD')
      CALL DFTDNS(WORK(KDMAT),WORK(KLAST),LWRK,0)
      CALL DFTINT(WORK(KDMAT),1,1,.FALSE.,WORK(KLAST),LWRK,
     &            DFTMOLGRD,WORK(KDMAT),ELE) 
      RETURN
      END
C
      SUBROUTINE DFTMOLGRD(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                     COORD,WGHT,DMAT)
C
C     T. Helgaker oct 2003
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
#include "dftmolgrad.h"
#include <inforb.h>
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          DMAT(NBAST,NBAST)
C
#include <dftinf.h>
#include <energy.h>
#include <symmet.h>
#include <shells.h>
c
      LOGICAL ACTIVE(0:NBAST)
      DIMENSION VXC(NBLEN), VXB(NBLEN), VX(5)
      DIMENSION TMP(NBLEN), KVALS(3,3)
      DATA ((KVALS(I,J), I = 1,3), J = 1,3) /1, 2, 3, 2, 4, 5, 3, 5, 6/
C
#include <ibtfun.h>
C 
      DO I = 0, NBAST
         ACTIVE(I) = .FALSE.
      END DO
      DO ISYM = 1, NSYM
         DO IBLA = 1, NBLCNT(ISYM)
         DO I = NBLOCKS(1,IBLA,ISYM), NBLOCKS(2,IBLA,ISYM) 
            ACTIVE(I) = .TRUE.
         END DO
         END DO
      END DO
C
C     Exchange-correlation contribution to molecular gradient 
C
      IF (DOGGA) THEN
         DO I = 1, NBLEN
            GRDNRM = SQRT(GRADA(1,I)**2 + GRADA(2,I)**2 + GRADA(3,I)**2)
            CALL DFTPTF0(RHOA(I),GRDNRM,WGHT(I),VX)
            VXC(I) = D2*VX(1)
            VXB(I) = D2*VX(2)/GRDNRM
         END DO
      ELSE
         DO I = 1, NBLEN
            CALL DFTPTF0(RHOA(I),D0,WGHT(I),VX)
            VXC(I) = D2*VX(1)
         END DO
      END IF
C
      DO IX = 1, 3
         DO ISYM = 1, NSYM
            IORBA = 0
            DO ISHELA = 1, KMAX
               ISCOOR = IPTCNT(3*(NCENT(ISHELA) - 1) + IX,0,1)
               DO ICOMPA = 1, KHKT(ISHELA)
                  IORBA = IORBA + 1
                  IA = IPTSYM(IORBA,ISYM-1)
                  KA = IPTSYM(IORBA,IBTXOR(ISYM-1,ISYMAX(IX,1)))
                  IF (ACTIVE(IA) .AND. KA.GT.0) THEN
                     IF (.NOT.DOGGA) THEN
                        DO I = 1, NBLEN
                           TMP(I) = D0
                        END DO
                        DO IBLB = 1, NBLCNT(ISYM)
                        DO IB=NBLOCKS(1,IBLB,ISYM),NBLOCKS(2,IBLB,ISYM) 
                           DO I = 1, NBLEN
                              TMP(I) = TMP(I) + GAO(I,IB,1)*DMAT(IB,IA)
                           END DO
                        END DO
                        END DO
                        FRC = D0
                        DO I = 1, NBLEN
                           FRC = FRC + VXC(I)*TMP(I)*GAO(I,KA,IX+1)
                        END DO
                     ELSE
                        FRC = D0
                        K1 = KVALS(1,IX) + 4
                        K2 = KVALS(2,IX) + 4
                        K3 = KVALS(3,IX) + 4
                        DO I = 1, NBLEN
                           GA  = GAO(I,KA,IX+1)
                           GAX = GRADA(1,I)*GA
                           GAY = GRADA(2,I)*GA
                           GAZ = GRADA(3,I)*GA
                           GA2 = GRADA(1,I)*GAO(I,KA,K1) 
     &                         + GRADA(2,I)*GAO(I,KA,K2)
     &                         + GRADA(3,I)*GAO(I,KA,K3)
                           GD = D0
                           GF = D0
                           DO IBLB = 1, NBLCNT(ISYM)
                              ISTRB = NBLOCKS(1,IBLB,ISYM)
                              IENDB = NBLOCKS(2,IBLB,ISYM) 
                              DO IB = ISTRB, IENDB
                                 GD = GD + DMAT(IB,IA)*GAO(I,IB,1)
                                 GF = GF + DMAT(IB,IA)*(GAO(I,IB,1)*GA2
     &                                                + GAO(I,IB,2)*GAX
     &                                                + GAO(I,IB,3)*GAY
     &                                                + GAO(I,IB,4)*GAZ)
                              END DO
                           END DO
                           FRC = FRC + VXC(I)*GD*GA + VXB(I)*GF
                        END DO
                     END IF
                     GRADFT(ISCOOR) = GRADFT(ISCOOR) - FRC 
                   END IF
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
