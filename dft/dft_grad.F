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
      CALL DFTINT(WORK(KDMAT),1,.TRUE.,.FALSE.,WORK(KLAST),LWRK,
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
C
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          DMAT(NBAST,NBAST)
C
#include <inforb.h>
#include <dftinf.h>
#include <energy.h>
#include <symmet.h>
#include <shells.h>
c
      LOGICAL ACTIVE(NBAST)
      DIMENSION VXC(NBLEN), VXB(NBLEN), VX(5)
      DIMENSION TMP(NBLEN)
C
#include <ibtfun.h>
C 
C     DO ISYM = 1, NSYM
C     DO IBLB = 1, NBLCNT(ISYM)
C     istr = NBLOCKS(1,IBLB,ISYM)
C     iend =  NBLOCKS(2,IBLB,ISYM) 
C     write (2,*) ' isym, iblb,istr,iend ',isym,iblb,istr,iend
C     write (2,*) ' 0 der ints '
C     call output(gao(1,1,1),1,nblen,istr,iend,nblen,nbast,1,2)
C     write (2,*) ' 1 der ints '
C     call output(gao(1,1,2),1,nblen,istr,iend,nblen,nbast,1,2)
C     call output(gao(1,1,3),1,nblen,istr,iend,nblen,nbast,1,2)
C     call output(gao(1,1,4),1,nblen,istr,iend,nblen,nbast,1,2)
C     write (2,*) ' 2 der ints '
C     call output(gao(1,1,5),1,nblen,istr,iend,nblen,nbast,1,2)
C     call output(gao(1,1,6),1,nblen,istr,iend,nblen,nbast,1,2)
C     call output(gao(1,1,7),1,nblen,istr,iend,nblen,nbast,1,2)
C     call output(gao(1,1,8),1,nblen,istr,iend,nblen,nbast,1,2)
C     call output(gao(1,1,9),1,nblen,istr,iend,nblen,nbast,1,2)
C     call output(gao(1,1,10),1,nblen,istr,iend,nblen,nbast,1,2)
C     end do
C     end do
      DO I = 1, NBAST
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
         IF (IX.EQ.1) THEN
            K1 = 1
            K2 = 2
            K3 = 3
         ELSE IF (IX.EQ.2) THEN
            K1 = 2
            K2 = 4
            K3 = 5
         ELSE
            K1 = 3
            K2 = 5
            K3 = 6
         END IF
C
         DO ISYM = 1, NSYM
            IRPAX = IBTXOR(ISYM-1,ISYMAX(IX,1))
            IORBA = 0
            DO ISHELA = 1, KMAX
               ISCOOR = IPTCNT(3*(NCENT(ISHELA) - 1) + IX,0,1)
               DO ICOMPA = 1, KHKT(ISHELA)
                  IORBA = IORBA + 1
                  IA = IPTSYM(IORBA,ISYM-1)
                  KA = IPTSYM(IORBA,IRPAX)
                  IF (KA.GT.0 .AND. ACTIVE(IA)) THEN
C
                  IF (.NOT.DOGGA) THEN
                     DO I = 1, NBLEN
                        TMP(I) = D0
                     END DO
                     DO IBLB = 1, NBLCNT(ISYM)
                     DO IB = NBLOCKS(1,IBLB,ISYM), NBLOCKS(2,IBLB,ISYM) 
                        DO I = 1, NBLEN
                           TMP(I) = TMP(I) + GAO(I,IB,1)*DMAT(IB,IA)
                        END DO
                     END DO
                     END DO
                     FRC = D0
                     DO I = 1, NBLEN
                        FRC = FRC + VXC(I)*TMP(I)*GAO(I,KA,IX+1)
                     END DO
                     GRADFT(ISCOOR) = GRADFT(ISCOOR) - FRC 
                  ELSE
                     FRC = D0
                     DO I = 1, NBLEN
                        GA  = GAO(I,KA,IX+1)
                        GAX = GRADA(1,I)*GA
                        GAY = GRADA(2,I)*GA
                        GAZ = GRADA(3,I)*GA
                        GA2 = GRADA(1,I)*GAO(I,KA,K1+4) 
     &                      + GRADA(2,I)*GAO(I,KA,K2+4)
     &                      + GRADA(3,I)*GAO(I,KA,K3+4)
                        GD = D0
                        GF = D0
                        DO IBLB = 1, NBLCNT(ISYM)
                           ISTRB = NBLOCKS(1,IBLB,ISYM)
                           IENDB = NBLOCKS(2,IBLB,ISYM) 
                           DO IB = ISTRB, IENDB
                              GD = GD + DMAT(IB,IA)*GAO(I,IB,1)
                              GF = GF + DMAT(IB,IA)*(GAO(I,IB,1)*GA2
     &                                             + GAO(I,IB,2)*GAX
     &                                             + GAO(I,IB,3)*GAY
     &                                             + GAO(I,IB,4)*GAZ)
                           END DO
                        END DO
                        FRC = FRC + VXC(I)*GD*GA + VXB(I)*GF
                     END DO
                     GRADFT(ISCOOR) = GRADFT(ISCOOR) - FRC 
                  END IF
C                    FRC = D0
C                    DO I = 1, NBLEN
C                       GD = D0 
C                       DO IBLB = 1, NBLCNT(ISYM)
C                          ISTRB = NBLOCKS(1,IBLB,ISYM)
C                          IENDB = NBLOCKS(2,IBLB,ISYM) 
C                          DO IB = ISTRB, IENDB
C                             GD = GD + GAO(I,IB,1)*DMAT(IB,IA)
C                          END DO
C                       END DO
C                       FRC = FRC + VXC(I)*GD*GAO(I,KA,IX+1)
C                    END DO
C                    GRADFT(ISCOOR) = GRADFT(ISCOOR) - FRC 
                  END IF
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
