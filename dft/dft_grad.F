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
      IF (NASHT.GT.0) THEN 
          CALL DFTMOLGRADAB(WORK,LWORK,1)
          RETURN 
      END IF    
C
      KDMAT = 1
      KLAST = KDMAT + NBAST*NBAST
      LWRK  = LWORK - KLAST +1
      IF(KLAST.GT.LWORK) CALL QUIT('NOMEM IN DFTGRAD')
      CALL DFTDNS(WORK(KDMAT),WORK(KLAST),LWRK,0)
      CALL KICK_SLAVES_GRAD(NBAST,WORK(KDMAT),IPRINT)
      CALL DFTINT(WORK(KDMAT),1,1,.FALSE.,WORK(KLAST),LWRK,
     &            DFTMOLGRD,WORK(KDMAT),ELE) 
      CALL GRADSLAVE_COLLECT(GRADFT,WORK(KLAST),LWRK)
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
C
      SUBROUTINE KICK_SLAVES_GRAD(NBAST,DMAT,IPRINT)
#if defined (VAR_MPI)
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
#include "infpar.h"
#include "mpif.h"
C
      DIMENSION DMAT(NBAST,NBAST)
      IF (MYNUM .EQ. MASTER) THEN
         IPRTYP = 6
         CALL MPI_BCAST(IPRTYP,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL KSMSYNC(DMAT)
      END IF
      RETURN
#endif
      END
#if defined (VAR_MPI)
      SUBROUTINE DFT_GRAD_SLAVE(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
#include "dftmolgrad.h"
#include "energy.h"
      DIMENSION WORK(LWORK)
      EXTERNAL DFTMOLGRD 
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      DOGGA = DFT_ISGGA()
      KDMAT  = 1
      KFREE  = KDMAT + N2BASX
      IF (KFREE .GT. LWORK)CALL STOPIT('DFT_GRAD_SLAVE',' ',KFREE,LWORK)
      LFREE = LWORK - KFREE + 1
      CALL DFTINTBCAST
      CALL KSMSYNC(WORK(KDMAT))
      CALL DZERO(GRADFT,MXCOOR)
      CALL DFTINT(WORK(KDMAT),1,1,.FALSE.,WORK(KFREE),LFREE,
     &            DFTMOLGRD,WORK(KDMAT),ELE) 
      CALL GRADSLAVE_COLLECT(GRADFT,WORK(KFREE),LFREE)
      RETURN
      END
#endif
      SUBROUTINE GRADSLAVE_COLLECT(GRADMOL,WORK,LWORK)
#if defined (VAR_MPI)
#include "implicit.h"
#include "mxcent.h"
#include "mpif.h"
      DIMENSION GRADMOL(MXCOOR), WORK(LWORK)
      CALL DCOPY(MXCOOR,GRADMOL,1,WORK(1),1)
      CALL MPI_Reduce(WORK,GRADMOL,MXCOOR,MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      RETURN
#endif
      END
C
C   DFT contribution to open-shell gradients 
C   Adaptation of "old" T. Helgaker code  
C   ZR 
C 
      SUBROUTINE UDFTMOLGRDAB(GAO,DMATA,DMATB,RHA,RHB,
     &           VRA,VRB,VZA,VZB,VZG)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
C
      LOGICAL DOGGA, DFT_ISGGA
      EXTERNAL DFT_ISGGA
      DIMENSION DMATA(NBAST,NBAST), DMATB(NBAST,NBAST)
      DIMENSION GAO(NBAST,*)
      DIMENSION RHA(3),RHB(3)
C
#include <inforb.h>
#include <energy.h>
#include <symmet.h>
#include <shells.h>
C
#include <ibtfun.h>
C
      DOGGA = DFT_ISGGA()  
      VRA = D2*VRA 
      VRB = D2*VRB
      GRDNRMA = DSQRT(RHA(1)**2 + RHA(2)**2 + RHA(3)**2)
      GRDNRMB = DSQRT(RHB(1)**2 + RHB(2)**2 + RHB(3)**2)
      VZA = D2*VZA/GRDNRMA 
      VZB = D2*VZB/GRDNRMB
      VZG = D2*VZG
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
         DO IREPA = 0, MAXREP
            ISTR = IBAS(IREPA+1) + 1
            IEND = IBAS(IREPA+1) + NBAS(IREPA+1)
            IRPAX = IBTXOR(IREPA,ISYMAX(IX,1))
            IORBA = 0 
            DO ISHELA = 1, KMAX
               ISCOOR = IPTCNT(3*(NCENT(ISHELA) - 1) + IX,0,1)
               DO ICOMPA = 1, KHKT(ISHELA)
                  IORBA = IORBA + 1
                  IA = IPTSYM(IORBA,IREPA)
                  KA = IPTSYM(IORBA,IRPAX)
                  IF (KA.GT.0) THEN 
                    IF (.NOT.DOGGA) THEN    
                        GDA = D0
                        GDB = D0
                        DO IB = ISTR, IEND
                           GDA = GDA + GAO(IB,1)*DMATA(IB,IA)
                           GDB = GDB + GAO(IB,1)*DMATB(IB,IA) 
                        END DO
                        FRC = VRA*GDA*GAO(KA,IX+1)
                        FRC = FRC + VRB*GDB*GAO(KA,IX+1)
                        
                    ELSE 
                        GA   = GAO(KA,IX+1)
                        GAAX = RHA(1)*GA
                        GAAY = RHA(2)*GA
                        GAAZ = RHA(3)*GA
                        GAA2 = RHA(1)*GAO(KA,K1+4) + RHA(2)*GAO(KA,K2+4)
     &                                        + RHA(3)*GAO(KA,K3+4)
                        GDA = D0
                        GFA = D0
                        GGA = D0
C                            
                        GA   = GAO(KA,IX+1)
                        GABX = RHB(1)*GA
                        GABY = RHB(2)*GA
                        GABZ = RHB(3)*GA
                        GAB2 = RHB(1)*GAO(KA,K1+4) + RHB(2)*GAO(KA,K2+4)
     &                                        + RHB(3)*GAO(KA,K3+4)
                        GDB = D0
                        GFB = D0
                        GGB = D0
                        DO IB = ISTR, IEND
                           GAB = GAO(IB,1)*GAB2 + GAO(IB,2)*GABX +
     &                           GAO(IB,3)*GABY + GAO(IB,4)*GABZ 
                           GAA = GAO(IB,1)*GAA2 + GAO(IB,2)*GAAX + 
     &                           GAO(IB,3)*GAAY + GAO(IB,4)*GAAZ
                           GDB = GDB + DMATB(IB,IA)*GAO(IB,1)
                           GDA = GDA + DMATA(IB,IA)*GAO(IB,1)
                           GFB = GFB + DMATB(IB,IA)*GAB 
                           GFA = GFA + DMATA(IB,IA)*GAA
                           GGB = GGB + DMATB(IB,IA)*GAA 
                           GGA = GGA + DMATA(IB,IA)*GAB     
                        END DO 
                        FRC = VRA*GDA*GA+VZA*GFA+VZG*GGA
                        FRC = FRC+VRB*GDB*GA+VZB*GFB+VZG*GGB   
                    END IF   
                    GRADFT(ISCOOR) = GRADFT(ISCOOR) - FRC
                  END IF
               END DO
            END DO
         END DO
      END DO     
      RETURN
      END

 

