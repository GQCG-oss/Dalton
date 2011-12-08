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
      CALL KICK_SLAVES_BRHS(NGEODRV,DOVB,NBAST,WORK(KDMAT),IPRINT)
      CALL DFTINT(WORK(KDMAT),1,NGEODRV,.TRUE.,WORK(KLAST),LWRK,
     &            DFTBRHS1,WORK(KFMAT),ELE)
      CALL DFT_BRHS_COLLECT(WORK(KFMAT),N2BASX,WORK(KLAST),LWRK)
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
C
C /* Deck dftbhes */
      SUBROUTINE DFTBHES(SUSCP,WORK,LWORK,IPRINT)

C     This subroutine calculates the static DFT contributions
C     to the magnetizability (second derivative of the energy
C     with respect to the magnetic field).
C     
C     D.J. Wilson, T. Helgaker  May 2004
C     
#include <implicit.h>
      DIMENSION SUSCP(3,3)
#include <priunit.h>
#include <inforb.h>
      DIMENSION WORK(LWORK)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
      EXTERNAL DFTBHES1
C
      IF (IPRINT.GE.5) CALL TITLER('Output from DFTBHES','*',103)
C     
      KDMAT  = 1 
      KSUSDF = KDMAT  +  N2BASX
      KLAST  = KSUSDF +  9 
      LWRK   = LWORK  -  KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTBHES',' ',KLAST,LWORK)
      CALL DFTDNS(WORK(KDMAT),WORK(KLAST),LWRK,0)
      CALL DZERO(WORK(KSUSDF),9)
      CALL KICK_SLAVES_SUSCEP(N2BASX,WORK(KDMAT),IPRINT)
      CALL DFTINT(WORK(KDMAT),1,0,.TRUE.,WORK(KLAST),LWRK,DFTBHES1,
     &            WORK(KDMAT),ELE)
      CALL DFT_SUSCEP_COLLECT(WORK(KSUSDF),WORK(KLAST),LWRK)
C     
C add dft contribution to the result (SUSCP=SUSDFT probably) and print
C     
      DO I = 1, 3
      DO J = 1, 3
         SUSCP(J,I) = WORK(KSUSDF + (I-1)*3 + J - 1)
      END DO
      END DO
      IF (IPRINT.GT.10) THEN
         CALL HEADER('SUSDFT contribution to magnetizability (au)',-1)
         CALL DSCAL(9,-D1,WORK(KSUSDF),1)
         CALL POLPRI(WORK(KSUSDF),'   ',-2)
         CALL DSCAL(9,-D1,WORK(KSUSDF),1)
      END IF
      RETURN
      END
C
C /* Deck dftbhes1 */
      SUBROUTINE DFTBHES1(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,RHOA,GRADA,
     &                    COORD,WGHT,BDATA)
C
C     D.J.Wilson May 04
C
#include <implicit.h>
#include <inforb.h>
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN), WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8), NBLOCKS(2,LDAIB,8),
     &          BDATA(9 + N2BASX)
#include <dftinf.h>
      KDMAT  = 1
      KSUSDF = KDMAT  + N2BASX
      KLAST  = KSUSDF + 9
      CALL DFTBHES2(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO(1,1,1),GAO(1,1,2),
     &              RHOA,GRADA,COORD,WGHT,BDATA(KDMAT),BDATA(KSUSDF))
      RETURN
      END
C
C /* Deck dftbhes2 */
      SUBROUTINE DFTBHES2(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,GAO1,
     &                    RHOA,GRADA,COORD,WGHT,DMAT,SUSDFT)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, DP25 = 0.25D0, DP5 = 0.5D0, D2 = 2.0D0)
#include <inforb.h>
#include <orgcom.h>
#include <symmet.h>
#include <shells.h>
#include <dftcom.h>
#include "dftbrhs.h"
C
      DIMENSION NBLOCKS(2,LDAIB,8), COORD(3,NBLEN), NBLCNT(8),
     &          WGHT(NBLEN), RHOA(NBLEN), GRADA(3,NBLEN),
     &          DMAT(NBAST,NBAST), GAO(NBLEN,NBAST),
     &          GAO1(NBLEN,NBAST,3)
      DIMENSION PHASE(3), RXRHO(3), SUSDFT(3,3)
      DIMENSION VXC(NBLEN), VXB(NBLEN), VX(9), DENS(5)
C
      LOGICAL ACTIVE(NBAST)
      LOGICAL DFT_ISGGA, DOGGA
      EXTERNAL DFT_ISGGA
      DOGGA = DFT_ISGGA()
C
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
      CALL DZERO(DENS,5)
      CALL DZERO(VX,9)
      IF (DOGGA) THEN
         DO I = 1, NBLEN  
            GRDNRM = SQRT(GRADA(1,I)**2+GRADA(2,I)**2+GRADA(3,I)**2)
            DENS(1) = DP5*RHOA(I)
            DENS(2) = DP5*RHOA(I)
            DENS(3) = DP5*GRDNRM
            DENS(4) = DP5*GRDNRM
            DENS(5) = DENS(3)*DENS(4)
            CALL DFTPOT1(VX,WGHT(I),DENS(1),.FALSE.)
            VXC(I)  = DP25*VX(1)
            VXB(I)  = DP5*(VX(2)/GRDNRM + VX(9))
         END DO
      ELSE
         DO I = 1, NBLEN
            DENS(1) = DP5*RHOA(I)
            DENS(2) = DP5*RHOA(I)
            CALL DFTPOT1(VX,WGHT(I),DENS(1),.FALSE.)
            VXC(I) = DP25*VX(1)
         END DO
      END IF
C 
      DO ISYM = 1, NSYM 
         IORBA = 0        
         DO ISHELA = 1, KMAX
            CORAX  = CENT(ISHELA,1,1)
            CORAY  = CENT(ISHELA,2,1)
            CORAZ  = CENT(ISHELA,3,1)
            DO ICOMPA = 1, KHKT(ISHELA)
               IORBA = IORBA + 1
               IA = IPTSYM(IORBA,ISYM-1)
               IF (IA.GT.0 .AND. ACTIVE(IA)) THEN
                  IORBB = 0
                  DO ISHELB = 1, KMAX
                     CORBX  = CENT(ISHELB,1,1)
                     CORBY  = CENT(ISHELB,2,1)
                     CORBZ  = CENT(ISHELB,3,1)
                     ABX = CORAX - CORBX
                     ABY = CORAY - CORBY
                     ABZ = CORAZ - CORBZ
                     DO ICOMPB = 1, KHKT(ISHELB)
                        IORBB = IORBB + 1
                        IB = IPTSYM(IORBB,ISYM-1)
                        IF (IB.GT.0 .AND. ACTIVE(IB)) THEN
                           DBA = DMAT(IB,IA)
                           DO I = 1, NBLEN
                              DG   = DBA*GAO(I,IA)
                              DGVC = VXC(I)*DG*GAO(I,IB)
                              RX = COORD(1,I) - ORIGIN(1)
                              RY = COORD(2,I) - ORIGIN(2)
                              RZ = COORD(3,I) - ORIGIN(3)
                              PHASE(1) = ABY*RZ - ABZ*RY
                              PHASE(2) = ABZ*RX - ABX*RZ
                              PHASE(3) = ABX*RY - ABY*RX
                              IF (DOGGA) THEN
                                 DGVB = VXB(I)*DG*GAO(I,IB)
                                 VGA  = VXB(I)*DG*
     &                                        (GRADA(1,I)*GAO1(I,IB,1)
     &                                       + GRADA(2,I)*GAO1(I,IB,2)
     &                                       + GRADA(3,I)*GAO1(I,IB,3))
                                 RXRHO(1) = ABY*GRADA(3,I)
     &                                                  - ABZ*GRADA(2,I)
                                 RXRHO(2) = ABZ*GRADA(1,I)
     &                                                  - ABX*GRADA(3,I)
                                 RXRHO(3) = ABX*GRADA(2,I)
     &                                                  - ABY*GRADA(1,I)
                              END IF
                              DO M = 1, 3
                              DO N = 1, 3
                                 SUSDFT(M,N) = SUSDFT(M,N)
     &                                       - PHASE(M)*PHASE(N)*DGVC
                                 IF (DOGGA) THEN
                                    SUSDFT(M,N) = SUSDFT(M,N)
     &                                   - RXRHO(M)*PHASE(N)*DGVB
     &                                   - PHASE(M)*PHASE(N)*VGA
                                 END IF
                              END DO
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE KICK_SLAVES_BRHS(NGEODRV,DOVB,NBAST,DMAT,IPRINT)
#if defined (VAR_MPI)
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
C
      LOGICAL DOVB
      DIMENSION DMAT(NBAST,NBAST)
C
      IF (MYNUM .EQ. MASTER) THEN
         IPRTYP = 7
         CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DOVB,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(NGEODRV,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL KSMSYNC(DMAT)
      END IF
      RETURN
#endif
      END
C
#if defined (VAR_MPI)
      SUBROUTINE DFT_BRHS_SLAVE(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
#include "dftbrhs.h"
C
      DIMENSION WORK(LWORK)
#include "dftcom.h"
#include "mpif.h"
      EXTERNAL DFTBRHS1
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      KDMAT = 1
      KFMAT = KDMAT + N2BASX
      KLAST = KFMAT + N2BASX*3
      LWRK  = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL QUIT('Out of memory in DFT_BRHS_SLAVE')
      CALL MPI_BCAST(DOVB,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NGEODRV,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL DFTINTBCAST
      CALL KSMSYNC(WORK(KDMAT))
      CALL DZERO(WORK(KFMAT),3*N2BASX)
      CALL DFTINT(WORK(KDMAT),1,NGEODRV,.TRUE.,WORK(KLAST),LWRK,
     &            DFTBRHS1,WORK(KFMAT),ELE)
      CALL DFT_BRHS_COLLECT(WORK(KFMAT),N2BASX,WORK(KLAST),LWRK)
      RETURN
      END
#endif
C
      SUBROUTINE DFT_BRHS_COLLECT(FMAT,N2BASX,WORK,LWORK)
#if defined (VAR_MPI)
#include "implicit.h"
#include "mxcent.h"
#include "mpif.h"
      DIMENSION FMAT(3*N2BASX), WORK(LWORK)
      CALL DCOPY(3*N2BASX,FMAT,1,WORK,1)
      CALL MPI_REDUCE(WORK,FMAT,3*N2BASX,MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      RETURN
#endif
      END
C
      SUBROUTINE KICK_SLAVES_SUSCEP(N2BASX,DMAT,IPRINT)
#if defined (VAR_MPI)
#include "implicit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
C
      DIMENSION DMAT(N2BASX)
      IF (MYNUM .EQ. MASTER) THEN
         IPRTYP = 10
         CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL KSMSYNC(DMAT)
      END IF
      RETURN
#endif
      END
C
#if defined (VAR_MPI)
      SUBROUTINE DFT_SUSCEP_SLAVE(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
C
      DIMENSION WORK(LWORK)
      EXTERNAL DFTBHES1
C
      KDMAT  = 1 
      KSUSDF = KDMAT  +  N2BASX
      KLAST  = KSUSDF +  9 
      LWRK   = LWORK  -  KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('DFT_SUSCEP_SLAVE',' ',
     &                                  KLAST,LWORK)
      CALL DFTINTBCAST
      CALL KSMSYNC(WORK(KDMAT))
      CALL DZERO(WORK(KSUSDF),9)
      CALL DFTINT(WORK(KDMAT),1,0,.TRUE.,WORK(KLAST),LWRK,DFTBHES1,
     &            WORK(KDMAT),ELE)
      CALL DFT_SUSCEP_COLLECT(WORK(KSUSDF),WORK(KLAST),LWRK)
      RETURN
      END
#endif
C
      SUBROUTINE DFT_SUSCEP_COLLECT(SUSCEP,WORK,LWORK)
#if defined (VAR_MPI)
#include "implicit.h"
#include "mpif.h"
      DIMENSION SUSCEP(3,3), WORK(LWORK)
      CALL DCOPY(9,SUSCEP,1,WORK,1)
      CALL MPI_REDUCE(WORK,SUSCEP,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                MPI_COMM_WORLD,IERR)
      RETURN
#endif
      END
