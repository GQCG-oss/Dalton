C /* Deck dftexp */
      SUBROUTINE DFTEXP(DMAT,GAO,GAO1,GAO2,VXC,VXB,RHX,RHY,RHZ,
     &                  CORPX,CORPY,CORPZ,DOGGA)
C
C     Exchange-correlation contribution to orbital-exponent gradient
C
C     T. Helgaker dec 02 
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
C
      PARAMETER (D0=0.D0, D1=1.D0, D2=2.D0, D3=3.D0, D4=4.D0, D7=7.D0)
C
      LOGICAL DOGGA 
      DIMENSION DMAT(NBAST,NBAST), 
     &          GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6)
      DIMENSION GA10(MXAQN), GB00VC(MXAQN), 
     &          GA11(MXAQN), GB00VB(MXAQN), GB01VB(MXAQN)
C
#include <inforb.h>
#include <energy.h>
#include <symmet.h>
#include <shells.h>
#include <onecom.h>
#include <lmns.h>
#include <nuclei.h>
#include <primit.h>
#include <sphtrm.h>
#include <orgcom.h>
#include <dftinf.h>
#include <dftcom.h>
#include <expopt.h>
C
#include <sdpre.h>
#include <ibtfun.h>
C
      DO IREPA = 0, MAXREP
         ISTR = IBAS(IREPA+1) + 1
         IEND = IBAS(IREPA+1) + NBAS(IREPA+1)
         IORBA = 0
         DO ISHELA = 1, KMAX
            NHKTA  = NHKT(ISHELA)
            KHKTA  = KHKT(ISHELA)
            KCKTA  = KCKT(ISHELA)
            SPHRA  = SPHR(ISHELA)
            NUMCFA = NUMCF(ISHELA)
            JSTA   = JSTRT(ISHELA)
            NUCA   = NUCO(ISHELA)
            MULA   = ISTBAO(ISHELA)
            DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULA).EQ.0) THEN
C
C              GB00VC, GB00VB, and GB01VB
C
               DO ICOMPA = 1, KHKTA 
                  FAC = PT(IBTAND(ISYMOP,
     &                     IBTXOR(IREPA,ISYMAO(NHKTA,ICOMPA))))
                  IA = IPTSYM(IORBA+ICOMPA,IREPA)
                  GD = D0 
                  DO IB = ISTR, IEND 
                     GD = GD + GAO(IB)*DMAT(IB,IA)
                  END DO
                  GB00VC(ICOMPA) = D2*FAC*VXC*GD
                  IF (DOGGA) THEN
                     G0 = D0 
                     G1 = D0 
                     DO IB = ISTR, IEND 
                        G0 = G0 + DMAT(IB,IA)*GAO(IB)
                        G1 = G1 + DMAT(IB,IA)*(RHX*GAO1(IB,1) 
     &                                       + RHY*GAO1(IB,2)
     &                                       + RHZ*GAO1(IB,3))
                     END DO
                     GB00VB(ICOMPA) = FAC*VXB*G0
                     GB01VB(ICOMPA) = FAC*VXB*G1
                  END IF
               END DO
C
               PAX=CORPX-PT(IBTAND(ISYMAX(1,1),ISYMOP))*CENT(ISHELA,1,1)
               PAY=CORPY-PT(IBTAND(ISYMAX(2,1),ISYMOP))*CENT(ISHELA,2,1)
               PAZ=CORPZ-PT(IBTAND(ISYMAX(3,1),ISYMOP))*CENT(ISHELA,3,1)
               PX2 = PAX**2
               PY2 = PAY**2
               PZ2 = PAZ**2
               PA2 = PX2 + PY2 + PZ2
               CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
C
               DO I = JSTA + 1, JSTA + NUCA
                  ALPHA = PRIEXP(I)
                  GA = PRICCF(I,NUMCFA)*DEXP(-ALPHA*PA2)
C
                  IF (SPHRA) CALL DZERO(GA10,KHKTA)
                  DO ICOMPA = 1, KCKTA
                     L = LVALUA(ICOMPA)
                     M = MVALUA(ICOMPA)
                     N = NVALUA(ICOMPA)
                     C1 = (D2*SDPRE(L+M+N) + D3)/(D4*ALPHA) - PA2
                     CI = (PAX**L)*(PAY**M)*(PAZ**N)*C1*GA
                     IF (SPHRA) THEN
                        IOFF = ISPADR(NHKTA) + (ICOMPA-1)*KHKTA - 1 
                        DO K = 1, KHKTA
                           GA10(K) = GA10(K) + CSP(IOFF+K)*CI
                        END DO
                     ELSE
                        GA10(ICOMPA) = CI 
                     END IF
                  END DO
                  DO ICOMPA = 1, KHKTA 
                     ALPGRD(I) = ALPGRD(I) + GA10(ICOMPA)*GB00VC(ICOMPA)
                  END DO
C
                  IF (DOGGA) THEN
                     IF (SPHRA) CALL DZERO(GA11,KHKTA)
                     DO ICOMPA = 1, KCKTA
                        L = LVALUA(ICOMPA)
                        M = MVALUA(ICOMPA)
                        N = NVALUA(ICOMPA)
                        TR = GA*((D2*SDPRE(L+M+N) + D3)/(D4*ALPHA)-PA2)
                        TS = -D2*(ALPHA*TR + GA)
                        IF (L.EQ.0) THEN
                           TX = TS*PAX*(PAY**M)*(PAZ**N)
                        ELSE
                           TX = (PX2*TS + SDPRE(L)*TR)
     &                         *(PAX**(L-1))*(PAY**M)*(PAZ**N)
                        END IF
                        IF (M.EQ.0) THEN
                           TY = TS*(PAX**L)*PAY*(PAZ**N)
                        ELSE
                           TY = (PY2*TS + SDPRE(M)*TR)
     &                         *(PAX**L)*(PAY**(M-1))*(PAZ**N)
                        END IF
                        IF (N.EQ.0) THEN
                           TZ = TS*(PAX**L)*(PAY**M)*PAZ
                        ELSE
                           TZ = (PZ2*TS + SDPRE(N)*TR)
     &                         *(PAX**L)*(PAY**M)*(PAZ**(N-1))
                        END IF
                        CI = RHX*TX + RHY*TY + RHZ*TZ
                        IF (SPHRA) THEN
                           IOFF = ISPADR(NHKTA) + (ICOMPA-1)*KHKTA - 1 
                           DO K = 1, KHKTA
                              GA11(K) = GA11(K) + CSP(IOFF+K)*CI
                           END DO
                        ELSE
                           GA11(ICOMPA) = CI 
                        END IF
                     END DO
                     DO ICOMPA = 1, KHKTA 
                        ALPGRD(I)=ALPGRD(I)+GA11(ICOMPA)*GB00VB(ICOMPA)
     &                                     +GA10(ICOMPA)*GB01VB(ICOMPA)
                     END DO
                  END IF
               END DO
            END IF
            END DO
            IORBA = IORBA + KHKTA
         END DO
      END DO
      RETURN
      END
