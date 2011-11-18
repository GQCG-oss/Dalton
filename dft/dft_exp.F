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
      SUBROUTINE DFTEXPGRAD(WORK,LFREE,IPRINT)
#include <implicit.h>
#include <dummy.h>
#include <inforb.h>
#include <priunit.h>
      EXTERNAL DEXPGC
      DIMENSION WORK(LFREE)
      KDMAT = 1
      KLAST = KDMAT + N2BASX
      LWRK  = LFREE - KLAST +1
      IF(KLAST.GT.LFREE) CALL QUIT('NOMEM IN DFTEXPGRAD')
      CALL DFTDNS(WORK(KDMAT),WORK(KLAST),LWRK,0)
      CALL DFTINT(WORK(KDMAT),1,0,.FALSE.,WORK(KLAST),LWRK,
     &            DEXPGC,WORK(KDMAT),ELE)
      IF(IPRINT.GE.0) WRITE(LUPRI,'(A,F20.14)')
     &     'Orbital-exponent gradient integration, electrons:', ELE
      END

C /* Deck dftexp */
      SUBROUTINE DEXPGC(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                  RHOA,GRADA,COORD,WGHT,DMAT)
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
#include <inforb.h>
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          DMAT(NBAST,NBAST)
      LOGICAL DFT_ISGGA,DOGGA
      EXTERNAL DFT_ISGGA
      DIMENSION VX(5), VXC(NBLEN), VXB(NBLEN)
C     
      DOGGA = DFT_ISGGA()
C     
      IF (DOGGA) THEN
         DO I = 1, NBLEN
            GRDNRM = SQRT(GRADA(1,I)**2 + GRADA(2,I)**2 + GRADA(3,I)**2)
            CALL DFTPTF0(RHOA(I),GRDNRM,WGHT(I),VX)
            VXC(I) = VX(1)
            VXB(I) = D2*VX(2)/GRDNRM
         END DO
      ELSE
         DO I = 1, NBLEN
            CALL DFTPTF0(RHOA(I),D0,WGHT(I),VX)
            VXC(I) = VX(1)
         END DO
      END IF
      DO I = 1, NBLEN
         CALL DFTEXP(DMAT,GAO(I,1,1),GAO(I,1,2),
     &               NBLEN,NBLCNT,NBLOCKS,LDAIB,
     &               VXC(I),VXB(I),
     &               GRADA(1,I),GRADA(2,I),GRADA(3,I),
     &               COORD(1,I),COORD(2,I),COORD(3,I),DOGGA)
      END DO
      END
      SUBROUTINE DFTEXP(DMAT,GAO,GAO1,NBLEN,NBLCNT,NBLOCKS,LDAIB,
     &                  VXC,VXB,RHX,RHY,RHZ,
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
#include <inforb.h>
C
      LOGICAL DOGGA 
      DIMENSION DMAT(NBAST,NBAST), 
     &          GAO(NBLEN,NBAST), GAO1(NBLEN,NBAST,3),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8)
      DIMENSION GA10(MXAQN), GB00VC(MXAQN), 
     &          GA11(MXAQN), GB00VB(MXAQN), GB01VB(MXAQN)
C
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
#include <ibtfun.h>
C
      DO IREPA = 0, MAXREP
      DO IBL = 1, NBLCNT(IREPA+1)
         ISTR  = NBLOCKS(1,IBL,IREPA+1)
         IEND  = NBLOCKS(2,IBL,IREPA+1)
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
                     GD = GD + GAO(1,IB)*DMAT(IB,IA)
                  END DO
                  GB00VC(ICOMPA) = D2*FAC*VXC*GD
                  IF (DOGGA) THEN
                     G0 = D0 
                     G1 = D0 
                     DO IB = ISTR, IEND 
                        G0 = G0 + DMAT(IB,IA)*GAO(1,IB)
                        G1 = G1 + DMAT(IB,IA)*(RHX*GAO1(1,IB,1) 
     &                                       + RHY*GAO1(1,IB,2)
     &                                       + RHZ*GAO1(1,IB,3))
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
                     C1 = (D2*(L+M+N) + D3)/(D4*ALPHA) - PA2
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
                        TR = GA*((D2*(L+M+N) + D3)/(D4*ALPHA)-PA2)
                        TS = -D2*(ALPHA*TR + GA)
                        IF (L.EQ.0) THEN
                           TX = TS*PAX*(PAY**M)*(PAZ**N)
                        ELSE
                           TX = (PX2*TS + L*TR)
     &                         *(PAX**(L-1))*(PAY**M)*(PAZ**N)
                        END IF
                        IF (M.EQ.0) THEN
                           TY = TS*(PAX**L)*PAY*(PAZ**N)
                        ELSE
                           TY = (PY2*TS + M*TR)
     &                         *(PAX**L)*(PAY**(M-1))*(PAZ**N)
                        END IF
                        IF (N.EQ.0) THEN
                           TZ = TS*(PAX**L)*(PAY**M)*PAZ
                        ELSE
                           TZ = (PZ2*TS + N*TR)
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
      END DO
      RETURN
      END
