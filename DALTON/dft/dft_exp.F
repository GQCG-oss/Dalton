!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
! FILE : dft_exp.F
!
! Subroutines for .EXPGRA option for DFT wave functions
!
cfrj-st : parallelized exp-grad spring 2013
cfrj  : SUBROUTINE KICK_SLAVES_EXPG, DFT_EXPG_SLAVE and EXPGSLAVE_COLLECT modelled after corresponding dft_grad routines
cfrj  : the latter is apparently not used
cfrj-end
      SUBROUTINE DFTEXPGRAD(WORK,LFREE,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "maxorb.h"
#include "priunit.h"
#include "inforb.h"
#include "expopt.h"
#include "dftmolgrad.h"

      DIMENSION WORK(LFREE)
      EXTERNAL DEXPGC
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
c
      DOGGA = DFT_ISGGA()
C      
c     CALL DZERO(tmpvec,MXPRIM)
c
      NDMAT = 1
      KDMAT = 1
      KLAST = KDMAT + N2BASX
      LWRK  = LFREE - KLAST +1
      IF(KLAST.GT.LFREE) CALL QUIT('NOMEM IN DFTEXPGRAD')
      CALL DFTDNS(WORK(KDMAT),WORK(KLAST),LWRK,0)
      CALL KICK_SLAVES_EXPG(NBAST,NDMAT,WORK(KDMAT),IPRINT)
      CALL DFTINT(WORK(KDMAT),1,0,.FALSE.,WORK(KLAST),LWRK,
     &            DEXPGC,WORK(KDMAT),ELE)
c     CALL EXPGSLAVE_COLLECT(tmpvec,WORK(KLAST),LWRK)
      IF(IPRINT.GE.0) WRITE(LUPRI,'(A,F20.14)')
     &     'Orbital-exponent gradient integration, electrons:', ELE
      call flush(lupri)
      RETURN
      END
C
      SUBROUTINE KICK_SLAVES_EXPG(NBAST,NDMAT,DMAT,IPRINT)
#if defined (VAR_MPI)
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
#include "infpar.h"
#include "mpif.h"
C
      DIMENSION DMAT(NBAST,NBAST,NDMAT)
      IF (MYNUM .EQ. MASTER) THEN
cfrj: define new IPRTYP for expgrad
         IPRTYP = 666
         CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL MPI_BCAST(NDMAT,1,my_MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(DMAT,NDMAT*NBAST*NBAST,MPI_DOUBLE_PRECISION,0,
     &                  MPI_COMM_WORLD,IERR)        
      END IF
      RETURN
#endif
      END
c
#if defined (VAR_MPI)
      SUBROUTINE DFT_EXPG_SLAVE(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
#include "dftmolgrad.h"
#include "energy.h"
#include "mpif.h"
#include "expopt.h"
      DIMENSION WORK(LWORK)
      EXTERNAL DEXPGC
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      DOGGA = DFT_ISGGA()
      KDMAT  = 1
      KFREE  = KDMAT + 2*N2BASX
      IF (KFREE .GT. LWORK)CALL STOPIT('DFT_EXPG_SLAVE',' ',KFREE,LWORK)
      LFREE = LWORK - KFREE + 1
      CALL DFTINTBCAST
      CALL MPI_BCAST(NDMAT,1,my_MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WORK(KDMAT),NDMAT*NBAST*NBAST,MPI_DOUBLE_PRECISION,
     &               0,MPI_COMM_WORLD,IERR)       
c     CALL DZERO(tmpvec,MXPRIM)
      CALL DFTINT(WORK(KDMAT),NDMAT,0,.FALSE.,WORK(KFREE),LFREE,
     &               DEXPGC,WORK(KDMAT),ELE)
c     CALL EXPGSLAVE_COLLECT(tmpvec,WORK(KFREE),LFREE)
      RETURN
      END
#endif
c
      SUBROUTINE EXPGSLAVE_COLLECT(GRADMOL,WORK,LWORK)
#if defined (VAR_MPI)
#include "implicit.h"
#include "mxcent.h"
#include "mpif.h"
#include "maxorb.h"
      DIMENSION GRADMOL(MXPRIM), WORK(LWORK)
      CALL DCOPY(MXPRIM,GRADMOL,1,WORK(1),1)
      CALL MPI_Reduce(WORK,GRADMOL,MXPRIM,MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      RETURN
#endif
      END

C /* Deck dftexp */
      SUBROUTINE DEXPGC(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                  RHOA,GRADA,DST,VFA,XCPOT,COORD,WGHT,DMAT)
#include "implicit.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
#include "inforb.h"
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),
     &          RHOA(NBLEN), GRADA(3,NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          DMAT(NBAST,NBAST), DST(NATOMS), VFA(NBLEN),
     &          XCPOT(NBLEN)
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
      RETURN
      END
      SUBROUTINE DFTEXP(DMAT,GAO,GAO1,NBLEN,NBLCNT,NBLOCKS,LDAIB,
     &                  VXC,VXB,RHX,RHY,RHZ,
     &                  CORPX,CORPY,CORPZ,DOGGA)
C
C     Exchange-correlation contribution to orbital-exponent gradient
C
C     T. Helgaker dec 02 
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
C
      PARAMETER (D0=0.D0, D1=1.D0, D2=2.D0, D3=3.D0, D4=4.D0, D7=7.D0)
C
#include "inforb.h"
C
      LOGICAL DOGGA 
      DIMENSION DMAT(NBAST,NBAST), 
     &          GAO(NBLEN,NBAST), GAO1(NBLEN,NBAST,3),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8)
      DIMENSION GA10(MXAQN), GB00VC(MXAQN), 
     &          GA11(MXAQN), GB00VB(MXAQN), GB01VB(MXAQN)
C
#include "energy.h"
#include "symmet.h"
#include "shells.h"
#include "onecom.h"
#include "lmns.h"
#include "nuclei.h"
#include "primit.h"
#include "sphtrm.h"
#include "orgcom.h"
#include "dftinf.h"
#include "dftcom.h"
#include "expopt.h"
C

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
            IF (IAND(ISYMOP,MULA).EQ.0) THEN
C
C              GB00VC, GB00VB, and GB01VB
C
               DO ICOMPA = 1, KHKTA 
                  FAC = PT(IAND(ISYMOP,
     &                     IEOR(IREPA,ISYMAO(NHKTA,ICOMPA))))
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
               PAX=CORPX-PT(IAND(ISYMAX(1,1),ISYMOP))*CENT(ISHELA,1,1)
               PAY=CORPY-PT(IAND(ISYMAX(2,1),ISYMOP))*CENT(ISHELA,2,1)
               PAZ=CORPZ-PT(IAND(ISYMAX(3,1),ISYMOP))*CENT(ISHELA,3,1)
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
C -- end of dft/dft_exp.F --
