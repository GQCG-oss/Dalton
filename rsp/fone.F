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
C  /* Deck fonedn */
      SUBROUTINE FONEDN(NSIM,FCOCO,FVOCO,FCOEX,FVOEX)
C 1)
C SCALE FCOCO AND FVOEX TO HAVE RIGHT FACTORS
C 2)
C ADD CONTRIBUTIONS FROM OCCUPIED SECONDARY TO FCOCO AND FVOCO
C (USE THAT COULOMB CONTRIBUTION IS SYMMETRIC) AND
C 3)
C SUM UP ONE INDEX TRANSFORMED DENSITY CONTRIBUTIONS TO FOCK
C MATRICES ( FCOEX  = FCOCO  + FCOEX )
C
#include <implicit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <wrkrsp.h>
#include <infrsp.h>
C
      DIMENSION FCOCO(NORBT,NORBT,*),FVOCO(NORBT,NORBT,*)
      DIMENSION FCOEX(NORBT,NORBT,*),FVOEX(NORBT,NORBT,*)
C
      PARAMETER ( D2 = 2.0D0  , DP5 = 0.5D0 )
C
      IF (NASHT.GT.0) CALL DSCAL(NORBT*NORBT*NSIM,DP5,FVOEX,1)
      IF (DODFT .OR. DOHFSRDFT .OR. DOMCSRDFT)
     &   CALL DSCAL(NORBT*NORBT*NSIM,HFXFAC,FCOEX,1)
      IF (TRPLET) RETURN
      IF (NISHT.GT.0) CALL DSCAL(NORBT*NORBT*NSIM,D2,FCOCO,1)
      DO 100 IPSYM = 1,NSYM
         IORBP = IORB(IPSYM) + 1
         NOCCP = NOCC(IPSYM)
         IQSYM = MULD2H(IPSYM,KSYMOP)
         IORBQ = IORB(IQSYM) + 1
         NOCCQ = NOCC(IQSYM)
         NORBQ = NORB(IQSYM)
         IF ((NOCCP.EQ.0).OR.((NORBQ-NOCCQ).EQ.0)) GO TO 100
         DO 400 ISIM = 1,NSIM
            DO 200 IP = IORBP,IORBP+NOCCP-1
               DO 300 IQ = IORBQ+NOCCQ,IORBQ+NORBQ-1
                  FCOCO(IP,IQ,ISIM) = FCOCO(IQ,IP,ISIM)
                  IF (NASHT.GT.0) FVOCO(IP,IQ,ISIM) = FVOCO(IQ,IP,ISIM)
 300           CONTINUE
 200        CONTINUE
 400     CONTINUE
 100  CONTINUE
      DO 500 IPSYM = 1,NSYM
         IORBP = IORB(IPSYM) + 1
         NORBP = NORB(IPSYM)
         IQSYM = MULD2H(IPSYM,KSYMOP)
         IORBQ = IORB(IQSYM) + 1
         NORBQ = NORB(IQSYM)
         IF ((NORBQ.EQ.0).OR.(NORBP.EQ.0)) GO TO 500
         DO 800 ISIM = 1,NSIM
            DO 600 IP = IORBP,IORBP+NORBP-1
               DO 700 IQ = IORBQ,IORBQ+NORBQ-1
                  FCOEX(IP,IQ,ISIM) = FCOCO(IP,IQ,ISIM)
     *                              + FCOEX(IP,IQ,ISIM)
                  IF (NASHT.GT.0) FVOEX(IP,IQ,ISIM) = FVOCO(IP,IQ,ISIM)
     *                              + FVOEX(IP,IQ,ISIM)
 700           CONTINUE
 600        CONTINUE
 800     CONTINUE
 500  CONTINUE
      RETURN
      END
C  /* Deck fonemu */
      SUBROUTINE FONEMU(NSIM,ICI1,IDI1,H2,
     *                  FCOEX,FVOEX,ZYMAT,DENA,DENB,WRK,LWRK)
C
C CALCULATE CONTRIBUTIONS TO INACTIVE AND ACTIVE FOCK MATRICES
C WHICH ORIGINATE FROM OCCUPIED-OCCUPIED COULOMN DISTRIBUTIONS
C (ONLY THE EXCHANGE PART OF THE ONE-INDEX TRANSFORMED
C  DENSITY CONTRIBUTE)
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) (JQ/PR)*ZYMAT(R,J)  (Q OCC)
C                          - SUM(R) (PJ/RQ)*ZYMAT(J,R)  (P OCC)
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) (XQ/PR)*DENA(R,X)   (Q OCC)
C                          - SUM(R) (PY/QR)*DENB(R,Y)   (P OCC)
C
#include <implicit.h>
C
C Used from common blocks:
C   INFORB :
C   INFIND : ISMO,IOBTYP
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <orbtypdef.h>
C
      DIMENSION H2(NORBT,*)
      DIMENSION FCOEX(NORBT,NORBT,*),FVOEX(NORBT,NORBT,*)
      DIMENSION ZYMAT(NORBT,NORBT,*)
      DIMENSION DENA(NORBT,NASHT,*),DENB(NORBT,NASHT,*),WRK(*)
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
      ELSE
         ICI = IDI1
         IDI = ICI1
      END IF
C
C     Find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need occupied-occupied distributions, i.e. itypcd .le. 3
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPCD .GT. 3) RETURN
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      IOFFC = IORB(ICSYM) + 1
      IOFFD = IORB(IDSYM) + 1
      ICDSYM = MULD2H(ICSYM,IDSYM)
      DO 100 ISIM = 1,NSIM
         DO 200 IRSYM = 1,NSYM
            IPSYM = MULD2H(IRSYM,ICDSYM)
            IOFFR = IORB(IRSYM) + 1
            IOFFP = IORB(IPSYM) + 1
            NORBR = NORB(IRSYM)
            NORBP = NORB(IPSYM)
            IF ( (NORBR.EQ.0) .OR. (NORBP.EQ.0) ) GO TO 200
C
            ICRSYM = MULD2H(ICSYM,IRSYM)
            IDRSYM = MULD2H(IDSYM,IRSYM)
            IF ( KSYMOP.EQ.IDRSYM ) THEN
C
               IF ((ITYPCD.EQ.1) .OR. (ITYPCD.EQ.2)) THEN
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) (JQ/PR)*ZYMAT(R,J)  (Q OCC)
C                                    QJ/PR
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       ZYMAT(IOFFR,IDI,ISIM),NORBT,1.D0,
     &                       FCOEX(IOFFP,ICI,ISIM),NORBT)
C
C                          - SUM(R) (PJ/RQ)*ZYMAT(J,R)  (P OCC)
C
                  CALL DGEMM('N','N',1,NORBP,NORBR,-1.D0,
     &                       ZYMAT(IDI,IOFFR,ISIM),NORBT,
     &                       H2(IOFFR,IOFFP),NORBT,1.D0,
     &                       FCOEX(ICI,IOFFP,ISIM),NORBT)
               END IF
            END IF
            IF (ICRSYM.EQ.KSYMOP) THEN
               IF ((ITYPCD.EQ.1) .AND. (ICI.NE.IDI)) THEN
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) (JQ/PR)*ZYMAT(R,J)  (Q OCC)
C                                    QJ/PR
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       ZYMAT(IOFFR,ICI,ISIM),NORBT,1.D0,
     &                       FCOEX(IOFFP,IDI,ISIM),NORBT)
C
C                          - SUM(R) (PJ/RQ)*ZYMAT(J,R)  (P OCC)
C
                  CALL DGEMM('N','N',1,NORBP,NORBR,-1.D0,
     &                       ZYMAT(ICI,IOFFR,ISIM),NORBT,
     &                       H2(IOFFR,IOFFP),NORBT,1.D0,
     &                       FCOEX(IDI,IOFFP,ISIM),NORBT)
               END IF
            END IF
            IF ( ICRSYM.EQ.KSYMOP ) THEN
               IF ((ITYPCD.EQ.2).OR.(ITYPCD.EQ.3)) THEN
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) (XQ/PR)*DENA(R,X)   (Q OCC)
C
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       DENA(IOFFR,NCIW,ISIM),NORBT,1.D0,
     &                       FVOEX(IOFFP,IDI,ISIM),NORBT)
C
C                          - SUM(R) (PY/QR)*DENB(R,Y)   (P OCC)
C                                    YP/QR
C
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       DENB(IOFFR,NCIW,ISIM),NORBT,0.D0,
     &                       WRK(1),NORBT)
                  DO 450 IQ = IOFFP,IOFFP+NORBP-1
                     FVOEX(IDI,IQ,ISIM) = FVOEX(IDI,IQ,ISIM)
     *                                       - WRK(IQ-IOFFP+1)
 450              CONTINUE
               END IF
            END IF
            IF (IDRSYM.EQ.KSYMOP) THEN
               IF ((ITYPCD.EQ.3).AND.(ICI.NE.IDI)) THEN
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) (XQ/PR)*DENA(R,X)   (Q OCC)
C
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       DENA(IOFFR,NDIW,ISIM),NORBT,1.D0,
     &                       FVOEX(IOFFP,ICI,ISIM),NORBT)
C
C                          - SUM(R) (PY/QR)*DENB(R,Y)   (P OCC)
C                                    YP/QR
C
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       DENB(IOFFR,NDIW,ISIM),NORBT,0.D0,
     &                       WRK(1),NORBT)
                  DO 451 IQ = IOFFP,IOFFP+NORBP-1
                     FVOEX(ICI,IQ,ISIM) = FVOEX(ICI,IQ,ISIM)
     *                                       - WRK(IQ-IOFFP+1)
 451              CONTINUE
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck fonedr */
      SUBROUTINE FONEDR(NSIM,ICI1,IDI1,H2D,FCOCO,FVOCO,
     *                  FCOEX,FVOEX,ZYMAT,DENA,DENB,WRK,LWRK)
C
C CALCULATE CONTRIBUTIONS TO INACTIVE AND ACTIVE FOCK MATRICES
C WHICH ORIGINATE FROM OCCUPIED-OCCUPIED DIRAC DISTRIBUTIONS
C
C F(CV)OEX CONTAIN EXCHANCE CONTRIBUTION
C F(CV)OMU CONTAIN COULOMN CONTRIBUTIONS
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) <JP/QR>*ZYMAT(R,J) P OCC Q SEC
C                          - SUM(R) <QJ/RP>*ZYMAT(J,R) Q OCC P SEC
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) <XP/QR>*DENA(R,X)  P OCC Q SEC
C                          - SUM(R) <QY/RP>*DENB(R,Y)  Q OCC P SEC
C
C  FCOCO(P,Q) = FCOCO(P,Q) + SUM(R) <QJ/PR>*(ZYMAT(J,R)-ZYMAT(R,J))
C
C  FVOCO(P,Q) = FVOCO(P,Q) + SUM(R) <QY/PR>*(DENB(R,Y)-DENA(R,Y))
C
#include <implicit.h>
C
C Used from common blocks:
C   INFORB :
C   INFIND : ISMO,IOBTYP
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <orbtypdef.h>
      DIMENSION H2D(NORBT,*)
      DIMENSION FCOCO(NORBT,NORBT,*),FVOCO(NORBT,NORBT,*)
      DIMENSION FCOEX(NORBT,NORBT,*),FVOEX(NORBT,NORBT,*)
      DIMENSION ZYMAT(NORBT,NORBT,*)
      DIMENSION DENA(NORBT,NASHT,*),DENB(NORBT,NASHT,*),WRK(*)
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      NDITR = 1
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
      ELSE
         ICI = IDI1
         IDI = ICI1
         CALL DGETRN(H2D,NORBT,NORBT)
         NDITR = -NDITR
      END IF
C
C     Find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need occupied-occupied distributions, i.e. itypcd .le. 3
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPCD .GT. 3) RETURN
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      IOFFC = IORB(ICSYM) + 1
      IOFFD = IORB(IDSYM) + 1
      ICDSYM = MULD2H(ICSYM,IDSYM)
      DO 100 ISIM = 1,NSIM
         DO 200 IRSYM = 1,NSYM
            IPSYM = MULD2H(IRSYM,ICDSYM)
            NORBR = NORB(IRSYM)
            NSSHP = NSSH(IPSYM)
            IF (NORBR.EQ.0) GO TO 200
            IOFFR = IORB(IRSYM) + 1
            IOFFP = IORB(IPSYM) + 1
            NORBP = NORB(IPSYM)
            IF (NORBP.EQ.0) GO TO 200
            NOCCP = NOCC(IPSYM)
            IDRSYM = MULD2H(IDSYM,IRSYM)
            ICRSYM = MULD2H(ICSYM,IRSYM)
            IF ( IDRSYM.EQ.KSYMOP) THEN
               IF ((.NOT.TRPLET).AND.(ITYPCD.LE.2)) THEN
C
C  FCOCO(P,Q) = FCOCO(P,Q) + SUM(R) <QJ/PR>*(ZYMAT(J,R)-ZYMAT(R,J))
C
                  DO 1100 IR=1,NORBR
                     WRK(IR) = ZYMAT(IDI,IOFFR-1+IR,ISIM)-
     *                         ZYMAT(IOFFR-1+IR,IDI,ISIM)
 1100             CONTINUE
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2D(IOFFP,IOFFR),NORBT,
     &                       WRK,NORBT,1.D0,
     &                       FCOCO(IOFFP,ICI,ISIM),NORBT)
               END IF
            END IF
            IF ( ICRSYM.EQ.KSYMOP) THEN
               IF ((.NOT.TRPLET).AND.(ITYPCD.GE.2)) THEN
C
C  FVOCO(P,Q) = FVOCO(P,Q) + SUM(R) <QY/PR>*(DENB(R,Y)-DENA(R,Y))
C                                    YQ/RP        R,Y       R,Y
C
                  DO 1200 IR=1,NORBR
                     WRK(IR) = DENB(IOFFR-1+IR,NCIW,ISIM)-
     *                         DENA(IOFFR-1+IR,NCIW,ISIM)
 1200             CONTINUE
                  CALL DGEMM('T','N',NORBP,1,NORBR,1.D0,
     &                       H2D(IOFFR,IOFFP),NORBT,
     &                       WRK,NORBT,1.D0,
     &                       FVOCO(IOFFP,IDI,ISIM),NORBT)
               END IF
            END IF
            IF (NSSHP.EQ.0) GO TO 200
            IF ( IDRSYM.EQ.KSYMOP) THEN
               IF ((ITYPCD.EQ.1) .OR. (ITYPCD.EQ.2)) THEN
C
C  FCOEX(P,Q) = FCOEX(P,Q) - SUM(R) <QJ/RP>*ZYMAT(J,R) Q OCC P SEC
C

                  DO 111 IR=1,NORBR
                     WRK(IR) =-ZYMAT(IDI,IOFFR-1+IR,ISIM)
 111              CONTINUE
                  CALL DGEMM('T','N',NSSHP,1,NORBR,1.D0,
     &                       H2D(IOFFR,IOFFP+NOCCP),NORBT,
     &                       WRK(1),NORBR,1.D0,
     &                       FCOEX(IOFFP+NOCCP,ICI,ISIM),NORBT)
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) <JP/QR>*ZYMAT(R,J) P OCC Q SEC
C                                    PJ RQ
C
                  CALL DGEMM('T','N',1,NSSHP,NORBR,1.D0,
     &                       ZYMAT(IOFFR,IDI,ISIM),NORBT,
     &                       H2D(IOFFR,IOFFP+NOCCP),NORBT,1.D0,
     &                       FCOEX(ICI,IOFFP+NOCCP,ISIM),NORBT)
               END IF
            END IF
            IF ( ICRSYM.EQ.KSYMOP) THEN
               IF ((ITYPCD.EQ.2).OR.(ITYPCD.EQ.3)) THEN
C
C  FVOEX(P,Q) = FVOEX(P,Q) - SUM(R) <QY/RP>*DENB(R,Y) Q OCC P SEC
C                                    YQ/PR       R,Y
                  CALL DGEMM('N','N',NSSHP,1,NORBR,-1.D0,
     &                       H2D(IOFFP+NOCCP,IOFFR),NORBT,
     &                       DENB(IOFFR,NCIW,ISIM),NORBT,1.D0,
     &                       FVOEX(IOFFP+NOCCP,IDI,ISIM),NORBT)
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) <XP/QR>*DENA(R,X)   (P OCC)
C
                  CALL DGEMM('N','N',NSSHP,1,NORBR,1.D0,
     &                       H2D(IOFFP+NOCCP,IOFFR),NORBT,
     &                       DENA(IOFFR,NCIW,ISIM),NORBT,0.D0,
     &                       WRK(1),NORBT)
                  DO 350 IQ = IOFFP+NOCCP,IOFFP+NORBP-1
                     IQ1 = IQ - IOFFP - NOCCP + 1
                     FVOEX(IDI,IQ,ISIM) = FVOEX(IDI,IQ,ISIM) + WRK(IQ1)
 350              CONTINUE
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
      IF (((ITYPCD.EQ.1).OR.(ITYPCD.EQ.3)) .AND. (ICI.NE.IDI) ) THEN
         CALL DGETRN(H2D,NORBT,NORBT)
         NDITR = -NDITR
         DO 110 ISIM = 1,NSIM
            DO 210 IRSYM = 1,NSYM
               IPSYM = MULD2H(IRSYM,ICDSYM)
               NORBR = NORB(IRSYM)
               NSSHP = NSSH(IPSYM)
               IF (NORBR.EQ.0) GO TO 210
               IOFFR = IORB(IRSYM) + 1
               IOFFP = IORB(IPSYM) + 1
               NORBP = NORB(IPSYM)
               IF (NORBP.EQ.0) GO TO 210
               NOCCP = NOCC(IPSYM)
C
               IDRSYM = MULD2H(IDSYM,IRSYM)
               ICRSYM = MULD2H(ICSYM,IRSYM)
               IF (ICRSYM.EQ.KSYMOP) THEN
                  IF ((.NOT.TRPLET).AND.(ITYPCD.EQ.1)) THEN
C
C  FCOCO(P,Q) = FCOCO(P,Q) + SUM(R) <QJ/PR>*(ZYMAT(J,R)-ZYMAT(R,J))
C
                     DO 1101 IR=1,NORBR
                        WRK(IR) = ZYMAT(ICI,IOFFR-1+IR,ISIM)-
     *                         ZYMAT(IOFFR-1+IR,ICI,ISIM)
 1101                CONTINUE
                     CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                          H2D(IOFFP,IOFFR),NORBT,
     &                          WRK,NORBT,1.D0,
     &                          FCOCO(IOFFP,IDI,ISIM),NORBT)
C
                  END IF
               END IF
               IF (IDRSYM.EQ.KSYMOP) THEN
                  IF ((.NOT.TRPLET).AND.(ITYPCD.EQ.3)) THEN
C
C  FVOCO(P,Q) = FVOCO(P,Q) + SUM(R) <QY/PR>*(DENB(R,Y)-DENA(R,Y))
C                                    YQ/RP        R,Y       R,Y
C
                     DO 1201 IR=1,NORBR
                        WRK(IR) = DENB(IOFFR-1+IR,NDIW,ISIM)-
     *                            DENA(IOFFR-1+IR,NDIW,ISIM)
 1201                CONTINUE
                     CALL DGEMM('T','N',NORBP,1,NORBR,1.D0,
     &                          H2D(IOFFR,IOFFP),NORBT,
     &                          WRK,NORBT,1.D0,
     &                          FVOCO(IOFFP,ICI,ISIM),NORBT)
C
                  END IF
               END IF
               IF (NSSHP.EQ.0) GO TO 210
               IF (ICRSYM.EQ.KSYMOP) THEN
                  IF (ITYPCD.EQ.1) THEN
C
C  FCOEX(P,Q) = FCOEX(P,Q) - SUM(R) <QJ/RP>*ZYMAT(J,R) Q OCC P SEC
C

                     DO 112 IR=1,NORBR
                        WRK(IR) =-ZYMAT(ICI,IOFFR-1+IR,ISIM)
 112                 CONTINUE
                     CALL DGEMM('T','N',NSSHP,1,NORBR,1.D0,
     &                          H2D(IOFFR,IOFFP+NOCCP),NORBT,
     &                          WRK(1),NORBR,1.D0,
     &                          FCOEX(IOFFP+NOCCP,IDI,ISIM),NORBT)
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) <JP/QR>*ZYMAT(R,J) P OCC Q SEC
C                                    PJ RQ
C
                     CALL DGEMM('T','N',1,NSSHP,NORBR,1.D0,
     &                          ZYMAT(IOFFR,ICI,ISIM),NORBT,
     &                          H2D(IOFFR,IOFFP+NOCCP),NORBT,1.D0,
     &                          FCOEX(IDI,IOFFP+NOCCP,ISIM),NORBT)
C
                  END IF
               END IF
               IF (IDRSYM.EQ.KSYMOP) THEN
                  IF (ITYPCD.EQ.3) THEN
C
C  FVOEX(P,Q) = FVOEX(P,Q) - SUM(R) <QY/RP>*DENB(R,Y) Q OCC P SEC
C                                    YQ/PR       R,Y
                     CALL DGEMM('N','N',NSSHP,1,NORBR,-1.D0,
     &                          H2D(IOFFP+NOCCP,IOFFR),NORBT,
     &                          DENB(IOFFR,NDIW,ISIM),NORBT,1.D0,
     &                          FVOEX(IOFFP+NOCCP,ICI,ISIM),NORBT)
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) <XP/QR>*DENA(R,X)   (P OCC)
C
                     CALL DGEMM('N','N',NSSHP,1,NORBR,1.D0,
     &                          H2D(IOFFP+NOCCP,IOFFR),NORBT,
     &                          DENA(IOFFR,NDIW,ISIM),NORBT,0.D0,
     &                          WRK(1),NORBT)
                     DO 353 IQ = IOFFP+NOCCP,IOFFP+NORBP-1
                        IQ1 = IQ - IOFFP - NOCCP + 1
                        FVOEX(ICI,IQ,ISIM) = FVOEX(ICI,IQ,ISIM)
     *                                      + WRK(IQ1)
 353                 CONTINUE
                  END IF
               END IF
C
 210        CONTINUE
 110     CONTINUE
      END IF
      IF ( NDITR.LT.0) CALL DGETRN(H2D,NORBT,NORBT)
      RETURN
      END
C  /* Deck ftddr */
      SUBROUTINE FTDDR(NSIM,ICI1,IDI1,FVTD,DVT,H2D)
C
C CALCULATE EXCHANGE CONTRIBUTIONS TO ACTIVE FOCK MATRIX WITH A
C TRANSITION DENSITY MATRIX
C
C  FV(P,Q) = FV(P,Q) -0.5*SUM(X,Y)*<XY/QP>*DVT(X,Y)
C THE SUM X,Y ARE TAKEN BY READING IN ALL DIRAC DISTRIBUTIONS
C
#include <implicit.h>
C
C Used from common blocks:
C   INFORB :
C   INFIND : ISMO,IOBTYP
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <orbtypdef.h>
C
      DIMENSION FVTD(NORBT,NORBT,*),DVT(NASHT,NASHT,*)
      PARAMETER ( DMP5 = -0.5D0 )
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      NDITR = 1
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
      ELSE
         ICI = IDI1
         IDI = ICI1
         CALL DGETRN(H2D,NORBT,NORBT)
         NDITR = -NDITR
      END IF
C
C     Find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need occupied-occupied distributions, i.e. itypcd .le. 3
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPCD .NE. 3) RETURN
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      ICDSYM = MULD2H(ICSYM,IDSYM)
      IF ( ICDSYM.EQ.KSYMOP) THEN
         DO 100 ISIM = 1,NSIM
            FAC = DMP5*DVT(NDIW,NCIW,ISIM)
            CALL DAXPY(NORBT*NORBT,FAC,H2D,1,FVTD(1,1,ISIM),1)
 100     CONTINUE
         IF (NCIW.NE.NDIW) THEN
            CALL DGETRN(H2D,NORBT,NORBT)
            NDITR = -NDITR
            DO 200 ISIM = 1,NSIM
               FAC = DMP5*DVT(NCIW,NDIW,ISIM)
               CALL DAXPY(NORBT*NORBT,FAC,H2D,1,FVTD(1,1,ISIM),1)
 200        CONTINUE
         END IF
      END IF
      IF ( NDITR.LT.0) CALL DGETRN(H2D,NORBT,NORBT)
      RETURN
      END
C  /* Deck ftdmu */
      SUBROUTINE FTDMU(NSIM,ICI1,IDI1,FVTD,DVT,H2)
C
C CALCULATE EXCHANGE CONTRIBUTIONS TO ACTIVE FOCK MATRIX WITH A
C TRANSITION DENSITY MATRIX
C
C  FV(P,Q) = FV(P,Q) + SUM(X,Y)*(PQ/XY)*DVT(X,Y)
C THE SUM X,Y ARE TAKEN BY READING IN ALL MULLIKEN DISTRIBUTIONS
C
#include <implicit.h>
C
C Used from common blocks:
C   INFORB :
C   INFIND : ISMO,IOBTYP
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <orbtypdef.h>
C
      DIMENSION FVTD(NORBT,NORBT,*),DVT(NASHT,NASHT,*)
      DIMENSION H2(NORBT,*)
C
      IF (TRPLET) RETURN
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      NDITR = 1
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
      ELSE
         ICI = IDI1
         IDI = ICI1
      END IF
C
C     Find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need occupied-occupied distributions, i.e. itypcd .le. 3
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPCD .NE. 3) RETURN
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      ICDSYM = MULD2H(ICSYM,IDSYM)
      IF ( ICDSYM.EQ.KSYMOP) THEN
         DO 100 ISIM = 1,NSIM
            IF (NCIW.NE.NDIW) THEN
               FAC = DVT(NCIW,NDIW,ISIM)+DVT(NDIW,NCIW,ISIM)
            ELSE
               FAC = DVT(NCIW,NDIW,ISIM)
            ENDIF
            CALL DAXPY(NORBT*NORBT,FAC,H2,1,FVTD(1,1,ISIM),1)
 100     CONTINUE
      END IF
      RETURN
      END
