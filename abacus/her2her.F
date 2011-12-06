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
C  /* Deck hr2drv */
      SUBROUTINE HR2DRV(HERINT,INDHER,COOR12,COOR34,EXP12,EXP34,FAC12,
     &                  FAC34,NTUV,WORK,LWORK,JMAX,MAXDER,NUABCD,IPQ0X,
     &                  IPQ0Y,IPQ0Z,NOINT,ONECEN,NUC1,NUC2,NUC12,NUC3,
     &                  NUC4,NUC34,THRESH,IPRINT,SIGNT,IODDHR)
C
C     T. Helgaker, Sep. 91
C
C     gnrinf.h needed for PANAS correction factor, and CHIVAL,ERFEXP,
C              SRINTS for partition of r_12 (MCDFT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      LOGICAL NOINT, ONECEN
      DIMENSION HERINT(NUABCD,NTUV), INDHER(*), IODDHR(*),
     &          COOR12(*), COOR34(*), EXP12(*), EXP34(*), FAC12(*),
     &          FAC34(*), SIGNT(*), WORK(LWORK)
#include <gnrinf.h>
#include <twosta.h>
      LRJ000 = (JMAX + 1)*NUABCD
      KRJ000 = 1
      KPQX   = KRJ000 + LRJ000
      KPQY   = KPQX  + NUABCD
      KPQZ   = KPQY  + NUABCD
      KLAST  = KPQZ  + NUABCD
      IF (KLAST .GT. LWORK) CALL STOPIT('HR2DRV',' ',KLAST,LWORK)
      MWRJ00 = MAX(MWRJ00,LRJ000)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
      LWRK   = LWORK - KLAST + 1
      CALL HR2DR1(HERINT,INDHER,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,
     &            WORK(KRJ000),WORK(KPQX),WORK(KPQY),WORK(KPQZ),
     &            WORK(KLAST),LWRK,NTUV,JMAX,MAXDER,NUABCD,IPQ0X,IPQ0Y,
     &            IPQ0Z,NOINT,ONECEN,NUC1,NUC2,NUC12,NUC3,NUC4,NUC34,
     &            THRESH,IPRINT,SIGNT,IODDHR,PANAS,CHIVAL,ERFEXP,SRINTS)
      LWTOT  = LWTOT - KLAST
      RETURN
      END
C  /* Deck hr2dr1 */
      SUBROUTINE HR2DR1(HERINT,INDHER,COOR12,COOR34,EXP12,EXP34,FAC12,
     &                  FAC34,RJ000,PQX,PQY,PQZ,WORK,LWORK,NTUV,JMAX,
     &                  MAXDER,NUABCD,IPQ0X,IPQ0Y,IPQ0Z,NOINT,ONECEN,
     &                  NUC1,NUC2,NUC12,NUC3,NUC4,NUC34,THRESH,IPRINT,
     &                  SIGNT,IODDHR,PANAS,CHIVAL,ERFEXP,SRINTS)
C     Modified for r12 integrals (WK/UniKA/15-11-2002).
#include <implicit.h>
#include <priunit.h>
#include <r12int.h>
#include <drw2el.h>
#include <maxaqn.h>
      PARAMETER (D1 = 1.D0)
      LOGICAL NOINT, ONECEN, ERFEXP, SRINTS
#include <twosta.h>
C     HERINT is declared 3-dimensional for r12 integrals (WK/UniKA/14-11-2002).
      DIMENSION HERINT(NUABCD,NTUV,2), INDHER(*), IODDHR(*),
     &          RJ000(NUABCD,*), WORK(LWORK),
     &          PQX(NUABCD), PQY(NUABCD), PQZ(NUABCD),
     &          COOR12(*), COOR34(*), EXP12(*), EXP34(*), FAC12(*),
     &          FAC34(*), SIGNT(*)
C
      IF (R12EIN) THEN
C
C        Gaussian-damped R12 integrals (WK/UniKA/15-11-2002).
C        ====================================================
C
         LRJ000 = NUABCD * (JMAX + 1)
         LHRINT = NUABCD * NTUV
         KFACNT = 1
         KALPHA = KFACNT 
         KHARGE = KALPHA + NUABCD
         KHARGF = KHARGE + NUABCD
         KHALPH = KHARGF + NUABCD
         KHBETA = KHALPH + NUABCD
         KHPRE1 = KHBETA + NUABCD
         KHPRE2 = KHPRE1 + NUABCD
         KHPRE3 = KHPRE2 + NUABCD
         KHEXPP = KHPRE3 + NUABCD
         KHEXPQ = KHEXPP + NUABCD
         KRO000 = KHEXPQ + NUABCD
         KHRIN3 = KRO000 + LRJ000
         KHRIN4 = KHRIN3 + LHRINT
         KHRIN5 = KHRIN4 + LHRINT
         LSTEIN = KHRIN5 + LHRINT
         LWKHER = LWORK - LSTEIN + 1
         IF (LWKHER .LE. 0) CALL STOPIT('HR2DR1',' ',LSTEIN,LWORK)
         IF (TKTIME) TIMSTR = SECOND()
         CALL R00G(RJ000,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,PQX,PQY,
     &             PQZ,JMAX,NOINT,NUABCD,NUC1,NUC2,NUC12,NUC3,NUC4,
     &             NUC34,THRESH,ONECEN,IPRINT,IPQ0X,IPQ0Y,IPQ0Z,SIGNT,
     &             WORK(KFACNT),WORK(KHEXPP),WORK(KHEXPQ))
         CALL ERIIFC(NTUV,NUABCD,IPQ0X,IPQ0Y,IPQ0Z,JMAX)
         CALL ERIHER(HERINT(1,1,1),HERINT(1,1,2),RJ000,PQX,INDHER,
     &              IODDHR,WORK(KFACNT),WORK(KALPHA),WORK(KHARGE),
     &              WORK(KHARGF),WORK(KHALPH),WORK(KHBETA),
     &              WORK(KHPRE1),WORK(KHPRE2),WORK(KHPRE3),
     &              WORK(KRO000),WORK(KHRIN3),WORK(KHRIN4),
     &              WORK(KHRIN5),WORK(KHEXPP),WORK(KHEXPQ),
     &              WORK(LSTEIN),LWKHER,IPRINT)
         IF (TKTIME) THEN
            TIMEND = SECOND()
            TIME = TIMEND - TIMSTR
            THERIX(JMAX) = THERIX(JMAX) + TIME
            THERI = THERI + TIME
         END IF
C
      ELSE  
C
C        Incomplete gamma function
C        =========================
C
         IF (TKTIME) TIMSTR = SECOND()
         CALL R000(RJ000,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,PQX,PQY,
     &             PQZ,JMAX,NOINT,NUABCD,NUC1,NUC2,NUC12,NUC3,NUC4,
     &             NUC34,THRESH,ONECEN,IPRINT,IPQ0X,IPQ0Y,IPQ0Z,SIGNT,
     &             PANAS,CHIVAL,ERFEXP,SRINTS,WORK,LWORK)
         IF (TKTIME) THEN
            TIMEND = SECOND()
            TIME = TIMEND - TIMSTR
            TR000X(JMAX) = TR000X(JMAX) + TIME
            TR000 = TR000 + TIME
            TIMSTR = TIMEND
         END IF
C
C        Hermite integrals
C        =================
C
         IF (.NOT.NOINT) THEN
C
C           Calculate integrals
C
            CALL HERI(HERINT,WORK,LWORK,RJ000,PQX,PQY,PQZ,INDHER,JMAX,
     &             MAXDER,NUABCD,NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
C
            IF (R12INT .OR. U12INT .OR. BPH2OO) THEN
C
C              Calculate Hermite integrals Q(t,u,v) over r12 from
C              R(t,u,v) integrals over 1/r12 (WK/UniKA/14-11-2002).
C
               KWKAMA = 1
               KWKAMB = KWKAMA + NUABCD
               KWKLST = KWKAMB + NUABCD
               IF (KWKLST .GT. LWORK) 
     &         CALL STOPIT('HR2DR1','ABEQ52',KWKLST,LWORK)
               CALL ABEQ52(HERINT(1,1,2),HERINT(1,1,1),WORK(KWKAMA),
     &                     WORK(KWKAMB),PQX,PQY,PQZ,INDHER,JMAX,
     &                     EXP12,EXP34,NUC12,NUC34,NUABCD,
     &                     NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
            END IF
            IF (TKTIME) THEN
               TIMEND = SECOND()
               TIME = TIMEND - TIMSTR
               THERIX(JMAX) = THERIX(JMAX) + TIME
               THERI = THERI + TIME
               TIMSTR = TIMEND
            END IF
         END IF
      END IF
C
      RETURN
      END
C  /* Deck r000 */
      SUBROUTINE R000(RJ000,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,PQX,
     &                PQY,PQZ,JMAX,NOINT,NUABCD,NUC1,NUC2,NUC12,NUC3,
     &                NUC4,NUC34,THRESH,ONECEN,IPRINT,IPQ0X,IPQ0Y,IPQ0Z,
     &                SIGNT,PANAS,CHIVAL,ERFEXP,SRINTS,WORK,LWORK)
C
C     TUH 84
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <twosta.h>
#include <r12int.h>
      LOGICAL ONECEN, NOINT, ERFEXP, SRINTS
      DIMENSION RJ000(*), PQX(*), PQY(*), PQZ(*), COOR12(*), COOR34(*),
     &          EXP12 (*), EXP34(*), FAC12 (*), FAC34(*), SIGNT(3),
     &          WORK(LWORK)
C-----------------------------------------------------------------------
C     =============================================
C     Parameters used for memory-allocations
C     =============================================
      LXJ000 = (JMAX+1)*NUABCD
      IF(ERFEXP) THEN
C      Allocate extra memory for exponential factors 
C      in modified two-electron operator. \jkp
         NTALPHX = 2
      ELSE
         NTALPHX = 1
      ENDIF
      IF (SRINTS) THEN
         NTALPHX = NTALPHX + 1
      END IF
C     =============================================
C     Allocations
C     =============================================
      KWVALU = 1
      KALPHA = KWVALU + NTALPHX*NUABCD
      KALPHJ = KALPHA + NTALPHX*NUABCD
      KINDAD = KALPHJ + NTALPHX*NUABCD
      KEJ000 = KINDAD + (NUABCD + 1)/IRAT
      IF(ERFEXP) THEN
         KSJ000 = KEJ000 + LXJ000
      ELSE
         KSJ000 = KEJ000
      ENDIF
      IF(SRINTS) THEN
         KLAST  = KSJ000 + LXJ000
      ELSE
         KLAST  = KSJ000
      ENDIF
      IF (KLAST .GT. LWORK) CALL STOPIT('R000',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C     =============================================
      CALL R0001(RJ000,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,PQX,PQY,
     &           PQZ,JMAX,NOINT,NUABCD,NUC1,NUC2,NUC12,NUC3,NUC4,NUC34,
     &           THRESH,ONECEN,IPRINT,IPQ0X,IPQ0Y,IPQ0Z,SIGNT,
     &           WORK(KWVALU),WORK(KALPHJ),WORK(KALPHA),WORK(KINDAD),
     &           PANAS,CHIVAL,ERFEXP,NTALPHX,WORK(KEJ000),WORK(KSJ000),
     &           SRINTS,WORK(KLAST),LWRK)
      RETURN
      END
C  /* Deck r0001 */
      SUBROUTINE R0001(RJ000,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,PQX,
     &                PQY,PQZ,JMAX,NOINT,NUABCD,NUC1,NUC2,NUC12,NUC3,
     &                NUC4,NUC34,THRESH,ONECEN,IPRINT,IPQ0X,IPQ0Y,IPQ0Z,
     &                SIGNT,WVALU,TALPHA,SCLFAC,INDADR,PANAS,CHIVAL,
     &                ERFEXP,NTALPHX,EJ000,SJ000,SRINTS,WORK,LWORK)
C
C     TUH 84
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <aovec.h>
#include <codata.h>
#include <drw2el.h>
      PARAMETER (D0=0.0D0,D1=1.0D0,D2=2.0D0,D3=3.0D0,DP25 = 0.25D0)
      LOGICAL ONECEN, NOINT, ERFEXP, SRINTS
      DIMENSION RJ000(NUABCD,0:JMAX),EJ000(NUABCD,0:JMAX),
     &          SJ000(NUABCD,0:JMAX),PQX(NUABCD),PQY(NUABCD),
     &          PQZ(NUABCD),COOR12(NUC1*NUC2,3),COOR34(NUC3*NUC4,3),
     &          EXP12 (*), EXP34(*), FAC12 (*), FAC34(*), SIGNT(3),
     &          WVALU(NUABCD,NTALPHX), TALPHA(NUABCD,NTALPHX), 
     &          INDADR(NUABCD),SCLFAC(NUABCD,NTALPHX),WORK(LWORK)
#include <twosta.h>
#include <dftcom.h>
#include <subdir.h>
C-----------------------------------------------------------------------
C
C     *****************
C     ***** RJ000 *****
C     *****************
C
C
C     Special case: One-center Integrals
C     ==================================
C
C     Note: There should be no testing for small integrals since
C     this may in the case of one-center integrals introduce
C     numerical instabilities for large exponents.
C
      NOINT = .FALSE.
      IF (ERFEXP) IEJ = 2
      IF (SRINTS.AND.ERFEXP) IEJ = 3
      IF (ERFEXP) CALL DZERO(EJ000(1,0),(JMAX+1)*NUABCD)
      IF (AD2DAR) DRWFAC=DP25*ALPHA2*DARFAC
      IF (DO2DAR) DRWFAC=DP25*ALPHA2
      IF (S4CENT) DRWFAC=-DP25/PI
      IF (ONECEN) THEN
         IODS = 0
         CALL DZERO(PQX,NUABCD)
         CALL DZERO(PQY,NUABCD)
         CALL DZERO(PQZ,NUABCD)
         IPQ0X = 1
         IPQ0Y = 1
         IPQ0Z = 1
         DO 100 IOD12 = 1, NUC12
            EXPP   = EXP12(IOD12)
            FAC12I = FAC12(IOD12)
            IF (DO2DAR .OR. S4CENT) THEN
C
C              Compute two-electron Dirac delta function
C              (Wim Klopper, University of Utrecht, March 4, 1999)
C              Merged into Dalton1.1 by A. Halkier 29/11/99.
C
               DO IOD34 = 1, NUC34
                  IODS   = IODS + 1
                  EXPQ   = EXP34(IOD34)
                  EXPPQ  = EXPP + EXPQ
                  EXPPQI = D1/EXPPQ 
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQI)
                  TALPHA(IODS,1) = - D2*EXPP*EXPQ*EXPPQI
                  SCLFAC(IODS,1) = FACTOR*TALPHA(IODS,1)
                  RJ000(IODS,0) = DRWFAC*SCLFAC(IODS,1)
               END DO
            ELSE IF (AD2DAR) THEN
C
C              Add two-electron Darwin integrals with perturbation
C              parameter DARFAC onto repulsion integrals.
C              (Wim Klopper, University of Utrecht, March 4, 1999)
C              Merged into Dalton1.1 by A. Halkier 29/11/99.
C
               DO IOD34 = 1, NUC34
                  IODS   = IODS + 1
                  EXPQ   = EXP34(IOD34)
                  EXPPQ  = EXPP + EXPQ
                  EXPPQI = D1/EXPPQ 
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQI)
                  TALPHA(IODS,1) = - D2*EXPP*EXPQ*EXPPQI
                  SCLFAC(IODS,1) = FACTOR*TALPHA(IODS,1)
                  RJ000(IODS,0) = FACTOR+DRWFAC*SCLFAC(IODS,1)
               END DO
            ELSE IF (HFXMU .NE. D0) THEN
               DO IOD34 = 1, NUC34
                  IODS   = IODS + 1
                  EXPQ   = EXP34(IOD34)
                  EXPPQ  = EXPP + EXPQ
                  EXPPQI = D1/EXPPQ 
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQI)
                  ALPHA  = EXPP*EXPQ*EXPPQI
                  BETA   = HFXMU**2/(ALPHA + HFXMU**2)
                  TALPHA(IODS,1) = - D2*ALPHA*BETA
                  RJ000(IODS,0) = SQRT(BETA)*FACTOR
               END DO
            ELSE IF (PANAS .EQ. D0 .AND. CHIVAL .EQ. -1.D0) THEN
C
C              Regular two-electron integrals
C
               DO IOD34 = 1, NUC34
                  IODS   = IODS + 1
                  EXPQ   = EXP34(IOD34)
                  EXPPQ  = EXPP + EXPQ
                  EXPPQI = D1/EXPPQ
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQI)
                  TALPHA(IODS,1) = - D2*EXPP*EXPQ*EXPPQI
                  RJ000(IODS,0) = FACTOR
               END DO
            ELSE IF (SRINTS) THEN
C
C              Short-range two-electron integrals (MCDFT)
C
               CHISQ     = CHIVAL**2
               DO IOD34 = 1, NUC34
                  IODS   = IODS + 1
c              --- Regular
                  EXPQ   = EXP34(IOD34)
                  EXPPQ  = EXPP + EXPQ
                  EXPPQI = D1/EXPPQ
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQI)
                  TALPHA(IODS,1) = - D2*EXPP*EXPQ*EXPPQI
                  RJ000(IODS,0) = FACTOR
c              --- Modified
                  EXPPI  = D1/EXPP
                  EXPQI  = D1/EXPQ
                  EXPPQ  = EXPPI*EXPQI
                  EXPPQI = EXPPI + EXPQI + D1/CHISQ
                  EXPPQII= D1/EXPPQI
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQ*EXPPQII)
                  TALPHA(IODS,2) = - D2*EXPPQII
                  SJ000(IODS,0) = FACTOR
                  IF(ERFEXP) THEN
                     AI = D3/CHISQ
                     EXPFAC  = EXPPI + EXPQI + AI
                     EXPFACI = D1/EXPFAC
                     TALPHA(IODS,IEJ) = - D2*EXPFACI
                     EJ000(IODS,0) = - FAC12I*FAC34(IOD34)*CHIVAL*
     &               SQRT((AI*EXPFACI)**3/(EXPP*EXPQ))
                  ENDIF
c              ---
               END DO
            ELSE IF (PANAS .EQ. D0) THEN
C
C              Long-range two-electron integrals (MCDFT)
C
               CHISQ     = CHIVAL**2
               DO IOD34 = 1, NUC34
                  IODS   = IODS + 1
                  EXPQ   = EXP34(IOD34)
                  EXPPI  = D1/EXPP
                  EXPQI  = D1/EXPQ
                  EXPPQ  = EXPPI*EXPQI
                  EXPPQI = EXPPI + EXPQI + D1/CHISQ
                  EXPPQII= D1/EXPPQI
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQ*EXPPQII)
                  TALPHA(IODS,1) = - D2*EXPPQII
                  RJ000(IODS,0) = FACTOR
                  IF(ERFEXP) THEN
                     AI = D3/CHISQ
                     EXPFAC  = EXPPI + EXPQI + AI
                     EXPFACI = D1/EXPFAC
                     TALPHA(IODS,IEJ) = - D2*EXPFACI
                     EJ000(IODS,0) = - FAC12I*FAC34(IOD34)*CHIVAL*
     &               SQRT((AI*EXPFACI)**3/(EXPP*EXPQ))
                  ENDIF
               END DO
            ELSE
C
C              Two-electron integrals regularized according to I.Panas
C
               DO IOD34 = 1, NUC34
                  IODS   = IODS + 1
                  EXPQ   = EXP34(IOD34)
                  EXPPQ  = EXPP + EXPQ
                  EXPMOD = EXPPQ/D2
                  EXPPQP = EXPMOD + SQRT(EXPMOD**2 + EXPP*EXPQ*PANAS)
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(D1/EXPPQP)
                  TALPHA(IODS,1) = - D2*EXPP*EXPQ/EXPPQP
                  RJ000(IODS,0) = FACTOR
               END DO
            END IF
  100    CONTINUE
         IF (JMAX .GT. 0) THEN
            IF (DO2DAR .OR. S4CENT) THEN
               DO J = 1, JMAX
                  DO I = 1, NUABCD
                     SCLFAC(I,1) = TALPHA(I,1)*SCLFAC(I,1)
                     RJ000(I,J) = DRWFAC*SCLFAC(I,1)
                  END DO
               END DO
            ELSE IF (AD2DAR) THEN
               DO J = 1, JMAX
                  FAC = D1/(2.0D0*J + 1.0D0)
                  DO I = 1, NUABCD
                     RJ000(I,J) = (FAC+DRWFAC*TALPHA(I,1))*SCLFAC(I,1)
                     SCLFAC(I,1) = TALPHA(I,1)*SCLFAC(I,1)
                  END DO
               END DO
            ELSE
               CALL DCOPY(NUABCD,TALPHA(1,1),1,SCLFAC(1,1),1)
               IF(SRINTS) CALL DCOPY(NUABCD,TALPHA(1,2),1,SCLFAC(1,2),1)
               IF(ERFEXP) 
     &                CALL DCOPY(NUABCD,TALPHA(1,IEJ),1,SCLFAC(1,IEJ),1)
               IF (SRINTS) THEN
C             ... Calculate SJ000 for short-range integrals
                  DO J = 1, JMAX
                     FAC = D1/(2.0D0*J + 1.0D0)
                     DO I = 1, NUABCD
                        SJ000(I,J) = FAC*SCLFAC(I,2)*SJ000(I,0)
                        SCLFAC(I,2) = TALPHA(I,2)*SCLFAC(I,2)
                     END DO
                  END DO
               END IF
               IF (ERFEXP) THEN
C             ... Calculate EJ000 for long-range integrals
                  DO J = 1, JMAX
                     DO I = 1, NUABCD
                        EJ000(I,J)  = EJ000(I,0)*SCLFAC(I,IEJ)
                        SCLFAC(I,IEJ) = TALPHA(I,IEJ)*SCLFAC(I,IEJ)
                     END DO
                  END DO
               END IF
C             ... Calculate regular two-electron integrals
               DO J = 1, JMAX
                  FAC = D1/(2.0D0*J + 1.0D0)
                  DO I = 1, NUABCD
                     RJ000(I,J) = FAC*SCLFAC(I,1)*RJ000(I,0)
                     SCLFAC(I,1) = TALPHA(I,1)*SCLFAC(I,1)
                  END DO
               END DO
            END IF
         END IF   ! JMAX .GT. 0
C
C     General case: Multicenter Integrals
C     ===================================
C
      ELSE
         IF (.NOT.DPATH1) THEN
            SGN12X = - SIGNT(1)
            SGN12Y = - SIGNT(2)
            SGN12Z = - SIGNT(3)
            SGN34X = - D1
            SGN34Y = - D1
            SGN34Z = - D1
         ELSE
            SGN12X = D1
            SGN12Y = D1
            SGN12Z = D1
            SGN34X = SIGNT(1)
            SGN34Y = SIGNT(2)
            SGN34Z = SIGNT(3)
         END IF
C
C    It turns out 64-bits IRIX compiler misinterprets this statement with
C       -O2 optimization. These fixes prevent it. K.Ruud-Jan96
#if defined (SYS_IRIX)
         IODS  = 1
         NODS  = 1
#else
         IODS  = 0
         NODS  = 0
#endif
         DO 300 IOD12 = 1, NUC12
            EXPP   = EXP12(IOD12)
            FAC12I = FAC12(IOD12)
            PX     = SGN12X*COOR12(IOD12,1)
            PY     = SGN12Y*COOR12(IOD12,2)
            PZ     = SGN12Z*COOR12(IOD12,3)
            IF (HFXMU .NE. D0) THEN
               DO IOD34 = 1, NUC34
#if !defined (SYS_IRIX)
                  IODS = IODS + 1
#endif
                  EXPQ   = EXP34(IOD34)
                  EXPPQI = D1/(EXPP + EXPQ)
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQI)
                  ALPHA = EXPP*EXPQ*EXPPQI
                  BETA = HFXMU**2/(ALPHA + HFXMU**2)
                  SCLFAC(IODS,1) = SQRT(BETA)*FACTOR
#if !defined (SYS_IRIX)
                  NODS = NODS + 1
#endif
                  TALPHA(IODS,1) = - D2*ALPHA*BETA
                  PQXI = PX - SGN34X*COOR34(IOD34,1)
                  PQYI = PY - SGN34Y*COOR34(IOD34,2)
                  PQZI = PZ - SGN34Z*COOR34(IOD34,3)
                  PQX(IODS) = PQXI
                  PQY(IODS) = PQYI
                  PQZ(IODS) = PQZI
                  WVALU(NODS,1)=ALPHA*BETA*
     &               (PQXI*PQXI+PQYI*PQYI+PQZI*PQZI)
                  INDADR(NODS) = IODS
#if defined (SYS_IRIX)
                  IODS = IODS + 1
                  NODS = NODS + 1
#endif
               END DO
            ELSE IF (PANAS .EQ. D0 .AND. CHIVAL .EQ. -1.D0) THEN
C
C              Regular two-electron integrals
C
               DO IOD34 = 1, NUC34
#if !defined (SYS_IRIX)
                  IODS = IODS + 1
#endif
                  EXPQ   = EXP34(IOD34)
                  EXPPQI = D1/(EXPP + EXPQ)
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQI)
                  SCLFAC(IODS,1) = FACTOR
#if !defined (SYS_IRIX)
                  NODS = NODS + 1
#endif
                  ALPHA = EXPP*EXPQ*EXPPQI
                  TALPHA(IODS,1) = - D2*ALPHA
                  PQXI = PX - SGN34X*COOR34(IOD34,1)
                  PQYI = PY - SGN34Y*COOR34(IOD34,2)
                  PQZI = PZ - SGN34Z*COOR34(IOD34,3)
                  PQX(IODS) = PQXI
                  PQY(IODS) = PQYI
                  PQZ(IODS) = PQZI
                  WVALU(NODS,1) = ALPHA*(PQXI*PQXI+PQYI*PQYI+PQZI*PQZI)
                  INDADR(NODS) = IODS
#if defined (SYS_IRIX)
                  IODS = IODS + 1
                  NODS = NODS + 1
#endif
               END DO
            ELSE IF (SRINTS) THEN
C
C              Short-range two-electron integrals (MCDFT)
C
               CHISQ     = CHIVAL**2
               DO IOD34 = 1, NUC34
#if !defined (SYS_IRIX)
                  IODS = IODS + 1
#endif
c                --- Regular
                  EXPQ   = EXP34(IOD34)
                  EXPPQI = D1/(EXPP + EXPQ)
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQI)
                  SCLFAC(IODS,1) = FACTOR
c                --- Modified
                  EXPPI  = D1/EXPP
                  EXPQI  = D1/EXPQ
                  EXPPQ  = EXPPI*EXPQI
                  EXPPQI2= EXPPI + EXPQI + D1/CHISQ
                  EXPPQII= D1/EXPPQI2
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQ*EXPPQII)
                  SCLFAC(IODS,2) = FACTOR
#if !defined (SYS_IRIX)
                  NODS = NODS + 1
#endif
                  ALPHA = EXPP*EXPQ*EXPPQI
                  TALPHA(IODS,1) = - D2*ALPHA
                  TALPHA(IODS,2) = - D2*EXPPQII
                  PQXI = PX - SGN34X*COOR34(IOD34,1)
                  PQYI = PY - SGN34Y*COOR34(IOD34,2)
                  PQZI = PZ - SGN34Z*COOR34(IOD34,3)
                  PQX(IODS) = PQXI
                  PQY(IODS) = PQYI
                  PQZ(IODS) = PQZI
                  PQXYZI    = PQXI*PQXI+PQYI*PQYI+PQZI*PQZI
                  IF(ERFEXP) THEN
                     AI = D3/CHISQ
                     EXPFAC  = EXPPI + EXPQI + AI
                     EXPFACI = D1/EXPFAC
                     FACTOR = FAC12I*FAC34(IOD34)*
     &               CHIVAL*SQRT(D1/(EXPP*EXPQ))*(SQRT(AI*EXPFACI))**3
                     EJ000(IODS,0) = - FACTOR * DEXP(- EXPFACI*PQXYZI)
                     TALPHA(IODS,IEJ) = - D2*EXPFACI
                  ENDIF
                  WVALU(NODS,1) = ALPHA*PQXYZI
                  WVALU(NODS,2) = EXPPQII*PQXYZI
                  INDADR(NODS) = IODS
#if defined (SYS_IRIX)
                  IODS = IODS + 1
                  NODS = NODS + 1
#endif
               END DO
            ELSE IF (PANAS .EQ. D0) THEN
C
C              Long-range two-electron integrals (MCSCF-DFT)
C
               CHISQ     = CHIVAL**2
               DO IOD34 = 1, NUC34
#if !defined (SYS_IRIX)
                  IODS = IODS + 1
#endif
                  EXPQ   = EXP34(IOD34)
                  EXPPI  = D1/EXPP
                  EXPQI  = D1/EXPQ
                  EXPPQ  = EXPPI*EXPQI
                  EXPPQI = EXPPI + EXPQI + D1/CHISQ
                  EXPPQII= D1/EXPPQI
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(EXPPQ*EXPPQII)
                  SCLFAC(IODS,1) = FACTOR
#if !defined (SYS_IRIX)
                  NODS = NODS + 1
#endif
                  TALPHA(IODS,1) = - D2*EXPPQII
                  PQXI = PX - SGN34X*COOR34(IOD34,1)
                  PQYI = PY - SGN34Y*COOR34(IOD34,2)
                  PQZI = PZ - SGN34Z*COOR34(IOD34,3)
                  PQX(IODS) = PQXI
                  PQY(IODS) = PQYI
                  PQZ(IODS) = PQZI
                  PQXYZI    = PQXI*PQXI+PQYI*PQYI+PQZI*PQZI
                  WVALU(NODS,1) = EXPPQII*PQXYZI
                  IF (ERFEXP) THEN
                     AI = D3/CHISQ
                     EXPFAC  = EXPPI + EXPQI + AI
                     EXPFACI = D1/EXPFAC
                     TALPHA(IODS,IEJ) = - D2*EXPFACI
                     FACTOR = FAC12I*FAC34(IOD34)*CHIVAL*
     &                        SQRT(D1/(EXPP*EXPQ))*(SQRT(AI*EXPFACI))**3
                     EJ000(IODS,0) = - FACTOR* DEXP(- EXPFACI*PQXYZI)
                  ENDIF
                  INDADR(NODS) = IODS
#if defined (SYS_IRIX)
                  IODS = IODS + 1
                  NODS = NODS + 1
#endif
               END DO
            ELSE
C
C              Two-electron integrals regularized according to I.Panas
C
               DO IOD34 = 1, NUC34
#if !defined (SYS_IRIX)
                  IODS = IODS + 1
#endif
                  EXPQ   = EXP34(IOD34)
                  EXPPQ  = EXPP + EXPQ
                  EXPMOD = EXPPQ/D2
                  EXPPQP = EXPMOD + SQRT(EXPMOD**2 + EXPP*EXPQ*PANAS)
                  FACTOR = FAC12I*FAC34(IOD34)*SQRT(D1/EXPPQP)
                  SCLFAC(IODS,1) = FACTOR
#if !defined (SYS_IRIX)
                  NODS = NODS + 1
#endif
                  ALPHA = EXPP*EXPQ/EXPPQP
                  TALPHA(IODS,1) = - D2*ALPHA
                  PQXI = PX - SGN34X*COOR34(IOD34,1)
                  PQYI = PY - SGN34Y*COOR34(IOD34,2)
                  PQZI = PZ - SGN34Z*COOR34(IOD34,3)
                  PQX(IODS) = PQXI
                  PQY(IODS) = PQYI
                  PQZ(IODS) = PQZI
                  WVALU(NODS,1) = ALPHA*(PQXI*PQXI+PQYI*PQYI+PQZI*PQZI)
                  INDADR(NODS) = IODS
#if defined (SYS_IRIX)
                  IODS = IODS + 1
                  NODS = NODS + 1
#endif
               END DO
            END IF
 300     CONTINUE
#if defined (SYS_IRIX)
         IODS = IODS - 1
         NODS = NODS - 1
#endif
         IF (.NOT.NOINT) THEN
            IPQ0X = 1
            IPQ0Y = 1
            IPQ0Z = 1
            IF (DASUM(NUABCD,PQX,1) .GT. THRESH) IPQ0X = 0
            IF (DASUM(NUABCD,PQY,1) .GT. THRESH) IPQ0Y = 0
            IF (DASUM(NUABCD,PQZ,1) .GT. THRESH) IPQ0Z = 0
C
            CALL DZERO(RJ000(1,0),(JMAX+1)*NUABCD)
            IF (SRINTS) CALL DZERO(SJ000(1,0),(JMAX+1)*NUABCD)
C
C           Calculate gamma function
C           ========================
C
C           Allocations
C
            KINDAD = 1
            KWVALS = KINDAD + (3*NODS + 1)/IRAT
            KFJW   = KWVALS + 3*NODS
            KREXPW = KFJW   + NODS*(JMAX + 1)
            KLAST  = KREXPW + NODS
            IF (KLAST .GT. LWORK) CALL STOPIT('R0001',' ',KLAST,LWORK)
            IF (.NOT. (DO2DAR .OR. S4CENT)) THEN
               CALL GETGAM(NODS,INDADR,WVALU(1,1),RJ000,JMAX,NUABCD,
     &                     WORK(KFJW),WORK(KINDAD),WORK(KWVALS),
     &                     WORK(KREXPW),IPRINT)
               IF (SRINTS)
     &            CALL GETGAM(NODS,INDADR,WVALU(1,2),SJ000,JMAX,NUABCD,
     &                        WORK(KFJW),WORK(KINDAD),WORK(KWVALS),
     &                        WORK(KREXPW),IPRINT)
               IF (ERFEXP)
     &            CALL DCOPY(NUABCD,TALPHA(1,IEJ),1,SCLFAC(1,IEJ),1)
            END IF
C
C           Scale gamma function
C
            IF (DO2DAR .OR. AD2DAR .OR. S4CENT) THEN
              DO I = 1, NUABCD
                WVALU(I,1) = DRWFAC*TALPHA(I,1)*EXP(-WVALU(I,1))
              END DO
             IF (DO2DAR .OR. S4CENT) THEN
              DO J = 0, JMAX
                DO I = 1, NUABCD
                  RJ000(I,J) = WVALU(I,1)
                END DO
              END DO 
             ELSE IF (AD2DAR) THEN
              DO J = 0, JMAX
                DO I = 1, NUABCD
                  RJ000(I,J) = RJ000(I,J) + WVALU(I,1)
                END DO
              END DO
             END IF
            ENDIF
            IF (SRINTS) THEN
C          ... Calculate SJ000 for short-range integrals
               DO J = 0, JMAX
                  DO I = 1, NUABCD
                     SJ000(I,J) = SCLFAC(I,2)*SJ000(I,J)
                     SCLFAC(I,2) = TALPHA(I,2)*SCLFAC(I,2)
                  END DO
               END DO
            END IF
            IF (ERFEXP) THEN
C          ... Calculate EJ000 for long-range integrals
               DO J = 1, JMAX
                  DO I = 1, NUABCD
                  EJ000(I,J)  = EJ000(I,0)*SCLFAC(I,IEJ)
                  SCLFAC(I,IEJ) = TALPHA(I,IEJ)*SCLFAC(I,IEJ)
                  END DO
               END DO
            END IF
C          ... Calculate regular two-electron integrals
            DO J = 0, JMAX
               DO I = 1, NUABCD
                  RJ000(I,J) = SCLFAC(I,1)*RJ000(I,J)
                  SCLFAC(I,1) = TALPHA(I,1)*SCLFAC(I,1)
               END DO
            END DO
         END IF
C
      END IF
C
C     For modified operators, R000(I,J) can now be assigned correctly
C
      IF(SRINTS) THEN
         IF (ERFEXP) 
     &      CALL DAXPY(NUABCD*(JMAX+1), D1,EJ000(1,0),1,SJ000(1,0),1)
            CALL DAXPY(NUABCD*(JMAX+1),-D1,SJ000(1,0),1,RJ000(1,0),1)
      ELSE
         IF (ERFEXP)
     &      CALL DAXPY(NUABCD*(JMAX+1),D1,EJ000(1,0),1,RJ000(1,0),1)
      ENDIF
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1010) JMAX
         WRITE (LUPRI, 1020) NUC12, NUC34
         WRITE (LUPRI, 1030) THRESH
         WRITE (LUPRI, 1040) ONECEN
         WRITE (LUPRI, 1045) NOINT
         WRITE (LUPRI, 1050) NUABCD
         IF (IPRINT .GT. 20) THEN
            WRITE (LUPRI, 1100)
            ISTART = 0
            DO J = 0, JMAX
               WRITE (LUPRI, 1110) J
               IADR = 0
               DO I = 1, NUC12
                  WRITE (LUPRI, 1120) I
                  WRITE (LUPRI, 1130) (RJ000(IADR + K,J), K = 1, NUC34)
                  IADR = IADR + NUC34
               END DO
            END DO
         END IF
      END IF
 1000 FORMAT (//,' ********** SUBROUTINE R000 **********')
 1010 FORMAT (//,'  JMAX    ',I7)
 1020 FORMAT (   '  NUC     ',2I7)
 1030 FORMAT (   '  THRESH  ',1P,D12.4)
 1040 FORMAT (   '  ONECEN  ',L7)
 1045 FORMAT (   '  NOINT   ',L7)
 1050 FORMAT (   '  NUABCD  ',I7)
 1100 FORMAT(//,' ***** RJ000 - INTEGRALS *****')
 1110 FORMAT( /,'  J     ',I7)
 1120 FORMAT( /,'  NUC12 ',I7,/)
 1130 FORMAT(1P,6D12.4)
      RETURN
      END
C  /* Deck getgam */
      SUBROUTINE GETGAM(NODS,INDADR,WVALU,RJ000,JMAX,NUABCD,FJWS,INDADS,
     &                  WVALS,REXPW,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <pi.h>
      PARAMETER (HALF = 0.5D0)
      PARAMETER (D1 = 1.D0, D10 = 10.D0)
      PARAMETER (D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.1D0)
      PARAMETER (COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0,
     &           COEF5 = - D1/120.D0, COEF6 = D1/720.D0)
      PARAMETER (GFAC0 =  D2*0.4999489092 D0,
     &           GFAC1 = -D2*0.2473631686 D0,
     &           GFAC2 =  D2*0.321180909  D0,
     &           GFAC3 = -D2*0.3811559346 D0)
      PARAMETER (SQRPIH = SQRTPI/D2)
      PARAMETER (PID4 = PI/D4, PID4I = D4/PI)
#include <gamcom.h>
      DIMENSION FJWS(NODS,0:JMAX), INDADR(*),
     &          WVALU(*), RJ000(NUABCD,0:JMAX), INDADS(NODS,3),
     &          WVALS(NODS,3), REXPW(NODS)
C-----------------------------------------------------------------------
C
      NODS1 = 0
      NODS2 = 0
      NODS3 = 0
      D2JP36 = 2.0D0*JMAX + 36.0D0
C NECgh980808 As sth suggested, this loop DO 100 can savely 
C             be vectorized.
#include <ivdep.h>
      DO 100 I = 1, NODS
         WVAL = WVALU(I)
         IF (WVAL .LT. D12) THEN
            NODS1 = NODS1 + 1
            INDADS(NODS1,1) = INDADR(I)
            WVALS(NODS1,1)  = WVAL
         ELSE IF (WVAL .LE. D2JP36) THEN
            NODS2 = NODS2 + 1
            INDADS(NODS2,2) = INDADR(I)
            WVALS(NODS2,2)  = WVAL
         ELSE
            NODS3 = NODS3 + 1
            INDADS(NODS3,3) = INDADR(I)
            WVALS(NODS3,3)  = WVAL
         END IF
  100 CONTINUE
C
C     WVAL < 12
C
      IF (NODS1 .GT. 0) THEN
         ISTRT0 = 1 + 121*JMAX
         DO 200 I = 1, NODS1
            WVAL = WVALS(I,1)
            IPNT = NINT(D10*WVAL)
            WDIF = WVAL - TENTH*IPNT
            ISTART = ISTRT0 + IPNT
            FJWS(I,JMAX) = (((((COEF6*TABFJW(ISTART + 726)*WDIF
     &                         + COEF5*TABFJW(ISTART + 605))*WDIF
     &                         + COEF4*TABFJW(ISTART + 484))*WDIF
     &                         + COEF3*TABFJW(ISTART + 363))*WDIF
     &                         + COEF2*TABFJW(ISTART + 242))*WDIF
     &                         - TABFJW(ISTART + 121))*WDIF
     &                         + TABFJW(ISTART)
  200    CONTINUE
         IF (JMAX .GT. 0) THEN
            DO 300 I = 1, NODS1
               REXPW(I) = HALF*EXP(-WVALS(I,1))
  300       CONTINUE
            DO 310 J = JMAX - 1, 0, -1
               FCT = D2/(2.0D0*J + 1.0D0)
               DO 320 I = 1, NODS1
                  FJWS(I,J) = FCT*(WVALS(I,1)*FJWS(I,J+1) + REXPW(I))
  320          CONTINUE
  310       CONTINUE
            DO 330 J = 1, JMAX
               DO 340 I = 1, NODS1
                  RJ000(INDADS(I,1),J) = FJWS(I,J)
  340          CONTINUE
  330       CONTINUE
         END IF
         DO 350 I = 1, NODS1
            RJ000(INDADS(I,1),0) = FJWS(I,0)
  350    CONTINUE
      END IF
C
C     Near asymptotic region
C
      IF (NODS2 .GT. 0) THEN
         DO 400 I = 1, NODS2
            WVAL       = WVALS(I,2)
            REXPW(I)   = HALF*EXP(-WVAL)
            WVALS(I,2) = D1/WVAL
  400    CONTINUE
         DO 410 I = 1, NODS2
            RWVAL = WVALS(I,2)
            GVAL = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
            FJWS(I,0) = SQRPIH*SQRT(RWVAL) - REXPW(I)*GVAL*RWVAL
  410    CONTINUE
         DO 420 I = 1, NODS2
            RJ000(INDADS(I,2),0) = FJWS(I,0)
  420    CONTINUE
         DO 430 J = 1, JMAX
            FCT = J - HALF
            DO 440 I = 1, NODS2
               FJWS(I,J) = (FCT*FJWS(I,J-1) - REXPW(I))*WVALS(I,2)
  440       CONTINUE
            DO 450 I = 1, NODS2
               RJ000(INDADS(I,2),J) = FJWS(I,J)
  450       CONTINUE
  430    CONTINUE
      END IF
C
C     Asymptotic region
C
      IF (NODS3 .GT. 0) THEN
         DO 500 I = 1, NODS3
            WVALS(I,3) = PID4/WVALS(I,3)
  500    CONTINUE
         DO 510 I = 1, NODS3
            FJWS(I,0) = SQRT(WVALS(I,3))
  510    CONTINUE
         DO 520 I = 1, NODS3
            RJ000(INDADS(I,3),0) = FJWS(I,0)
  520    CONTINUE
         DO 530 J = 1, JMAX
            FACTOR = PID4I*(J - HALF)
            DO 540 I = 1, NODS3
               FJWS(I,J) = FACTOR*FJWS(I,J-1)*WVALS(I,3)
  540       CONTINUE
            DO 550 I = 1, NODS3
               RJ000(INDADS(I,3),J) = FJWS(I,J)
  550       CONTINUE
  530    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck heri */
      SUBROUTINE HERI(HERINT,WORK,LWORK,RJ000,PQX,PQY,PQZ,INDHER,JMAX,
     &                MAXDER,NUABCD,NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,
     &                IPRINT)
C
C     tuh fall 1984
C
C     Modified Jul 28 88 to avoid multiplying zero with
C     undetermined numbers - tuh
C
C     Modified for triangular addressing March 92 - tuh
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      INTEGER T, U, V, TUV
      DIMENSION WORK(LWORK), RJ000(NUABCD,0:JMAX), IODDHR(*),
     &          PQX(NUABCD), PQY(NUABCD), PQZ(NUABCD),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), HERINT(NUABCD,NTUV)
#include <twosta.h>
#include <hertop.h>
C
C     The final R(T,U,V) integrals are arranged as follows:
C
C     R(000)
C     R(100) R(010) R(001)
C     R(200) R(110) R(101) R(020) R(011) R(002)
C     R(300) R(210) R(201) R(120) R(111) R(102) R(030) R(021) R(012)
C                                                             R(003)
C     Special case JMAX = 0
C     =====================
C
      IF (JMAX .EQ. 0) THEN
         CALL DCOPY(NUABCD,RJ000,1,HERINT,1)
      ELSE
C
C        Allocate work space
C        ===================
C
         KHRWRK = 1
         KLAST  = KHRWRK + NTUV*NUABCD
         IF (KLAST .GT. LWORK) CALL STOPIT('HERI',' ',KLAST,LWORK)
         MWHRSQ = MAX(MWHRSQ,NTUV*NUABCD)
         LWTOT  = LWTOT + KLAST
         MWTOT  = MAX(MWTOT,LWTOT)
C
C        Recursion loop for Hermite integrals
C        ====================================
C
         IPQX = IPQ0X + 1
         IPQY = IPQ0Y + 1
         IPQZ = IPQ0Z + 1
         CALL DZERO(HERINT,NUABCD*NTUV)
         DO 200 JVAL = 1, JMAX
            IF (MOD(JMAX-JVAL,2).EQ.0) THEN
               CALL HRECUR(HERINT,WORK(KHRWRK),JVAL,RJ000,PQX,PQY,PQZ,
     &                     INDHER,JMAX,MAXDER,NUABCD,NTUV,IPQX,IPQY,
     &                     IPQZ)
            ELSE
               CALL HRECUR(WORK(KHRWRK),HERINT,JVAL,RJ000,PQX,PQY,PQZ,
     &                     INDHER,JMAX,MAXDER,NUABCD,NTUV,IPQX,IPQY,
     &                     IPQZ)
            END IF
  200    CONTINUE
         LWTOT  = LWTOT - KLAST
      END IF
C
C     Print section
C     =============
C
      IF (IPRINT .GE. 10) THEN
         CALL TITLER('Output from HERI','*',103)
         WRITE (LUPRI,'(2X,A,I5)') 'JMAX  ', JMAX
         WRITE (LUPRI,'(2X,A,I5)') 'NUABCD', NUABCD
         WRITE (LUPRI,'(2X,A,I5)') 'NTUV  ', NTUV
         IF (IPRINT .GE. 20) THEN
            CALL HEADER('Hermite integrals R(t,u,v)',1)
            DO 300 J = 0, JMAX
              DO 320 T = J, 0, -1
                DO 330 U = J - T, 0, -1
                  V = J - T - U
                  TUV = INDHER(T,U,V)
                  IF (IODDHR(TUV) .EQ. 0) THEN
                    WRITE (LUPRI,'(2X,3(A,I1),A,2X,5F12.8/,
     &                                                 (12X,5F12.8))')
     &              'R(',T,',',U,',',V,')', (HERINT(I,TUV),I=1,NUABCD)
                  WRITE (LUPRI,'()')
                  END IF
  330           CONTINUE
  320         CONTINUE
  300      CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck hrecur */
      SUBROUTINE HRECUR(CUR,OLD,JVAL,RJ000,PQX,PQY,PQZ,INDHER,JMAX,
     &                  MAXDER,NUABCD,NTUV,IPQX,IPQY,IPQZ)
C
#include <implicit.h>
      INTEGER T, U, V, TUV
      LOGICAL PQXGT0, PQYGT0, PQZGT0
      DIMENSION CUR(NUABCD,NTUV), OLD(NUABCD,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          PQX(NUABCD), PQY(NUABCD), PQZ(NUABCD),
     &          RJ000(NUABCD,0:JMAX)
#include <doxyz.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      PQXGT0 = IPQX .EQ. 1
      PQYGT0 = IPQY .EQ. 1
      PQZGT0 = IPQZ .EQ. 1
C
C     JVAL = 1
C     ========
C
      IF (JVAL .EQ. 1) THEN
         CALL DCOPY(NUABCD,RJ000(1,JMAX-1),1,CUR(1,1),1)
         IF (PQXGT0) THEN
            DO 110 I = 1, NUABCD
               CUR(I,2) = PQX(I)*RJ000(I,JMAX)
  110       CONTINUE
         END IF
         IF (PQYGT0) THEN
            DO 120 I = 1, NUABCD
               CUR(I,3) = PQY(I)*RJ000(I,JMAX)
  120       CONTINUE
         END IF
         IF (PQZGT0) THEN
            DO 130 I = 1, NUABCD
               CUR(I,4) = PQZ(I)*RJ000(I,JMAX)
  130       CONTINUE
         END IF
C
C     JVAL > 1
C     ========
C
      ELSE
         MAXT   = JMAX
         MAXU   = JMAX
         MAXV   = JMAX
         IF (.NOT.DOX) MAXT = JMAX - MAXDER
         IF (.NOT.DOY) MAXU = JMAX - MAXDER
         IF (.NOT.DOZ) MAXV = JMAX - MAXDER
C
C        R(0,0,0)
C
         CALL DCOPY(NUABCD,RJ000(1,JMAX-JVAL),1,CUR,1)
C
C        R(T,0,0)
C
         IF (PQXGT0) THEN
            DO 200 I = 1, NUABCD
               CUR(I,2) = PQX(I)*OLD(I,1)
  200       CONTINUE
            DO 300 T = 2, MIN(MAXT,JVAL)
               TMIN1 = T - 1.0D0
               TUV   = INDHER(T  ,0,0)
               M1T   = INDHER(T-1,0,0)
               M2T   = INDHER(T-2,0,0)
               DO 310 I = 1, NUABCD
                  CUR(I,TUV) = PQX(I)*OLD(I,M1T) + TMIN1*OLD(I,M2T)
  310          CONTINUE
  300       CONTINUE
         ELSE
            DO 400 T = 2, MIN(MAXT,JVAL), 2
               TMIN1 = T - 1.0D0
               TUV   = INDHER(T  ,0,0)
               M2T   = INDHER(T-2,0,0)
               DO 410 I = 1, NUABCD
                  CUR(I,TUV) = TMIN1*OLD(I,M2T)
  410          CONTINUE
  400       CONTINUE
         END IF
C
C        R(T,U,0)
C
         IF (PQYGT0) THEN
            DO 500 T = 0, MIN(MAXT,JVAL - 1), IPQX
               TUV = INDHER(T,1,0)
               M1U = INDHER(T,0,0)
               DO 510 I = 1, NUABCD
                  CUR(I,TUV) = PQY(I)*OLD(I,M1U)
  510          CONTINUE
  500       CONTINUE
            DO 600 U = 2, MIN(MAXU,JVAL)
               UMIN1  = U - 1.0D0
               DO 610 T = 0, MIN(MAXT,JVAL - U), IPQX
                  TUV = INDHER(T,U  ,0)
                  M1U = INDHER(T,U-1,0)
                  M2U = INDHER(T,U-2,0)
                  DO 620 I = 1, NUABCD
                     CUR(I,TUV) = PQY(I)*OLD(I,M1U) + UMIN1*OLD(I,M2U)
  620             CONTINUE
  610          CONTINUE
  600       CONTINUE
         ELSE
            DO 700 U = 2, MIN(MAXU,JVAL), 2
               UMIN1  = U - 1.0D0
               DO 710 T = 0, MIN(MAXT,JVAL - U), IPQX
                  TUV = INDHER(T,U  ,0)
                  M2U = INDHER(T,U-2,0)
                  DO 720 I = 1, NUABCD
                     CUR(I,TUV) = UMIN1*OLD(I,M2U)
  720             CONTINUE
  710          CONTINUE
  700       CONTINUE
         END IF
C
C        R(T,U,V)
C
         IF (PQZGT0) THEN
            IUMAX  = JVAL - 1
            DO 800 U = 0, MIN(MAXU,IUMAX), IPQY
               DO 810 T = 0, MIN(MAXT,IUMAX - U), IPQX
                  TUV = INDHER(T,U,1)
                  M1V = INDHER(T,U,0)
                  DO 820 I = 1, NUABCD
                     CUR(I,TUV) = PQZ(I)*OLD(I,M1V)
  820             CONTINUE
  810          CONTINUE
  800       CONTINUE
            DO 900 V = 2, MIN(MAXV,JVAL)
               VMIN1  =  V - 1.0D0
               IUMAX  = JVAL - V
               DO 910 U = 0, MIN(MAXU,IUMAX), IPQY
                  DO 920 T = 0, MIN(MAXT,IUMAX - U), IPQX
                     TUV = INDHER(T,U,V  )
                     M1V = INDHER(T,U,V-1)
                     M2V = INDHER(T,U,V-2)
                     DO 930 I = 1, NUABCD
                        CUR(I,TUV) = PQZ(I)*OLD(I,M1V)+VMIN1*OLD(I,M2V)
  930                CONTINUE
  920             CONTINUE
  910          CONTINUE
  900       CONTINUE
         ELSE
            DO 1000 V = 2, MIN(MAXV,JVAL), 2
               VMIN1  = V - 1.0D0
               IUMAX  = JVAL - V
               DO 1010 U = 0, MIN(MAXU,IUMAX), IPQY
                  DO 1020 T = 0, MIN(MAXT,IUMAX - U), IPQX
                     TUV = INDHER(T,U,V  )
                     M2V = INDHER(T,U,V-2)
                     DO 1030 I = 1, NUABCD
                        CUR(I,TUV) = VMIN1*OLD(I,M2V)
 1030                CONTINUE
 1020             CONTINUE
 1010          CONTINUE
 1000       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck herswp */
      SUBROUTINE HERSWP(HERINT,NTUV,NUABCD,WORK,LWORK,JMAX,NUCAB,
     &                  NUCCD,IPRINT)
C
C     TUH 87
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION HERINT(NUABCD,NTUV), WORK(LWORK)
      NHRINT = NTUV*NUABCD
      IF (NHRINT .GT. LWORK) CALL STOPIT('HERSWP',' ',LWORK,NHRINT)
      ISTR10 = 1
      ISTR2  = 1
      DO 100 I = 1, NUCAB
         ISTR1 = ISTR10
         DO 200 J = 1, NUCCD
            CALL DCOPY(NTUV,HERINT(ISTR2,1),NUABCD,WORK(ISTR1),NUABCD)
            ISTR1 = ISTR1 + NUCAB
            ISTR2 = ISTR2 + 1
  200    CONTINUE
         ISTR10 = ISTR10 + 1
  100 CONTINUE
      CALL DCOPY(NHRINT,WORK,1,HERINT,1)
      IF (IPRINT .GE. 25) THEN
         CALL HEADER('OUTPUT FROM HERSWP',-1)
         DO 500 ITUV = 1, NTUV
            WRITE (LUPRI,'(/,A,I5/)') ' NR ',ITUV
            WRITE (LUPRI,'(6F12.8)') (HERINT(I,ITUV),I=1,NUABCD)
  500    CONTINUE
      END IF
      RETURN
      END
C  /* Deck r000x */
#if defined (VAR_OLDGAM)
      SUBROUTINE R000X(RJ000,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,PQX,
     &                 PQY,PQZ,JMAX,NOINT,NUABCD,NUC1,NUC2,NUC12,NUC3,
     &                 NUC4,NUC34,THRESH,ONECEN,IPRINT,IPQ0X,IPQ0Y,
     &                 IPQ0Z)
C
C     TUH 84
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0, D2 = 2.D0)
      LOGICAL ONECEN, NOINT
      DIMENSION RJ000(NUABCD,0:JMAX),
     &          PQX(NUABCD), PQY(NUABCD), PQZ(NUABCD),
     &          COOR12(NUC1*NUC2,3), COOR34(NUC3*NUC4,3),
     &          EXP12 (NUC1*NUC2),   EXP34 (NUC3*NUC4),
     &          FAC12 (NUC1*NUC2),   FAC34 (NUC3*NUC4)
#include <gamcom.h>
#include <twosta.h>
#include <subdir.h>
      IF (.NOT.DPATH1) THEN
         SIGN   = - D1
      ELSE
         SIGN   = D1
      END IF
      JMAX0 = JMAX
      NOINT = .TRUE.
      IPQ0X = 1
      IPQ0Y = 1
      IPQ0Z = 1
C
C     *****************
C     ***** RJ000 *****
C     *****************
C
C     Special case: One-center Integrals
C     ==================================
C
      IF (ONECEN) THEN
         IODS = 0
         DO 200 IOD12 = 1, NUC12
            EXPP   = EXP12(IOD12)
            FAC12I = FAC12(IOD12)
            DO 210 IOD34 = 1, NUC34
               IODS   = IODS + 1
               EXPQ   = EXP34(IOD34)
               FAC34I = FAC34(IOD34)
               EXPPQI = EXPP + EXPQ
               FACTOR = FAC12I*FAC34I/SQRT(EXPPQI)
               IF (ABS(FACTOR) .GT. THRESH) THEN
                  NOINT = .FALSE.
                  ALPHA = EXPP*EXPQ/EXPPQI
                  PQX(IODS) = D0
                  PQY(IODS) = D0
                  PQZ(IODS) = D0
                  TALPHA = ALPHA + ALPHA
                  FJ0INV = D1
                  DO 220 I = 0, JMAX
                     RJ000(IODS,I) = FACTOR/FJ0INV
                     FACTOR = - TALPHA*FACTOR
                     FJ0INV = FJ0INV + D2
  220             CONTINUE
               ELSE
                  DO 225 I = 0 ,JMAX
                     RJ000(IODS,I) = D0
  225             CONTINUE
               END IF
  210       CONTINUE
  200    CONTINUE
C
C     General case: Multicenter Integrals
C     ===================================
C
      ELSE
         IODS = 0
         DO 300 IOD12 = 1, NUC12
            EXPP   = EXP12(IOD12)
            FAC12I = FAC12(IOD12)
            PX     = COOR12(IOD12,1)
            PY     = COOR12(IOD12,2)
            PZ     = COOR12(IOD12,3)
            DO 310 IOD34 = 1, NUC34
               IODS   = IODS + 1
               EXPQ   = EXP34(IOD34)
               FAC34I = FAC34(IOD34)
               EXPPQI = EXPP + EXPQ
               FACTOR = FAC12I*FAC34I/SQRT(EXPPQI)
               IF (ABS(FACTOR) .GT. THRESH) THEN
                  NOINT = .FALSE.
                  ALPHA = EXPP*EXPQ/EXPPQI
                  PQXI = PX - COOR34(IOD34,1)
                  PQYI = PY - COOR34(IOD34,2)
                  PQZI = PZ - COOR34(IOD34,3)
                  IF (ABS(PQXI) .GT. THRESH) IPQ0X = 0
                  IF (ABS(PQYI) .GT. THRESH) IPQ0Y = 0
                  IF (ABS(PQZI) .GT. THRESH) IPQ0Z = 0
                  PQX(IODS) = SIGN*PQXI
                  PQY(IODS) = SIGN*PQYI
                  PQZ(IODS) = SIGN*PQZI
                  WVAL = ALPHA*(PQXI*PQXI + PQYI*PQYI + PQZI*PQZI)
                  CALL GAMFUN
                  TALPHA = ALPHA + ALPHA
                  DO 320 I = 0, JMAX
                     RJ000(IODS,I) = FACTOR*FJW(I)
                     FACTOR = - TALPHA*FACTOR
  320             CONTINUE
               ELSE
                  DO 325 I = 0 ,JMAX
                     RJ000(IODS,I) = D0
  325             CONTINUE
               END IF
  310       CONTINUE
  300    CONTINUE
      END IF
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1010) JMAX
         WRITE (LUPRI, 1020) NUC12, NUC34
         WRITE (LUPRI, 1030) THRESH
         WRITE (LUPRI, 1040) ONECEN
         WRITE (LUPRI, 1045) NOINT
         WRITE (LUPRI, 1050) NUABCD
         IF (IPRINT .GT. 20) THEN
            WRITE (LUPRI, 1100)
            ISTART = 0
            DO 2000 J = 0, JMAX
               WRITE (LUPRI, 1110) J
               IADR = 0
               DO 2100 I = 1, NUC12
                  WRITE (LUPRI, 1120) I
                  WRITE (LUPRI, 1130) (RJ000(IADR + K,J), K = 1, NUC34)
                  IADR = IADR + NUC34
 2100          CONTINUE
 2000       CONTINUE
         END IF
      END IF
 1000 FORMAT (//,' ********** SUBROUTINE R000 **********')
 1010 FORMAT (//,'  JMAX    ',I7)
 1020 FORMAT (   '  NUC     ',2I7)
 1030 FORMAT (   '  THRESH  ',1P,D12.4)
 1040 FORMAT (   '  ONECEN  ',L7)
 1045 FORMAT (   '  NOINT   ',L7)
 1050 FORMAT (   '  NUABCD  ',I7)
 1100 FORMAT(//,' ***** RJ000 - INTEGRALS *****')
 1110 FORMAT( /,'  J     ',I7)
 1120 FORMAT( /,'  NUC12 ',I7,/)
 1130 FORMAT(1P,6D12.4)
      RETURN
      END
#endif
