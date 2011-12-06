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
C  /* Deck odcvec */
      SUBROUTINE ODCVEC(COOR12,EXP12,FAC12,CONT1,CONT2,JMAX1,JMAX2,
     &                  NSET1,NSET2,NUC1,NUC2,NUCT12,NORB1,NORB2,NPCO1,
     &                  NPCO2,NUCS12,JSTR1,JSTR2,TCON12,TPR12,GEN12,
     &                  ITYPE,THRESH,MAXDER,MUL1,MUL2,NODC12,NORT12,
     &                  NIND12,NPNT12,NRED12,KHKT1,KHKT2,EXPECT,DIRFCK,
     &                  WORK,LWORK,RPRI12,RCNT12,IXP12,IPRINT)
C
C     TUH Apr 11 1988
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <iratdef.h>
      LOGICAL TPR12, TCON12, GEN12, EXPECT, DIRFCK, RPRI12, RCNT12
      DIMENSION NPCO1(*), NPCO2(*), JSTR1(*), JSTR2(*), COOR12(*),
     &          EXP12(*), FAC12(*), CONT1(*), CONT2(*), NUCS12(*),
     &          NUCT12(*), NORT12(*), NIND12(*), NPNT12(*), IXP12(*),
     &          WORK(LWORK)
#include <twosta.h>
      IF (TKTIME) TIMSTR = SECOND()
      IF (GEN12) THEN
         KPOINT = 1
         KREDP1 = KPOINT + (2*NUC1*NUC2*NODC12 + 1)/IRAT
         KREDP2 = KREDP1 + (NUC1*NODC12 + 1)/IRAT
         KREDC1 = KREDP2 + (NUC2*NODC12 + 1)/IRAT
         KREDC2 = KREDC1 + (NORB1 + 1)/IRAT
         KREDCC = KREDC2 + (NORB2 + 1)/IRAT
         KCNT1  = KREDCC + (NORB1*NORB2 + 1)/IRAT
         KCNT2  = KCNT1  + NUC1*NORB1
         KCMX1  = KCNT2  + NUC2*NORB2
         KCMX2  = KCMX1  + NUC1
         KFCPRM = KCMX2  + NUC2
         KFACCP = KFCPRM + NUC1*NUC2*NODC12
         KFACNT = KFACCP + NORB1*NUC2
         KFTOR1 = KFACNT + NORB1*NORB2
         KFTOR2 = KFTOR1 + (NUC1 + 1)/IRAT
         KLAST  = KFTOR2 + (NUC2 + 1)/IRAT
         IF (KLAST .GT. LWORK) CALL STOPIT('ODCVEC',' ',KLAST,LWORK)
         CALL ODCGEN(COOR12,EXP12,FAC12,CONT1,CONT2,NUC1,NUC2,NORB1,
     &               NORB2,NSET1,NSET2,NPCO1,NPCO2,JSTR1,JSTR2,NUCT12,
     &               TCON12,TPR12,THRESH,ITYPE,MUL1,MUL2,NODC12,NORT12,
     &               NIND12,NPNT12,NRED12,KHKT1,KHKT2,EXPECT,DIRFCK,
     &               WORK(KPOINT),WORK(KREDP1),WORK(KREDP2),
     &               WORK(KREDC1),WORK(KREDC2),WORK(KREDCC),
     &               WORK(KCNT1),WORK(KCNT2),WORK(KCMX1),WORK(KCMX2),
     &               WORK(KFCPRM),WORK(KFACCP),WORK(KFACNT),
     &               WORK(KFTOR1),WORK(KFTOR2),RPRI12,RCNT12,IPRINT)
      ELSE
         CALL ODCSEG(COOR12,EXP12,FAC12,NUC1,NUC2,NORB1,NORB2,NSET1,
     &               NSET2,NPCO1,NPCO2,NUCS12,JSTR1,JSTR2,NUCT12,TCON12,
     &               TPR12,THRESH,ITYPE,MUL1,MUL2,NODC12,NORT12,
     &               NIND12,KHKT1,KHKT2,EXPECT,DIRFCK,IXP12,IPRINT)
      END IF
      IF (TKTIME) THEN
         TIMEND = SECOND()
         TODCVE = TODCVE + TIMEND - TIMSTR
         TIMSTR = TIMEND
      END IF
      RETURN
      END
C  /* Deck odcgen */
      SUBROUTINE ODCGEN(COOR12,EXP12,FAC12,CONT1,CONT2,NUC1,NUC2,
     &                  NORB1,NORB2,NSET1,NSET2,NPCO1,NPCO2,
     &                  JSTR1,JSTR2,NUCT12,TCON12,TPR12,THRESH,ITYPE,
     &                  MUL1,MUL2,NODC12,NORT12,NIND12,NPNT12,
     &                  NRED12,KHKT1,KHKT2,EXPECT,DIRFCK,NPOINT,NREDP1,
     &                  NREDP2,NREDC1,NREDC2,NREDCC,CNT1,CNT2,CFMAX1,
     &                  CFMAX2,FACPRM,FACCP,FACCNT,NFTOR1,NFTOR2,RPRI12,
     &                  RCNT12,IPRINT)
C
C     TUH 84
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <pi.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
      LOGICAL TCON12, TPR12, ALLSYM, EXPECT, DIRFCK, RPRI12, RCNT12,
     &        DOTHIS
#include <primit.h>
#include <symmet.h>
      COMMON /MAXOLD/ SRFMAX,OLDMAX
      DIMENSION NPCO1(NSET1,2,0:NODC12), NPCO2(NSET2,2,0:NODC12),
     &          JSTR1(*), JSTR2(*),
     &          COOR12(NUC1*NUC2,3,3,NODC12),
     &          EXP12(NUC1*NUC2,3,NODC12),
     &          FAC12(NUC1*NUC2,NODC12), NUCT12(NODC12),
     &          NPRI12(8), NIND12(NORB1*NORB2,2), NORT12(8),
     &          CONT1(NORB1*NUC1,2,NODC12), CONT2(NORB2*NUC2,2,NODC12),
     &          NPOINT(NUC1*NUC2,2,NODC12), NPNT12(NUC1*NUC2,2,NODC12),
     &          NREDP1(NUC1,NODC12), NREDP2(NUC2,NODC12),
     &          NREDC1(NORB1), NREDC2(NORB2),
     &          CNT1(NORB1,NUC1), CNT2(NUC2,NORB2),
     &          CFMAX1(NUC1), CFMAX2(NUC2),
     &          FACPRM(NUC1,NUC2,NODC12), FACCP(NORB1,NUC2),
     &          FACCNT(NORB1,NORB2),
     &          XFAC(8), YFAC(8), ZFAC(8),
     &          NREDCC(NORB1,NORB2), NRED12(NORB1*NORB2),
     &          NFTOR1(NUC1), NFTOR2(NUC2)
C

#include <ibtfun.h>
      XAND(I) = PT(IBTAND(ISYMAX(1,1),I))
      YAND(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IBTAND(ISYMAX(3,1),I))
C
      CALL DZERO(FACPRM,NUC1*NUC2*NODC12)
C
C     ***********************************************
C     ***** Set up overlap distribution vectors *****
C     ***********************************************
C
      IF (NSET1.NE.1 .OR. NSET2.NE.1) THEN
         WRITE (LUPRI,'(1X,A)')     ' NSET must be 1 in this version.'
         WRITE (LUPRI,'(1X,A,2I5)') ' NSET1, NSET2: ',NSET1,NSET2
         CALL QUIT('Illegal value of NSET in ODCGEN.')
      END IF
C
      IF (ITYPE .EQ. 12) THEN
         IELCTR = 1
         SRFMAX = D0
         THRSH  = THRESH
      ELSE
         IELCTR = 2
C
C     This adjustment of screening thresholds seems to lead to stalling and
C     odd convergence behaviour. Disabled. kr-230106
C     hjaaj-071017: reenabled; I think the problem was SRFMAX was not saved
C
         THRSH  = THRESH/OLDMAX
Chj      THRSH  = THRESH
      END IF
C
C     XFAC, YFAC and ZFAC
C     ===================
C
      ISYM = 0
      DO 100 ISYMOP = 0, MAXOPR
      IF (IBTAND(ISYMOP,IBTOR(MUL1,MUL2)) .EQ. 0) THEN
         ISYM = ISYM + 1
         XFAC(ISYM) = XAND(ISYMOP)
         YFAC(ISYM) = YAND(ISYMOP)
         ZFAC(ISYM) = ZAND(ISYMOP)
      END IF
  100 CONTINUE
C
C     Contraction matrices (no reductions)
C     ====================================
C
C hj aug99: I note that screening is disabled (DOTHIS always true)
C  and a lot of superfluous work is therefore done. My guess is that
C  this has been done because of symmetry breaking when a primitive
C  has been included sometimes and sometimes not. Therefore I change
C  the test to an absolute CFMAX test instead of for each primitive,
C  this must be safe, and reactivate screening.
C  - inserted CFMX1 and CFMX2 and CFMAX = CFMX1*CFMX2
C
      CFMX1 = D0
      DO 200 I = 1, NUC1
         II = JSTR1(1) + I
         CFMAX = D0
         DO 210 J = 1, NORB1
            CCF = PRICCF(II,J)
            CNT1(J,I) = CCF
            CFMAX = MAX(CFMAX,ABS(CCF))
  210    CONTINUE
Chj      CFMAX1(I) = CFMAX
         CFMX1 = MAX(CFMX1,CFMAX)
  200 CONTINUE
      CFMX2 = D0
      DO 205 I = 1, NUC2
         II = JSTR2(1) + I
         CFMAX = D0
         DO 215 J = 1, NORB2
            CCF = PRICCF(II,J)
            CNT2(I,J) = CCF
            CFMAX = MAX(CFMAX,ABS(CCF))
  215    CONTINUE
Chj      CFMAX2(I) = CFMAX
         CFMX2 = MAX(CFMX2,CFMAX)
  205 CONTINUE
      CFMAX = CFMX1*CFMX2
C
C     Primitive vectors
C     =================
C
      CALL IZERO(NUCT12,NODC12)
      CALL IZERO(NPRI12,NODC12)
      CALL IZERO(NREDP1,NUC1*NODC12)
      CALL IZERO(NREDP2,NUC2*NODC12)
      CALL IZERO(NPOINT,NUC1*NUC2*NODC12)
      DO 300 IPRM2 = 1, NUC2
         IPRIM2 = JSTR2(1) + IPRM2
         CRX20  = PRICRX(IPRIM2)
         CRY20  = PRICRY(IPRIM2)
         CRZ20  = PRICRZ(IPRIM2)
         EXP2   = PRIEXP(IPRIM2)
Chj      CFMX2  = CFMAX2(IPRM2)
         DO 310 IPRM1 = 1, NUC1
            IPRIM1 = JSTR1(1) + IPRM1
            CRX1   = PRICRX(IPRIM1)
            CRY1   = PRICRY(IPRIM1)
            CRZ1   = PRICRZ(IPRIM1)
            EXP1   = PRIEXP(IPRIM1)
Chj         CFMAX  = CFMAX1(IPRM1)*CFMX2
            EXPP   = EXP1 + EXP2
            EXPPI  = D1/EXPP
            EXPQ   = EXP1*EXP2*EXPPI
            PREXP  = R2PI52*EXPPI
C
C           screening test
C
C           DOTHIS = .FALSE.
C dothis set to .true. 
C Must not be futher changed without extensive
C testing on large molecules! tuh 08.12.05
C
            DOTHIS = .TRUE.
            DO 315 ISYM = 1, NODC12
Chj: next line moved from 320 loop to here
               NPRI12(ISYM) = NPRI12(ISYM) + 1
               CRX2   = XFAC(ISYM)*CRX20
               CRY2   = YFAC(ISYM)*CRY20
               CRZ2   = ZFAC(ISYM)*CRZ20
               DIFX   = CRX1 - CRX2
               DIFY   = CRY1 - CRY2
               DIFZ   = CRZ1 - CRZ2
               DIST   = DIFX*DIFX + DIFY*DIFY + DIFZ*DIFZ
               FAC12I = PREXP*EXP(-EXPQ*DIST)
               ABSFAC = ABS(CFMAX*FAC12I)
               IF (ABSFAC .GT. THRSH) THEN
                  SRFMAX = MAX(SRFMAX,ABSFAC)
                  DOTHIS = .TRUE.
               END IF
  315       CONTINUE
Chj aug99: inserted next line
            IF (.NOT. DOTHIS) GO TO 310
C
            DO 320 ISYM = 1, NODC12
               CRX2   = XFAC(ISYM)*CRX20
               CRY2   = YFAC(ISYM)*CRY20
               CRZ2   = ZFAC(ISYM)*CRZ20
               DIFX   = CRX1 - CRX2
               DIFY   = CRY1 - CRY2
               DIFZ   = CRZ1 - CRZ2
               DIST   = DIFX*DIFX + DIFY*DIFY + DIFZ*DIFZ
Chj            NPRI12(ISYM) = NPRI12(ISYM) + 1
               FAC12I = PREXP*EXP(-EXPQ*DIST)
               IPRI12 = NPRI12(ISYM)
Chj            ABSFAC = ABS(CFMAX*FAC12I)
Chj            IF (DOTHIS) THEN
                  NUCT12(ISYM) = NUCT12(ISYM) + 1
                  NREDP1(IPRM1,ISYM) = 1
                  NREDP2(IPRM2,ISYM) = 1
Chj               SRFMAX = MAX(SRFMAX,ABSFAC)
                  EXP1PI = EXP1*EXPPI
                  EXP2PI = EXP2*EXPPI
                  CORPX  = EXP1PI*CRX1 + EXP2PI*CRX2
                  CORPY  = EXP1PI*CRY1 + EXP2PI*CRY2
                  CORPZ  = EXP1PI*CRZ1 + EXP2PI*CRZ2
                  IPRIM  = NUCT12(ISYM)
                  NPOINT(IPRIM,1,ISYM)   = IPRM1
                  NPOINT(IPRIM,2,ISYM)   = IPRM2
                  EXP12 (IPRIM,1,  ISYM) = EXPP
                  EXP12 (IPRIM,2,  ISYM) = D2*EXP1
                  EXP12 (IPRIM,3,  ISYM) = D2*EXP2
                  FAC12 (IPRIM,    ISYM) = FAC12I
                  COOR12(IPRIM,1,1,ISYM) = CORPX
                  COOR12(IPRIM,2,1,ISYM) = CORPY
                  COOR12(IPRIM,3,1,ISYM) = CORPZ
                  COOR12(IPRIM,1,2,ISYM) = CORPX - CRX1
                  COOR12(IPRIM,2,2,ISYM) = CORPY - CRY1
                  COOR12(IPRIM,3,2,ISYM) = CORPZ - CRZ1
                  COOR12(IPRIM,1,3,ISYM) = CORPX - CRX2
                  COOR12(IPRIM,2,3,ISYM) = CORPY - CRY2
                  COOR12(IPRIM,3,3,ISYM) = CORPZ - CRZ2
Chj            END IF
               FACPRM(IPRM1,IPRM2,ISYM) = FAC12I
  320       CONTINUE
  310    CONTINUE
  300 CONTINUE
      IF (ISUM(NODC12,NUCT12,1) .EQ. 0) RETURN
      IF (ITYPE .EQ. 12) OLDMAX = SRFMAX
C
C     NPNT12
C     ======
C
      RPRI12 = .FALSE.
      DO 400 ISYM = 1, NODC12
         IF (NUCT12(ISYM) .LT. NPRI12(ISYM)) RPRI12 = .TRUE.
         NUCR1 = ISUM(NUC1,NREDP1(1,ISYM),1)
         NUCR2 = ISUM(NUC2,NREDP2(1,ISYM),1)
C
         CALL IZERO(NFTOR1,NUC1)
         CALL IZERO(NFTOR2,NUC2)
C
         IRED = 0
         DO 410 I = 1, NUC1
            IF (NREDP1(I,ISYM) .EQ. 1) THEN
               IRED = IRED + 1
               NFTOR1(I) = IRED
            END IF
  410    CONTINUE
C
         IRED = 0
         DO 420 I = 1, NUC2
            IF (NREDP2(I,ISYM) .EQ. 1) THEN
               IRED = IRED + 1
               NFTOR2(I) = IRED
            END IF
  420    CONTINUE
C
         DO 430 IJ = 1, NUCT12(ISYM)
            I = NFTOR1(NPOINT(IJ,1,ISYM))
            J = NFTOR2(NPOINT(IJ,2,ISYM))
            NPNT12(IJ,1,ISYM) = (J - 1)*NUCR1 + I
            NPNT12(IJ,2,ISYM) = (I - 1)*NUCR2 + J
  430    CONTINUE
         NPCO1(1,1,ISYM) = NUCR1
         NPCO2(1,1,ISYM) = NUCR2
  400 CONTINUE
C
C     Screening of contracted orbitals
C     ================================
C
      CALL IZERO(NREDC1,NORB1)
      CALL IZERO(NREDC2,NORB2)
      CALL IZERO(NREDCC,NORB1*NORB2)
      DO 500 ISYM = 1, NODC12
         CALL MXMA(CNT1,1,NORB1,FACPRM(1,1,ISYM),1,NUC1,FACCP,1,NORB1,
     &             NORB1,NUC1,NUC2)
         CALL MXMA(FACCP,1,NORB1,CNT2,1,NUC2,FACCNT,1,NORB1,
     &             NORB1,NUC2,NORB2)
ckr         CALL DGEMM('N','N',NORB1,NUC2,NUC1,D1,CNT1,NORB1,
ckr     &        FACPRM(1,1,ISYM),NUC1,D0,FACCP,NORB1)
ckr         CALL DGEMM('N','N',NORB1,NORB2,NUC2,D1,FACCP,NORB1,CNT2,NUC2,
ckr     &        D0,FACCNT,NORB1)
         DO 510 J = 1, NORB2
         DO 510 I = 1, NORB1
c            IF (ABS(FACCNT(I,J)) .GT. THRSH) THEN
               NREDC1(I)   = 1
               NREDC2(J)   = 1
               NREDCC(I,J) = 1
c            END IF
  510    CONTINUE
  500 CONTINUE
      NORR1 = ISUM(NORB1,NREDC1,1)
      NORR2 = ISUM(NORB2,NREDC2,1)
      IF ((NORR1.EQ.0) .OR. (NORR2.EQ.0)) THEN
         CALL IZERO(NUCT12,NODC12)
         RETURN
      END IF
C
C     NRED12
C     ======
C
      IJRED = 0
      IRED  = 0
      DO 600 I = 1, NORB1
         IF (NREDC1(I) .EQ. 1) IRED = IRED + 1
         JRED = 0
         JEND = NORB2
         IF (TCON12) JEND = I
         DO 610 J = 1, JEND
            IF (NREDC2(J) .EQ. 1) JRED = JRED + 1
            IF (NREDCC(I,J) .EQ. 1) THEN
               IJRED = IJRED + 1
               NRED12(IJRED) = (IRED - 1)*NORR2 + JRED
            END IF
  610    CONTINUE
  600 CONTINUE
C
C     NIND12
C     ======
C
      IJRED = 0
      DO 700 I = 1, NORB1
         JEND = NORB2
         IF (TCON12) JEND = I
         DO 710 J = 1, JEND
            IF (NREDCC(I,J) .EQ. 1) THEN
               IJRED = IJRED + 1
               NIND12(IJRED,1) = I
               NIND12(IJRED,2) = J
            END IF
  710    CONTINUE
  700 CONTINUE
C
      IJALL = NORB1*NORB2
      IF (TCON12) IJALL = NORB1*(NORB1 + 1)/2
      CALL IZERO(NORT12,NODC12)
      NORT12(1) = IJRED
      RCNT12    = IJRED .LT. IJALL
C
C     Reduced contraction matrices
C     ============================
C
      DO 800 ISYM = 1, NODC12
         NPCO1(1,2,ISYM) = NORR1
         NPRIU  = NPCO1(1,1,ISYM)
         IADR1  = 1
         IADR21 = 1
         DO 820 I = 1, NUC1
         IF (NREDP1(I,ISYM).EQ.1) THEN
            IADR2 = IADR21
            DO 830 J = 1, NORB1
            IF (NREDC1(J) .EQ. 1) THEN
               CCF = CNT1(J,I)
               CONT1(IADR1,1,ISYM) = CCF
               CONT1(IADR2,2,ISYM) = CCF
               IADR1 = IADR1 + 1
               IADR2 = IADR2 + NPRIU
            END IF
  830       CONTINUE
            IADR21 = IADR21 + 1
         END IF
  820    CONTINUE
  800 CONTINUE
C
      DO 805 ISYM = 1, NODC12
         NPCO2(1,2,ISYM) = NORR2
         NPRIU  = NPCO2(1,1,ISYM)
         IADR1  = 1
         IADR21 = 1
         DO 825 I = 1, NUC2
         IF (NREDP2(I,ISYM).EQ.1) THEN
            IADR2 = IADR21
            DO 835 J = 1, NORB2
            IF (NREDC2(J) .EQ. 1) THEN
               CCF = CNT2(I,J)
               CONT2(IADR1,1,ISYM) = CCF
               CONT2(IADR2,2,ISYM) = CCF
               IADR1 = IADR1 + 1
               IADR2 = IADR2 + NPRIU
            END IF
  835       CONTINUE
            IADR21 = IADR21 + 1
         END IF
  825    CONTINUE
  805 CONTINUE
C
      DO 900 ISYM = 1, NODC12
         CALL REDSTA(0,NUCT12(ISYM),NPRI12(ISYM),NORT12(1),IJALL)
  900 CONTINUE
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         WRITE(LUPRI, 2000) IELCTR
         WRITE(LUPRI,'(A,2I5)')' NUC   ', NUC1, NUC2
         WRITE(LUPRI,'(A,3I5)')' NCONT ', NORB1, NORB2, IJRED
         WRITE(LUPRI,'(A,2I5)')' NSET  ', NSET1, NSET2
         WRITE(LUPRI,'(A,2L5)')' RPRI12,RCNT12', RPRI12, RCNT12
         WRITE(LUPRI,'(A,8I5)')' NUC12 ', (NUCT12(I),I=1,NODC12)
         WRITE(LUPRI,'(A,8I5)')' NPCO1(*,1,0)',(NPCO1(I,1,0),I=1,NODC12)
         WRITE(LUPRI,'(A,8I5)')' NPCO1(*,2,0)',(NPCO1(I,1,0),I=1,NODC12)
         WRITE(LUPRI,'(A,8I5)')' NPCO2(*,1,0)',(NPCO2(I,1,0),I=1,NODC12)
         WRITE(LUPRI,'(A,8I5)')' NPCO2(*,2,0)',(NPCO2(I,1,0),I=1,NODC12)
         WRITE(LUPRI,'(A,12I5/,(7X,12I5))')
     &      ' NRED12',(NRED12(I),I=1,IJRED)
         DO 1000 ISYM = 1, NODC12
            WRITE(LUPRI,'(/,A,I5,/)') ' Symmetry (compressed):',ISYM
            WRITE(LUPRI,'(A,12I5/,(8X,12I5))') ' NPNT12 - 1',
     &           (NPNT12(I,1,ISYM),I=1,NUCT12(ISYM))
            WRITE(LUPRI,'(A,12I5/,(8X,12I5))') ' NPNT12 - 2',
     &           (NPNT12(I,2,ISYM),I=1,NUCT12(ISYM))
            IOFF = 1
            DO 1010 I = 1, NSET1
               WRITE (LUPRI, '(/A,I2,A)')
     &            ' Contraction matrices for set ',I,'.'
               NPRM1 = NPCO1(I,1,ISYM)
               NRC1  = NPCO1(I,2,ISYM)
               WRITE (LUPRI, '(A,2I5)')
     &          ' Number of prim. and contr. functions:',NPRM1,NRC1
               CALL HEADER('CONT1',-1)
               CALL OUTPUT(CONT1(IOFF,1,ISYM),1,NRC1,1,NPRM1,
     &                     NRC1,NPRM1,1,LUPRI)
               CALL HEADER('CONT1 transposed',-1)
               CALL OUTPUT(CONT1(IOFF,2,ISYM),1,NPRM1,1,NRC1,NPRM1,
     &                     NRC1,1,LUPRI)
               IOFF = IOFF + NPRM1*NRC1
 1010       CONTINUE
            IOFF = 1
            DO 1020 I = 1, NSET2
               WRITE (LUPRI, '(/A,I2,A)')
     &                   ' Contraction matrices for set ',I,'.'
               NPRM2 = NPCO2(I,1,ISYM)
               NRC2  = NPCO2(I,2,ISYM)
               WRITE (LUPRI, '(A,2I5)')
     &          ' Number of prim. and contr. functions:',NPRM2,NRC2
               CALL HEADER('CONT2',-1)
               CALL OUTPUT(CONT2(IOFF,1,ISYM),1,NRC2,1,NPRM2,
     &                     NRC2,NPRM2,1,LUPRI)
               CALL HEADER('CONT2 transposed',-1)
               CALL OUTPUT(CONT2(IOFF,2,ISYM),1,NPRM2,1,NRC2,NPRM2,
     &                     NRC2,1,LUPRI)
               IOFF = IOFF + NPRM2*NRC2
 1020       CONTINUE
 1000    CONTINUE
      END IF
 2000 FORMAT(//,' <<<<<<<<<< SUBROUTINE ODCGEN >>>>>>>>>>',
     *       //,' Overlap distribution arrays for electron ',I1,/)
      RETURN
      END
C  /* Deck odcseg */
      SUBROUTINE ODCSEG(COOR12,EXP12,FAC12,NUC1,NUC2,NORB1,NORB2,
     &                  NSET1,NSET2,NPCO1,NPCO2,NUCS12,JSTR1,JSTR2,
     &                  NUCT12,TCON12,TPR12,THRESH,ITYPE,MUL1,MUL2,
     &                  NODC12,NORT12,NIND12,KHKT1,KHKT2,EXPECT,DIRFCK,
     &                  IXP12,IPRINT)
C
C     TUH 84
C
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <pi.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
      LOGICAL TCON12, TPR12, DIAG12, ALLSYM, EXPECT, DIRFCK, DOTHIS
#include <primit.h>
#include <symmet.h>
#include <expopt.h>
      COMMON /MAXOLD/ SRFMAX,OLDMAX
      DIMENSION NPCO1(NSET1,2), NPCO2(NSET2,2), JSTR1(*), JSTR2(*),
     &          COOR12(NUC1*NUC2,3,3,NODC12),
     &          EXP12(NUC1*NUC2,3,NODC12),
     &          FAC12(NUC1*NUC2,NODC12), NUCT12(NODC12),
     &          NUCS12(NORB1*NORB2,NODC12), NPRI12(8),
     &          NIND12(NORB1*NORB2,2), NORT12(8),
     &          IXP12(NUC1*NUC2,2,NODC12)
C
#include <ibtfun.h>
      XAND(I) = PT(IBTAND(ISYMAX(1,1),I))
      YAND(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IBTAND(ISYMAX(3,1),I))
C
C     ***********************************************
C     ***** Set up overlap distribution vectors *****
C     ***********************************************
C
      IF (ITYPE .EQ. 12) THEN
         IELCTR = 1
         SRFMAX = D0
         THRSH  = THRESH
      ELSE
         IELCTR = 2
C
C     This adjustment of screening thresholds seems to lead to stalling and
C     odd convergence behaviour. Disabled. kr-230106
C     hjaaj-071017: reenabled; I think the problem was SRFMAX was not saved
C
         THRSH  = THRESH/OLDMAX
Chj      THRSH  = THRESH
      END IF
Chj aug99: find max relevant isymop
C          test is trivially true for ISYMOP .eq. 1
      MAXOP12 = 0
      DO ISYMOP = 1,MAXOPR
         IF (IBTAND(ISYMOP,IBTOR(MUL1,MUL2)) .EQ. 0) MAXOP12 = ISYMOP
      END DO
C
C     Run over contracted orbitals
C     ============================
C
      ISET12 = 0
      CALL IZERO(NUCT12,NODC12)
      CALL IZERO(NORT12,NODC12)
      CALL IZERO(NUCS12,NORB1*NORB2*NODC12)
      DO 100 ISET1 = 1, NSET1
         ISTR1 = JSTR1(ISET1) + 1
         IEND1 = JSTR1(ISET1) + NPCO1(ISET1,1)
         IORB1 = - NPCO1(ISET1,2)
         IF (TCON12) THEN
            JEND2 = ISET1
         ELSE
            JEND2 = NSET2
         END IF
         DO 200 ISET2 = 1, JEND2
            ISTR2 = JSTR2(ISET2) + 1
            IEND2 = JSTR2(ISET2) + NPCO2(ISET2,1)
            IORB2 = - NPCO2(ISET2,2)
            CALL IZERO(NPRI12,NODC12)
C
C           Run over primitives
C           ===================
C
C NECgh980908 We do integral screening for MAXOPR=0 in a bit more efficent
C way. Save 75% of odcseg time on NEC, but probably also a bit on other
C systems, but for this we have to duplicate the DO 300 loop.
C hj-aug99: bugfix (ISYMOP not defined) + moved ISYMOP test outside loops
C           (ISYMOP test is maybe superfluous for MAXOPR .eq. 0 ???)
C           + optimized further by using MAXOP12 instead of MAXOPR
C             and removing redundant code for MAXOPR.eq.0
C
          IF (MAXOP12 .EQ. 0) THEN
            NSYMOP = 1
            DO 301 IPRIM1 = ISTR1, IEND1
               CCF1 = PRICCF(IPRIM1,IORB1)
               EXP1 = PRIEXP(IPRIM1)
               CRX1 = PRICRX(IPRIM1)
               CRY1 = PRICRY(IPRIM1)
               CRZ1 = PRICRZ(IPRIM1)
               DO 401 IPRIM2 = ISTR2, IEND2
                  CCF2 = PRICCF(IPRIM2,IORB2)
                  EXP2 = PRIEXP(IPRIM2)
                  CRX2 = PRICRX(IPRIM2)
                  CRY2 = PRICRY(IPRIM2)
                  CRZ2 = PRICRZ(IPRIM2)
C
                  EXPP  = EXP1 + EXP2
                  EXPPI = D1/EXPP
                  CCF12 = CCF1*CCF2*R2PI52*EXPPI
                  EXPQ  = EXP1*EXP2*EXPPI
C
                  DIFX = CRX1 - CRX2
                  DIFY = CRY1 - CRY2
                  DIFZ = CRZ1 - CRZ2
                  DIST12 = DIFX*DIFX + DIFY*DIFY + DIFZ*DIFZ
                  FAC12I = CCF12*EXP(-EXPQ*DIST12)
                  ABSFAC = ABS(FAC12I)
C                 screening test:
                  IF (ABSFAC .GT. THRSH) THEN
                        NPRI12(NSYMOP) = NPRI12(NSYMOP) + 1
                        NUCT12(NSYMOP) = NUCT12(NSYMOP) + 1
                        SRFMAX = MAX(SRFMAX,ABSFAC)
                        EXP1PI = EXP1*EXPPI
                        EXP2PI = EXP2*EXPPI
                        CORPX  = EXP1PI*CRX1 + EXP2PI*CRX2
                        CORPY  = EXP1PI*CRY1 + EXP2PI*CRY2
                        CORPZ  = EXP1PI*CRZ1 + EXP2PI*CRZ2
                        IPRIM  = NUCT12(NSYMOP)
                        EXP12 (IPRIM,1,NSYMOP)   = EXPP
                        EXP12 (IPRIM,2,NSYMOP)   = D2*EXP1
                        EXP12 (IPRIM,3,NSYMOP)   = D2*EXP2
                        IF (EXPGRA) THEN
                           IXP12(ISET12+1,1,NSYMOP) = IPRIM1 
                           IXP12(ISET12+1,2,NSYMOP) = IPRIM2 
                        END IF
                        FAC12 (IPRIM,NSYMOP)     = FAC12I
                        COOR12(IPRIM,1,1,NSYMOP) = CORPX
                        COOR12(IPRIM,2,1,NSYMOP) = CORPY
                        COOR12(IPRIM,3,1,NSYMOP) = CORPZ
                        COOR12(IPRIM,1,2,NSYMOP) = CORPX - CRX1
                        COOR12(IPRIM,2,2,NSYMOP) = CORPY - CRY1
                        COOR12(IPRIM,3,2,NSYMOP) = CORPZ - CRZ1
                        COOR12(IPRIM,1,3,NSYMOP) = CORPX - CRX2
                        COOR12(IPRIM,2,3,NSYMOP) = CORPY - CRY2
                        COOR12(IPRIM,3,3,NSYMOP) = CORPZ - CRZ2
                  END IF
  401          CONTINUE
  301       CONTINUE
C NECgh980908 Here comes the original code for MAXOPR .GT. 0
          ELSE
            DO 300 IPRIM1 = ISTR1, IEND1
               CCF1 = PRICCF(IPRIM1,IORB1)
               EXP1 = PRIEXP(IPRIM1)
               CRX1 = PRICRX(IPRIM1)
               CRY1 = PRICRY(IPRIM1)
               CRZ1 = PRICRZ(IPRIM1)
               DO 400 IPRIM2 = ISTR2, IEND2
                  CCF2  = PRICCF(IPRIM2,IORB2)
                  EXP2  = PRIEXP(IPRIM2)
                  CRX20 = PRICRX(IPRIM2)
                  CRY20 = PRICRY(IPRIM2)
                  CRZ20 = PRICRZ(IPRIM2)
C
                  EXPP  = EXP1 + EXP2
                  EXPPI = D1/EXPP
                  CCF12 = CCF1*CCF2*R2PI52*EXPPI
                  EXPQ  = EXP1*EXP2*EXPPI
C
C                 screening test:
C
Chj               DOTHIS = .FALSE.
                  DO 510 ISYMOP = 0,MAXOP12
                  IF (IBTAND(ISYMOP,IBTOR(MUL1,MUL2)) .EQ. 0) THEN
                     CRX2 = XAND(ISYMOP)*CRX20
                     CRY2 = YAND(ISYMOP)*CRY20
                     CRZ2 = ZAND(ISYMOP)*CRZ20
                     DIFX = CRX1 - CRX2
                     DIFY = CRY1 - CRY2
                     DIFZ = CRZ1 - CRZ2
                     DIST12 = DIFX*DIFX + DIFY*DIFY + DIFZ*DIFZ
                     FAC12I = CCF12*EXP(-EXPQ*DIST12)
                     ABSFAC = ABS(FAC12I)
Chj                  IF (ABSFAC .GT. THRSH) DOTHIS = .TRUE.
                     IF (ABSFAC .GT. THRSH) GO TO 511
                  END IF
  510             CONTINUE
Chj aug99: do not calculate integral, next IPRIM2
C (DOTHIS test was inside 500 loop)
                  GO TO 400
C
  511             NSYMOP = 0
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
C                 Cray fails in some cases without this directive
CDIR$ NOVECTOR
#endif
                  DO 500 ISYMOP = 0,MAXOP12
                  IF (IBTAND(ISYMOP,IBTOR(MUL1,MUL2)) .EQ. 0) THEN
                     NSYMOP = NSYMOP + 1
                     CRX2 = XAND(ISYMOP)*CRX20
                     CRY2 = YAND(ISYMOP)*CRY20
                     CRZ2 = ZAND(ISYMOP)*CRZ20
                     DIFX = CRX1 - CRX2
                     DIFY = CRY1 - CRY2
                     DIFZ = CRZ1 - CRZ2
                     DIST12 = DIFX*DIFX + DIFY*DIFY + DIFZ*DIFZ
                     FAC12I = CCF12*EXP(-EXPQ*DIST12)
                     ABSFAC = ABS(FAC12I)
Chj                  IF (DOTHIS) THEN ! we now only enter here if needed
                        NPRI12(NSYMOP) = NPRI12(NSYMOP) + 1
                        NUCT12(NSYMOP) = NUCT12(NSYMOP) + 1
                        SRFMAX = MAX(SRFMAX,ABSFAC)
                        EXP1PI = EXP1*EXPPI
                        EXP2PI = EXP2*EXPPI
                        CORPX  = EXP1PI*CRX1 + EXP2PI*CRX2
                        CORPY  = EXP1PI*CRY1 + EXP2PI*CRY2
                        CORPZ  = EXP1PI*CRZ1 + EXP2PI*CRZ2
                        IPRIM  = NUCT12(NSYMOP)
                        EXP12 (IPRIM,1,NSYMOP)   = EXPP
                        EXP12 (IPRIM,2,NSYMOP)   = D2*EXP1
                        EXP12 (IPRIM,3,NSYMOP)   = D2*EXP2
                        IF (EXPGRA) THEN
                           IXP12(ISET12+1,1,NSYMOP) = IPRIM1 
                           IXP12(ISET12+1,2,NSYMOP) = IPRIM2 
                        END IF
                        FAC12 (IPRIM,NSYMOP)     = FAC12I
                        COOR12(IPRIM,1,1,NSYMOP) = CORPX
                        COOR12(IPRIM,2,1,NSYMOP) = CORPY
                        COOR12(IPRIM,3,1,NSYMOP) = CORPZ
                        COOR12(IPRIM,1,2,NSYMOP) = CORPX - CRX1
                        COOR12(IPRIM,2,2,NSYMOP) = CORPY - CRY1
                        COOR12(IPRIM,3,2,NSYMOP) = CORPZ - CRZ1
                        COOR12(IPRIM,1,3,NSYMOP) = CORPX - CRX2
                        COOR12(IPRIM,2,3,NSYMOP) = CORPY - CRY2
                        COOR12(IPRIM,3,3,NSYMOP) = CORPZ - CRZ2
Chj                  END IF
                  END IF
  500          CONTINUE
  400          CONTINUE
  300       CONTINUE
          END IF
C NECgh980908 End of duplicate code for MAXOPR .GT. 0
            IF (ISUM(NODC12,NPRI12,1) .GT. 0) THEN
               ISET12 = ISET12 + 1
               CALL ICOPY(NODC12,NPRI12,1,NUCS12(ISET12,1),NORB1*NORB2)
               NIND12(ISET12,1) = ISET1
               NIND12(ISET12,2) = ISET2
            END IF
  200    CONTINUE
  100 CONTINUE
      NORT12(1) = ISET12
      IF (ITYPE .EQ. 12) OLDMAX = SRFMAX
C
      IJALL = NORB1*NORB2
      IF (TCON12) IJALL = NORB1*(NORB1 + 1)/2
      DO 600 I = 1, NODC12
         NTOT = NUC1*NUC2
         CALL REDSTA(0,NUCT12(I),NTOT,ISET12,IJALL)
  600 CONTINUE
C
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI, 1000) IELCTR
         WRITE (LUPRI, '(A,2I5)')   ' NUC   ', NUC1, NUC2
         WRITE (LUPRI, '(A,2I5)')   ' NCONT ', NORB1, NORB2
         WRITE (LUPRI, '(A,2I5)')   ' NSET  ', NSET1, NSET2
         WRITE (LUPRI, '(A,8I5)')  ' NUC12 ', (NUCT12(I),I=1,NODC12)
      END IF
 1000 FORMAT(//,' <<<<<<<<<< SUBROUTINE ODCSEG >>>>>>>>>>',
     *       //,' Overlap distribution arrays for electron ',I1,/)
      RETURN
      END
C  /* Deck odcoef */
      SUBROUTINE ODCOEF(COEF12,COOR12,EXP12,WORK,LWORK,JMAX1,JMAX2,
     &                  NHKT1,NHKT2,NSET1,NSET2,NUC1,NUC2,NUC12,MXUC12,
     &                  NORB1,NORB2,NPCO1,NPCO2,NUCS12,JSTR1,JSTR2,
     &                  SIGN1X,SIGN1Y,SIGN1Z,SIGN2X,SIGN2Y,SIGN2Z,COR1X,
     &                  COR1Y,COR1Z,COR2X,COR2Y,COR2Z,SAM12X,SAM12Y,
     &                  SAM12Z,I0X,I0Y,I0Z,ONECEN,DO1,DO2,BIGVEC,UNDIFF,
     &                  LONDON,SPNORB,DIA2SO,ZFS2EL,ITYPE,THRESH,
     &                  MAXDER,IPRINT)
C
C     TUH Apr 11 1988
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      LOGICAL SAM12X, SAM12Y, SAM12Z, ONECEN, DO1, DO2, BIGVEC, UNDIFF,
     &        SPNORB, DIA2SO, LONDON, DIA2UP, EL2CD, ZFS2EL
      DIMENSION COEF12(*), NPCO1(*), NPCO2(*), JSTR1(*), JSTR2(*),
     &          COOR12(*), EXP12(*), NUCS12(*), ORIGO(3), WORK(LWORK)
#include <orgcom.h>
#include <twosta.h>
C
      DIA2UP = DIA2SO .AND. ITYPE.EQ.12
      EL2CD  = DIA2SO .AND. ITYPE.EQ.34
      IF (TKTIME) TIMSTR = SECOND()
C
      IF (LONDON) THEN
         ORIGO(1) = ORIGIN(1)
         ORIGO(2) = ORIGIN(2)
         ORIGO(3) = ORIGIN(3)
      ELSE IF (DIA2UP) THEN
         ORIGO(1) = GAGORG(1)
         ORIGO(2) = GAGORG(2)
         ORIGO(3) = GAGORG(3)
      END IF
C
      IF (BIGVEC) THEN
         CALL CORDIF(NSET1,NSET2,THRESH,SAM12X,SAM12Y,SAM12Z,IPRINT,
     &               NPCO1,NPCO2,JSTR1,JSTR2)
         IF (LONDON) THEN
            CALL QUIT('BIGVEC not allowed with London orbitals.')
         END IF
      ELSE
         DIFX  = SIGN1X*COR1X - SIGN2X*COR2X
         DIFY  = SIGN1Y*COR1Y - SIGN2Y*COR2Y
         DIFZ  = SIGN1Z*COR1Z - SIGN2Z*COR2Z
         SAM12X = ABS(DIFX) .LT. THRESH
         SAM12Y = ABS(DIFY) .LT. THRESH
         SAM12Z = ABS(DIFZ) .LT. THRESH
         IF (LONDON.OR.DIA2UP) THEN
            SAM12X = SAM12X .AND. ABS(COR1X - ORIGO(1)) .LT. THRESH
            SAM12Y = SAM12Y .AND. ABS(COR1Y - ORIGO(2)) .LT. THRESH
            SAM12Z = SAM12Z .AND. ABS(COR1Z - ORIGO(3)) .LT. THRESH
         END IF
      END IF
      I0X = 0
      I0Y = 0
      I0Z = 0
      IF (SAM12X) I0X = 1
      IF (SAM12Y) I0Y = 1
      IF (SAM12Z) I0Z = 1
C
C     Work space allocation
C
      KHEXPI = 1
      KCOOR  = KHEXPI + NUC12
      IF (LONDON.OR.DIA2UP) THEN
         KLAST = KCOOR + MXUC12
      ELSE
         KLAST = KCOOR
      END IF
      LWRK   = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ODCOEF',' ',KLAST,LWORK)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
      CALL EXCOEF(COEF12,JMAX1,JMAX2,COOR12,WORK(KCOOR),EXP12,
     &            WORK(KHEXPI),ORIGO,BIGVEC,ONECEN,UNDIFF,
     &            NHKT1,NHKT2,NUC1,
     &            NUC2,NUC12,MXUC12,I0X,I0Y,I0Z,DIFX,DIFY,DIFZ,THRESH,
     &            ITYPE,MAXDER,LONDON,SPNORB,DIA2UP,EL2CD,ZFS2EL,
     &            DO1,DO2,IPRINT,WORK(KLAST),LWRK)
      LWTOT  = LWTOT - KLAST
      IF (TKTIME) TEXCOE = TEXCOE + SECOND() - TIMSTR
      RETURN
      END
C  /* Deck excoef */
      SUBROUTINE EXCOEF(COEF12,JMAX1,JMAX2,COOR12,COOR1,EXP12,HEXPPI,
     &                  ORIGO,
     &                  BIGVEC,ONECEN,UNDIFF,NHKT1,NHKT2,NUC1,NUC2,
     &                  NUC12,MXUC12,I120X,I120Y,I120Z,DIF12X,DIF12Y,
     &                  DIF12Z,THRESH,ITYPE,MAXDER,LONDON,SPNORB,DIA2UP,
     &                  EL2CD,ZFS2EL,DERECC,DERECD,IPRINT,WORK,LWORK)
C
C     TUH 84
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      LOGICAL YEQX, ZEQX, ZEQY, ONECEN, DER1, DER2, DER12,
     &        DERECC, DERECD, BIGVEC, UNDIFF, SPNORB, DIA2UP, EL2CD,
     &        LONDON, ZFS2EL
      DIMENSION COEF12(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,*),
     &          COOR12(NUC1*NUC2,3,3), HEXPPI(*), EXP12(NUC1*NUC2,3),
     &          COOR1(MXUC12),  ORIGO(3), WORK(LWORK)
#include <twosta.h>
#include <expopt.h>
#include <drw2el.h>
#include <r12int.h>
C
      LA = NHKT1 - 1
      LB = NHKT2 - 1
C     Space allocation for R12 method (WK/UniKA/15-11-2002).
      IF (U12INT) THEN
         KAOVRP = 1
         KBOVRP = KAOVRP + NUC12
         KLAST  = KBOVRP + NUC12
         IF (KLAST .GT. LWORK) CALL STOPIT('EXCOEF',' ',KLAST,LWORK)
         MWTOT  = MAX(MWTOT,LWTOT + KLAST)
      END IF
      IF (EXPGRA) THEN
         MAXDIF = 2
      ELSE IF (SPNORB .OR. BPH2OO) THEN
         MAXDIF = 1
      ELSE IF (FINDPT.OR.DPTINT) THEN
         MAXDIF = 2
      ELSE IF (DIA2UP .OR. EL2CD) THEN
         MAXDIF = 1
      ELSE IF (LONDON) THEN 
         MAXDIF = MAXDER
      ELSE IF (ZFS2EL) THEN
         IF (ITYPE .EQ. 12) THEN
            MAXDIF = MAXDER
         ELSE
            MAXDIF = 0
         END IF
      ELSE IF (ONECEN) THEN
         MAXDIF = 0
      ELSE
         MAXDIF = MAXDER
      END IF
      MAX1 = LA + MAXDIF
      MAX2 = LB + MAXDIF
      IF (LONDON .OR. DIA2UP) MAX2 = LB
      IF (ITYPE .EQ. 12) THEN
         SIGN  = D1
         DER1  = .TRUE.
         DER2  = .TRUE.
         DER12 = .TRUE.
      ELSE
         SIGN  = - D1
         DER1  = DERECC
         DER2  = DERECD
         DER12 = DER1 .AND. DER2
      END IF
      YEQX = (ABS(DIF12X-DIF12Y).LT.THRESH).AND..NOT.(BIGVEC.OR.LONDON)
     &       .AND. (I120X .EQ. I120Y)
      ZEQX = (ABS(DIF12X-DIF12Z).LT.THRESH).AND..NOT.(BIGVEC.OR.LONDON)
     &       .AND. (I120X .EQ. I120Z)
      ZEQY = (ABS(DIF12Y-DIF12Z).LT.THRESH).AND..NOT.(BIGVEC.OR.LONDON)
     &       .AND. (I120Y .EQ. I120Z)
      IF (DIA2UP) THEN
         YEQX = .FALSE.
         ZEQX = .FALSE.
         ZEQY = .FALSE.
      END IF
      LNG  = MXUC12*(JMAX1+JMAX2+1)*(JMAX1+1)*(JMAX2+1)
C
C     *************************
C     ***** Print section *****
C     *************************
C
      IF (IPRINT .GE. 10) THEN
         IELECT = 1
         IF (ITYPE .EQ. 34) IELECT = 2
         WRITE (LUPRI, 1000) IELECT
         WRITE (LUPRI, '(2X,A,2I5)')   ' NHKT1/2   ',NHKT1, NHKT2
         WRITE (LUPRI, '(2X,A, I5)')   ' NUC12     ',NUC12
         WRITE (LUPRI, '(2X,A, I5)')   ' MXUC12    ',MXUC12
         WRITE (LUPRI, '(2X,A,F12.6)') ' SIGN      ',SIGN
         WRITE (LUPRI, '(2X,A, L5)')   ' YEQX      ',YEQX
         WRITE (LUPRI, '(2X,A, L5)')   ' ZEQX      ',ZEQX
         WRITE (LUPRI, '(2X,A, L5)')   ' ZEQY      ',ZEQY
      END IF
C
C     *************************************************************
C     ***** Calculate undifferentiated expansion coefficients *****
C     *************************************************************
C
      IF ((.NOT. (UNDIFF .AND. (LA+LB .EQ. 0))) .OR. FINDPT
     *    .OR. DPTINT .OR. U12INT .OR. BPH2OO) THEN
C
         DO 100 I = 1, NUC12
            HEXPPI(I) = SIGN*DP5/EXP12(I,1)
  100    CONTINUE
C 
C        Calculate (a-b)/(a+b) factor (WK/UniKA/15-11-2002).
         IF (U12INT)
     &   CALL AMINBP(COEF12(1,0,0,0,1,3),NUC12,EXP12(1,2),EXP12(1,3),
     &               HEXPPI,SIGN,IPRINT)
C
C        x component
C        ===========
C
         CALL TWOODC(COEF12(1,0,0,0,1,1),JMAX1,JMAX2,NUC12,MXUC12,MAX1,
     &               MAX2,SIGN,I120X,COOR12(1,1,2),COOR12(1,1,3),HEXPPI,
     &               'EX00',IPRINT)
         IF (U12INT)
     &   CALL WKEQ30(1,COEF12(1,0,0,0,1,1),COEF12(1,0,0,0,1,2),
     &               COEF12(1,0,0,0,3,1),
     &               JMAX1,JMAX2,NUC12,MAX1,MAX2,SIGN,
     &               I120X,COOR12(1,1,2),COOR12(1,1,3),EXP12(1,2),
     &               EXP12(1,3),WORK(KAOVRP),WORK(KBOVRP),HEXPPI,
     &               'EX01',IPRINT)
c        argument list changed! C. Villani, Uni-Ka, Spring 2005
C
C        y component
C        ===========
C
         IF (YEQX) THEN
            CALL DCOPY(LNG,COEF12(1,0,0,0,1,1),1,COEF12(1,0,0,0,2,1),1)
            IF (U12INT)
     &      CALL DCOPY(LNG,COEF12(1,0,0,0,1,2),1,COEF12(1,0,0,0,2,2),1)
         ELSE
            CALL TWOODC(COEF12(1,0,0,0,2,1),JMAX1,JMAX2,NUC12,MXUC12,
     &                  MAX1,MAX2,SIGN,I120Y,COOR12(1,2,2),
     &                  COOR12(1,2,3),HEXPPI,'EY00',IPRINT)
            IF (U12INT)
     &      CALL WKEQ30(1,COEF12(1,0,0,0,2,1),COEF12(1,0,0,0,2,2),
     &                  COEF12(1,0,0,0,3,1),
     &                  JMAX1,JMAX2,NUC12,MAX1,MAX2,SIGN,
     &                  I120Y,COOR12(1,2,2),COOR12(1,2,3),EXP12(1,2),
     &                  EXP12(1,3),WORK(KAOVRP),WORK(KBOVRP),HEXPPI,
     &                  'EY01',IPRINT)
c           argument list changed! C. Villani, Uni-Ka, Spring 2005
         END IF
C
C        z component
C        ===========
C
         IF (ZEQX) THEN
            CALL DCOPY(LNG,COEF12(1,0,0,0,1,1),1,COEF12(1,0,0,0,3,1),1)
            IF (U12INT)
     &      CALL DCOPY(LNG,COEF12(1,0,0,0,1,2),1,COEF12(1,0,0,0,3,2),1)
         ELSE IF (ZEQY) THEN
            CALL DCOPY(LNG,COEF12(1,0,0,0,2,1),1,COEF12(1,0,0,0,3,1),1)
            IF (U12INT)
     &      CALL DCOPY(LNG,COEF12(1,0,0,0,2,2),1,COEF12(1,0,0,0,3,2),1)
         ELSE
            CALL TWOODC(COEF12(1,0,0,0,3,1),JMAX1,JMAX2,NUC12,MXUC12,
     &                  MAX1,MAX2,SIGN,I120Z,COOR12(1,3,2),
     &                  COOR12(1,3,3),HEXPPI,'EZ00',IPRINT)
            IF (U12INT)
     &      CALL WKEQ30(1,COEF12(1,0,0,0,3,1),COEF12(1,0,0,0,3,2),
     &                  COEF12(1,0,0,0,3,1),
     &                  JMAX1,JMAX2,NUC12,MAX1,MAX2,SIGN,
     &                  I120Z,COOR12(1,3,2),COOR12(1,3,3),EXP12(1,2),
     &                  EXP12(1,3),WORK(KAOVRP),WORK(KBOVRP),HEXPPI,
     &                  'EZ01',IPRINT)
c          argument list changed! C. Villani, Uni-Ka, Spring 2005
         END IF
      END IF
C
C     ******************************************************
C     ***** Expansion coefficients for London orbitals *****
C     ******************************************************
C
      IF (LONDON .OR. DIA2UP) THEN
         CALL DERLON(COEF12,COOR12,COOR1,ORIGO,NUC1,NUC2,NUC12,MXUC12,
     &               JMAX1,JMAX2,LA,LB,I120X,I120Y,I120Z,MAXDIF,IPRINT)
      ELSE IF (EL2CD) THEN
         CALL DER1AB(COEF12(1,0,0,0,1,1),COEF12(1,0,0,0,1,2),
     &               NUC12,MXUC12,JMAX1,JMAX2,
     &               EXP12(1,2),EXP12(1,3),LA,LB,'EX  ',I120X,IPRINT)
         CALL DER1AB(COEF12(1,0,0,0,2,1),COEF12(1,0,0,0,2,2),
     &               NUC12,MXUC12,JMAX1,JMAX2,
     &               EXP12(1,2),EXP12(1,3),LA,LB,'EY  ',I120Y,IPRINT)
         CALL DER1AB(COEF12(1,0,0,0,3,1),COEF12(1,0,0,0,3,2),
     &               NUC12,MXUC12,JMAX1,JMAX2,
     &               EXP12(1,2),EXP12(1,3),LA,LB,'EZ  ',I120Z,IPRINT)
C
C     ********************************************************
C     ***** Expansion coefficients for exponent gradient *****
C     ********************************************************
C
      ELSE IF (EXPGRA) THEN
         CALL DEREXP(COEF12(1,0,0,0,1,1),
     &               COEF12(1,0,0,0,1,2),COEF12(1,0,0,0,1,3),
     &               NUC12,MXUC12,JMAX1,JMAX2,
     &               EXP12(1,2),EXP12(1,3),LA,LB,'EX  ',I120X,IPRINT)
         CALL DEREXP(COEF12(1,0,0,0,2,1),
     &               COEF12(1,0,0,0,2,2),COEF12(1,0,0,0,2,3),
     &               NUC12,MXUC12,JMAX1,JMAX2,
     &               EXP12(1,2),EXP12(1,3),LA,LB,'EX  ',I120Y,IPRINT)
         CALL DEREXP(COEF12(1,0,0,0,3,1),
     &               COEF12(1,0,0,0,3,2),COEF12(1,0,0,0,3,3),
     &               NUC12,MXUC12,JMAX1,JMAX2,
     &               EXP12(1,2),EXP12(1,3),LA,LB,'EX  ',I120Z,IPRINT)
      ELSE
C
C
C        **************************************************************
C        **** Calculate first derivatives of overlap distributions ****
C        **************************************************************
C
         IF (MAXDIF .GT. 0) THEN
C
C           x component
C           ===========
C
            CALL DERONE(COEF12(1,0,0,0,1,1),COEF12(1,0,0,0,1,2),
     &                  COEF12(1,0,0,0,1,3),NUC12,MXUC12,JMAX1,JMAX2,
     &                  EXP12(1,2),EXP12(1,3),DER1,DER2,LA,LB,I120X,
     &                  'EX10','EX01',IPRINT)
C
C           y component
C           ===========
C
            IF (YEQX) THEN
               IF (DER1) CALL DCOPY(LNG,COEF12(1,0,0,0,1,2),1,
     &                                  COEF12(1,0,0,0,2,2),1)
               IF (DER2) CALL DCOPY(LNG,COEF12(1,0,0,0,1,3),1,
     &                                  COEF12(1,0,0,0,2,3),1)
            ELSE
               CALL DERONE(COEF12(1,0,0,0,2,1),COEF12(1,0,0,0,2,2),
     &                     COEF12(1,0,0,0,2,3),NUC12,MXUC12,JMAX1,JMAX2,
     &                     EXP12(1,2),EXP12(1,3),DER1,DER2,LA,LB,I120Y,
     &                     'EY10','EY01',IPRINT)
            END IF
C
C           z components
C           ============
C
            IF (ZEQX) THEN
               IF (DER1) CALL DCOPY(LNG,COEF12(1,0,0,0,1,2),1,
     &                                  COEF12(1,0,0,0,3,2),1)
               IF (DER2) CALL DCOPY(LNG,COEF12(1,0,0,0,1,3),1,
     &                                  COEF12(1,0,0,0,3,3),1)
            ELSE IF (ZEQY) THEN
               IF (DER1) CALL DCOPY(LNG,COEF12(1,0,0,0,2,2),1,
     &                                  COEF12(1,0,0,0,3,2),1)
               IF (DER2) CALL DCOPY(LNG,COEF12(1,0,0,0,2,3),1,
     &                                  COEF12(1,0,0,0,3,3),1)
            ELSE
               CALL DERONE(COEF12(1,0,0,0,3,1),COEF12(1,0,0,0,3,2),
     &                     COEF12(1,0,0,0,3,3),NUC12,MXUC12,JMAX1,JMAX2,
     &                     EXP12(1,2),EXP12(1,3),DER1,DER2,LA,LB,I120Z,
     &                     'EZ10','EZ01',IPRINT)
            END IF
         END IF
C
C        *************************************************************
C        *** Calculate second derivatives of overlap distributions ***
C        *************************************************************
C
         IF (MAXDIF .EQ. 2) THEN
C
C            x components
C            ============
C
               CALL DERTWO(COEF12(1,0,0,0,1,1),COEF12(1,0,0,0,1,4),
     &                     COEF12(1,0,0,0,1,5),COEF12(1,0,0,0,1,6),
     &                     NUC12,MXUC12,JMAX1,JMAX2,EXP12(1,2),
     &                     EXP12(1,3),DER1,DER2,DER12,LA,LB,I120X,'X',
     &                     IPRINT)
C
C           y components
C           ============
C
            IF (YEQX) THEN
               IF (DER1)  CALL DCOPY(LNG,COEF12(1,0,0,0,1,4),1,
     &                                   COEF12(1,0,0,0,2,4),1)
               IF (DER12) CALL DCOPY(LNG,COEF12(1,0,0,0,1,5),1,
     &                                   COEF12(1,0,0,0,2,5),1)
               IF (DER2)  CALL DCOPY(LNG,COEF12(1,0,0,0,1,6),1,
     &                                   COEF12(1,0,0,0,2,6),1)
            ELSE
               CALL DERTWO(COEF12(1,0,0,0,2,1),COEF12(1,0,0,0,2,4),
     &                     COEF12(1,0,0,0,2,5),COEF12(1,0,0,0,2,6),
     &                     NUC12,MXUC12,JMAX1,JMAX2,EXP12(1,2),
     &                     EXP12(1,3),DER1,DER2,DER12,LA,LB,I120Y,'Y',
     &                     IPRINT)
            END IF
C
C           z components
C           ============
C
            IF (ZEQX) THEN
               IF (DER1)  CALL DCOPY(LNG,COEF12(1,0,0,0,1,4),1,
     &                                   COEF12(1,0,0,0,3,4),1)
               IF (DER12) CALL DCOPY(LNG,COEF12(1,0,0,0,1,5),1,
     &                                   COEF12(1,0,0,0,3,5),1)
               IF (DER2)  CALL DCOPY(LNG,COEF12(1,0,0,0,1,6),1,
     &                                   COEF12(1,0,0,0,3,6),1)
            ELSE IF (ZEQY) THEN
               IF (DER1)  CALL DCOPY(LNG,COEF12(1,0,0,0,2,4),1,
     &                                   COEF12(1,0,0,0,3,4),1)
               IF (DER12) CALL DCOPY(LNG,COEF12(1,0,0,0,2,5),1,
     &                                   COEF12(1,0,0,0,3,5),1)
               IF (DER2)  CALL DCOPY(LNG,COEF12(1,0,0,0,2,6),1,
     &                                   COEF12(1,0,0,0,3,6),1)
            ELSE
               CALL DERTWO(COEF12(1,0,0,0,3,1),COEF12(1,0,0,0,3,4),
     &                     COEF12(1,0,0,0,3,5),COEF12(1,0,0,0,3,6),
     &                     NUC12,MXUC12,JMAX1,JMAX2,EXP12(1,2),
     &                     EXP12(1,3),DER1,DER2,DER12,LA,LB,I120Z,'Z',
     &                     IPRINT)
            END IF
         END IF
      END IF
      RETURN
 1000 FORMAT(//,' <<<<<<<<<< SUBROUTINE EXCOEF >>>>>>>>>>',
     &       //,' Overlap distribution coefficients for electron ',I1,/)
      END
C  /* Deck twoodc */
      SUBROUTINE TWOODC(COEF12,JMAX1,JMAX2,NUC12,MXUC12,MAX1,MAX2,SIGN,
     &                  I120,DIFPA,DIFPB,HEXPPI,WORD,IPRINT)
C
C     TUH 84
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 =0.00 D00, D1 =1.00 D00, D2 =2.00 D00)
      INTEGER T, A, B, AB
      CHARACTER WORD*4
      DIMENSION DIFPA(*), DIFPB(*), HEXPPI(*)
      DIMENSION COEF12(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2)
#include <ibtfun.h>
C
C     ****************************
C     ********** AB > 0 **********
C     ****************************
C
      IF (I120 .EQ. 0) THEN
C
C        ***** RUN OVER A *****
C
         DO 100 A = 0, MAX1
C
C           ***** E(0,0) *****
C
            IF (A .EQ. 0) THEN
               DO 200 I = 1, NUC12
                  COEF12(I,0,0,0) = D1
  200          CONTINUE
C
C           ***** E(1,0) *****
C
            ELSE IF (A .EQ. 1) THEN
               DO 300 I = 1, NUC12
                  COEF12(I,0,1,0) = DIFPA(I)
                  COEF12(I,1,1,0) = HEXPPI(I)
  300          CONTINUE
C
C           ***** E(2,0) *****
C
            ELSE IF (A .EQ. 2) THEN
               DO 400 I = 1, NUC12
                  DIFPAI = DIFPA(I)
                  EXPPIH = HEXPPI(I)
                  COEF12(I,0,2,0) = DIFPAI*DIFPAI + SIGN*EXPPIH
                  COEF12(I,1,2,0) = D2*DIFPAI*EXPPIH
                  COEF12(I,2,2,0) = EXPPIH*EXPPIH
  400          CONTINUE
C
C            ***** E(A,0) *****
C
            ELSE
               DO 500 I = 1, NUC12
                  DIFPAI = DIFPA(I)
                  EXPPIH = HEXPPI(I)
                  COEF12(I,  0,A,0) = DIFPAI*COEF12(I,  0,A-1,0)
     &                                + SIGN*COEF12(I,  1,A-1,0)
                  COEF12(I,A-1,A,0) = EXPPIH*COEF12(I,A-2,A-1,0)
     &                              + DIFPAI*COEF12(I,A-1,A-1,0)
                  COEF12(I,  A,A,0) = EXPPIH*COEF12(I,A-1,A-1,0)

  500          CONTINUE
               DO 510 T = 1, A - 2
                  T1 = SIGN*(T + 1.0D0)
                  DO 520 I = 1, NUC12
                     COEF12(I,T,A,0) = HEXPPI(I)*COEF12(I,T-1,A-1,0)
     *                                + DIFPA(I)*COEF12(I,  T,A-1,0)
     *                                      + T1*COEF12(I,T+1,A-1,0)
  520             CONTINUE
  510          CONTINUE
            END IF
C
C           ***** RUN OVER B *****
C
            DO 600 B = 1, MAX2
               AB = A + B
C
C              ***** E(0,1) *****
C
               IF (AB .EQ. 1) THEN
                  DO 700 I = 1, NUC12
                     COEF12(I,0,0,1) = DIFPB(I)
                     COEF12(I,1,0,1) = HEXPPI(I)
  700             CONTINUE
               ELSE IF (AB .EQ. 2) THEN
C
C                 ***** E(0,2) *****
C
                  IF (A .EQ. 0) THEN
                     DO 800 I = 1, NUC12
                        DIFPBI = DIFPB(I)
                        EXPPIH = HEXPPI(I)
                        COEF12(I,0,0,2) = DIFPBI*DIFPBI + SIGN*EXPPIH
                        COEF12(I,1,0,2) = D2*DIFPBI*EXPPIH
                        COEF12(I,2,0,2) = EXPPIH*EXPPIH
  800                CONTINUE
C
C                 ***** E(1,1) *****
C
                  ELSE
                     DO 810 I = 1, NUC12
                        DIFPAI = DIFPA(I)
                        DIFPBI = DIFPB(I)
                        EXPPIH = HEXPPI(I)
                        COEF12(I,0,1,1) = DIFPAI*DIFPBI + SIGN*EXPPIH
                        COEF12(I,1,1,1) = (DIFPAI + DIFPBI)*EXPPIH
                        COEF12(I,2,1,1) = EXPPIH*EXPPIH
  810                CONTINUE
                  END IF
C
C              ***** E(A,B) *****
C
               ELSE
                  DO 900 I = 1, NUC12
                     DIFPBI = DIFPB(I)
                     EXPPIH = HEXPPI(I)
                     COEF12(I,   0,A,B) = DIFPBI*COEF12(I,   0,A,B-1)
     &                                    + SIGN*COEF12(I,   1,A,B-1)
                     COEF12(I,AB-1,A,B) = EXPPIH*COEF12(I,AB-2,A,B-1)
     &                                  + DIFPBI*COEF12(I,AB-1,A,B-1)
                     COEF12(I,  AB,A,B) = EXPPIH*COEF12(I,AB-1,A,B-1)
  900             CONTINUE
                  DO 910 T = 1, AB - 2
                     T1 = SIGN*(T + 1.0D0)
                     DO 920 I = 1, NUC12
                        COEF12(I,T,A,B)=HEXPPI(I)*COEF12(I,T-1,A,B-1)
     &                                 + DIFPB(I)*COEF12(I,  T,A,B-1)
     &                                       + T1*COEF12(I,T+1,A,B-1)
  920               CONTINUE
  910             CONTINUE
               END IF
  600       CONTINUE
  100    CONTINUE
C
C     ****************************
C     ********** AB = 0 **********
C     ****************************
C
      ELSE
C
C        ***** RUN OVER A *****
C
         DO 105 A = 0, MAX1
C
C           ***** E(0,0) *****
C
            IF (A .EQ. 0) THEN
               DO 205 I = 1, NUC12
                  COEF12(I,0,0,0) = D1
  205          CONTINUE
C
C           ***** E(1,0) *****
C
            ELSE IF (A .EQ. 1) THEN
               DO 305 I = 1, NUC12
                  COEF12(I,1,1,0) = HEXPPI(I)
  305          CONTINUE
C
C           ***** E(2,0) *****
C
            ELSE IF (A .EQ. 2) THEN
               DO 405 I = 1, NUC12
                  EXPPIH = HEXPPI(I)
                  COEF12(I,0,2,0) = SIGN*EXPPIH
                  COEF12(I,2,2,0) = EXPPIH*EXPPIH
  405          CONTINUE
C
C            ***** E(A,0) *****
C
            ELSE
               IF (IBTAND(1,A) .EQ. 0) THEN
                  DO 505 I = 1, NUC12
                     COEF12(I,0,A,0) =      SIGN*COEF12(I,  1,A-1,0)
                     COEF12(I,A,A,0) = HEXPPI(I)*COEF12(I,A-1,A-1,0)
  505             CONTINUE
               ELSE
                  DO 506 I = 1, NUC12
                     COEF12(I,A,A,0) = HEXPPI(I)*COEF12(I,A-1,A-1,0)
  506             CONTINUE
               END IF
               DO 515 T = 2 - IBTAND(1,A), A - 2, 2
                  T1 = SIGN*(T + 1.0D0)
                  DO 525 I = 1, NUC12
                     COEF12(I,T,A,0) = HEXPPI(I)*COEF12(I,T-1,A-1,0)
     &                                      + T1*COEF12(I,T+1,A-1,0)
  525             CONTINUE
  515          CONTINUE
            END IF
C
C           ***** RUN OVER B *****
C
            DO 605 B = 1, MAX2
               AB = A + B
C
C              ***** E(0,1) *****
C
               IF (AB .EQ. 1) THEN
                  DO 705 I = 1, NUC12
                     COEF12(I,1,0,1) = HEXPPI(I)
  705             CONTINUE
C
C              ***** E(1,1) AND E(0,2) *****
C
               ELSE IF (AB .EQ. 2) THEN
                  DO 805 I = 1, NUC12
                     EXPPIH = HEXPPI(I)
                     COEF12(I,0,A,B) = SIGN*EXPPIH
                     COEF12(I,2,A,B) = EXPPIH*EXPPIH
  805             CONTINUE
C
C              ***** E(A,B) *****
C
               ELSE
                  IF (IBTAND(1,AB) .EQ. 0) THEN
                     DO 905 I = 1, NUC12
                       COEF12(I, 0,A,B) =      SIGN*COEF12(I,   1,A,B-1)
                       COEF12(I,AB,A,B) = HEXPPI(I)*COEF12(I,AB-1,A,B-1)
  905                CONTINUE
                  ELSE
                     DO 906 I = 1, NUC12
                       COEF12(I,AB,A,B) = HEXPPI(I)*COEF12(I,AB-1,A,B-1)
  906                CONTINUE
                  END IF
                  DO 915 T = 2 - IBTAND(1,AB), AB - 2, 2
                     T1 = SIGN*(T + 1.0D0)
                     DO 925 I = 1, NUC12
                        COEF12(I,T,A,B) = HEXPPI(I)*COEF12(I,T-1,A,B-1)
     &                                         + T1*COEF12(I,T+1,A,B-1)
  925                CONTINUE
  915             CONTINUE
               END IF
  605       CONTINUE
  105    CONTINUE
      END IF
      IF (IPRINT .LT. 10) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      WRITE (LUPRI, 1000)
      WRITE (LUPRI, 1010) MAX1, MAX2
      WRITE (LUPRI, 1030) I120
      IF (IPRINT .LT. 20) RETURN
      DO 2000 A = 0, MAX1
         DO 2100 B = 0, MAX2
            DO 2200 T = 0, A + B
               IF (IBTAND(A + B - T,I120) .EQ. 0) THEN
                  WRITE (LUPRI, 1100) WORD, A, B, T
                  WRITE (LUPRI, 1130) (COEF12(I,T,A,B), I = 1, NUC12)
               END IF
 2200       CONTINUE
 2100    CONTINUE
 2000 CONTINUE
      RETURN
 1000 FORMAT (/,'  <<<<<<<<<< SUBROUTINE TWOODC >>>>>>>>>>',/)
 1010 FORMAT ('  MAX1/B:     ',2I7)
 1030 FORMAT ('  I120:     ',I7)
 1100 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
 1130 FORMAT(1X,6F12.8)
      END
C  /* Deck derone */
      SUBROUTINE DERONE(UCOEF,ACOEF,BCOEF,NUC12,MXUC12,JMAX1,JMAX2,AXP,
     &                  BXP,DER1,DER2,MAX1,MAX2,I120,WORDA,WORDB,IPRINT)
C
C     TUH 84
C     Rewritten Oct 91, tuh
C
#include <implicit.h>
#include <priunit.h>
      CHARACTER WORDA*4, WORDB*4
      LOGICAL DER1, DER2
      INTEGER A, B, AB, T
      DIMENSION UCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          ACOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          BCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          AXP(*), BXP(*)
#include <ibtfun.h>
C
      DO 100 A = 0, MAX1
         DO 200 B = 0, MAX2
            AB = A + B
            DO 300 T = IBTAND(AB + 1,I120), AB + 1, I120 + 1
C
C              Differentiation with respect to A
C              =================================
C
               IF (DER1) THEN
                  IF ((A .EQ. 0) .OR. (T .GE. AB)) THEN
                     DO 400 I = 1, NUC12
                        ACOEF(I,T,A,B) = AXP(I)*UCOEF(I,T,A+1,B)
  400                CONTINUE
                  ELSE
                     DO 410 I = 1, NUC12
                        ACOEF(I,T,A,B) = AXP(I)*UCOEF(I,T,A+1,B)
     &                                      - A*UCOEF(I,T,A-1,B)
  410                CONTINUE
                  END IF
               END IF
C
C              Differentiation with respect to B
C              =================================
C
               IF (DER2) THEN
                  IF ((B .EQ. 0) .OR. (T .GE. AB)) THEN
                     DO 500 I = 1, NUC12
                        BCOEF(I,T,A,B) = BXP(I)*UCOEF(I,T,A,B+1)
  500                CONTINUE
                  ELSE
                     DO 510 I = 1, NUC12
                        BCOEF(I,T,A,B) = BXP(I)*UCOEF(I,T,A,B+1)
     &                                      - B*UCOEF(I,T,A,B-1)
  510                CONTINUE
                  END IF
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER
     &       ('First-derivative expansion coefficients in DERONE',-1)
         WRITE (LUPRI,'(1X,A,2L5)') ' DER1/2 ', DER1, DER2
         WRITE (LUPRI,'(1X,A,2I5)') ' MAX1/2 ', MAX1, MAX2
         WRITE (LUPRI,'(1X,A, I5)') ' NUC12  ', NUC12
         DO 600 A = 0, MAX1
            DO 610 B = 0, MAX2
               DO 620 T = 0, A + B + 1
                  IF (IBTAND(A + B + 1 - T,I120) .EQ. 0) THEN
                    WRITE (LUPRI,1100) WORDA, A, B, T
                    WRITE(LUPRI,'(1X,6F12.8)')(ACOEF(I,T,A,B),I=1,NUC12)
                    WRITE (LUPRI,1100) WORDB, A, B, T
                    WRITE(LUPRI,'(1X,6F12.8)')(BCOEF(I,T,A,B),I=1,NUC12)
                  END IF
  620          CONTINUE
  610       CONTINUE
  600    CONTINUE
      END IF
      RETURN
 1100 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
      END
C  /* Deck dertwo */
      SUBROUTINE DERTWO(UCOEF,ACOEF,XCOEF,BCOEF,NUC12,MXUC12,JMAX1,
     &                  JMAX2,AXP,BXP,DER1,DER2,DER12,MAX1,MAX2,I120,
     &                  WORD,IPRINT)
C
C     Oct 91 tuh
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, AB, T
      CHARACTER WORD*1
      LOGICAL DER1, DER2, DER12, AGT0, BGT0
      DIMENSION UCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          ACOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          XCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          BCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          AXP(*), BXP(*)
#include <ibtfun.h>
C
      DO 100 A = 0, MAX1
         AGT0 = A .GT. 0
         XA_TRI = A*(A-1.0D0)
         DO 200 B = 0, MAX2
            BGT0 = B .GT. 0
            XB_TRI = B*(B-1.0D0)
            AB = A + B
            DO 300 T = IBTAND(AB,I120), AB + 2, I120 + 1
C
C              Differentiation with respect to A
C              =================================
C
               IF (DER1) THEN
                  IF ((T .LE. AB - 2) .AND. (A .GT. 1)) THEN
                     DO 400 I = 1, NUC12
                       ACOEF(I,T,A,B) = AXP(I)*(AXP(I)*UCOEF(I,T,A+2,B)
     &                               - (2.0D0*A+1.0D0)*UCOEF(I,T,A,  B))
     &                                + XA_TRI        *UCOEF(I,T,A-2,B)
  400                CONTINUE
                  ELSE IF (T .LE. AB) THEN
                     DO 410 I = 1, NUC12
                       ACOEF(I,T,A,B) = AXP(I)*(AXP(I)*UCOEF(I,T,A+2,B)
     &                               - (2.0D0*A+1.0D0)*UCOEF(I,T,A,  B))
  410                CONTINUE
                  ELSE
                     DO 420 I = 1, NUC12
                       ACOEF(I,T,A,B) = AXP(I)*AXP(I)*UCOEF(I,T,A+2,B)
  420                CONTINUE
                  END IF
               END IF
C
C              Differentiation with respect to B
C              =================================
C
               IF (DER2) THEN
                  IF ((T .LE. AB - 2) .AND. (B .GT. 1)) THEN
                     DO 500 I = 1, NUC12
                       BCOEF(I,T,A,B) = BXP(I)*(BXP(I)*UCOEF(I,T,A,B+2)
     &                               - (2.0D0*B+1.0D0)*UCOEF(I,T,A,B  ))
     &                                + XB_TRI        *UCOEF(I,T,A,B-2)
  500                CONTINUE
                  ELSE IF (T .LE. AB) THEN
                     DO 510 I = 1, NUC12
                       BCOEF(I,T,A,B) = BXP(I)*(BXP(I)*UCOEF(I,T,A,B+2)
     &                               - (2.0D0*B+1.0D0)*UCOEF(I,T,A,B  ))
  510                CONTINUE
                  ELSE
                     DO 520 I = 1, NUC12
                       BCOEF(I,T,A,B) = BXP(I)*BXP(I)*UCOEF(I,T,A,B+2)
  520                CONTINUE
                  END IF
               END IF
C
C              Differentiation with respect to A and B
C              =======================================
C
               IF (DER12) THEN
                  IF (AGT0 .AND. BGT0 .AND. T.LE.AB-2) THEN
                     DO 600 I = 1, NUC12
                       XCOEF(I,T,A,B) = AXP(I)*BXP(I)*UCOEF(I,T,A+1,B+1)
     &                              - AXP(I)*B*UCOEF(I,T,A+1,B-1)
     &                              - BXP(I)*A*UCOEF(I,T,A-1,B+1)
     &                              +    (A*B)*UCOEF(I,T,A-1,B-1)
  600                CONTINUE
                  ELSE IF (AGT0 .AND. BGT0 .AND. T.LE.AB) THEN
                     DO 610 I = 1, NUC12
                       XCOEF(I,T,A,B) = AXP(I)*BXP(I)*UCOEF(I,T,A+1,B+1)
     &                                - AXP(I)*B     *UCOEF(I,T,A+1,B-1)
     &                                - A     *BXP(I)*UCOEF(I,T,A-1,B+1)
  610                CONTINUE
                  ELSE IF (AGT0 .AND. T.LE.AB) THEN
                     DO 620 I = 1, NUC12
                       XCOEF(I,T,A,B) = AXP(I)*BXP(I)*UCOEF(I,T,A+1,B+1)
     &                                - A     *BXP(I)*UCOEF(I,T,A-1,B+1)
  620                CONTINUE
                  ELSE IF (BGT0 .AND. T.LE.AB) THEN
                     DO 630 I = 1, NUC12
                       XCOEF(I,T,A,B) = AXP(I)*BXP(I)*UCOEF(I,T,A+1,B+1)
     &                                - AXP(I)*B     *UCOEF(I,T,A+1,B-1)
  630                CONTINUE
                  ELSE
                     DO 640 I = 1, NUC12
                       XCOEF(I,T,A,B) = AXP(I)*BXP(I)*UCOEF(I,T,A+1,B+1)
  640                CONTINUE
                  END IF
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
C
C     Print section
C     =============
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER
     &       ('Second-derivative expansion coefficients in DERTWO',-1)
         WRITE (LUPRI,'(1X,2A)')    ' Direction:', WORD
         WRITE (LUPRI,'(1X,A,3L5)') ' DER    ', DER1, DER12, DER2
         WRITE (LUPRI,'(1X,A,2I5)') ' MAX1/2 ', MAX1, MAX2
         WRITE (LUPRI,'(1X,A, I5)') ' NUC12  ', NUC12
         DO 700 A = 0, MAX1
            DO 710 B = 0, MAX2
               DO 720 T = 0, A + B + 2
                  IF (IBTAND(A + B + 2 - T,I120) .EQ. 0) THEN
                     WRITE (LUPRI,1100) 'E20',WORD, A, B, T,
     &                                  'E11',WORD, A, B, T,
     &                                  'E02',WORD, A, B, T
                     WRITE (LUPRI,'(11X,45A)') ('-',I=1,45)
                     DO 730 I = 1, NUC12
                        WRITE (LUPRI,'(10X,3(F12.8,5X))')ACOEF(I,T,A,B),
     &                                                   XCOEF(I,T,A,B),
     &                                                   BCOEF(I,T,A,B)
  730                CONTINUE
                  END IF
  720          CONTINUE
  710       CONTINUE
  700    CONTINUE
      END IF
      RETURN
 1100 FORMAT (//,11X,3(A3,A1,'(',I1,',',I1,';',I1,')',6X))
      END
C  /* Deck derlon */
      SUBROUTINE DERLON(COEF12,COOR12,COOR1,ORIGO,
     &                  NUC1,NUC2,NUC12,MXUC12,
     &                  JMAX1,JMAX2,MAX1,MAX2,I120X,I120Y,I120Z,MAXDIF,
     &                  IPRINT)
C
C     tuh Sep 4 92
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      INTEGER A, B, ABM, T, X, IODD(3)
      DIMENSION COEF12(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,*),
     &          COOR12(NUC1*NUC2,3,3), COOR1(MXUC12), ORIGO(3)
#include <abainf.h>
#include <chrxyz.h>
#include <orgcom.h>
#include <ibtfun.h>
C
      IODD(1) = I120X
      IODD(2) = I120Y
      IODD(3) = I120Z
      MADD = 0
      IF (MAXDIF .EQ. 2) MADD = 1
      DO 100 X = 1, 3
         DO 110 I = 1, NUC12
            COOR1(I) = COOR12(I,X,1) - COOR12(I,X,2) - ORIGO(X)
  110    CONTINUE
C
C        First derivative
C
         DO 310 A = 0, MAX1 + MADD
            DO 410 B = 0, MAX2
               ABM = A + B + 1
               DO 510 T = IBTAND(ABM,IODD(X)), ABM, IODD(X) + 1
                  IF (T.LT.ABM .AND. IODD(X).EQ.0) THEN
                     DO 610 I = 1, NUC12
                        COEF12(I,T,A,B,X,2) = COEF12(I,T,A+1,B,X,1)
     &                             + COOR1(I)*COEF12(I,T,A,  B,X,1)
  610                CONTINUE
                  ELSE
                     DO 615 I = 1, NUC12
                        COEF12(I,T,A,B,X,2) = COEF12(I,T,A+1,B,X,1)
  615                CONTINUE
                  END IF
  510          CONTINUE
  410       CONTINUE
  310    CONTINUE
C
C        Second derivative
C
         IF (MAXDIF .EQ. 2) THEN
            DO 320 A = 0, MAX1
               DO 420 B = 0, MAX2
                  ABM = A + B + 2
                  DO 520 T = IBTAND(ABM,IODD(X)), ABM, IODD(X) + 1
                     IF (T.LT.ABM .AND. IODD(X).EQ.0) THEN
                        DO 620 I = 1, NUC12
                           COEF12(I,T,A,B,X,4) = COEF12(I,T,A+1,B,X,2)
     &                                + COOR1(I)*COEF12(I,T,A,  B,X,2)
  620                   CONTINUE
                     ELSE
                        DO 625 I = 1, NUC12
                           COEF12(I,T,A,B,X,4) = COEF12(I,T,A+1,B,X,2)
  625                   CONTINUE
                     END IF
  520             CONTINUE
  420          CONTINUE
  320       CONTINUE
         END IF
  100 CONTINUE
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('London expansion coefficients in DERLON',-1)
         WRITE (LUPRI,'(1X,A,2I5)') ' MAX1/2 ', MAX1, MAX2
         WRITE (LUPRI,'(1X,A, I5)') ' NUC12  ', NUC12
         WRITE (LUPRI,'(1X,A, I5)') ' MAXDIF ', MAXDIF
         DO 700 M = 1, MAXDIF
         DO 700 X = 1, 3
            IF (M .EQ. 1) THEN
               CALL HEADER('First derivative London coefficients ('
     &                      //CHRXYZ(-X)//' direction)',-1)
            ELSE
               CALL HEADER('Second derivative London coefficients ('
     &                      //CHRXYZ(-X)//' direction)',-1)
            END IF
            IF (M .EQ. 1) THEN
               MM = 2
            ELSE
               MM = 4
            END IF
            DO 800 A = 0, MAX1
            DO 800 B = 0, MAX2
            DO 800 T = 0, A + B + M
               WRITE (LUPRI,1000) 'E', A, B, T
               IF (IBTAND(A + B + M - T,IODD(X)) .EQ. 0) THEN
                  WRITE(LUPRI,'(1X,6F12.8)')
     &               (COEF12(I,T,A,B,X,MM),I=1,NUC12)
               END IF
  800       CONTINUE
  700    CONTINUE
      END IF
      RETURN
 1000 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
      END
C  /* Deck der1ab */
      SUBROUTINE DER1AB(UCOEF,ABCOEF,NUC12,MXUC12,JMAX1,JMAX2,AXP,
     &                  BXP,MAX1,MAX2,WORD,I120,IPRINT)
C
C     TUH 84
C     Rewritten Oct 91, tuh
C
#include <implicit.h>
#include <priunit.h>
      CHARACTER WORD*4
      INTEGER A, B, AB, T
      DIMENSION UCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          ABCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          AXP(*), BXP(*)
#include <ibtfun.h>
C
      DO 100 A = 0, MAX1
         DO 200 B = 0, MAX2
            AB = A + B
            DO 300 T = IBTAND(AB + 1,I120), AB + 1, I120 + 1
               IF ((A .EQ. 0) .OR. (T .GE. AB)) THEN
                  DO 400 I = 1, NUC12
                     ABCOEF(I,T,A,B) = AXP(I)*UCOEF(I,T,A+1,B)
  400             CONTINUE
               ELSE
                  DO 410 I = 1, NUC12
                     ABCOEF(I,T,A,B) = AXP(I)*UCOEF(I,T,A+1,B)
     &                               - A     *UCOEF(I,T,A-1,B)
  410             CONTINUE
               END IF
               IF ((B .EQ. 0) .OR. (T .GE. AB)) THEN
                  DO 500 I = 1, NUC12
                     ABCOEF(I,T,A,B) = ABCOEF(I,T,A,B) 
     &                               + BXP(I)*UCOEF(I,T,A,B+1)
  500             CONTINUE
               ELSE
                  DO 510 I = 1, NUC12
                     ABCOEF(I,T,A,B) = ABCOEF(I,T,A,B) 
     &                               + BXP(I)*UCOEF(I,T,A,B+1)
     &                               - B     *UCOEF(I,T,A,B-1)
  510             CONTINUE
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER
     &       ('First-derivative expansion coefficients in DER1AB',-1)
         WRITE (LUPRI,'(1X,A,2I5)') ' MAX1/2 ', MAX1, MAX2
         WRITE (LUPRI,'(1X,A, I5)') ' NUC12  ', NUC12
         DO 600 A = 0, MAX1
            DO 610 B = 0, MAX2
               DO 620 T = 0, A + B + 1
                  IF (IBTAND(A + B + 1 - T,I120) .EQ. 0) THEN
                   WRITE (LUPRI,1100) WORD, A, B, T
                   WRITE(LUPRI,'(1X,6F12.8)')(ABCOEF(I,T,A,B),I=1,NUC12)
                  END IF
  620          CONTINUE
  610       CONTINUE
  600    CONTINUE
      END IF
      RETURN
 1100 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
      END
C  /* Deck derexp */
      SUBROUTINE DEREXP(UCOEF,ACOEF,BCOEF,NUC12,MXUC12,JMAX1,JMAX2,AXP,
     &                  BXP,MAX1,MAX2,WORD,I120,IPRINT)
C
C     TUH 92
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP25 =0.25 D00, DP5 = 0.50 D0, D2 = 2.0 D0)
      CHARACTER WORD*4
      INTEGER A, B, AB, T
      DIMENSION UCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          ACOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          BCOEF(MXUC12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2),
     &          AXP(*), BXP(*)
#include <ibtfun.h>
C
      DO A = 0, MAX1
      DO B = 0, MAX2
         DO T = 0, A + B 
         IF (IBTAND(A + B - T,I120) .EQ. 0) THEN
            DO I = 1, NUC12
               AINV = D2*AXP(I)**(-1)
               BINV = D2*BXP(I)**(-1)
               ACOEF(I,T,A,B) = (DP5*A + DP25)*AINV*UCOEF(I,T,A,B)
               BCOEF(I,T,A,B) = (DP5*B + DP25)*BINV*UCOEF(I,T,A,B)
            END DO
         END IF
         END DO
      END DO
      END DO
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Scaled expansion coefficients in DEREXP',-1)
         WRITE (LUPRI,'(1X,A,2I5)') ' MAX1/2 ', MAX1, MAX2
         WRITE (LUPRI,'(1X,A, I5)') ' NUC12  ', NUC12
         DO A = 0, MAX1
         DO B = 0, MAX2
            DO T = 0, A + B
            IF (IBTAND(A + B - T,I120) .EQ. 0) THEN
               WRITE (LUPRI,1100) WORD, A, B, T
               WRITE(LUPRI,'(1X,6F12.8)')(ACOEF(I,T,A,B),I=1,NUC12)
               WRITE(LUPRI,'(1X,6F12.8)')(BCOEF(I,T,A,B),I=1,NUC12)
            END IF
            END DO
         END DO
         END DO
      END IF
      RETURN
 1100 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
      END
C  /* Deck excsgn */
      SUBROUTINE EXCSGN(COEFCD,JMAXC,JMAXD,NHKTC,NHKTD,NUCCD,MXUCCD,
     &                  ICD0X,ICD0Y,ICD0Z,MAXDER,LONDON,SPNORB,DIA2SO,
     &                  EXPGRA,ISYMT,IPRINT)
C
C     tuh March 92
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D1 = 1.0D0)
      INTEGER X, T, P
      LOGICAL DOCOMP(3), SPNORB, DIA2SO, LONDON, EXPGRA
      DIMENSION COEFCD(MXUCCD,0:JMAXC+JMAXD,0:JMAXC,0:JMAXD,3,*),
     &          SGNOLD(3), SGNNEW(3)
#include <symmet.h>
#include <drw2el.h>
#include <r12int.h>
#include <chrnos.h>
#include <chrxyz.h>
      SAVE SGNOLD, DOCOMP, MAXDIF
#include <ibtfun.h>
C
      XAND(I) = PT(IBTAND(ISYMAX(1,1),I))
      YAND(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IBTAND(ISYMAX(3,1),I))
C
      IF (ISYMT .EQ. 0) THEN
         SGNOLD(1) = D1
         SGNOLD(2) = D1
         SGNOLD(3) = D1
         DOCOMP(1) = ICD0X .EQ. 0
         DOCOMP(2) = ICD0Y .EQ. 0
         DOCOMP(3) = ICD0Z .EQ. 0
         IF (SPNORB .OR. DIA2SO .OR. BPH2OO) THEN
            MAXDIF = 1
         ELSE IF (FINDPT .OR. DPTINT) THEN
            MAXDIF = 2
         ELSE
            MAXDIF = MAXDER
         END IF
      ELSE
         SGNNEW(1) = XAND(ISYMT)
         SGNNEW(2) = YAND(ISYMT)
         SGNNEW(3) = ZAND(ISYMT)
         DO 100 X = 1, 3
         IF (SGNNEW(X) .NE. SGNOLD(X)) THEN
            IF (DOCOMP(X)) THEN
C
C              London orbitals
C
               IF (EXPGRA) THEN
                  DO I = 0, NHKTC + 1
                  DO J = 0, NHKTD + 1
                     DO T = 0, I + J
                     IF (MOD(I+J-T,2).EQ.1) THEN
                        DO P = 1, NUCCD
                           COEFCD(P,T,I,J,X,1) = -COEFCD(P,T,I,J,X,1)
                        END DO
                     END IF
                     END DO
                  END DO
                  END DO
C
                  DO I = 0, NHKTC - 1
                  DO J = 0, NHKTD - 1
                     DO T = 0, I + J
                     IF (MOD(I+J-T,2).EQ.1) THEN
                        DO P = 1, NUCCD
                           COEFCD(P,T,I,J,X,2) = -COEFCD(P,T,I,J,X,2)
                           COEFCD(P,T,I,J,X,3) = -COEFCD(P,T,I,J,X,3)
                        END DO
                     END IF
                     END DO
                  END DO
                  END DO
               ELSE IF (LONDON) THEN
                  DO 300 M = 0, MAXDIF
                     IF (M .EQ. 0) MM = 1
                     IF (M .EQ. 1) MM = 2
                     IF (M .EQ. 2) MM = 4
                     DO 310 I = 0, NHKTC - 1
                     DO 310 J = 0, NHKTD - 1
                        DO 320 T = 0, I + J + M
                        IF (MOD(I+J+M-T,2).EQ.1) THEN
                           DO 330 P = 1, NUCCD
                              COEFCD(P,T,I,J,X,MM)=-COEFCD(P,T,I,J,X,MM)
  330                      CONTINUE
                        END IF
  320                   CONTINUE
  310                CONTINUE
  300             CONTINUE
C
C              diamagnetic spin-orbit
C
               ELSE IF (DIA2SO) THEN
                  DO M = 0, MAXDIF
                     DO I = 0, NHKTC - 1
                     DO J = 0, NHKTD - 1
                        DO T = 0, I + J + M
                        IF (MOD(I+J+M-T,2).EQ.1) THEN
                           DO P = 1, NUCCD
                            COEFCD(P,T,I,J,X,M+1)=-COEFCD(P,T,I,J,X,M+1)
                           END DO
                        END IF
                        END DO
                     END DO
                     END DO
                  END DO
C
C              and all others
C
               ELSE
                  DO 400 I = 0, NHKTC - 1
                  DO 400 J = 0, NHKTD - 1
                     DO 410 T = 0, I + J
                     IF (MOD(I+J+T,2).EQ.1) THEN
                        DO 420 P = 1, NUCCD
                           COEFCD(P,T,I,J,X,1) = -COEFCD(P,T,I,J,X,1)
  420                   CONTINUE
                     END IF
  410                CONTINUE
  400             CONTINUE
C                 For [r12,t1+T2] integrals (WK/UniKA/15-11-2002).
                  IF (U12INT) THEN
                     DO 530 I = 0, NHKTC - 1
                     DO 530 J = 0, NHKTD - 1
                        DO 540 T = 0, I + J
                        IF (MOD(I+J+1-T,2).EQ.1) THEN
                           DO 550 P = 1, NUCCD
                              COEFCD(P,T,I,J,X,2) = -COEFCD(P,T,I,J,X,2)
  550                      CONTINUE
                        END IF
  540                   CONTINUE
  530                CONTINUE
                  END IF
                  IF (MAXDIF .GT. 0) THEN
                     DO 500 I = 0, NHKTC - 1
                     DO 500 J = 0, NHKTD - 1
                        DO 510 T = 0, I + J + 1
                        IF (MOD(I+J+1-T,2).EQ.1) THEN
                           DO 520 P = 1, NUCCD
                              COEFCD(P,T,I,J,X,2) = -COEFCD(P,T,I,J,X,2)
                              COEFCD(P,T,I,J,X,3) = -COEFCD(P,T,I,J,X,3)
  520                      CONTINUE
                        END IF
  510                   CONTINUE
  500                CONTINUE
                  END IF
                  IF (MAXDIF .GT. 1) THEN
                     DO 600 I = 0, NHKTC - 1
                     DO 600 J = 0, NHKTD - 1
                        DO 610 T = 0, I + J + 2
                        IF (MOD(I+J+2-T,2).EQ.1) THEN
                           DO 620 P = 1, NUCCD
                              COEFCD(P,T,I,J,X,4) = -COEFCD(P,T,I,J,X,4)
                              COEFCD(P,T,I,J,X,5) = -COEFCD(P,T,I,J,X,5)
                              COEFCD(P,T,I,J,X,6) = -COEFCD(P,T,I,J,X,6)
  620                      CONTINUE
                        END IF
  610                   CONTINUE
  600                CONTINUE
                  END IF
               END IF
            END IF
            SGNOLD(X) = SGNNEW(X)
         END IF
  100    CONTINUE
      END IF
C
C     Print Section
C     =============
C
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from EXCSGN','*',103)
         IDERIV = 0
         DO 700 ITYPE = 1, (MAXDIF + 1)*(MAXDIF + 2)/2
            IF (ITYPE .EQ. 2) IDERIV = IDERIV + 1
            IF (ITYPE .EQ. 4) IDERIV = IDERIV + 1
            IF (LONDON .AND. (ITYPE .NE. 1 .OR.
     &                        ITYPE .NE. 2 .OR.
     &                        ITYPE .NE. 4)) GO TO 700
            DO 800 X = 1, 3
               IF (X.EQ.1) ICD0 = ICD0X
               IF (X.EQ.2) ICD0 = ICD0Y
               IF (X.EQ.3) ICD0 = ICD0Z
               CALL HEADER('Expansion coefficients for direction '
     &                     //CHRXYZ(X)//' and differentiation level '
     &                     //CHRNOS(IDERIV)//':',-1)
               DO 900 I = 0, NHKTC - 1
                  DO 910 J = 0, NHKTD - 1
                     DO 920 T = 0, I + J + IDERIV
                     IF (IBTAND(I + J + IDERIV - T,ICD0) .EQ. 0) THEN
                        WRITE (LUPRI,1000) 'E', I, J, T
                        WRITE (LUPRI,'(1X,6F12.8)')
     &                     (COEFCD(P,T,I,J,X,ITYPE),P=1,NUCCD)
                     END IF
  920                CONTINUE
  910             CONTINUE
  900          CONTINUE
  800       CONTINUE
  700    CONTINUE
      END IF
      RETURN
 1000 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
      END
