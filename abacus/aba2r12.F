C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
C  /* Deck eriher */
      SUBROUTINE ABEQ52(Q,R,W,A,PQX,PQY,PQZ,INDHER,JMAX,
     &                  EXP12,EXP34,NUC12,NUC34,NUABCD,
     &                  NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
C
C     Written by Wim Klopper (University of Karlsruhe, 14 November 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DM2 = -2.0D0)
      DIMENSION Q(*), R(*), W(*), A(*), EXP12(NUC12), EXP34(NUC34),
     &          IODDHR(*), INDHER(*), PQX(*), PQY(*), PQZ(*)
      IODS = 0
      DO  IOD12 = 1, NUC12
         EXPP = EXP12(IOD12)
         DO  IOD34 = 1, NUC34
            IODS = IODS + 1
            EXPQ = EXP34(IOD34)
            W(IODS) = DM2*EXPP*EXPQ/(EXPP + EXPQ)
         END DO
      END DO
      CALL WKEQ52(Q,R,W,A,PQX,PQY,PQZ,INDHER,JMAX,NUABCD,
     &            NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
      RETURN
      END
C  /* Deck r00G */
      SUBROUTINE R00G(RJ000,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,PQX,
     &               PQY,PQZ,JMAX,NOINT,NUABCD,NUC1,NUC2,NUC12,NUC3,
     &               NUC4,NUC34,THRESH,ONECEN,IPRINT,IPQ0X,IPQ0Y,IPQ0Z,
     &               SIGNT,FACINT,HEXPP,HEXPQ)
C     Copy of R0001 for use with R12EIN (WK/UniKA/20-11-2002).
#include <implicit.h>
#include <priunit.h>
#include <subdir.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0, D2 = 2.D0, DP25 = 0.25D0)
      LOGICAL ONECEN, NOINT
      DIMENSION RJ000(NUABCD,0:JMAX),
     &          PQX(NUABCD), PQY(NUABCD), PQZ(NUABCD),
     &          COOR12(NUC1*NUC2,3), COOR34(NUC3*NUC4,3),
     &          EXP12(*), EXP34(*), FAC12 (*), FAC34(*), SIGNT(3),
     &          FACINT(*), HEXPP(*) ,HEXPQ(*)
      NOINT = .FALSE.
      IF (ONECEN) THEN
         CALL DZERO(PQX,NUABCD)
         CALL DZERO(PQY,NUABCD)
         CALL DZERO(PQZ,NUABCD)
         IPQ0X = 1
         IPQ0Y = 1
         IPQ0Z = 1
         IODS = NUABCD
         NODS = NUABCD
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
C	It turns out 64-bits IRIX compiler misinterprets this statement with
C       -O2 optimization. These fixes prevent it. K.Ruud-Jan96
#if defined (SYS_IRIX)
         IODS  = 1
         NODS  = 1
#else
         IODS  = 0
         NODS  = 0
#endif
         DO 300 IOD12 = 1, NUC12
            PX     = SGN12X*COOR12(IOD12,1)
            PY     = SGN12Y*COOR12(IOD12,2)
            PZ     = SGN12Z*COOR12(IOD12,3)
            DO 310 IOD34 = 1, NUC34
#if !defined (SYS_IRIX)
               IODS = IODS + 1
               NODS = NODS + 1
#endif
               PQXI = PX - SGN34X*COOR34(IOD34,1)
               PQYI = PY - SGN34Y*COOR34(IOD34,2)
               PQZI = PZ - SGN34Z*COOR34(IOD34,3)
               PQX(IODS) = PQXI
               PQY(IODS) = PQYI
               PQZ(IODS) = PQZI
#if defined (SYS_IRIX)
               IODS = IODS + 1
               NODS = NODS + 1
#endif
  310       CONTINUE
  300    CONTINUE
#if defined (SYS_IRIX)
	 IODS = IODS - 1
	 NODS = NODS - 1
#endif
         IPQ0X = 1
         IPQ0Y = 1
         IPQ0Z = 1
         IF (DASUM(NUABCD,PQX,1) .GT. THRESH) IPQ0X = 0
         IF (DASUM(NUABCD,PQY,1) .GT. THRESH) IPQ0Y = 0
         IF (DASUM(NUABCD,PQZ,1) .GT. THRESH) IPQ0Z = 0
      END IF
C
      IJ = 0
      DO IOD12 = 1, NUC12
         FAC = FAC12(IOD12)
         DO IOD34 = 1, NUC34
            IJ = IJ + 1
            FACINT(IJ) = FAC * FAC34(IOD34)
            HEXPP(IJ) = EXP12(IOD12)
            HEXPQ(IJ) = EXP34(IOD34)
         END DO
      END DO
      RETURN
      END
C  /* Deck r12wrt */
      SUBROUTINE R12WRT(BUF,LBUF,ICOUNT,ITYPE,INDA,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <ibtpar.h>
      DIMENSION BUF(LBUF,5)
#include <twosta.h>
#include <r12int.h>
#include <inftap.h>
#include <nuclei.h>
#include <symmet.h>
      INTEGER JBUF(5), MBUF
      SAVE JBUF, MBUF
#if defined (VAR_SPLITFILES)
#include <2gbdef.h>
      CHARACTER*80 FNNAME, FNNM2
#include <chrnos.h>
#endif
#include <ibtfun.h>
#include <sdpre.h>
C
       IF (ITYPE .EQ. -1) THEN
         IDUM = 0
         I = 1
         CALL GPOPEN(LUR12(I),'AOXYZ','UNKNOWN',' ',' ',IDUM,.FALSE.)
         CALL NEWLAB('AOXYZINT',LUR12(I),LUPRI)
         IF (V12INT) THEN
            I = I + 1
            CALL GPOPEN(LUR12(I),'AOV12','UNKNOWN',' ',' ',IDUM,.FALSE.)
            CALL NEWLAB('AOV12INT',LUR12(I),LUPRI)
         END IF
         IF (R12INT) THEN
            I = I + 1
            CALL GPOPEN(LUR12(I),'AOR12','UNKNOWN',' ',' ',IDUM,.FALSE.)
            CALL NEWLAB('AOR12INT',LUR12(I),LUPRI)
         END IF
         IF (U12INT) THEN
            I = I + 1
            CALL GPOPEN(LUR12(I),'AOU12','UNKNOWN',' ',' ',IDUM,.FALSE.)
            CALL NEWLAB('AOU12INT',LUR12(I),LUPRI)
         END IF
         IF (U21INT) THEN
            I = I + 1
            CALL GPOPEN(LUR12(I),'AOU21','UNKNOWN',' ',' ',IDUM,.FALSE.)
            CALL NEWLAB('AOU21INT',LUR12(I),LUPRI)
         END IF
         IF (INDA .NE. 0) THEN
            MBUF = LBUF*2 + 4/IRAT
         ELSE
            MBUF = LBUF*2 + 2/IRAT
         END IF
      END IF
C
      DO 100 KR12 = 1, NOPP12 + 1
       IF (ITYPE .EQ. -1) THEN
         ICOUNT = 0
         JBUF(KR12) = MBUF + 8
       ELSE IF (ITYPE .EQ. 0) THEN
#if defined (VAR_SPLITFILES)
         IF (JBUF(KR12) + MBUF .GT. I2GB) THEN
            INQUIRE(UNIT=LUR12(KR12),NAME=FNNAME)
            LN = 1
 10         CONTINUE
            IF (FNNAME(LN:LN) .NE. ' ') THEN
               LN = LN + 1
               GOTO 10
            END IF
            LN = LN - 1
            CALL GPCLOSE(LUR12(KR12),'KEEP')
            I = LN - 1
            IF (FNNAME(I:I) .NE. '-') THEN
               FNNM2 = FNNAME(1:LN)//'-0'
               LN = LN + 2
            ELSE
               READ(FNNAME(LN:),'(I1)') INUM
               INUM = INUM + 1
               IF (INUM .GT. 9) THEN
                  WRITE (LUPRI,'(/A)') ' DALTON needs to split a '//
     &                 ' file more than 11 times.',
     &                 ' This is currently not supported'
                  CALL QUIT('Too many splittings of a file')
               END IF
               FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
            END IF
            CALL GPOPEN(LUR12(KR12),FNNM2(1:LN),'UNKNOWN',' ',' ',
     &           IDUMMY,.FALSE.)
            JBUF(KR12) = MBUF
         END IF
#endif
         IF (INDA .NE. 0) WRITE (LUR12(KR12)) INDA
         WRITE (LUR12(KR12)) (BUF(L,KR12),L=1,LBUF),ICOUNT
         JBUF(KR12) = JBUF(KR12) + MBUF
         IF (KR12 .EQ. NOPP12 + 1) ICOUNT = 0
       ELSE
         IF (INDA .NE. 0) THEN
            WRITE (LUR12(KR12)) -INDA
            WRITE (LUR12(KR12)) (BUF(L,KR12),L=1,LBUF),ICOUNT
            JBUF(KR12) = JBUF(KR12) + MBUF
         ELSE
            IF (ICOUNT .GT. 0) THEN
               WRITE (LUR12(KR12)) (BUF(L,KR12),L=1,LBUF),ICOUNT
               JBUF(KR12) = JBUF(KR12) + MBUF
            END IF
            WRITE (LUR12(KR12)) (BUF(L,KR12),L=1,LBUF),-1
         END IF
         CALL GPCLOSE(LUR12(KR12),'KEEP')
       END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck rn2out */
      SUBROUTINE RN2OUT(SO,NSOINT,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,
     &                  THRESH,NINDAB,NINDCD,IPRINT)
C
C     Version of UN2OUT for R12 integrals (WK/26-11-2002).
C
#include <implicit.h>
#include <priunit.h>
#include <r12int.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
      PARAMETER (LBUF = 600)
      LOGICAL DCMPAB, DCMPCD, DCMPAC, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, IAEQIC, IALTIC, IPNTLG(3,*), NOTEST,
     &        GTTHRS
      DIMENSION SO(NSOINT,*), BUF(LBUF,5), IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
#include <nuclei.h>
#include <twocom.h>
#include <symmet.h>
      SAVE BUF, XBUF, ICOUNT
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine RN2OUT',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL?  ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB? ', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR? ', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
         WRITE (LUPRI,'(2X,A,2L5)') 'DIAGAB/CD', DIAGAB, DIAGCD
         WRITE (LUPRI,'(2X,A,2L5)') 'TCONAB/CD', TCONAB, TCONCD
         WRITE (LUPRI,'(2X,A,2L5)') 'SHAEQB/CD', SHAEQB, SHCEQD
         WRITE (LUPRI,'(2X,A, L5)') 'SHABAB', SHABAB
      END IF
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) CALL R12WRT(BUF,LBUF,ICOUNT,-1,0,IPRINT)
C
      ISOFF  = 0
      NBUFCL = 0
      NSTART = ICOUNT
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      DO 100 I = 1, NINTS
         NSTRNA = IPNTNO(1,I)
         NSTRNB = IPNTNO(2,I)
         NSTRNC = IPNTNO(3,I)
         NSTRND = IPNTNO(4,I)
         IREPA  = IPNTRP(1,I)
         IREPB  = IPNTRP(2,I)
         IREPC  = IPNTRP(3,I)
         IREPD  = IBTXOR(IBTXOR(IREPA,IREPB),IREPC)
         IF (NOTEST) THEN
               INT = 0
               DO 200 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  DO 210 ICD = 1, NORBCD
                     INT = INT + 1
                     IC = KHKTC*(NINDCD(ICD,1) - 1)
                     ID = KHKTD*(NINDCD(ICD,2) - 1)
                     INDC = IPTSYM(NSTRNC + IC,IREPC)
                     INDD = IPTSYM(NSTRND + ID,IREPD)
                     CALL LAB64U(SO(ISOFF+INT,1),NSOINT,
     &                           INDA,INDB,INDC,INDD,XABCD,
     &                           THRESH,GTTHRS,IPRINT)
                     IF (GTTHRS) THEN
                        ICOUNT = ICOUNT + 1
                        BUF(ICOUNT,1) = XABCD
                        DO IOPP = 1, NOPP12
                           BUF(ICOUNT,IOPP+1) = SO(ISOFF+INT,IOPP)
                        END DO
                        IF (ICOUNT .EQ. LBUF) THEN
                           NBUFCL = NBUFCL + 1
                           CALL R12WRT(BUF,LBUF,ICOUNT,0,0,IPRINT)
                        END IF
                     END IF
  210             CONTINUE
  200          CONTINUE
         ELSE
            DCMPAB = IPNTLG(1,I)
            DCMPCD = IPNTLG(2,I)
            DCABAB = IPNTLG(3,I)
            DRALTB = IREPA .LT. IREPB
            DRCLTD = IREPC .LT. IREPD
            DRABAB = DCABAB .AND. IREPA.EQ.IREPC .AND. IREPB.EQ.IREPD
            INT = 0
            DO 300 IAB = 1, NORBAB
               IA = KHKTA*(NINDAB(IAB,1) - 1)
               IB = KHKTB*(NINDAB(IAB,2) - 1)
               IF (DCMPAB) THEN
                  IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                     INT = INT + NORBCD
                     GO TO 300
                  END IF
               END IF
               INDA = IPTSYM(NSTRNA + IA,IREPA)
               INDB = IPTSYM(NSTRNB + IB,IREPB)
               DO 310 ICD = 1,NORBCD
                  IC = KHKTC*(NINDCD(ICD,1) - 1)
                  ID = KHKTD*(NINDCD(ICD,2) - 1)
                  INT = INT + 1
                  IF (DCMPCD ) THEN
                     IF (ID.GT.IC) GO TO 310
                     IF (DRCLTD .AND. ID.EQ.IC) GO TO 310
                  END IF
                  IF (DRABAB) THEN
                     IF (IA.LT.IC.OR.(IA.EQ.IC.AND.IB.LT.ID)) GOTO 310
                  END IF
                  INDC = IPTSYM(NSTRNC + IC,IREPC)
                  INDD = IPTSYM(NSTRND + ID,IREPD)
                  CALL LAB64U(SO(ISOFF+INT,1),NSOINT,
     &                        INDA,INDB,INDC,INDD,XABCD,
     &                        THRESH,GTTHRS,IPRINT)
                  IF (GTTHRS) THEN
                     ICOUNT = ICOUNT + 1
                     BUF(ICOUNT,1) = XABCD
                     DO IOPP = 1, NOPP12
                        BUF(ICOUNT,IOPP+1) = SO(ISOFF+INT,IOPP)
                     END DO
                     IF (ICOUNT. EQ. LBUF) THEN
                        CALL R12WRT(BUF,LBUF,ICOUNT,0,0,IPRINT)
                        NBUFCL = NBUFCL + 1
                     END IF
                  END IF
  310          CONTINUE
  300       CONTINUE
         END IF
         ISOFF = ISOFF + NOABCD
  100 CONTINUE
      NGINT = LBUF*NBUFCL + ICOUNT - NSTART
      CALL DELSTA(0,NGINT)
C
C     *************************************
C     ***** Last call to empty buffer *****
C     *************************************
C
      IF (LAST) CALL R12WRT(BUF,LBUF,ICOUNT,1,0,IPRINT)
      RETURN
      END
