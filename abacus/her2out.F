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
C  /* Deck un2out */
      SUBROUTINE UN2OUT(SO,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,THRESH,
     &                  NINDAB,NINDCD,IPRINT)
C
C     Write out blocks of symmetry integrals, eliminating duplicates
C
C                                          880412   PRT
C
C     Some low-brain work has been done by TUH
C
C     Rewritten to allow for triangular looping and to eliminate
C     all duplicates by testing  880601 TUH
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
      PARAMETER (LBUF = 600)
      LOGICAL DCMPAB, DCMPCD, DCMPAC, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, IAEQIC, IALTIC, IPNTLG(3,*), NOTEST
      DIMENSION SO(*), BUF(LBUF), IBUF(LBUF,2),
     &          IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
#include <nuclei.h>
#include <twocom.h>
#include <symmet.h>
      SAVE BUF, IBUF, ICOUNT
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine UN2OUT',-1)
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
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,-1,NBITS,0,
     &                       IPRINT)
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
            IF (NIBUF .EQ. 1) THEN
               INT = 0
               DO 200 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
                  DO 210 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1)=MAX(INDAB,INDCD)*IBIT2+INDAB+INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                       END IF
                     END IF
  210             CONTINUE
  200          CONTINUE
            ELSE
               INT = 0
               DO 205 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
                  DO 215 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1) = MAX(INDAB,INDCD)
                       IBUF(ICOUNT,2) = MIN(INDAB,INDCD)
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                       END IF
                     END IF
  215             CONTINUE
  205          CONTINUE
            END IF
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
               INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
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
                  SOINT = SO(ISOFF+INT)
                  IF (ABS(SOINT) .GT. THRESH) THEN
                     IF (NIBUF .EQ. 1) THEN
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1)=MAX(INDAB,INDCD)*IBIT2+INDAB+INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                          NBUFCL = NBUFCL + 1
                       END IF
                     ELSE
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1) = MAX(INDAB,INDCD)
                       IBUF(ICOUNT,2) = MIN(INDAB,INDCD)
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                          NBUFCL = NBUFCL + 1
                       END IF
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
      IF (LAST) THEN
         WRITE (LUPRI,'(/A,1P,D10.2)') ' Threshold for neglecting '//
     &      'two-electron integrals:',THRESH
         CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,1,NBITS,0,IPRINT)
      END IF
      RETURN
      END
C  /* Deck un2wrt */
      SUBROUTINE UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,ITYPE,NBITS,INDA,
     &                  IPRINT)
C
C     Write undifferentiated 2-electron integrals to disk file AOTWOINT
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <ibtpar.h>
      DIMENSION BUF(LBUF), IBUF(LBUF,NIBUF)
#include <twosta.h>
#include <inftap.h>
#include <nuclei.h>
#include <symmet.h>
      SAVE NBUF, JBUF
#if defined (VAR_SPLITFILES)
#include <2gbdef.h>
      CHARACTER*80 FNNAME, FNNM2
#include <chrnos.h>
#endif
#include <ibtfun.h>
C
      IF (ITYPE .EQ. -1) THEN
         CALL REWSPL(LUINTA)
         CALL NEWLAB('BASINFO ',LUINTA,LUPRI)
         WRITE (LUINTA) MAXREP+1,(NAOS(I),I=1,8),LBUF,NIBUF,NBITS
         CALL NEWLAB('BASTWOEL',LUINTA,LUPRI)
         JLEN = 608
         ICOUNT = 0
         NBUF = 0
         JBUF = 608
      ELSE IF (ITYPE .EQ. 0) THEN
#if defined (VAR_SPLITFILES)
C
         IF ((JBUF + LBUF*(IRAT + NIBUF) + 1) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
            INQUIRE(UNIT=LUINTA,NAME=FNNAME)
            LN = 1
 10         CONTINUE
            IF (FNNAME(LN:LN) .NE. ' ') THEN
               LN = LN + 1
               GOTO 10
            END IF
            LN = LN - 1
            CALL GPCLOSE(LUINTA,'KEEP')
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
            CALL GPOPEN(LUINTA,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &           .FALSE.)
            JBUF = 0
         END IF
         IF (INDA .NE. 0) THEN
            WRITE (LUINTA) INDA
            JBUF = JBUF + 1
         END IF
         WRITE (LUINTA) BUF,IBUF,ICOUNT
         JBUF = JBUF + LBUF*(IRAT + NIBUF) + 1
#else
         IF (INDA .NE. 0) WRITE (LUINTA) INDA
         WRITE (LUINTA) BUF,IBUF,ICOUNT
#endif
         NBUF = NBUF + 1
C
         IF (IPRINT .GE. 6) THEN
            WRITE (LUPRI,'(2X,A,I5,A/)')
     &         'Integral buffer #',NBUF,' has been written.'
            IBIT1 = 2**NBITS - 1
            DO 100 INT = 1, ICOUNT
               IF (NIBUF .EQ. 1) THEN
                  I = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
                  J = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
                  K = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
                  L = IBTAND(       IBUF(INT,1),         IBIT1)
               ELSE
                  I = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
                  J = IBTAND(       IBUF(INT,1),         IBIT1)
                  K = IBTAND(IBTSHR(IBUF(INT,2),  NBITS),IBIT1)
                  L = IBTAND(       IBUF(INT,2),         IBIT1)
               END IF
               WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8)')
     &                      ' ## ', I, J, K, L, BUF(INT)
  100       CONTINUE
         END IF
         ICOUNT = 0
      ELSE
         IF (INDA .NE. 0) THEN
            WRITE (LUINTA) -INDA
            WRITE (LUINTA) BUF,IBUF,ICOUNT
         ELSE
            IF (ICOUNT .GT. 0) WRITE (LUINTA) BUF,IBUF,ICOUNT
            WRITE (LUINTA) BUF,IBUF,-1
         END IF
C
         IF (IPRINT .GE. 6) THEN
            IF (ICOUNT .GT. 0) THEN
               WRITE (LUPRI,'(2X,A,I5,A/)')
     &            'Integral buffer #',NBUF+1,' has been written.'
               IBIT1 = 2**NBITS - 1
               DO 200 INT = 1, ICOUNT
                  IF (NIBUF .EQ. 1) THEN
                     I = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
                     J = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
                     K = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
                     L = IBTAND(       IBUF(INT,1),         IBIT1)
                  ELSE
                     I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
                     J = IBTAND(       IBUF(INT,1),       IBIT1)
                     K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
                     L = IBTAND(       IBUF(INT,2),       IBIT1)
                  END IF
                  WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8)')
     &                         ' ## ', I, J, K, L, BUF(INT)
  200          CONTINUE
            END IF
         END IF
C
C        Statistics
C
         IF (INDA .EQ. 0) THEN
C           N2WRIT = LBUF*NBUF + ICOUNT
C           hjaaj: do all intermediates with integer*8 N2WRIT
C                  to avoid integer*4 overflows for big basis sets
            N2WRIT = LBUF
            N2WRIT = N2WRIT*NBUF
            N2WRIT = N2WRIT + ICOUNT
C
            IF (ICOUNT.GT.0 .AND. INDA.NE.0) THEN
               NBUF = NBUF + 2
            ELSE
               NBUF = NBUF + 1
            END IF
            IF (IRAT .EQ. 1) LWORD = 8
            IF (IRAT .EQ. 2) LWORD = 4
Chjaaj:     FBYTES is calculed in four lines, to avoid integer overflows.
            FBYTES = LBUF*IRAT + NIBUF*LBUF + 1
            FBYTES = NBUF*FBYTES
            FBYTES = LWORD*FBYTES
            FBYTES = FBYTES / (1024.D0**2)
            FNALL  = (NBASIS*(NBASIS + 1))/2
            FNALL  = (FNALL*(FNALL + 1.0D0))*0.5D0
            PERCNT = N2WRIT
            PERCNT = 100.0D0*PERCNT / FNALL
            WRITE (LUPRI,'(A,I12,A,F5.1,A,/A,F12.3/)')
     &         ' Number of two-electron integrals written:',N2WRIT,
     &         ' (',PERCNT,'% )',
     &         ' Megabytes written:                       ',FBYTES
         END IF
         CALL GPCLOSE(LUINTA,'KEEP')
      END IF
      RETURN
      END
C  /* Deck fckout */
      SUBROUTINE FCKOUT(FMAT,DMAT,NDMAT,SO,IPNTNO,IPNTRP,IPNTLG,
     &                  FIRST,LAST,THRESH,NINDAB,NINDCD,IPRINT)
C
C     Write out blocks of symmetry integrals, eliminating duplicates
C
C                                          880412   PRT
C
C     Some low-brain work has been done by TUH
C
C     Rewritten to allow for triangular looping and to eliminate
C     all duplicates by testing  880601 TUH
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
      LOGICAL DCMPAB, DCMPCD, DCMPAC, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, IAEQIC, IALTIC, IPNTLG(3,*), NOTEST
      DIMENSION SO(*), BUF(600), IBUF(4,600), IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          FMAT(*), DMAT(*)
#include <twocom.h>
#include <symmet.h>
      SAVE BUF, IBUF, ICOUNT
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine FCKOUT',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL? ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB?', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR?', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
         WRITE (LUPRI,'(2X,A,2L5)') 'DIAGAB/CD', DIAGAB, DIAGCD
         WRITE (LUPRI,'(2X,A,2L5)') 'TCONAB/CD', TCONAB, TCONCD
         WRITE (LUPRI,'(2X,A,2L5)') 'SHAEQB/CD', SHAEQB, SHCEQD
         WRITE (LUPRI,'(2X,A, L5)') 'SHABAB', SHABAB
      END IF
C
      IF (FIRST) ICOUNT = 0
      NSTART = ICOUNT
      NBUFCL = 0
      ISOFF = 0
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
                  SOINT = SO(ISOFF+INT)
                  IF (ABS(SOINT) .GT. THRESH) THEN
                     IC = KHKTC*(NINDCD(ICD,1) - 1)
                     ID = KHKTD*(NINDCD(ICD,2) - 1)
                     INDC = IPTSYM(NSTRNC + IC,IREPC)
                     INDD = IPTSYM(NSTRND + ID,IREPD)
                     ICOUNT = ICOUNT + 1
                     BUF   (ICOUNT) = SOINT
                     IBUF(1,ICOUNT) = INDA
                     IBUF(2,ICOUNT) = INDB
                     IBUF(3,ICOUNT) = INDC
                     IBUF(4,ICOUNT) = INDD
                     IF (ICOUNT.EQ.600) THEN
                        NBUFCL = NBUFCL + 1
                        CALL FCKDIR(FMAT,DMAT,NDMAT,BUF,IBUF,ICOUNT,
     &                              IPRINT)
                        ICOUNT = 0
                     END IF
                  END IF
  210          CONTINUE
  200       CONTINUE
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
                  SOINT = SO(ISOFF+INT)
                  IF (ABS(SOINT) .GT. THRESH) THEN
                     INDC = IPTSYM(NSTRNC + IC,IREPC)
                     INDD = IPTSYM(NSTRND + ID,IREPD)
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT)   = SOINT
                     IBUF(1,ICOUNT) = INDA
                     IBUF(2,ICOUNT) = INDB
                     IBUF(3,ICOUNT) = INDC
                     IBUF(4,ICOUNT) = INDD
                     IF (ICOUNT.EQ.600) THEN
                        CALL FCKDIR(FMAT,DMAT,NDMAT,BUF,IBUF,ICOUNT,
     &                              IPRINT)
                        ICOUNT = 0
                        NBUFCL = NBUFCL + 1
                     END IF
                  END IF
  310          CONTINUE
  300       CONTINUE
         END IF
         ISOFF = ISOFF + NOABCD
  100 CONTINUE
      NGINT = 600*NBUFCL + ICOUNT - NSTART
      CALL DELSTA(0,NGINT)
C
C     *************************************
C     ***** Last call to empty buffer *****
C     *************************************
C
      IF (LAST .AND. (ICOUNT.GT.0)) THEN
         CALL FCKDIR(FMAT,DMAT,NDMAT,BUF,IBUF,ICOUNT,IPRINT)
      END IF
      RETURN
      END
C  /* Deck fckdir */
      SUBROUTINE FCKDIR(FMAT,DMAT,NDMAT,BUF,IBUF,LENGTH,IPRINT)
C
C     Henrik Koch and Trygve Helgaker 18-NOV-1991.
C
C     This subroutine adds derivative two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.50D00)
      INTEGER P, Q, R, S
#include <inforb.h>
      DIMENSION FMAT(NBAST,NBAST,NDMAT), DMAT(NBAST,NBAST,NDMAT),
     &          BUF(LENGTH), IBUF(4,LENGTH)
#include <ibtfun.h>
C
      DO 100 INT = 1, LENGTH
         DINT = BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DP5*DINT
         IF (R.EQ.S)              DINT = DP5*DINT
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DP5*DINT
         EINT = DP5*DINT
         DO 200 I = 1, NDMAT
            FMAT(P,Q,I) = FMAT(P,Q,I) + DINT*(DMAT(R,S,I)+DMAT(S,R,I))
            FMAT(Q,P,I) = FMAT(Q,P,I) + DINT*(DMAT(R,S,I)+DMAT(S,R,I))
            FMAT(R,S,I) = FMAT(R,S,I) + DINT*(DMAT(P,Q,I)+DMAT(Q,P,I))
            FMAT(S,R,I) = FMAT(S,R,I) + DINT*(DMAT(P,Q,I)+DMAT(Q,P,I))
            FMAT(P,R,I) = FMAT(P,R,I) - EINT*DMAT(S,Q,I)
            FMAT(R,P,I) = FMAT(R,P,I) - EINT*DMAT(Q,S,I)
            FMAT(P,S,I) = FMAT(P,S,I) - EINT*DMAT(R,Q,I)
            FMAT(S,P,I) = FMAT(S,P,I) - EINT*DMAT(Q,R,I)
            FMAT(Q,R,I) = FMAT(Q,R,I) - EINT*DMAT(S,P,I)
            FMAT(R,Q,I) = FMAT(R,Q,I) - EINT*DMAT(P,S,I)
            FMAT(Q,S,I) = FMAT(Q,S,I) - EINT*DMAT(R,P,I)
            FMAT(S,Q,I) = FMAT(S,Q,I) - EINT*DMAT(P,R,I)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck dr2out */
      SUBROUTINE DR2OUT(SO,IATOM,MULE,FIRST,LAST,SPNORB,SQ12EL,THRESH,
     &                  IPNTNO,IPNTRP,IPNTLG,NINDAB,NINDCD,IPRINT)
C
C     tuh Sep 92
C
C     Write derivative integrals (IATOM .gt. 0), magnetic integrals
C     (IATOM .eq. 0), or spin-orbit integrals (SPNORB .eqv. .TRUE.)
C     to file (LU2DER).
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
      LOGICAL NOTEST, DCMPAB, DCMPCD, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, SPNORB, SQ12EL, IPNTLG(3,NINTMX,*)
      DIMENSION IPNTNO(4,NINTMX,*), IPNTRP(3,NINTMX,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
      DIMENSION SO(*), BUF(600), IBUF(600)
#include <nuclei.h>
#include <twocom.h>
#include <dorps.h>
#include <symmet.h>
#include <doxyz.h>
      SAVE BUF, IBUF, ICOUNT, NBUFCL
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine DR2OUT',-1)
C
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
         CALL QUIT('More than 255 basis functions is not'//
     &   ' implemented for 32-bit integers, use 64-bit integers')
!        hjaaj Feb 2011: the problem is that NIBUF=2 is not implemented
!        here, see DS2OU1 for how to implement if you have the time ...
      END IF
#endif
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
      IF (FIRST) CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,-1,
     &                       SPNORB,IPRINT)
C
      ISOFF = 0
      DO 100 ICOOR = 1, 3
      IF (DOXYZ(ICOOR)) THEN
         DO 200 IREPE = 0,MAXREP
            IF (SPNORB) THEN
             IRPXYZ=IBTXOR(ISYMAX(1,1),IBTXOR(ISYMAX(2,1),ISYMAX(3,1)))
             IF (IREPE .NE. IBTXOR(ISYMAX(ICOOR,1),IRPXYZ)) GOTO 200
            ELSE
             IF(.NOT.DOREPS(IREPE)) GO TO 200
             IF(IBTAND(MULE,IBTXOR(IREPE,ISYMAX(ICOOR,1))).NE.0)GOTO 200
            END IF
            CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,2,SPNORB,
     &                  IPRINT)
            NBUFCL = 0
            IREPX  = IPTREP(IREPE,2)
            NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
            DO 300 I = 1, NINTSR(IREPX)
               NSTRNA = IPNTNO(1,I,IREPX)
               NSTRNB = IPNTNO(2,I,IREPX)
               NSTRNC = IPNTNO(3,I,IREPX)
               NSTRND = IPNTNO(4,I,IREPX)
               IREPA  = IPNTRP(1,I,IREPX)
               IREPB  = IPNTRP(2,I,IREPX)
               IREPC  = IPNTRP(3,I,IREPX)
               IREPD  = IBTXOR(IBTXOR(IBTXOR(IREPA,IREPB),IREPC),IREPE)
               IF (NOTEST) THEN
                  INT = 0
                  DO 400 IAB = 1, NORBAB
                     IA = KHKTA*(NINDAB(IAB,1) - 1)
                     IB = KHKTB*(NINDAB(IAB,2) - 1)
                     INDA = IPTSYM(NSTRNA + IA,IREPA)
                     INDB = IPTSYM(NSTRNB + IB,IREPB)
                     INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
                     DO 410 ICD = 1, NORBCD
                        INT = INT + 1
                        SOINT = SO(ISOFF+INT)
                        IF (ABS(SOINT) .GT. THRESH) THEN
                           IC = KHKTC*(NINDCD(ICD,1) - 1)
                           ID = KHKTD*(NINDCD(ICD,2) - 1)
                           INDC = IPTSYM(NSTRNC + IC,IREPC)
                           INDD = IPTSYM(NSTRND + ID,IREPD)
                           IF (SQ12EL) THEN
                              IF (INDC .GT. INDD) THEN
                                 INDCD  = INDC*(IBIT1+1) + INDD
                              ELSE IF (SPNORB) THEN
                                 SOINT = - SOINT
                                 INDCD  = INDD*(IBIT1+1) + INDC
                              ELSE
                                 INDCD  = INDD*(IBIT1+1) + INDC
                              END IF
                              IABCD = INDAB*(IBIT2+1) + INDCD
                           ELSE
                              INDCD = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                              IABCD = MAX(INDAB,INDCD)*IBIT2+INDAB+INDCD
                           END IF
                           ICOUNT = ICOUNT + 1
                           BUF (ICOUNT) = SOINT
                           IBUF(ICOUNT) = IABCD
                           IF (ICOUNT.EQ.600) THEN
                              NBUFCL = NBUFCL + 1
                              CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                    IREPE,0,SPNORB,IPRINT)
                           END IF
                        END IF
  410                CONTINUE
  400             CONTINUE
               ELSE
                  DCMPAB = IPNTLG(1,I,IREPX)
                  DCMPCD = IPNTLG(2,I,IREPX)
                  DCABAB = IPNTLG(3,I,IREPX)
                  DRALTB = IREPA .LT. IREPB
                  DRCLTD = IREPC .LT. IREPD
                  DRABAB = DCABAB.AND.IREPA.EQ.IREPC.AND.IREPB.EQ.IREPD
                  INT = 0
                  DO 500 IAB = 1, NORBAB
                     IA = KHKTA*(NINDAB(IAB,1) - 1)
                     IB = KHKTB*(NINDAB(IAB,2) - 1)
                     IF (DCMPAB) THEN
                        IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                           INT = INT + NORBCD
                           GO TO 500
                        END IF
                     END IF
                     INDA = IPTSYM(NSTRNA + IA,IREPA)
                     INDB = IPTSYM(NSTRNB + IB,IREPB)
                     INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
                     DO 510 ICD = 1,NORBCD
                        IC = KHKTC*(NINDCD(ICD,1) - 1)
                        ID = KHKTD*(NINDCD(ICD,2) - 1)
                        INT = INT + 1
                        IF (DCMPCD ) THEN
                           IF (ID.GT.IC) GO TO 510
                           IF (DRCLTD .AND. ID.EQ.IC) GO TO 510
                        END IF
                        IF (DRABAB) THEN
                           IF (IA.LT.IC.OR.(IA.EQ.IC.AND.IB.LT.ID))
     &                        GOTO 510
                        END IF
                        SOINT = SO(ISOFF+INT)
                        IF (ABS(SOINT) .GT. THRESH) THEN
                           INDC = IPTSYM(NSTRNC + IC,IREPC)
                           INDD = IPTSYM(NSTRND + ID,IREPD)
                           IF (SQ12EL) THEN
                              IF (INDC .GT. INDD) THEN
                                 INDCD  = INDC*(IBIT1+1) + INDD
                              ELSE IF (SPNORB) THEN
                                 SOINT = - SOINT
                                 INDCD  = INDD*(IBIT1+1) + INDC
                              ELSE
                                 INDCD  = INDD*(IBIT1+1) + INDC
                              END IF
                              IABCD = INDAB*(IBIT2+1) + INDCD
                           ELSE
                              INDCD = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                              IABCD = MAX(INDAB,INDCD)*IBIT2+INDAB+INDCD
                           END IF
                           ICOUNT = ICOUNT + 1
                           BUF (ICOUNT) = SOINT
                           IBUF(ICOUNT) = IABCD
                           IF (ICOUNT.EQ.600) THEN
                              NBUFCL = NBUFCL + 1
                              CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                    IREPE,0,SPNORB,IPRINT)
                           END IF
                        END IF
  510                CONTINUE
  500             CONTINUE
               END IF
               ISOFF = ISOFF + NOABCD
  300       CONTINUE
  200    CONTINUE
      END IF
  100 CONTINUE
      IF (LAST) CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,1,SPNORB,
     &                      IPRINT)
      RETURN
      END
C  /* Deck dr2wrt */
      SUBROUTINE DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,ITYPE,SPNORB,
     &                  IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
      PARAMETER (TEN14 = 1.1D14, LBUF = 600)
      LOGICAL SPNORB
      DIMENSION BUF(600), IBUF(600)
#include <inftap.h>
#include <nodint.h>
#include <symmet.h>
#include <chrxyz.h>
#if defined (VAR_SPLITFILES)
#include <2gbdef.h>
      CHARACTER*80 FNNAME, FNNM2
#include <chrnos.h>
#endif
      SAVE NBUF, NFLAGS, IRPOLD, ICROLD, ICTOLD, NBFOLD, JBUF
C
#include <ibtfun.h>
C
C     Initialization
C     ==============
C
      IF (ITYPE .EQ. -1) THEN
         IF (SPNORB) THEN
            CALL GPOPEN(LU2DER,'AO2SOINT','UNKNOWN',' ','UNFORMATTED',
     *                  IDUMMY,.FALSE.)
            CALL REWSPL(LU2DER)
            CALL NEWLAB('AO2SOINT',LU2DER,LUPRI)
         ELSE IF (IATOM .EQ. 0) THEN
            CALL GPOPEN (LU2DER,'AO2MGINT','UNKNOWN',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)   
            CALL REWSPL(LU2DER)
            CALL NEWLAB('AO2MGINT',LU2DER,LUPRI)
         ELSE
            IF (LU2DER .LE. 0) THEN
               CALL QUIT('illegal unit number LU2DER in DR2WRT')
            ENDIF
            CALL REWSPL(LU2DER)
         END IF
C
         JBUF   = 256
         ICOUNT =  0
         ICTOLD =  0
         NBUF   =  0
         NBFOLD =  0
         NFLAGS =  0
         IRPOLD = -1
         ICROLD = -1
         DO 100 IREP = 0, MAXREP
            NDSINT(1,IREP) = 0
            NDSINT(2,IREP) = 0
            NDSINT(3,IREP) = 0
  100    CONTINUE
C
C     Write new buffer
C     ================
C
      ELSE IF (ITYPE .EQ. 0) THEN
C
#if defined (VAR_SPLITFILES)
         IF ((JBUF + LBUF*(IRAT + 1) + 1) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
            INQUIRE(UNIT=LU2DER,NAME=FNNAME)
            LN = 1
 12         CONTINUE
            IF (FNNAME(LN:LN) .NE. ' ') THEN
               LN = LN + 1
               GOTO 12
            END IF
            LN = LN - 1
            CALL GPCLOSE(LU2DER,'KEEP')
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
            CALL GPOPEN(LU2DER,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &           .FALSE.)
            JBUF = 0
         END IF
         WRITE (LU2DER) BUF,IBUF,ICOUNT
         JBUF = JBUF + LBUF*(IRAT + 1) + 1
#else
         WRITE (LU2DER) BUF, IBUF, ICOUNT
#endif
         IF (IPRINT .GE. 6) CALL DR2PRI(BUF,IBUF,ICOUNT,NBUF)
         NBUF   = NBUF + 1
         ICOUNT = 0
C
C     Insert tag for new coordinate or symmetry
C     =========================================
C
      ELSE IF (ITYPE .EQ. 2) THEN
         IF (IREPE.NE.IRPOLD .OR. ICOOR.NE.ICROLD) THEN
            INTS = (NBUF - NBFOLD)*600 + ICOUNT - ICTOLD
            IF (INTS .GT. 0) NDSINT(ICROLD,IRPOLD) =
     &                       NDSINT(ICROLD,IRPOLD) + INTS
            ICOUNT = ICOUNT + 1
            BUF(ICOUNT)  = TEN14
#ifdef VAR_INT64
            IBUF(ICOUNT) = IATOM*2**48+IREPE*2**32+ICOOR*2**16
#else
            IBUF(ICOUNT) = IATOM*2**24+IREPE*2**16+ICOOR*2**8
#endif
            ICTOLD = ICOUNT
            NBFOLD = NBUF
            IRPOLD = IREPE
            ICROLD = ICOOR
            NFLAGS = NFLAGS + 1
            IF (ICOUNT .EQ. 600) THEN
#if defined (VAR_SPLITFILES)
C
               IF ((JBUF + LBUF*(IRAT + 1) + 1) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
                  INQUIRE(UNIT=LU2DER,NAME=FNNAME)
                  LN = 1
 11               CONTINUE
                  IF (FNNAME(LN:LN) .NE. ' ') THEN
                     LN = LN + 1
                     GOTO 11
                  END IF
                  LN = LN - 1
                  CALL GPCLOSE(LU2DER,'KEEP')
                  I = LN - 1
                  IF (FNNAME(I:I) .NE. '-') THEN
                     FNNM2 = FNNAME(1:LN)//'-0'
                     LN = LN + 2
                  ELSE
                     READ(FNNAME(LN:),'(I1)') INUM
                     INUM = INUM + 1
                     IF (INUM .GT. 9) THEN
                        WRITE (LUPRI,'(/A)') ' DALTON needs to split '//
     &                 ' a file more than 11 times.',
     &                 ' This is currently not supported'
                        CALL QUIT('Too many splittings of a file')
                     END IF
                     FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
                  END IF
                  CALL GPOPEN(LU2DER,FNNM2(1:LN),'UNKNOWN',' ',' ',
     &                        IDUMMY,.FALSE.)
                  JBUF = 0
               END IF
               WRITE (LU2DER) BUF,IBUF,ICOUNT
               JBUF = JBUF + LBUF*(IRAT + 1) + 1
#else 
               WRITE (LU2DER) BUF,IBUF,ICOUNT
#endif
               IF (IPRINT.GE.6) CALL DR2PRI(BUF,IBUF,ICOUNT,NBUF)
               NBUF   = NBUF + 1
               ICOUNT = 0
            END IF
         END IF
C
C     Empty final buffer and do statistics
C     ====================================
C
      ELSE
         INTS = (NBUF - NBFOLD)*600 + ICOUNT - ICTOLD
         IF (INTS .GT. 0) NDSINT(ICROLD,IRPOLD) =
     &                    NDSINT(ICROLD,IRPOLD) + INTS
#if defined (VAR_SPLITFILES)
C
         IF ((JBUF + 2*LBUF*(IRAT + 1) + 1) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
            INQUIRE(UNIT=LU2DER,NAME=FNNAME)
            LN = 1
 10         CONTINUE
            IF (FNNAME(LN:LN) .NE. ' ') THEN
               LN = LN + 1
               GOTO 10
            END IF
            LN = LN - 1
            CALL GPCLOSE(LU2DER,'KEEP')
            I = LN - 1
            IF (FNNAME(I:I) .NE. '-') THEN
               FNNM2 = FNNAME(1:LN)//'-0'
               LN = LN + 2
            ELSE
               READ(FNNAME(LN:),'(I1)') INUM
               INUM = INUM + 1
               IF (INUM .GT. 9) THEN
                  WRITE (LUPRI,'(/A)') ' DALTON needs to split '//
     &                 ' a file more than 11 times.',
     &                 ' This is currently not supported'
                  CALL QUIT('Too many splittings of a file')
               END IF
               FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
            END IF
            CALL GPOPEN(LU2DER,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &           .FALSE.)
            JBUF = 0
         END IF
         IF (ICOUNT .GT. 0) THEN
            WRITE (LU2DER) BUF, IBUF, ICOUNT
            IF (IPRINT .GE. 6) CALL DR2PRI(BUF,IBUF,ICOUNT,NBUF)
         END IF
         WRITE (LU2DER) BUF, IBUF, -1
         JBUF = JBUF + 2*LBUF*(IRAT + 1) + 1
#else
         IF (ICOUNT .GT. 0) THEN
            WRITE (LU2DER) BUF, IBUF, ICOUNT
            IF (IPRINT .GE. 6) CALL DR2PRI(BUF,IBUF,ICOUNT,NBUF)
         END IF
         WRITE (LU2DER) BUF, IBUF, -1
#endif
         END FILE LU2DER
         IF (SPNORB) CALL GPCLOSE(LU2DER,'KEEP')
         NTOTAL = 600*NBUF + ICOUNT - NFLAGS
         NOINTS = NTOTAL
       IF (IPRINT .GE. 2) THEN
         IF (SPNORB) THEN
            CALL AROUND('Number of written 2-el. spin-orbit integrals')
         ELSE
            CALL AROUND('Number of written derivative integrals')
         END IF
         CALL HEADER('Symmetry           x         y         z',1)
         DO 500 IREP = 0, MAXREP
            IF (NDSINT(1,IREP)+NDSINT(2,IREP)+NDSINT(3,IREP).GT.0) THEN
               WRITE (LUPRI,'(2X,I5,5X,3I10)')
     &               IREP + 1, (NDSINT(I,IREP),I=1,3)
            END IF
 500     CONTINUE
         IF (SPNORB) THEN
            WRITE (LUPRI,'(/1X,I10,A,/,1X,I10,A//)')
     &         NTOTAL,' spin-orbit two-electron integrals and ',
     &         NFLAGS,' flags have been written on disk.'
         ELSE
            WRITE (LUPRI,'(/1X,I10,A,/,1X,I10,A//)')
     &         NTOTAL,' differentiated two-electron integrals and ',
     &         NFLAGS,' flags have been written on disk.'
         END IF
       END IF
       RETURN
      END IF
C
      RETURN
      END
C  /* Deck dr2pri */
      SUBROUTINE DR2PRI(BUF,IBUF,ICOUNT,NBUF)
#include <implicit.h>
#include <priunit.h>
#include <ibtpar.h>
      DIMENSION BUF(600), IBUF(600)
      SAVE IATOM, ICOOR, IREPE
#include <chrxyz.h>
#include <ibtfun.h>
      WRITE(LUPRI,'(2X,A,I5,A/)')
     &   'Integral buffer #',NBUF + 1,' has been written.'
      DO 100 INT = 1, ICOUNT
         I = IBTAND(IBTSHR(IBUF(INT),24),IBT08)
         J = IBTAND(IBTSHR(IBUF(INT),16),IBT08)
         K = IBTAND(IBTSHR(IBUF(INT),8), IBT08)
         L = IBTAND(       IBUF(INT),    IBT08)
         IF (L .EQ. 0) THEN
            ICOOR = IBTAND(IBTSHR(IBUF(INT), 8),IBT08)
            IREPE = IBTAND(IBTSHR(IBUF(INT),16),IBT08)
            IATOM = IBTAND(IBTSHR(IBUF(INT),24),IBT08)
         ELSE
            WRITE (LUPRI,'(10X,A,2X,4I4,5X,I2,A,I2,5X,1P,D16.8)')
     &         ' ## ',I,J,K,L,IATOM,CHRXYZ(ICOOR),IREPE,BUF(INT)
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck mg2out */
      SUBROUTINE MG2OUT(SO,FIRST,LAST,THRESH,IPNTNO,IPNTRP,IPNTLG,
     &                  NINDAB,NINDCD,IPRINT)
C
C     tuh Sep 92
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
      LOGICAL NOTEST, DCMPAB, DCMPCD, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, SPNORB, IPNTLG(3,NINTMX,*)
      DIMENSION IPNTNO(4,NINTMX,*), IPNTRP(3,NINTMX,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
      DIMENSION SO(*), BUF(600), IBUF(600)
#include <nuclei.h>
#include <twocom.h>
#include <dorps.h>
#include <symmet.h>
#include <doxyz.h>
      SAVE BUF, IBUF, ICOUNT, NBUFCL
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine MG2OUT',-1)

#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
         CALL QUIT('More than 255 basis functions and mag.der. not'//
     &   ' implemented for 32-bit integers, use 64-bit integers')
!        hjaaj Feb 2011: the problem is that NIBUF=2 is not implemented
!        here, see DS2OU1 for how to implement if you have the time ...
      END IF
#endif
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
      IATOM = 0
      SPNORB = .FALSE.
      IF (FIRST) CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,-1,
     &                       SPNORB,IPRINT)
C
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      NBUFCL = 0
      ISOFF0  = 0
      DO 100 ICOOR = 1, 3
         IREPE = ISYMAX(ICOOR,2)
         CALL DR2WRT(BUF,IBUF,ICOUNT,0,ICOOR,IREPE,2,.FALSE.,IPRINT)
         IREPX = IPTREP(IREPE,2)
         NINTX = NINTSR(IREPX)
         ISOFF1 = ISOFF0
         ISOFF2 = ISOFF0 + NOABCD*NINTX
         DO 200 I = 1, NINTX
            NSTRNA = IPNTNO(1,I,IREPX)
            NSTRNB = IPNTNO(2,I,IREPX)
            NSTRNC = IPNTNO(3,I,IREPX)
            NSTRND = IPNTNO(4,I,IREPX)
            IREPA  = IPNTRP(1,I,IREPX)
            IREPB  = IPNTRP(2,I,IREPX)
            IREPC  = IPNTRP(3,I,IREPX)
            IREPD  = IBTXOR(IBTXOR(IBTXOR(IREPA,IREPB),IREPC),IREPE)
            IF (NOTEST) THEN
               INT = 0
               DO 400 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  DO 410 ICD = 1, NORBCD
                     IC = KHKTC*(NINDCD(ICD,1) - 1)
                     ID = KHKTD*(NINDCD(ICD,2) - 1)
                     INDC = IPTSYM(NSTRNC + IC,IREPC)
                     INDD = IPTSYM(NSTRND + ID,IREPD)
                     INT = INT + 1
                     SOINT = SO(ISOFF1 + INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                        INDAB = INDA*(IBIT1+1) + INDB
                        INDCD = INDC*(IBIT1+1) + INDD
                        ICOUNT = ICOUNT + 1
                        BUF (ICOUNT) = SOINT
                        IBUF(ICOUNT) = INDAB*(IBIT2+1) + INDCD
                        IF (ICOUNT.EQ.600) THEN
                           NBUFCL = NBUFCL + 1
                           CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                 IREPE,0,SPNORB,IPRINT)
                        END IF
                     END IF
                     SOINT = SO(ISOFF2 + INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                        INDAB = INDA*(IBIT1+1) + INDB
                        INDCD = INDD*(IBIT1+1) + INDC
                        ICOUNT = ICOUNT + 1
                        BUF (ICOUNT) = SOINT
                        IBUF(ICOUNT) = INDAB*(IBIT2+1) + INDCD
                        IF (ICOUNT.EQ.600) THEN
                           NBUFCL = NBUFCL + 1
                           CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                 IREPE,0,SPNORB,IPRINT)
                        END IF
                     END IF
  410             CONTINUE
  400          CONTINUE
            ELSE
               DCMPAB = IPNTLG(1,I,IREPX)
               DCMPCD = IPNTLG(2,I,IREPX)
               DCABAB = IPNTLG(3,I,IREPX)
               DRALTB = IREPA .LT. IREPB
               DRCLTD = IREPC .LT. IREPD
               DRABAB = DCABAB.AND.IREPA.EQ.IREPC.AND.IREPB.EQ.IREPD
               INT = 0
               DO 500 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  IF (DCMPAB) THEN
                     IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                        INT = INT + NORBCD
                        GO TO 500
                     END IF
                  END IF
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = MAX(INDA,INDB)*(IBIT1+1) + INDA + INDB
                  DO 510 ICD = 1,NORBCD
                     IC   = KHKTC*(NINDCD(ICD,1) - 1)
                     ID   = KHKTD*(NINDCD(ICD,2) - 1)
                     INDC = IPTSYM(NSTRNC + IC,IREPC)
                     INDD = IPTSYM(NSTRND + ID,IREPD)
                     INT = INT + 1
                     IF (DCMPCD ) THEN
                        IF (ID.GT.IC) GO TO 510
                        IF (DRCLTD .AND. ID.EQ.IC) GO TO 510
                     END IF
                     IF (DRABAB) THEN
                        IF (IA.LT.IC.OR.(IA.EQ.IC.AND.IB.LT.ID))
     &                     GOTO 510
                     END IF
                     SOINT = SO(ISOFF1 + INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                        INDAB = INDA*(IBIT1+1) + INDB
                        INDCD = INDC*(IBIT1+1) + INDD
                        ICOUNT = ICOUNT + 1
                        BUF (ICOUNT) = SOINT
                        IBUF(ICOUNT) = INDAB*(IBIT2+1) + INDCD
                        IF (ICOUNT.EQ.600) THEN
                           NBUFCL = NBUFCL + 1
                           CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                 IREPE,0,SPNORB,IPRINT)
                        END IF
                     END IF
                     SOINT = SO(ISOFF2 + INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                     IF (INDA .EQ. INDB) GO TO 510
                     IF (INDC .EQ. INDD) GO TO 510
                        INDAB = INDA*(IBIT1+1) + INDB
                        INDDC = INDD*(IBIT1+1) + INDC
                        ICOUNT = ICOUNT + 1
                        BUF (ICOUNT) = SOINT
                        IBUF(ICOUNT) = INDAB*(IBIT2+1) + INDDC
                        IF (ICOUNT.EQ.600) THEN
                           NBUFCL = NBUFCL + 1
                           CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,
     &                                 IREPE,0,SPNORB,IPRINT)
                        END IF
                     END IF
  510             CONTINUE
  500          CONTINUE
            END IF
            ISOFF1 = ISOFF1 + NOABCD
            ISOFF2 = ISOFF2 + NOABCD
  200    CONTINUE
         ISOFF0 = ISOFF0 + 2*NINTX*NOABCD
  100 CONTINUE
      IF (LAST) CALL DR2WRT(BUF,IBUF,ICOUNT,IATOM,ICOOR,IREPE,1,SPNORB,
     &                      IPRINT)
      RETURN
      END
C  /* Deck ds2out */
      SUBROUTINE DS2OUT(SO,WRKBUF,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,
     &                  THRESH,NINDAB,NINDCD,IORBSH,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      LOGICAL FIRST, LAST, IPNTLG(*)
      DIMENSION SO(*), WRKBUF(*), IPNTNO(*), IPNTRP(*), NINDAB(*),
     &          NINDCD(*), IORBSH(*)
#include <disbuf.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
C
C---------------------------------
C     Call sort and write routine.
C---------------------------------
C
Cholesky
      IF (CHOINT) THRESH = 1.0D-40
Cholesky
C
      CALL DS2OU1(SO,WRKBUF(KDSBF),WRKBUF(KDSIBF),WRKBUF(KDSNCT),
     &            WRKBUF(KDSORB),WRKBUF(KORBDS),IPNTNO,IPNTRP,IPNTLG,
     &            FIRST,LAST,THRESH,NINDAB,NINDCD,LDSBUF,NDIST,IORBSH,
     &            IPRINT)
C
      RETURN
      END
C  /* Deck ds2ou1 */
      SUBROUTINE DS2OU1(SO,BUF,IBUF,NCOUNT,IDSORB,IORBDS,IPNTNO,IPNTRP,
     &                  IPNTLG,FIRST,LAST,THRESH,NINDAB,NINDCD,LDSBUF,
     &                  NDIST,IORBSH,IPRINT)
C
C     Write out blocks of symmetry integrals, eliminating duplicates
C     Used when "ADISTR" distributions (all gabcd for fixed a)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
#include <eribuf.h>
#include <nuclei.h>
C
      INTEGER*8 NWRIT
      SAVE      NWRIT
      LOGICAL DCMPCD, DRCLTD, FIRST, LAST, IPNTLG(3,*), NOTEST,
     &        DOINDX
      DIMENSION SO(*), BUF(LBUF,NDIST), IBUF(LBUF,2,NDIST),
     &          IPNTNO(4,*), IPNTRP(3,*), NCOUNT(NDIST),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          IDSORB(NDIST), IORBDS(NBASIS), IORBSH(*)
#include <twocom.h>
#include <symmet.h>
#include <drw2el.h>
C
Cholesky
#include <ccdeco.h>
Cholesky
C
#include <ibtfun.h>
C
      IF (LBUF .NE. LDSBUF) THEN
        WRITE (LUPRI,*) 'LBUF  :',LBUF
        WRITE (LUPRI,*) 'LDSBUF:',LDSBUF
        CALL QUIT('Error in DS2OU1')
      END IF
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine DS2OUT',-1)
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
         WRITE (LUPRI,'(2X,A, I5)') 'NDIST ', NDIST
      END IF
C
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) THEN
         CALL IZERO(NCOUNT,NDIST)
         CALL UN2WRN(BUF,IBUF,ICOUNT,-1,0,IPRINT)
C         CALL UN2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,-1,NBITS,0,IPRINT)
         DOINDX = .TRUE.
         CALL AINDEX(ISHELA,NAINTS,IDSORB,DOINDX,IORBSH,IPRINT)
         DO 50 IDIST = 1, NDIST
             IORBDS(IDSORB(IDIST)) = IDIST
   50    CONTINUE
         NWRIT = 0
      END IF
C
      ISOFF  = 0
      NBUFCL = 0
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
            IF (NIBUF .EQ. 1) THEN
               INT = 0
               DO 200 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
Cholesky
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
Cholesky
                  INDAB = INDA*(IBIT1 + 1) + IPTSYM(NSTRNB + IB,IREPB)
                  IDIST = IORBDS(INDA)
Cholesky
                  IF (DIACAL .AND. CHOINT) THEN
                     IDIST = 1
                  ENDIF
Cholesky
                  DO 210 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
Cholesky
                       IF (DIACAL .AND. CHOINT) THEN
                          IF (INDA .NE.  INDC) GOTO 210
                          IF (INDB .NE.  INDD) GOTO 210
                       ENDIF
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
Cholesky
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       IF (BPH2OO.AND.INDD.GT.INDC) THEN 
                          BUF (ICOUNT,IDIST) = -SOINT
                       ELSE
                          BUF (ICOUNT,IDIST) = SOINT
                       END IF
                       IBUF(ICOUNT,1,IDIST) = INDAB*(IBIT2 + 1) + INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
Cholesky
                          INDX = INDA
                          IF (DIACAL .AND. CHOINT) THEN
                             INDX = 1
                          ENDIF
Cholesky
                          CALL UN2WRN(BUF(1,IDIST),IBUF(1,1,IDIST),
     &                                ICOUNT,0,INDX,IPRINT)
Cho  &                                ICOUNT,0,INDA,IPRINT)
C                          CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,
C     &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
  210             CONTINUE
  200          CONTINUE
            ELSE
               INT = 0
               DO 205 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = INDA*(IBIT1 + 1) + INDB
                  IDIST = IORBDS(INDA)
Cholesky
                  IF (DIACAL .AND. CHOINT) THEN
C                    IF (INDB .GT.  INDA) GOTO 205
                     IDIST = 1
                  ENDIF
Cholesky
                  DO 215 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
Cholesky
                       IF (DIACAL .AND. CHOINT) THEN
                          IF (INDA .NE.  INDC) GOTO 215
                          IF (INDB .NE.  INDD) GOTO 215
                       ENDIF
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
Cholesky
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       IF (BPH2OO.AND.INDD.GT.INDC) THEN 
                          BUF (ICOUNT,IDIST) = -SOINT
                       ELSE
                          BUF (ICOUNT,IDIST) = SOINT
                       END IF
                       IBUF(ICOUNT,1,IDIST) = INDAB
                       IBUF(ICOUNT,2,IDIST) = INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
Cholesky
                          INDX = INDA
                          IF (DIACAL .AND. CHOINT) THEN
                             INDX = 1
                          ENDIF
Cholesky
                          CALL UN2WRN(BUF(1,IDIST),IBUF(1,1,IDIST),
     &                                ICOUNT,0,INDX,IPRINT)
Cho  &                                ICOUNT,0,INDA,IPRINT)
C                          CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,
C     &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
  215             CONTINUE
  205          CONTINUE
            END IF
         ELSE
            DCMPCD = IPNTLG(2,I)
            DRCLTD = IREPC .LT. IREPD
            INT = 0
            DO 300 IAB = 1, NORBAB
               IA = KHKTA*(NINDAB(IAB,1) - 1)
               IB = KHKTB*(NINDAB(IAB,2) - 1)
               INDA = IPTSYM(NSTRNA + IA,IREPA)
               INDB = IPTSYM(NSTRNB + IB,IREPB)
               INDAB = INDA*(IBIT1 + 1) + INDB
               IDIST = IORBDS(INDA)
Cholesky
               IF (DIACAL .AND. CHOINT) THEN
C                 IF (INDB .GT. INDA) GOTO 300
                  IDIST = 1
               ENDIF
Cholesky
               DO 310 ICD = 1,NORBCD
                  IC = KHKTC*(NINDCD(ICD,1) - 1)
                  ID = KHKTD*(NINDCD(ICD,2) - 1)
                  INT = INT + 1
                  IF (DCMPCD ) THEN
                     IF (ID.GT.IC) GO TO 310
                     IF (DRCLTD .AND. ID.EQ.IC) GO TO 310
                  END IF
                  SOINT = SO(ISOFF+INT)
                  IF (ABS(SOINT) .GT. THRESH) THEN
                     IF (NIBUF .EQ. 1) THEN
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
Cholesky
                       IF (DIACAL .AND. CHOINT) THEN
                          IF (INDA .NE. INDC) GOTO 310
                          IF (INDB .NE. INDD) GOTO 310
                       ENDIF
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       IF (BPH2OO.AND.INDD.GT.INDC) THEN 
                          BUF (ICOUNT,IDIST) = -SOINT
                       ELSE
                          BUF (ICOUNT,IDIST) = SOINT
                       END IF
                       IBUF(ICOUNT,1,IDIST) = INDAB*(IBIT2 + 1) + INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
Cholesky
                       INDX = INDA
                       IF (DIACAL .AND. CHOINT) THEN
                          INDX = 1
                       ENDIF
Cholesky
                          CALL UN2WRN(BUF(1,IDIST),IBUF(1,1,IDIST),
     &                                ICOUNT,0,INDX,IPRINT)
Cho  &                                ICOUNT,0,INDA,IPRINT)
C                          CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,
C     &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NBUFCL = NBUFCL + 1
                          NCOUNT(IDIST) = 0
                       END IF
                     ELSE
Cholesky               NCOUNT(IDIST) = NCOUNT(IDIST) + 1
Cholesky               ICOUNT = NCOUNT(IDIST)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
Cholesky
                       IF (DIACAL .AND. CHOINT) THEN
                          IF (INDA .NE. INDC) GOTO 310
                          IF (INDB .NE. INDD) GOTO 310
                       ENDIF
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       IF (BPH2OO.AND.INDD.GT.INDC) THEN 
                          BUF (ICOUNT,IDIST) = -SOINT
                       ELSE
                          BUF (ICOUNT,IDIST) = SOINT
                       END IF
                       IBUF(ICOUNT,1,IDIST) = INDAB
                       IBUF(ICOUNT,2,IDIST) = INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
Cholesky
                       INDX = INDA
                       IF (DIACAL .AND. CHOINT) THEN
                          INDX = 1
                       ENDIF
Cholesky
                          CALL UN2WRN(BUF(1,IDIST),IBUF(1,1,IDIST),
     &                                ICOUNT,0,INDX,IPRINT)
Cho  &                                ICOUNT,0,INDA,IPRINT)
C                          CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,
C     &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NBUFCL = NBUFCL + 1
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
                  END IF
  310          CONTINUE
  300       CONTINUE
         END IF
         ISOFF = ISOFF + NOABCD
  100 CONTINUE
      NWRIT = NWRIT + LBUF*NBUFCL
C
C     *************************************
C     ***** Last call to empty buffer *****
C     *************************************
C
      IF (LAST) THEN
Cholesky
         NUMD = NDIST
         IF (DIACAL .AND. CHOINT) THEN
            NUMD = 1
         ENDIF
Cholesky
         DO 400 IDIST = 1, NDIST
            NWRIT = NWRIT + NCOUNT(IDIST)
Cholesky
            INDX = IDSORB(IDIST)
            IF (DIACAL .AND. CHOINT) THEN
               INDX = 1
            ENDIF
            CALL UN2WRN(BUF(1,IDIST),IBUF(1,1,IDIST),
     &                  NCOUNT(IDIST),1,INDX,IPRINT)
Cholesky
Chol        CALL UN2WRN(BUF(1,IDIST),IBUF(1,1,IDIST),
Chol &                  NCOUNT(IDIST),1,IDSORB(IDIST),IPRINT)
C           CALL UN2WRT(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,NIBUF,
C     &                 NCOUNT(IDIST),1,NBITS,IDSORB(IDIST),IPRINT)
  400    CONTINUE
C        HJAaJ Sep08: calculation of PERCNT rewritten, to avoid integer*4 overflows
C        NALL   = NDIST*(NBASIS + 1)*(NBASIS**2)/2
         FNBASIS= NBASIS
         FNALL  = NDIST*(FNBASIS+1.0D0)*FNBASIS*FNBASIS*0.5D0
         PERCNT = NWRIT
         PERCNT = 100.0D0*PERCNT / FNALL
         IF ((.NOT. CHOINT) .AND. (IPRINT.GE.0)) 
     &         WRITE (LUPRI,'(/A,I10,A,F4.1,A)')
     &         ' Number of two-electron spin-orbit integrals written:',
     &         NWRIT,' (',PERCNT,'%)'
      END IF
C
      RETURN
      END
C  /* Deck un2wrn */
      SUBROUTINE UN2WRN(BUF,IBUF,ICOUNT,ITYPE,INDA,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <ibtpar.h>
#include <eribuf.h>
      DIMENSION BUF(LBUF), IBUF(LBUF,NIBUF)
#include <drw2el.h>
#include <twosta.h>
#include <inftap.h>
#include <nuclei.h>
#include <symmet.h>
#include <ibtfun.h>
C
      NBUFX(0) = NBUFX(0) + 1
C
      IF (ITYPE .EQ. -1) THEN
         REWIND LUINTR
         NBUFX(0) = NBUFX(0) - 1
c        CALL NEWLAB('BASINFO ',LUINTA,LUPRI)
c        WRITE (LUINTA) MAXREP+1,(NAOS(I),I=1,8),LBUF,NIBUF,NBITS
c        CALL NEWLAB('BASTWOEL',LUINTA,LUPRI)
c        ICOUNT = 0
c        NBUFX(0) = 0
      ELSE IF (ITYPE .EQ. 0) THEN
         IF (INDA .NE. 0) WRITE (LUINTR) INDA
         WRITE (LUAORC(0),REC=NBUFX(0)) BUF,IBUF,ICOUNT
C
         IF (IPRINT .GE. 6) THEN
c        if (bph2oo) then
            WRITE (LUPRI,'(2X,A,I5,A/)')
     &         'Integral buffer #',NBUFX(0),' has been written.'
            IBIT1 = 2**NBITS - 1
            DO 100 INT = 1, ICOUNT
               IF (NIBUF .EQ. 1) THEN
                  I = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
                  J = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
                  K = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
                  L = IBTAND(       IBUF(INT,1),         IBIT1)
               ELSE
                  I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
                  J = IBTAND(       IBUF(INT,1),       IBIT1)
                  K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
                  L = IBTAND(       IBUF(INT,2),       IBIT1)
               END IF
               WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8)')
     &                      ' ## ', I, J, K, L, BUF(INT)
  100       CONTINUE
         END IF
         ICOUNT = 0
      ELSE
         IF (INDA .NE. 0) THEN
            WRITE (LUINTR) INDA
            WRITE (LUAORC(0),REC=NBUFX(0)) BUF,IBUF,ICOUNT
         ELSE
            CALL QUIT('Error in UN2WRN')
         END IF
C
         IF (IPRINT .GE. 6) THEN
c        if (bph2oo) then
            IF (ICOUNT .GT. 0) THEN
               WRITE (LUPRI,'(2X,A,I5,A/)')
     &            'Integral buffer #',NBUFX(0),' has been written.'
               IBIT1 = 2**NBITS - 1
               DO 200 INT = 1, ICOUNT
                  IF (NIBUF .EQ. 1) THEN
                     I = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
                     J = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
                     K = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
                     L = IBTAND(       IBUF(INT,1),         IBIT1)
                  ELSE
                     I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
                     J = IBTAND(       IBUF(INT,1),       IBIT1)
                     K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
                     L = IBTAND(       IBUF(INT,2),       IBIT1)
                  END IF
                  WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8)')
     &                         ' ## ', I, J, K, L, BUF(INT)
  200          CONTINUE
            END IF
         END IF
C
C        Statistics
C
         IF (INDA .EQ. 0) THEN
C           N2WRIT = LBUF*NBUFX(0) + ICOUNT
C           hjaaj: do all intermediates with integer*8 N2WRIT
C                  to avoid integer*4 overflows
            N2WRIT = LBUF
            N2WRIT = N2WRIT*NBUFX(0)
            N2WRIT = N2WRIT + ICOUNT
C
            IF (ICOUNT.GT.0 .AND. INDA.NE.0) THEN
               NBUFX(0) = NBUFX(0) + 2
            ELSE
               NBUFX(0) = NBUFX(0) + 1
            END IF
            IF (IRAT .EQ. 1) LWORD = 8
            IF (IRAT .EQ. 2) LWORD = 4
            FMBYTES = LWORD*(LBUF*IRAT + NIBUF*LBUF + 1)
            FMBYTES = NBUFX(0) * FMBYTES / (1024.D0**2)
            FNALL  = (NBASIS*(NBASIS + 1))/2
            FNALL  = (FNALL*(FNALL + 1))/2
            PERCNT = N2WRIT
            PERCNT = 100.0D0 * PERCNT / FNALL
            WRITE (LUPRI,'(/A,I10,A,F4.1,A/A,F10.3//)')
     &         ' Number of two-electron integrals written:',N2WRIT,
     &         ' (',PERCNT,'%)',
     &         ' Megabytes written:                       ',FMBYTES
         END IF
      END IF
      RETURN
      END
C  /* Deck sr2out */
      SUBROUTINE SR2OUT(SO,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,THRESH,
     &                  NINDAB,NINDCD,IPRINT)
C
C     Adjusted for writing to LUSRINT for short-range 2-elektron
C     integrals (MCDFT)   Feb. 2003 JKP
C
C     Based on UN2OUT
C        Write out blocks of symmetry integrals, eliminating duplicates
C                                   880412   PRT
C        Rewritten to allow for triangular looping and to eliminate
C        all duplicates by testing  880601 TUH
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
      PARAMETER (LBUF = 600)
      LOGICAL DCMPAB, DCMPCD, DCMPAC, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, IAEQIC, IALTIC, IPNTLG(3,*), NOTEST
      DIMENSION SO(*), BUF(LBUF), IBUF(LBUF,2),
     &          IPNTNO(4,*), IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
#include <nuclei.h>
#include <twocom.h>
#include <symmet.h>
      SAVE BUF, IBUF, ICOUNT
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine SR2OUT',-1)
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
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) CALL SR2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,-1,NBITS,0,
     &                       IPRINT)
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
            IF (NIBUF .EQ. 1) THEN
               INT = 0
               DO 200 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
                  DO 210 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1)=MAX(INDAB,INDCD)*IBIT2+INDAB+INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL SR2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                       END IF
                     END IF
  210             CONTINUE
  200          CONTINUE
            ELSE
               INT = 0
               DO 205 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
                  DO 215 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT = SO(ISOFF+INT)
                     IF (ABS(SOINT) .GT. THRESH) THEN
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1) = MAX(INDAB,INDCD)
                       IBUF(ICOUNT,2) = MIN(INDAB,INDCD)
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL SR2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                       END IF
                     END IF
  215             CONTINUE
  205          CONTINUE
            END IF
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
               INDAB = MAX(INDA,INDB)*IBIT1 + INDA + INDB
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
                  SOINT = SO(ISOFF+INT)
                  IF (ABS(SOINT) .GT. THRESH) THEN
                     IF (NIBUF .EQ. 1) THEN
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1)=MAX(INDAB,INDCD)*IBIT2+INDAB+INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL SR2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                          NBUFCL = NBUFCL + 1
                       END IF
                     ELSE
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       ICOUNT = ICOUNT + 1
                       BUF (ICOUNT) = SOINT
                       IBUF(ICOUNT,1) = MAX(INDAB,INDCD)
                       IBUF(ICOUNT,2) = MIN(INDAB,INDCD)
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL SR2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,0,
     &                                NBITS,0,IPRINT)
                          NBUFCL = NBUFCL + 1
                       END IF
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
      IF (LAST) THEN
         WRITE (LUPRI,'(A,1P,D10.2)') ' Threshold for neglecting '//
     &      'two-electron short-range integrals:',THRESH
         CALL SR2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,1,NBITS,0,IPRINT)
      END IF
      RETURN
      END
C  /* Deck sr2wrt */
      SUBROUTINE SR2WRT(BUF,IBUF,LBUF,NIBUF,ICOUNT,ITYPE,NBITS,INDA,
     &                  IPRINT)
C
C     Write short-range integrals to disk file SRTWOINT
C     (based on UN2WRT)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <ibtpar.h>
      DIMENSION BUF(LBUF), IBUF(LBUF,NIBUF)
#include <twosta.h>
#include <inftap.h>
#include <nuclei.h>
#include <symmet.h>
      SAVE NBUF, JBUF
#if defined (VAR_SPLITFILES)
#include <2gbdef.h>
      CHARACTER*80 FNNAME, FNNM2
#include <chrnos.h>
#endif
#include <ibtfun.h>
C
      IF (ITYPE .EQ. -1) THEN
C        ... First call
         CALL REWSPL(LUSRINT)
         CALL NEWLAB('BASINFO ',LUSRINT,LUPRI)
         WRITE (LUSRINT) MAXREP+1,(NAOS(I),I=1,8),LBUF,NIBUF,NBITS
         CALL NEWLAB('BASTWOEL',LUSRINT,LUPRI)
         JLEN = 608
         ICOUNT = 0
         NBUF = 0
         JBUF = 608
      ELSE IF (ITYPE .EQ. 0) THEN
C        ... between first and last call
#if defined (VAR_SPLITFILES)
C
         IF ((JBUF + LBUF*(IRAT + NIBUF) + 1) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
            INQUIRE(UNIT=LUSRINT,NAME=FNNAME)
            LN = 1
 10         CONTINUE
            IF (FNNAME(LN:LN) .NE. ' ') THEN
               LN = LN + 1
               GOTO 10
            END IF
            LN = LN - 1
            CALL GPCLOSE(LUSRINT,'KEEP')
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
            CALL GPOPEN(LUSRINT,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &           .FALSE.)
            JBUF = 0
         END IF
         IF (INDA .NE. 0) THEN
            WRITE (LUSRINT) INDA
            JBUF = JBUF + 1
         END IF
         WRITE (LUSRINT) BUF,IBUF,ICOUNT
         JBUF = JBUF + LBUF*(IRAT + NIBUF) + 1
#else
         IF (INDA .NE. 0) WRITE (LUSRINT) INDA
         WRITE (LUSRINT) BUF,IBUF,ICOUNT
#endif
         NBUF = NBUF + 1
C
         IF (IPRINT .GE. 6) THEN
            WRITE (LUPRI,'(2X,A,I5,A/)')
     &         'Integral buffer #',NBUF,' has been written.'
            IBIT1 = 2**NBITS - 1
            DO 100 INT = 1, ICOUNT
               IF (NIBUF .EQ. 1) THEN
                  I = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
                  J = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
                  K = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
                  L = IBTAND(       IBUF(INT,1),         IBIT1)
               ELSE
                  I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
                  J = IBTAND(       IBUF(INT,1),       IBIT1)
                  K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
                  L = IBTAND(       IBUF(INT,2),       IBIT1)
               END IF
               WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8)')
     &                      ' ## ', I, J, K, L, BUF(INT)
  100       CONTINUE
         END IF
         ICOUNT = 0
      ELSE
C        ... ITYPE .gt. 0: last call
         IF (INDA .NE. 0) THEN
            WRITE (LUSRINT) -INDA
            WRITE (LUSRINT) BUF,IBUF,ICOUNT
         ELSE
            IF (ICOUNT .GT. 0) WRITE (LUSRINT) BUF,IBUF,ICOUNT
            WRITE (LUSRINT) BUF,IBUF,-1
         END IF
C
         IF (IPRINT .GE. 6) THEN
            IF (ICOUNT .GT. 0) THEN
               WRITE (LUPRI,'(2X,A,I5,A/)')
     &            'Integral buffer #',NBUF+1,' has been written.'
               IBIT1 = 2**NBITS - 1
               DO 200 INT = 1, ICOUNT
                  IF (NIBUF .EQ. 1) THEN
                     I = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
                     J = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
                     K = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
                     L = IBTAND(       IBUF(INT,1),         IBIT1)
                  ELSE
                     I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
                     J = IBTAND(       IBUF(INT,1),       IBIT1)
                     K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
                     L = IBTAND(       IBUF(INT,2),       IBIT1)
                  END IF
                  WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8)')
     &                         ' ## ', I, J, K, L, BUF(INT)
  200          CONTINUE
            END IF
         END IF
C
C        Statistics
C
         IF (INDA .EQ. 0) THEN
C           N2WRIT = LBUF*NBUF + ICOUNT
C           hjaaj: do all intermediates with integer*8 N2WRIT
C                  to avoid integer*4 overflows for big basis sets
            N2WRIT = LBUF
            N2WRIT = N2WRIT*NBUF
            N2WRIT = N2WRIT + ICOUNT
C
            IF (ICOUNT.GT.0 .AND. INDA.NE.0) THEN
               NBUF = NBUF + 2
            ELSE
               NBUF = NBUF + 1
            END IF
            IF (IRAT .EQ. 1) LWORD = 8
            IF (IRAT .EQ. 2) LWORD = 4
            FBYTES = LWORD*(LBUF*IRAT + NIBUF*LBUF + 1)
            FBYTES = NBUF*FBYTES
            FBYTES = FBYTES / (1024.D0**2)
            FNALL  = (NBASIS*(NBASIS + 1))/2
            FNALL  = (FNALL*(FNALL + 1.0D0))*0.5D0
            PERCNT = N2WRIT
            PERCNT = 100.0D0* PERCNT / FNALL
            WRITE (LUPRI,'(/A,I12,A,F5.1,A/A,F12.3/)')
     &         ' Number of two-electron short-range integrals written:',
     &         N2WRIT,' (',PERCNT,'% )',
     &         ' Megabytes written:                       ',FBYTES
         END IF
         CALL GPCLOSE(LUSRINT,'KEEP')
      END IF
      RETURN
      END
