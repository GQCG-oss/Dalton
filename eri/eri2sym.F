C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
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
C
C$Id: eri2sym.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
C  /* Deck erisym */
      SUBROUTINE ERISYM(AO,SO,IODDCC,IPNTUV,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      DIMENSION SO(NCCX,KHKTAB,KHKTCD),
     &          AO(NCCX,KHKTAB,KHKTCD),
     &          RMAT(64,0:7), SMAT(64,0:7), TMAT(64,0:7),
     &          IPNTUV(KC2MAX,2,2), IODDCC(NRTOP)
#include <ericom.h>
#include <symmet.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERI2SYM',-1)
         WRITE (LUPRI,'(1X,A,4I5)')
     &         ' NHKTA,... ',NHKTA,NHKTB,NHKTC,NHKTD
         WRITE (LUPRI,'(1X,A,4I5)')
     &         ' ISTBLA,...',ISTBLA,ISTBLB,ISTBLC,ISTBLD
         WRITE (LUPRI,'(1X,A,4I5)')
     &         ' MLTPA,... ',MLTPA,MLTPB,MLTPC,MLTPD
         WRITE (LUPRI,'(1X,A,3I5)')
     &         ' ISTBLR,...',ISTBLR,ISTBLS,ISTBLT
         WRITE (LUPRI,'(1X,A,3I5)')
     &         ' MLTPR,... ', MLTPR,MLTPS,MLTPT
         WRITE (LUPRI,'(1X,A,2I5)')
     &         ' NCCS, NCCX', NCCS, NCCX
         IF (IPRINT .GT. 15) THEN
            CALL HEADER('AO and SO integrals in ERI2SYM',-1)
         END IF
      END IF
C
      CALL ERITRA(RMAT,MLTPR,ISTBLR)
      CALL ERITRA(SMAT,MLTPS,ISTBLS)
      CALL ERITRA(TMAT,MLTPT,ISTBLT)
C
C     Run over components
C     ===================
C
      MAXB = KHKTB
      MAXD = KHKTD
C
      ICMPAB = 0
      DO 300 ICOMPA = 1, KHKTA
         IF (TKMPAB) MAXB = ICOMPA
      DO 300 ICOMPB = 1, MAXB
         ICMPAB = ICMPAB + 1
C
         ICMPCD = 0
      DO 300 ICOMPC = 1, KHKTC
         IF (TKMPCD) MAXD = ICOMPC
      DO 300 ICOMPD = 1, MAXD
         ICMPCD = ICMPCD + 1
C
         IODDAB = IODDCC(IPNTUV(ICMPAB,2,1))
         IODDCD = IODDCC(IPNTUV(ICMPCD,2,2))
         IF (IODDAB .EQ. IODDCD) THEN
C
C           Cartesian phase factors
C           =======================
C
            IPARB  = ISYMAO(NHKTB,ICOMPB)
            IPARC  = ISYMAO(NHKTC,ICOMPC)
            IPARD  = ISYMAO(NHKTD,ICOMPD)
C
            CALL DGEMM('T','N',MLTPS*MLTPT*NCCS,MLTPR,MLTPR,
     &                 D1,AO(1,ICMPAB,ICMPCD),MLTPR,
     &                 RMAT(1,IPARB),MLTPR,
     &                 D0,SO(1,ICMPAB,ICMPCD),MLTPS*MLTPT*NCCS)
            CALL DGEMM('T','N',MLTPT*NCCS*MLTPR,MLTPS,MLTPS,
     &                 D1,SO(1,ICMPAB,ICMPCD),MLTPS,
     &                 SMAT(1,IPARD),MLTPS,
     &                 D0,AO(1,ICMPAB,ICMPCD),MLTPT*NCCS*MLTPR)
            CALL DGEMM('T','N',NCCS*MLTPR*MLTPS,MLTPT,MLTPT,
     &                 D1,AO(1,ICMPAB,ICMPCD),MLTPT,
     &                 TMAT(1,IBTXOR(IPARC,IPARD)),MLTPT,
     &                 D0,SO(1,ICMPAB,ICMPCD),NCCS*MLTPR*MLTPS)
C
C           Print
C           =====
C
            IF (IPRINT .GT. 15) THEN
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' AO integrals for components ',ICMPAB,ICMPCD
               CALL OUTPUT(AO(1,ICMPAB,ICMPCD),1,NCCS,1,MLTPX,
     &                     NCCS,MLTPX,1,LUPRI)
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' SO integrals for components ',ICMPAB,ICMPCD
               CALL OUTPUT(SO(1,ICMPAB,ICMPCD),1,NCCS,1,MLTPX,
     &                     NCCS,MLTPX,1,LUPRI)
            END IF
C
         END IF
C
  300 CONTINUE
C
      RETURN
      END
C  /* Deck eritra */
      SUBROUTINE ERITRA(RMAT,MLTPX,ISTBLX)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION RMAT(64,0:7)
#include <symmet.h>
#include <ibtfun.h>
C
      IRED = 0
      DO 100 I = 0, MAXOPR
      IF (IBTAND(I,ISTBLX) .EQ. 0) THEN
         IRED = IRED + 1
         JRED = 0
         DO 200 J = 0, MAXOPR
         IF (IBTAND(J,ISTBLX) .EQ. 0) THEN
            JRED = JRED + 1
            JI = (IRED - 1)*MLTPX + JRED
            DO 300 K = 0, MAXREP
               RMAT(JI,K) = PT(IBTAND(J,IBTXOR(I,K)))
  300       CONTINUE
         END IF
  200    CONTINUE
      END IF
  100 CONTINUE
C
      RETURN
      END
