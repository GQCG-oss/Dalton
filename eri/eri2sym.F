C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
C  /* Deck erisym */
      SUBROUTINE ERISYM(AO,SO,IATOM,ICRB,IREPE,IODDCC,IPNTUV,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      INTEGER R, S, T
      DIMENSION AO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD),
     &          SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD),
     &          RMAT(64,0:7), SMAT(64,0:7), TMAT(64,0:7),
     &          IPNTUV(KC2MAX,0:NRDER,2), IODDCC(NRTOP)
      DIMENSION TMP(NCCS,MLTPR,MLTPS,MLTPT)
#include <ericom.h>
#include <maxorb.h>
#include <symmet.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERISYM',-1)
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
      END IF
C
      CALL ERITRA(RMAT,MLTPR,ISTBLR)
      CALL ERITRA(SMAT,MLTPS,ISTBLS)
      CALL ERITRA(TMAT,MLTPT,ISTBLT)
C
      IPARE = 0
      IF (GDER.AND.ICRB.GT.0) IPARE = IBTXOR(IREPE,ISYMAX(ICRB,1))
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
         IODDAB = IODDCC(IPNTUV(ICMPAB,ICRB,1))
C
         ICMPCD = 0
      DO 300 ICOMPC = 1, KHKTC
         IF (TKMPCD) MAXD = ICOMPC
      DO 300 ICOMPD = 1, MAXD
         ICMPCD = ICMPCD + 1
C
         IF (IODDAB.EQ. IODDCC(IPNTUV(ICMPCD,0,2))) THEN
            IF (IPRINT .GT. 15) THEN
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' AO integrals for components ',ICMPAB,ICMPCD
               CALL OUTPUT(AO(1,1,1,1,ICMPAB,ICMPCD),1,NCCS,1,MLTPX,
     &                     NCCS,MLTPX,1,LUPRI)
            END IF
C
C           Cartesian phase factors
C           =======================
C
            IPARB  = ISYMAO(NHKTB,ICOMPB)
            IPARC  = ISYMAO(NHKTC,ICOMPC)
            IPARD  = ISYMAO(NHKTD,ICOMPD)
            IF (GDER) THEN
               IF (IATOM.EQ.2) IPARB = IBTXOR(IPARB,IPARE)
               IF (IATOM.EQ.3) IPARC = IBTXOR(IPARC,IPARE)
               IF (IATOM.EQ.4) IPARD = IBTXOR(IPARD,IPARE)
            END IF
C
C           Transform T
C
            CALL DGEMM('N','N',NCCS*MLTPR*MLTPS,MLTPT,MLTPT,
     &                 D1,AO(1,1,1,1,ICMPAB,ICMPCD),NCCS*MLTPR*MLTPS,
     &                 TMAT(1,IBTXOR(IPARC,IPARD)),MLTPT,
     &                 D0,SO(1,1,1,1,ICMPAB,ICMPCD),NCCS*MLTPR*MLTPS)
C
C           Transform S
C
            DO 400 T = 1, MLTPT
               CALL DGEMM('N','N',NCCS*MLTPR,MLTPS,MLTPS,
     &                    D1,SO(1,1,1,T,ICMPAB,ICMPCD),NCCS*MLTPR,
     &                    SMAT(1,IPARD),MLTPS,
     &                    D0,TMP(1,1,1,T),NCCS*MLTPR)
  400       CONTINUE
C
C           Transform R
C
            DO 500 T = 1, MLTPT
            DO 500 S = 1, MLTPS
               CALL DGEMM('N','N',NCCS,MLTPR,MLTPR,
     &                    D1,TMP(1,1,S,T),NCCS,
     &                    RMAT(1,IPARB),MLTPR,
     &                    D0,SO(1,1,S,T,ICMPAB,ICMPCD),NCCS)
  500       CONTINUE
C
C#if defined (SYS_CRAY)
C            CALL SGEMM('T','N',MLTPS*MLTPT*NCCS,MLTPR,MLTPR,
C     &                 D1,AO(1,ICMPAB,ICMPCD),MLTPR,
C     &                 RMAT(1,IPARB),MLTPR,
C     &                 D0,SO(1,ICMPAB,ICMPCD),MLTPS*MLTPT*NCCS)
C            CALL SGEMM('T','N',MLTPT*NCCS*MLTPR,MLTPS,MLTPS,
C     &                 D1,SO(1,ICMPAB,ICMPCD),MLTPS,
C     &                 SMAT(1,IPARD),MLTPS,
C     &                 D0,AO(1,ICMPAB,ICMPCD),MLTPT*NCCS*MLTPR)
C            CALL SGEMM('T','N',NCCS*MLTPR*MLTPS,MLTPT,MLTPT,
C     &                 D1,AO(1,ICMPAB,ICMPCD),MLTPT,
C     &                 TMAT(1,IBTXOR(IPARC,IPARD)),MLTPT,
C     &                 D0,SO(1,ICMPAB,ICMPCD),NCCS*MLTPR*MLTPS)
C#endif
C#if defined (SYS_DEC)
C            CALL QUIT('This will require ESSL (eri2sym)')
C#endif
C#if !defined (SYS_CRAY) && !defined (SYS_DEC)
C            CALL DGEMM('T','N',MLTPS*MLTPT*NCCS,MLTPR,MLTPR,
C     &                 D1,AO(1,ICMPAB,ICMPCD),MLTPR,
C     &                 RMAT(1,IPARB),MLTPR,
C     &                 D0,SO(1,ICMPAB,ICMPCD),MLTPS*MLTPT*NCCS)
C            CALL DGEMM('T','N',MLTPT*NCCS*MLTPR,MLTPS,MLTPS,
C     &                 D1,SO(1,ICMPAB,ICMPCD),MLTPS,
C     &                 SMAT(1,IPARD),MLTPS,
C     &                 D0,AO(1,ICMPAB,ICMPCD),MLTPT*NCCS*MLTPR)
C            CALL DGEMM('T','N',NCCS*MLTPR*MLTPS,MLTPT,MLTPT,
C     &                 D1,AO(1,ICMPAB,ICMPCD),MLTPT,
C     &                 TMAT(1,IBTXOR(IPARC,IPARD)),MLTPT,
C     &                 D0,SO(1,ICMPAB,ICMPCD),NCCS*MLTPR*MLTPS)
C#endif
C
C           Print
C           =====
C
            IF (IPRINT .GT. 15) THEN
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' SO integrals for components ',ICMPAB,ICMPCD
               CALL OUTPUT(SO(1,1,1,1,ICMPAB,ICMPCD),1,NCCS,1,MLTPX,
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
#include <mxcent.h>
      DIMENSION RMAT(64,0:7)
#include <maxorb.h>
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
