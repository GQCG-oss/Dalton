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
C$Id: eri2sph.F,v 1.2 2001-02-12 18:17:55 vebjornb Exp $
C
#include <single.h>
C
C  /* Deck cr1sph */
      SUBROUTINE CR1SPH(HCCONT,HCINT,FCINT,HCSINT,FCSINT,CSQ1,CSQ2,
     &                  IODDHC,IPNTUV,IODDKC,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER TUV
      LOGICAL   FCSINT(NTUV34,KHKT2),
     &          FCINT (NTUV34,KHKT12)
      DIMENSION HCCONT(NCCPP,NTUV34,KCKT12),
     &          HCSINT(NCCPP,NTUV34,KHKT2),
     &          HCINT (NCCPP,NTUV34,KHKT12),
     &          CSQ1(KHKT1,KCKT1),
     &          CSQ2(KHKT2,KCKT2),
     &          IODDHC(NRTOP),
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          IODDKC(KC2MAX)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR1SPH','*',103)
C
      DO 10 TUV = 1, NTUV34
      DO 10 IKMP12 = 1, KHKT12
         FCINT(TUV,IKMP12) = .TRUE.
   10 CONTINUE
C
C     Transformation of both indices
C     ==============================
C
      IF (SPHR1 .AND. SPHR2) THEN
         DO 100 ICOMP1 = 1, KCKT1
C
C           First half transformation:
C
            DO 110 TUV = 1, NTUV34
            DO 110 IKOMP2 = 1, KHKT2
               FCSINT(TUV,IKOMP2) = .TRUE.
  110       CONTINUE
C
            DO 200 ICOMP2 = 1, KCKT2
               ICMP12 = (ICOMP1 - 1)*KCKT2 + ICOMP2
               IODD12 = IODDKC(ICMP12)
               DO 210 IKOMP2 = 1, KHKT2
                  SPHFAC = CSQ2(IKOMP2,ICOMP2)
                  IF (ABS(SPHFAC) .GT. D0) THEN
                     DO 220 TUV = 1, NTUV34
                     IF (IODDHC(TUV) .EQ. IODD12) THEN
                        IF (FCSINT(TUV,IKOMP2)) THEN
                           FCSINT(TUV,IKOMP2) = .FALSE.
                           DO 230 I = 1, NCCPP
                              HCSINT(I,TUV,IKOMP2)=
     &                                      SPHFAC*HCCONT(I,TUV,ICMP12)
  230                      CONTINUE
                        ELSE
                           DO 235 I = 1, NCCPP
                              HCSINT(I,TUV,IKOMP2)=HCSINT(I,TUV,IKOMP2)
     &                                    + SPHFAC*HCCONT(I,TUV,ICMP12)
  235                      CONTINUE
                        END IF
                     END IF
  220                CONTINUE
                  END IF
  210          CONTINUE
  200       CONTINUE
C
C           Second half transformation:
C
            IKMP12 = 0
            DO 300 IKOMP1 = 1, KHKT1
               SPHFAC = CSQ1(IKOMP1,ICOMP1)
               IF (ABS(SPHFAC) .GT. D0) THEN
                  MAX2 = KHKT2
                  IF (TKMP12) MAX2 = IKOMP1
                  DO 310 IKOMP2 = 1, MAX2
                     IKMP12 = IKMP12 + 1
                     IODD12 = IODDHC(IPNTUV(IKMP12,0,IELCT1))
                     DO 320 TUV = 1, NTUV34
                     IF (IODDHC(TUV) .EQ. IODD12) THEN
                        IF (FCINT(TUV,IKMP12)) THEN
                           FCINT(TUV,IKMP12) = .FALSE.
                           DO 330 I = 1, NCCPP
                              HCINT(I,TUV,IKMP12) =
     &                                      SPHFAC*HCSINT(I,TUV,IKOMP2)
  330                      CONTINUE
                        ELSE
                           DO 335 I = 1, NCCPP
                              HCINT(I,TUV,IKMP12) = HCINT(I,TUV,IKMP12)
     &                                    + SPHFAC*HCSINT(I,TUV,IKOMP2)
  335                      CONTINUE
                        END IF
                     END IF
  320                CONTINUE
  310             CONTINUE
               ELSE IF (TKMP12) THEN
                  IKMP12 = IKMP12 + IKOMP1
               ELSE
                  IKMP12 = IKMP12 + KHKT2
               END IF
  300       CONTINUE
C
  100    CONTINUE
C
C     Transformation of first index only
C     ==================================
C
      ELSE IF (SPHR1) THEN
         DO 400 ICOMP1 = 1, KCKT1
         DO 400 IKOMP1 = 1, KHKT1
            SPHFAC = CSQ1(IKOMP1,ICOMP1)
            IF (ABS(SPHFAC) .GT. D0) THEN
            DO 410 IKOMP2 = 1, KHKT2
               ICMP12 = (ICOMP1 - 1)*KHKT2 + IKOMP2
               IKMP12 = (IKOMP1 - 1)*KHKT2 + IKOMP2
               IODD12 = IODDHC(IPNTUV(IKMP12,0,IELCT1))
               DO 420 TUV = 1, NTUV34
               IF (IODDHC(TUV) .EQ. IODD12) THEN
                  IF (FCINT(TUV,IKMP12)) THEN
                     FCINT(TUV,IKMP12) = .FALSE.
                     DO 430 I = 1, NCCPP
                        HCINT(I,TUV,IKMP12) =
     &                                SPHFAC*HCCONT(I,TUV,ICMP12)
  430                CONTINUE
                  ELSE
                     DO 435 I = 1, NCCPP
                        HCINT(I,TUV,IKMP12) = HCINT(I,TUV,IKMP12)
     &                              + SPHFAC*HCCONT(I,TUV,ICMP12)
  435                CONTINUE
                  END IF
               END IF
  420          CONTINUE
  410       CONTINUE
            END IF
  400    CONTINUE
C
C     Transformation of second index only
C     ===================================
C
      ELSE
         DO 500 ICOMP2 = 1, KCKT2
         DO 500 IKOMP2 = 1, KHKT2
            SPHFAC = CSQ2(IKOMP2,ICOMP2)
            IF (ABS(SPHFAC) .GT. D0) THEN
            DO 510 IKOMP1 = 1, KHKT1
               ICMP12 = (IKOMP1 - 1)*KCKT2 + ICOMP2
               IKMP12 = (IKOMP1 - 1)*KHKT2 + IKOMP2
               IODD12 = IODDHC(IPNTUV(IKMP12,0,IELCT1))
               DO 520 TUV = 1, NTUV34
               IF (IODDHC(TUV) .EQ. IODD12) THEN
                  IF (FCINT(TUV,IKMP12)) THEN
                     FCINT(TUV,IKMP12) = .FALSE.
                     DO 530 I = 1, NCCPP
                        HCINT(I,TUV,IKMP12) =
     &                                SPHFAC*HCCONT(I,TUV,ICMP12)
  530                CONTINUE
                  ELSE
                     DO 535 I = 1, NCCPP
                        HCINT(I,TUV,IKMP12) = HCINT(I,TUV,IKMP12)
     &                              + SPHFAC*HCCONT(I,TUV,ICMP12)
  535                CONTINUE
                  END IF
               END IF
  520          CONTINUE
  510       CONTINUE
            END IF 
  500    CONTINUE 
      END IF
      RETURN
      END
C  /* Deck cr2sph */
      SUBROUTINE CR2SPH(CCONT,CSINT,FCSINT,AOINT,FAOINT,FLEVEL,
     &                  CSQ3,CSQ4,CSA3,CSA4,CSV3,CSV4,BDER3,BDER4,
     &                  NCENTR,IODDCC,IPNTUV,IODDKC,INDX,IODX,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      LOGICAL   FCSINT(KHKT12,KHKT4),
     &          FAOINT(KHKT12,KHKT34),
     &          FLEVEL, BDER3, BDER4
      DIMENSION CCONT(NCCCC,KHKT12,KCKT34),
     &          CSINT(NCCCC,KHKT12,KHKT4),
     &          AOINT(NCCT,KHKTAB,KHKTCD),
     &          CSQ3(NCSQ1,NCSQ2), CSQ4(NCSQ1,NCSQ2),
     &          CSA3(KHKT3,KCKT3), CSA4(KHKT4,KCKT4),
     &          CSV3(NCCCC,KHKT3,KCKT3), CSV4(NCCCC,KHKT4,KCKT4),
     &          IODDCC(NRTOP), NCENTR(NCCT,2),
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          IODDKC(KC2MAX)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (BDER3.OR.BDER4) THEN
         CALL CR2SVC(CSQ3,CSA3,CSV3,NCENTR(1,1),NCENTR(1,2),
     &               NCCCC,KHKT3,KCKT3)
         CALL CR2SPX(CCONT,CSINT,FCSINT,AOINT,FAOINT,FLEVEL,
     &               CSA3,CSQ4,CSV3,CSV4,BDER3,BDER4,
     &               IODDCC,IPNTUV,IODDKC,INDX,IODX,IPRINT)
      ELSE 
         CALL CR2SPX(CCONT,CSINT,FCSINT,AOINT,FAOINT,FLEVEL,
     &               CSQ3,CSQ4,CSV3,CSV4,BDER3,BDER4,
     &               IODDCC,IPNTUV,IODDKC,INDX,IODX,IPRINT)
      END IF
      RETURN
      END
C  /* Deck cr2svc */
      SUBROUTINE CR2SVC(CSQ,CSA,CSV,NCNT1,NCNT2,NCCCC,KHKTX,KCKTX)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION CSQ(NCSQ1,NCSQ2), CSA(KHKTX,KCKTX), 
     &          CSV(NCCCC,KHKTX,KCKTX), NCNT1(*), NCNT2(*) 
#include <ericom.h>
#include <nuclei.h>
      IADR = 1
      DO 100 ICMP = 1, KCKTX
      DO 100 IKMP = 1, KHKTX
         CSQM = D0
         DO 200 I = 1, NCCCC
            CSQI = CSQ(IADR,NCNT2(I)) - CSQ(IADR,NCNT1(I))
            CSQM = MAX(CSQM,ABS(CSQI)) 
            CSV(I,IKMP,ICMP) = CSQI 
  200    CONTINUE
         CSA(IKMP,ICMP) = CSQM
         IADR = IADR + 1
  100 CONTINUE
      RETURN
      END
C  /* Deck cr2sph */
      SUBROUTINE CR2SPX(CCONT,CSINT,FCSINT,AOINT,FAOINT,FLEVEL,
     &                  CSQ3,CSQ4,CSV3,CSV4,BDER3,BDER4,
     &                  IODDCC,IPNTUV,IODDKC,INDX,IODX,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      LOGICAL   FCSINT(KHKT12,KHKT4),
     &          FAOINT(KHKT12,KHKT34),
     &          FLEVEL, PATH1, BDER3, BDER4
      DIMENSION CCONT(NCCCC,KHKT12,KCKT34),
     &          CSINT(NCCCC,KHKT12,KHKT4),
     &          AOINT(NCCT,KHKTAB,KHKTCD),
     &          CSQ3(KHKT3,KCKT3), CSQ4(KHKT4,KCKT4),
     &          CSV3(NCCCC,KHKT3,KCKT3), CSV4(NCCCC,KHKT4,KCKT4),
     &          IODDCC(NRTOP),
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          IODDKC(KC2MAX)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      PATH1 = IPATH .EQ. 1
C
C     IF (FLEVEL) THEN
         DO 10 IKMP12 = 1, KHKT12
         DO 10 IKMP34 = 1, KHKT34
            FAOINT(IKMP12,IKMP34) = .TRUE.
            FAOINT(IKMP12,IKMP34) = .FALSE.
   10    CONTINUE
C     END IF
C
C     Transformation of both indices
C     ==============================
C
      CALL DZERO(CSINT,NCCCC*KHKT12*KHKT4)
      IF (SPHR3 .AND. SPHR4) THEN
         ICMP34 = 0
         DO 100 ICOMP3 = 1,KCKT3
C
C           First half transformation
C
            DO 110 IKMP12 = 1, KHKT12
            DO 110 IKOMP4 = 1, KHKT4
               FCSINT(IKMP12,IKOMP4) = .TRUE.
  110       CONTINUE
            DO 200 ICOMP4 = 1, KCKT4
               ICMP34 = ICMP34 + 1
               IODD34 = IODDKC(ICMP34)
               DO 210 IKOMP4 = 1, KHKT4
                  SPHFAC = CSQ4(IKOMP4,ICOMP4)
                  IF (ABS(SPHFAC).GT.D0) THEN
                     DO 220 IKMP12 = 1, KHKT12
                     IF(IODDCC(IPNTUV(IKMP12,0,IELCT1)).EQ.IODD34)THEN
                        IF (FCSINT(IKMP12,IKOMP4)) THEN
                           FCSINT(IKMP12,IKOMP4) = .FALSE.
                           IF (.NOT.BDER4) THEN
                              DO 230 I = 1, NCCCC
                                 CSINT(I,IKMP12,IKOMP4) =
     &                              SPHFAC*CCONT(I,IKMP12,ICMP34)
  230                         CONTINUE
                           ELSE
                              DO 231 I = 1, NCCCC
                                 CSINT(I,IKMP12,IKOMP4) = 
     &                               CSV4(I,IKOMP4,ICOMP4)
     &                                   *CCONT(I,IKMP12,ICMP34)
  231                         CONTINUE
                           END IF
                        ELSE
                           IF (.NOT.BDER4) THEN
                              DO 235 I = 1, NCCCC
                                 CSINT(I,IKMP12,IKOMP4) =
     &                                 CSINT(I,IKMP12,IKOMP4)
     &                                 + SPHFAC*CCONT(I,IKMP12,ICMP34)
  235                         CONTINUE
                           ELSE
                              DO 236 I = 1, NCCCC
                                 CSINT(I,IKMP12,IKOMP4) =
     &                                 CSINT(I,IKMP12,IKOMP4)
     &                                    + CSV4(I,IKOMP4,ICOMP4)
     &                                         *CCONT(I,IKMP12,ICMP34)
  236                         CONTINUE
                           END IF 
                        END IF
                     END IF
  220                CONTINUE
                  END IF
  210          CONTINUE
  200       CONTINUE
C
C           Second half transformation
C
            IKMP34 = 0
            DO 300 IKOMP3 = 1, KHKT3
               SPHFAC = CSQ3(IKOMP3,ICOMP3)
               IF (ABS(SPHFAC) .GT. D0) THEN
                  MAX4 = KHKT4
                  IF (TKMP34) MAX4 = IKOMP3
                  DO 310 IKOMP4 = 1, MAX4
                     IKMP34 = IKMP34 + 1
C                    IODD34 = IODDCC(IPNTUV(IKMP34,INDX,IELCT2))
                     IODD34 = IODDCC(IPNTUV(IKMP34,0,IELCT2))
                     IODD34 = IBTXOR(IODD34,IODX)
                     DO 320 IKMP12 = 1, KHKT12
                     IF (IODDCC(IPNTUV(IKMP12,0,IELCT1)).EQ.IODD34) THEN
                        IF (PATH1) THEN
                          IKMPAB = IKMP12
                          IKMPCD = IKMP34
                        ELSE
                          IKMPAB = IKMP34
                          IKMPCD = IKMP12
                        END IF
                        IF (FAOINT(IKMP12,IKMP34)) THEN
                           FAOINT(IKMP12,IKMP34) = .FALSE.
                           IF (.NOT.BDER3) THEN
                              DO 330 I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                   SPHFAC*CSINT(I,IKMP12,IKOMP4)
  330                         CONTINUE
                           ELSE
                              DO 331 I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                              CSV3(I,IKOMP3,ICOMP3)
     &                                 *CSINT(I,IKMP12,IKOMP4)
  331                         CONTINUE
                           END IF
                        ELSE
                           IF (.NOT.BDER3) THEN
                              DO 335 I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                 AOINT(I,IKMPAB,IKMPCD)
     &                                 + SPHFAC*CSINT(I,IKMP12,IKOMP4)
  335                         CONTINUE
                           ELSE
                              DO 336 I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                              AOINT(I,IKMPAB,IKMPCD)
     &                                 + CSV3(I,IKOMP3,ICOMP3)
     &                                      *CSINT(I,IKMP12,IKOMP4)
  336                         CONTINUE
                           END IF
                        END IF
                     END IF
  320                CONTINUE
  310             CONTINUE
               ELSE IF (TKMP34) THEN
                  IKMP34 = IKMP34 + IKOMP3
               ELSE
                  IKMP34 = IKMP34 + KHKT4
               END IF
  300       CONTINUE
  100    CONTINUE
C
C     Transformation of first index only
C     ==================================
C
      ELSE IF (SPHR3) THEN
         DO 400 ICOMP3 = 1, KCKT3
         DO 400 IKOMP3 = 1, KHKT3
            SPHFAC = CSQ3(IKOMP3,ICOMP3)
            IF (ABS(SPHFAC) .GT. D0) THEN
               DO 420 IKOMP4 = 1, KHKT4
                  ICMP34 = (ICOMP3 - 1)*KHKT4 + IKOMP4
                  IKMP34 = (IKOMP3 - 1)*KHKT4 + IKOMP4
                  IODD34 = IODDCC(IPNTUV(IKMP34,INDX,IELCT2))
                  IODD34 = IBTXOR(IODD34,IODX)
                  DO 430 IKMP12 = 1, KHKT12
                  IF (IODDCC(IPNTUV(IKMP12,0,IELCT1)).EQ.IODD34) THEN 
                     IF (PATH1) THEN
                       IKMPAB = IKMP12
                       IKMPCD = IKMP34
                     ELSE
                       IKMPAB = IKMP34
                       IKMPCD = IKMP12
                     END IF
                     IF (FAOINT(IKMP12,IKMP34)) THEN
                        FAOINT(IKMP12,IKMP34) = .FALSE.
                        IF (.NOT.BDER3) THEN
                           DO 440 I = 1, NCCCC
                              AOINT(I,IKMPAB,IKMPCD) =
     &                                     SPHFAC*CCONT(I,IKMP12,ICMP34)
  440                      CONTINUE
                        ELSE
                           DO 441 I = 1, NCCCC
                              AOINT(I,IKMPAB,IKMPCD) =
     &                           CSV3(I,IKOMP3,ICOMP3)
     &                              *CCONT(I,IKMP12,ICMP34)
  441                      CONTINUE
                        END IF
                     ELSE
                        IF (.NOT.BDER3) THEN
                           DO 445 I = 1, NCCCC
                              AOINT(I,IKMPAB,IKMPCD) =
     &                                            AOINT(I,IKMPAB,IKMPCD)
     &                                   + SPHFAC*CCONT(I,IKMP12,ICMP34)
  445                      CONTINUE
                        ELSE
                           DO 446 I = 1, NCCCC
                              AOINT(I,IKMPAB,IKMPCD) =
     &                           AOINT(I,IKMPAB,IKMPCD)
     &                              + CSV3(I,IKOMP3,ICOMP3)
     &                                 *CCONT(I,IKMP12,ICMP34)
  446                      CONTINUE
                        END IF
                     END IF
                  END IF
  430             CONTINUE
  420          CONTINUE
            END IF
  400    CONTINUE
C
C     Transformation of second index only
C     ===================================
C
      ELSE
         DO 500 IKOMP3 = 1, KHKT3
         DO 500 ICOMP4 = 1, KCKT4
            ICMP34 = (IKOMP3 - 1)*KCKT4 + ICOMP4
               DO 510 IKOMP4 = 1, KHKT4
                  SPHFAC = CSQ4(IKOMP4,ICOMP4)
                  IF (ABS(SPHFAC).GT.D0) THEN
                     IKMP34 = (IKOMP3 - 1)*KHKT4 + IKOMP4
                     IODD34 = IODDCC(IPNTUV(IKMP34,INDX,IELCT2))
                     IODD34 = IBTXOR(IODD34,IODX)
                     DO 520 IKMP12 = 1, KHKT12
                     IF (IODDCC(IPNTUV(IKMP12,0,IELCT1)).EQ.IODD34) THEN
                        IF (PATH1) THEN
                          IKMPAB = IKMP12
                          IKMPCD = IKMP34
                        ELSE
                          IKMPAB = IKMP34
                          IKMPCD = IKMP12
                        END IF
                        IF (FAOINT(IKMP12,IKMP34)) THEN
                           FAOINT(IKMP12,IKMP34) = .FALSE.
                           IF (.NOT.BDER4) THEN
                              DO 530 I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                   SPHFAC*CCONT(I,IKMP12,ICMP34)
  530                         CONTINUE
                           ELSE
                              DO 531 I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                              CSV4(I,IKOMP4,ICOMP4)
     &                                 *CCONT(I,IKMP12,ICMP34)
  531                         CONTINUE
                           END IF 
                        ELSE
                           IF (.NOT.BDER4) THEN
                              DO 535 I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                          AOINT(I,IKMPAB,IKMPCD)
     &                                 + SPHFAC*CCONT(I,IKMP12,ICMP34)
  535                         CONTINUE
                           ELSE
                              DO 536 I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                          AOINT(I,IKMPAB,IKMPCD)
     &                              + CSV4(I,IKOMP4,ICOMP4)
     &                                  *CCONT(I,IKMP12,ICMP34)
  536                         CONTINUE
                           END IF
                        END IF
                     END IF
  520                CONTINUE
                  END IF
  510          CONTINUE
  500    CONTINUE
      END IF
      RETURN
      END
C  /* Deck eribun */
      SUBROUTINE ERIBUN(AOINT,AOINT0,COORAO,IODDCC,IPNTUV,R,LWORK,
     &                  IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.5D0)
      DIMENSION AOINT(NCCT,KHKTAB,KHKTCD,6), AOINT0(NCCT,KHKTAB,KHKTCD), 
     &          COORAO(NCCT,3,4),
     &          R(NCCT,6), IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          ORG(3)

#include <ericom.h>
#include <orgcom.h>
#include <eriao.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from ERIBUN','*',103)
      ORG(1) = ORIGIN(1)
      ORG(2) = ORIGIN(2)
      ORG(3) = ORIGIN(3)
C
C     IF (JMAX.EQ.1 .OR. LWORK.LT.6*NCCT) THEN
      IF (.FALSE.) THEN
         DO 100 I = 1, NCCT
            AX = COORAO(I,1,1) - ORG(1)
            AY = COORAO(I,2,1) - ORG(2)
            AZ = COORAO(I,3,1) - ORG(3)
            BX = COORAO(I,1,2) - ORG(1) 
            BY = COORAO(I,2,2) - ORG(2)
            BZ = COORAO(I,3,2) - ORG(3)
            CX = COORAO(I,1,3) - ORG(1)
            CY = COORAO(I,2,3) - ORG(2)
            CZ = COORAO(I,3,3) - ORG(3)
            DX = COORAO(I,1,4) - ORG(1)
            DY = COORAO(I,2,4) - ORG(2)
            DZ = COORAO(I,3,4) - ORG(3)
            R(I,1) = DP5*(AY*BZ - AZ*BY)
            R(I,2) = DP5*(AZ*BX - AX*BZ)
            R(I,3) = DP5*(AX*BY - AY*BX)
            R(I,4) = DP5*(CY*DZ - CZ*DY)
            R(I,5) = DP5*(CZ*DX - CX*DZ)
            R(I,6) = DP5*(CX*DY - CY*DX)
  100    CONTINUE
         DO 200 IAB = 1, KHKTAB 
         DO 200 ICD = 1, KHKTCD 
         IF (IODDCC(IPNTUV(IAB,0,1)).EQ.IODDCC(IPNTUV(ICD,0,2))) THEN
            DO 300 I = 1, NCCT
               AO = AOINT0(I,IAB,ICD)
               AOINT(I,IAB,ICD,1) = AO*R(I,1)
               AOINT(I,IAB,ICD,2) = AO*R(I,2)
               AOINT(I,IAB,ICD,3) = AO*R(I,3)
               AOINT(I,IAB,ICD,4) = AO*R(I,4)
               AOINT(I,IAB,ICD,5) = AO*R(I,5)
               AOINT(I,IAB,ICD,6) = AO*R(I,6)
  300       CONTINUE
         END IF
  200    CONTINUE
      ELSE
         DO 400 IAB = 1, KHKTAB 
         DO 400 ICD = 1, KHKTCD 
         IF (IODDCC(IPNTUV(IAB,0,1)).EQ.IODDCC(IPNTUV(ICD,0,2))) THEN
            DO 500 I = 1, NCCT
               AO = DP5*AOINT0(I,IAB,ICD)
               AX = COORAO(I,1,1) - ORG(1) 
               AY = COORAO(I,2,1) - ORG(2) 
               AZ = COORAO(I,3,1) - ORG(3) 
               BX = COORAO(I,1,2) - ORG(1)
               BY = COORAO(I,2,2) - ORG(2)
               BZ = COORAO(I,3,2) - ORG(3)
               CX = COORAO(I,1,3) - ORG(1) 
               CY = COORAO(I,2,3) - ORG(2) 
               CZ = COORAO(I,3,3) - ORG(3) 
               DX = COORAO(I,1,4) - ORG(1)
               DY = COORAO(I,2,4) - ORG(2)
               DZ = COORAO(I,3,4) - ORG(3)
               AOINT(I,IAB,ICD,1) = (AY*BZ - AZ*BY)*AO 
               AOINT(I,IAB,ICD,2) = (AZ*BX - AX*BZ)*AO 
               AOINT(I,IAB,ICD,3) = (AX*BY - AY*BX)*AO 
               AOINT(I,IAB,ICD,4) = (CY*DZ - CZ*DY)*AO 
               AOINT(I,IAB,ICD,5) = (CZ*DX - CX*DZ)*AO 
               AOINT(I,IAB,ICD,6) = (CX*DY - CY*DX)*AO 
  500       CONTINUE
         END IF
  400    CONTINUE
      END IF
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('Final spherical integrals - ERIBUN',-1)
         DO 900 K = 1, 6
            INDX = MOD(K-1,3)+1
            DO 910 I = 1, KHKTAB
            DO 910 J = 1, KHKTCD
C           IF (IODDCC(IPNTUV(I,0,1)) .EQ. 
C    &          IODDCC(IPNTUV(J,INDX,2))) THEN
               WRITE (LUPRI, '(/1X,A,3I3)') ' K,ICMPAB, ICMPCD ', K,I,J
               CALL OUTPUT(AOINT(1,I,J,K),1,1,1,NCCCC,1,NCCCC,1,LUPRI)
C           END IF
  910       CONTINUE
  900    CONTINUE
      END IF
      RETURN
      END
