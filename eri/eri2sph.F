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
C$Id: eri2sph.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
C  /* Deck cr1sph */
      SUBROUTINE CR1SPH(HCCONT,HCINT,FCINT,HCSINT,FCSINT,CSP1,CSP2,
     &                  IODDHC,IPNTUV,IPRINT)
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
     &          CSP1(KHKT1,KCKT1),
     &          CSP2(KHKT2,KCKT2),
     &          IODDHC(NRTOP),
     &          IPNTUV(KC2MAX,2,2)
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
               IODD12 = IODDHC(IPNTUV(ICMP12,1,1))
               DO 210 IKOMP2 = 1, KHKT2
                  SPHFAC = CSP2(IKOMP2,ICOMP2)
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
               SPHFAC = CSP1(IKOMP1,ICOMP1)
               IF (ABS(SPHFAC) .GT. D0) THEN
                  MAX2 = KHKT2
                  IF (TKMP12) MAX2 = IKOMP1
                  DO 310 IKOMP2 = 1, MAX2
                     IKMP12 = IKMP12 + 1
                     IODD12 = IODDHC(IPNTUV(IKMP12,2,1))
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
            SPHFAC = CSP1(IKOMP1,ICOMP1)
            IF (ABS(SPHFAC) .GT. D0) THEN
            DO 410 IKOMP2 = 1, KHKT2
               ICMP12 = (ICOMP1 - 1)*KHKT2 + IKOMP2
               IKMP12 = (IKOMP1 - 1)*KHKT2 + IKOMP2
               IODD12 = IODDHC(IPNTUV(IKMP12,2,1))
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
            SPHFAC = CSP2(IKOMP2,ICOMP2)
            IF (ABS(SPHFAC) .GT. D0) THEN
            DO 510 IKOMP1 = 1, KHKT1
               ICMP12 = (IKOMP1 - 1)*KCKT2 + ICOMP2
               IKMP12 = (IKOMP1 - 1)*KHKT2 + IKOMP2
               IODD12 = IODDHC(IPNTUV(IKMP12,2,1))
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
      SUBROUTINE CR2SPH(CCONT,CSINT,FCSINT,AOINT,FAOINT,CSP3,CSP4,
     &                  IODDCC,IPNTUV,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      LOGICAL   FCSINT(KHKT12,KHKT4),
     &          FAOINT(KHKT12,KHKT34)
      DIMENSION CCONT(NCCCC,KHKT12,KCKT34),
     &          CSINT(NCCCC,KHKT12,KHKT4),
     &          AOINT(NCCCC,KHKT12,KHKT34),
     &          CSP3(KHKT3,KCKT3), CSP4(KHKT4,KCKT4),
     &          IODDCC(NRTOP),
     &          IPNTUV(KC2MAX,2,2)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
C
      DO 10 IKMP12 = 1, KHKT12
      DO 10 IKMP34 = 1, KHKT34
         FAOINT(IKMP12,IKMP34) = .TRUE.
   10 CONTINUE
C
C     Transformation of both indices
C     ==============================
C
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
               IODD34 = IODDCC(IPNTUV(ICMP34,1,2))
               DO 210 IKOMP4 = 1, KHKT4
                  SPHFAC = CSP4(IKOMP4,ICOMP4)
                  IF (ABS(SPHFAC).GT.D0) THEN
                     DO 220 IKMP12 = 1, KHKT12
                     IF (IODDCC(IPNTUV(IKMP12,2,1)) .EQ. IODD34) THEN
                        IF (FCSINT(IKMP12,IKOMP4)) THEN
                           FCSINT(IKMP12,IKOMP4) = .FALSE.
                           DO 230 I = 1, NCCCC
                              CSINT(I,IKMP12,IKOMP4) =
     &                                SPHFAC*CCONT(I,IKMP12,ICMP34)
  230                      CONTINUE
                        ELSE
                           DO 235 I = 1, NCCCC
                              CSINT(I,IKMP12,IKOMP4) =
     &                              CSINT(I,IKMP12,IKOMP4)
     &                              + SPHFAC*CCONT(I,IKMP12,ICMP34)
  235                      CONTINUE
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
               SPHFAC = CSP3(IKOMP3,ICOMP3)
               IF (ABS(SPHFAC) .GT. D0) THEN
                  MAX4 = KHKT4
                  IF (TKMP34) MAX4 = IKOMP3
                  DO 310 IKOMP4 = 1, MAX4
                     IKMP34 = IKMP34 + 1
                     IODD34 = IODDCC(IPNTUV(IKMP34,2,2))
                     DO 320 IKMP12 = 1, KHKT12
                     IF (IODDCC(IPNTUV(IKMP12,2,1)) .EQ. IODD34) THEN
                        IF (FAOINT(IKMP12,IKMP34)) THEN
                           FAOINT(IKMP12,IKMP34) = .FALSE.
                           DO 330 I = 1, NCCCC
                              AOINT(I,IKMP12,IKMP34) =
     &                                SPHFAC*CSINT(I,IKMP12,IKOMP4)
  330                      CONTINUE
                        ELSE
                           DO 335 I = 1, NCCCC
                              AOINT(I,IKMP12,IKMP34) =
     &                              AOINT(I,IKMP12,IKMP34)
     &                              + SPHFAC*CSINT(I,IKMP12,IKOMP4)
  335                      CONTINUE
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
C     Transformatio of first index only
C     =================================
C
      ELSE IF (SPHR3) THEN
         DO 400 ICOMP3 = 1, KCKT3
         DO 400 IKOMP3 = 1, KHKT3
            SPHFAC = CSP3(IKOMP3,ICOMP3)
            IF (ABS(SPHFAC) .GT. D0) THEN
               DO 420 IKOMP4 = 1, KHKT4
                  ICMP34 = (ICOMP3 - 1)*KHKT4 + IKOMP4
                  IKMP34 = (IKOMP3 - 1)*KHKT4 + IKOMP4
                  IODD34 = IODDCC(IPNTUV(IKMP34,2,2))
                  DO 430 IKMP12 = 1, KHKT12
                  IF (IODDCC(IPNTUV(IKMP12,2,1)) .EQ. IODD34) THEN
                     IF (FAOINT(IKMP12,IKMP34)) THEN
                        FAOINT(IKMP12,IKMP34) = .FALSE.
                        DO 440 I = 1, NCCCC
                           AOINT(I,IKMP12,IKMP34) =
     &                                     SPHFAC*CCONT(I,IKMP12,ICMP34)
  440                   CONTINUE
                     ELSE
                        DO 445 I = 1, NCCCC
                           AOINT(I,IKMP12,IKMP34) =
     &                                            AOINT(I,IKMP12,IKMP34)
     &                                   + SPHFAC*CCONT(I,IKMP12,ICMP34)
  445                   CONTINUE
                     END IF
                  END IF
  430             CONTINUE
  420          CONTINUE
            END IF
  400    CONTINUE
C
C     Transformatio of second index only
C     ==================================
C
      ELSE
         DO 500 IKOMP3 = 1, KHKT3
         DO 500 ICOMP4 = 1, KCKT4
            ICMP34 = (IKOMP3 - 1)*KCKT4 + ICOMP4
               DO 510 IKOMP4 = 1, KHKT4
                  SPHFAC = CSP4(IKOMP4,ICOMP4)
                  IF (ABS(SPHFAC).GT.D0) THEN
                     IKMP34 = (IKOMP3 - 1)*KHKT4 + IKOMP4
                     IODD34 = IODDCC(IPNTUV(IKMP34,2,2))
                     DO 520 IKMP12 = 1, KHKT12
                     IF (IODDCC(IPNTUV(IKMP12,2,1)) .EQ. IODD34) THEN
                        IF (FAOINT(IKMP12,IKMP34)) THEN
                           FAOINT(IKMP12,IKMP34) = .FALSE.
                           DO 530 I = 1, NCCCC
                              AOINT(I,IKMP12,IKMP34) =
     &                                   SPHFAC*CCONT(I,IKMP12,ICMP34)
  530                      CONTINUE
                        ELSE
                           DO 535 I = 1, NCCCC
                              AOINT(I,IKMP12,IKMP34) =
     &                                          AOINT(I,IKMP12,IKMP34)
     &                                 + SPHFAC*CCONT(I,IKMP12,ICMP34)
  535                      CONTINUE
                        END IF
                     END IF
  520                CONTINUE
                  END IF
  510          CONTINUE
  500    CONTINUE
      END IF
      RETURN
      END
