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
                     DO TUV = 1, NTUV34
                     IF (IODDHC(TUV) .EQ. IODD12) THEN
                        IF (FCSINT(TUV,IKOMP2)) THEN
                           FCSINT(TUV,IKOMP2) = .FALSE.
                           DO I = 1, NCCPP
                              HCSINT(I,TUV,IKOMP2)=
     &                                      SPHFAC*HCCONT(I,TUV,ICMP12)
                           END DO
                        ELSE
                           DO I = 1, NCCPP
                              HCSINT(I,TUV,IKOMP2)=HCSINT(I,TUV,IKOMP2)
     &                                    + SPHFAC*HCCONT(I,TUV,ICMP12)
                           END DO
                        END IF
                     END IF
                     END DO
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
                     DO TUV = 1, NTUV34
                     IF (IODDHC(TUV) .EQ. IODD12) THEN
                        IF (FCINT(TUV,IKMP12)) THEN
                           FCINT(TUV,IKMP12) = .FALSE.
                           DO I = 1, NCCPP
                              HCINT(I,TUV,IKMP12) =
     &                                      SPHFAC*HCSINT(I,TUV,IKOMP2)
                           END DO
                        ELSE
                           DO I = 1, NCCPP
                              HCINT(I,TUV,IKMP12) = HCINT(I,TUV,IKMP12)
     &                                    + SPHFAC*HCSINT(I,TUV,IKOMP2)
                           END DO
                        END IF
                     END IF
                     END DO
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
               DO TUV = 1, NTUV34
               IF (IODDHC(TUV) .EQ. IODD12) THEN
                  IF (FCINT(TUV,IKMP12)) THEN
                     FCINT(TUV,IKMP12) = .FALSE.
                     DO I = 1, NCCPP
                        HCINT(I,TUV,IKMP12) =
     &                                SPHFAC*HCCONT(I,TUV,ICMP12)
                     END DO
                  ELSE
                     DO I = 1, NCCPP
                        HCINT(I,TUV,IKMP12) = HCINT(I,TUV,IKMP12)
     &                              + SPHFAC*HCCONT(I,TUV,ICMP12)
                     END DO
                  END IF
               END IF
               END DO
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
               DO TUV = 1, NTUV34
               IF (IODDHC(TUV) .EQ. IODD12) THEN
                  IF (FCINT(TUV,IKMP12)) THEN
                     FCINT(TUV,IKMP12) = .FALSE.
                     DO I = 1, NCCPP
                        HCINT(I,TUV,IKMP12) =
     &                                SPHFAC*HCCONT(I,TUV,ICMP12)
                     END DO
                  ELSE
                     DO I = 1, NCCPP
                        HCINT(I,TUV,IKMP12) = HCINT(I,TUV,IKMP12)
     &                              + SPHFAC*HCCONT(I,TUV,ICMP12)
                     END DO
                  END IF
               END IF
               END DO
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
     &          IODDCC(NRTOP), NCENTR(NPQBCX,2),
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          IODDKC(KC2MAX)
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (BDER3.OR.BDER4) THEN
         CALL CR2SVC(CSQ3,CSA3,CSV3,NCENTR(1,1),NCENTR(1,2),
     &               KHKT3,KCKT3)
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
      SUBROUTINE CR2SVC(CSQ,CSA,CSV,NCNT1,NCNT2,KHKTX,KCKTX)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION CSQ(NCSQ1,NCSQ2), CSA(KHKTX,KCKTX), 
     &          CSV(NCTFPQ,NPQBCX,KHKTX,KCKTX), 
     &          NCNT1(NPQBCX), NCNT2(NPQBCX)
#include <ericom.h>
#include <nuclei.h>
      IADR = 1
      DO ICMP = 1, KCKTX
      DO IKMP = 1, KHKTX
         CSQM = D0
         DO I = 1, NPQBCX
            CSQI = CSQ(IADR,NCNT2(I)) - CSQ(IADR,NCNT1(I))
            CSQM = MAX(CSQM,ABS(CSQI)) 
            DO J = 1, NCTFPQ
               CSV(J,I,IKMP,ICMP) = CSQI 
            END DO
         END DO
         CSA(IKMP,ICMP) = CSQM
         IADR = IADR + 1
      END DO
      END DO
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
c     IF (FLEVEL) THEN
         DO IKMP12 = 1, KHKT12
         DO IKMP34 = 1, KHKT34
c           FAOINT(IKMP12,IKMP34) = .TRUE.
            FAOINT(IKMP12,IKMP34) = .FALSE.
         END DO
         END DO
c     END IF
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
            DO IKMP12 = 1, KHKT12
            DO IKOMP4 = 1, KHKT4
               FCSINT(IKMP12,IKOMP4) = .TRUE.
            END DO
            END DO
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
                              DO I = 1, NCCCC
                                 CSINT(I,IKMP12,IKOMP4) =
     &                              SPHFAC*CCONT(I,IKMP12,ICMP34)
                              END DO
                           ELSE
                              DO I = 1, NCCCC
                                 CSINT(I,IKMP12,IKOMP4) = 
     &                               CSV4(I,IKOMP4,ICOMP4)
     &                                   *CCONT(I,IKMP12,ICMP34)
                              END DO
                           END IF
                        ELSE
                           IF (.NOT.BDER4) THEN
                              DO I = 1, NCCCC
                                 CSINT(I,IKMP12,IKOMP4) =
     &                                 CSINT(I,IKMP12,IKOMP4)
     &                                 + SPHFAC*CCONT(I,IKMP12,ICMP34)
                              END DO
                           ELSE
                              DO I = 1, NCCCC
                                 CSINT(I,IKMP12,IKOMP4) =
     &                                 CSINT(I,IKMP12,IKOMP4)
     &                                    + CSV4(I,IKOMP4,ICOMP4)
     &                                         *CCONT(I,IKMP12,ICMP34)
                              END DO
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
                              DO I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                   SPHFAC*CSINT(I,IKMP12,IKOMP4)
                              END DO
                           ELSE
                              DO I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                              CSV3(I,IKOMP3,ICOMP3)
     &                                 *CSINT(I,IKMP12,IKOMP4)
                              END DO
                           END IF
                        ELSE
                           IF (.NOT.BDER3) THEN
                              DO I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                 AOINT(I,IKMPAB,IKMPCD)
     &                                 + SPHFAC*CSINT(I,IKMP12,IKOMP4)
                              END DO
                           ELSE
                              DO I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                              AOINT(I,IKMPAB,IKMPCD)
     &                                 + CSV3(I,IKOMP3,ICOMP3)
     &                                      *CSINT(I,IKMP12,IKOMP4)
                              END DO
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
                           DO I = 1, NCCCC
                              AOINT(I,IKMPAB,IKMPCD) =
     &                                     SPHFAC*CCONT(I,IKMP12,ICMP34)
                           END DO
                        ELSE
                           DO I = 1, NCCCC
                              AOINT(I,IKMPAB,IKMPCD) =
     &                           CSV3(I,IKOMP3,ICOMP3)
     &                              *CCONT(I,IKMP12,ICMP34)
                           END DO
                        END IF
                     ELSE
                        IF (.NOT.BDER3) THEN
                           DO I = 1, NCCCC
                              AOINT(I,IKMPAB,IKMPCD) =
     &                                            AOINT(I,IKMPAB,IKMPCD)
     &                                   + SPHFAC*CCONT(I,IKMP12,ICMP34)
                           END DO
                        ELSE
                           DO I = 1, NCCCC
                              AOINT(I,IKMPAB,IKMPCD) =
     &                           AOINT(I,IKMPAB,IKMPCD)
     &                              + CSV3(I,IKOMP3,ICOMP3)
     &                                 *CCONT(I,IKMP12,ICMP34)
                           END DO
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
                              DO I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                   SPHFAC*CCONT(I,IKMP12,ICMP34)
                              END DO
                           ELSE
                              DO I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                              CSV4(I,IKOMP4,ICOMP4)
     &                                 *CCONT(I,IKMP12,ICMP34)
                              END DO
                           END IF 
                        ELSE
                           IF (.NOT.BDER4) THEN
                              DO I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                          AOINT(I,IKMPAB,IKMPCD)
     &                                 + SPHFAC*CCONT(I,IKMP12,ICMP34)
                              END DO
                           ELSE
                              DO I = 1, NCCCC
                                 AOINT(I,IKMPAB,IKMPCD) =
     &                                          AOINT(I,IKMPAB,IKMPCD)
     &                              + CSV4(I,IKOMP4,ICOMP4)
     &                                  *CCONT(I,IKMP12,ICMP34)
                              END DO
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
      SUBROUTINE ERIBUN(AOINT,AOINT0,COORAO,IODDCC,IPNTUV,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.5D0)
      DIMENSION AOINT(NCCT,KHKTAB,KHKTCD,6), AOINT0(NCCT,KHKTAB,KHKTCD), 
     &          COORAO(NPQBCX,3,4),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2)

#include <ericom.h>
#include <orgcom.h>
#include <eriao.h>
#include <hertop.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from ERIBUN','*',103)
C
      DO IAB = 1, KHKTAB 
      DO ICD = 1, KHKTCD 
      IF (IODDCC(IPNTUV(IAB,0,1)).EQ.IODDCC(IPNTUV(ICD,0,2))) THEN
         IF (GCONAB .OR. GCONCD) THEN
            INT1 = 0
            DO N = 1, NPQBCX
               AX = COORAO(N,1,1) - ORIGIN(1) 
               AY = COORAO(N,2,1) - ORIGIN(2) 
               AZ = COORAO(N,3,1) - ORIGIN(3) 
               BX = COORAO(N,1,2) - ORIGIN(1)
               BY = COORAO(N,2,2) - ORIGIN(2)
               BZ = COORAO(N,3,2) - ORIGIN(3)
               CX = COORAO(N,1,3) - ORIGIN(1) 
               CY = COORAO(N,2,3) - ORIGIN(2) 
               CZ = COORAO(N,3,3) - ORIGIN(3) 
               DX = COORAO(N,1,4) - ORIGIN(1)
               DY = COORAO(N,2,4) - ORIGIN(2)
               DZ = COORAO(N,3,4) - ORIGIN(3)
               ABX = AY*BZ - AZ*BY
               ABY = AZ*BX - AX*BZ
               ABZ = AX*BY - AY*BX
               CDX = CY*DZ - CZ*DY
               CDY = CZ*DX - CX*DZ
               CDZ = CX*DY - CY*DX
               DO I = 1, NCTFA
               DO J = 1, NCTFB
               DO K = 1, NCTFC
               DO L = 1, NCTFD
                  INT1 = INT1 + 1
                  INT2 = (I - 1)*NCTFB + J + (K - 1)*NCTFAB*NCTFD
     &                 + (L - 1)*NCTFAB + (N - 1)*NCTFPQ
                  AO = DP5*AOINT0(INT1,IAB,ICD)
                  AOINT(INT2,IAB,ICD,1) = ABX*AO 
                  AOINT(INT2,IAB,ICD,2) = ABY*AO 
                  AOINT(INT2,IAB,ICD,3) = ABZ*AO 
                  AOINT(INT1,IAB,ICD,4) = CDX*AO 
                  AOINT(INT1,IAB,ICD,5) = CDY*AO 
                  AOINT(INT1,IAB,ICD,6) = CDZ*AO 
               END DO
               END DO
               END DO
               END DO
            END DO
          ELSE
            DO I = 1, NCCX
               AO = DP5*AOINT0(I,IAB,ICD)
               AX = COORAO(I,1,1) - ORIGIN(1) 
               AY = COORAO(I,2,1) - ORIGIN(2) 
               AZ = COORAO(I,3,1) - ORIGIN(3) 
               BX = COORAO(I,1,2) - ORIGIN(1)
               BY = COORAO(I,2,2) - ORIGIN(2)
               BZ = COORAO(I,3,2) - ORIGIN(3)
               CX = COORAO(I,1,3) - ORIGIN(1) 
               CY = COORAO(I,2,3) - ORIGIN(2) 
               CZ = COORAO(I,3,3) - ORIGIN(3) 
               DX = COORAO(I,1,4) - ORIGIN(1)
               DY = COORAO(I,2,4) - ORIGIN(2)
               DZ = COORAO(I,3,4) - ORIGIN(3)
               AOINT(I,IAB,ICD,1) = (AY*BZ - AZ*BY)*AO 
               AOINT(I,IAB,ICD,2) = (AZ*BX - AX*BZ)*AO 
               AOINT(I,IAB,ICD,3) = (AX*BY - AY*BX)*AO 
               AOINT(I,IAB,ICD,4) = (CY*DZ - CZ*DY)*AO 
               AOINT(I,IAB,ICD,5) = (CZ*DX - CX*DZ)*AO 
               AOINT(I,IAB,ICD,6) = (CX*DY - CY*DX)*AO 
            END DO
          END IF
      END IF
      END DO
      END DO
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('Final spherical integrals - ERIBUN',-1)
         DO K = 1, 6
            DO I = 1, KHKTAB
            DO J = 1, KHKTCD
            IF (IODDCC(IPNTUV(I,0,1)) .EQ. 
     &          IODDCC(IPNTUV(J,0,2))) THEN
               WRITE (LUPRI, '(/1X,A,3I3)') ' K,ICMPAB, ICMPCD ', K,I,J
               CALL OUTPUT(AOINT(1,I,J,K),1,1,1,NCCCC,1,NCCCC,1,LUPRI)
            END IF
            END DO
            END DO
         END DO
      END IF
      RETURN
      END
